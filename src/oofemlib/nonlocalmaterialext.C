/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "gausspoint.h"
#include "floatarray.h"
#include "element.h"
#include "timestep.h"
#include "integrationrule.h"
#include "nonlocalmaterialext.h"
#include "material.h"
#include "spatiallocalizer.h"
#include "domain.h"
#include "nonlocalbarrier.h"
#include "mathfem.h"
#include "dynamicinputrecord.h"

#ifdef __PARALLEL_MODE
 #include "parallel.h"
#endif

#include <list>

namespace oofem {
// flag forcing the inclusion of all elements with volume inside support of weight function.
// This forces inclusion of all integration points of these elements, even if weight is zero
// If not defined (default) only integration points with nonzero weight are included.
// #define NMEI_USE_ALL_ELEMENTS_IN_SUPPORT


// constructor
NonlocalMaterialExtensionInterface :: NonlocalMaterialExtensionInterface(Domain *d)  : Interface()
{
    domain = d;
    regionMap.resize( d->giveNumberOfRegions() ); /*lastUpdatedStateCounter = 0;*/
    if ( this->hasBoundedSupport() ) {
        permanentNonlocTableFlag = true;
    } else {
        permanentNonlocTableFlag = false;
    }

    cl = 0.;
    suprad = 0.;
    mm = 0.;
    weightFun = WFT_Unknown;
    scaling = ST_Unknown;
    averagedVar = AVT_Unknown;

    cl0 = 0.;
    beta = 0.;
    zeta = 0.;
    nlvar = NLVT_Standard;

    px = 0.;
}

void
NonlocalMaterialExtensionInterface :: updateDomainBeforeNonlocAverage(TimeStep *tStep)
{
    Domain *d = this->giveDomain();

    if ( d->giveNonlocalUpdateStateCounter() == tStep->giveSolutionStateCounter() ) {
        return; // already updated
    }

    OOFEM_LOG_DEBUG ("Updating Before NonlocAverage\n");
    for ( auto &elem : d->giveElements() ) {
        elem->updateBeforeNonlocalAverage(tStep);
    }

    // mark last update counter to prevent multiple updates
    d->setNonlocalUpdateStateCounter( tStep->giveSolutionStateCounter() );
}

void
NonlocalMaterialExtensionInterface :: buildNonlocalPointTable(GaussPoint *gp)
{
    double elemVolume, integrationVolume = 0.;

    NonlocalMaterialStatusExtensionInterface *statusExt =
        static_cast< NonlocalMaterialStatusExtensionInterface * >( gp->giveMaterialStatus()->
                                                                  giveInterface(NonlocalMaterialStatusExtensionInterfaceType) );
    if ( !statusExt ) {
        OOFEM_ERROR("local material status encountered");
    }

    if ( !statusExt->giveIntegrationDomainList()->empty() ) {
        return;                                                  // already done
    }

    auto iList = statusExt->giveIntegrationDomainList();

    FloatArray gpCoords, jGpCoords, shiftedGpCoords;
    if ( gp->giveElement()->computeGlobalCoordinates( gpCoords, gp->giveNaturalCoordinates() ) == 0 ) {
        OOFEM_ERROR("computeGlobalCoordinates of target failed");
    }

    // If nonlocal variation is set to the distance-based approach, a new nonlocal radius
    // is calculated as a function of the distance from the Gauss point to the nonlocal boundaries
    if ( nlvar == NLVT_DistanceBasedLinear || nlvar == NLVT_DistanceBasedExponential ) {
        //      cl=cl0;
        cl = giveDistanceBasedInteractionRadius(gpCoords);
        suprad = evaluateSupportRadius();
    }

    // If the mesh represents a periodic cell, nonlocal interaction is considered not only for the real neighbors
    // but also for their periodic images, shifted by +px or -px in the x-direction. In the implementation,
    // instead of shifting the potential neighbors, we shift the receiver point gp. In the non-periodic case (typical),
    // px=0 and the following loop is executed only once.  

    int nx = 0; // typical case
    if ( px > 0. ) nx = 1; // periodicity taken into account

    for ( int ix = -nx; ix <= nx; ix++ ) { // loop over periodic images shifted in x-direction
        SpatialLocalizer :: elementContainerType elemSet;
        shiftedGpCoords = gpCoords;
        shiftedGpCoords.at(1) += ix*px;

    // ask domain spatial localizer for list of elements with IP within this zone
#ifdef NMEI_USE_ALL_ELEMENTS_IN_SUPPORT
        this->giveDomain()->giveSpatialLocalizer()->giveAllElementsWithNodesWithinBox(elemSet, shiftedGpCoords, suprad);
        // insert element containing given gp
        elemSet.insert( gp->giveElement()->giveNumber() );
#else
        this->giveDomain()->giveSpatialLocalizer()->giveAllElementsWithIpWithinBox_EvenIfEmpty(elemSet, shiftedGpCoords, suprad);
#endif
        // initialize iList
        iList->reserve(elemSet.giveSize());
        for ( auto elindx: elemSet ) {
            Element *ielem = this->giveDomain()->giveElement(elindx);
            if ( regionMap.at( ielem->giveRegionNumber() ) == 0 ) {
                for ( auto &jGp: *ielem->giveDefaultIntegrationRulePtr() ) {
                    if ( ielem->computeGlobalCoordinates( jGpCoords, jGp->giveNaturalCoordinates() ) ) {
                        double weight = this->computeWeightFunction(shiftedGpCoords, jGpCoords);

                        //manipulate weights for a special averaging of strain (OFF by default)
                        this->manipulateWeight(weight, gp, jGp);

                        this->applyBarrierConstraints(shiftedGpCoords, jGpCoords, weight);
#ifdef NMEI_USE_ALL_ELEMENTS_IN_SUPPORT
                        if ( 1 ) {
#else
                        if ( weight > 0. ) {
#endif
                            localIntegrationRecord ir;
                            ir.nearGp = jGp;  // store gp
                            elemVolume = weight * jGp->giveElement()->computeVolumeAround(jGp);
                            ir.weight = elemVolume; // store gp weight
                            iList->push_back(ir); // store own copy in list
                            integrationVolume += elemVolume;
                        }
                    } else {
                        OOFEM_ERROR("computeGlobalCoordinates of target failed");
                    }
                }
            }
        } // loop over elements
        iList->shrink_to_fit();
    }

    statusExt->setIntegrationScale(integrationVolume); // store scaling factor
}

void
NonlocalMaterialExtensionInterface :: rebuildNonlocalPointTable(GaussPoint *gp, IntArray *contributingElems)
{
    double weight, elemVolume, integrationVolume = 0.;

    NonlocalMaterialStatusExtensionInterface *statusExt =
        static_cast< NonlocalMaterialStatusExtensionInterface * >( gp->giveMaterialStatus()->
                                                                  giveInterface(NonlocalMaterialStatusExtensionInterfaceType) );

    if ( !statusExt ) {
        OOFEM_ERROR("local material status encountered");
    }

    auto iList = statusExt->giveIntegrationDomainList();
    iList->clear();

    if ( contributingElems == NULL ) {
        // no element table provided, use standard method
        this->buildNonlocalPointTable(gp);
    } else {
        FloatArray gpCoords, jGpCoords;
        int _size = contributingElems->giveSize();
        if ( gp->giveElement()->computeGlobalCoordinates( gpCoords, gp->giveNaturalCoordinates() ) == 0 ) {
            OOFEM_ERROR("computeGlobalCoordinates of target failed");
        }

        //If nonlocal variation is set to the distance-based approach calculates  new nonlocal radius
        // based on the distance from the nonlocal boundaries
        if ( nlvar == NLVT_DistanceBasedLinear || nlvar == NLVT_DistanceBasedExponential ) {
            cl = cl0;
            cl = giveDistanceBasedInteractionRadius(gpCoords);
            suprad = evaluateSupportRadius();
        }

        // initialize iList
        iList->reserve(_size);
        for ( int _e = 1; _e <= _size; _e++ ) {
            Element *ielem = this->giveDomain()->giveElement( contributingElems->at(_e) );
            if ( regionMap.at( ielem->giveRegionNumber() ) == 0 ) {
                for ( auto &jGp:* ielem->giveDefaultIntegrationRulePtr() ) {
                    if ( ielem->computeGlobalCoordinates( jGpCoords, jGp->giveNaturalCoordinates() ) ) {
                        weight = this->computeWeightFunction(gpCoords, jGpCoords);

                        //manipulate weights for a special averaging of strain (OFF by default)
                        this->manipulateWeight(weight, gp, jGp);

                        this->applyBarrierConstraints(gpCoords, jGpCoords, weight);
#ifdef NMEI_USE_ALL_ELEMENTS_IN_SUPPORT
                        if ( 1 ) {
#else
                        if ( weight > 0. ) {
#endif
                            localIntegrationRecord ir;
                            ir.nearGp = jGp;     // store gp
                            elemVolume = weight * jGp->giveElement()->computeVolumeAround(jGp);
                            ir.weight = elemVolume; // store gp weight
                            iList->push_back(ir); // store own copy in list
                            integrationVolume += elemVolume;
                        }
                    } else {
                        OOFEM_ERROR("computeGlobalCoordinates of target failed");
                    }
                }
            }
        } // loop over elements

        statusExt->setIntegrationScale(integrationVolume); // remember scaling factor
#ifdef __PARALLEL_MODE
 #ifdef __VERBOSE_PARALLEL
        fprintf( stderr, "%d(%d):", gp->giveElement()->giveGlobalNumber(), gp->giveNumber() );
        for ( auto &lir: *iList ) {
            fprintf(stderr, "%d,%d(%e)", lir.nearGp->giveElement()->giveGlobalNumber(), lir.nearGp->giveNumber(), lir.weight);
        }

        fprintf(stderr, "\n");
 #endif
#endif
    }
}


std :: vector< localIntegrationRecord > *
NonlocalMaterialExtensionInterface :: giveIPIntegrationList(GaussPoint *gp)
{
    NonlocalMaterialStatusExtensionInterface *statusExt =
        static_cast< NonlocalMaterialStatusExtensionInterface * >( gp->giveMaterialStatus()->
                                                                  giveInterface(NonlocalMaterialStatusExtensionInterfaceType) );

    if ( !statusExt ) {
        OOFEM_ERROR("local material status encountered");
    }

    if ( statusExt->giveIntegrationDomainList()->empty() ) {
        this->buildNonlocalPointTable(gp);
    }

    return statusExt->giveIntegrationDomainList();
}

void
NonlocalMaterialExtensionInterface :: endIPNonlocalAverage(GaussPoint *gp)
{
    NonlocalMaterialStatusExtensionInterface *statusExt =
        static_cast< NonlocalMaterialStatusExtensionInterface * >( gp->giveMaterialStatus()->
                                                                  giveInterface(NonlocalMaterialStatusExtensionInterfaceType) );

    if ( !statusExt ) {
        OOFEM_ERROR("local material status encountered");
    }

    if ( ( !this->hasBoundedSupport() ) || ( !permanentNonlocTableFlag ) ) {
        statusExt->clear();
    }
}

double
NonlocalMaterialExtensionInterface :: computeWeightFunction(double distance)
{
    if ( weightFun == WFT_UniformOverElement ) { // uniform function over one element
        return 1.;
    }

    if ( distance > suprad || distance < 0. ) {
        return 0.;
    }

    double aux = distance / this->cl;
    double iwf = giveIntegralOfWeightFunction( this->domain->giveNumberOfSpatialDimensions() );

    switch ( weightFun ) {
    case WFT_Bell: // Bell shaped function (quartic spline)
        aux = ( 1. - aux * aux );
        return aux * aux / iwf;

    case WFT_Gauss: // Gauss function
        return exp(-aux * aux) / iwf;

    case WFT_Green: // Function corresponding in 1D to Green's function of Helmholtz equation (implicit gradient model)
        //printf("%14g %14g\n",distance,exp(-aux) / iwf);
        return exp(-aux) / iwf;

    case WFT_Green_21: // Green function reduced from 2D to 1D
    {
        /*
         * if (this->domain->giveNumberOfSpatialDimensions() != 1){
         * OOFEM_ERROR("this type of weight function can be used for a 1D problem only");
         * }
         */
        iwf = giveIntegralOfWeightFunction(2); // indeed
        double x = distance;
        double y = 0.;
        double r = sqrt(x * x + y * y);
        double sum = exp(-r / this->cl);
        double h = this->cl / 10.; // 10 could later be replaced by an optional parameter
        do {
            y += h;
            r = sqrt(x * x + y * y);
            sum += 2. * exp(-r / this->cl);
        } while ( r <= suprad );
        //printf("%14g %14g\n",distance,sum * h / iwf);
        return sum * h / iwf;
    }

    case WFT_Uniform: // uniform function over an interaction distance
        return 1. / iwf;

    default:
        OOFEM_WARNING("unknown type of weight function %d", weightFun);
        return 0.0;
    }
}

double
NonlocalMaterialExtensionInterface :: computeWeightFunction(const FloatArray &src, const FloatArray &coord)
{
    return computeWeightFunction( src.distance(coord) );
}

double
NonlocalMaterialExtensionInterface :: giveIntegralOfWeightFunction(const int spatial_dimension)
{
    const double pi = M_PI;
    switch ( weightFun ) {
    case WFT_Bell:
        switch ( spatial_dimension ) {
        case 1: return cl * 16. / 15.;

        case 2: return cl * cl * pi / 3.;

        case 3: return cl * cl * cl * pi * 8. / 105.;

        default: return 1.;
        }

    case WFT_Gauss:
        switch ( spatial_dimension ) {
        case 1: return cl *sqrt(pi);

        case 2: return cl * cl * pi;

        case 3: return cl *cl *cl *sqrt(pi *pi *pi) / 4.;

        default: return 1.;
        }

    case WFT_Green:
    case WFT_Green_21:
        switch ( spatial_dimension ) {
        case 1: return cl * 2.;

        case 2: return cl * cl * 2. * pi;

        case 3: return cl * cl * cl * 8. * pi;

        default: return 1.;
        }

    case WFT_Uniform:
        switch ( spatial_dimension ) {
        case 1: return cl * 2.;

        case 2: return cl * cl * pi;

        case 3: return cl * cl * cl * ( 4. / 3. ) * pi;

        default: return 1.;
        }

    default: return 1.;
    }
}

double
NonlocalMaterialExtensionInterface :: maxValueOfWeightFunction()
{
    double iwf = giveIntegralOfWeightFunction( this->domain->giveNumberOfSpatialDimensions() );
    return 1. / iwf;
}

double
NonlocalMaterialExtensionInterface :: evaluateSupportRadius()
{
    switch ( weightFun ) {
    case WFT_Bell:  return cl;

    case WFT_Gauss: return 2.5 * cl;

    case WFT_Green: case WFT_Green_21: return 6. * cl;

    case WFT_Uniform:  return cl;

    case WFT_UniformOverElement:  return 0.0; // to make sure that only Gauss points of the same element will be considered as neighbors

    default:        return cl;
    }
}


IRResultType
NonlocalMaterialExtensionInterface :: initializeFrom(InputRecord *ir)
{
    IRResultType result;              // Required by IR_GIVE_FIELD macro

    if ( ir->hasField(_IFT_NonlocalMaterialExtensionInterface_regionmap) ) {
        IR_GIVE_FIELD(ir, regionMap, _IFT_NonlocalMaterialExtensionInterface_regionmap);
        if ( regionMap.giveSize() != this->giveDomain()->giveNumberOfRegions() ) {
            OOFEM_ERROR("regionMap size mismatch");
        }
    } else {
        regionMap.zero();
    }

    if ( this->hasBoundedSupport() ) {
        permanentNonlocTableFlag = true;
    } else {
        permanentNonlocTableFlag = false;
    }
    IR_GIVE_OPTIONAL_FIELD(ir, this->permanentNonlocTableFlag, _IFT_NonlocalMaterialExtensionInterface_permanentNonlocTableFlag);

    // read the characteristic length
    IR_GIVE_FIELD(ir, cl, _IFT_NonlocalMaterialExtensionInterface_r);
    if ( cl < 0.0 ) {
        cl = 0.0;
    }

    // read the type of weight function
    int val = WFT_Bell;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_NonlocalMaterialExtensionInterface_wft);
    this->weightFun = ( WeightFunctionType ) val;

    // this is introduced for compatibility of input format with previous versions
    // ("averagingtype 1" in the input means that the weight function
    //   should be uniform over an element,
    //   which can now be described as "wft 5")
    int avt = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, avt, _IFT_NonlocalMaterialExtensionInterface_averagingtype);
    if ( avt == 1 ) {
        weightFun = WFT_UniformOverElement;
    }

    // evaluate the support radius based on type of weight function and characteristic length
    suprad = evaluateSupportRadius();

    // read the optional parameter for overnonlocal formulation
    mm = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, mm, _IFT_NonlocalMaterialExtensionInterface_m);

    // read the type of scaling
    val = ST_Standard;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_NonlocalMaterialExtensionInterface_scalingtype);
    this->scaling = ( ScalingType ) val;

    // read the type of averaged variable
    val = AVT_EqStrain;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_NonlocalMaterialExtensionInterface_averagedquantity);
    this->averagedVar = ( AveragedVarType ) val;

    // read the nonlocal variation type (default is zero)
    cl0 = cl;
    int nlvariation = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, nlvariation, _IFT_NonlocalMaterialExtensionInterface_nonlocalvariation);
    if ( nlvariation == 1 ) {
        nlvar = NLVT_DistanceBasedLinear;
        IR_GIVE_FIELD(ir, beta, _IFT_NonlocalMaterialExtensionInterface_beta);
        IR_GIVE_FIELD(ir, zeta, _IFT_NonlocalMaterialExtensionInterface_zeta);
    } else if ( nlvariation == 2 ) {
        nlvar = NLVT_StressBased;
        IR_GIVE_FIELD(ir, beta, _IFT_NonlocalMaterialExtensionInterface_beta);
    } else if ( nlvariation == 3 ) {
        nlvar = NLVT_DistanceBasedExponential;
        IR_GIVE_FIELD(ir, beta, _IFT_NonlocalMaterialExtensionInterface_beta);
        IR_GIVE_FIELD(ir, zeta, _IFT_NonlocalMaterialExtensionInterface_zeta);
    }

    // read the periodic shift (default is zero)
    px = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, px, _IFT_NonlocalMaterialExtensionInterface_px);

    return IRRT_OK;
}


void NonlocalMaterialExtensionInterface :: giveInputRecord(DynamicInputRecord &input)
{
    if ( regionMap.giveSize() ) {
        input.setField(this->regionMap, _IFT_NonlocalMaterialExtensionInterface_regionmap);
    }
    input.setField(this->permanentNonlocTableFlag, _IFT_NonlocalMaterialExtensionInterface_permanentNonlocTableFlag);
    input.setField(this->cl, _IFT_NonlocalMaterialExtensionInterface_r);
    input.setField(this->weightFun, _IFT_NonlocalMaterialExtensionInterface_wft);
    input.setField(this->mm, _IFT_NonlocalMaterialExtensionInterface_m);
    input.setField(this->scaling, _IFT_NonlocalMaterialExtensionInterface_scalingtype);
    input.setField(this->averagedVar, _IFT_NonlocalMaterialExtensionInterface_averagedquantity);
    input.setField(this->nlvar, _IFT_NonlocalMaterialExtensionInterface_nonlocalvariation);

    if ( nlvar == NLVT_DistanceBasedLinear ||  nlvar == NLVT_DistanceBasedExponential ) {
        input.setField(this->beta, _IFT_NonlocalMaterialExtensionInterface_beta);
        input.setField(this->zeta, _IFT_NonlocalMaterialExtensionInterface_zeta);
    } else if ( nlvar == NLVT_StressBased ) {
        input.setField(this->beta, _IFT_NonlocalMaterialExtensionInterface_beta);
    }
}



//int NonlocalMaterialExtension :: giveElementRegion (Element* element) {return element->giveCrossSection()->giveNumber();}
//int NonlocalMaterialExtension :: giveNumberOfRegions () {return domain->giveNumberOfCrossSectionModels();}

/*
 * bool
 * NonlocalMaterialExtensionInterface::isBarrierActivated (const FloatArray& c1, const FloatArray& c2) const
 * {
 * int ib, nbarrier = domain->giveNumberOfNonlocalBarriers();
 *
 * for (ib=1; ib<=nbarrier; ib++) {
 * if (domain->giveNonlocalBarrier(ib)->isActivated(c1,c2))
 * return true;
 * }
 *
 * return false;
 * }
 */

void
NonlocalMaterialExtensionInterface :: applyBarrierConstraints(const FloatArray &gpCoords, const FloatArray &jGpCoords, double &weight)
{
    int ib, nbarrier = domain->giveNumberOfNonlocalBarriers();
    bool shieldFlag = false;

    for ( ib = 1; ib <= nbarrier; ib++ ) {
        domain->giveNonlocalBarrier(ib)->applyConstraint(gpCoords, jGpCoords, weight, shieldFlag, this);
        if ( shieldFlag ) {
            weight = 0.0;
            return;
        }
    }
}

void
NonlocalMaterialExtensionInterface :: manipulateWeight(double &weight, GaussPoint *gp, GaussPoint *jGp)
{
    Element *ielem = jGp->giveElement();
    IntegrationRule *iRule = ielem->giveDefaultIntegrationRulePtr();

    if ( ielem->giveMaterial()->hasProperty(AVERAGING_TYPE, jGp) ) {
        if ( ielem->giveMaterial()->give(AVERAGING_TYPE, jGp) == 1 ) {
            weight = 1. / ( iRule->giveNumberOfIntegrationPoints() ); //assign the same weights over the whole element
        }
    }
}


double
NonlocalMaterialExtensionInterface :: giveDistanceBasedInteractionRadius(const FloatArray &gpCoords)
{
    double distance = 1.e10; // Initially distance from the boundary is set to the maximum value
    double temp;
    int ib, nbarrier = domain->giveNumberOfNonlocalBarriers();
    for ( ib = 1; ib <= nbarrier; ib++ ) { //Loop over all the nonlocal barriers to find minimum distance from the boundary
        temp = domain->giveNonlocalBarrier(ib)->calculateMinimumDistanceFromBoundary(gpCoords);
        if ( distance > temp ) { //Check to find minimum distance from boundary from all nonlocal boundaries
            distance = temp;
        }
    }

    //Calculate interaction radius based on the minimum distance from the nonlocal boundaries
    double newradius = 0.0;
    if ( nlvar == NLVT_DistanceBasedLinear ) {
        if ( distance < zeta * cl0 ) {
            newradius = ( 1. - beta ) / ( zeta * cl0 ) * distance + beta;
        } else   {
            newradius = 1.;
        }
    } else if ( nlvar == NLVT_DistanceBasedExponential ) {
        newradius = 1. - ( 1. - beta ) * exp( -distance / ( zeta * cl0 ) );
    }
    return newradius * cl0;
}
////////////////////////////////////////////////////////////////////////////////////////////////////

NonlocalMaterialStatusExtensionInterface :: NonlocalMaterialStatusExtensionInterface() : Interface(), integrationDomainList()
{
    integrationScale = 0.;
}

NonlocalMaterialStatusExtensionInterface :: ~NonlocalMaterialStatusExtensionInterface()
{
    ;
}
} // end namespace oofem
