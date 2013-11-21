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
}

void
NonlocalMaterialExtensionInterface :: updateDomainBeforeNonlocAverage(TimeStep *atTime)
{
    int nelem, i;
    Domain *d = this->giveDomain();
    nelem = d->giveNumberOfElements();

    if ( d->giveNonlocalUpdateStateCounter() == atTime->giveSolutionStateCounter() ) {
        return; // already updated
    }

    for ( i = 1; i <= nelem; i++ ) {
        d->giveElement(i)->updateBeforeNonlocalAverage(atTime);
    }

    // mark last update counter to prevent multiple updates
    d->setNonlocalUpdateStateCounter( atTime->giveSolutionStateCounter() );
}

void
NonlocalMaterialExtensionInterface :: buildNonlocalPointTable(GaussPoint *gp)
{
    double weight, elemVolume, integrationVolume = 0.;

    NonlocalMaterialStatusExtensionInterface *statusExt =
        static_cast< NonlocalMaterialStatusExtensionInterface * >( gp->giveMaterialStatus()->
        giveInterface(NonlocalMaterialStatusExtensionInterfaceType) );
    std::list< localIntegrationRecord > *iList;

    Element *ielem;
    GaussPoint *jGp;
    IntegrationRule *iRule;

    if ( !statusExt ) {
        OOFEM_ERROR("NonlocalMaterialExtensionInterface::buildNonlocalPointTable : local material status encountered");
    }

    // test for bounded support - if no bounded support, the nonlocal point table is
    // big vasting of space.
    /*
     * if ( !this->hasBoundedSupport() ) {
     *  return;
     * }
     */

    if ( !statusExt->giveIntegrationDomainList()->empty() ) {
        return;                                                  // already done
    }

    iList = statusExt->giveIntegrationDomainList();

    FloatArray gpCoords, jGpCoords;
    SpatialLocalizer :: elementContainerType elemSet;
    if ( gp->giveElement()->computeGlobalCoordinates( gpCoords, * ( gp->giveCoordinates() ) ) == 0 ) {
        OOFEM_ERROR("NonlocalMaterialExtensionInterface::buildNonlocalPointTable: computeGlobalCoordinates of target failed");
    }

    //If nonlocal variation is set to the distance-based approach, a new nonlocal radius
    // is calculated as a function of the distance from the gausspoint to the nonlocal boundaries
    if ( nlvar == NLVT_DistanceBased ) {
        //      cl=cl0;
        cl = giveDistanceBasedInteractionRadius(gpCoords);
        suprad = evaluateSupportRadius();
    }


    // ask domain spatial localizer for list of elements with IP within this zone
#ifdef NMEI_USE_ALL_ELEMENTS_IN_SUPPORT
    this->giveDomain()->giveSpatialLocalizer()->giveAllElementsWithNodesWithinBox(elemSet, gpCoords, suprad);
    // insert element containing given gp
    elemSet.insert( gp->giveElement()->giveNumber() );
#else
    this->giveDomain()->giveSpatialLocalizer()->giveAllElementsWithIpWithinBox(elemSet, gpCoords, suprad);
#endif
    // initialize iList


    SpatialLocalizer :: elementContainerType :: iterator pos;

    for ( pos = elemSet.begin(); pos !=  elemSet.end(); ++pos ) {
        ielem = this->giveDomain()->giveElement(* pos);
        if ( regionMap.at( ielem->giveRegionNumber() ) == 0 ) {
            iRule = ielem->giveDefaultIntegrationRulePtr();
            for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
                jGp = iRule->getIntegrationPoint(j);
                if ( ielem->computeGlobalCoordinates( jGpCoords, * ( jGp->giveCoordinates() ) ) ) {
                    weight = this->computeWeightFunction(gpCoords, jGpCoords);

                    //manipulate weights for a special averaging of strain (OFF by default)
                    this->manipulateWeight(weight, gp, jGp);

                    /*
                     * if ((weight > NonlocalMaterialZeroWeight) && (!this->isBarrierActivated(gpCoords, jGpCoords))) {
                     */
                    this->applyBarrierConstraints(gpCoords, jGpCoords, weight);
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
                    OOFEM_ERROR("NonlocalMaterialExtensionInterface::buildNonlocalPointTable: computeGlobalCoordinates of target failed");
                }
            }
        }
    } // loop over elements

    statusExt->setIntegrationScale(integrationVolume); // remember scaling factor

    /*
     * // Old implementation without spatial localizer
     *
     * FloatArray jGpCoords;
     * for (i=1; i<=nelem; i++) {
     * ielem = this->giveDomain()->giveElement(i);
     * if (regionMap.at(ielem->giveRegionNumber()) == 0) {
     * iRule = ielem->giveDefaultIntegrationRulePtr ();
     * for (j=0 ; j < iRule->giveNumberOfIntegrationPoints() ; j++) {
     * jGp = iRule->getIntegrationPoint(j) ;
     * if (ielem->computeGlobalCoordinates (jGpCoords, *(jGp->giveCoordinates()))) {
     *   weight = this->computeWeightFunction (gpCoords, jGpCoords);
     *   if (weight > NonlocalMaterialZeroWeight) {
     *    localIntegrationRecord ir;
     *    ir.nearGp = jGp;                   // store gp
     *    elemVolume = weight * jGp->giveElement()->computeVolumeAround (jGp);
     *    ir.weight = elemVolume;            // store gp weight
     *    iList->append(ir); // store own copy in list
     *    integrationVolume += elemVolume;
     *   }
     * } else _error ("buildNonlocalPointTable: computeGlobalCoordinates failed");
     * }
     * }
     * } // loop over elements
     * statusExt->setIntegrationScale (integrationVolume); // remember scaling factor
     */
}

void
NonlocalMaterialExtensionInterface :: rebuildNonlocalPointTable(GaussPoint *gp, IntArray *contributingElems)
{
    double weight, elemVolume, integrationVolume = 0.;

    NonlocalMaterialStatusExtensionInterface *statusExt =
        static_cast< NonlocalMaterialStatusExtensionInterface * >( gp->giveMaterialStatus()->
        giveInterface(NonlocalMaterialStatusExtensionInterfaceType) );
    std::list< localIntegrationRecord > *iList;

    Element *ielem;
    GaussPoint *jGp;
    IntegrationRule *iRule;

    if ( !statusExt ) {
        OOFEM_ERROR("NonlocalMaterialExtensionInterface::buildNonlocalPointTable : local material status encountered");
    }

    // test for bounded support - if no bounded support, the nonlocal point table is
    // big vasting of space.
    /*
     * if ( !this->hasBoundedSupport() ) {
     *  return;
     * }
     */

    iList = statusExt->giveIntegrationDomainList();
    iList->clear();

    if ( contributingElems == NULL ) {
        // no element table provided, use standard method
        this->buildNonlocalPointTable(gp);
    } else {
        FloatArray gpCoords, jGpCoords;
        int _size = contributingElems->giveSize();
        if ( gp->giveElement()->computeGlobalCoordinates( gpCoords, * ( gp->giveCoordinates() ) ) == 0 ) {
            OOFEM_ERROR("NonlocalMaterialExtensionInterface::buildNonlocalPointTable: computeGlobalCoordinates of target failed");
        }

        //If nonlocal variation is set to the distance-based approach calculates  new nonlocal radius
        // based on the distance from the nonlocal boundaries
        if ( nlvar == NLVT_DistanceBased ) {
            cl = cl0;
            cl = giveDistanceBasedInteractionRadius(gpCoords);
            suprad = evaluateSupportRadius();
        }

        // initialize iList
        for ( int _e = 1; _e <= _size; _e++ ) {
            ielem = this->giveDomain()->giveElement( contributingElems->at(_e) );
            if ( regionMap.at( ielem->giveRegionNumber() ) == 0 ) {
                iRule = ielem->giveDefaultIntegrationRulePtr();
                for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
                    jGp = iRule->getIntegrationPoint(j);
                    if ( ielem->computeGlobalCoordinates( jGpCoords, * ( jGp->giveCoordinates() ) ) ) {
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
                        OOFEM_ERROR("NonlocalMaterialExtensionInterface::buildNonlocalPointTable: computeGlobalCoordinates of target failed");
                    }
                }
            }
        } // loop over elements

        statusExt->setIntegrationScale(integrationVolume); // remember scaling factor
#ifdef __PARALLEL_MODE
 #ifdef __VERBOSE_PARALLEL
        std::list< localIntegrationRecord > :: iterator pos;
        fprintf( stderr, "%d(%d):", gp->giveElement()->giveGlobalNumber(), gp->giveNumber() );
        for ( pos = iList->begin(); pos != iList->end(); ++pos ) {
            fprintf(stderr, "%d,%d(%e)", pos->nearGp->giveElement()->giveGlobalNumber(), pos->nearGp->giveNumber(), pos->weight);
        }

        fprintf(stderr, "\n");
 #endif
#endif
    }
}


std::list< localIntegrationRecord > *
NonlocalMaterialExtensionInterface :: giveIPIntegrationList(GaussPoint *gp)
{
    NonlocalMaterialStatusExtensionInterface *statusExt =
        static_cast< NonlocalMaterialStatusExtensionInterface * >( gp->giveMaterialStatus()->
        giveInterface(NonlocalMaterialStatusExtensionInterfaceType) );

    if ( !statusExt ) {
        OOFEM_ERROR("NonlocalMaterialExtensionInterface::givIPIntegrationList : local material status encountered");
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
        OOFEM_ERROR("NonlocalMaterialExtensionInterface::givIPIntegrationList : local material status encountered");
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
    double iwf = giveIntegralOfWeightFunction(this->domain->giveNumberOfSpatialDimensions()); 

    switch ( weightFun ) {
    case WFT_Bell: // Bell shaped function (quartic spline)
        aux = ( 1. - aux * aux );
        return aux * aux / iwf;

    case WFT_Gauss: // Gauss function
        return exp(-aux * aux) / iwf;

    case WFT_Green: // Function corresponding in 1D to Green's function of Helmholtz equation (implicit gradient model)
        return exp(-aux) / iwf;

    case WFT_Green_21: // Green function reduced from 2D to 1D
      {
	if (this->domain->giveNumberOfSpatialDimensions() != 1){
	  OOFEM_ERROR("NonlocalMaterialExtensionInterface :: computeWeightFunction - this type of weight function can be used for a 1D problem only\n");
	}
	iwf = giveIntegralOfWeightFunction(2); // indeed
	double x = distance;
	double y = 0.;
	double r = sqrt(x*x + y*y);
	double sum = exp(-r / this->cl);
	double h = suprad / 100.; // 100 should later be replaced by an optional parameter
	do {
	  y += h;
	  r = sqrt(x*x + y*y);
	  sum += 2.*exp(-r / this->cl);
	}
	while (r <= suprad);
        return sum / iwf;
      }

    case WFT_Uniform: // uniform function over an interaction distance
        return 1. / iwf;

    default:
        OOFEM_WARNING2("NonlocalMaterialExtensionInterface :: computeWeightFunction - unknown type of weight function %d", weightFun);
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
    const double pi = 3.1415926;
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
        case 1: return cl * sqrt(pi);

        case 2: return cl * cl * pi;

        case 3: return cl * cl * cl * sqrt(pi * pi * pi) / 4.;

        default: return 1.;
        }

    case WFT_Green:
        switch ( spatial_dimension ) {
        case 1: return cl * 2.;

        case 2: return cl * cl * 2. * pi;

        case 3: return cl * cl * cl * 2. * pi;

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
    double iwf = giveIntegralOfWeightFunction(this->domain->giveNumberOfSpatialDimensions());
    return 1. / iwf;
}

double
NonlocalMaterialExtensionInterface :: evaluateSupportRadius()
{
    switch ( weightFun ) {
    case WFT_Bell:  return cl;

    case WFT_Gauss: return 2.5 * cl;

    case WFT_Green: return 6. * cl;

    case WFT_Uniform:  return cl;

    case WFT_UniformOverElement:  return 0.0; // to make sure that only Gauss points of the same element will be considered as neighbors

    default:        return cl;
    }
}


IRResultType
NonlocalMaterialExtensionInterface :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;              // Required by IR_GIVE_FIELD macro

    if ( ir->hasField(_IFT_NonlocalMaterialExtensionInterface_regionmap) ) {
        IR_GIVE_FIELD(ir, regionMap, _IFT_NonlocalMaterialExtensionInterface_regionmap);
        if ( regionMap.giveSize() != this->giveDomain()->giveNumberOfRegions() ) {
            OOFEM_ERROR("NonlocalMaterialExtensionInterface::instanciateFrom: regionMap size mismatch");
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
    this->weightFun = (WeightFunctionType)val;

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
    this->scaling = (ScalingType)val;

    // read the type of averaged variable
    val = AVT_EqStrain;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_NonlocalMaterialExtensionInterface_averagedquantity);
    this->averagedVar = (AveragedVarType)val;

    //Read the nonlocal variation type (default is zero)
    cl0 = cl;
    int nlvariation = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, nlvariation, _IFT_NonlocalMaterialExtensionInterface_nonlocalvariation);
    if ( nlvariation == 1 ) {
        nlvar = NLVT_DistanceBased;
        IR_GIVE_FIELD(ir, beta, _IFT_NonlocalMaterialExtensionInterface_beta);
        IR_GIVE_FIELD(ir, zeta, _IFT_NonlocalMaterialExtensionInterface_zeta);
    } else if ( nlvariation == 2 ) {
        nlvar = NLVT_StressBased;
        IR_GIVE_FIELD(ir, beta, _IFT_NonlocalMaterialExtensionInterface_beta);
    }

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
    if ( nlvar == NLVT_DistanceBased ) {
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
            weight = 1. / ( iRule->giveNumberOfIntegrationPoints() ); //asign the same weights over the whole element
        }
    }
}


double
NonlocalMaterialExtensionInterface :: giveDistanceBasedInteractionRadius(const FloatArray &gpCoords)
{
    double distance = zeta * cl0; // Initially distance from the boundary is set to the maximum value
    double temp;
    int ib, nbarrier = domain->giveNumberOfNonlocalBarriers();
    for ( ib = 1; ib <= nbarrier; ib++ ) { //Loop over all the nonlocal barriers to find minimum distance from the boundary
        temp = domain->giveNonlocalBarrier(ib)->calculateMinimumDistanceFromBoundary(gpCoords, zeta * cl0);
        if ( distance > temp ) { //Check to find minimum distance from boundary from all nonlocal boundaries
            distance = temp;
        }
    }

    //Calculate interaction radius based on the minimum distance from the nonlocal boundaries
    return ( ( 1 - beta ) / ( zeta * cl0 ) * distance + beta ) * cl0;
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
