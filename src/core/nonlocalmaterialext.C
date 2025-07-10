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

#ifdef __MPI_PARALLEL_MODE
 #include "parallel.h"
#endif
#ifdef _OPENMP
#include <omp.h>
#endif
#include <list>

namespace oofem {
// flag forcing the inclusion of all elements with volume inside support of weight function.
// This forces inclusion of all integration points of these elements, even if weight is zero
// If not defined (default) only integration points with nonzero weight are included.
// #define NMEI_USE_ALL_ELEMENTS_IN_SUPPORT

#ifdef _OPENMP
  omp_lock_t NonlocalMaterialExtensionInterface::updateDomainBeforeNonlocAverageLock;
#endif

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

    Rf = 0.;
    exponent = 1.;
    averType = 0;

    gridSize = 0;
    initDiag = 0.;
    order = 1;
    centDiff = 2;

#ifdef _OPENMP
    omp_init_lock(&NonlocalMaterialExtensionInterface::updateDomainBeforeNonlocAverageLock);
#endif
}

void
NonlocalMaterialExtensionInterface :: updateDomainBeforeNonlocAverage(TimeStep *tStep) const
{
    Domain *d = this->domain;

    if ( d->giveNonlocalUpdateStateCounter() == tStep->giveSolutionStateCounter() ) {
        return; // already updated
    }
    #ifdef _OPENMP
    omp_set_lock(&NonlocalMaterialExtensionInterface::updateDomainBeforeNonlocAverageLock); // if not initialized yet; one thread can proceed with init; others have to wait until init completed
    if ( d->giveNonlocalUpdateStateCounter() == tStep->giveSolutionStateCounter() ) {
      omp_unset_lock(&NonlocalMaterialExtensionInterface::updateDomainBeforeNonlocAverageLock);
        return; // already updated
    }
#endif   

    OOFEM_LOG_DEBUG("Updating Before NonlocAverage\n");
    for ( auto &elem : d->giveElements() ) {
        elem->updateBeforeNonlocalAverage(tStep);
    }

    // mark last update counter to prevent multiple updates
    d->setNonlocalUpdateStateCounter( tStep->giveSolutionStateCounter() );
#ifdef _OPENMP
    omp_unset_lock(&NonlocalMaterialExtensionInterface::updateDomainBeforeNonlocAverageLock);
#endif
}

void
NonlocalMaterialExtensionInterface :: buildNonlocalPointTable(GaussPoint *gp) const
{
    double elemVolume, integrationVolume = 0.;
    double cl=this->cl, suprad;  // bp: local to be thread safe

    NonlocalMaterialStatusExtensionInterface *statusExt =
        static_cast< NonlocalMaterialStatusExtensionInterface * >( gp->giveMaterialStatus()->
                                                                   giveInterface(NonlocalMaterialStatusExtensionInterfaceType) );
    if ( !statusExt ) {
        OOFEM_ERROR("local material status encountered");
    }

    if ( !statusExt->giveIntegrationDomainList()->empty() ) {
        return;                                                  // already done
    }

    // Compute the volume around the Gauss point and store it in the nonlocal material status
    // (it will be used by modifyNonlocalWeightFunctionAround)
    elemVolume = gp->giveElement()->computeVolumeAround(gp);
    statusExt->setVolumeAround(elemVolume);

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
      suprad = evaluateSupportRadius(cl);
    } else {
      suprad = this->suprad;
    }

    // If the mesh represents a periodic cell, nonlocal interaction is considered not only for the real neighbors
    // but also for their periodic images, shifted by +px or -px in the x-direction. In the implementation,
    // instead of shifting the potential neighbors, we shift the receiver point gp. In the non-periodic case (typical),
    // px=0 and the following loop is executed only once.

    int nx = 0; // typical case
    if ( px > 0. ) {
        nx = 1;            // periodicity taken into account
    }
    for ( int ix = -nx; ix <= nx; ix++ ) { // loop over periodic images shifted in x-direction
        SpatialLocalizer :: elementContainerType elemSet;
        shiftedGpCoords = gpCoords;
        shiftedGpCoords.at(1) += ix * px;

        // ask domain spatial localizer for list of elements with IP within this zone
#ifdef NMEI_USE_ALL_ELEMENTS_IN_SUPPORT
        this->domain->giveSpatialLocalizer()->giveAllElementsWithNodesWithinBox(elemSet, shiftedGpCoords, suprad);
        // insert element containing given gp
        elemSet.insert( gp->giveElement()->giveNumber() );
#else
        this->domain->giveSpatialLocalizer()->giveAllElementsWithIpWithinBox_EvenIfEmpty(elemSet, shiftedGpCoords, suprad);
#endif
        // initialize iList
        iList->reserve(elemSet.giveSize());
        for ( auto elindx : elemSet ) {
            Element *ielem = this->domain->giveElement(elindx);
            if ( regionMap.at( ielem->giveRegionNumber() ) == 0 ) {
                for ( auto &jGp : *ielem->giveDefaultIntegrationRulePtr() ) {
                    if ( ielem->computeGlobalCoordinates( jGpCoords, jGp->giveNaturalCoordinates() ) ) {
                      double weight = this->computeWeightFunction(cl, shiftedGpCoords, jGpCoords);

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
  NonlocalMaterialExtensionInterface :: rebuildNonlocalPointTable(GaussPoint *gp, IntArray *contributingElems) const
{
    double weight, elemVolume, integrationVolume = 0.;
    double cl = this->cl;
    
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
        }

        // initialize iList
        iList->reserve(_size);
        for ( int _e = 1; _e <= _size; _e++ ) {
            Element *ielem = const_cast<NonlocalMaterialExtensionInterface*>(this)->giveDomain()->giveElement( contributingElems->at(_e) );
            if ( regionMap.at( ielem->giveRegionNumber() ) == 0 ) {
                for ( auto &jGp : *ielem->giveDefaultIntegrationRulePtr() ) {
                    if ( ielem->computeGlobalCoordinates( jGpCoords, jGp->giveNaturalCoordinates() ) ) {
                      weight = this->computeWeightFunction(cl, gpCoords, jGpCoords);

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
#ifdef __MPI_PARALLEL_MODE
 #ifdef __VERBOSE_PARALLEL
        fprintf( stderr, "%d(%d):", gp->giveElement()->giveGlobalNumber(), gp->giveNumber() );
        for ( auto &lir : *iList ) {
            fprintf(stderr, "%d,%d(%e)", lir.nearGp->giveElement()->giveGlobalNumber(), lir.nearGp->giveNumber(), lir.weight);
        }

        fprintf(stderr, "\n");
 #endif
#endif
    }
}

// This is the method used by eikonal nonlocal models to adjust the nonlocal interaction
// depending on the evolution of some internal variables such as damage
void
  NonlocalMaterialExtensionInterface :: modifyNonlocalWeightFunctionAround(GaussPoint *gp) const
{
    Element *elem = gp->giveElement();
    FloatArray coords;
    // Grid on which the eikonal equation will be solved (used by eikonal nonlocal models)
    Grid grid(2 *gridSize + 1, 2 *gridSize + 1);
    // Auxiliary matrix to store minimum distances of grid points from Gauss points
    FloatMatrix minDist2(2 *gridSize + 1, 2 *gridSize + 1);
    
    elem->computeGlobalCoordinates( coords, gp->giveNaturalCoordinates() );
    int dim = coords.giveSize();
    if ( dim == 1 ) {
        modifyNonlocalWeightFunction_1D_Around(gp);
        return;
    }

#ifdef DEBUG
    if ( dim != 2 ) {
        OOFEM_ERROR("NonlocalMaterialExtensionInterface :: modifyNonlocalWeightFunctionAround is implemented only for 1D and 2D problems, sorry.\n");
    }
#endif

    // Compute the "speed" at each grid node, depending on damage or a similar variable (note that "speed" will be deleted by grid)
    FloatMatrix &speed = grid.givePrescribedField();
    // This is a simple initialization that leads to standard Euclidean distance
    /*
     * for (i=1; i<=2*gridSize+1; i++)
     * for (j=1; j<=2*gridSize+1; j++)pos
     * speed->at(i,j) = 1.;
     */

    auto *list = this->giveIPIntegrationList(gp);
    Element *ngpElem;
    FloatArray ngpCoords(2);

    // This is the proper initialization for a truly eikonal damage model

    // Loop over neighboring Gauss points - transfer of damage to the grid "speed"
    double xgp = coords.at(1);
    double ygp = coords.at(2);
    int ipos = 0;
    double xngp, yngp, damgp;
    for ( auto lip : *list ) {
        ipos++;
        // Determine the global coordinates of the neighboring Gauss point
        ngpElem = ( lip.nearGp )->giveElement();
        ngpElem->computeGlobalCoordinates( ngpCoords, lip.nearGp->giveNaturalCoordinates() );
        // Get damage at the neighboring Gauss point
        xngp = mapToGridCoord(ngpCoords.at(1), xgp);
        yngp = mapToGridCoord(ngpCoords.at(2), ygp);
        damgp = this->giveNonlocalMetricModifierAt(lip.nearGp);
        // For the first neighbor, set damage as initial guess
        if ( ipos == 1 ) {
            for ( int i = 1; i <= 2 * gridSize + 1; i++ ) {
                for ( int j = 1; j <= 2 * gridSize + 1; j++ ) {
                    minDist2.at(i, j) = dist2FromGridNode(xngp, yngp, j, i);
                    speed.at(i, j) = damgp;
                }
            }
            // For the other neighbors, check whether distance is smaller and update damage
        } else {
            for ( int i = 1; i <= 2 * gridSize + 1; i++ ) {
                for ( int j = 1; j <= 2 * gridSize + 1; j++ ) {
                    double dist2 = dist2FromGridNode(xngp, yngp, j, i);
                    if ( dist2 < minDist2.at(i, j) ) {
                        minDist2.at(i, j) = dist2;
                        speed.at(i, j) = damgp;
                    }
                }
            }
        }
    }

    // Transform damage to speed
    for ( int i = 1; i <= 2 * gridSize + 1; i++ ) {
        for ( int j = 1; j <= 2 * gridSize + 1; j++ ) {
          speed.at(i, j) = 1. / computeDistanceModifier( cl, speed.at(i, j) );
        }
    }

    // Initialize grid for distance evaluation
    grid.unFreeze();

    // Set the details of the method that should be used by the grid
    grid.setMethod(order, initDiag, centDiff);

    // Set zero distance at the grid center
    FloatMatrix center(1, 2);
    center.at(1, 1) = center.at(1, 2) = ( double ) gridSize + 1;
    grid.setZeroValues(& center);

    // The fast marching method will be invoked implicitly, by asking for the solution

    double wsum = 0.;
    // Loop over neighboring Gauss points - evaluation of new weights
    for ( auto lip : *list ) {
        // Determine the global coordinates of a neighboring Gauss point
        ngpElem = ( lip.nearGp )->giveElement();
        ngpElem->computeGlobalCoordinates( ngpCoords, lip.nearGp->giveNaturalCoordinates() );
        // Find the nearest node of the grid (transformation from physical coordinates to grid coordinates)
        int j = mapToGridPoint( ngpCoords.at(1), coords.at(1) );
        if ( j < 1 ) {
            j = 1;
        } else if ( j > 2 * gridSize + 1 ) {
            j = 2 * gridSize + 1;
        }
        int i = mapToGridPoint( ngpCoords.at(2), coords.at(2) );
        if ( i < 1 ) {
            i = 1;
        } else if ( i > 2 * gridSize + 1 ) {
            i = 2 * gridSize + 1;
        }

        // Get solution value from the nearest grid point
        double distance = ( suprad / gridSize ) * grid.giveSolutionValueAt(i, j);
        if ( distance < 0. ) {
            printf("Warning\n");
        }
        // Compute and store the weight function for that distance
        // double volumeAround = ngpElem->computeVolumeAround(lip.nearGp); // old style
        // More efficient:
        NonlocalMaterialStatusExtensionInterface *statusExt =
            static_cast< NonlocalMaterialStatusExtensionInterface * >( lip.nearGp->giveMaterialStatus()->
                                                                       giveInterface(NonlocalMaterialStatusExtensionInterfaceType) );
        double volumeAround = statusExt->giveVolumeAround();
        double w = computeWeightFunction(cl, distance) * volumeAround;
        lip.weight = w;
        wsum += w;
    }

    // update the stored sum of nonlocal interaction weights (to be used for potential rescaling)
    NonlocalMaterialStatusExtensionInterface *statusExt =
        static_cast< NonlocalMaterialStatusExtensionInterface * >( gp->giveMaterialStatus()->
                                                                   giveInterface(NonlocalMaterialStatusExtensionInterfaceType) );
    statusExt->setIntegrationScale(wsum);
}

// Simple algorithm, limited to 1D, can be used for comparison
void
NonlocalMaterialExtensionInterface :: modifyNonlocalWeightFunction_1D_Around(GaussPoint *gp) const
{
    auto *list = this->giveIPIntegrationList(gp);
    std :: vector< localIntegrationRecord > :: iterator postarget;

    // find the current Gauss point (target) in the list of it neighbors
    for ( auto pos = list->begin(); pos != list->end(); ++pos ) {
        if ( pos->nearGp == gp ) {
            postarget = pos;
        }
    }

    Element *elem = gp->giveElement();
    FloatArray coords;
    elem->computeGlobalCoordinates( coords, gp->giveNaturalCoordinates() );
    double xtarget = coords.at(1);

    double wsum = 0., xprev, damageprev = 0.;

    // process the list from the target to the end
    double distance = 0.; // distance modified by damage
    xprev = xtarget;
    for ( auto pos = postarget; pos != list->end(); ++pos ) {
        Element *nearElem = ( pos->nearGp )->giveElement();
        nearElem->computeGlobalCoordinates( coords, pos->nearGp->giveNaturalCoordinates() );
        double x = coords.at(1);
        double damage = this->giveNonlocalMetricModifierAt(pos->nearGp);
        /*
         * nonlocStatus = static_cast< IDNLMaterialStatus * >( this->giveStatus(pos->nearGp) );
         * damage = nonlocStatus->giveTempDamage();
         * if ( damage == 0. ) {
         *  damage = nonlocStatus->giveDamage();
         * }
         */

        if ( pos != postarget ) {
            distance += computeModifiedLength(x - xprev, damage, damageprev);
            //( x - xprev ) * 0.5 * ( computeDistanceModifier(damage) + computeDistanceModifier(damageprev) );
        }

        double w = computeWeightFunction(cl, distance) * nearElem->computeVolumeAround(pos->nearGp);
        pos->weight = w;
        wsum += w;
        xprev = x;
        damageprev = damage;
    }

    // process the list from the target to the beginning
    distance = 0.;
    for ( auto pos = postarget; pos != list->begin(); --pos ) {
        Element *nearElem = ( pos->nearGp )->giveElement();
        nearElem->computeGlobalCoordinates( coords, pos->nearGp->giveNaturalCoordinates() );
        double x = coords.at(1);
        double damage = this->giveNonlocalMetricModifierAt(pos->nearGp);
        /*
         * nonlocStatus = static_cast< IDNLMaterialStatus * >( this->giveStatus(pos->nearGp) );
         * damage = nonlocStatus->giveTempDamage();
         * if ( damage == 0. ) {
         *  damage = nonlocStatus->giveDamage();
         * }
         */

        if ( pos != postarget ) {
            distance += computeModifiedLength(xprev - x, damage, damageprev);
            //distance += ( xprev - x ) * 0.5 * ( computeDistanceModifier(damage) + computeDistanceModifier(damageprev) );
            double w = computeWeightFunction(cl, distance) * nearElem->computeVolumeAround(pos->nearGp);
            pos->weight = w;
            wsum += w;
        }

        xprev = x;
        damageprev = damage;
    }

    // the beginning must be treated separately
    auto pos = list->begin();
    if ( pos != postarget ) {
        Element *nearElem = ( pos->nearGp )->giveElement();
        nearElem->computeGlobalCoordinates( coords, pos->nearGp->giveNaturalCoordinates() );
        double x = coords.at(1);
        double damage = this->giveNonlocalMetricModifierAt(pos->nearGp);
        /*
         * nonlocStatus = static_cast< IDNLMaterialStatus * >( this->giveStatus(pos->nearGp) );
         * damage = nonlocStatus->giveTempDamage();
         * if ( damage == 0. ) {
         *  damage = nonlocStatus->giveDamage();
         * }
         */

        distance += computeModifiedLength(xprev - x, damage, damageprev);
        //distance += ( xprev - x ) * 0.5 * ( computeDistanceModifier(damage) + computeDistanceModifier(damageprev) );
        double w = computeWeightFunction(cl, distance) * nearElem->computeVolumeAround(pos->nearGp);
        pos->weight = w;
        wsum += w;
    }

    // update the stored sum of nonlocal interaction weights (to be used for potential rescaling)
    NonlocalMaterialStatusExtensionInterface *statusExt =
        static_cast< NonlocalMaterialStatusExtensionInterface * >( gp->giveMaterialStatus()->
                                                                   giveInterface(NonlocalMaterialStatusExtensionInterfaceType) );
    statusExt->setIntegrationScale(wsum);
}

double
NonlocalMaterialExtensionInterface :: computeModifiedLength(double length, double dam1, double dam2) const
{
    if ( averType == 6 ) { // different (improved) integration scheme
      return length * 2. / ( 1. / computeDistanceModifier(cl, dam1) + 1. / computeDistanceModifier(cl, dam2) );
    } else { // standard integration scheme
        //printf("%g %g %g %g %g %g\n",dam1,dam2,computeDistanceModifier(dam1),computeDistanceModifier(dam2),length,length * 0.5 * ( computeDistanceModifier(dam1) + computeDistanceModifier(dam2) ));
      return length * 0.5 * ( computeDistanceModifier(cl, dam1) + computeDistanceModifier(cl, dam2) );
    }
}

double
NonlocalMaterialExtensionInterface :: computeDistanceModifier(double cl, double damage) const
{
    switch ( averType ) {
    case 2: return 1. / ( Rf / cl + ( 1. - Rf / cl ) * pow(1. - damage, exponent) );

    case 3: if ( damage == 0. ) {
            return 1.;
    } else {
            return 1. / ( 1. - ( 1. - Rf / cl ) * pow(damage, exponent) );
    }

    case 4: return 1. / pow(Rf / cl, damage);

    case 5: return ( 2. * cl ) / ( cl + Rf + ( cl - Rf ) * cos(M_PI * damage) );

    case 6: return 1. / sqrt(1. - damage);

    default: return 1.;
    }
}

std :: vector< localIntegrationRecord > *
NonlocalMaterialExtensionInterface :: giveIPIntegrationList(GaussPoint *gp) const
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
NonlocalMaterialExtensionInterface :: endIPNonlocalAverage(GaussPoint *gp) const
{
    auto statusExt = static_cast< NonlocalMaterialStatusExtensionInterface * >( gp->giveMaterialStatus()->
                                                                   giveInterface(NonlocalMaterialStatusExtensionInterfaceType) );

    if ( !statusExt ) {
        OOFEM_ERROR("local material status encountered");
    }

    if ( ( !this->hasBoundedSupport() ) || ( !permanentNonlocTableFlag ) ) {
        statusExt->clear();
    }
}

double
NonlocalMaterialExtensionInterface :: computeWeightFunction(double cl, double distance) const
{
    if ( weightFun == WFT_UniformOverElement ) { // uniform function over one element
        return 1.;
    }

    if ( distance > suprad || distance < 0. ) {
        return 0.;
    }

    double aux = distance / cl;
    double iwf = giveIntegralOfWeightFunction( cl, this->domain->giveNumberOfSpatialDimensions() );

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
        iwf = giveIntegralOfWeightFunction(cl, 2); // indeed
        double x = distance;
        double y = 0.;
        double r = sqrt(x * x + y * y);
        double sum = exp(-r / cl);
        double h = cl / 10.; // 10 could later be replaced by an optional parameter
        do {
            y += h;
            r = sqrt(x * x + y * y);
            sum += 2. * exp(-r / cl);
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
  NonlocalMaterialExtensionInterface :: computeWeightFunction(const double cl, const FloatArray &src, const FloatArray &coord) const
{
  return computeWeightFunction( cl, distance(src, coord) );
}

double
NonlocalMaterialExtensionInterface :: giveIntegralOfWeightFunction(double cl, const int spatial_dimension) const
{
    const double pi = M_PI;
    switch ( weightFun ) {
    case WFT_Bell:
        switch ( spatial_dimension ) {
        case 1: return cl * 16. / 15.;

        case 2: return cl * cl * pi / 3.;

        case 3: return cl * cl * cl * pi * 32. / 105.;

        default: return 1.;
        }

    case WFT_Gauss:
        switch ( spatial_dimension ) {
        case 1: return cl * sqrt(pi);

        case 2: return cl * cl * pi;

        case 3: return cl * cl * cl * sqrt(pi * pi * pi);

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
  double iwf = giveIntegralOfWeightFunction( cl, this->domain->giveNumberOfSpatialDimensions() );
    return 1. / iwf;
}

double
NonlocalMaterialExtensionInterface :: evaluateSupportRadius(double cl) const
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


void
NonlocalMaterialExtensionInterface :: initializeFrom(InputRecord &ir)
{
    if ( ir.hasField(_IFT_NonlocalMaterialExtensionInterface_regionmap) ) {
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
    suprad = evaluateSupportRadius(cl);

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

    // read parameters used by models with variable characteristic length (eikonal nonlocal models)
    averType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, averType, _IFT_NonlocalMaterialExtensionInterface_averagingtype);
    if ( averType == 2 ) {
        exponent = 0.5; // default value for averaging type 2
    }

    if ( averType == 3 ) {
        exponent = 1.; // default value for averaging type 3
    }

    if ( averType == 2 || averType == 3 ) {
        IR_GIVE_OPTIONAL_FIELD(ir, exponent, _IFT_NonlocalMaterialExtensionInterface_exp);
    }

    if ( averType >= 2 && averType <= 5 ) {
        IR_GIVE_OPTIONAL_FIELD(ir, Rf, _IFT_NonlocalMaterialExtensionInterface_rf);
    }

    if ( averType >= 2 && averType <= 6 ) { // eikonal models
        gridSize = 10; // default value
        IR_GIVE_OPTIONAL_FIELD(ir, gridSize, _IFT_NonlocalMaterialExtensionInterface_gridsize);
        order = 1; // default value
        IR_GIVE_OPTIONAL_FIELD(ir, order, _IFT_NonlocalMaterialExtensionInterface_order);
        initDiag = 0.; // default value
        IR_GIVE_OPTIONAL_FIELD(ir, initDiag, _IFT_NonlocalMaterialExtensionInterface_initdiag);
        centDiff = 2; // default value
        IR_GIVE_OPTIONAL_FIELD(ir, centDiff, _IFT_NonlocalMaterialExtensionInterface_centdiff);
    }
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

    input.setField(this->averType, _IFT_NonlocalMaterialExtensionInterface_averagingtype);
    if ( averType == 2 || averType == 3 ) {
        input.setField(this->exponent, _IFT_NonlocalMaterialExtensionInterface_exp);
    }
    if ( averType >= 2 && averType <= 5 ) {
        input.setField(this->Rf, _IFT_NonlocalMaterialExtensionInterface_rf);
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
NonlocalMaterialExtensionInterface :: applyBarrierConstraints(const FloatArray &gpCoords, const FloatArray &jGpCoords, double &weight) const
{
    int ib, nbarrier = domain->giveNumberOfNonlocalBarriers();
    bool shieldFlag = false;

    for ( ib = 1; ib <= nbarrier; ib++ ) {
      domain->giveNonlocalBarrier(ib)->applyConstraint(cl, gpCoords, jGpCoords, weight, shieldFlag, *this);
        if ( shieldFlag ) {
            weight = 0.0;
            return;
        }
    }
}

void
NonlocalMaterialExtensionInterface :: manipulateWeight(double &weight, GaussPoint *gp, GaussPoint *jGp) const
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
NonlocalMaterialExtensionInterface :: giveDistanceBasedInteractionRadius(const FloatArray &gpCoords) const
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
        } else {
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
    volumeAround = 0.;
}

NonlocalMaterialStatusExtensionInterface :: ~NonlocalMaterialStatusExtensionInterface()
{
    ;
}
} // end namespace oofem
