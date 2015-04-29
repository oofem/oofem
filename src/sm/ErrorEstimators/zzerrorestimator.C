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

#include "../sm/ErrorEstimators/zzerrorestimator.h"
#include "../sm/Elements/structuralelement.h"
#include "../sm/Materials/structuralmaterial.h"
#include "domain.h"
#include "dofmanager.h"
#include "element.h"
#include "gausspoint.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "mathfem.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "timestep.h"
#include "integrationrule.h"
#include "feinterpol.h"
#include "connectivitytable.h"
#include "errorestimatortype.h"
#include "classfactory.h"
#include "engngm.h"
#include "parallelcontext.h"

#include <vector>

namespace oofem {
REGISTER_ErrorEstimator(ZZErrorEstimator, EET_ZZEE);

#ifdef EXPERIMENT
FloatArray sNorms;
#endif


int
ZZErrorEstimator :: estimateError(EE_ErrorMode mode, TimeStep *tStep)
{
    int nelems = this->domain->giveNumberOfElements();
    ZZErrorEstimatorInterface *interface;
    double sNorm;
    InternalStateType type = IStype;

    if ( mode == temporaryEM ) {
        type = IST_StressTensorTemp;                  // OOFEM_ERROR("temporaryEM mode not supported");
    }

    if ( this->stateCounter == tStep->giveSolutionStateCounter() ) {
        return 1;
    }

    NodalRecoveryModel *oldSmoother, *rm = NULL;
    if ( this->nodalRecoveryType == ZZRecovery ) {
        rm = new ZZNodalRecoveryModel(this->domain);
    } else if ( this->nodalRecoveryType == SPRRecovery ) {
        rm = new SPRNodalRecoveryModel(this->domain);
    } else {
        OOFEM_ERROR("unknown nodal recovery type");
    }

    // first set the domain Smoother to suitable one, keep old one to be recovered
    oldSmoother = this->domain->giveSmoother();
    this->domain->setSmoother(rm, 0); // do not delete old one

    // create a new set containing all elements
    Set elemSet(0, this->domain);
    elemSet.addAllElements();
    // recover nodal values on entire elemSet (domain)
    rm->recoverValues(elemSet, type, tStep);

#ifdef ZZErrorEstimator_ElementResultCashed
    this->eNorms.resize(nelems);
 #ifdef EXPERIMENT
    sNorms.resize(nelems);
 #endif
#else
    double eNorm;
#endif

    this->globalENorm = this->globalSNorm = 0.0;
    // loop over domain's elements
    for ( int ielem = 1; ielem <= nelems; ielem++ ) {
        if ( this->skipRegion( this->domain->giveElement(ielem)->giveRegionNumber() ) ) {
            continue;
        }

        interface = static_cast< ZZErrorEstimatorInterface * >( this->domain->giveElement(ielem)->giveInterface(ZZErrorEstimatorInterfaceType) );
        if ( interface == NULL ) {
            OOFEM_ERROR("Element has no ZZ error estimator interface defined");
        }

#ifdef ZZErrorEstimator_ElementResultCashed
        interface->ZZErrorEstimatorI_computeElementContributions(eNorms.at(ielem), sNorm, this->normType, type, tStep);
        this->globalENorm += eNorms.at(ielem) * eNorms.at(ielem);
 #ifdef EXPERIMENT
        sNorms.at(ielem) = sNorm;
 #endif
#else
        interface->ZZErrorEstimatorI_computeElementContributions(eNorm, sNorm, this->normType, type, tStep);
        this->globalENorm += eNorm * eNorm;
#endif
        this->globalSNorm += sNorm * sNorm;
    }

    FloatArray gnorms;
    ParallelContext *parallel_context = this->domain->giveEngngModel()->giveParallelContext(this->domain->giveNumber());
    parallel_context->accumulate({this->globalENorm, this->globalSNorm}, gnorms);
    this->globalENorm = gnorms [ 0 ];
    this->globalSNorm = gnorms [ 1 ];
    

    // recover the stored smoother
    this->domain->setSmoother(oldSmoother); //delete old one (the rm)

    // report the error estimate
    double pe = sqrt( this->globalENorm / ( this->globalENorm + this->globalSNorm ) );
    OOFEM_LOG_RELEVANT("Relative stress error estimate: %5.2f%%\n", pe * 100.0);

    this->globalENorm = sqrt(this->globalENorm);
    this->globalSNorm = sqrt(this->globalSNorm);
    this->globalErrorEstimate = pe;


    this->stateCounter = tStep->giveSolutionStateCounter();
    return 1;
}

double
ZZErrorEstimator :: giveElementError(EE_ErrorType type, Element *elem, TimeStep *tStep)
{
    if ( this->skipRegion( elem->giveRegionNumber() ) ) {
        return 0.0;
    }

    this->estimateError(equilibratedEM, tStep);

#ifdef ZZErrorEstimator_ElementResultCashed
    if ( type == internalStressET ) {
        return this->eNorms.at( elem->giveNumber() );
    }

#else
    ZZErrorEstimatorInterface *interface;
    if ( type == internalStressET ) {
        double e, s;
        interface = static_cast< ZZErrorEstimatorInterface * >( elem->giveInterface(ZZErrorEstimatorInterfaceType) );
        if ( interface ) {
            interface->ZZErrorEstimatorI_computeElementContributions(e, s, this->normType, tStep);
            return e;
        } else {
            return 0.0;
        }
    }

#endif
    return 0.0;
}

double
ZZErrorEstimator :: giveValue(EE_ValueType type, TimeStep *tStep)
{
    this->estimateError(equilibratedEM, tStep);
    if ( type == globalErrorEEV ) {
        return this->globalENorm;
    } else if ( type == globalNormEEV ) {
        return this->globalSNorm;
    } else if ( type == relativeErrorEstimateEEV ) {
        return this->globalErrorEstimate;
    } else {
        return 0.0;
    }
}

RemeshingCriteria *
ZZErrorEstimator :: giveRemeshingCrit()
{
    if ( !this->rc ) {
        this->rc.reset( new ZZRemeshingCriteria(1, this) );
    }

    return this->rc.get();
}


IRResultType
ZZErrorEstimator :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    int n;

    ErrorEstimator :: initializeFrom(ir);
    n = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, n, _IFT_ZZErrorEstimator_normtype);
    if ( n == 1 ) {
        this->normType = EnergyNorm;
    } else {
        this->normType = L2Norm; // default
    }

    n = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, n, _IFT_ZZErrorEstimator_recoverytype);
    if ( n == 1 ) {
        this->nodalRecoveryType = SPRRecovery;
    } else {
        this->nodalRecoveryType = ZZRecovery; // default
    }

    return this->giveRemeshingCrit()->initializeFrom(ir);
}

void
ZZErrorEstimatorInterface :: ZZErrorEstimatorI_computeElementContributions(double &eNorm, double &sNorm,
                                                                           ZZErrorEstimator :: NormType norm,
                                                                           InternalStateType type,
                                                                           TimeStep *tStep)
{
    int nDofMans;
    FEInterpolation *interpol = element->giveInterpolation();
    const FloatArray *recoveredStress;
    FloatArray sig, lsig, diff, ldiff, n;
    FloatMatrix nodalRecoveredStreses;

    nDofMans = element->giveNumberOfDofManagers();
    // assemble nodal recovered stresses
    for ( int i = 1; i <= element->giveNumberOfNodes(); i++ ) {
        element->giveDomain()->giveSmoother()->giveNodalVector( recoveredStress,
                                                            element->giveDofManager(i)->giveNumber() );
        if ( i == 1 ) {
            nodalRecoveredStreses.resize( nDofMans, recoveredStress->giveSize() );
        }
        for ( int j = 1; j <= recoveredStress->giveSize(); j++ ) {
            nodalRecoveredStreses.at(i, j) = recoveredStress->at(j);
        }
    }
    /* Note: The recovered stresses should be in global coordinate system. This is important for shells, for example, to make
     * sure that forces and moments in the same directions are averaged. For elements where local and global coordina systems
     * are the same this does not matter.
     */

    eNorm = sNorm = 0.0;

    // compute  the e-norm and s-norm
    if ( norm == ZZErrorEstimator :: L2Norm ) {
        for ( GaussPoint *gp: *this->ZZErrorEstimatorI_giveIntegrationRule() ) {
            double dV = element->computeVolumeAround(gp);
            interpol->evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(element) );

            diff.beTProductOf(nodalRecoveredStreses, n);

            element->giveIPValue(sig, gp, type, tStep);
            /* the internal stress difference is in global coordinate system */
            diff.subtract(sig);

            eNorm += diff.computeSquaredNorm() * dV;
            sNorm += sig.computeSquaredNorm() * dV;
        }
    } else if ( norm == ZZErrorEstimator :: EnergyNorm ) {
        FloatArray help, ldiff_reduced, lsig_reduced;
        FloatMatrix D, DInv;
        StructuralElement *selem = static_cast< StructuralElement * >(element);

        for ( GaussPoint *gp: *this->ZZErrorEstimatorI_giveIntegrationRule() ) {
            double dV = element->computeVolumeAround(gp);
            interpol->evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(element) );
            selem->computeConstitutiveMatrixAt(D, TangentStiffness, gp, tStep);
            DInv.beInverseOf(D);

            diff.beTProductOf(nodalRecoveredStreses, n);

            element->giveIPValue(sig, gp, type, tStep); // returns full value now
            diff.subtract(sig);
            /* the internal stress difference is in global coordinate system */
            /* needs to be transformed into local system to compute associated energy */
            this->ZZErrorEstimatorI_computeLocalStress(ldiff, diff);
            StructuralMaterial :: giveReducedSymVectorForm( ldiff_reduced, ldiff, gp->giveMaterialMode() );

            help.beProductOf(DInv, ldiff_reduced);
            eNorm += ldiff_reduced.dotProduct(help) * dV;
            this->ZZErrorEstimatorI_computeLocalStress(lsig, sig);
            StructuralMaterial :: giveReducedSymVectorForm( lsig_reduced, lsig, gp->giveMaterialMode() );
            help.beProductOf(DInv, lsig_reduced);
            sNorm += lsig_reduced.dotProduct(help) * dV;
        }
    } else {
        OOFEM_ERROR("unsupported norm type");
    }

    eNorm = sqrt(eNorm);
    sNorm = sqrt(sNorm);
}


ZZRemeshingCriteria :: ZZRemeshingCriteria(int n, ErrorEstimator *e) : RemeshingCriteria(n, e)
{
    this->mode = stressBased;
    this->stateCounter = 0;
}

double
ZZRemeshingCriteria :: giveRequiredDofManDensity(int num, TimeStep *tStep, int relative)
{
    double size;

    this->estimateMeshDensities(tStep);
    size = this->nodalDensities.at(num);
    if ( size >= 0 ) {
        size = max(minElemSize, size);
        if ( relative ) {
            return size / this->giveDofManDensity(num);
        } else {
            return size;
        }
    } else {
        // size negative -> marks undetermined size
        return size;
    }
}


RemeshingStrategy
ZZRemeshingCriteria :: giveRemeshingStrategy(TimeStep *tStep)
{
    this->estimateMeshDensities(tStep);
    return this->remeshingStrategy;
}

int
ZZRemeshingCriteria :: estimateMeshDensities(TimeStep *tStep)
{
    int nelem, nnode, elemPolyOrder, ielemNodes;
    double globValNorm = 0.0, globValErrorNorm = 0.0, elemErrLimit, eerror, iratio, currDensity, elemSize;
    EE_ErrorType errorType = indicatorET;
    double pe, coeff = 2.0;

    if ( stateCounter == tStep->giveSolutionStateCounter() ) {
        return 1;
    }

    nelem = this->domain->giveNumberOfElements();
    nnode = this->domain->giveNumberOfDofManagers();

    //std::vector<char> nodalDensities(nnode);
    this->nodalDensities.resize(nnode);
    std :: vector< char > dofManInitFlag(nnode);
    for ( int i = 0; i < nnode; i++ ) {
        dofManInitFlag [ i ] = 0;
    }

    // compute element error limit based on equally distribution error principle
    if ( mode == stressBased ) {
        globValNorm      = this->ee->giveValue(globalNormEEV, tStep);
        globValErrorNorm = this->ee->giveValue(globalErrorEEV, tStep);
        errorType = internalStressET;
    } else {
        OOFEM_ERROR("unknown mode");
    }


    // test if solution is allowable
    pe = sqrt( globValErrorNorm * globValErrorNorm / ( globValErrorNorm * globValErrorNorm + globValNorm * globValNorm ) );
    if ( pe <= this->requiredError ) {
        this->remeshingStrategy = NoRemeshing_RS;
        // added by bp
        /*
         * for (i=1; i<=nnode; i++) nodalDensities.at(i) = this->giveDofManDensity (i);
         * stateCounter = tStep->giveSolutionStateCounter();
         * return 1;  // remove to force computing densities
         */
    } else {
        this->remeshingStrategy = RemeshingFromPreviousState_RS;
    }

    elemErrLimit = sqrt( ( globValNorm * globValNorm + globValErrorNorm * globValErrorNorm ) / nelem ) *
    this->requiredError * coeff;

    for ( auto &elem : domain->giveElements() ) {

        if ( this->ee->skipRegion( elem->giveRegionNumber() ) ) {
            continue;
        }

        eerror = this->ee->giveElementError(errorType, elem.get(), tStep);
        iratio = eerror / elemErrLimit;
        if ( fabs(iratio) < 1.e-3 ) {
            continue;
        }

        if ( iratio > 1.0 ) {
            this->remeshingStrategy = RemeshingFromPreviousState_RS;
        }

        if ( iratio < 0.5 ) {
            iratio = 0.5;
        }

        //  if (iratio > 5.0)iratio = 5.0;

        currDensity = elem->computeMeanSize();
        elemPolyOrder = elem->giveInterpolation()->giveInterpolationOrder();
        elemSize = currDensity / pow(iratio, 1.0 / elemPolyOrder);

        ielemNodes = elem->giveNumberOfDofManagers();
        for ( int j = 1; j <= ielemNodes; j++ ) {
            int jnode = elem->giveDofManager(j)->giveNumber();
            if ( dofManInitFlag [ jnode - 1 ] ) {
                this->nodalDensities.at(jnode) = min(this->nodalDensities.at(jnode), elemSize);
            } else {
                this->nodalDensities.at(jnode) = elemSize;
                dofManInitFlag [ jnode - 1 ] = 1;
            }
        }
    }

    // init non-initialized nodes -> those in skip regions
    for ( int i = 0; i < nnode; i++ ) {
        if ( dofManInitFlag [ i ] == 0 ) {
            this->nodalDensities.at(i + 1) = this->giveDofManDensity(i + 1);
        }
    }

    // remember time stamp
    stateCounter = tStep->giveSolutionStateCounter();
    return 1;
}

IRResultType
ZZRemeshingCriteria :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->requiredError, _IFT_ZZRemeshingCriteria_requirederror);
    IR_GIVE_FIELD(ir, this->minElemSize, _IFT_ZZRemeshingCriteria_minelemsize);

    return IRRT_OK;
}


double
ZZRemeshingCriteria :: giveDofManDensity(int num)
{
    int isize;
    bool init = false;
    ConnectivityTable *ct = domain->giveConnectivityTable();
    const IntArray *con;
    double density;

    con = ct->giveDofManConnectivityArray(num);
    isize = con->giveSize();

#if 0
    // Minimum density
    for ( int i = 1; i <= isize; i++ ) {
        Element *ielem = domain->giveElement( con->at(i) );
        if (i==1)
            density = ielem->computeMeanSize();
        else
            density = min(density, ielem->computeMeanSize());
    }
#endif

    // Average density

    density = 0.0;
    for ( int i = 1; i <= isize; i++ ) {
        Element *ielem = domain->giveElement( con->at(i) );
        init = true;
        density += ielem->computeMeanSize();
    }
    if ( init ) {
        density /= isize;
    } else {
        // the nodal mesh density could not be determined
        density = -1;
    }

    return density;
}
} // end namespace oofem
