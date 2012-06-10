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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "zzerrorestimator.h"
#include "domain.h"
#include "dofmanager.h"
#include "element.h"
#include "gausspnt.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "mathfem.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"
#include "timestep.h"
#include "structuralmaterial.h"
#include "structuralms.h"
#include "integrationrule.h"
#include "conTable.h"
#ifndef __MAKEDEPEND
 #include <vector>
#endif

namespace oofem {
#ifdef EXPERIMENT
FloatArray sNorms;
#endif


int
ZZErrorEstimator :: estimateError(EE_ErrorMode mode, TimeStep *tStep)
{
    int ielem, nelems = this->domain->giveNumberOfElements();
    ZZErrorEstimatorInterface *interface;
    double sNorm;
    InternalStateType type = IST_StressTensor;

    if ( mode == temporaryEM ) {
        type = IST_StressTensorTemp;                  // _error ("estimateError: temporaryEM mode not supported");
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
        _error("estimateError: unknown nodal recovery type");
    }

    // first set the domain Smoother to suitable one, keep old one to be recovered
    oldSmoother = this->domain->giveSmoother();
    this->domain->setSmoother(rm, 0); // do not delete old one

    // recover nodal values
    rm->recoverValues(type, tStep);

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
    for ( ielem = 1; ielem <= nelems; ielem++ ) {
        if ( this->skipRegion( this->domain->giveElement(ielem)->giveRegionNumber() ) ) {
            continue;
        }

        interface = ( ZZErrorEstimatorInterface * ) this->domain->giveElement(ielem)->giveInterface(ZZErrorEstimatorInterfaceType);
        if ( interface == NULL ) {
            _error("estimateError: Element has no ZZ error estimator interface defined");
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

#ifdef __PARALLEL_MODE
    // compute global ENorm and SNorm by summing up comtributions on all partitions
    double lnorms [ 2 ] = {
        this->globalENorm, this->globalSNorm
    };
    double gnorms [ 2 ] = {
        0.0, 0.0
    };
    MPI_Allreduce(lnorms, gnorms, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    this->globalENorm = gnorms [ 0 ];
    this->globalSNorm = gnorms [ 1 ];
#endif

    // recover the stored smoother
    this->domain->setSmoother(oldSmoother); //delete old one (the rm)

    // report the error estimate
    double pe = sqrt( this->globalENorm / ( this->globalENorm + this->globalSNorm ) );
    OOFEM_LOG_RELEVANT("Relative stress error estimate: %5.2f%%\n", pe * 100.0);

    this->globalENorm = sqrt(this->globalENorm);
    this->globalSNorm = sqrt(this->globalSNorm);


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
        interface = ( ZZErrorEstimatorInterface * ) elem->giveInterface(ZZErrorEstimatorInterfaceType);
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
    }
    // return sqrt (this->globalENorm/(this->globalENorm + this->globalSNorm));
    else if ( type == globalNormEEV ) {
        return this->globalSNorm;
    } else {
        return 0.0;
    }
}

RemeshingCriteria *
ZZErrorEstimator :: giveRemeshingCrit()
{
    if ( this->rc ) {
        return this->rc;
    }

    return ( this->rc = new ZZRemeshingCriteria(1, this) );
}


IRResultType
ZZErrorEstimator :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    int n;

    ErrorEstimator :: initializeFrom(ir);
    n = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, n, IFT_ZZErrorEstimator_normtype, "normtype"); // Macro
    if ( n == 1 ) {
        this->normType = EnergyNorm;
    } else {
        this->normType = L2Norm; // default
    }

    n = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, n, IFT_ZZErrorEstimator_recoverytype, "recoverytype"); // Macro
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
    int i, j, k, size, nDofMans;
    Element *elem = this->ZZErrorEstimatorI_giveElement();
    IntegrationRule *iRule = elem->giveDefaultIntegrationRulePtr();
    const FloatArray *recoveredStress;
    FloatArray sig, diff;
    FloatMatrix nodalRecoveredStreses, n;
    GaussPoint *gp;
    double dV;

    nDofMans = elem->giveNumberOfDofManagers();
    size = ( ( StructuralMaterial * ) elem->giveMaterial() )->
           giveSizeOfReducedStressStrainVector( iRule->getIntegrationPoint(0)->giveMaterialMode() );
    nodalRecoveredStreses.resize(nDofMans, size);
    // assemble nodal recovered stresses
    for ( i = 1; i <= elem->giveNumberOfNodes(); i++ ) {
        elem->giveDomain()->giveSmoother()->giveNodalVector( recoveredStress, elem->giveDofManager(i)->giveNumber(),
                                                            elem->giveRegionNumber() );
        for ( j = 1; j <= size; j++ ) {
            nodalRecoveredStreses.at(i, j) = recoveredStress->at(j);
        }
    }

    eNorm = sNorm = 0.0;
    diff.resize(size);
    diff.zero();

    // compute  the e-norm and s-norm
    if ( norm == ZZErrorEstimator :: L2Norm ) {
        for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            gp  = iRule->getIntegrationPoint(i);
            dV  = elem->computeVolumeAround(gp);
            this->ZZErrorEstimatorI_computeEstimatedStressInterpolationMtrx(n, gp, IST_StressTensor);
            for ( j = 1; j <= size; j++ ) {
                for ( k = 1; k <= nDofMans; k++ ) {
                    diff.at(j) += n.at(1, k) * nodalRecoveredStreses.at(k, j);
                }
            }

            if ( type == IST_StressTensor ) {
                sig = ( ( StructuralMaterialStatus * ) elem->giveMaterial()->giveStatus(gp) )->giveStressVector();
            } else if ( type == IST_StressTensorTemp ) {
                sig = ( ( StructuralMaterialStatus * ) elem->giveMaterial()->giveStatus(gp) )->giveTempStressVector();
            } else {
                OOFEM_ERROR("ZZErrorEstimatorI_computeElementContributions: unsuported InternalStateType");
            }

            diff.subtract(sig);

            eNorm += diff.computeSquaredNorm() * dV;
            sNorm += sig.computeSquaredNorm() * dV;
        }
    } else if ( norm == ZZErrorEstimator :: EnergyNorm ) {
        FloatArray help;
        FloatMatrix DInv;

        for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            gp  = iRule->getIntegrationPoint(i);
            dV  = elem->computeVolumeAround(gp);
            this->ZZErrorEstimatorI_computeEstimatedStressInterpolationMtrx(n, gp, IST_StressTensor);
            ( ( StructuralMaterial * ) elem->giveMaterial() )
            ->giveCharacteristicComplianceMatrix(DInv, ReducedForm, TangentStiffness,
                                                 gp, tStep);
            for ( j = 1; j <= size; j++ ) {
                for ( k = 1; k <= nDofMans; k++ ) {
                    diff.at(j) += n.at(1, k) * nodalRecoveredStreses.at(k, j);
                }
            }

            if ( type == IST_StressTensor ) {
                sig = ( ( StructuralMaterialStatus * ) elem->giveMaterial()->giveStatus(gp) )->giveStressVector();
            } else if ( type == IST_StressTensorTemp ) {
                sig = ( ( StructuralMaterialStatus * ) elem->giveMaterial()->giveStatus(gp) )->giveTempStressVector();
            } else {
                OOFEM_ERROR("ZZErrorEstimatorI_computeElementContributions: unsuported InternalStateType");
            }

            diff.subtract(sig);

            help.beProductOf(DInv, diff);
            eNorm += diff.dotProduct(help) * dV;
            help.beProductOf(DInv, sig);
            sNorm += sig.dotProduct(help)   * dV;
        }
    } else {
        OOFEM_ERROR("ZZErrorEstimatorInterface::ZZErrorEstimatorI_computeElementContributions unsupported norm type");
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
    size = max(minElemSize, size);
    if ( relative ) {
        return size / this->giveDofManDensity(num);
    } else {
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
    int i, j, nelem, nnode, jnode, elemPolyOrder, ielemNodes;
    double globValNorm = 0.0, globValErrorNorm = 0.0, elemErrLimit, eerror, iratio, currDensity, elemSize;
    Element *ielem;
    EE_ErrorType errorType = indicatorET;
    ZZRemeshingCriteriaInterface *interface;
    double pe, coeff = 2.0;

    if ( stateCounter == tStep->giveSolutionStateCounter() ) {
        return 1;
    }

    nelem = this->domain->giveNumberOfElements();
    nnode = this->domain->giveNumberOfDofManagers();

    //std::vector<char> nodalDensities(nnode);
    this->nodalDensities.resize(nnode);
    std :: vector< char >dofManInitFlag(nnode);
    for ( i = 0; i < nnode; i++ ) {
        dofManInitFlag [ i ] = 0;
    }

    // compute element error limit based on equally distribution error principle
    if ( mode == stressBased ) {
        globValNorm      = this->ee->giveValue(globalNormEEV, tStep);
        globValErrorNorm = this->ee->giveValue(globalErrorEEV, tStep);
        errorType = internalStressET;
    } else {
        _error("estimateMeshDensities: unknown mode");
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

    for ( i = 1; i <= nelem; i++ ) {
        ielem = domain->giveElement(i);

        if ( this->ee->skipRegion( ielem->giveRegionNumber() ) ) {
            continue;
        }

        interface = ( ZZRemeshingCriteriaInterface * )
                    domain->giveElement(i)->giveInterface(ZZRemeshingCriteriaInterfaceType);
        if ( !interface ) {
            _error("estimateMeshDensities: element does not support ZZRemeshingCriteriaInterface");
        }

        eerror = this->ee->giveElementError(errorType, ielem, tStep);
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

        currDensity = interface->ZZRemeshingCriteriaI_giveCharacteristicSize();
        elemPolyOrder = interface->ZZRemeshingCriteriaI_givePolynOrder();
        elemSize = currDensity / __OOFEM_POW(iratio, 1.0 / elemPolyOrder);

        ielemNodes = ielem->giveNumberOfDofManagers();
        for ( j = 1; j <= ielemNodes; j++ ) {
            jnode = ielem->giveDofManager(j)->giveNumber();
            if ( dofManInitFlag [ jnode - 1 ] ) {
                this->nodalDensities.at(jnode) = min(this->nodalDensities.at(jnode), elemSize);
            } else {
                this->nodalDensities.at(jnode) = elemSize;
                dofManInitFlag [ jnode - 1 ] = 1;
            }
        }
    }

    // init non-initialized nodes -> those in skip regions
    for ( i = 0; i < nnode; i++ ) {
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
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->requiredError, IFT_ZZRemeshingCriteria_requirederror, "requirederror"); // Macro
    IR_GIVE_FIELD(ir, this->minElemSize, IFT_ZZRemeshingCriteria_minelemsize, "minelemsize"); // Macro

    return IRRT_OK;
}


double
ZZRemeshingCriteria :: giveDofManDensity(int num)
{
    int i, isize;
    ConnectivityTable *ct = domain->giveConnectivityTable();
    const IntArray *con;
    ZZRemeshingCriteriaInterface *interface;
    double density;

    con = ct->giveDofManConnectivityArray(num);
    isize = con->giveSize();

    /*
     * // Minimum density
     *
     * for (i=1; i<=isize; i++) {
     * interface = (ZZRemeshingCriteriaInterface*)
     * domain->giveElement(con->at(i))->giveInterface (ZZRemeshingCriteriaInterfaceType);
     * if (!interface) {
     * _error ("giveDofManDensity: element does not support ZZRemeshingCriteriaInterface");
     * }
     * if (i==1) density = interface->ZZRemeshingCriteriaI_giveCharacteristicSize ();
     * else density = min (density, interface->ZZRemeshingCriteriaI_giveCharacteristicSize ());
     * }
     */

    // Average density

    density = 0.0;
    for ( i = 1; i <= isize; i++ ) {
        interface = ( ZZRemeshingCriteriaInterface * )
                    domain->giveElement( con->at(i) )->giveInterface(ZZRemeshingCriteriaInterfaceType);
        if ( !interface ) {
            _error("giveDofManDensity: element does not support ZZRemeshingCriteriaInterface");
        }

        density += interface->ZZRemeshingCriteriaI_giveCharacteristicSize();
    }

    density /= isize;

    return density;
}
} // end namespace oofem
