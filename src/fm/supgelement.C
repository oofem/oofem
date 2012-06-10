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

#include "supgelement.h"
#include "domain.h"
#include "load.h"
#include "gausspnt.h"
#include "intarray.h"
#include "flotarry.h"
#include "flotmtrx.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

namespace oofem {
SUPGElement :: SUPGElement(int n, Domain *aDomain) :
    FMElement(n, aDomain)
{ }


SUPGElement :: ~SUPGElement()
// Destructor.
{ }

IRResultType
SUPGElement :: initializeFrom(InputRecord *ir)
{
    //const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    //IRResultType result;                   // Required by IR_GIVE_FIELD macro

    FMElement :: initializeFrom(ir);
    this->computeGaussPoints();
    return IRRT_OK;
}



void
SUPGElement :: giveCharacteristicMatrix(FloatMatrix &answer,
                                         CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver according to mtrx
//
{
    if ( mtrx == AccelerationTerm_MB ) {
        this->computeAccelerationTerm_MB(answer, tStep);
    } else if ( mtrx == AdvectionDerivativeTerm_MB ) {
        this->computeAdvectionDerivativeTerm_MB(answer, tStep);
    } else if ( mtrx == DiffusionDerivativeTerm_MB ) {
        this->computeDiffusionDerivativeTerm_MB(answer, TangentStiffness, tStep);
    } else if ( mtrx == SecantDiffusionDerivativeTerm_MB ) {
        this->computeDiffusionDerivativeTerm_MB(answer, SecantStiffness, tStep);
    } else if ( mtrx == TangentDiffusionDerivativeTerm_MB ) {
        this->computeDiffusionDerivativeTerm_MB(answer, TangentStiffness, tStep);
    } else if ( mtrx == InitialDiffusionDerivativeTerm_MB ) {
        this->computeDiffusionDerivativeTerm_MB(answer, ElasticStiffness, tStep);
    } else if ( mtrx == PressureTerm_MB ) {
        this->computePressureTerm_MB(answer, tStep);
    } else if ( mtrx == LinearAdvectionTerm_MC ) {
        this->computeLinearAdvectionTerm_MC(answer, tStep);
    } else if ( mtrx == AdvectionDerivativeTerm_MC ) {
        this->computeAdvectionDerivativeTerm_MC(answer, tStep);
    } else if ( mtrx == DiffusionDerivativeTerm_MC ) {
        this->computeDiffusionDerivativeTerm_MC(answer, tStep);
    } else if ( mtrx == AccelerationTerm_MC ) {
        this->computeAccelerationTerm_MC(answer, tStep);
    } else if ( mtrx == PressureTerm_MC ) {
        this->computePressureTerm_MC(answer, tStep);
    } else if ( mtrx == BCLhsTerm_MB ) {
      this->computeBCLhsTerm_MB(answer, tStep);
    } else if ( mtrx == BCLhsPressureTerm_MB ) {
      this->computeBCLhsPressureTerm_MB(answer, tStep);
    } else if ( mtrx == LSICStabilizationTerm_MB ) {
      this->computeLSICStabilizationTerm_MB(answer, tStep);
    } else if ( mtrx == StiffnessMatrix) {
        // support for stokes solver
        IntArray vloc, ploc;
        FloatMatrix h;
        int size = this->computeNumberOfDofs(EID_MomentumBalance_ConservationEquation);
        this->giveLocalVelocityDofMap(vloc);
        this->giveLocalPressureDofMap(ploc);
        answer.resize(size,size); answer.zero();
        //this->computeAdvectionDerivativeTerm_MB(h, tStep);  answer.assemble(h, vloc);
        this->computeDiffusionDerivativeTerm_MB(h, TangentStiffness, tStep); answer.assemble(h, vloc);
        this->computePressureTerm_MB(h, tStep); answer.assemble(h, vloc, ploc);
        this->computeLinearAdvectionTerm_MC(h, tStep); answer.assemble(h, ploc, vloc);
        this->computeAdvectionDerivativeTerm_MC(h, tStep); answer.assemble(h, ploc, vloc);
        this->computeDiffusionDerivativeTerm_MC(h, tStep); answer.assemble(h, ploc, vloc);
        this->computePressureTerm_MC(h, tStep); answer.assemble(h, ploc);
        this->computeBCLhsTerm_MB(h, tStep); answer.assemble(h, vloc);
        this->computeBCLhsPressureTerm_MB(h, tStep); answer.assemble(h, vloc, ploc);
        //this->computeLSICStabilizationTerm_MB(h, tStep); answer.assemble(h, vloc);
    } else {
        _error("giveCharacteristicMatrix: Unknown Type of characteristic mtrx.");
    }
}


void
SUPGElement :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
                                         TimeStep *tStep)
//
// returns characteristics vector of receiver according to requested type
//
{
    if ( mtrx == AdvectionTerm_MB ) {
        this->computeAdvectionTerm_MB(answer, tStep);
    } else if ( mtrx == DiffusionTerm_MB ) {
        this->computeDiffusionTerm_MB(answer, tStep);
    } else if ( mtrx == AdvectionTerm_MC ) {
        this->computeAdvectionTerm_MC(answer, tStep);
    } else if ( mtrx == BCRhsTerm_MB ) {
        this->computeBCRhsTerm_MB(answer, tStep);
    } else if ( mtrx == BCRhsTerm_MC ) {
        this->computeBCRhsTerm_MC(answer, tStep);
    } else if ( mtrx == DiffusionTerm_MC ) {
        this->computeDiffusionTerm_MC(answer, tStep);
    } else if ( mtrx == ExternalForcesVector) {
        // stokes flow
        IntArray vloc, ploc;
        FloatArray h;
        int size = this->computeNumberOfDofs(EID_MomentumBalance_ConservationEquation);
        this->giveLocalVelocityDofMap(vloc);
        this->giveLocalPressureDofMap(ploc);
        answer.resize(size); answer.zero();
        this->computeBCRhsTerm_MB(h, tStep); answer.assemble(h, vloc);
        this->computeBCRhsTerm_MC(h, tStep); answer.assemble(h, ploc);
    } else if ( mtrx == InternalForcesVector) {
        // stokes flow
        IntArray vloc, ploc;
        FloatArray h;
        int size = this->computeNumberOfDofs(EID_MomentumBalance_ConservationEquation);
        this->giveLocalVelocityDofMap(vloc);
        this->giveLocalPressureDofMap(ploc);
        answer.resize(size); answer.zero();
        //this->computeAdvectionTerm_MB(h, tStep); answer.assemble(h, vloc);
        this->computeAdvectionTerm_MC(h, tStep); answer.assemble(h, ploc);
        this->computeDiffusionTerm_MB(h, tStep); answer.assemble(h, vloc);
        this->computeDiffusionTerm_MC(h, tStep); answer.assemble(h, ploc);

        FloatMatrix m1;
        FloatArray v,p;
        // add lsic stabilization term
        //this->giveCharacteristicMatrix(m1, LSICStabilizationTerm_MB, tStep);
        //m1.times( lscale / ( dscale * uscale * uscale ) );
        this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, v);
        //h.beProductOf(m1, v);
        //answer.assemble(h, vloc);
        this->giveCharacteristicMatrix(m1, LinearAdvectionTerm_MC, tStep);
        //m1.times( 1. / ( dscale * uscale ) );
        h.beProductOf(m1, v);
        answer.assemble(h,ploc);

        // add pressure term
        this->giveCharacteristicMatrix(m1, PressureTerm_MB, tStep);
        this->computeVectorOf(EID_ConservationEquation, VM_Total, tStep, v);
        h.beProductOf(m1, v);
        answer.assemble(h, vloc);

        // pressure term
        this->giveCharacteristicMatrix(m1, PressureTerm_MC, tStep);
        this->computeVectorOf(EID_ConservationEquation, VM_Total, tStep, v);
        h.beProductOf(m1, v);
        answer.assemble(h, ploc);

    } else {
        _error("giveCharacteristicVector: Unknown Type of characteristic mtrx.");
    }
}





void
SUPGElement :: computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray eps;

    // compute deviatoric strain
    this->computeDeviatoricStrain(eps, gp, tStep);
    // call material to compute stress
    ( ( FluidDynamicMaterial * ) this->giveMaterial() )->computeDeviatoricStressVector(answer, gp, eps, tStep);
}



void
SUPGElement :: computeBCLhsTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    bcType boundarytype;
    int i, n, side;
    int nLoads = 0;
    int undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    Load *load;
    //bcType loadtype;
    FloatMatrix helpMatrix;
    // loop over boundary load array
    helpMatrix.resize(undofs, undofs);
    helpMatrix.zero();

    answer.resize(undofs, undofs);
    answer.zero();

    nLoads    = this->giveBoundaryLoadArray()->giveSize() / 2;

    if ( nLoads ) {
        for ( i = 1; i <= nLoads; i++ ) {
            n     = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
            side    = boundaryLoadArray.at(i * 2);
            load  = dynamic_cast< Load * >( domain->giveLoad(n) );
            boundarytype = load->giveType();
            if ( boundarytype == SlipWithFriction ) {
                this->computeSlipWithFrictionBCTerm_MB(helpMatrix, ( Load * ) load, side, atTime);
            } else if ( boundarytype == PenetrationWithResistance ) {
                this->computePenetrationWithResistanceBCTerm_MB(helpMatrix, ( Load * ) load, side, atTime);
            } else {
                helpMatrix.resize(undofs, undofs);
                helpMatrix.zero();
                // _error("computeForceLoadVector : unsupported load type class");
            }

            answer.add(helpMatrix);
        }
    }
}


void
SUPGElement :: computeBCLhsPressureTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    bcType boundarytype;
    int i, n, side;
    int nLoads = 0;
    int undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);
    Load *load;
    //bcType loadtype;
    FloatMatrix helpMatrix;
    // loop over boundary load array
    helpMatrix.resize(undofs, pndofs);
    helpMatrix.zero();

    answer.resize(undofs, pndofs);
    answer.zero();

    nLoads    = this->giveBoundaryLoadArray()->giveSize() / 2;

    if ( nLoads ) {
        for ( i = 1; i <= nLoads; i++ ) {
            n     = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
            side    = boundaryLoadArray.at(i * 2);
            load  = dynamic_cast< Load * >( domain->giveLoad(n) );
            boundarytype = load->giveType();
            if ( boundarytype == OutFlowBC ) {
                this->computeOutFlowBCTerm_MB(helpMatrix, side, atTime);
            } else {
                helpMatrix.resize(undofs, pndofs);
                helpMatrix.zero();
                //_warning("computeForceLoadVector : unsupported load type class");
            }

            answer.add(helpMatrix);
        }
    }
}




double
SUPGElement :: giveCharacteristicValue(CharType mtrx, TimeStep *tStep)
{
    if ( mtrx == CriticalTimeStep ) {
        return this->computeCriticalTimeStep(tStep);
    } else {
        _error("giveCharacteristicValue: Unknown Type of characteristic mtrx.");
    }

    return 0.0;
}



int
SUPGElement :: checkConsistency()
//
// check internal consistency
// mainly tests, whether material and crossSection data
// are safe for conversion to "Structural" versions
//
{
    int result = 1;
    /*
     * if (!this->giveMaterial()->testMaterialExtension(Material_TransportCapability)) {
     * _warning("checkConsistency : material without support for transport problems");
     * result =0;
     * }
     */
    return result;
}

void
SUPGElement :: updateInternalState(TimeStep *stepN)
{
    int i, j;
    IntegrationRule *iRule;
    FloatArray stress;

    // force updating strains & stresses
    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            computeDeviatoricStress(stress, iRule->getIntegrationPoint(j), stepN);
        }
    }
}

void
SUPGElement :: printOutputAt(FILE *file, TimeStep *stepN)
// Performs end-of-step operations.
{
    int i;

#ifdef __PARALLEL_MODE
    fprintf( file, "element %d [%8d] :\n", this->giveNumber(), this->giveGlobalNumber() );
#else
    fprintf(file, "element %d :\n", number);
#endif

    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        integrationRulesArray [ i ]->printOutputAt(file, stepN);
    }
}


#ifdef __OOFEG
int
SUPGElement :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                       int node, TimeStep *atTime)
{
    int indx = 1;
    Node *n = this->giveNode(node);

    if ( type == IST_Velocity ) {
        answer.resize( this->giveSpatialDimension() );
        int dofindx;
        if ( ( dofindx = n->findDofWithDofId(V_u) ) ) {
            answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(EID_MomentumBalance, VM_Total, atTime);
        }

        if ( ( dofindx = n->findDofWithDofId(V_v) ) ) {
            answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(EID_MomentumBalance, VM_Total, atTime);
        }

        if ( ( dofindx = n->findDofWithDofId(V_w) ) ) {
            answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(EID_MomentumBalance, VM_Total, atTime);
        }

        return 1;
    } else if ( type == IST_Pressure ) {
        int dofindx;
        if ( ( dofindx = n->findDofWithDofId(P_f) ) ) {
            answer.resize(1);
            answer.at(1) = n->giveDof(dofindx)->giveUnknown(EID_ConservationEquation, VM_Total, atTime);
            return 1;
        } else {
            return 0;
        }
    } else {
        return Element :: giveInternalStateAtNode(answer, type, mode, node, atTime);
    }
}

#endif

int
SUPGElement::giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_DeviatoricStrain ) {
        this->computeDeviatoricStrain (answer, gp, tStep);
        return 1;
    } else if ( type == IST_DeviatoricStress ) {
        this->computeDeviatoricStress (answer, gp, tStep);
        return 1;
    } else {
        return FMElement::giveIPValue(answer, gp, type, tStep);
    }
}

int
SUPGElement::giveIPValueSize(InternalStateType type, GaussPoint *gp)
{
    MaterialMode mmode = gp->giveMaterialMode();
    if ( ( type == IST_DeviatoricStrain ) || (type == IST_DeviatoricStress) ) {
        if ( mmode == _2dFlow ) {
            return 3;
        } else if (mmode == _2dAxiFlow ) {
            return 4;
        } else if (mmode == _3dFlow ) {
            return 6;
        } else {
            OOFEM_ERROR ("SUPGElement::giveIPValueSize: material mode not supported");
            return 0;
        }
    } else {
        return FMElement::giveIPValueSize(type, gp);
    }
}

InternalStateValueType
SUPGElement::giveIPValueType(InternalStateType type)
{
    if ( ( type == IST_DeviatoricStrain ) || (type == IST_DeviatoricStress ) ) {
        return ISVT_TENSOR_S3;
    } else {
        return FMElement::giveIPValueType(type);
    }
}

int
SUPGElement::giveIntVarCompFullIndx(IntArray &answer, InternalStateType type)
{
    MaterialMode mmode = this->giveMaterialMode();
    if ( ( type == IST_DeviatoricStrain ) || ( type == IST_DeviatoricStress ) ) {
        if ( mmode == _2dFlow ) {
            answer.resize(6);
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(6) = 3;
            return 1;
        } else if (mmode == _2dAxiFlow ) {
            answer.resize(6);
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(3) = 3;
            answer.at(6) = 4;
            return 1;
        } else if (mmode == _3dFlow ) {
            answer.resize(6);
            for (int i=1; i<=6; i++) answer.at(i) = i;
            return 1;
        } else {
            OOFEM_ERROR ("FluidDynamicMaterial :: giveIntVarCompFullIndx: material mode not supported");
            return 0;
        }
    } else {
        return FMElement::giveIntVarCompFullIndx(answer, type);
    }
}


#if 0
void
SUPGElement::computeVectorOfPrescribed (EquationID ut, ValueModeType type, TimeStep* stepN, FloatArray& answer)
{
    double scale;
    Element::computeVectorOfPrescribed (ut, type, stepN, answer);
    if (domain->giveEngngModel()->giveEquationScalingFlag()) {
    if (ut == EID_MomentumBalance) {
        scale = domain->giveEngngModel()->giveVariableScale(VST_Velocity);
    } else if (ut == EID_ConservationEquation) {
        scale = domain->giveEngngModel()->giveVariableScale(VST_Pressure);
    } else scale = 1.0;
        answer.times (1.0/scale);
    }
}
#endif
} // end namespace oofem
