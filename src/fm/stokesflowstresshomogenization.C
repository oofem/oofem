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
 *               Copyright (C) 1993 - 2010   Borek Patzak
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

#include "stokesflowstresshomogenization.h"
#include "prescribedgradient.h"
#include "inputrecord.h"
#include "usrdefsub.h"

namespace oofem {
StokesFlowStressHomogenization :: StokesFlowStressHomogenization(int i, EngngModel *_master) : StokesFlow(i, _master)
{
    // Homogenization related variables
    this->ht = HT_None;
    this->K_cf = NULL;
    this->K_cc = NULL;
    this->C = NULL;
}

StokesFlowStressHomogenization :: ~StokesFlowStressHomogenization()
{
    if ( this->K_cf ) {
        delete this->K_cf;
    }

    if ( this->K_cc ) {
        delete this->K_cc;
    }

    if ( this->C ) {
        delete this->C;
    }
}

bool StokesFlowStressHomogenization :: activateHomogenizationMode(HomogenizationType ht) {
    this->ht = ht;
    if ( ht == HT_StressDirichlet ) {
        this->updateC();
        return true;
    }

    return false;
}

bool StokesFlowStressHomogenization :: computeMacroStress(FloatArray &answer, const FloatArray &input, TimeStep *tStep)
{
    if ( this->ht == HT_StressDirichlet ) {
        PrescribedGradient *pt = dynamic_cast< PrescribedGradient * >( this->giveDomain(1)->giveBc(1) );
        if ( pt == NULL ) {
            OOFEM_ERROR("Wrong boundary conditions");
        }

        pt->setPrescribedGradientVoigt(input);
        this->solveYourselfAt(tStep);
        this->computeMacroStressFromDirichlet(answer, tStep);
        return true;
    } else  {
        OOFEM_ERROR("Unknown homogenization type");
    }

    return false;
}

void StokesFlowStressHomogenization :: computeMacroTangent(FloatMatrix &answer, TimeStep *tStep) {
    if ( this->ht == HT_StressDirichlet ) {
        this->computeMacroStressTangentFromDirichlet(answer, tStep);
    } else {
        OOFEM_ERROR("Unknown homogenization type");
    }
}

void StokesFlowStressHomogenization :: updateC()
// This is written in a very general way, even if this class only supports incompressible flow.
// Remove C when remeshing and it'll be recreated.
// v_prescribed = C.d = (x-xbar).d;
// C = [x 0 y]
//     [0 y x]
//     [ ... ] in 2D, voigt form [d_11, d_22, d_12]
// C = [x 0 0 y z 0]
//     [0 y 0 x 0 z]
//     [0 0 z 0 x y]
//     [ ......... ] in 3D, voigt form [d_11, d_22, d_33, d_23, d_13, d_12]
{
    Domain *domain = this->giveDomain(1);
    int nNodes = domain->giveNumberOfDofManagers();

    int dofType [ 3 ] = {
        0, 0, 0
    };
    domainType d = domain->giveDomainType();
    if ( d == _2dIncompressibleFlow || d == _3dIncompressibleFlow ) {
        dofType [ 0 ] = V_u;
        dofType [ 1 ] = V_v;
        dofType [ 2 ] = V_w;
    } else if ( d == _3dMode || d == _2dPlaneStressMode || d == _PlaneStrainMode || d == _2dTrussMode )     {
        dofType [ 0 ] = D_u;
        dofType [ 1 ] = D_v;
        dofType [ 2 ] = D_w;
    } else   {
        OOFEM_ERROR("Unsupported domain type to construct C matrix");
    }

    int nsd = domain->giveNumberOfSpatialDimensions();
    int npeq = this->giveNumberOfPrescribedDomainEquations(1, EID_MomentumBalance);
    if ( !C ) {
        C = new FloatMatrix(npeq, nsd *( nsd + 1 ) / 2);
    }

    C->zero();

    PrescribedGradient *pt = dynamic_cast< PrescribedGradient * >( domain->giveBc(1) );
    if ( !pt ) {
        OOFEM_ERROR("Incorrect boundary condition, should be prescribed gradient.");
    }

    FloatArray &cCoords = pt->giveCenterCoordinate();
    double xbar = cCoords.at(1), ybar = cCoords.at(2), zbar;
    if ( nsd == 3 ) {
        zbar = cCoords.at(3);
    }

    for ( int i = 1; i <= nNodes; i++ ) {
        Node *n = domain->giveNode(i);
        FloatArray *coords = n->giveCoordinates();
        Dof *d1 = n->giveDofWithID(dofType [ 0 ]);
        Dof *d2 = n->giveDofWithID(dofType [ 1 ]);
        int k1 = d1->__givePrescribedEquationNumber();
        int k2 = d2->__givePrescribedEquationNumber();
        if ( nsd == 2 ) {
            if ( k1 ) {
                C->at(k1, 1) = coords->at(1) - xbar;
                C->at(k1, 3) = coords->at(2) - ybar;
            }

            if ( k2 ) {
                C->at(k2, 2) = coords->at(2) - ybar;
                C->at(k2, 3) = coords->at(1) - xbar;
            }
        } else   { // nsd == 3
            OOFEM_ERROR("COMPLETELY UNTESTED!");
            Dof *d3 = n->giveDofWithID(dofType [ 2 ]);
            int k3 = d3->__givePrescribedEquationNumber();

            if ( k1 ) {
                C->at(k1, 1) = coords->at(1) - xbar;
                C->at(k1, 4) = coords->at(2) - ybar;
                C->at(k1, 5) = coords->at(3) - zbar;
            }

            if ( k2 ) {
                C->at(k2, 2) = coords->at(2) - ybar;
                C->at(k2, 4) = coords->at(1) - xbar;
                C->at(k2, 6) = coords->at(3) - zbar;
            }

            if ( k3 ) {
                C->at(k3, 3) = coords->at(3) - zbar;
                C->at(k3, 5) = coords->at(1) - xbar;
                C->at(k3, 6) = coords->at(2) - ybar;
            }
        }
    }

    C->printYourself();
}

void StokesFlowStressHomogenization :: computeMacroStressFromDirichlet(FloatArray &answer, TimeStep *tStep)
// stress = C'.R_c
{
    int npeq = this->giveNumberOfPrescribedEquations(EID_MomentumBalance_ConservationEquation);
    FloatArray R_c(npeq), R_ext(npeq);
    R_c.zero();
    R_ext.zero();
    this->assembleVectorFromElements( R_c, tStep, EID_MomentumBalance_ConservationEquation, NodalInternalForcesVector, VM_Total,
                                     EModelDefaultPrescribedEquationNumbering(), this->giveDomain(1) );
    this->assembleVectorFromElements( R_ext, tStep, EID_MomentumBalance_ConservationEquation, LoadVector, VM_Total,
                                     EModelDefaultPrescribedEquationNumbering(), this->giveDomain(1) );
    R_c.subtract(R_ext);
    answer.beTProductOf(* this->C, R_c);
}

void StokesFlowStressHomogenization :: computeMacroStressTangentFromDirichlet(FloatMatrix &tangent, TimeStep *tStep)
// a = [a_c; a_f];
// K.a = [R_c,0];
// [K_cc, K_cf; K_fc, K_ff].[a_c;a_f] = [R_c; 0];
// a_c = d.[x-x_b] = [x-x_b].d = C.d
// E = C'.(K_cc - K_cf.K_ff^(-1).K_fc).C
//   = C'.(K_cc.C - K_cf.(K_ff^(-1).(K_fc.C)))
//   = C'.(K_cc.C - K_cf.a)
//   = C'.X
{
    int neq = this->giveNumberOfEquations(EID_MomentumBalance_ConservationEquation);

    // The matrices at the last step (if the mesh is changed, equations renumbered, these should be removed.)
    if ( !this->stiffnessMatrix ) {
        OOFEM_ERROR("Trying to compute the tangent before solving the actual problem!");
    }

    if ( !this->K_cc ) {
        this->K_cc = CreateUsrDefSparseMtrx(sparseMtrxType);
        this->K_cc->buildInternalStructure( this, 1, EID_MomentumBalance_ConservationEquation,
                                           EModelDefaultPrescribedEquationNumbering(),
                                           EModelDefaultPrescribedEquationNumbering() );
    }

    if ( !this->K_cf ) {
        this->K_cf = CreateUsrDefSparseMtrx(sparseMtrxType);
        this->K_cf->buildInternalStructure( this, 1, EID_MomentumBalance_ConservationEquation,
                                           EModelDefaultPrescribedEquationNumbering(),
                                           EModelDefaultEquationNumbering() );
    }

    this->K_cc->zero();
    this->assemble( this->K_cc, tStep, EID_MomentumBalance_ConservationEquation,
                   StiffnessMatrix, EModelDefaultPrescribedEquationNumbering(),
                   this->giveDomain(1) );

    this->K_cf->zero();
    this->assemble( this->K_cf, tStep, EID_MomentumBalance_ConservationEquation,
                   StiffnessMatrix, EModelDefaultPrescribedEquationNumbering(),
                   EModelDefaultEquationNumbering(), this->giveDomain(1) );

    int dim = this->C->giveNumberOfColumns();
    FloatMatrix X, K_cfa, K_fcC, a(neq, dim);

    this->K_cf->timesT(* this->C, K_fcC);

    // This could be done in a more efficient way.
    //nMethod->solve(this->stiffnessMatrix, &K_fcC, a);
    // I propose to add this routine. Default implementation would be to solve for each column.
    // But since we can only solve for vectors, we have to copy it up.
    SparseLinearSystemNM *linMethod = CreateUsrDefSparseLinSolver(ST_Petsc, 1, this->giveDomain(1), this);
    FloatArray ai(neq), K_fcCi(neq);
    for ( int i = 0; i < dim; i++ ) {
        for ( int j = 0; j < neq; j++ ) {
            K_fcCi(j) = K_fcC(j, i);
        }

        linMethod->solve(this->stiffnessMatrix, & K_fcCi, & ai);
        for ( int j = 0; j < neq; j++ ) {
            a(j, i) = ai(j);
        }
    }

    delete(linMethod);
    this->K_cc->times(* this->C, X);
    this->K_cf->times(a, K_cfa);
    X.subtract(K_cfa);
    tangent.beTProductOf(* this->C, X);
}

IRResultType StokesFlowStressHomogenization :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";
    IRResultType result;
    int val;

    StokesFlow :: initializeFrom(ir);

    val = ( int ) ST_Petsc;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_StokesFlow_lstype, "lstype");
    this->solverType = ( LinSystSolverType ) val;

    val = ( int ) HT_None;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_StokesFlowStressHomogenization_homogenizationtype, "homogenizationtype");
    this->ht = ( HomogenizationType ) val;

    return IRRT_OK;
}
} // end namespace oofem
