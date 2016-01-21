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

#include "PrescribedGenStrainShell7shell7base.h"
#include "dofiditem.h"
#include "dofmanager.h"
#include "dof.h"
#include "valuemodetype.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "function.h"
#include "engngm.h"
#include "set.h"
#include "node.h"
#include "element.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "feinterpol.h"

#include "sparsemtrx.h"
#include "sparselinsystemnm.h"

namespace oofem {
REGISTER_BoundaryCondition(PrescribedGenStrainShell7);

double PrescribedGenStrainShell7 :: give(Dof *dof, ValueModeType mode, double time)
{
    DofIDItem id = dof->giveDofID();
    FloatArray *coords = dof->giveDofManager()->giveCoordinates();

    if ( coords->giveSize() != this->centerCoord.giveSize() ) {
        OOFEM_ERROR("PrescribedGenStrainShell7 :: give - Size of coordinate system different from center coordinate in b.c.");
    }

    double factor = 0;
    if ( mode == VM_Total ) {
        factor = this->giveTimeFunction()->evaluateAtTime(time);
    } else if ( mode == VM_Velocity ) {
        factor = this->giveTimeFunction()->evaluateVelocityAtTime(time);
    } else if ( mode == VM_Acceleration ) {
        factor = this->giveTimeFunction()->evaluateAccelerationAtTime(time);
    } else {
        OOFEM_ERROR("Should not be called for value mode type then total, velocity, or acceleration.");
    }

    // Reminder: u_i = F_ij . (x_j - xb_j) = d_ij . dx_j
    FloatArray dx;
    dx.beDifferenceOf(* coords, this->centerCoord);

    FloatArray u;
    u.beProductOf(gradient, dx);
    u.times( factor );

    switch ( id ) {
    case D_u:
    case V_u:
        return u.at(1);

    case D_v:
    case V_v:
        return u.at(2);

    case D_w:
    case V_w:
        return u.at(3);

    default:
        return 0.0;
    }
}

void PrescribedGenStrainShell7 :: setPrescribedGenStrainShell7Voigt(const FloatArray &t)
{
    int n = t.giveSize();
    if ( n == 3 ) { // Then 2D
        this->gradient.resize(2, 2);
        this->gradient.at(1, 1) = t.at(1);
        this->gradient.at(2, 2) = t.at(2);
        this->gradient.at(1, 2) = this->gradient.at(2, 1) = t.at(3);
    } else if ( n == 6 ) { // Then 3D
        this->gradient.resize(3, 3);
        this->gradient.at(1, 1) = t.at(1);
        this->gradient.at(2, 2) = t.at(2);
        this->gradient.at(3, 3) = t.at(3);
        // In voigt form, assuming the use of gamma_12 instead of eps_12
        this->gradient.at(1, 2) = this->gradient.at(2, 1) = t.at(6) * 0.5;
        this->gradient.at(1, 3) = this->gradient.at(3, 1) = t.at(5) * 0.5;
        this->gradient.at(2, 3) = this->gradient.at(3, 2) = t.at(4) * 0.5;
    } else {
        OOFEM_ERROR("setPrescribedTensorVoigt: Tensor is in strange voigt format. Should be 3 or 6. Use setPrescribedTensor directly if needed.");
    }
}

void PrescribedGenStrainShell7 :: updateCoefficientMatrix(FloatMatrix &C)
// This is written in a very general way, supporting both fm and sm problems.
// v_prescribed = C.d = (x-xbar).d;
// C = [x 0 y]
//     [0 y x]
//     [ ... ] in 2D, voigt form [d_11, d_22, d_12]
// C = [x 0 0 y z 0]
//     [0 y 0 x 0 z]
//     [0 0 z 0 x y]
//     [ ......... ] in 3D, voigt form [d_11, d_22, d_33, d_23, d_13, d_12]
{
    Domain *domain = this->giveDomain();
    int nNodes = domain->giveNumberOfDofManagers();

    int nsd = domain->giveNumberOfSpatialDimensions();
    int npeq = domain->giveEngngModel()->giveNumberOfDomainEquations( domain->giveNumber(), EModelDefaultPrescribedEquationNumbering() );
    C.resize(npeq, nsd * ( nsd + 1 ) / 2);
    C.zero();

    FloatArray &cCoords = this->giveCenterCoordinate();
    double xbar = cCoords.at(1), ybar = cCoords.at(2), zbar = 0.0;
    if ( nsd == 3 ) {
        zbar = cCoords.at(3);
    }

    for ( int i = 1; i <= nNodes; i++ ) {
        Node *n = domain->giveNode(i);
        FloatArray *coords = n->giveCoordinates();
        Dof *d1 = n->giveDofWithID( this->dofs(0) );
        Dof *d2 = n->giveDofWithID( this->dofs(1) );
        int k1 = d1->__givePrescribedEquationNumber();
        int k2 = d2->__givePrescribedEquationNumber();
        if ( nsd == 2 ) {
            if ( k1 ) {
                C.at(k1, 1) = coords->at(1) - xbar;
                C.at(k1, 3) = coords->at(2) - ybar;
            }

            if ( k2 ) {
                C.at(k2, 2) = coords->at(2) - ybar;
                C.at(k2, 3) = coords->at(1) - xbar;
            }
        } else { // nsd == 3
            OOFEM_ERROR("PrescribedGenStrainShell7 :: updateCoefficientMatrix - 3D Not tested yet!");
            Dof *d3 = n->giveDofWithID( this->dofs(2) );
            int k3 = d3->__givePrescribedEquationNumber();

            if ( k1 ) {
                C.at(k1, 1) = coords->at(1) - xbar;
                C.at(k1, 4) = coords->at(2) - ybar;
                C.at(k1, 5) = coords->at(3) - zbar;
            }
            if ( k2 ) {
                C.at(k2, 2) = coords->at(2) - ybar;
                C.at(k2, 4) = coords->at(1) - xbar;
                C.at(k2, 6) = coords->at(3) - zbar;
            }
            if ( k3 ) {
                C.at(k3, 3) = coords->at(3) - zbar;
                C.at(k3, 5) = coords->at(1) - xbar;
                C.at(k3, 6) = coords->at(2) - ybar;
            }
        }
    }
}


double PrescribedGenStrainShell7 :: domainSize()
{
    int nsd = this->domain->giveNumberOfSpatialDimensions();
    double domain_size = 0.0;
    // This requires the boundary to be consistent and ordered correctly.
    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);
        FEInterpolation *fei = e->giveInterpolation();
        domain_size += fei->evalNXIntegral( boundary, FEIElementGeometryWrapper(e) );
    }
    return domain_size / nsd;
}


void PrescribedGenStrainShell7 :: computeField(FloatArray &sigma, EquationID eid, TimeStep *tStep)
{
    EngngModel *emodel = this->domain->giveEngngModel();
    int npeq = emodel->giveNumberOfDomainEquations( this->giveDomain()->giveNumber(), EModelDefaultPrescribedEquationNumbering() );
    FloatArray R_c(npeq), R_ext(npeq);

    R_c.zero();
    R_ext.zero();
    emodel->assembleVector( R_c, tStep, eid, InternalForceAssembler(), VM_Total,
                            EModelDefaultPrescribedEquationNumbering(), this->giveDomain() );
    emodel->assembleVector( R_ext, tStep, eid, ExternalForceAssembler(), VM_Total,
                            EModelDefaultPrescribedEquationNumbering(), this->giveDomain() );
    R_c.subtract(R_ext);

    // Condense it;
    FloatMatrix C;
    this->updateCoefficientMatrix(C);
    sigma.beTProductOf(C, R_c);
    sigma.times( 1. / this->domainSize() );
}


void PrescribedGenStrainShell7 :: computeTangent(FloatMatrix &tangent, EquationID eid, TimeStep *tStep)
// a = [a_c; a_f];
// K.a = [R_c,0];
// [K_cc, K_cf; K_fc, K_ff].[a_c;a_f] = [R_c; 0];
// a_c = d.[x-x_b] = [x-x_b].d = C.d
// E = C'.(K_cc - K_cf.K_ff^(-1).K_fc).C
//   = C'.(K_cc.C - K_cf.(K_ff^(-1).(K_fc.C)))
//   = C'.(K_cc.C - K_cf.a)
//   = C'.X
{
    // Fetch some information from the engineering model
    EngngModel *rve = this->giveDomain()->giveEngngModel();
    ///@todo Get this from engineering model
    std :: unique_ptr< SparseLinearSystemNM > solver( classFactory.createSparseLinSolver( ST_Petsc, this->domain, this->domain->giveEngngModel() ) ); // = rve->giveLinearSolver();
    SparseMtrxType stype = solver->giveRecommendedMatrix(true);
    EModelDefaultEquationNumbering fnum;
    EModelDefaultPrescribedEquationNumbering pnum;

    // Set up and assemble tangent FE-matrix which will make up the sensitivity analysis for the macroscopic material tangent.
    std :: unique_ptr< SparseMtrx > Kff( classFactory.createSparseMtrx(stype) );
    std :: unique_ptr< SparseMtrx > Kfp( classFactory.createSparseMtrx(stype) );
    std :: unique_ptr< SparseMtrx > Kpf( classFactory.createSparseMtrx(stype) );
    std :: unique_ptr< SparseMtrx > Kpp( classFactory.createSparseMtrx(stype) );
    if ( !Kff ) {
        OOFEM_ERROR("MixedGradientPressureBC :: computeTangents - Couldn't create sparse matrix of type %d\n", stype);
    }
    Kff->buildInternalStructure(rve, 1, eid, fnum);
    Kfp->buildInternalStructure(rve, 1, eid, fnum, pnum);
    Kpf->buildInternalStructure(rve, 1, eid, pnum, fnum);
    Kpp->buildInternalStructure(rve, 1, eid, pnum);
    rve->assemble(*Kff, tStep, eid, TangentStiffnessMatrix, fnum, this->domain);
    rve->assemble(*Kfp, tStep, eid, TangentStiffnessMatrix, fnum, pnum, this->domain);
    rve->assemble(*Kpf, tStep, eid, TangentStiffnessMatrix, pnum, fnum, this->domain);
    rve->assemble(*Kpp, tStep, eid, TangentStiffnessMatrix, pnum, this->domain);

    FloatMatrix C, X, Kpfa, KfpC, a;

    this->updateCoefficientMatrix(C);
    Kpf->timesT(C, KfpC);
    solver->solve(*Kff, KfpC, a);
    Kpp->times(C, X);
    Kpf->times(a, Kpfa);
    X.subtract(Kpfa);
    tangent.beTProductOf(C, X);
    tangent.times( 1. / this->domainSize() );
}


IRResultType PrescribedGenStrainShell7 :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->gradient, _IFT_PrescribedGenStrainShell7_gradient);

    this->centerCoord.resize( this->gradient.giveNumberOfColumns() );
    this->centerCoord.zero();
    IR_GIVE_OPTIONAL_FIELD(ir, this->centerCoord, _IFT_PrescribedGenStrainShell7_centercoords)

    return GeneralBoundaryCondition :: initializeFrom(ir);
}


void PrescribedGenStrainShell7 :: giveInputRecord(DynamicInputRecord &input)
{
    BoundaryCondition :: giveInputRecord(input);
    input.setField(this->gradient, _IFT_PrescribedGenStrainShell7_gradient);
    input.setField(this->centerCoord, _IFT_PrescribedGenStrainShell7_centercoords);
}
} // end namespace oofem
