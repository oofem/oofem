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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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
#include "sparsemtrx.h"

#include "line2boundaryelement.h"

#include "tr21stokes.h"

namespace oofem {
StokesFlowStressHomogenization :: StokesFlowStressHomogenization(int i, EngngModel *_master) : StokesFlow(i, _master)
{
    this->K_cf = NULL;
    this->K_cc = NULL;
    this->linNumericalMethod = NULL;
    this->C.beEmptyMtrx();
}

StokesFlowStressHomogenization :: ~StokesFlowStressHomogenization()
{
    if ( this->K_cf ) {
        delete this->K_cf;
    }
    if ( this->K_cc ) {
        delete this->K_cc;
    }
    if ( this->linNumericalMethod ) {
        delete this->linNumericalMethod;
    }
}

SparseLinearSystemNM *StokesFlowStressHomogenization :: giveLinearNumericalMethod()
{
    if (!this->linNumericalMethod) {
        this->linNumericalMethod = CreateUsrDefSparseLinSolver(ST_Petsc, 1, this->giveDomain(1), this);
    }
    return this->linNumericalMethod;
}

PrescribedGradient *StokesFlowStressHomogenization :: giveDirichletBC()
{
    PrescribedGradient *pt = dynamic_cast< PrescribedGradient * >( this->giveDomain(1)->giveBc(1) );
    if ( pt == NULL ) {
        OOFEM_ERROR("StokesFlowStressHomogenization :: computeMacroStress - Wrong boundary conditions");
    }
    return pt;
}

bool StokesFlowStressHomogenization :: computeMacroStress(FloatArray &answer, const FloatArray &input, TimeStep *tStep)
{
    this->giveDirichletBC()->setPrescribedGradientVoigt(input);
    this->solveYourselfAt(tStep);
    this->computeMacroStressFromDirichlet(answer, tStep);
    return true;
}

void StokesFlowStressHomogenization :: computeMacroTangent(FloatMatrix &answer, TimeStep *tStep)
{
    this->computeMacroStressTangentFromDirichlet(answer, tStep);
}

double StokesFlowStressHomogenization :: computeSize(int di)
{
    // This could be based on the topology description as well.
    Domain *d = this->giveDomain(di);
    int nsd = d->giveNumberOfSpatialDimensions();
    double domain_size = 0.0;

    // This requires the boundary to be consistent and ordered correctly.
    for (int i = 1; i <= d->giveNumberOfElements(); ++i) {
        //BoundaryElement *e = dynamic_cast< BoundaryElement* >(d->giveElement(i));
        Line2BoundaryElement *e = dynamic_cast< Line2BoundaryElement* >(d->giveElement(i));
        if (e) {
            domain_size += e->computeNXIntegral();
        }
    }
    return domain_size/nsd;
}

int StokesFlowStressHomogenization :: forceEquationNumbering(int di)
{
    int rs = StokesFlow :: forceEquationNumbering(di);
    // Then mesh is replaced, all the stored structures are trashed.
    if ( this->K_cf ) {
        delete this->K_cf;
        this->K_cf = NULL;
    }
    if ( this->K_cc ) {
        delete this->K_cc;
        this->K_cc = NULL;
    }
    this->C.beEmptyMtrx();
    return rs;
}

void StokesFlowStressHomogenization :: updateYourself(TimeStep *tStep)
{
    StokesFlow :: updateYourself(tStep);
    this->C.beEmptyMtrx();
}

void StokesFlowStressHomogenization :: computeMacroStressFromDirichlet(FloatArray &answer, TimeStep *tStep)
// stress = C'.R_c
{
    double volume;
    int npeq = this->giveNumberOfPrescribedEquations(EID_MomentumBalance_ConservationEquation);
    FloatArray R_c(npeq), R_ext(npeq);
    R_c.zero();
    R_ext.zero();
    this->assembleVectorFromElements( R_c, tStep, EID_MomentumBalance_ConservationEquation, NodalInternalForcesVector, VM_Total,
                                     EModelDefaultPrescribedEquationNumbering(), this->giveDomain(1) );
    this->assembleVectorFromElements( R_ext, tStep, EID_MomentumBalance_ConservationEquation, LoadVector, VM_Total,
                                     EModelDefaultPrescribedEquationNumbering(), this->giveDomain(1) );
    R_c.subtract(R_ext);
    if (!this->C.isNotEmpty()) {
        this->giveDirichletBC()->updateCoefficientMatrix(this->C);
    }
    answer.beTProductOf(this->C, R_c);
    volume = this->computeSize(1);
    answer.times(1/volume);

    // DEBUGGING;
//    FloatArray stress(3); stress.zero();
//    for (int i = 1; i <= this->giveDomain(1)->giveNumberOfElements(); ++i) {
//        Tr21Stokes *e = dynamic_cast<Tr21Stokes*>(this->giveDomain(1)->giveElement(i));
//        if (e) e->computeStressIntegral(stress, tStep);
//    }
//    stress.times(1/volume);
//    answer.printYourself();
//    stress.printYourself();
//    OOFEM_ERROR("DEBUG QUIT");
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
    double volume;
    // The matrices at the last step (if the mesh is changed, equations renumbered, these should be removed.)
    if (!this->C.isNotEmpty()) {
        this->giveDirichletBC()->updateCoefficientMatrix(this->C);
    }
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

    FloatMatrix X, K_cfa, K_fcC, a;
    this->K_cf->timesT(this->C, K_fcC);
    this->giveLinearNumericalMethod()->solve(this->stiffnessMatrix, K_fcC, a);
    this->K_cc->times(this->C, X);
    this->K_cf->times(a, K_cfa);
    X.subtract(K_cfa);
    tangent.beTProductOf(this->C, X);
    volume = this->computeSize(1);
    tangent.times(1/volume);
}

} // end namespace oofem
