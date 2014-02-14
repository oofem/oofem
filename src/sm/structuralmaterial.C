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

#include "structuralmaterial.h"
#include "domain.h"
#include "verbose.h"
#include "structuralms.h"
#include "structuralelement.h"
#include "nlstructuralelement.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "engngm.h"
#include "fieldmanager.h"
#include "dynamicinputrecord.h"
#include "eleminterpmapperinterface.h"

namespace oofem {
int
StructuralMaterial :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    return mode == _3dMat || mode == _PlaneStress || mode == _PlaneStrain  || mode == _1dMat ||
           mode == _PlateLayer || mode == _2dBeamLayer || mode == _Fiber;
}


void
StructuralMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    ///@todo Move this to StructuralCrossSection ?
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _3dMat ) {
        this->giveRealStressVector_3d(answer, gp, reducedStrain, tStep);
    } else if ( mode == _PlaneStrain ) {
        this->giveRealStressVector_PlaneStrain(answer, gp, reducedStrain, tStep);
    } else if ( mode == _PlaneStress ) {
        this->giveRealStressVector_PlaneStress(answer, gp, reducedStrain, tStep);
    } else if ( mode == _1dMat ) {
        this->giveRealStressVector_1d(answer, gp, reducedStrain, tStep);
    } else if ( mode == _2dBeamLayer ) {
        this->giveRealStressVector_2dBeamLayer(answer, gp, reducedStrain, tStep);
    } else if ( mode == _PlateLayer ) {
        this->giveRealStressVector_PlateLayer(answer, gp, reducedStrain, tStep);
    } else if ( mode == _Fiber ) {
        this->giveRealStressVector_Fiber(answer, gp, reducedStrain, tStep);
    }
}


void
StructuralMaterial :: giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    OOFEM_ERROR2( "%s :: giveRealStressVector_3d - 3d mode not supported", this->giveClassName() );
}


void
StructuralMaterial :: giveRealStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    FloatArray vE, vS;
    StructuralMaterial :: giveFullSymVectorForm(vE, reducedStrain, _PlaneStrain);
    this->giveRealStressVector_3d(vS, gp, vE, tStep);
    StructuralMaterial :: giveReducedSymVectorForm(answer, vS, _PlaneStrain);
}


void
StructuralMaterial :: giveRealStressVector_StressControl(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, const IntArray &strainControl, TimeStep *tStep)
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    IntArray stressControl;
    FloatArray vE, increment_vE, vS, reducedvS;
    FloatMatrix tangent, reducedTangent;
    // Iterate to find full vE.
    // Compute the negated the array of control since we need stressControl as well;

    stressControl.resize( 9 - strainControl.giveSize() );
    for ( int i = 1, j = 1; i <= 9; i++ ) {
        if ( !strainControl.contains(i) ) {
            stressControl.at(j++) = i;
        }
    }

    // Initial guess;
    vE = status->giveStrainVector();
    for ( int i = 1; i <= strainControl.giveSize(); ++i ) {
        vE.at( strainControl.at(i) ) = reducedStrain.at(i);
    }

    // Iterate to find full vE.
    for ( int k = 0; k < 100; k++ ) { // Allow for a generous 100 iterations.
        this->giveRealStressVector_3d(vS, gp, vE, tStep);
        vS.printYourself();
        reducedvS.beSubArrayOf(vS, stressControl);
        if ( reducedvS.computeNorm() < 1e-6 ) { ///@todo We need a tolerance here!
            StructuralMaterial :: giveReducedVectorForm(answer, vS, _1dMat);
            return;
        }

        this->give3dMaterialStiffnessMatrix(tangent, TangentStiffness, gp, tStep);
        reducedTangent.beSubMatrixOf(tangent, stressControl, stressControl);
        reducedTangent.solveForRhs(reducedvS, increment_vE);
        vE.assemble(increment_vE, stressControl);
    }

    OOFEM_WARNING("StructuralMaterial :: giveRealStressVector_StressControl - Iteration did not converge");
    answer.clear();
}


void
StructuralMaterial :: giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    IntArray strainControl;
    StructuralMaterial :: giveVoigtSymVectorMask(strainControl, _PlaneStress);
    this->giveRealStressVector_StressControl(answer, gp, reducedStrain, strainControl, tStep);
}


void
StructuralMaterial :: giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    IntArray strainControl;
    StructuralMaterial :: giveVoigtSymVectorMask(strainControl, _1dMat);
    this->giveRealStressVector_StressControl(answer, gp, reducedStrain, strainControl, tStep);
}


void
StructuralMaterial :: giveRealStressVector_2dBeamLayer(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    IntArray strainControl;
    StructuralMaterial :: giveVoigtSymVectorMask(strainControl, _2dBeamLayer);
    this->giveRealStressVector_StressControl(answer, gp, reducedStrain, strainControl, tStep);
}


void
StructuralMaterial :: giveRealStressVector_PlateLayer(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    IntArray strainControl;
    StructuralMaterial :: giveVoigtSymVectorMask(strainControl, _PlateLayer);
    this->giveRealStressVector_StressControl(answer, gp, reducedStrain, strainControl, tStep);
}


void
StructuralMaterial :: giveRealStressVector_Fiber(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    IntArray strainControl;
    StructuralMaterial :: giveVoigtSymVectorMask(strainControl, _Fiber);
    this->giveRealStressVector_StressControl(answer, gp, reducedStrain, strainControl, tStep);
}


void
StructuralMaterial :: giveFirstPKStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &vF, TimeStep *tStep)
{
    // Default implementation used if this method is not overloaded by the particular material model.
    // 1) Compute Green-Lagrange strain and call standard method for small strains.
    // 2) Treat stress as second Piola-Kirchhoff stress and convert to first Piola-Kirchhoff stress.
    // 3) Set state variables F, P

    FloatArray vE, vS;
    FloatMatrix F, E;
    F.beMatrixForm(vF);
    E.beTProductOf(F, F);
    E.at(1, 1) -= 1.0;
    E.at(2, 2) -= 1.0;
    E.at(3, 3) -= 1.0;
    E.times(0.5);
    vE.beSymVectorFormOfStrain(E);      // 6

    ///@todo Have this function:
    this->giveRealStressVector_3d(vS, gp, vE, tStep);
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    // Compute first PK stress from second PK stress
    FloatMatrix P, S;
    S.beMatrixForm(vS);
    P.beProductOf(F, S);
    answer.beVectorForm(P);

    status->letTempPVectorBe(answer);
    status->letTempFVectorBe(vF);
}


void
StructuralMaterial :: giveFirstPKStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedvF, TimeStep *tStep)
{
    FloatArray vF, vP;
    StructuralMaterial :: giveFullVectorFormF(vF, reducedvF, _PlaneStrain);
    this->giveFirstPKStressVector_3d(vP, gp, vF, tStep);
    StructuralMaterial :: giveReducedVectorForm(answer, vP, _PlaneStrain);
}


void
StructuralMaterial :: giveFirstPKStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedvF, TimeStep *tStep)
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    IntArray F_control, P_control; // Determines which components are controlled by F and P resp.
    FloatArray vF, increment_vF, vP, vP_control;
    FloatMatrix tangent, tangent_Pcontrol;
    // Iterate to find full vF.
    StructuralMaterial :: giveVoigtVectorMask(F_control, _PlaneStress);
    // Compute the negated the array of control since we need P_control as well;
    P_control.resize( 9 - F_control.giveSize() );
    for ( int i = 1, j = 1; i <= 9; i++ ) {
        if ( !F_control.contains(i) ) {
            P_control.at(j++) = i;
        }
    }

    // Initial guess;
    vF = status->giveFVector();
    for ( int i = 1; i <= F_control.giveSize(); ++i ) {
        vF.at( F_control.at(i) ) = reducedvF.at(i);
    }

    // Iterate to find full vF.
    for ( int k = 0; k < 100; k++ ) { // Allow for a generous 100 iterations.
        this->giveFirstPKStressVector_3d(vP, gp, vF, tStep);
        vP_control.beSubArrayOf(vP, P_control);
        if ( vP_control.computeNorm() < 1e-6 ) { ///@todo We need a tolerance here!
            StructuralMaterial :: giveReducedVectorForm(answer, vP, _1dMat);
            return;
        }

        this->give3dMaterialStiffnessMatrix_dPdF(tangent, TangentStiffness, gp, tStep);
        tangent_Pcontrol.beSubMatrixOf(tangent, P_control, P_control);
        tangent_Pcontrol.solveForRhs(vP_control, increment_vF);
        vF.assemble(increment_vF, P_control);
    }

    OOFEM_WARNING("StructuralMaterial :: giveFirstPKStressVector_PlaneStress - Iteration did not converge");
    answer.clear();
}


void
StructuralMaterial :: giveFirstPKStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedvF, TimeStep *tStep)
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    IntArray P_control; // Determines which components are controlled by P resp.
    FloatArray old_vF, vF, increment_vF, vP, vP_control;
    FloatMatrix tangent, tangent_Pcontrol;
    // Compute the negated the array of control since we need P_control as well;
    P_control.resize(8);
    for ( int i = 1; i <= 8; i++ ) {
        P_control.at(i) = i + 1;
    }

    // Initial guess;
    vF = status->giveFVector();
    vF.at(1) = reducedvF.at(1);
    // Iterate to find full vF.
    for ( int k = 0; k < 100; k++ ) { // Allow for a generous 100 iterations.
        this->giveFirstPKStressVector_3d(vP, gp, vF, tStep);
        vP_control.beSubArrayOf(vP, P_control);
        if ( vP_control.computeNorm() < 1e-6 ) { ///@todo We need a tolerance here!
            StructuralMaterial :: giveReducedVectorForm(answer, vP, _1dMat);
            return;
        }

        this->give3dMaterialStiffnessMatrix_dPdF(tangent, TangentStiffness, gp, tStep);
        tangent_Pcontrol.beSubMatrixOf(tangent, P_control, P_control);
        tangent_Pcontrol.solveForRhs(vP_control, increment_vF);
        vF.assemble(increment_vF, P_control);
    }

    OOFEM_WARNING("StructuralMaterial :: giveFirstPKStressVector_1d - Iteration did not converge");
    answer.clear();
}


void
StructuralMaterial :: convert_dSdE_2_dPdF(FloatMatrix &answer, const FloatMatrix &C, FloatArray &S, FloatArray &F, MaterialMode matMode)
{
    // Converts the reduced dSdE-stiffness to reduced dPdF-sitiffness for different MaterialModes
    // Performs the following operation dPdF = I_ik * S_jl + F_im F_kn C_mjnl,
    // See for example: G.A. Holzapfel, Nonlinear Solid Mechanics: A Continuum Approach for
    // Engineering, 2000, ISBN-10: 0471823198.

    if ( matMode == _3dMat ) {
        //Save terms associated with H = [du/dx, dv/dy, dw/dz, dv/dz, du/dz, du/dy, dw/dy, dw/dx, dv/dx]

        answer.resize(9, 9);
        answer(0, 0) = F(0) * C(0, 0) * F(0) + F(0) * C(0, 5) * F(5) + F(0) * C(0, 4) * F(4) + F(5) * C(5, 0) * F(0) + F(5) * C(5, 5) * F(5) + F(5) * C(5, 4) * F(4) + F(4) * C(4, 0) * F(0) + F(4) * C(4, 5) * F(5) + F(4) * C(4, 4) * F(4) + S(0);
        answer(0, 1) = F(0) * C(0, 5) * F(8) + F(0) * C(0, 1) * F(1) + F(0) * C(0, 3) * F(3) + F(5) * C(5, 5) * F(8) + F(5) * C(5, 1) * F(1) + F(5) * C(5, 3) * F(3) + F(4) * C(4, 5) * F(8) + F(4) * C(4, 1) * F(1) + F(4) * C(4, 3) * F(3) + 0.0;
        answer(0, 2) = F(0) * C(0, 4) * F(7) + F(0) * C(0, 3) * F(6) + F(0) * C(0, 2) * F(2) + F(5) * C(5, 4) * F(7) + F(5) * C(5, 3) * F(6) + F(5) * C(5, 2) * F(2) + F(4) * C(4, 4) * F(7) + F(4) * C(4, 3) * F(6) + F(4) * C(4, 2) * F(2) + 0.0;
        answer(0, 3) = F(0) * C(0, 4) * F(8) + F(0) * C(0, 3) * F(1) + F(0) * C(0, 2) * F(3) + F(5) * C(5, 4) * F(8) + F(5) * C(5, 3) * F(1) + F(5) * C(5, 2) * F(3) + F(4) * C(4, 4) * F(8) + F(4) * C(4, 3) * F(1) + F(4) * C(4, 2) * F(3) + 0.0;
        answer(0, 4) = F(0) * C(0, 4) * F(0) + F(0) * C(0, 3) * F(5) + F(0) * C(0, 2) * F(4) + F(5) * C(5, 4) * F(0) + F(5) * C(5, 3) * F(5) + F(5) * C(5, 2) * F(4) + F(4) * C(4, 4) * F(0) + F(4) * C(4, 3) * F(5) + F(4) * C(4, 2) * F(4) + S(4);
        answer(0, 5) = F(0) * C(0, 5) * F(0) + F(0) * C(0, 1) * F(5) + F(0) * C(0, 3) * F(4) + F(5) * C(5, 5) * F(0) + F(5) * C(5, 1) * F(5) + F(5) * C(5, 3) * F(4) + F(4) * C(4, 5) * F(0) + F(4) * C(4, 1) * F(5) + F(4) * C(4, 3) * F(4) + S(5);
        answer(0, 6) = F(0) * C(0, 5) * F(7) + F(0) * C(0, 1) * F(6) + F(0) * C(0, 3) * F(2) + F(5) * C(5, 5) * F(7) + F(5) * C(5, 1) * F(6) + F(5) * C(5, 3) * F(2) + F(4) * C(4, 5) * F(7) + F(4) * C(4, 1) * F(6) + F(4) * C(4, 3) * F(2) + 0.0;
        answer(0, 7) = F(0) * C(0, 0) * F(7) + F(0) * C(0, 5) * F(6) + F(0) * C(0, 4) * F(2) + F(5) * C(5, 0) * F(7) + F(5) * C(5, 5) * F(6) + F(5) * C(5, 4) * F(2) + F(4) * C(4, 0) * F(7) + F(4) * C(4, 5) * F(6) + F(4) * C(4, 4) * F(2) + 0.0;
        answer(0, 8) = F(0) * C(0, 0) * F(8) + F(0) * C(0, 5) * F(1) + F(0) * C(0, 4) * F(3) + F(5) * C(5, 0) * F(8) + F(5) * C(5, 5) * F(1) + F(5) * C(5, 4) * F(3) + F(4) * C(4, 0) * F(8) + F(4) * C(4, 5) * F(1) + F(4) * C(4, 4) * F(3) + 0.0;
        answer(1, 0) = F(8) * C(5, 0) * F(0) + F(8) * C(5, 5) * F(5) + F(8) * C(5, 4) * F(4) + F(1) * C(1, 0) * F(0) + F(1) * C(1, 5) * F(5) + F(1) * C(1, 4) * F(4) + F(3) * C(3, 0) * F(0) + F(3) * C(3, 5) * F(5) + F(3) * C(3, 4) * F(4) + 0.0;
        answer(1, 1) = F(8) * C(5, 5) * F(8) + F(8) * C(5, 1) * F(1) + F(8) * C(5, 3) * F(3) + F(1) * C(1, 5) * F(8) + F(1) * C(1, 1) * F(1) + F(1) * C(1, 3) * F(3) + F(3) * C(3, 5) * F(8) + F(3) * C(3, 1) * F(1) + F(3) * C(3, 3) * F(3) + S(1);
        answer(1, 2) = F(8) * C(5, 4) * F(7) + F(8) * C(5, 3) * F(6) + F(8) * C(5, 2) * F(2) + F(1) * C(1, 4) * F(7) + F(1) * C(1, 3) * F(6) + F(1) * C(1, 2) * F(2) + F(3) * C(3, 4) * F(7) + F(3) * C(3, 3) * F(6) + F(3) * C(3, 2) * F(2) + 0.0;
        answer(1, 3) = F(8) * C(5, 4) * F(8) + F(8) * C(5, 3) * F(1) + F(8) * C(5, 2) * F(3) + F(1) * C(1, 4) * F(8) + F(1) * C(1, 3) * F(1) + F(1) * C(1, 2) * F(3) + F(3) * C(3, 4) * F(8) + F(3) * C(3, 3) * F(1) + F(3) * C(3, 2) * F(3) + S(3);
        answer(1, 4) = F(8) * C(5, 4) * F(0) + F(8) * C(5, 3) * F(5) + F(8) * C(5, 2) * F(4) + F(1) * C(1, 4) * F(0) + F(1) * C(1, 3) * F(5) + F(1) * C(1, 2) * F(4) + F(3) * C(3, 4) * F(0) + F(3) * C(3, 3) * F(5) + F(3) * C(3, 2) * F(4) + 0.0;
        answer(1, 5) = F(8) * C(5, 5) * F(0) + F(8) * C(5, 1) * F(5) + F(8) * C(5, 3) * F(4) + F(1) * C(1, 5) * F(0) + F(1) * C(1, 1) * F(5) + F(1) * C(1, 3) * F(4) + F(3) * C(3, 5) * F(0) + F(3) * C(3, 1) * F(5) + F(3) * C(3, 3) * F(4) + 0.0;
        answer(1, 6) = F(8) * C(5, 5) * F(7) + F(8) * C(5, 1) * F(6) + F(8) * C(5, 3) * F(2) + F(1) * C(1, 5) * F(7) + F(1) * C(1, 1) * F(6) + F(1) * C(1, 3) * F(2) + F(3) * C(3, 5) * F(7) + F(3) * C(3, 1) * F(6) + F(3) * C(3, 3) * F(2) + 0.0;
        answer(1, 7) = F(8) * C(5, 0) * F(7) + F(8) * C(5, 5) * F(6) + F(8) * C(5, 4) * F(2) + F(1) * C(1, 0) * F(7) + F(1) * C(1, 5) * F(6) + F(1) * C(1, 4) * F(2) + F(3) * C(3, 0) * F(7) + F(3) * C(3, 5) * F(6) + F(3) * C(3, 4) * F(2) + 0.0;
        answer(1, 8) = F(8) * C(5, 0) * F(8) + F(8) * C(5, 5) * F(1) + F(8) * C(5, 4) * F(3) + F(1) * C(1, 0) * F(8) + F(1) * C(1, 5) * F(1) + F(1) * C(1, 4) * F(3) + F(3) * C(3, 0) * F(8) + F(3) * C(3, 5) * F(1) + F(3) * C(3, 4) * F(3) + S(5);
        answer(2, 0) = F(7) * C(4, 0) * F(0) + F(7) * C(4, 5) * F(5) + F(7) * C(4, 4) * F(4) + F(6) * C(3, 0) * F(0) + F(6) * C(3, 5) * F(5) + F(6) * C(3, 4) * F(4) + F(2) * C(2, 0) * F(0) + F(2) * C(2, 5) * F(5) + F(2) * C(2, 4) * F(4) + 0.0;
        answer(2, 1) = F(7) * C(4, 5) * F(8) + F(7) * C(4, 1) * F(1) + F(7) * C(4, 3) * F(3) + F(6) * C(3, 5) * F(8) + F(6) * C(3, 1) * F(1) + F(6) * C(3, 3) * F(3) + F(2) * C(2, 5) * F(8) + F(2) * C(2, 1) * F(1) + F(2) * C(2, 3) * F(3) + 0.0;
        answer(2, 2) = F(7) * C(4, 4) * F(7) + F(7) * C(4, 3) * F(6) + F(7) * C(4, 2) * F(2) + F(6) * C(3, 4) * F(7) + F(6) * C(3, 3) * F(6) + F(6) * C(3, 2) * F(2) + F(2) * C(2, 4) * F(7) + F(2) * C(2, 3) * F(6) + F(2) * C(2, 2) * F(2) + S(2);
        answer(2, 3) = F(7) * C(4, 4) * F(8) + F(7) * C(4, 3) * F(1) + F(7) * C(4, 2) * F(3) + F(6) * C(3, 4) * F(8) + F(6) * C(3, 3) * F(1) + F(6) * C(3, 2) * F(3) + F(2) * C(2, 4) * F(8) + F(2) * C(2, 3) * F(1) + F(2) * C(2, 2) * F(3) + 0.0;
        answer(2, 4) = F(7) * C(4, 4) * F(0) + F(7) * C(4, 3) * F(5) + F(7) * C(4, 2) * F(4) + F(6) * C(3, 4) * F(0) + F(6) * C(3, 3) * F(5) + F(6) * C(3, 2) * F(4) + F(2) * C(2, 4) * F(0) + F(2) * C(2, 3) * F(5) + F(2) * C(2, 2) * F(4) + 0.0;
        answer(2, 5) = F(7) * C(4, 5) * F(0) + F(7) * C(4, 1) * F(5) + F(7) * C(4, 3) * F(4) + F(6) * C(3, 5) * F(0) + F(6) * C(3, 1) * F(5) + F(6) * C(3, 3) * F(4) + F(2) * C(2, 5) * F(0) + F(2) * C(2, 1) * F(5) + F(2) * C(2, 3) * F(4) + 0.0;
        answer(2, 6) = F(7) * C(4, 5) * F(7) + F(7) * C(4, 1) * F(6) + F(7) * C(4, 3) * F(2) + F(6) * C(3, 5) * F(7) + F(6) * C(3, 1) * F(6) + F(6) * C(3, 3) * F(2) + F(2) * C(2, 5) * F(7) + F(2) * C(2, 1) * F(6) + F(2) * C(2, 3) * F(2) + S(3);
        answer(2, 7) = F(7) * C(4, 0) * F(7) + F(7) * C(4, 5) * F(6) + F(7) * C(4, 4) * F(2) + F(6) * C(3, 0) * F(7) + F(6) * C(3, 5) * F(6) + F(6) * C(3, 4) * F(2) + F(2) * C(2, 0) * F(7) + F(2) * C(2, 5) * F(6) + F(2) * C(2, 4) * F(2) + S(4);
        answer(2, 8) = F(7) * C(4, 0) * F(8) + F(7) * C(4, 5) * F(1) + F(7) * C(4, 4) * F(3) + F(6) * C(3, 0) * F(8) + F(6) * C(3, 5) * F(1) + F(6) * C(3, 4) * F(3) + F(2) * C(2, 0) * F(8) + F(2) * C(2, 5) * F(1) + F(2) * C(2, 4) * F(3) + 0.0;
        answer(3, 0) = F(8) * C(4, 0) * F(0) + F(8) * C(4, 5) * F(5) + F(8) * C(4, 4) * F(4) + F(1) * C(3, 0) * F(0) + F(1) * C(3, 5) * F(5) + F(1) * C(3, 4) * F(4) + F(3) * C(2, 0) * F(0) + F(3) * C(2, 5) * F(5) + F(3) * C(2, 4) * F(4) + 0.0;
        answer(3, 1) = F(8) * C(4, 5) * F(8) + F(8) * C(4, 1) * F(1) + F(8) * C(4, 3) * F(3) + F(1) * C(3, 5) * F(8) + F(1) * C(3, 1) * F(1) + F(1) * C(3, 3) * F(3) + F(3) * C(2, 5) * F(8) + F(3) * C(2, 1) * F(1) + F(3) * C(2, 3) * F(3) + S(3);
        answer(3, 2) = F(8) * C(4, 4) * F(7) + F(8) * C(4, 3) * F(6) + F(8) * C(4, 2) * F(2) + F(1) * C(3, 4) * F(7) + F(1) * C(3, 3) * F(6) + F(1) * C(3, 2) * F(2) + F(3) * C(2, 4) * F(7) + F(3) * C(2, 3) * F(6) + F(3) * C(2, 2) * F(2) + 0.0;
        answer(3, 3) = F(8) * C(4, 4) * F(8) + F(8) * C(4, 3) * F(1) + F(8) * C(4, 2) * F(3) + F(1) * C(3, 4) * F(8) + F(1) * C(3, 3) * F(1) + F(1) * C(3, 2) * F(3) + F(3) * C(2, 4) * F(8) + F(3) * C(2, 3) * F(1) + F(3) * C(2, 2) * F(3) + S(2);
        answer(3, 4) = F(8) * C(4, 4) * F(0) + F(8) * C(4, 3) * F(5) + F(8) * C(4, 2) * F(4) + F(1) * C(3, 4) * F(0) + F(1) * C(3, 3) * F(5) + F(1) * C(3, 2) * F(4) + F(3) * C(2, 4) * F(0) + F(3) * C(2, 3) * F(5) + F(3) * C(2, 2) * F(4) + 0.0;
        answer(3, 5) = F(8) * C(4, 5) * F(0) + F(8) * C(4, 1) * F(5) + F(8) * C(4, 3) * F(4) + F(1) * C(3, 5) * F(0) + F(1) * C(3, 1) * F(5) + F(1) * C(3, 3) * F(4) + F(3) * C(2, 5) * F(0) + F(3) * C(2, 1) * F(5) + F(3) * C(2, 3) * F(4) + 0.0;
        answer(3, 6) = F(8) * C(4, 5) * F(7) + F(8) * C(4, 1) * F(6) + F(8) * C(4, 3) * F(2) + F(1) * C(3, 5) * F(7) + F(1) * C(3, 1) * F(6) + F(1) * C(3, 3) * F(2) + F(3) * C(2, 5) * F(7) + F(3) * C(2, 1) * F(6) + F(3) * C(2, 3) * F(2) + 0.0;
        answer(3, 7) = F(8) * C(4, 0) * F(7) + F(8) * C(4, 5) * F(6) + F(8) * C(4, 4) * F(2) + F(1) * C(3, 0) * F(7) + F(1) * C(3, 5) * F(6) + F(1) * C(3, 4) * F(2) + F(3) * C(2, 0) * F(7) + F(3) * C(2, 5) * F(6) + F(3) * C(2, 4) * F(2) + 0.0;
        answer(3, 8) = F(8) * C(4, 0) * F(8) + F(8) * C(4, 5) * F(1) + F(8) * C(4, 4) * F(3) + F(1) * C(3, 0) * F(8) + F(1) * C(3, 5) * F(1) + F(1) * C(3, 4) * F(3) + F(3) * C(2, 0) * F(8) + F(3) * C(2, 5) * F(1) + F(3) * C(2, 4) * F(3) + S(4);
        answer(4, 0) = F(0) * C(4, 0) * F(0) + F(0) * C(4, 5) * F(5) + F(0) * C(4, 4) * F(4) + F(5) * C(3, 0) * F(0) + F(5) * C(3, 5) * F(5) + F(5) * C(3, 4) * F(4) + F(4) * C(2, 0) * F(0) + F(4) * C(2, 5) * F(5) + F(4) * C(2, 4) * F(4) + S(4);
        answer(4, 1) = F(0) * C(4, 5) * F(8) + F(0) * C(4, 1) * F(1) + F(0) * C(4, 3) * F(3) + F(5) * C(3, 5) * F(8) + F(5) * C(3, 1) * F(1) + F(5) * C(3, 3) * F(3) + F(4) * C(2, 5) * F(8) + F(4) * C(2, 1) * F(1) + F(4) * C(2, 3) * F(3) + 0.0;
        answer(4, 2) = F(0) * C(4, 4) * F(7) + F(0) * C(4, 3) * F(6) + F(0) * C(4, 2) * F(2) + F(5) * C(3, 4) * F(7) + F(5) * C(3, 3) * F(6) + F(5) * C(3, 2) * F(2) + F(4) * C(2, 4) * F(7) + F(4) * C(2, 3) * F(6) + F(4) * C(2, 2) * F(2) + 0.0;
        answer(4, 3) = F(0) * C(4, 4) * F(8) + F(0) * C(4, 3) * F(1) + F(0) * C(4, 2) * F(3) + F(5) * C(3, 4) * F(8) + F(5) * C(3, 3) * F(1) + F(5) * C(3, 2) * F(3) + F(4) * C(2, 4) * F(8) + F(4) * C(2, 3) * F(1) + F(4) * C(2, 2) * F(3) + 0.0;
        answer(4, 4) = F(0) * C(4, 4) * F(0) + F(0) * C(4, 3) * F(5) + F(0) * C(4, 2) * F(4) + F(5) * C(3, 4) * F(0) + F(5) * C(3, 3) * F(5) + F(5) * C(3, 2) * F(4) + F(4) * C(2, 4) * F(0) + F(4) * C(2, 3) * F(5) + F(4) * C(2, 2) * F(4) + S(2);
        answer(4, 5) = F(0) * C(4, 5) * F(0) + F(0) * C(4, 1) * F(5) + F(0) * C(4, 3) * F(4) + F(5) * C(3, 5) * F(0) + F(5) * C(3, 1) * F(5) + F(5) * C(3, 3) * F(4) + F(4) * C(2, 5) * F(0) + F(4) * C(2, 1) * F(5) + F(4) * C(2, 3) * F(4) + S(3);
        answer(4, 6) = F(0) * C(4, 5) * F(7) + F(0) * C(4, 1) * F(6) + F(0) * C(4, 3) * F(2) + F(5) * C(3, 5) * F(7) + F(5) * C(3, 1) * F(6) + F(5) * C(3, 3) * F(2) + F(4) * C(2, 5) * F(7) + F(4) * C(2, 1) * F(6) + F(4) * C(2, 3) * F(2) + 0.0;
        answer(4, 7) = F(0) * C(4, 0) * F(7) + F(0) * C(4, 5) * F(6) + F(0) * C(4, 4) * F(2) + F(5) * C(3, 0) * F(7) + F(5) * C(3, 5) * F(6) + F(5) * C(3, 4) * F(2) + F(4) * C(2, 0) * F(7) + F(4) * C(2, 5) * F(6) + F(4) * C(2, 4) * F(2) + 0.0;
        answer(4, 8) = F(0) * C(4, 0) * F(8) + F(0) * C(4, 5) * F(1) + F(0) * C(4, 4) * F(3) + F(5) * C(3, 0) * F(8) + F(5) * C(3, 5) * F(1) + F(5) * C(3, 4) * F(3) + F(4) * C(2, 0) * F(8) + F(4) * C(2, 5) * F(1) + F(4) * C(2, 4) * F(3) + 0.0;
        answer(5, 0) = F(0) * C(5, 0) * F(0) + F(0) * C(5, 5) * F(5) + F(0) * C(5, 4) * F(4) + F(5) * C(1, 0) * F(0) + F(5) * C(1, 5) * F(5) + F(5) * C(1, 4) * F(4) + F(4) * C(3, 0) * F(0) + F(4) * C(3, 5) * F(5) + F(4) * C(3, 4) * F(4) + S(5);
        answer(5, 1) = F(0) * C(5, 5) * F(8) + F(0) * C(5, 1) * F(1) + F(0) * C(5, 3) * F(3) + F(5) * C(1, 5) * F(8) + F(5) * C(1, 1) * F(1) + F(5) * C(1, 3) * F(3) + F(4) * C(3, 5) * F(8) + F(4) * C(3, 1) * F(1) + F(4) * C(3, 3) * F(3) + 0.0;
        answer(5, 2) = F(0) * C(5, 4) * F(7) + F(0) * C(5, 3) * F(6) + F(0) * C(5, 2) * F(2) + F(5) * C(1, 4) * F(7) + F(5) * C(1, 3) * F(6) + F(5) * C(1, 2) * F(2) + F(4) * C(3, 4) * F(7) + F(4) * C(3, 3) * F(6) + F(4) * C(3, 2) * F(2) + 0.0;
        answer(5, 3) = F(0) * C(5, 4) * F(8) + F(0) * C(5, 3) * F(1) + F(0) * C(5, 2) * F(3) + F(5) * C(1, 4) * F(8) + F(5) * C(1, 3) * F(1) + F(5) * C(1, 2) * F(3) + F(4) * C(3, 4) * F(8) + F(4) * C(3, 3) * F(1) + F(4) * C(3, 2) * F(3) + 0.0;
        answer(5, 4) = F(0) * C(5, 4) * F(0) + F(0) * C(5, 3) * F(5) + F(0) * C(5, 2) * F(4) + F(5) * C(1, 4) * F(0) + F(5) * C(1, 3) * F(5) + F(5) * C(1, 2) * F(4) + F(4) * C(3, 4) * F(0) + F(4) * C(3, 3) * F(5) + F(4) * C(3, 2) * F(4) + S(3);
        answer(5, 5) = F(0) * C(5, 5) * F(0) + F(0) * C(5, 1) * F(5) + F(0) * C(5, 3) * F(4) + F(5) * C(1, 5) * F(0) + F(5) * C(1, 1) * F(5) + F(5) * C(1, 3) * F(4) + F(4) * C(3, 5) * F(0) + F(4) * C(3, 1) * F(5) + F(4) * C(3, 3) * F(4) + S(1);
        answer(5, 6) = F(0) * C(5, 5) * F(7) + F(0) * C(5, 1) * F(6) + F(0) * C(5, 3) * F(2) + F(5) * C(1, 5) * F(7) + F(5) * C(1, 1) * F(6) + F(5) * C(1, 3) * F(2) + F(4) * C(3, 5) * F(7) + F(4) * C(3, 1) * F(6) + F(4) * C(3, 3) * F(2) + 0.0;
        answer(5, 7) = F(0) * C(5, 0) * F(7) + F(0) * C(5, 5) * F(6) + F(0) * C(5, 4) * F(2) + F(5) * C(1, 0) * F(7) + F(5) * C(1, 5) * F(6) + F(5) * C(1, 4) * F(2) + F(4) * C(3, 0) * F(7) + F(4) * C(3, 5) * F(6) + F(4) * C(3, 4) * F(2) + 0.0;
        answer(5, 8) = F(0) * C(5, 0) * F(8) + F(0) * C(5, 5) * F(1) + F(0) * C(5, 4) * F(3) + F(5) * C(1, 0) * F(8) + F(5) * C(1, 5) * F(1) + F(5) * C(1, 4) * F(3) + F(4) * C(3, 0) * F(8) + F(4) * C(3, 5) * F(1) + F(4) * C(3, 4) * F(3) + 0.0;
        answer(6, 0) = F(7) * C(5, 0) * F(0) + F(7) * C(5, 5) * F(5) + F(7) * C(5, 4) * F(4) + F(6) * C(1, 0) * F(0) + F(6) * C(1, 5) * F(5) + F(6) * C(1, 4) * F(4) + F(2) * C(3, 0) * F(0) + F(2) * C(3, 5) * F(5) + F(2) * C(3, 4) * F(4) + 0.0;
        answer(6, 1) = F(7) * C(5, 5) * F(8) + F(7) * C(5, 1) * F(1) + F(7) * C(5, 3) * F(3) + F(6) * C(1, 5) * F(8) + F(6) * C(1, 1) * F(1) + F(6) * C(1, 3) * F(3) + F(2) * C(3, 5) * F(8) + F(2) * C(3, 1) * F(1) + F(2) * C(3, 3) * F(3) + 0.0;
        answer(6, 2) = F(7) * C(5, 4) * F(7) + F(7) * C(5, 3) * F(6) + F(7) * C(5, 2) * F(2) + F(6) * C(1, 4) * F(7) + F(6) * C(1, 3) * F(6) + F(6) * C(1, 2) * F(2) + F(2) * C(3, 4) * F(7) + F(2) * C(3, 3) * F(6) + F(2) * C(3, 2) * F(2) + S(3);
        answer(6, 3) = F(7) * C(5, 4) * F(8) + F(7) * C(5, 3) * F(1) + F(7) * C(5, 2) * F(3) + F(6) * C(1, 4) * F(8) + F(6) * C(1, 3) * F(1) + F(6) * C(1, 2) * F(3) + F(2) * C(3, 4) * F(8) + F(2) * C(3, 3) * F(1) + F(2) * C(3, 2) * F(3) + 0.0;
        answer(6, 4) = F(7) * C(5, 4) * F(0) + F(7) * C(5, 3) * F(5) + F(7) * C(5, 2) * F(4) + F(6) * C(1, 4) * F(0) + F(6) * C(1, 3) * F(5) + F(6) * C(1, 2) * F(4) + F(2) * C(3, 4) * F(0) + F(2) * C(3, 3) * F(5) + F(2) * C(3, 2) * F(4) + 0.0;
        answer(6, 5) = F(7) * C(5, 5) * F(0) + F(7) * C(5, 1) * F(5) + F(7) * C(5, 3) * F(4) + F(6) * C(1, 5) * F(0) + F(6) * C(1, 1) * F(5) + F(6) * C(1, 3) * F(4) + F(2) * C(3, 5) * F(0) + F(2) * C(3, 1) * F(5) + F(2) * C(3, 3) * F(4) + 0.0;
        answer(6, 6) = F(7) * C(5, 5) * F(7) + F(7) * C(5, 1) * F(6) + F(7) * C(5, 3) * F(2) + F(6) * C(1, 5) * F(7) + F(6) * C(1, 1) * F(6) + F(6) * C(1, 3) * F(2) + F(2) * C(3, 5) * F(7) + F(2) * C(3, 1) * F(6) + F(2) * C(3, 3) * F(2) + S(1);
        answer(6, 7) = F(7) * C(5, 0) * F(7) + F(7) * C(5, 5) * F(6) + F(7) * C(5, 4) * F(2) + F(6) * C(1, 0) * F(7) + F(6) * C(1, 5) * F(6) + F(6) * C(1, 4) * F(2) + F(2) * C(3, 0) * F(7) + F(2) * C(3, 5) * F(6) + F(2) * C(3, 4) * F(2) + S(5);
        answer(6, 8) = F(7) * C(5, 0) * F(8) + F(7) * C(5, 5) * F(1) + F(7) * C(5, 4) * F(3) + F(6) * C(1, 0) * F(8) + F(6) * C(1, 5) * F(1) + F(6) * C(1, 4) * F(3) + F(2) * C(3, 0) * F(8) + F(2) * C(3, 5) * F(1) + F(2) * C(3, 4) * F(3) + 0.0;
        answer(7, 0) = F(7) * C(0, 0) * F(0) + F(7) * C(0, 5) * F(5) + F(7) * C(0, 4) * F(4) + F(6) * C(5, 0) * F(0) + F(6) * C(5, 5) * F(5) + F(6) * C(5, 4) * F(4) + F(2) * C(4, 0) * F(0) + F(2) * C(4, 5) * F(5) + F(2) * C(4, 4) * F(4) + 0.0;
        answer(7, 1) = F(7) * C(0, 5) * F(8) + F(7) * C(0, 1) * F(1) + F(7) * C(0, 3) * F(3) + F(6) * C(5, 5) * F(8) + F(6) * C(5, 1) * F(1) + F(6) * C(5, 3) * F(3) + F(2) * C(4, 5) * F(8) + F(2) * C(4, 1) * F(1) + F(2) * C(4, 3) * F(3) + 0.0;
        answer(7, 2) = F(7) * C(0, 4) * F(7) + F(7) * C(0, 3) * F(6) + F(7) * C(0, 2) * F(2) + F(6) * C(5, 4) * F(7) + F(6) * C(5, 3) * F(6) + F(6) * C(5, 2) * F(2) + F(2) * C(4, 4) * F(7) + F(2) * C(4, 3) * F(6) + F(2) * C(4, 2) * F(2) + S(4);
        answer(7, 3) = F(7) * C(0, 4) * F(8) + F(7) * C(0, 3) * F(1) + F(7) * C(0, 2) * F(3) + F(6) * C(5, 4) * F(8) + F(6) * C(5, 3) * F(1) + F(6) * C(5, 2) * F(3) + F(2) * C(4, 4) * F(8) + F(2) * C(4, 3) * F(1) + F(2) * C(4, 2) * F(3) + 0.0;
        answer(7, 4) = F(7) * C(0, 4) * F(0) + F(7) * C(0, 3) * F(5) + F(7) * C(0, 2) * F(4) + F(6) * C(5, 4) * F(0) + F(6) * C(5, 3) * F(5) + F(6) * C(5, 2) * F(4) + F(2) * C(4, 4) * F(0) + F(2) * C(4, 3) * F(5) + F(2) * C(4, 2) * F(4) + 0.0;
        answer(7, 5) = F(7) * C(0, 5) * F(0) + F(7) * C(0, 1) * F(5) + F(7) * C(0, 3) * F(4) + F(6) * C(5, 5) * F(0) + F(6) * C(5, 1) * F(5) + F(6) * C(5, 3) * F(4) + F(2) * C(4, 5) * F(0) + F(2) * C(4, 1) * F(5) + F(2) * C(4, 3) * F(4) + 0.0;
        answer(7, 6) = F(7) * C(0, 5) * F(7) + F(7) * C(0, 1) * F(6) + F(7) * C(0, 3) * F(2) + F(6) * C(5, 5) * F(7) + F(6) * C(5, 1) * F(6) + F(6) * C(5, 3) * F(2) + F(2) * C(4, 5) * F(7) + F(2) * C(4, 1) * F(6) + F(2) * C(4, 3) * F(2) + S(5);
        answer(7, 7) = F(7) * C(0, 0) * F(7) + F(7) * C(0, 5) * F(6) + F(7) * C(0, 4) * F(2) + F(6) * C(5, 0) * F(7) + F(6) * C(5, 5) * F(6) + F(6) * C(5, 4) * F(2) + F(2) * C(4, 0) * F(7) + F(2) * C(4, 5) * F(6) + F(2) * C(4, 4) * F(2) + S(0);
        answer(7, 8) = F(7) * C(0, 0) * F(8) + F(7) * C(0, 5) * F(1) + F(7) * C(0, 4) * F(3) + F(6) * C(5, 0) * F(8) + F(6) * C(5, 5) * F(1) + F(6) * C(5, 4) * F(3) + F(2) * C(4, 0) * F(8) + F(2) * C(4, 5) * F(1) + F(2) * C(4, 4) * F(3) + 0.0;
        answer(8, 0) = F(8) * C(0, 0) * F(0) + F(8) * C(0, 5) * F(5) + F(8) * C(0, 4) * F(4) + F(1) * C(5, 0) * F(0) + F(1) * C(5, 5) * F(5) + F(1) * C(5, 4) * F(4) + F(3) * C(4, 0) * F(0) + F(3) * C(4, 5) * F(5) + F(3) * C(4, 4) * F(4) + 0.0;
        answer(8, 1) = F(8) * C(0, 5) * F(8) + F(8) * C(0, 1) * F(1) + F(8) * C(0, 3) * F(3) + F(1) * C(5, 5) * F(8) + F(1) * C(5, 1) * F(1) + F(1) * C(5, 3) * F(3) + F(3) * C(4, 5) * F(8) + F(3) * C(4, 1) * F(1) + F(3) * C(4, 3) * F(3) + S(5);
        answer(8, 2) = F(8) * C(0, 4) * F(7) + F(8) * C(0, 3) * F(6) + F(8) * C(0, 2) * F(2) + F(1) * C(5, 4) * F(7) + F(1) * C(5, 3) * F(6) + F(1) * C(5, 2) * F(2) + F(3) * C(4, 4) * F(7) + F(3) * C(4, 3) * F(6) + F(3) * C(4, 2) * F(2) + 0.0;
        answer(8, 3) = F(8) * C(0, 4) * F(8) + F(8) * C(0, 3) * F(1) + F(8) * C(0, 2) * F(3) + F(1) * C(5, 4) * F(8) + F(1) * C(5, 3) * F(1) + F(1) * C(5, 2) * F(3) + F(3) * C(4, 4) * F(8) + F(3) * C(4, 3) * F(1) + F(3) * C(4, 2) * F(3) + S(4);
        answer(8, 4) = F(8) * C(0, 4) * F(0) + F(8) * C(0, 3) * F(5) + F(8) * C(0, 2) * F(4) + F(1) * C(5, 4) * F(0) + F(1) * C(5, 3) * F(5) + F(1) * C(5, 2) * F(4) + F(3) * C(4, 4) * F(0) + F(3) * C(4, 3) * F(5) + F(3) * C(4, 2) * F(4) + 0.0;
        answer(8, 5) = F(8) * C(0, 5) * F(0) + F(8) * C(0, 1) * F(5) + F(8) * C(0, 3) * F(4) + F(1) * C(5, 5) * F(0) + F(1) * C(5, 1) * F(5) + F(1) * C(5, 3) * F(4) + F(3) * C(4, 5) * F(0) + F(3) * C(4, 1) * F(5) + F(3) * C(4, 3) * F(4) + 0.0;
        answer(8, 6) = F(8) * C(0, 5) * F(7) + F(8) * C(0, 1) * F(6) + F(8) * C(0, 3) * F(2) + F(1) * C(5, 5) * F(7) + F(1) * C(5, 1) * F(6) + F(1) * C(5, 3) * F(2) + F(3) * C(4, 5) * F(7) + F(3) * C(4, 1) * F(6) + F(3) * C(4, 3) * F(2) + 0.0;
        answer(8, 7) = F(8) * C(0, 0) * F(7) + F(8) * C(0, 5) * F(6) + F(8) * C(0, 4) * F(2) + F(1) * C(5, 0) * F(7) + F(1) * C(5, 5) * F(6) + F(1) * C(5, 4) * F(2) + F(3) * C(4, 0) * F(7) + F(3) * C(4, 5) * F(6) + F(3) * C(4, 4) * F(2) + 0.0;
        answer(8, 8) = F(8) * C(0, 0) * F(8) + F(8) * C(0, 5) * F(1) + F(8) * C(0, 4) * F(3) + F(1) * C(5, 0) * F(8) + F(1) * C(5, 5) * F(1) + F(1) * C(5, 4) * F(3) + F(3) * C(4, 0) * F(8) + F(3) * C(4, 5) * F(1) + F(3) * C(4, 4) * F(3) + S(0);
    } else if ( matMode == _PlaneStress ) {
        // Save terms associated with H = [du/dx dv/dy du/dy dv/dx]

        answer.resize(4, 4);
        answer(0, 0) = F(0) * C(0, 0) * F(0) + F(0) * C(0, 2) * F(2) + F(2) * C(2, 0) * F(0) + F(2) * C(2, 2) * F(2) + S(0);
        answer(0, 1) = F(0) * C(0, 2) * F(3) + F(0) * C(0, 1) * F(1) + F(2) * C(2, 2) * F(3) + F(2) * C(2, 1) * F(1) + 0.0;
        answer(0, 2) = F(0) * C(0, 2) * F(0) + F(0) * C(0, 1) * F(2) + F(2) * C(2, 2) * F(0) + F(2) * C(2, 1) * F(2) + S(2);
        answer(0, 3) = F(0) * C(0, 0) * F(3) + F(0) * C(0, 2) * F(1) + F(2) * C(2, 0) * F(3) + F(2) * C(2, 2) * F(1) + 0.0;
        answer(1, 0) = F(3) * C(2, 0) * F(0) + F(3) * C(2, 2) * F(2) + F(1) * C(1, 0) * F(0) + F(1) * C(1, 2) * F(2) + 0.0;
        answer(1, 1) = F(3) * C(2, 2) * F(3) + F(3) * C(2, 1) * F(1) + F(1) * C(1, 2) * F(3) + F(1) * C(1, 1) * F(1) + S(1);
        answer(1, 2) = F(3) * C(2, 2) * F(0) + F(3) * C(2, 1) * F(2) + F(1) * C(1, 2) * F(0) + F(1) * C(1, 1) * F(2) + 0.0;
        answer(1, 3) = F(3) * C(2, 0) * F(3) + F(3) * C(2, 2) * F(1) + F(1) * C(1, 0) * F(3) + F(1) * C(1, 2) * F(1) + S(2);
        answer(2, 0) = F(0) * C(2, 0) * F(0) + F(0) * C(2, 2) * F(2) + F(2) * C(1, 0) * F(0) + F(2) * C(1, 2) * F(2) + S(2);
        answer(2, 1) = F(0) * C(2, 2) * F(3) + F(0) * C(2, 1) * F(1) + F(2) * C(1, 2) * F(3) + F(2) * C(1, 1) * F(1) + 0.0;
        answer(2, 2) = F(0) * C(2, 2) * F(0) + F(0) * C(2, 1) * F(2) + F(2) * C(1, 2) * F(0) + F(2) * C(1, 1) * F(2) + S(1);
        answer(2, 3) = F(0) * C(2, 0) * F(3) + F(0) * C(2, 2) * F(1) + F(2) * C(1, 0) * F(3) + F(2) * C(1, 2) * F(1) + 0.0;
        answer(3, 0) = F(3) * C(0, 0) * F(0) + F(3) * C(0, 2) * F(2) + F(1) * C(2, 0) * F(0) + F(1) * C(2, 2) * F(2) + 0.0;
        answer(3, 1) = F(3) * C(0, 2) * F(3) + F(3) * C(0, 1) * F(1) + F(1) * C(2, 2) * F(3) + F(1) * C(2, 1) * F(1) + S(2);
        answer(3, 2) = F(3) * C(0, 2) * F(0) + F(3) * C(0, 1) * F(2) + F(1) * C(2, 2) * F(0) + F(1) * C(2, 1) * F(2) + 0.0;
        answer(3, 3) = F(3) * C(0, 0) * F(3) + F(3) * C(0, 2) * F(1) + F(1) * C(2, 0) * F(3) + F(1) * C(2, 2) * F(1) + S(0);
    } else if ( matMode == _PlaneStrain ) {
        //Save terms associated with H = [du/dx, dv/dy, dw/dz, du/dy, dv/dx] //@todo not fully checked

        answer.resize(5, 5);
        answer(0, 0) = F(0) * C(0, 0) * F(0) + F(0) * C(0, 3) * F(3) + F(3) * C(3, 0) * F(0) + F(3) * C(3, 3) * F(3) + S(0);
        answer(0, 1) = F(0) * C(0, 3) * F(4) + F(0) * C(0, 1) * F(1) + F(3) * C(3, 3) * F(4) + F(3) * C(3, 1) * F(1) + 0.0;
        answer(0, 2) = F(0) * C(0, 2) * F(2) + F(3) * C(3, 2) * F(2) + 0.0;
        answer(0, 3) = F(0) * C(0, 3) * F(0) + F(0) * C(0, 1) * F(3) + F(3) * C(3, 3) * F(0) + F(3) * C(3, 1) * F(3) + S(3);
        answer(0, 4) = F(0) * C(0, 0) * F(4) + F(0) * C(0, 3) * F(1) + F(3) * C(3, 0) * F(4) + F(3) * C(3, 3) * F(1) + 0.0;
        answer(1, 0) = F(4) * C(3, 0) * F(0) + F(4) * C(3, 3) * F(3) + F(1) * C(1, 0) * F(0) + F(1) * C(1, 3) * F(3) + 0.0;
        answer(1, 1) = F(4) * C(3, 3) * F(4) + F(4) * C(3, 1) * F(1) + F(1) * C(1, 3) * F(4) + F(1) * C(1, 1) * F(1) + S(1);
        answer(1, 2) = F(4) * C(3, 2) * F(2) + F(1) * C(1, 2) * F(2) + 0.0;
        answer(1, 3) = F(4) * C(3, 3) * F(0) + F(4) * C(3, 1) * F(3) + F(1) * C(1, 3) * F(0) + F(1) * C(1, 1) * F(3) + 0.0;
        answer(1, 4) = F(4) * C(3, 0) * F(4) + F(4) * C(3, 3) * F(1) + F(1) * C(1, 0) * F(4) + F(1) * C(1, 3) * F(1) + S(3);
        answer(2, 0) = F(2) * C(2, 0) * F(0) + F(2) * C(2, 3) * F(3) + 0.0;
        answer(2, 1) = F(2) * C(2, 3) * F(4) + F(2) * C(2, 1) * F(1) + 0.0;
        answer(2, 2) = F(2) * C(2, 2) * F(2) + S(2);
        answer(2, 3) = F(2) * C(2, 3) * F(0) + F(2) * C(2, 1) * F(3) + 0.0;
        answer(2, 4) = F(2) * C(2, 0) * F(4) + F(2) * C(2, 3) * F(1) + 0.0;
        answer(3, 0) = F(0) * C(3, 0) * F(0) + F(0) * C(3, 3) * F(3) + F(3) * C(1, 0) * F(0) + F(3) * C(1, 3) * F(3) + S(3);
        answer(3, 1) = F(0) * C(3, 3) * F(4) + F(0) * C(3, 1) * F(1) + F(3) * C(1, 3) * F(4) + F(3) * C(1, 1) * F(1) + 0.0;
        answer(3, 2) = F(0) * C(3, 2) * F(2) + F(3) * C(1, 2) * F(2) + 0.0;
        answer(3, 3) = F(0) * C(3, 3) * F(0) + F(0) * C(3, 1) * F(3) + F(3) * C(1, 3) * F(0) + F(3) * C(1, 1) * F(3) + S(1);
        answer(3, 4) = F(0) * C(3, 0) * F(4) + F(0) * C(3, 3) * F(1) + F(3) * C(1, 0) * F(4) + F(3) * C(1, 3) * F(1) + 0.0;
        answer(4, 0) = F(4) * C(0, 0) * F(0) + F(4) * C(0, 3) * F(3) + F(1) * C(3, 0) * F(0) + F(1) * C(3, 3) * F(3) + 0.0;
        answer(4, 1) = F(4) * C(0, 3) * F(4) + F(4) * C(0, 1) * F(1) + F(1) * C(3, 3) * F(4) + F(1) * C(3, 1) * F(1) + S(3);
        answer(4, 2) = F(4) * C(0, 2) * F(2) + F(1) * C(3, 2) * F(2) + 0.0;
        answer(4, 3) = F(4) * C(0, 3) * F(0) + F(4) * C(0, 1) * F(3) + F(1) * C(3, 3) * F(0) + F(1) * C(3, 1) * F(3) + 0.0;
        answer(4, 4) = F(4) * C(0, 0) * F(4) + F(4) * C(0, 3) * F(1) + F(1) * C(3, 0) * F(4) + F(1) * C(3, 3) * F(1) + S(0);
    } else if ( matMode == _1dMat ) {
        //Save terms associated with H = [du/dx]
        /// @todo is this really correct??

        answer.resize(1, 1);
        answer(0, 0) = F(0) * C(0, 0) * F(0) + S(0);
    }
}

void
StructuralMaterial :: give_dPdF_from(const FloatMatrix &dSdE, FloatMatrix &answer, GaussPoint *gp)
{
    // Default implementation for converting dSdE to dPdF. This includes updating the
    // state variables of P and F.
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    FloatArray reducedvF, reducedvP, reducedvS;
    const FloatArray &vF = status->giveTempFVector();
    const FloatArray &vP = status->giveTempPVector();
    const FloatArray &vS = status->giveTempStressVector();

    MaterialMode matMode = gp->giveMaterialMode();
    ///@todo This is unnecessary and inefficient, just a hack for now (it will be changed when everything is in working order!)
    StructuralMaterial :: giveReducedVectorForm(reducedvF, vF, matMode);
    StructuralMaterial :: giveReducedVectorForm(reducedvP, vP, matMode);
    StructuralMaterial :: giveReducedSymVectorForm(reducedvS, vS, matMode);
    //this->convert_P_2_S(reducedvS, reducedvP, reducedvF, matMode);
    this->convert_dSdE_2_dPdF(answer, dSdE, reducedvS, reducedvF, matMode);
}


void
StructuralMaterial :: giveStiffnessMatrix(FloatMatrix &answer,
                                          MatResponseMode rMode,
                                          GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMat:
        this->give3dMaterialStiffnessMatrix(answer, rMode, gp, tStep);
        break;
    case _PlaneStress:
        this->givePlaneStressStiffMtrx(answer, rMode, gp, tStep);
        break;
    case _PlaneStrain:
        this->givePlaneStrainStiffMtrx(answer, rMode, gp, tStep);
        break;
    case _1dMat:
        this->give1dStressStiffMtrx(answer, rMode, gp, tStep);
        break;

    case _PlateLayer:
        this->givePlateLayerStiffMtrx(answer, rMode, gp, tStep);
        break;
    case _2dBeamLayer:
        this->give2dBeamLayerStiffMtrx(answer, rMode, gp, tStep);
        break;
    case _Fiber:
        this->giveFiberStiffMtrx(answer, rMode, gp, tStep);
        break;
    default:
        OOFEM_ERROR2( "StructuralMaterial :: giveStiffnessMatrix : unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}


void
StructuralMaterial :: give3dMaterialStiffnessMatrix_dPdF(FloatMatrix &answer,
                                                         MatResponseMode mode,
                                                         GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix dSdE;
    this->give3dMaterialStiffnessMatrix(dSdE, mode, gp, tStep);
    this->give_dPdF_from(dSdE, answer, gp);
}


void
StructuralMaterial :: givePlaneStressStiffMtrx_dPdF(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix dSdE;
    this->givePlaneStressStiffMtrx(dSdE, mode, gp, tStep);
    this->give_dPdF_from(dSdE, answer, gp);
}


void
StructuralMaterial :: givePlaneStrainStiffMtrx_dPdF(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix dSdE;
    this->givePlaneStrainStiffMtrx(dSdE, mode, gp, tStep);
    this->give_dPdF_from(dSdE, answer, gp);
}


void
StructuralMaterial :: give1dStressStiffMtrx_dPdF(FloatMatrix &answer,
                                                 MatResponseMode mode,
                                                 GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix dSdE;
    this->give1dStressStiffMtrx(dSdE, mode, gp, tStep);
    this->give_dPdF_from(dSdE, answer, gp);
}


void
StructuralMaterial :: give3dMaterialStiffnessMatrix_dCde(FloatMatrix &answer,
                                                         MatResponseMode mode,
                                                         GaussPoint *gp, TimeStep *tStep)
{
    ///@todo what should be default implementaiton?
    OOFEM_ERROR("StructuralMaterial :: give3dMaterialStiffnessMatrix_dCde - There is no default implementation");
}


void
StructuralMaterial :: givePlaneStressStiffMtrx_dCde(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("StructuralMaterial :: givePlaneStressStiffMtrx_dCde - There is no default implementation");
}


void
StructuralMaterial :: givePlaneStrainStiffMtrx_dCde(FloatMatrix &answer,
                                                    MatResponseMode mode,
                                                    GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("StructuralMaterial :: givePlaneStrainStiffMtrx_dCde - There is no default implementation");
}


void
StructuralMaterial :: give1dStressStiffMtrx_dCde(FloatMatrix &answer,
                                                 MatResponseMode mode,
                                                 GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("StructuralMaterial :: give1dStressStiffMtrx_dCde - There is no default implementation");
}


void
StructuralMaterial :: convert_P_2_S(FloatArray &answer, const FloatArray &reducedvP, const FloatArray &reducedvF, MaterialMode matMode)
{
    // Converts first Piola-Kirchhoff stress to second Piola-Kirchhoff stress: S = inv(F)*P
    // Output size will be according to MaterialMode

    FloatArray vF, vP;
    StructuralMaterial :: giveFullVectorFormF(vF, reducedvF, matMode); // 9 components
    StructuralMaterial :: giveFullVectorForm(vP, reducedvP, matMode);
    FloatMatrix F, P, S, invF;
    F.beMatrixForm(vF);
    P.beMatrixForm(vP);
    invF.beInverseOf(F);
    S.beProductOf(invF, P);
    FloatArray vS;
    vS.beSymVectorForm(S); // 6 components

    StructuralMaterial :: giveReducedSymVectorForm(answer, vS, matMode); // convert back to reduced size
}


void
StructuralMaterial :: convert_S_2_P(FloatArray &answer, const FloatArray &reducedvS, const FloatArray &reducedvF, MaterialMode matMode)
{
    // Converts second Piola-Kirchhoff stress to first Piola-Kirchhoff stress: P = F*S
    // Output size will be according to MaterialMode

    FloatArray vF, vS, vP;
    StructuralMaterial :: giveFullVectorFormF(vF, reducedvF, matMode);   // 9 components
    StructuralMaterial :: giveFullSymVectorForm(vS, reducedvS, matMode); // 6 components
    FloatMatrix F, P, S;
    F.beMatrixForm(vF);
    S.beMatrixForm(vS);
    P.beProductOf(F, S);
    vP.beVectorForm(P);
    StructuralMaterial :: giveReducedVectorForm(answer, vP, matMode);   // convert back to reduced size
}


void
StructuralMaterial :: giveStressDependentPartOfStrainVector(FloatArray &answer, GaussPoint *gp,
                                                            const FloatArray &reducedStrainVector,
                                                            TimeStep *tStep, ValueModeType mode)
{
    /*
     * This functions subtract from reducedStrainVector its stress independent part
     * caused by temperature, shrinkage and possibly by other phenomena.
     */
    FloatArray epsilonTemperature;

    answer = reducedStrainVector;
    this->computeStressIndependentStrainVector(epsilonTemperature, gp, tStep, mode);
    if ( epsilonTemperature.giveSize() ) {
        answer.subtract(epsilonTemperature);
    }
}


int
StructuralMaterial :: giveSizeOfVoigtSymVector(MaterialMode mode)
{
    IntArray indx;
    StructuralMaterial :: giveVoigtSymVectorMask(indx, mode);
    return indx.giveSize();
}


int
StructuralMaterial :: giveSizeOfVoigtVector(MaterialMode mode)
{
    IntArray indx;
    StructuralMaterial :: giveVoigtVectorMask(indx, mode);
    return indx.giveSize();
}


void
StructuralMaterial :: giveInvertedVoigtVectorMask(IntArray &answer, MaterialMode mmode)
{
    IntArray mask;
    answer.resize( StructuralMaterial :: giveVoigtSymVectorMask(mask, mmode) );
    answer.zero();
    for ( int i = 1; i <= mask.giveSize(); i++ ) {
        answer.at( mask.at(i) ) = i;
    }
}


int
StructuralMaterial :: giveVoigtSymVectorMask(IntArray &answer, MaterialMode mmode)
{
    // The same as giveVoigtVectorMask but returns a mask corresponding to a symmetric
    // second order tensor.
    //
    // Returns a mask of the vector indices corresponding to components in a symmetric
    // second order tensor of some stress/strain/deformation measure that performs work.
    // Thus, components corresponding to imposed zero stress (e.g. plane stress etc.) are
    // not included. On the other hand, if zero strain components are imposed( e.g. plane
    // strain etc.) this condition must be taken into account in geometrical relations.
    // Therefore, these corresponding components are included in the reduced vector.
    // Which components to include are given by the particular MaterialMode.

    switch ( mmode ) {
    case _3dMat:
    case _3dMicroplane:
        answer.resize(6);
        for ( int i = 1; i <= 6; i++ ) {
            answer.at(i) = i;
        }

        return 6;

    case _PlaneStress:
        answer.resize(3);
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 6;
        return 6;

    case _PlaneStrain:
        answer.resize(4);
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 3;
        answer.at(4) = 6;
        return 6;

    case _1dMat:
        answer.resize(1);
        answer.at(1) = 1;
        return 6;

    case _PlateLayer:
        answer.resize(5);
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 4;
        answer.at(4) = 5;
        answer.at(5) = 6;
        return 6;

    case _2dBeamLayer:
        answer.resize(2);
        answer.at(1) = 1;
        answer.at(2) = 5;
        return 6;

    case _Fiber:
        answer.resize(3);
        answer.at(1) = 1;
        answer.at(2) = 5;
        answer.at(3) = 6;
        return 6;

    case _2dPlate:
        answer.resize(5);
        answer.at(1) = 4;
        answer.at(2) = 5;
        answer.at(3) = 6;
        answer.at(4) = 7;
        answer.at(5) = 8;
        return 8;

    case _2dBeam:
        answer.resize(3);
        answer.at(1) = 1;
        answer.at(2) = 4;
        answer.at(3) = 7;
        return 8;

    case _3dBeam: ///@todo This isn't actually fixed yet. Should be made compatible with 3dShell and 2dBeam
        answer.resize(6);
#if 1
        for ( int i = 1; i <= 6; i++ ) {
            answer.at(i) = i;
        }

        return 6;

#else
        answer.at(1) = 1;
        answer.at(2) = 5;
        answer.at(3) = 6;
        answer.at(4) = 7;
        answer.at(5) = 8;
        answer.at(6) = 9;
        return 12;

#endif
    case _3dShell:
        answer.resize(8);
        for ( int i = 1; i <= 8; i++ ) {
            answer.at(i) = i;
        }

        return 8;

    case _PlaneStressRot:
        answer.resize(4);
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 6;
        answer.at(4) = 7;
        return 7;

    case _1dInterface:
        answer.setValues(1, 1);
        return 1;

    case _2dInterface:
        answer.setValues(2, 1, 2);
        return 2;

    case _3dInterface:
        answer.setValues(3, 1, 2, 3);
        return 3;

    case _2dLattice:
        answer.setValues(3, 1, 2, 3);
        return 3;

    case _Unknown:
        answer.clear();
        return 0;

    default:
        OOFEM_ERROR2( "StructuralMaterial :: giveVoigtSymVectorMask : unknown mode (%s)", __MaterialModeToString(mmode) );
        return 0;
    }
}


int
StructuralMaterial :: giveVoigtVectorMask(IntArray &answer, MaterialMode mmode)
{
    // Returns a mask of the vector indices corresponding to components in a general
    // (non-symmetric) second order tensor of some stress/strain/deformation measure that
    // performs work. Thus, components corresponding to imposed zero stress (e.g. plane
    // stress etc.) are not included. On the other hand, if zero strain components are
    // imposed( e.g. plane strain etc.) this condition must be taken into account in
    // geometrical relations. Therefore, these corresponding components are included in
    // the reduced vector. Which components to include are given by the particular MaterialMode.
    //
    /// @todo add additional modes if they relevant.

    switch ( mmode ) {
    case _3dMat:
        answer.resize(9);
        for ( int i = 1; i <= 9; i++ ) {
            answer.at(i) = i;
        }

        return 9;

    case _PlaneStress:
        answer.resize(4);
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 6;
        answer.at(4) = 9;
        return 9;

    case _PlaneStrain:
        answer.resize(5);
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 3;
        answer.at(4) = 6;
        answer.at(5) = 9;
        return 9;

    case _1dMat:
        answer.resize(1);
        answer.at(1) = 1;
        return 9;

    default:
        //OOFEM_ERROR2( "StructuralMaterial :: giveVoigtVectorMask: unknown mode (%s)", __MaterialModeToString(mmode) );
        return 0;
    }
}


void
StructuralMaterial :: givePlaneStressStiffMtrx(FloatMatrix &answer,
                                               MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep)
//
// returns Mat stiffness for PlaneStress
//
{
    FloatMatrix m3d, invMat3d, invAnswer;

    this->give3dMaterialStiffnessMatrix(m3d, mode, gp, tStep);

    invMat3d.beInverseOf(m3d);

    invAnswer.resize(3, 3);
    //invAnswer.beSubMatrixOf(invMat3d, indx, indx);

    invAnswer.at(1, 1) = invMat3d.at(1, 1);
    invAnswer.at(1, 2) = invMat3d.at(1, 2);
    invAnswer.at(1, 3) = invMat3d.at(1, 6);

    invAnswer.at(2, 1) = invMat3d.at(2, 1);
    invAnswer.at(2, 2) = invMat3d.at(2, 2);
    invAnswer.at(2, 3) = invMat3d.at(2, 6);

    invAnswer.at(3, 1) = invMat3d.at(6, 1);
    invAnswer.at(3, 2) = invMat3d.at(6, 2);
    invAnswer.at(3, 3) = invMat3d.at(6, 6);

    answer.beInverseOf(invAnswer);
}

void
StructuralMaterial :: givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                               MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep)
//
// return material stiffness matrix for PlaneStrain mode
//
{
    FloatMatrix m3d;

    this->give3dMaterialStiffnessMatrix(m3d, mode, gp, tStep);

    answer.resize(4, 4);
    answer.zero();
    //answer.beSubMatrixOf(m3d, indx, indx);

    answer.at(1, 1) = m3d.at(1, 1);
    answer.at(1, 2) = m3d.at(1, 2);
    answer.at(1, 4) = m3d.at(1, 6);

    answer.at(2, 1) = m3d.at(2, 1);
    answer.at(2, 2) = m3d.at(2, 2);
    answer.at(2, 4) = m3d.at(2, 6);

    answer.at(3, 1) = m3d.at(3, 1);
    answer.at(3, 2) = m3d.at(3, 2);
    answer.at(3, 4) = m3d.at(3, 6);

    answer.at(4, 1) = m3d.at(6, 1);
    answer.at(4, 2) = m3d.at(6, 2);
    answer.at(4, 4) = m3d.at(6, 6);
}

void
StructuralMaterial :: give1dStressStiffMtrx(FloatMatrix &answer,
                                            MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *tStep)
//
// return material stiffness matrix for 1d stress strain mode
//
{
    FloatMatrix m3d, invMat3d;
    double val11;

    this->give3dMaterialStiffnessMatrix(m3d, mode, gp, tStep);

    invMat3d.beInverseOf(m3d);
    val11 = invMat3d.at(1, 1);
    answer.resize(1, 1);
    answer.at(1, 1) = 1. / val11;
}


void
StructuralMaterial :: give2dBeamLayerStiffMtrx(FloatMatrix &answer,
                                               MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep)
//
// return material stiffness matrix for2dBeamLayer mode
//
{
    FloatMatrix m3d, invMat3d, invMatLayer(2, 2);

    this->give3dMaterialStiffnessMatrix(m3d, mode, gp, tStep);

    invMat3d.beInverseOf(m3d);

    invMatLayer.at(1, 1) = invMat3d.at(1, 1);
    invMatLayer.at(1, 2) = invMat3d.at(1, 5);
    invMatLayer.at(2, 1) = invMat3d.at(5, 1);
    invMatLayer.at(2, 2) = invMat3d.at(5, 5);

    answer.beInverseOf(invMatLayer);
}


void
StructuralMaterial :: givePlateLayerStiffMtrx(FloatMatrix &answer,
                                              MatResponseMode mode,
                                              GaussPoint *gp,
                                              TimeStep *tStep)
//
// return material stiffness matrix for 2dPlateLayer
//
{
    FloatMatrix m3d, invMat3d, invMatLayer(5, 5);

    this->give3dMaterialStiffnessMatrix(m3d, mode, gp, tStep);

    invMat3d.beInverseOf(m3d);
    //invMatLayer.beSubMatrixOf(invMat3d, indx, indx);

    for ( int i = 1; i <= 2; i++ ) {
        for ( int j = 1; j <= 2; j++ ) {
            invMatLayer.at(i, j) = invMat3d.at(i, j);
        }
    }

    for ( int i = 4; i <= 6; i++ ) {
        for ( int j = 4; j <= 6; j++ ) {
            invMatLayer.at(i - 1, j - 1) = invMat3d.at(i, j);
        }
    }

    for ( int i = 1; i <= 2; i++ ) {
        for ( int j = 4; j <= 6; j++ ) {
            invMatLayer.at(i, j - 1) = invMat3d.at(i, j);
            invMatLayer.at(j - 1, i) = invMat3d.at(j, i);
        }
    }

    answer.beInverseOf(invMatLayer);
}

void
StructuralMaterial :: giveFiberStiffMtrx(FloatMatrix &answer,
                                         MatResponseMode mode,
                                         GaussPoint *gp,
                                         TimeStep *tStep)
//
// return material stiffness matrix for 2dPlateLayer
//
{
    FloatMatrix m3d, invMat3d, invMatLayer(3, 3);

    this->give3dMaterialStiffnessMatrix(m3d, mode, gp, tStep);

    invMat3d.beInverseOf(m3d);

    invMatLayer.at(1, 1) = invMat3d.at(1, 1);
    invMatLayer.at(1, 2) = invMat3d.at(1, 5);
    invMatLayer.at(1, 3) = invMat3d.at(1, 6);
    invMatLayer.at(2, 1) = invMat3d.at(5, 1);
    invMatLayer.at(2, 2) = invMat3d.at(5, 5);
    invMatLayer.at(2, 3) = invMat3d.at(5, 6);
    invMatLayer.at(3, 1) = invMat3d.at(6, 1);
    invMatLayer.at(3, 2) = invMat3d.at(6, 5);
    invMatLayer.at(3, 3) = invMat3d.at(6, 6);

    answer.beInverseOf(invMatLayer);
}


void
StructuralMaterial :: computePrincipalValues(FloatArray &answer, const FloatArray &s, stressStrainPrincMode mode)
//
// This function computes the principal values of strains or stresses.
// Strains/stresses are stored in vector form in array s.
// Engineering notation is used.
//
// Problem size (3D/2D) is recognized automatically according to the
// vector size.
// If size = 6 -> 3D problem, then array s contains:
//                            {Sxx,Syy,Szz,Syz,Szx,Sxy} if mode = principal_stress
//                            {Exx,Eyy,Ezz,GMyz,GMzx,GMxy} if mode = principal_strain
// if size = 3 -> 2D problem, then array s contains:
//                            {Sxx,Syy,Sxy} if mode = principal_stress
//                            {Exx,Eyy,GMxy} if mode = principal_strain
//
// if size = 4 -> 2D problem (with normal out-of-plane component), then array s contains:
//                            {Sxx,Syy,Szz,Sxy} if mode = principal_stress
//                            {Exx,Eyy,Ezz,GMxy} if mode = principal_strain
//
// Return Values:
//
//    array answer -> principal strains or stresses
//
{
    int size = s.giveSize();
    if ( !( ( size == 3 ) || ( size == 4 ) || ( size == 6 ) ) ) {
        OOFEM_ERROR("StructuralMaterial :: ComputePrincipalValues: Vector size mismatch");
    }

    double swap;
    int nonzeroFlag = 0;
    if ( ( size == 3 ) || ( size == 4 ) ) {
        // 2D problem
        double ast, dst, D = 0.0;
        answer.resize(size - 1);

        for ( int i = 1; i <= size; i++ ) {
            if ( fabs( s.at(i) ) > 1.e-20 ) {
                nonzeroFlag = 1;
            }
        }

        if ( nonzeroFlag == 0 ) {
            answer.zero();
            return;
        }

        ast = s.at(1) + s.at(2);
        dst = s.at(1) - s.at(2);
        if ( mode == principal_strain ) {
            D = dst * dst + s.at(size) * s.at(size);
        } else if ( mode == principal_stress ) {
            D = dst * dst + 4.0 * s.at(size) * s.at(size);
        } else {
            OOFEM_ERROR("StructuralMaterial :: ComputePrincipalValues: not supported");
        }

        if ( D < 0. ) {
            OOFEM_ERROR("StructuralMaterial :: ComputePrincipalValues: Imaginary roots ");
        }

        D = sqrt(D);
        answer.at(1) = 0.5 * ( ast - D );
        answer.at(2) = 0.5 * ( ast + D );
        if ( size == 4 ) {
            answer.at(3) = s.at(3);
        }
    } else {
        // 3D problem
        double I1 = 0.0, I2 = 0.0, I3 = 0.0, help, s1, s2, s3;

        for ( int i = 1; i <= size; i++ ) {
            if ( fabs( s.at(i) ) > 1.e-20 ) {
                nonzeroFlag = 1;
            }
        }

        answer.resize(3);
        answer.zero();
        if ( nonzeroFlag == 0 ) {
            return;
        }



        if ( mode == principal_stress ) {
            I1 = s.at(1) + s.at(2) + s.at(3);
            I2 = s.at(1) * s.at(2) + s.at(2) * s.at(3) + s.at(3) * s.at(1) -
            ( s.at(4) * s.at(4) + s.at(5) * s.at(5) + s.at(6) * s.at(6) );
            I3 = s.at(1) * s.at(2) * s.at(3) + 2. * s.at(4) * s.at(5) * s.at(6) -
            ( s.at(1) * s.at(4) * s.at(4) + s.at(2) * s.at(5) * s.at(5) +
             s.at(3) * s.at(6) * s.at(6) );
        } else if ( mode == principal_deviatoricstress ) {
            help = ( s.at(1) + s.at(2) + s.at(3) ) / 3.0;
            I1 = 0.;
            I2 = -( 1. / 6. ) * ( ( s.at(1) - s.at(2) ) * ( s.at(1) - s.at(2) ) + ( s.at(2) - s.at(3) ) * ( s.at(2) - s.at(3) ) +
                                 ( s.at(3) - s.at(1) ) * ( s.at(3) - s.at(1) ) ) - s.at(4) * s.at(4) - s.at(5) * s.at(5) -
            s.at(6) * s.at(6);
            I3 = ( s.at(1) - help ) * ( s.at(2) - help ) * ( s.at(3) - help ) + 2. * s.at(4) * s.at(5) * s.at(6) -
            s.at(5) * s.at(5) * ( s.at(2) - help ) - s.at(4) * s.at(4) * ( s.at(1) - help ) -
            s.at(6) * s.at(6) * ( s.at(3) - help );
        } else if ( mode == principal_strain ) {
            I1 = s.at(1) + s.at(2) + s.at(3);
            I2 = s.at(1) * s.at(2) + s.at(2) * s.at(3) + s.at(3) * s.at(1) -
            0.25 * ( s.at(4) * s.at(4) + s.at(5) * s.at(5) + s.at(6) * s.at(6) );
            I3 = s.at(1) * s.at(2) * s.at(3) +
            0.25 * ( s.at(4) * s.at(5) * s.at(6) - s.at(1) * s.at(4) * s.at(4) -
                    s.at(2) * s.at(5) * s.at(5) - s.at(3) * s.at(6) * s.at(6) );
        } else {
            OOFEM_ERROR("StructuralMaterial :: ComputePrincipalValues: not supported");
        }

        /*
         * Call cubic3r to ensure that all three real eigenvalues will be found, because we have symmetric tensor.
         * This allows to overcome various rounding errors when solving general cubic equation.
         */
        int n;
        cubic3r( ( double ) -1., I1, -I2, I3, & s1, & s2, & s3, & n );

        if ( n > 0 ) {
            answer.at(1) = s1;
        }

        if ( n > 1 ) {
            answer.at(2) = s2;
        }

        if ( n > 2 ) {
            answer.at(3) = s3;
        }
    }

    //sort the results
    for ( int i = 1; i < answer.giveSize(); i++ ) {
        for ( int j = 1; j < answer.giveSize(); j++ ) {
            if ( answer.at(j + 1) > answer.at(j) ) {
                swap = answer.at(j + 1);
                answer.at(j + 1) = answer.at(j);
                answer.at(j) = swap;
            }
        }
    }
}

void
StructuralMaterial :: computePrincipalValDir(FloatArray &answer, FloatMatrix &dir,
                                             const FloatArray &s,
                                             stressStrainPrincMode mode)
//
// This function computes the principal values & directions corresponding to principal values
// of strains or streses.
// strains/streses are stored in vector form in array s.
// Engineering notation is used.
//
// Problem size (3D/2D) is recognized automatically according to
// vector size.
// If size = 6 -> 3D problem, then array s contains:
//                            {Sxx,Syy,Szz,Syz,Szx,Sxy} if mode = principal_stress
//                            {Exx,Eyy,Ezz,GMyz,GMzx,GMxy} if mode = principal_strain
// if size = 3 -> 2D problem, then array s contains:
//                            {Sxx,Syy,Sxy} if mode = principal_stress
//                            {Exx,Eyy,GMxy} if mode = principal_strain
//
// mode      - principal strains
//           - principal stress
//
// Input Values:
// mode
// s
//
// Return Values:
//
// matrix dir -> principal directions of strains or stresses
// array answer -> principal strains or stresses
//
{
    FloatMatrix ss;
    FloatArray sp;
    double swap;
    int nval, size = s.giveSize();
    int nonzeroFlag = 0;

    // printf ("size is %d\n",size);
    if ( !( ( size == 3 ) || ( size == 4 ) || ( size == 6 ) ) ) {
        OOFEM_ERROR("StructuralMaterial :: computePrincipalValDir: Vector size mismatch");
    }

    if ( ( size == 3 ) || ( size == 4 ) ) {
        // 2D problem
        ss.resize(2, 2);
        answer.resize(2);

        for ( int i = 1; i <= size; i++ ) {
            if ( fabs( s.at(i) ) > 1.e-20 ) {
                nonzeroFlag = 1;
            }
        }

        if ( nonzeroFlag == 0 ) {
            answer.zero();
            ss.zero();
            return;
        }

        ss.at(1, 1) = s.at(1);
        ss.at(2, 2) = s.at(2);

        if ( mode == principal_strain ) {
            ss.at(1, 2) = ss.at(2, 1) = 0.5 * s.at(size);
        } else if ( mode == principal_stress ) {
            ss.at(1, 2) = ss.at(2, 1) = s.at(size);
        } else {
            OOFEM_ERROR("StructuralMaterial :: computePrincipalValDir: not supported");
        }
    } else {
        // 3D problem
        double help;
        ss.resize(3, 3);
        answer.resize(3);

        for ( int i = 1; i <= size; i++ ) {
            if ( fabs( s.at(i) ) > 1.e-20 ) {
                nonzeroFlag = 1;
            }
        }

        if ( nonzeroFlag == 0 ) {
            answer.zero();
            ss.zero();
            return;
        }

        if ( mode == principal_stress ) {
            ss.at(1, 1) = s.at(1);
            ss.at(2, 2) = s.at(2);
            ss.at(3, 3) = s.at(3);
            ss.at(1, 2) = ss.at(2, 1) = s.at(6);
            ss.at(1, 3) = ss.at(3, 1) = s.at(5);
            ss.at(2, 3) = ss.at(3, 2) = s.at(4);
        } else if ( mode == principal_deviatoricstress ) {
            help = ( s.at(1) + s.at(2) + s.at(3) ) / 3.0;
            ss.at(1, 1) = s.at(1) - help;
            ss.at(2, 2) = s.at(2) - help;
            ss.at(3, 3) = s.at(3) - help;
            ss.at(1, 2) = ss.at(2, 1) = s.at(6);
            ss.at(1, 3) = ss.at(3, 1) = s.at(5);
            ss.at(2, 3) = ss.at(3, 2) = s.at(4);
        } else if ( mode == principal_strain ) {
            ss.at(1, 1) = s.at(1);
            ss.at(2, 2) = s.at(2);
            ss.at(3, 3) = s.at(3);
            ss.at(1, 2) = ss.at(2, 1) = 0.5 * s.at(6);
            ss.at(1, 3) = ss.at(3, 1) = 0.5 * s.at(5);
            ss.at(2, 3) = ss.at(3, 2) = 0.5 * s.at(4);
        } else {
            OOFEM_ERROR("StructuralMaterial :: computePrincipalDirection: not supported");
        }
    }

#if 0
    ss.Jacobi(& answer, & dir, & i);
#else
    ss.jaco_(answer, dir, 10);
#endif
    // sort results
    nval = 2;
    if ( size == 6 ) {
        nval = 3;
    }

    for ( int ii = 1; ii < nval; ii++ ) {
        for ( int jj = 1; jj < nval; jj++ ) {
            if ( answer.at(jj + 1) > answer.at(jj) ) {
                // swap eigen values and eigen vectors
                swap = answer.at(jj + 1);
                answer.at(jj + 1) = answer.at(jj);
                answer.at(jj) = swap;
                for ( int kk = 1; kk <= nval; kk++ ) {
                    swap = dir.at(kk, jj + 1);
                    dir.at(kk, jj + 1) = dir.at(kk, jj);
                    dir.at(kk, jj) = swap;
                }
            }
        }
    }
}


double
StructuralMaterial :: computeVonMisesStress(const FloatArray *currentStress)
{
    double J2;
    double v1, v2, v3;

    if ( currentStress == NULL || currentStress->giveSize() != 6 ) {
        return 0.0;
    }

    v1 = ( ( currentStress->at(1) - currentStress->at(2) ) * ( currentStress->at(1) - currentStress->at(2) ) );
    v2 = ( ( currentStress->at(2) - currentStress->at(3) ) * ( currentStress->at(2) - currentStress->at(3) ) );
    v3 = ( ( currentStress->at(3) - currentStress->at(1) ) * ( currentStress->at(3) - currentStress->at(1) ) );

    J2 = ( 1. / 6. ) * ( v1 + v2 + v3 ) + currentStress->at(4) * currentStress->at(4) +
    currentStress->at(5) * currentStress->at(5) + currentStress->at(6) * currentStress->at(6);

    return sqrt(3 * J2);
}


void
StructuralMaterial :: giveStrainVectorTranformationMtrx(FloatMatrix &answer,
                                                        const FloatMatrix &base,
                                                        bool transpose)
//
// returns transformation matrix for 3d - strains to another system of axes,
// given by base.
// In base (FloatMatrix[3,3]) there are on each column stored vectors of
// coordinate system to which we do transformation.
//
// If transpose == 1 we transpose base matrix before transforming
//
{
    FloatMatrix t;
    answer.resize(6, 6);
    answer.zero();

    if ( transpose ) {
        t.beTranspositionOf(base);
    } else {
        t = base;
    }

    answer.at(1, 1) = t.at(1, 1) * t.at(1, 1);
    answer.at(1, 2) = t.at(2, 1) * t.at(2, 1);
    answer.at(1, 3) = t.at(3, 1) * t.at(3, 1);
    answer.at(1, 4) = t.at(2, 1) * t.at(3, 1);
    answer.at(1, 5) = t.at(1, 1) * t.at(3, 1);
    answer.at(1, 6) = t.at(1, 1) * t.at(2, 1);

    answer.at(2, 1) = t.at(1, 2) * t.at(1, 2);
    answer.at(2, 2) = t.at(2, 2) * t.at(2, 2);
    answer.at(2, 3) = t.at(3, 2) * t.at(3, 2);
    answer.at(2, 4) = t.at(2, 2) * t.at(3, 2);
    answer.at(2, 5) = t.at(1, 2) * t.at(3, 2);
    answer.at(2, 6) = t.at(1, 2) * t.at(2, 2);

    answer.at(3, 1) = t.at(1, 3) * t.at(1, 3);
    answer.at(3, 2) = t.at(2, 3) * t.at(2, 3);
    answer.at(3, 3) = t.at(3, 3) * t.at(3, 3);
    answer.at(3, 4) = t.at(2, 3) * t.at(3, 3);
    answer.at(3, 5) = t.at(1, 3) * t.at(3, 3);
    answer.at(3, 6) = t.at(1, 3) * t.at(2, 3);

    answer.at(4, 1) = 2.0 * t.at(1, 2) * t.at(1, 3);
    answer.at(4, 2) = 2.0 * t.at(2, 2) * t.at(2, 3);
    answer.at(4, 3) = 2.0 * t.at(3, 2) * t.at(3, 3);
    answer.at(4, 4) = ( t.at(2, 2) * t.at(3, 3) + t.at(3, 2) * t.at(2, 3) );
    answer.at(4, 5) = ( t.at(1, 2) * t.at(3, 3) + t.at(3, 2) * t.at(1, 3) );
    answer.at(4, 6) = ( t.at(1, 2) * t.at(2, 3) + t.at(2, 2) * t.at(1, 3) );

    answer.at(5, 1) = 2.0 * t.at(1, 1) * t.at(1, 3);
    answer.at(5, 2) = 2.0 * t.at(2, 1) * t.at(2, 3);
    answer.at(5, 3) = 2.0 * t.at(3, 1) * t.at(3, 3);
    answer.at(5, 4) = ( t.at(2, 1) * t.at(3, 3) + t.at(3, 1) * t.at(2, 3) );
    answer.at(5, 5) = ( t.at(1, 1) * t.at(3, 3) + t.at(3, 1) * t.at(1, 3) );
    answer.at(5, 6) = ( t.at(1, 1) * t.at(2, 3) + t.at(2, 1) * t.at(1, 3) );

    answer.at(6, 1) = 2.0 * t.at(1, 1) * t.at(1, 2);
    answer.at(6, 2) = 2.0 * t.at(2, 1) * t.at(2, 2);
    answer.at(6, 3) = 2.0 * t.at(3, 1) * t.at(3, 2);
    answer.at(6, 4) = ( t.at(2, 1) * t.at(3, 2) + t.at(3, 1) * t.at(2, 2) );
    answer.at(6, 5) = ( t.at(1, 1) * t.at(3, 2) + t.at(3, 1) * t.at(1, 2) );
    answer.at(6, 6) = ( t.at(1, 1) * t.at(2, 2) + t.at(2, 1) * t.at(1, 2) );
}


void
StructuralMaterial :: giveStressVectorTranformationMtrx(FloatMatrix &answer,
                                                        const FloatMatrix &base,
                                                        bool transpose)
//
// returns transformation matrix for 3d - stress to another system of axes,
// given by base.
// In base (FloatMatrix[3,3]) there are on each column stored vectors of
// coordinate system to which we do transformation.
//
// If transpose == 1 we transpose base matrix before transforming
//
{
    FloatMatrix t;
    answer.resize(6, 6);
    answer.zero();

    if ( transpose ) {
        t.beTranspositionOf(base);
    } else {
        t = base;
    }

    answer.at(1, 1) = t.at(1, 1) * t.at(1, 1);
    answer.at(1, 2) = t.at(2, 1) * t.at(2, 1);
    answer.at(1, 3) = t.at(3, 1) * t.at(3, 1);
    answer.at(1, 4) = 2.0 * t.at(2, 1) * t.at(3, 1);
    answer.at(1, 5) = 2.0 * t.at(1, 1) * t.at(3, 1);
    answer.at(1, 6) = 2.0 * t.at(1, 1) * t.at(2, 1);

    answer.at(2, 1) = t.at(1, 2) * t.at(1, 2);
    answer.at(2, 2) = t.at(2, 2) * t.at(2, 2);
    answer.at(2, 3) = t.at(3, 2) * t.at(3, 2);
    answer.at(2, 4) = 2.0 * t.at(2, 2) * t.at(3, 2);
    answer.at(2, 5) = 2.0 * t.at(1, 2) * t.at(3, 2);
    answer.at(2, 6) = 2.0 * t.at(1, 2) * t.at(2, 2);

    answer.at(3, 1) = t.at(1, 3) * t.at(1, 3);
    answer.at(3, 2) = t.at(2, 3) * t.at(2, 3);
    answer.at(3, 3) = t.at(3, 3) * t.at(3, 3);
    answer.at(3, 4) = 2.0 * t.at(2, 3) * t.at(3, 3);
    answer.at(3, 5) = 2.0 * t.at(1, 3) * t.at(3, 3);
    answer.at(3, 6) = 2.0 * t.at(1, 3) * t.at(2, 3);

    answer.at(4, 1) = t.at(1, 2) * t.at(1, 3);
    answer.at(4, 2) = t.at(2, 2) * t.at(2, 3);
    answer.at(4, 3) = t.at(3, 2) * t.at(3, 3);
    answer.at(4, 4) = ( t.at(2, 2) * t.at(3, 3) + t.at(3, 2) * t.at(2, 3) );
    answer.at(4, 5) = ( t.at(1, 2) * t.at(3, 3) + t.at(3, 2) * t.at(1, 3) );
    answer.at(4, 6) = ( t.at(1, 2) * t.at(2, 3) + t.at(2, 2) * t.at(1, 3) );

    answer.at(5, 1) = t.at(1, 1) * t.at(1, 3);
    answer.at(5, 2) = t.at(2, 1) * t.at(2, 3);
    answer.at(5, 3) = t.at(3, 1) * t.at(3, 3);
    answer.at(5, 4) = ( t.at(2, 1) * t.at(3, 3) + t.at(3, 1) * t.at(2, 3) );
    answer.at(5, 5) = ( t.at(1, 1) * t.at(3, 3) + t.at(3, 1) * t.at(1, 3) );
    answer.at(5, 6) = ( t.at(1, 1) * t.at(2, 3) + t.at(2, 1) * t.at(1, 3) );

    answer.at(6, 1) = t.at(1, 1) * t.at(1, 2);
    answer.at(6, 2) = t.at(2, 1) * t.at(2, 2);
    answer.at(6, 3) = t.at(3, 1) * t.at(3, 2);
    answer.at(6, 4) = ( t.at(2, 1) * t.at(3, 2) + t.at(3, 1) * t.at(2, 2) );
    answer.at(6, 5) = ( t.at(1, 1) * t.at(3, 2) + t.at(3, 1) * t.at(1, 2) );
    answer.at(6, 6) = ( t.at(1, 1) * t.at(2, 2) + t.at(2, 1) * t.at(1, 2) );
}


void
StructuralMaterial :: givePlaneStressVectorTranformationMtrx(FloatMatrix &answer,
                                                             const FloatMatrix &base,
                                                             bool transpose)
//
// returns transformation matrix for 2d - stress to another system of axes,
// given by base.
// In base (FloatMatrix[2,2]) there are on each column stored vectors of
// coordinate system to which we do transformation.
//
// If transpose == 1 we transpose base matrix before transforming
//
{
    FloatMatrix t;
    answer.resize(3, 3);
    answer.zero();

    if ( transpose ) {
        t.beTranspositionOf(base);
    } else {
        t = base;
    }

    answer.at(1, 1) = t.at(1, 1) * t.at(1, 1);
    answer.at(1, 2) = t.at(2, 1) * t.at(2, 1);
    answer.at(1, 3) = 2.0 * t.at(1, 1) * t.at(2, 1);

    answer.at(2, 1) = t.at(1, 2) * t.at(1, 2);
    answer.at(2, 2) = t.at(2, 2) * t.at(2, 2);
    answer.at(2, 3) = 2.0 * t.at(1, 2) * t.at(2, 2);

    answer.at(3, 1) = t.at(1, 1) * t.at(1, 2);
    answer.at(3, 2) = t.at(2, 1) * t.at(2, 2);
    answer.at(3, 3) = t.at(1, 1) * t.at(2, 2) + t.at(2, 1) * t.at(1, 2);
}


void
StructuralMaterial :: transformStrainVectorTo(FloatArray &answer, const FloatMatrix &base,
                                              const FloatArray &strainVector, bool transpose)
//
// performs transformation of 3d-strain vector to another system of axes,
// given by base.
// In base (FloatMatrix[3,3]) there are on each column stored vectors of
// coordinate system to which we do transformation. These vectors must
// be expressed in the same coordinate system as strainVector
//
// If transpose == 1 we transpose base matrix before transforming
{
    FloatMatrix tt;

    StructuralMaterial :: giveStrainVectorTranformationMtrx(tt, base, transpose);
    answer.beProductOf(tt, strainVector);
}


void
StructuralMaterial :: transformStressVectorTo(FloatArray &answer, const FloatMatrix &base,
                                              const FloatArray &stressVector, bool transpose)
//
//
// performs transformation of 3d-stress vector to another system of axes,
// given by base.
// In base (FloatMatrix[3,3]) there are on each column stored vectors of
// coordinate system to which we do transformation. These vectors must
// be expressed in the same coordinate system as strainVector
// If transpose == 1 we transpose base matrix before transforming
//

{
    FloatMatrix tt;

    StructuralMaterial :: giveStressVectorTranformationMtrx(tt, base, transpose);
    answer.beProductOf(tt, stressVector);
}


void
StructuralMaterial :: sortPrincDirAndValCloseTo(FloatArray *pVal, FloatMatrix *pDir,
                                                FloatMatrix *toPDir)
//
// this method sorts newly computed principal values (pVal) and
// corresponding principal directions (pDir) to be closed to
// some (often previous) principal directions (toPDir).
//
// remark : pDir and toPDir should have eigen vectors stored in columns
// and normalized.
//
{
    int maxJ = 0, size;
    double cosine, maxCosine, swap;

#ifdef DEBUG
    if ( ( !pDir->isSquare() ) || ( !toPDir->isSquare() ) ) {
        OOFEM_ERROR("StructuralMaterial :: sortPrincDirandValCloseTo - Not square matrix");
    }

    if ( pDir->giveNumberOfRows() != toPDir->giveNumberOfRows() ) {
        OOFEM_ERROR("StructuralMaterial :: sortPrincDirandValCloseTo - Incompatible matrices");
    }

    if ( pDir->giveNumberOfRows() != pVal->giveSize() ) {
        OOFEM_ERROR("StructuralMaterial :: sortPrincDirandValCloseTo - Incompatible pVal Array size");
    }

#endif

    //
    // compute cosine matrix, where member i,j is cosine of angle
    // between toPDir i th eigen vector and j th pDir eigen vector
    //
    // sort pVal and pDir
    size = pDir->giveNumberOfRows();
    for ( int i = 1; i <= size - 1; i++ ) {
        // find closest pDir vector to toPDir i-th vector
        maxCosine = 0.0;
        for ( int j = i; j <= size; j++ ) {
            cosine = 0.;
            for ( int k = 1; k <= size; k++ ) {
                cosine += toPDir->at(k, i) * pDir->at(k, j);
            }

            cosine = fabs(cosine);
            if ( cosine > maxCosine ) {
                maxJ = j;
                maxCosine = cosine;
            }
        }

        // swap entries
        if ( maxJ != i ) {
            // swap eigenVectors and values
            swap = pVal->at(maxJ);
            pVal->at(maxJ) = pVal->at(i);
            pVal->at(i) = swap;
            for ( int k = 1; k <= size; k++ ) {
                swap = pDir->at(k, maxJ);
                pDir->at(k, maxJ) = pDir->at(k, i);
                pDir->at(k, i) = swap;
            }
        }
    }
}


int
StructuralMaterial :: setIPValue(const FloatArray &value, GaussPoint *gp, InternalStateType type)
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_StressTensor ) {
        status->letStressVectorBe(value);
        return 1;
    } else if ( type == IST_StrainTensor ) {
        status->letStrainVectorBe(value);
        return 1;
    } else if ( type == IST_StressTensorTemp ) {
        status->letTempStressVectorBe(value);
        return 1;
    } else if ( type == IST_StrainTensorTemp ) {
        status->letTempStrainVectorBe(value);
        return 1;
    } else {
        return 0;
    }
}


int
StructuralMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_StressTensor ) {
        StructuralMaterial :: giveFullSymVectorForm( answer, status->giveStressVector(), gp->giveMaterialMode() );
        return 1;
    } else if ( type == IST_StressTensor_Reduced ) {
        answer = status->giveStressVector();
        return 1;
    } else if ( type == IST_vonMisesStress ) {
        ///@todo What about the stress meassure in large deformations here? The internal state type should specify "Cauchy" or something.
        answer.resize(1);
        answer.at(1) = this->computeVonMisesStress( & status->giveStressVector() );
        return 1;
    } else if ( type == IST_StrainTensor ) {
        ///@todo Fill in correct full form values here! This just adds zeros!
        StructuralMaterial :: giveFullSymVectorForm( answer, status->giveStrainVector(), gp->giveMaterialMode() );
        if ( gp->giveMaterialMode() == _PlaneStress ) {
            double Nxy = this->give(NYxy, gp);
            double Nxz = this->give(NYxz, gp);
            double Nyz = this->give(NYyz, gp);
            double Nyx = Nxy * this->give(Ey, gp) / this->give(Ex, gp);
            answer.at(3) = ( -( Nxz + Nxy * Nyz ) * answer.at(1) - ( Nyz + Nxz * Nyx ) * answer.at(2) ) / ( 1. - Nxy * Nyx );
        }
        return 1;
    } else if ( type == IST_StrainTensor_Reduced ) {
        ///@todo Fill in correct full form values here! This just adds zeros!
        answer = status->giveStrainVector();
        return 1;
    } else if ( type == IST_StressTensorTemp ) {
        ///@todo Fill in correct full form values here! This just adds zeros!
        StructuralMaterial :: giveFullSymVectorForm( answer, status->giveTempStressVector(), gp->giveMaterialMode() );
        return 1;
    } else if ( type == IST_StrainTensorTemp ) {
        ///@todo Fill in correct full form values here! This just adds zeros!
        StructuralMaterial :: giveFullSymVectorForm( answer, status->giveTempStrainVector(), gp->giveMaterialMode() );
        return 1;
    } else if ( type == IST_PrincipalStressTensor || type == IST_PrincipalStressTempTensor ) {
        FloatArray s;

        if ( type == IST_PrincipalStressTensor ) {
            ///@todo Fill in correct full form values here! This just adds zeros!
            StructuralMaterial :: giveFullSymVectorForm( s, status->giveStressVector(), gp->giveMaterialMode() );
        } else {
            ///@todo Fill in correct full form values here! This just adds zeros!
            StructuralMaterial :: giveFullSymVectorForm( s, status->giveTempStressVector(), gp->giveMaterialMode() );
        }

        this->computePrincipalValues(answer, s, principal_stress);
        return 1;
    } else if ( type == IST_PrincipalStrainTensor || type == IST_PrincipalStrainTempTensor ) {
        FloatArray s;

        if ( type == IST_PrincipalStrainTensor ) {
            ///@todo Fill in correct full form values here! This just adds zeros!
            StructuralMaterial :: giveFullSymVectorForm( s, status->giveStrainVector(), gp->giveMaterialMode() );
        } else {
            ///@todo Fill in correct full form values here! This just adds zeros!
            StructuralMaterial :: giveFullSymVectorForm( s, status->giveTempStrainVector(), gp->giveMaterialMode() );
        }

        this->computePrincipalValues(answer, s, principal_strain);
        return 1;
    } else if ( type == IST_Temperature ) {
        /* add external source, if provided, such as staggered analysis */
        FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();
        FM_FieldPtr tf;
        int err;
        if ( ( tf = fm->giveField(FT_Temperature) ) ) {
            // temperature field registered
            FloatArray gcoords, et2;
            static_cast< StructuralElement * >( gp->giveElement() )->computeGlobalCoordinates( gcoords, * gp->giveCoordinates() );
            if ( ( err = tf->evaluateAt(answer, gcoords, VM_Total, tStep) ) ) {
                OOFEM_ERROR3("StructuralMaterial :: giveIPValue: tf->evaluateAt failed, element %d, error code %d", gp->giveElement()->giveNumber(), err);
            }
        } else {
            answer.resize(1);
            answer.zero();
        }

        return 1;
    } else if ( type == IST_CylindricalStressTensor || type == IST_CylindricalStrainTensor ) {
        FloatArray gc, val = status->giveStressVector();
        FloatMatrix base(3, 3);
        static_cast< StructuralElement * >( gp->giveElement() )->computeGlobalCoordinates( gc, * gp->giveCoordinates() );
        double l = sqrt( gc.at(1) * gc.at(1) + gc.at(2) * gc.at(2) );
        if ( l > 1.e-4 ) {
            base.at(1, 1) = gc.at(1) / l;
            base.at(2, 1) = gc.at(2) / l;
            base.at(3, 1) = 0.0;

            base.at(1, 2) = -1.0 * base.at(2, 1);
            base.at(2, 2) = base.at(1, 1);
            base.at(3, 2) = 0.0;

            base.at(1, 3) = 0.0;
            base.at(2, 3) = 0.0;
            base.at(3, 3) = 1.0;

            if ( type == IST_CylindricalStressTensor ) {
                this->transformStressVectorTo(answer, base, val, false);
            } else {
                this->transformStrainVectorTo(answer, base, val, false);
            }
        } else {
            answer = val;
        }

        return 1;
    } else if ( type == IST_DeformationGradientTensor ) {
        answer = status->giveFVector();
        return 1;
    } else if ( type == IST_FirstPKStressTensor ) {
        answer = status->givePVector();
        return 1;
    } else {
        return Material :: giveIPValue(answer, gp, type, tStep);
    }
    return 0;
}

void
StructuralMaterial :: computeStressIndependentStrainVector(FloatArray &answer,
                                                           GaussPoint *gp, TimeStep *tStep, ValueModeType mode)
{
    FloatArray et, e0, eigenstrain;
    FloatMatrix GCS;
    MaterialMode matmode = gp->giveMaterialMode();
    Element *elem = gp->giveElement();
    StructuralElement *selem = dynamic_cast< StructuralElement * >( gp->giveElement() );

    answer.clear();

    if ( tStep->giveIntrinsicTime() < this->castingTime ) {
        answer.zero();
        return;
    }

    //sum up all prescribed temperatures over an element
    //elem->computeResultingIPTemperatureAt(et, tStep, gp, mode);
    if ( selem ) {
        selem->computeResultingIPTemperatureAt(et, tStep, gp, mode);        // HUHU
    }

    //sum up all prescribed eigenstrain over an element
    if ( selem ) {
        selem->computeResultingIPEigenstrainAt(eigenstrain, tStep, gp, mode);
    }

    if ( eigenstrain.giveSize() != 0 && eigenstrain.giveSize() != giveSizeOfVoigtSymVector(matmode) ) {
        OOFEM_ERROR5( "StructuralMaterial :: Number of given eigenstrain components %d is different than required %d by material mode %s, element %d", eigenstrain.giveSize(), giveSizeOfVoigtSymVector(matmode), __MaterialModeToString(matmode), elem->giveNumber() );
    }

    /* add external source, if provided */
    FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();
    FM_FieldPtr tf;

    if ( ( tf = fm->giveField(FT_Temperature) ) ) {
        // temperature field registered
        FloatArray gcoords, et2;
        int err;
        elem->computeGlobalCoordinates( gcoords, * gp->giveCoordinates() );
        if ( ( err = tf->evaluateAt(et2, gcoords, mode, tStep) ) ) {
            OOFEM_ERROR3("StructuralMaterial :: computeStressIndependentStrainVector: tf->evaluateAt failed, element %d, error code %d", elem->giveNumber(), err);
        }

        if ( et2.isNotEmpty() ) {
            if ( et.isEmpty() ) {
                et = et2;
            } else {
                et.at(1) += et2.at(1);
            }
        }
    }


    if ( et.giveSize() ) { //found temperature boundary conditions or prescribed field
        FloatArray fullAnswer;

        this->giveThermalDilatationVector(e0, gp, tStep);

        if ( e0.giveSize() ) {
            fullAnswer = e0;
            if ( mode == VM_Total ) {
                fullAnswer.times(et.at(1) - this->referenceTemperature);
            } else {
                fullAnswer.times( et.at(1) );
            }

            StructuralMaterial :: giveReducedSymVectorForm( answer, fullAnswer, gp->giveMaterialMode() );
        }
    }

    //join temperature and eigenstrain vectors, compare vector sizes
    if ( answer.giveSize() ) {
        if ( eigenstrain.giveSize() ) {
            if ( answer.giveSize() != eigenstrain.giveSize() ) {
                OOFEM_ERROR4( "StructuralMaterial :: Vector of temperature strains has the size %d which is different with the size of eigenstrain vector %d, element %d", answer.giveSize(), eigenstrain.giveSize(), elem->giveNumber() );
            }

            answer.add(eigenstrain);
        }
    } else {
        if ( eigenstrain.giveSize() ) {
            answer = eigenstrain;
        }
    }
}


void
StructuralMaterial :: giveFullSymVectorForm(FloatArray &answer, const FloatArray &vec, MaterialMode matMode)
{
    if(vec.giveSize() == 6) {
        // If we use default 3D implementation to treat e.g. plane strain.
        answer = vec;
    }
    else {
        IntArray indx;
        answer.resize( StructuralMaterial :: giveVoigtSymVectorMask(indx, matMode) );
        answer.zero();
        answer.assemble(vec, indx);
    }
}


void
StructuralMaterial :: giveFullVectorForm(FloatArray &answer, const FloatArray &vec, MaterialMode matMode)
{
    IntArray indx;
    answer.resize( StructuralMaterial :: giveVoigtVectorMask(indx, matMode) );
    answer.zero();
    answer.assemble(vec, indx);
}


void
StructuralMaterial :: giveFullVectorFormF(FloatArray &answer, const FloatArray &vec, MaterialMode matMode)
{
    IntArray indx;
    answer.resize(9);
    answer.at(1) = answer.at(2) = answer.at(3) = 1.0;   // set diagonal terms

    StructuralMaterial :: giveVoigtVectorMask(indx, matMode);
    for ( int i = 1; i <= indx.giveSize(); i++ ) {
        answer.at( indx.at(i) ) = vec.at(i);
    }
}


void
StructuralMaterial :: giveReducedVectorForm(FloatArray &answer, const FloatArray &vec, MaterialMode matMode)
{
    IntArray indx;
    StructuralMaterial :: giveVoigtVectorMask(indx, matMode);
    answer.resize( indx.giveSize() );
    for ( int i = 1; i <= indx.giveSize(); i++ ) {
        answer.at(i) = vec.at( indx.at(i) );
    }
}


void
StructuralMaterial :: giveReducedSymVectorForm(FloatArray &answer, const FloatArray &vec, MaterialMode matMode)
{
    IntArray indx;
    StructuralMaterial :: giveVoigtSymVectorMask(indx, matMode);

    if(indx.giveSize() == vec.giveSize()) {
        answer = vec;
    }
    else {
        answer.resize( indx.giveSize() );
        for ( int i = 1; i <= indx.giveSize(); i++ ) {
            answer.at(i) = vec.at( indx.at(i) );
        }
    }
}


void
StructuralMaterial :: giveFullSymMatrixForm(FloatMatrix &answer, const FloatMatrix &red, MaterialMode matMode)
{
    IntArray indx;
    int size = StructuralMaterial :: giveVoigtSymVectorMask(indx, matMode);
    answer.resize(size, size);
    answer.zero();
    answer.assemble(red, indx, indx);
}


void
StructuralMaterial :: giveReducedMatrixForm(FloatMatrix &answer, const FloatMatrix &full, MaterialMode matMode)
{
    IntArray indx;
    StructuralMaterial :: giveVoigtVectorMask(indx, matMode);
    answer.beSubMatrixOf(full, indx, indx);
}


void
StructuralMaterial :: giveReducedSymMatrixForm(FloatMatrix &answer, const FloatMatrix &full, MaterialMode matMode)
{
    IntArray indx;
    StructuralMaterial :: giveVoigtSymVectorMask(indx, matMode);
    answer.beSubMatrixOf(full, indx, indx);
}



IRResultType
StructuralMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

#  ifdef VERBOSE
    // VERBOSE_PRINT1 ("Instanciating material ",this->giveNumber())
#  endif
    this->Material :: initializeFrom(ir);

    referenceTemperature = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, referenceTemperature, _IFT_StructuralMaterial_referencetemperature);

    return IRRT_OK;
}


void
StructuralMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    Material :: giveInputRecord(input);
    input.setField(this->referenceTemperature, _IFT_StructuralMaterial_referencetemperature);
}
} // end namespace oofem
