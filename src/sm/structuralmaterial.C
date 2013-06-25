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

#include "structuralmaterial.h"
#include "structuralcrosssection.h"
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

namespace oofem {
int
StructuralMaterial :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports given mode
//
{
    return mode == _3dMat        || mode == _PlaneStress ||
           mode == _PlaneStrain  || mode == _1dMat       ||
           mode == _2dPlateLayer || mode == _2dBeamLayer ||
           mode == _3dShellLayer || mode == _1dFiber;
}


void 
StructuralMaterial :: giveFirstPKStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
        const FloatArray &reducedvF, TimeStep *tStep) 
{
    // Default implementation used if this method is not overloaded by the particular material model.
    // 1) Compute Green-Lagrange strain and call standard method for small strains. 
    // 2) Treat stress as second Piola-Kirchhoff stress and convert to first Piola-Kirchhoff stress.
    // 3) Set state variables F, P 
    
    FloatArray reducedvE, reducedvS;
    FloatArray vF, vE;
    FloatMatrix F, E;
    StructuralMaterial :: giveFullVectorFormF(vF, reducedvF, gp->giveMaterialMode()); // 9
    F.beMatrixForm(vF);
    StructuralMaterial :: computeGreenLagrangeStrain(E,F);  // 3x3
    vE.beReducedVectorFormOfStrain(E);      // 6
    StructuralMaterial :: giveReducedSymVectorForm(reducedvE, vE, gp->giveMaterialMode()); //reduced
    
    this->giveRealStressVector(reducedvS, form, gp, reducedvE, tStep); // Treat stress obtained as second PK stress
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    
    // Compute first PK stress from second PK stress
    this->convert_S_2_P(answer, reducedvS, reducedvF, gp->giveMaterialMode());

    status->letTempPVectorBe(answer);
    status->letTempFVectorBe(reducedvF);
}


void 
StructuralMaterial :: computeGreenLagrangeStrain(FloatMatrix &answer, FloatMatrix &F)
{
    // Computes the Green-Lagrange strain tensor, E = 0.5*( C - I ), on matrix form. 
    // The size of the output matrix will be the same as the input matrix.

    answer.beTProductOf(F, F);    // C - Right Caucy-Green deformation tensor

    if ( answer.giveNumberOfRows() == 3 ) {
        answer.at(1, 1) -= 1.0;
        answer.at(2, 2) -= 1.0;
        answer.at(3, 3) -= 1.0;
    } else if ( answer.giveNumberOfRows() == 2 ) {
        answer.at(1, 1) -= 1.0;
        answer.at(2, 2) -= 1.0;
    } else if ( answer.giveNumberOfRows() == 1 ) {
        answer.at(1, 1) -= 1.0;
    } else {
        OOFEM_ERROR2("StructuralMaterial :: computeGreenLagrangeStrain - wrong size of input matrix (num rows = %d)", answer.giveNumberOfRows());
    }
    answer.times(0.5);
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

        answer.resize(9,9);
        answer(0,0) = F(0)*C(0,0)*F(0) + F(0)*C(0,5)*F(5) + F(0)*C(0,4)*F(4) + F(5)*C(5,0)*F(0) + F(5)*C(5,5)*F(5) + F(5)*C(5,4)*F(4) + F(4)*C(4,0)*F(0) + F(4)*C(4,5)*F(5) + F(4)*C(4,4)*F(4) + S(0); 
        answer(0,1) = F(0)*C(0,5)*F(8) + F(0)*C(0,1)*F(1) + F(0)*C(0,3)*F(3) + F(5)*C(5,5)*F(8) + F(5)*C(5,1)*F(1) + F(5)*C(5,3)*F(3) + F(4)*C(4,5)*F(8) + F(4)*C(4,1)*F(1) + F(4)*C(4,3)*F(3) + 0.0; 
        answer(0,2) = F(0)*C(0,4)*F(7) + F(0)*C(0,3)*F(6) + F(0)*C(0,2)*F(2) + F(5)*C(5,4)*F(7) + F(5)*C(5,3)*F(6) + F(5)*C(5,2)*F(2) + F(4)*C(4,4)*F(7) + F(4)*C(4,3)*F(6) + F(4)*C(4,2)*F(2) + 0.0; 
        answer(0,3) = F(0)*C(0,4)*F(8) + F(0)*C(0,3)*F(1) + F(0)*C(0,2)*F(3) + F(5)*C(5,4)*F(8) + F(5)*C(5,3)*F(1) + F(5)*C(5,2)*F(3) + F(4)*C(4,4)*F(8) + F(4)*C(4,3)*F(1) + F(4)*C(4,2)*F(3) + 0.0; 
        answer(0,4) = F(0)*C(0,4)*F(0) + F(0)*C(0,3)*F(5) + F(0)*C(0,2)*F(4) + F(5)*C(5,4)*F(0) + F(5)*C(5,3)*F(5) + F(5)*C(5,2)*F(4) + F(4)*C(4,4)*F(0) + F(4)*C(4,3)*F(5) + F(4)*C(4,2)*F(4) + S(4); 
        answer(0,5) = F(0)*C(0,5)*F(0) + F(0)*C(0,1)*F(5) + F(0)*C(0,3)*F(4) + F(5)*C(5,5)*F(0) + F(5)*C(5,1)*F(5) + F(5)*C(5,3)*F(4) + F(4)*C(4,5)*F(0) + F(4)*C(4,1)*F(5) + F(4)*C(4,3)*F(4) + S(5); 
        answer(0,6) = F(0)*C(0,5)*F(7) + F(0)*C(0,1)*F(6) + F(0)*C(0,3)*F(2) + F(5)*C(5,5)*F(7) + F(5)*C(5,1)*F(6) + F(5)*C(5,3)*F(2) + F(4)*C(4,5)*F(7) + F(4)*C(4,1)*F(6) + F(4)*C(4,3)*F(2) + 0.0; 
        answer(0,7) = F(0)*C(0,0)*F(7) + F(0)*C(0,5)*F(6) + F(0)*C(0,4)*F(2) + F(5)*C(5,0)*F(7) + F(5)*C(5,5)*F(6) + F(5)*C(5,4)*F(2) + F(4)*C(4,0)*F(7) + F(4)*C(4,5)*F(6) + F(4)*C(4,4)*F(2) + 0.0; 
        answer(0,8) = F(0)*C(0,0)*F(8) + F(0)*C(0,5)*F(1) + F(0)*C(0,4)*F(3) + F(5)*C(5,0)*F(8) + F(5)*C(5,5)*F(1) + F(5)*C(5,4)*F(3) + F(4)*C(4,0)*F(8) + F(4)*C(4,5)*F(1) + F(4)*C(4,4)*F(3) + 0.0; 
        answer(1,0) = F(8)*C(5,0)*F(0) + F(8)*C(5,5)*F(5) + F(8)*C(5,4)*F(4) + F(1)*C(1,0)*F(0) + F(1)*C(1,5)*F(5) + F(1)*C(1,4)*F(4) + F(3)*C(3,0)*F(0) + F(3)*C(3,5)*F(5) + F(3)*C(3,4)*F(4) + 0.0; 
        answer(1,1) = F(8)*C(5,5)*F(8) + F(8)*C(5,1)*F(1) + F(8)*C(5,3)*F(3) + F(1)*C(1,5)*F(8) + F(1)*C(1,1)*F(1) + F(1)*C(1,3)*F(3) + F(3)*C(3,5)*F(8) + F(3)*C(3,1)*F(1) + F(3)*C(3,3)*F(3) + S(1); 
        answer(1,2) = F(8)*C(5,4)*F(7) + F(8)*C(5,3)*F(6) + F(8)*C(5,2)*F(2) + F(1)*C(1,4)*F(7) + F(1)*C(1,3)*F(6) + F(1)*C(1,2)*F(2) + F(3)*C(3,4)*F(7) + F(3)*C(3,3)*F(6) + F(3)*C(3,2)*F(2) + 0.0; 
        answer(1,3) = F(8)*C(5,4)*F(8) + F(8)*C(5,3)*F(1) + F(8)*C(5,2)*F(3) + F(1)*C(1,4)*F(8) + F(1)*C(1,3)*F(1) + F(1)*C(1,2)*F(3) + F(3)*C(3,4)*F(8) + F(3)*C(3,3)*F(1) + F(3)*C(3,2)*F(3) + S(3); 
        answer(1,4) = F(8)*C(5,4)*F(0) + F(8)*C(5,3)*F(5) + F(8)*C(5,2)*F(4) + F(1)*C(1,4)*F(0) + F(1)*C(1,3)*F(5) + F(1)*C(1,2)*F(4) + F(3)*C(3,4)*F(0) + F(3)*C(3,3)*F(5) + F(3)*C(3,2)*F(4) + 0.0; 
        answer(1,5) = F(8)*C(5,5)*F(0) + F(8)*C(5,1)*F(5) + F(8)*C(5,3)*F(4) + F(1)*C(1,5)*F(0) + F(1)*C(1,1)*F(5) + F(1)*C(1,3)*F(4) + F(3)*C(3,5)*F(0) + F(3)*C(3,1)*F(5) + F(3)*C(3,3)*F(4) + 0.0; 
        answer(1,6) = F(8)*C(5,5)*F(7) + F(8)*C(5,1)*F(6) + F(8)*C(5,3)*F(2) + F(1)*C(1,5)*F(7) + F(1)*C(1,1)*F(6) + F(1)*C(1,3)*F(2) + F(3)*C(3,5)*F(7) + F(3)*C(3,1)*F(6) + F(3)*C(3,3)*F(2) + 0.0; 
        answer(1,7) = F(8)*C(5,0)*F(7) + F(8)*C(5,5)*F(6) + F(8)*C(5,4)*F(2) + F(1)*C(1,0)*F(7) + F(1)*C(1,5)*F(6) + F(1)*C(1,4)*F(2) + F(3)*C(3,0)*F(7) + F(3)*C(3,5)*F(6) + F(3)*C(3,4)*F(2) + 0.0; 
        answer(1,8) = F(8)*C(5,0)*F(8) + F(8)*C(5,5)*F(1) + F(8)*C(5,4)*F(3) + F(1)*C(1,0)*F(8) + F(1)*C(1,5)*F(1) + F(1)*C(1,4)*F(3) + F(3)*C(3,0)*F(8) + F(3)*C(3,5)*F(1) + F(3)*C(3,4)*F(3) + S(5); 
        answer(2,0) = F(7)*C(4,0)*F(0) + F(7)*C(4,5)*F(5) + F(7)*C(4,4)*F(4) + F(6)*C(3,0)*F(0) + F(6)*C(3,5)*F(5) + F(6)*C(3,4)*F(4) + F(2)*C(2,0)*F(0) + F(2)*C(2,5)*F(5) + F(2)*C(2,4)*F(4) + 0.0; 
        answer(2,1) = F(7)*C(4,5)*F(8) + F(7)*C(4,1)*F(1) + F(7)*C(4,3)*F(3) + F(6)*C(3,5)*F(8) + F(6)*C(3,1)*F(1) + F(6)*C(3,3)*F(3) + F(2)*C(2,5)*F(8) + F(2)*C(2,1)*F(1) + F(2)*C(2,3)*F(3) + 0.0; 
        answer(2,2) = F(7)*C(4,4)*F(7) + F(7)*C(4,3)*F(6) + F(7)*C(4,2)*F(2) + F(6)*C(3,4)*F(7) + F(6)*C(3,3)*F(6) + F(6)*C(3,2)*F(2) + F(2)*C(2,4)*F(7) + F(2)*C(2,3)*F(6) + F(2)*C(2,2)*F(2) + S(2); 
        answer(2,3) = F(7)*C(4,4)*F(8) + F(7)*C(4,3)*F(1) + F(7)*C(4,2)*F(3) + F(6)*C(3,4)*F(8) + F(6)*C(3,3)*F(1) + F(6)*C(3,2)*F(3) + F(2)*C(2,4)*F(8) + F(2)*C(2,3)*F(1) + F(2)*C(2,2)*F(3) + 0.0; 
        answer(2,4) = F(7)*C(4,4)*F(0) + F(7)*C(4,3)*F(5) + F(7)*C(4,2)*F(4) + F(6)*C(3,4)*F(0) + F(6)*C(3,3)*F(5) + F(6)*C(3,2)*F(4) + F(2)*C(2,4)*F(0) + F(2)*C(2,3)*F(5) + F(2)*C(2,2)*F(4) + 0.0; 
        answer(2,5) = F(7)*C(4,5)*F(0) + F(7)*C(4,1)*F(5) + F(7)*C(4,3)*F(4) + F(6)*C(3,5)*F(0) + F(6)*C(3,1)*F(5) + F(6)*C(3,3)*F(4) + F(2)*C(2,5)*F(0) + F(2)*C(2,1)*F(5) + F(2)*C(2,3)*F(4) + 0.0; 
        answer(2,6) = F(7)*C(4,5)*F(7) + F(7)*C(4,1)*F(6) + F(7)*C(4,3)*F(2) + F(6)*C(3,5)*F(7) + F(6)*C(3,1)*F(6) + F(6)*C(3,3)*F(2) + F(2)*C(2,5)*F(7) + F(2)*C(2,1)*F(6) + F(2)*C(2,3)*F(2) + S(3); 
        answer(2,7) = F(7)*C(4,0)*F(7) + F(7)*C(4,5)*F(6) + F(7)*C(4,4)*F(2) + F(6)*C(3,0)*F(7) + F(6)*C(3,5)*F(6) + F(6)*C(3,4)*F(2) + F(2)*C(2,0)*F(7) + F(2)*C(2,5)*F(6) + F(2)*C(2,4)*F(2) + S(4); 
        answer(2,8) = F(7)*C(4,0)*F(8) + F(7)*C(4,5)*F(1) + F(7)*C(4,4)*F(3) + F(6)*C(3,0)*F(8) + F(6)*C(3,5)*F(1) + F(6)*C(3,4)*F(3) + F(2)*C(2,0)*F(8) + F(2)*C(2,5)*F(1) + F(2)*C(2,4)*F(3) + 0.0; 
        answer(3,0) = F(8)*C(4,0)*F(0) + F(8)*C(4,5)*F(5) + F(8)*C(4,4)*F(4) + F(1)*C(3,0)*F(0) + F(1)*C(3,5)*F(5) + F(1)*C(3,4)*F(4) + F(3)*C(2,0)*F(0) + F(3)*C(2,5)*F(5) + F(3)*C(2,4)*F(4) + 0.0; 
        answer(3,1) = F(8)*C(4,5)*F(8) + F(8)*C(4,1)*F(1) + F(8)*C(4,3)*F(3) + F(1)*C(3,5)*F(8) + F(1)*C(3,1)*F(1) + F(1)*C(3,3)*F(3) + F(3)*C(2,5)*F(8) + F(3)*C(2,1)*F(1) + F(3)*C(2,3)*F(3) + S(3); 
        answer(3,2) = F(8)*C(4,4)*F(7) + F(8)*C(4,3)*F(6) + F(8)*C(4,2)*F(2) + F(1)*C(3,4)*F(7) + F(1)*C(3,3)*F(6) + F(1)*C(3,2)*F(2) + F(3)*C(2,4)*F(7) + F(3)*C(2,3)*F(6) + F(3)*C(2,2)*F(2) + 0.0; 
        answer(3,3) = F(8)*C(4,4)*F(8) + F(8)*C(4,3)*F(1) + F(8)*C(4,2)*F(3) + F(1)*C(3,4)*F(8) + F(1)*C(3,3)*F(1) + F(1)*C(3,2)*F(3) + F(3)*C(2,4)*F(8) + F(3)*C(2,3)*F(1) + F(3)*C(2,2)*F(3) + S(2); 
        answer(3,4) = F(8)*C(4,4)*F(0) + F(8)*C(4,3)*F(5) + F(8)*C(4,2)*F(4) + F(1)*C(3,4)*F(0) + F(1)*C(3,3)*F(5) + F(1)*C(3,2)*F(4) + F(3)*C(2,4)*F(0) + F(3)*C(2,3)*F(5) + F(3)*C(2,2)*F(4) + 0.0; 
        answer(3,5) = F(8)*C(4,5)*F(0) + F(8)*C(4,1)*F(5) + F(8)*C(4,3)*F(4) + F(1)*C(3,5)*F(0) + F(1)*C(3,1)*F(5) + F(1)*C(3,3)*F(4) + F(3)*C(2,5)*F(0) + F(3)*C(2,1)*F(5) + F(3)*C(2,3)*F(4) + 0.0; 
        answer(3,6) = F(8)*C(4,5)*F(7) + F(8)*C(4,1)*F(6) + F(8)*C(4,3)*F(2) + F(1)*C(3,5)*F(7) + F(1)*C(3,1)*F(6) + F(1)*C(3,3)*F(2) + F(3)*C(2,5)*F(7) + F(3)*C(2,1)*F(6) + F(3)*C(2,3)*F(2) + 0.0; 
        answer(3,7) = F(8)*C(4,0)*F(7) + F(8)*C(4,5)*F(6) + F(8)*C(4,4)*F(2) + F(1)*C(3,0)*F(7) + F(1)*C(3,5)*F(6) + F(1)*C(3,4)*F(2) + F(3)*C(2,0)*F(7) + F(3)*C(2,5)*F(6) + F(3)*C(2,4)*F(2) + 0.0; 
        answer(3,8) = F(8)*C(4,0)*F(8) + F(8)*C(4,5)*F(1) + F(8)*C(4,4)*F(3) + F(1)*C(3,0)*F(8) + F(1)*C(3,5)*F(1) + F(1)*C(3,4)*F(3) + F(3)*C(2,0)*F(8) + F(3)*C(2,5)*F(1) + F(3)*C(2,4)*F(3) + S(4); 
        answer(4,0) = F(0)*C(4,0)*F(0) + F(0)*C(4,5)*F(5) + F(0)*C(4,4)*F(4) + F(5)*C(3,0)*F(0) + F(5)*C(3,5)*F(5) + F(5)*C(3,4)*F(4) + F(4)*C(2,0)*F(0) + F(4)*C(2,5)*F(5) + F(4)*C(2,4)*F(4) + S(4); 
        answer(4,1) = F(0)*C(4,5)*F(8) + F(0)*C(4,1)*F(1) + F(0)*C(4,3)*F(3) + F(5)*C(3,5)*F(8) + F(5)*C(3,1)*F(1) + F(5)*C(3,3)*F(3) + F(4)*C(2,5)*F(8) + F(4)*C(2,1)*F(1) + F(4)*C(2,3)*F(3) + 0.0; 
        answer(4,2) = F(0)*C(4,4)*F(7) + F(0)*C(4,3)*F(6) + F(0)*C(4,2)*F(2) + F(5)*C(3,4)*F(7) + F(5)*C(3,3)*F(6) + F(5)*C(3,2)*F(2) + F(4)*C(2,4)*F(7) + F(4)*C(2,3)*F(6) + F(4)*C(2,2)*F(2) + 0.0; 
        answer(4,3) = F(0)*C(4,4)*F(8) + F(0)*C(4,3)*F(1) + F(0)*C(4,2)*F(3) + F(5)*C(3,4)*F(8) + F(5)*C(3,3)*F(1) + F(5)*C(3,2)*F(3) + F(4)*C(2,4)*F(8) + F(4)*C(2,3)*F(1) + F(4)*C(2,2)*F(3) + 0.0; 
        answer(4,4) = F(0)*C(4,4)*F(0) + F(0)*C(4,3)*F(5) + F(0)*C(4,2)*F(4) + F(5)*C(3,4)*F(0) + F(5)*C(3,3)*F(5) + F(5)*C(3,2)*F(4) + F(4)*C(2,4)*F(0) + F(4)*C(2,3)*F(5) + F(4)*C(2,2)*F(4) + S(2); 
        answer(4,5) = F(0)*C(4,5)*F(0) + F(0)*C(4,1)*F(5) + F(0)*C(4,3)*F(4) + F(5)*C(3,5)*F(0) + F(5)*C(3,1)*F(5) + F(5)*C(3,3)*F(4) + F(4)*C(2,5)*F(0) + F(4)*C(2,1)*F(5) + F(4)*C(2,3)*F(4) + S(3); 
        answer(4,6) = F(0)*C(4,5)*F(7) + F(0)*C(4,1)*F(6) + F(0)*C(4,3)*F(2) + F(5)*C(3,5)*F(7) + F(5)*C(3,1)*F(6) + F(5)*C(3,3)*F(2) + F(4)*C(2,5)*F(7) + F(4)*C(2,1)*F(6) + F(4)*C(2,3)*F(2) + 0.0; 
        answer(4,7) = F(0)*C(4,0)*F(7) + F(0)*C(4,5)*F(6) + F(0)*C(4,4)*F(2) + F(5)*C(3,0)*F(7) + F(5)*C(3,5)*F(6) + F(5)*C(3,4)*F(2) + F(4)*C(2,0)*F(7) + F(4)*C(2,5)*F(6) + F(4)*C(2,4)*F(2) + 0.0; 
        answer(4,8) = F(0)*C(4,0)*F(8) + F(0)*C(4,5)*F(1) + F(0)*C(4,4)*F(3) + F(5)*C(3,0)*F(8) + F(5)*C(3,5)*F(1) + F(5)*C(3,4)*F(3) + F(4)*C(2,0)*F(8) + F(4)*C(2,5)*F(1) + F(4)*C(2,4)*F(3) + 0.0; 
        answer(5,0) = F(0)*C(5,0)*F(0) + F(0)*C(5,5)*F(5) + F(0)*C(5,4)*F(4) + F(5)*C(1,0)*F(0) + F(5)*C(1,5)*F(5) + F(5)*C(1,4)*F(4) + F(4)*C(3,0)*F(0) + F(4)*C(3,5)*F(5) + F(4)*C(3,4)*F(4) + S(5); 
        answer(5,1) = F(0)*C(5,5)*F(8) + F(0)*C(5,1)*F(1) + F(0)*C(5,3)*F(3) + F(5)*C(1,5)*F(8) + F(5)*C(1,1)*F(1) + F(5)*C(1,3)*F(3) + F(4)*C(3,5)*F(8) + F(4)*C(3,1)*F(1) + F(4)*C(3,3)*F(3) + 0.0; 
        answer(5,2) = F(0)*C(5,4)*F(7) + F(0)*C(5,3)*F(6) + F(0)*C(5,2)*F(2) + F(5)*C(1,4)*F(7) + F(5)*C(1,3)*F(6) + F(5)*C(1,2)*F(2) + F(4)*C(3,4)*F(7) + F(4)*C(3,3)*F(6) + F(4)*C(3,2)*F(2) + 0.0; 
        answer(5,3) = F(0)*C(5,4)*F(8) + F(0)*C(5,3)*F(1) + F(0)*C(5,2)*F(3) + F(5)*C(1,4)*F(8) + F(5)*C(1,3)*F(1) + F(5)*C(1,2)*F(3) + F(4)*C(3,4)*F(8) + F(4)*C(3,3)*F(1) + F(4)*C(3,2)*F(3) + 0.0; 
        answer(5,4) = F(0)*C(5,4)*F(0) + F(0)*C(5,3)*F(5) + F(0)*C(5,2)*F(4) + F(5)*C(1,4)*F(0) + F(5)*C(1,3)*F(5) + F(5)*C(1,2)*F(4) + F(4)*C(3,4)*F(0) + F(4)*C(3,3)*F(5) + F(4)*C(3,2)*F(4) + S(3); 
        answer(5,5) = F(0)*C(5,5)*F(0) + F(0)*C(5,1)*F(5) + F(0)*C(5,3)*F(4) + F(5)*C(1,5)*F(0) + F(5)*C(1,1)*F(5) + F(5)*C(1,3)*F(4) + F(4)*C(3,5)*F(0) + F(4)*C(3,1)*F(5) + F(4)*C(3,3)*F(4) + S(1); 
        answer(5,6) = F(0)*C(5,5)*F(7) + F(0)*C(5,1)*F(6) + F(0)*C(5,3)*F(2) + F(5)*C(1,5)*F(7) + F(5)*C(1,1)*F(6) + F(5)*C(1,3)*F(2) + F(4)*C(3,5)*F(7) + F(4)*C(3,1)*F(6) + F(4)*C(3,3)*F(2) + 0.0; 
        answer(5,7) = F(0)*C(5,0)*F(7) + F(0)*C(5,5)*F(6) + F(0)*C(5,4)*F(2) + F(5)*C(1,0)*F(7) + F(5)*C(1,5)*F(6) + F(5)*C(1,4)*F(2) + F(4)*C(3,0)*F(7) + F(4)*C(3,5)*F(6) + F(4)*C(3,4)*F(2) + 0.0; 
        answer(5,8) = F(0)*C(5,0)*F(8) + F(0)*C(5,5)*F(1) + F(0)*C(5,4)*F(3) + F(5)*C(1,0)*F(8) + F(5)*C(1,5)*F(1) + F(5)*C(1,4)*F(3) + F(4)*C(3,0)*F(8) + F(4)*C(3,5)*F(1) + F(4)*C(3,4)*F(3) + 0.0; 
        answer(6,0) = F(7)*C(5,0)*F(0) + F(7)*C(5,5)*F(5) + F(7)*C(5,4)*F(4) + F(6)*C(1,0)*F(0) + F(6)*C(1,5)*F(5) + F(6)*C(1,4)*F(4) + F(2)*C(3,0)*F(0) + F(2)*C(3,5)*F(5) + F(2)*C(3,4)*F(4) + 0.0; 
        answer(6,1) = F(7)*C(5,5)*F(8) + F(7)*C(5,1)*F(1) + F(7)*C(5,3)*F(3) + F(6)*C(1,5)*F(8) + F(6)*C(1,1)*F(1) + F(6)*C(1,3)*F(3) + F(2)*C(3,5)*F(8) + F(2)*C(3,1)*F(1) + F(2)*C(3,3)*F(3) + 0.0; 
        answer(6,2) = F(7)*C(5,4)*F(7) + F(7)*C(5,3)*F(6) + F(7)*C(5,2)*F(2) + F(6)*C(1,4)*F(7) + F(6)*C(1,3)*F(6) + F(6)*C(1,2)*F(2) + F(2)*C(3,4)*F(7) + F(2)*C(3,3)*F(6) + F(2)*C(3,2)*F(2) + S(3); 
        answer(6,3) = F(7)*C(5,4)*F(8) + F(7)*C(5,3)*F(1) + F(7)*C(5,2)*F(3) + F(6)*C(1,4)*F(8) + F(6)*C(1,3)*F(1) + F(6)*C(1,2)*F(3) + F(2)*C(3,4)*F(8) + F(2)*C(3,3)*F(1) + F(2)*C(3,2)*F(3) + 0.0; 
        answer(6,4) = F(7)*C(5,4)*F(0) + F(7)*C(5,3)*F(5) + F(7)*C(5,2)*F(4) + F(6)*C(1,4)*F(0) + F(6)*C(1,3)*F(5) + F(6)*C(1,2)*F(4) + F(2)*C(3,4)*F(0) + F(2)*C(3,3)*F(5) + F(2)*C(3,2)*F(4) + 0.0; 
        answer(6,5) = F(7)*C(5,5)*F(0) + F(7)*C(5,1)*F(5) + F(7)*C(5,3)*F(4) + F(6)*C(1,5)*F(0) + F(6)*C(1,1)*F(5) + F(6)*C(1,3)*F(4) + F(2)*C(3,5)*F(0) + F(2)*C(3,1)*F(5) + F(2)*C(3,3)*F(4) + 0.0; 
        answer(6,6) = F(7)*C(5,5)*F(7) + F(7)*C(5,1)*F(6) + F(7)*C(5,3)*F(2) + F(6)*C(1,5)*F(7) + F(6)*C(1,1)*F(6) + F(6)*C(1,3)*F(2) + F(2)*C(3,5)*F(7) + F(2)*C(3,1)*F(6) + F(2)*C(3,3)*F(2) + S(1); 
        answer(6,7) = F(7)*C(5,0)*F(7) + F(7)*C(5,5)*F(6) + F(7)*C(5,4)*F(2) + F(6)*C(1,0)*F(7) + F(6)*C(1,5)*F(6) + F(6)*C(1,4)*F(2) + F(2)*C(3,0)*F(7) + F(2)*C(3,5)*F(6) + F(2)*C(3,4)*F(2) + S(5); 
        answer(6,8) = F(7)*C(5,0)*F(8) + F(7)*C(5,5)*F(1) + F(7)*C(5,4)*F(3) + F(6)*C(1,0)*F(8) + F(6)*C(1,5)*F(1) + F(6)*C(1,4)*F(3) + F(2)*C(3,0)*F(8) + F(2)*C(3,5)*F(1) + F(2)*C(3,4)*F(3) + 0.0; 
        answer(7,0) = F(7)*C(0,0)*F(0) + F(7)*C(0,5)*F(5) + F(7)*C(0,4)*F(4) + F(6)*C(5,0)*F(0) + F(6)*C(5,5)*F(5) + F(6)*C(5,4)*F(4) + F(2)*C(4,0)*F(0) + F(2)*C(4,5)*F(5) + F(2)*C(4,4)*F(4) + 0.0; 
        answer(7,1) = F(7)*C(0,5)*F(8) + F(7)*C(0,1)*F(1) + F(7)*C(0,3)*F(3) + F(6)*C(5,5)*F(8) + F(6)*C(5,1)*F(1) + F(6)*C(5,3)*F(3) + F(2)*C(4,5)*F(8) + F(2)*C(4,1)*F(1) + F(2)*C(4,3)*F(3) + 0.0; 
        answer(7,2) = F(7)*C(0,4)*F(7) + F(7)*C(0,3)*F(6) + F(7)*C(0,2)*F(2) + F(6)*C(5,4)*F(7) + F(6)*C(5,3)*F(6) + F(6)*C(5,2)*F(2) + F(2)*C(4,4)*F(7) + F(2)*C(4,3)*F(6) + F(2)*C(4,2)*F(2) + S(4); 
        answer(7,3) = F(7)*C(0,4)*F(8) + F(7)*C(0,3)*F(1) + F(7)*C(0,2)*F(3) + F(6)*C(5,4)*F(8) + F(6)*C(5,3)*F(1) + F(6)*C(5,2)*F(3) + F(2)*C(4,4)*F(8) + F(2)*C(4,3)*F(1) + F(2)*C(4,2)*F(3) + 0.0; 
        answer(7,4) = F(7)*C(0,4)*F(0) + F(7)*C(0,3)*F(5) + F(7)*C(0,2)*F(4) + F(6)*C(5,4)*F(0) + F(6)*C(5,3)*F(5) + F(6)*C(5,2)*F(4) + F(2)*C(4,4)*F(0) + F(2)*C(4,3)*F(5) + F(2)*C(4,2)*F(4) + 0.0; 
        answer(7,5) = F(7)*C(0,5)*F(0) + F(7)*C(0,1)*F(5) + F(7)*C(0,3)*F(4) + F(6)*C(5,5)*F(0) + F(6)*C(5,1)*F(5) + F(6)*C(5,3)*F(4) + F(2)*C(4,5)*F(0) + F(2)*C(4,1)*F(5) + F(2)*C(4,3)*F(4) + 0.0; 
        answer(7,6) = F(7)*C(0,5)*F(7) + F(7)*C(0,1)*F(6) + F(7)*C(0,3)*F(2) + F(6)*C(5,5)*F(7) + F(6)*C(5,1)*F(6) + F(6)*C(5,3)*F(2) + F(2)*C(4,5)*F(7) + F(2)*C(4,1)*F(6) + F(2)*C(4,3)*F(2) + S(5); 
        answer(7,7) = F(7)*C(0,0)*F(7) + F(7)*C(0,5)*F(6) + F(7)*C(0,4)*F(2) + F(6)*C(5,0)*F(7) + F(6)*C(5,5)*F(6) + F(6)*C(5,4)*F(2) + F(2)*C(4,0)*F(7) + F(2)*C(4,5)*F(6) + F(2)*C(4,4)*F(2) + S(0); 
        answer(7,8) = F(7)*C(0,0)*F(8) + F(7)*C(0,5)*F(1) + F(7)*C(0,4)*F(3) + F(6)*C(5,0)*F(8) + F(6)*C(5,5)*F(1) + F(6)*C(5,4)*F(3) + F(2)*C(4,0)*F(8) + F(2)*C(4,5)*F(1) + F(2)*C(4,4)*F(3) + 0.0; 
        answer(8,0) = F(8)*C(0,0)*F(0) + F(8)*C(0,5)*F(5) + F(8)*C(0,4)*F(4) + F(1)*C(5,0)*F(0) + F(1)*C(5,5)*F(5) + F(1)*C(5,4)*F(4) + F(3)*C(4,0)*F(0) + F(3)*C(4,5)*F(5) + F(3)*C(4,4)*F(4) + 0.0; 
        answer(8,1) = F(8)*C(0,5)*F(8) + F(8)*C(0,1)*F(1) + F(8)*C(0,3)*F(3) + F(1)*C(5,5)*F(8) + F(1)*C(5,1)*F(1) + F(1)*C(5,3)*F(3) + F(3)*C(4,5)*F(8) + F(3)*C(4,1)*F(1) + F(3)*C(4,3)*F(3) + S(5); 
        answer(8,2) = F(8)*C(0,4)*F(7) + F(8)*C(0,3)*F(6) + F(8)*C(0,2)*F(2) + F(1)*C(5,4)*F(7) + F(1)*C(5,3)*F(6) + F(1)*C(5,2)*F(2) + F(3)*C(4,4)*F(7) + F(3)*C(4,3)*F(6) + F(3)*C(4,2)*F(2) + 0.0; 
        answer(8,3) = F(8)*C(0,4)*F(8) + F(8)*C(0,3)*F(1) + F(8)*C(0,2)*F(3) + F(1)*C(5,4)*F(8) + F(1)*C(5,3)*F(1) + F(1)*C(5,2)*F(3) + F(3)*C(4,4)*F(8) + F(3)*C(4,3)*F(1) + F(3)*C(4,2)*F(3) + S(4); 
        answer(8,4) = F(8)*C(0,4)*F(0) + F(8)*C(0,3)*F(5) + F(8)*C(0,2)*F(4) + F(1)*C(5,4)*F(0) + F(1)*C(5,3)*F(5) + F(1)*C(5,2)*F(4) + F(3)*C(4,4)*F(0) + F(3)*C(4,3)*F(5) + F(3)*C(4,2)*F(4) + 0.0; 
        answer(8,5) = F(8)*C(0,5)*F(0) + F(8)*C(0,1)*F(5) + F(8)*C(0,3)*F(4) + F(1)*C(5,5)*F(0) + F(1)*C(5,1)*F(5) + F(1)*C(5,3)*F(4) + F(3)*C(4,5)*F(0) + F(3)*C(4,1)*F(5) + F(3)*C(4,3)*F(4) + 0.0; 
        answer(8,6) = F(8)*C(0,5)*F(7) + F(8)*C(0,1)*F(6) + F(8)*C(0,3)*F(2) + F(1)*C(5,5)*F(7) + F(1)*C(5,1)*F(6) + F(1)*C(5,3)*F(2) + F(3)*C(4,5)*F(7) + F(3)*C(4,1)*F(6) + F(3)*C(4,3)*F(2) + 0.0; 
        answer(8,7) = F(8)*C(0,0)*F(7) + F(8)*C(0,5)*F(6) + F(8)*C(0,4)*F(2) + F(1)*C(5,0)*F(7) + F(1)*C(5,5)*F(6) + F(1)*C(5,4)*F(2) + F(3)*C(4,0)*F(7) + F(3)*C(4,5)*F(6) + F(3)*C(4,4)*F(2) + 0.0; 
        answer(8,8) = F(8)*C(0,0)*F(8) + F(8)*C(0,5)*F(1) + F(8)*C(0,4)*F(3) + F(1)*C(5,0)*F(8) + F(1)*C(5,5)*F(1) + F(1)*C(5,4)*F(3) + F(3)*C(4,0)*F(8) + F(3)*C(4,5)*F(1) + F(3)*C(4,4)*F(3) + S(0);


    } else if ( matMode == _PlaneStress ) {
        // Save terms associated with H = [du/dx dv/dy du/dy dv/dx]

        answer.resize(4,4);
        answer(0,0) = F(0)*C(0,0)*F(0) + F(0)*C(0,2)*F(2) + F(2)*C(2,0)*F(0) + F(2)*C(2,2)*F(2) + S(0); 
        answer(0,1) = F(0)*C(0,2)*F(3) + F(0)*C(0,1)*F(1) + F(2)*C(2,2)*F(3) + F(2)*C(2,1)*F(1) + 0.0; 
        answer(0,2) = F(0)*C(0,2)*F(0) + F(0)*C(0,1)*F(2) + F(2)*C(2,2)*F(0) + F(2)*C(2,1)*F(2) + S(2); 
        answer(0,3) = F(0)*C(0,0)*F(3) + F(0)*C(0,2)*F(1) + F(2)*C(2,0)*F(3) + F(2)*C(2,2)*F(1) + 0.0; 
        answer(1,0) = F(3)*C(2,0)*F(0) + F(3)*C(2,2)*F(2) + F(1)*C(1,0)*F(0) + F(1)*C(1,2)*F(2) + 0.0; 
        answer(1,1) = F(3)*C(2,2)*F(3) + F(3)*C(2,1)*F(1) + F(1)*C(1,2)*F(3) + F(1)*C(1,1)*F(1) + S(1); 
        answer(1,2) = F(3)*C(2,2)*F(0) + F(3)*C(2,1)*F(2) + F(1)*C(1,2)*F(0) + F(1)*C(1,1)*F(2) + 0.0; 
        answer(1,3) = F(3)*C(2,0)*F(3) + F(3)*C(2,2)*F(1) + F(1)*C(1,0)*F(3) + F(1)*C(1,2)*F(1) + S(2); 
        answer(2,0) = F(0)*C(2,0)*F(0) + F(0)*C(2,2)*F(2) + F(2)*C(1,0)*F(0) + F(2)*C(1,2)*F(2) + S(2); 
        answer(2,1) = F(0)*C(2,2)*F(3) + F(0)*C(2,1)*F(1) + F(2)*C(1,2)*F(3) + F(2)*C(1,1)*F(1) + 0.0; 
        answer(2,2) = F(0)*C(2,2)*F(0) + F(0)*C(2,1)*F(2) + F(2)*C(1,2)*F(0) + F(2)*C(1,1)*F(2) + S(1); 
        answer(2,3) = F(0)*C(2,0)*F(3) + F(0)*C(2,2)*F(1) + F(2)*C(1,0)*F(3) + F(2)*C(1,2)*F(1) + 0.0; 
        answer(3,0) = F(3)*C(0,0)*F(0) + F(3)*C(0,2)*F(2) + F(1)*C(2,0)*F(0) + F(1)*C(2,2)*F(2) + 0.0; 
        answer(3,1) = F(3)*C(0,2)*F(3) + F(3)*C(0,1)*F(1) + F(1)*C(2,2)*F(3) + F(1)*C(2,1)*F(1) + S(2); 
        answer(3,2) = F(3)*C(0,2)*F(0) + F(3)*C(0,1)*F(2) + F(1)*C(2,2)*F(0) + F(1)*C(2,1)*F(2) + 0.0; 
        answer(3,3) = F(3)*C(0,0)*F(3) + F(3)*C(0,2)*F(1) + F(1)*C(2,0)*F(3) + F(1)*C(2,2)*F(1) + S(0); 

    } else if ( matMode == _PlaneStrain ) {
        //Save terms associated with H = [du/dx, dv/dy, dw/dz, du/dy, dv/dx] //@todo not fully checked

        answer.resize(5,5);
        answer(0,0) = F(0)*C(0,0)*F(0) + F(0)*C(0,3)*F(3) + F(3)*C(3,0)*F(0) + F(3)*C(3,3)*F(3) + S(0); 
        answer(0,1) = F(0)*C(0,3)*F(4) + F(0)*C(0,1)*F(1) + F(3)*C(3,3)*F(4) + F(3)*C(3,1)*F(1) + 0.0; 
        answer(0,2) = F(0)*C(0,2)*F(2) + F(3)*C(3,2)*F(2) + 0.0; 
        answer(0,3) = F(0)*C(0,3)*F(0) + F(0)*C(0,1)*F(3) + F(3)*C(3,3)*F(0) + F(3)*C(3,1)*F(3) + S(3); 
        answer(0,4) = F(0)*C(0,0)*F(4) + F(0)*C(0,3)*F(1) + F(3)*C(3,0)*F(4) + F(3)*C(3,3)*F(1) + 0.0; 
        answer(1,0) = F(4)*C(3,0)*F(0) + F(4)*C(3,3)*F(3) + F(1)*C(1,0)*F(0) + F(1)*C(1,3)*F(3) + 0.0; 
        answer(1,1) = F(4)*C(3,3)*F(4) + F(4)*C(3,1)*F(1) + F(1)*C(1,3)*F(4) + F(1)*C(1,1)*F(1) + S(1); 
        answer(1,2) = F(4)*C(3,2)*F(2) + F(1)*C(1,2)*F(2) + 0.0; 
        answer(1,3) = F(4)*C(3,3)*F(0) + F(4)*C(3,1)*F(3) + F(1)*C(1,3)*F(0) + F(1)*C(1,1)*F(3) + 0.0; 
        answer(1,4) = F(4)*C(3,0)*F(4) + F(4)*C(3,3)*F(1) + F(1)*C(1,0)*F(4) + F(1)*C(1,3)*F(1) + S(3); 
        answer(2,0) = F(2)*C(2,0)*F(0) + F(2)*C(2,3)*F(3) + 0.0; 
        answer(2,1) = F(2)*C(2,3)*F(4) + F(2)*C(2,1)*F(1) + 0.0; 
        answer(2,2) = F(2)*C(2,2)*F(2) + S(2); 
        answer(2,3) = F(2)*C(2,3)*F(0) + F(2)*C(2,1)*F(3) + 0.0; 
        answer(2,4) = F(2)*C(2,0)*F(4) + F(2)*C(2,3)*F(1) + 0.0; 
        answer(3,0) = F(0)*C(3,0)*F(0) + F(0)*C(3,3)*F(3) + F(3)*C(1,0)*F(0) + F(3)*C(1,3)*F(3) + S(3); 
        answer(3,1) = F(0)*C(3,3)*F(4) + F(0)*C(3,1)*F(1) + F(3)*C(1,3)*F(4) + F(3)*C(1,1)*F(1) + 0.0; 
        answer(3,2) = F(0)*C(3,2)*F(2) + F(3)*C(1,2)*F(2) + 0.0; 
        answer(3,3) = F(0)*C(3,3)*F(0) + F(0)*C(3,1)*F(3) + F(3)*C(1,3)*F(0) + F(3)*C(1,1)*F(3) + S(1); 
        answer(3,4) = F(0)*C(3,0)*F(4) + F(0)*C(3,3)*F(1) + F(3)*C(1,0)*F(4) + F(3)*C(1,3)*F(1) + 0.0; 
        answer(4,0) = F(4)*C(0,0)*F(0) + F(4)*C(0,3)*F(3) + F(1)*C(3,0)*F(0) + F(1)*C(3,3)*F(3) + 0.0; 
        answer(4,1) = F(4)*C(0,3)*F(4) + F(4)*C(0,1)*F(1) + F(1)*C(3,3)*F(4) + F(1)*C(3,1)*F(1) + S(3); 
        answer(4,2) = F(4)*C(0,2)*F(2) + F(1)*C(3,2)*F(2) + 0.0; 
        answer(4,3) = F(4)*C(0,3)*F(0) + F(4)*C(0,1)*F(3) + F(1)*C(3,3)*F(0) + F(1)*C(3,1)*F(3) + 0.0; 
        answer(4,4) = F(4)*C(0,0)*F(4) + F(4)*C(0,3)*F(1) + F(1)*C(3,0)*F(4) + F(1)*C(3,3)*F(1) + S(0);    

    } else if ( matMode == _1dMat ) {
        //Save terms associated with H = [du/dx] 
        /// @todo is this really correct??

        answer.resize(1,1);
        answer(0,0) = F(0)*C(0,0)*F(0) + S(0);
    }
}

void 
StructuralMaterial :: give_dPdF_from(const FloatMatrix &dSdE, FloatMatrix &answer, GaussPoint *gp)
{
    // Default implementation for converting dSdE to dPdF. This includes updating the 
    // state variables of P and F. 
    
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    FloatArray reducedvF, reducedvP, reducedvS, vP;
    reducedvF = status->giveTempFVector();
    reducedvP = status->giveTempPVector();   
    
    MaterialMode matMode = gp->giveMaterialMode();
    this->convert_P_2_S(reducedvS, reducedvP, reducedvF, matMode);
    this->convert_dSdE_2_dPdF(answer, dSdE, reducedvS, reducedvF, matMode);
}


void
StructuralMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                               MatResponseForm form, MatResponseMode rMode,
                                               GaussPoint *gp, TimeStep *atTime)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMat:
    case _3dMatGrad:
        this->give3dMaterialStiffnessMatrix(answer, rMode, gp, atTime);
        break;
    case _PlaneStress:
    case _PlaneStressGrad:
        this->givePlaneStressStiffMtrx(answer, form, rMode, gp, atTime);
        break;
    case _PlaneStrain:
    case _PlaneStrainGrad:
        this->givePlaneStrainStiffMtrx(answer, form, rMode, gp, atTime);
        break;
    case _1dMat:
    case _1dMatGrad:
        this->give1dStressStiffMtrx(answer, form, rMode, gp, atTime);
        break;
    case _2dPlateLayer:
        this->give2dPlateLayerStiffMtrx(answer, form, rMode, gp, atTime);
        break;
    case _3dShellLayer:
        this->give3dShellLayerStiffMtrx(answer, form, rMode, gp, atTime);
        break;
    case _2dBeamLayer:
        this->give2dBeamLayerStiffMtrx(answer, form, rMode, gp, atTime);
        break;
    case _1dFiber:
        this->give1dFiberStiffMtrx(answer, form, rMode, gp, atTime);
        break;
    default:
        OOFEM_ERROR2( "StructuralMaterial :: giveCharacteristicMatrix : unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}


void 
StructuralMaterial :: giveStiffnessMatrix_dPdF(FloatMatrix &answer, MatResponseForm form, MatResponseMode rMode,
                                               GaussPoint *gp, TimeStep *tStep)
{
    // Returns the stiffness matrix dPdF of the reciever according to MatResponseMode 
    // and MaterialMode

    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMat:
    case _3dMatGrad:
        this->give3dMaterialStiffnessMatrix_dPdF(answer, rMode, gp, tStep);
        break;
    case _PlaneStress:
    case _PlaneStressGrad:
        this->givePlaneStressStiffMtrx_dPdF(answer, form, rMode, gp, tStep);
        break;
    case _PlaneStrain:
    case _PlaneStrainGrad:
        this->givePlaneStrainStiffMtrx_dPdF(answer, form, rMode, gp, tStep);
        break;
    case _1dMat:
    case _1dMatGrad:
        this->give1dStressStiffMtrx_dPdF(answer, form, rMode, gp, tStep);
        break;
    case _2dPlateLayer:
        this->give2dPlateLayerStiffMtrx_dPdF(answer, form, rMode, gp, tStep);
        break;
    case _3dShellLayer:
        this->give3dShellLayerStiffMtrx_dPdF(answer, form, rMode, gp, tStep);
        break;
    case _2dBeamLayer:
        this->give2dBeamLayerStiffMtrx_dPdF(answer, form, rMode, gp, tStep);
        break;
    case _1dFiber:
        this->give1dFiberStiffMtrx_dPdF(answer, form, rMode, gp, tStep);
        break;
    default:
        OOFEM_ERROR2( "StructuralMaterial :: giveCharacteristicMatrix : unknown mode (%s)", __MaterialModeToString(mMode) );
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
                                       MatResponseForm form, MatResponseMode mode,
                                       GaussPoint *gp, TimeStep *tStep)
{ 
    FloatMatrix dSdE;
    this->givePlaneStressStiffMtrx(dSdE, form, mode, gp, tStep); 
    this->give_dPdF_from(dSdE, answer, gp);
}


void
StructuralMaterial :: givePlaneStrainStiffMtrx_dPdF(FloatMatrix &answer,
                                       MatResponseForm form, MatResponseMode mode,
                                       GaussPoint *gp, TimeStep *tStep)
{ 
    FloatMatrix dSdE;
    this->givePlaneStrainStiffMtrx(dSdE, form, mode, gp, tStep); 
    this->give_dPdF_from(dSdE, answer, gp);
}


void
StructuralMaterial :: give1dStressStiffMtrx_dPdF(FloatMatrix &answer,
                                       MatResponseForm form, MatResponseMode mode,
                                       GaussPoint *gp, TimeStep *tStep)
{ 
    FloatMatrix dSdE;
    this->give1dStressStiffMtrx(dSdE, form, mode, gp, tStep); 
    this->give_dPdF_from(dSdE, answer, gp);
}

void
StructuralMaterial :: give2dPlateLayerStiffMtrx_dPdF(FloatMatrix &answer,
                                       MatResponseForm form, MatResponseMode mode,
                                       GaussPoint *gp, TimeStep *tStep)
{ 
    FloatMatrix dSdE;
    this->give2dPlateLayerStiffMtrx(dSdE, form, mode, gp, tStep); 
    this->give_dPdF_from(dSdE, answer, gp);
}


void
StructuralMaterial :: give3dShellLayerStiffMtrx_dPdF(FloatMatrix &answer,
                                       MatResponseForm form, MatResponseMode mode,
                                       GaussPoint *gp, TimeStep *tStep)
{ 
    FloatMatrix dSdE;
    this->give3dShellLayerStiffMtrx(dSdE, form, mode, gp, tStep); 
    this->give_dPdF_from(dSdE, answer, gp);
}


void
StructuralMaterial :: give2dBeamLayerStiffMtrx_dPdF(FloatMatrix &answer,
                                       MatResponseForm form, MatResponseMode mode,
                                       GaussPoint *gp, TimeStep *tStep)
{ 
    FloatMatrix dSdE;
    this->give2dBeamLayerStiffMtrx(dSdE, form, mode, gp, tStep); 
    this->give_dPdF_from(dSdE, answer, gp);
}


void
StructuralMaterial :: give1dFiberStiffMtrx_dPdF(FloatMatrix &answer,
                                       MatResponseForm form, MatResponseMode mode,
                                       GaussPoint *gp, TimeStep *tStep)
{ 
    FloatMatrix dSdE;
    this->give1dFiberStiffMtrx(dSdE, form, mode, gp, tStep); 
    this->give_dPdF_from(dSdE, answer, gp);
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
    S.beProductOf(invF,P);
    FloatArray vS;
    vS.beReducedVectorForm(S); // 6 components 
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
    FloatMatrix F, P, S, invF;
    F.beMatrixForm(vF);
    S.beMatrixForm(vS);
    P.beProductOf(F,S);
    vP.beFullVectorForm(P); 
    StructuralMaterial :: giveReducedVectorForm(answer, vP, matMode);   // convert back to reduced size
}


void
StructuralMaterial :: giveIdentityVector(FloatArray &answer, MaterialMode matMode)
{
    // Create identity tensor on Voigt form according to MaterialMode.

    ///@todo This is a hack. We should just store 3D tensors always and be done with it.
    int size;
    if ( matMode == _3dBeam ) {
        size = 9;
    } else {
        size = StructuralMaterial :: giveSizeOfVoigtVector( matMode );
    }
    answer.resize(size);

    if ( size == 9 || size == 5 ) {
        answer.at(1) = answer.at(2) = answer.at(3) = 1.0;
    } else if ( size == 4 ) {
        answer.at(1) = answer.at(2) = 1.0;
    } else if (size == 1 ) {
        answer.at(1) = 1.0;
    } else {
        answer.resize(0);
    }
}

void
StructuralMaterial :: giveCharacteristicComplianceMatrix(FloatMatrix &answer,
                                                         MatResponseForm form, MatResponseMode rMode,
                                                         GaussPoint *gp,
                                                         TimeStep *atTime)
//
// Returns characteristic material compliance matrix of the receiver
// works for positive definite associated stiffnesses only
//
{
    FloatMatrix redInvAnswer, redAnswer;
    IntArray mask;

    this->giveCharacteristicMatrix(redInvAnswer, ReducedForm, rMode, gp, atTime);
    redAnswer.beInverseOf(redInvAnswer);

    if ( form == FullForm ) {
        int size = StructuralMaterial :: giveVoigtSymVectorMask(mask, gp->giveMaterialMode());
        answer.resize(size, size);
        answer.zero();
        answer.assemble(redAnswer,mask,mask);
    } else if ( form == ReducedForm ) {
        answer = redAnswer;
    } else {
        OOFEM_ERROR("StructuralMaterial :: giveCharacteristicComplianceMatrix - unsupported form mode");
    }
}


void
StructuralMaterial ::  reduceStiffMtrx3d(FloatMatrix &answer, MatResponseForm form, GaussPoint *gp,
                                         FloatMatrix &stiffMtrx3d) const
//
// Returns characteristic material stiffness matrix of the receiver
// reduced to corresponding mode obtained from gp.
{
    MaterialMode mode = gp->giveMaterialMode();
    switch ( mode ) {
    case _3dMat:
        answer = stiffMtrx3d;
        break;
    case _PlaneStress:
        this->reduceToPlaneStressStiffMtrx(answer, form, gp, stiffMtrx3d);
        break;
    case _PlaneStrain:
        this->reduceToPlaneStrainStiffMtrx(answer, form, gp, stiffMtrx3d);
        break;
    case _1dMat:
        this->reduceTo1dStressStiffMtrx(answer, form, gp, stiffMtrx3d);
        break;
    case _2dPlateLayer:
        this->reduceTo2dPlateLayerStiffMtrx(answer, form, gp, stiffMtrx3d);
        break;
    case _3dShellLayer:
        this->reduceTo3dShellLayerStiffMtrx(answer, form, gp, stiffMtrx3d);
        break;
    case _2dBeamLayer:
        this->reduceTo2dBeamLayerStiffMtrx(answer, form, gp, stiffMtrx3d);
        break;
    case _1dFiber:
        this->reduceTo1dFiberStiffMtrx(answer, form, gp, stiffMtrx3d);
        break;
    default:
        OOFEM_ERROR2( "StructuralMaterial :: reduceStiffMtrx3d : unknown mode (%s)", __MaterialModeToString(mode) );
    }
}


void
StructuralMaterial :: reduceComplMtrx3d(FloatMatrix &answer, MatResponseForm form, GaussPoint *gp,
                                        FloatMatrix &complMtrx3d) const
//
// Returns characteristic material compliance matrix of the receiver
// reduced to corresponding mode obtained from gp.
{
    MaterialMode mode = gp->giveMaterialMode();
    switch ( mode ) {
    case _3dMat:
        answer = complMtrx3d;
        break;
    case _PlaneStress:
        this->reduceToPlaneStressComplMtrx(answer, form, gp, complMtrx3d);
        break;
    case _PlaneStrain:
        this->reduceToPlaneStrainComplMtrx(answer, form, gp, complMtrx3d);
        break;
    case _1dMat:
        this->reduceTo1dStressComplMtrx(answer, form, gp, complMtrx3d);
        break;
    case _2dPlateLayer:
        this->reduceTo2dPlateLayerComplMtrx(answer, form, gp, complMtrx3d);
        break;
    case _3dShellLayer:
        this->reduceTo3dShellLayerComplMtrx(answer, form, gp, complMtrx3d);
        break;
    case _2dBeamLayer:
        this->reduceTo2dBeamLayerComplMtrx(answer, form, gp, complMtrx3d);
        break;
    case _1dFiber:
        this->reduceTo1dFiberComplMtrx(answer, form, gp, complMtrx3d);
        break;
    default:
        OOFEM_ERROR2( "StructuralMaterial :: reduceComplMtrx3d : unknown mode (%s)", __MaterialModeToString(mode) );
    }
}


void
StructuralMaterial :: giveStressDependentPartOfStrainVector(FloatArray &answer, GaussPoint *gp,
                                                            const FloatArray &reducedStrainVector,
                                                            TimeStep *stepN, ValueModeType mode)
{
    /*
     * This functions subtract from reducedStrainVector its stress independent part
     * caused by temperature, shrinkage and possibly by other phenomena.
     */
    FloatArray epsilonTemperature;
    StructuralCrossSection *cs = static_cast< StructuralCrossSection * >( gp->giveElement()->giveCrossSection() );

    answer = reducedStrainVector;
    cs->computeStressIndependentStrainVector(epsilonTemperature, gp, stepN, mode);
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
StructuralMaterial :: giveStressStrainMask(IntArray &answer, MatResponseForm form,
                                           MaterialMode mmode)
//
// this function returns mask of reduced(if form == ReducedForm)
// or Full(if form==FullForm) stressStrain vector in full or
// reduced StressStrainVector
// according to stressStrain mode of given gp.
//
//
// mask has size of reduced or full StressStrain Vector and  i-th component
// is index to full or reduced StressStrainVector where corresponding
// stressStrain resides.
//
// Reduced form is sub-vector (of stress or strain components),
// where components corresponding to imposed zero stress (plane stress,...)
// are not included. On the other hand, if zero strain component is imposed
// (Plane strain, ..) this condition must be taken into account in geometrical
// relations, and corresponding component is included in reduced vector.
//
{
    if ( form == ReducedForm ) {
        StructuralMaterial :: giveVoigtSymVectorMask(answer, mmode);
    } else if ( form == FullForm ) {
        IntArray mask;
        answer.resize( StructuralMaterial :: giveVoigtSymVectorMask(mask, mmode) );
        answer.zero();
        for ( int i = 1; i <= mask.giveSize(); i++ ) {
            answer.at(mask.at(i)) = i;
        }
    } else {
        OOFEM_ERROR("StructuralMaterial :: giveStressStrainMask : unknown form mode");
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
    case _2dPlateLayer:
    case _3dShellLayer:
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
    case _1dFiber:
        answer.resize(3);
        answer.at(1) = 1;
        answer.at(2) = 5;
        answer.at(3) = 6;
        return 6;
    case _3dMatGrad:
        answer.resize(7);
        for ( int i = 1; i <= 7; i++ ) {
            answer.at(i) = i;
        }
        return 7;
    case _PlaneStressGrad:
        answer.resize(4);
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 6;
        answer.at(4) = 7;
        return 7;
    case _PlaneStrainGrad:
        answer.resize(5);
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 3;
        answer.at(4) = 6;
        answer.at(5) = 7;
        return 7;
    case _1dMatGrad:
        answer.resize(2);
        answer.at(1) = 1;
        answer.at(2) = 7;
        return 7;
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
        answer.resize(0);
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

// Stiffness reduction methods
#if 1
void
StructuralMaterial :: reduceToPlaneStressStiffMtrx(FloatMatrix &answer,
                                                   MatResponseForm form, GaussPoint *gp,
                                                   FloatMatrix &stiffMtrx3d) const
//
// returns receiver's 2dPlaneStressMtrx constructed from
// stiffMtrx3d (general 3dMatrialStiffnessMatrix)
// (2dPlaneStres ==> sigma_z = tau_xz = tau_yz = 0.)
// This method works for general 3d stiff matrix
//
{
    FloatMatrix inv3d, invAnswer, reducedAnswer;

    // check if stiffMtrx is proper
    if ( ( stiffMtrx3d.isSquare() ) && ( stiffMtrx3d.giveNumberOfRows() == 6 ) ) {
        inv3d.beInverseOf(stiffMtrx3d);

        invAnswer.resize(3, 3);

        invAnswer.at(1, 1) = inv3d.at(1, 1);
        invAnswer.at(1, 2) = inv3d.at(1, 2);
        invAnswer.at(1, 3) = inv3d.at(1, 6);

        invAnswer.at(2, 1) = inv3d.at(2, 1);
        invAnswer.at(2, 2) = inv3d.at(2, 2);
        invAnswer.at(2, 3) = inv3d.at(2, 6);

        invAnswer.at(3, 1) = inv3d.at(6, 1);
        invAnswer.at(3, 2) = inv3d.at(6, 2);
        invAnswer.at(3, 3) = inv3d.at(6, 6);

        reducedAnswer.beInverseOf(invAnswer);

        if ( form == ReducedForm ) {
            answer = reducedAnswer;
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = reducedAnswer.at(1, 1);
            answer.at(1, 2) = reducedAnswer.at(1, 2);
            answer.at(1, 6) = reducedAnswer.at(1, 3);
            answer.at(2, 1) = reducedAnswer.at(2, 1);
            answer.at(2, 2) = reducedAnswer.at(2, 2);
            answer.at(2, 6) = reducedAnswer.at(2, 3);
            answer.at(6, 6) = reducedAnswer.at(3, 3);
            answer.at(6, 1) = reducedAnswer.at(3, 1);
            answer.at(6, 2) = reducedAnswer.at(3, 2);
        }

    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceToPlaneStressStiffMtrx : stiffMtrx size mismatch");
    }
}


void
StructuralMaterial :: reduceToPlaneStrainStiffMtrx(FloatMatrix &answer,
                                                   MatResponseForm form, GaussPoint *gp,
                                                   FloatMatrix &stiffMtrx3d) const
//
// returns receiver's 2dPlaneStrainMtrx constructed from
// general 3dMatrialStiffnessMatrix
// (2dPlaneStrain ==> eps_z = gamma_xz = gamma_yz = 0.)
// but we take ez, SigmaZ into account.
//
{
    // check if stiffMtrx is proper
    if ( ( stiffMtrx3d.isSquare() ) && ( stiffMtrx3d.giveNumberOfRows() == 6 ) ) {
        if ( form == ReducedForm ) {
            answer.resize(4, 4);
            answer.zero();

            answer.at(1, 1) = stiffMtrx3d.at(1, 1);
            answer.at(1, 2) = stiffMtrx3d.at(1, 2);
            answer.at(1, 4) = stiffMtrx3d.at(1, 6);

            answer.at(2, 1) = stiffMtrx3d.at(2, 1);
            answer.at(2, 2) = stiffMtrx3d.at(2, 2);
            answer.at(2, 4) = stiffMtrx3d.at(2, 6);

            answer.at(3, 1) = stiffMtrx3d.at(3, 1);
            answer.at(3, 2) = stiffMtrx3d.at(3, 2);
            answer.at(3, 4) = stiffMtrx3d.at(3, 6);

            answer.at(4, 1) = stiffMtrx3d.at(6, 1);
            answer.at(4, 2) = stiffMtrx3d.at(6, 2);
            answer.at(4, 4) = stiffMtrx3d.at(6, 6);
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = stiffMtrx3d.at(1, 1);
            answer.at(1, 2) = stiffMtrx3d.at(1, 2);
            answer.at(1, 6) = stiffMtrx3d.at(1, 6);

            answer.at(2, 1) = stiffMtrx3d.at(2, 1);
            answer.at(2, 2) = stiffMtrx3d.at(2, 2);
            answer.at(2, 6) = stiffMtrx3d.at(2, 6);

            answer.at(3, 1) = stiffMtrx3d.at(3, 1);
            answer.at(3, 2) = stiffMtrx3d.at(3, 2);
            answer.at(3, 6) = stiffMtrx3d.at(3, 6);

            answer.at(6, 6) = stiffMtrx3d.at(3, 3);
            answer.at(6, 1) = stiffMtrx3d.at(6, 1);
            answer.at(6, 2) = stiffMtrx3d.at(6, 2);
        }

    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceToPlaneStrainStiffMtrx :: stiffMtrx size mismatch");
    }
}


void
StructuralMaterial :: reduceTo1dStressStiffMtrx(FloatMatrix &answer,
                                                MatResponseForm form, GaussPoint *gp,
                                                FloatMatrix &stiffMtrx3d) const
//
//
// returns receiver's 1dMaterialStiffnessMAtrix constructed from
// general 3dMatrialStiffnessMatrix
// (1d case ==> sigma_y = sigma_z = tau_yz = tau_zx = tau_xy  = 0.)
// This method works only if 3dMateriallStiffnessMatrix
// has two 3x3 independent blocks
{
    FloatMatrix m3d11, inv3d;
    double val11;

    // check if stiffMtrx is proper
    if ( ( stiffMtrx3d.isSquare() ) && ( stiffMtrx3d.giveNumberOfRows() == 6 ) ) {
        inv3d.beInverseOf(stiffMtrx3d);
        val11 = inv3d.at(1, 1);

        if ( form == ReducedForm ) {
            answer.resize(1, 1);

            answer.at(1, 1) = 1. / val11;
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = 1. / val11;
        }

    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceTo1dStressStiffMtrx:: stiffMtrx3d size mismatch");
    }
}


void
StructuralMaterial :: reduceTo2dPlateLayerStiffMtrx(FloatMatrix &answer,
                                                    MatResponseForm form,
                                                    GaussPoint *gp,
                                                    FloatMatrix &stiffMtrx3d) const
//
// return material stiffness matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
{
    MaterialMode mode = gp->giveMaterialMode();
    FloatMatrix invMat3d, invMatLayer(5, 5), matLayer;

    if ( !( ( mode == _2dPlateLayer ) || ( mode == _3dShellLayer ) ) ) {
        _error("ReduceTo2dPlateLayerStiffMtrx : unsupported mode");
    }


    // check if stiffMtrx is proper
    if ( ( stiffMtrx3d.isSquare() ) && ( stiffMtrx3d.giveNumberOfRows() == 6 ) ) {
        invMat3d.beInverseOf(stiffMtrx3d);

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

        matLayer.beInverseOf(invMatLayer);

        if ( form == ReducedForm ) {
            answer = matLayer;
        } else {
            answer.resize(6, 6);
            answer.zero();

            for ( int i = 1; i <= 2; i++ ) {
                for ( int j = 1; j <= 2; j++ ) {
                    answer.at(i, j) = matLayer.at(i, j);
                }
            }

            for ( int i = 4; i <= 6; i++ ) {
                for ( int j = 4; j <= 6; j++ ) {
                    answer.at(i, j) = matLayer.at(i - 1, j - 1);
                }
            }

            for ( int i = 1; i <= 2; i++ ) {
                for ( int j = 4; j <= 6; j++ ) {
                    answer.at(i, j) = matLayer.at(i, j - 1);
                    answer.at(j, i) = matLayer.at(j - 1, i);
                }
            }
        }
    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceTo2dPlateLayerStiffMtrx : stiffMtrx size mismatch");
    }
}


void
StructuralMaterial :: reduceTo3dShellLayerStiffMtrx(FloatMatrix &answer,
                                                    MatResponseForm form,
                                                    GaussPoint *gp,
                                                    FloatMatrix &stiffMtrx3d) const
//
// return material stiffness matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}

{
    this->reduceTo2dPlateLayerStiffMtrx(answer, form, gp, stiffMtrx3d);
}


void
StructuralMaterial :: reduceTo2dBeamLayerStiffMtrx(FloatMatrix &answer,
                                                   MatResponseForm form,
                                                   GaussPoint *gp,
                                                   FloatMatrix &stiffMtrx3d) const
//
// return material stiffness matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
{
    MaterialMode mode = gp->giveMaterialMode();
    FloatMatrix invMat3d, invMatLayer(2, 2), matLayer;

    if ( mode != _2dBeamLayer ) {
        OOFEM_ERROR("StructuralMaterial :: ReduceTo2dBeamLayerStiffMtrx : unsupported mode");
    }

    if ( ( stiffMtrx3d.isSquare() ) && ( stiffMtrx3d.giveNumberOfRows() == 6 ) ) {
        invMat3d.beInverseOf(stiffMtrx3d);

        invMatLayer.at(1, 1) = invMat3d.at(1, 1);
        invMatLayer.at(1, 2) = invMat3d.at(1, 5);
        invMatLayer.at(2, 1) = invMat3d.at(5, 1);
        invMatLayer.at(2, 2) = invMat3d.at(5, 5);

        matLayer.beInverseOf(invMatLayer);

        if ( form == ReducedForm ) {
            answer = matLayer;
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = matLayer.at(1, 1);
            answer.at(1, 5) = matLayer.at(1, 2);
            answer.at(5, 1) = matLayer.at(2, 1);
            answer.at(5, 5) = matLayer.at(2, 2);
        }
    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceTo2dBeamLayerStiffMtrx: stiffMtrx3d size mismatch");
    }
}


void
StructuralMaterial :: reduceTo1dFiberStiffMtrx(FloatMatrix &answer,
                                               MatResponseForm form,
                                               GaussPoint *gp,
                                               FloatMatrix &stiffMtrx3d) const
//
// return material stiffness matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
{
    MaterialMode mode = gp->giveMaterialMode();
    FloatMatrix invMat3d, invMatLayer(3, 3), matLayer;

    if ( mode != _1dFiber ) {
        OOFEM_ERROR("StructuralMaterial :: reduceTo1dFiberStiffMtrx : unsupported mode");
    }

    if ( ( stiffMtrx3d.isSquare() ) && ( stiffMtrx3d.giveNumberOfRows() == 6 ) ) {
        invMat3d.beInverseOf(stiffMtrx3d);

        invMatLayer.at(1, 1) = invMat3d.at(1, 1);
        invMatLayer.at(1, 2) = invMat3d.at(1, 5);
        invMatLayer.at(1, 3) = invMat3d.at(1, 6);
        invMatLayer.at(2, 1) = invMat3d.at(5, 1);
        invMatLayer.at(2, 2) = invMat3d.at(5, 5);
        invMatLayer.at(2, 3) = invMat3d.at(5, 6);
        invMatLayer.at(3, 1) = invMat3d.at(6, 1);
        invMatLayer.at(3, 2) = invMat3d.at(6, 5);
        invMatLayer.at(3, 3) = invMat3d.at(6, 6);

        matLayer.beInverseOf(invMatLayer);

        if ( form == ReducedForm ) {
            answer = matLayer;
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = matLayer.at(1, 1);
            answer.at(1, 5) = matLayer.at(1, 2);
            answer.at(1, 6) = matLayer.at(1, 3);
            answer.at(5, 1) = matLayer.at(2, 1);
            answer.at(5, 5) = matLayer.at(2, 2);
            answer.at(5, 6) = matLayer.at(2, 3);
            answer.at(6, 1) = matLayer.at(3, 1);
            answer.at(6, 5) = matLayer.at(3, 2);
            answer.at(6, 6) = matLayer.at(3, 3);
        }
    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceTo1dFiberStiffMtrx: stiffMtrx3d size mismatch");
    }
}
#endif


// Compliance reduction methods
#if 1
void
StructuralMaterial :: reduceToPlaneStressComplMtrx(FloatMatrix &answer,
                                                   MatResponseForm form, GaussPoint *gp,
                                                   FloatMatrix &complMtrx3d) const
//
// returns receiver's 2dPlaneComplMtrx constructed from
// complMtrx3d (general 3dMatrialComplianceMatrix)
// (2dPlaneStres ==> sigma_z = tau_xz = tau_yz = 0.)
// This method works for general 3d compl matrix
//
{
    // check if complMtrx is proper
    if ( ( complMtrx3d.isSquare() ) && ( complMtrx3d.giveNumberOfRows() == 6 ) ) {
        if ( form == ReducedForm ) {
            answer.resize(3, 3);
            answer.zero();

            answer.at(1, 1) = complMtrx3d.at(1, 1);
            answer.at(1, 2) = complMtrx3d.at(1, 2);
            answer.at(1, 3) = complMtrx3d.at(1, 6);
            answer.at(2, 1) = complMtrx3d.at(2, 1);
            answer.at(2, 2) = complMtrx3d.at(2, 2);
            answer.at(2, 3) = complMtrx3d.at(2, 6);
            answer.at(3, 3) = complMtrx3d.at(6, 6);
            answer.at(3, 1) = complMtrx3d.at(6, 1);
            answer.at(3, 2) = complMtrx3d.at(6, 2);
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = complMtrx3d.at(1, 1);
            answer.at(1, 2) = complMtrx3d.at(1, 2);
            answer.at(1, 6) = complMtrx3d.at(1, 6);
            answer.at(2, 1) = complMtrx3d.at(2, 1);
            answer.at(2, 2) = complMtrx3d.at(2, 2);
            answer.at(2, 6) = complMtrx3d.at(2, 6);
            answer.at(6, 6) = complMtrx3d.at(6, 6);
            answer.at(6, 1) = complMtrx3d.at(6, 1);
            answer.at(6, 2) = complMtrx3d.at(6, 2);
        }
    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceToPlaneStressComplMtrx : complMtrx size mismatch");
    }
}


void
StructuralMaterial :: reduceToPlaneStrainComplMtrx(FloatMatrix &answer,
                                                   MatResponseForm form, GaussPoint *gp,
                                                   FloatMatrix &complMtrx3d) const
//
// returns receiver's 2dPlaneStrainMtrx constructed from
// general 3dMatrialComplianceMatrix
// (2dPlaneStrain ==> eps_z = gamma_xz = gamma_yz = 0.)
//
{
    FloatMatrix inv3d, invAnswer(3, 3), reducedAnswer;

    // check if complMtrx is proper
    if ( ( complMtrx3d.isSquare() ) && ( complMtrx3d.giveNumberOfRows() == 6 ) ) {
        inv3d.beInverseOf(complMtrx3d);

        invAnswer.at(1, 1) = inv3d.at(1, 1);
        invAnswer.at(1, 2) = inv3d.at(1, 2);
        invAnswer.at(1, 3) = inv3d.at(1, 6);

        invAnswer.at(2, 1) = inv3d.at(2, 1);
        invAnswer.at(2, 2) = inv3d.at(2, 2);
        invAnswer.at(2, 3) = inv3d.at(2, 6);

        invAnswer.at(3, 1) = inv3d.at(6, 1);
        invAnswer.at(3, 2) = inv3d.at(6, 2);
        invAnswer.at(3, 3) = inv3d.at(6, 6);

        reducedAnswer.beInverseOf(invAnswer);

        if ( form == ReducedForm ) {
            answer.resize(4, 4);
            answer.zero();

            for ( int i = 1; i <= 2; i++ ) {
                for ( int j = 1; j <= 2; j++ ) {
                    answer.at(i, j) = reducedAnswer.at(i, j);
                }
            }

            answer.at(1, 4) = reducedAnswer.at(1, 3);
            answer.at(2, 4) = reducedAnswer.at(2, 3);
            answer.at(4, 1) = reducedAnswer.at(3, 1);
            answer.at(4, 2) = reducedAnswer.at(3, 2);
            answer.at(4, 4) = reducedAnswer.at(3, 3);
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = reducedAnswer.at(1, 1);
            answer.at(1, 2) = reducedAnswer.at(1, 2);
            answer.at(1, 6) = reducedAnswer.at(1, 3);
            answer.at(2, 1) = reducedAnswer.at(2, 1);
            answer.at(2, 2) = reducedAnswer.at(2, 2);
            answer.at(2, 6) = reducedAnswer.at(2, 3);
            answer.at(6, 6) = reducedAnswer.at(3, 3);
            answer.at(6, 1) = reducedAnswer.at(3, 1);
            answer.at(6, 2) = reducedAnswer.at(3, 2);
        }

    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceToPlaneStrainComplMtrx :: complMtrx size mismatch");
    }
}


void
StructuralMaterial :: reduceTo1dStressComplMtrx(FloatMatrix &answer,
                                                MatResponseForm form, GaussPoint *gp,
                                                FloatMatrix &complMtrx3d) const
//
//
// returns receiver's 1dMaterialComplianceMAtrix constructed from
// general 3dMatrialComplianceMatrix
// (1d case ==> sigma_y = sigma_z = tau_yz = tau_zx = tau_xy  = 0.)
{
    // check if complMtrx is proper
    if ( ( complMtrx3d.isSquare() ) && ( complMtrx3d.giveNumberOfRows() == 6 ) ) {
        if ( form == ReducedForm ) {
            answer.resize(1, 1);

            answer.at(1, 1) = complMtrx3d.at(1, 1);
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = complMtrx3d.at(1, 1);
        }

        return;
    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceTo1dStressComplMtrx:: complMtrx3d size mismatch");
    }
}



void
StructuralMaterial :: reduceTo2dPlateLayerComplMtrx(FloatMatrix &answer,
                                                    MatResponseForm form,
                                                    GaussPoint *gp,
                                                    FloatMatrix &complMtrx3d) const
//
// return material compliance matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
{
    MaterialMode mode = gp->giveMaterialMode();

    if ( !( ( mode == _2dPlateLayer ) || ( mode == _3dShellLayer ) ) ) {
        OOFEM_ERROR("StructuralMaterial :: reduceTo2dPlateLayerComplMtrx : unsupported mode");
    }


    // check if complMtrx is proper
    if ( ( complMtrx3d.isSquare() ) && ( complMtrx3d.giveNumberOfRows() == 6 ) ) {
        if ( form == ReducedForm ) {
            answer.resize(5, 5);
            answer.zero();

            for ( int i = 1; i <= 2; i++ ) {
                for ( int j = 1; j <= 2; j++ ) {
                    answer.at(i, j) = complMtrx3d.at(i, j);
                }
            }

            for ( int i = 4; i <= 6; i++ ) {
                for ( int j = 4; j <= 6; j++ ) {
                    answer.at(i - 1, j - 1) = complMtrx3d.at(i, j);
                }
            }

            for ( int i = 1; i <= 2; i++ ) {
                for ( int j = 4; j <= 6; j++ ) {
                    answer.at(i, j - 1) = complMtrx3d.at(i, j);
                    answer.at(j - 1, i) = complMtrx3d.at(j, i);
                }
            }
        } else {
            answer.resize(6, 6);
            answer.zero();

            for ( int i = 1; i <= 2; i++ ) {
                for ( int j = 1; j <= 2; j++ ) {
                    answer.at(i, j) = complMtrx3d.at(i, j);
                }
            }

            for ( int i = 4; i <= 6; i++ ) {
                for ( int j = 4; j <= 6; j++ ) {
                    answer.at(i, j) = complMtrx3d.at(i, j);
                }
            }

            for ( int i = 1; i <= 2; i++ ) {
                for ( int j = 4; j <= 6; j++ ) {
                    answer.at(i, j) = complMtrx3d.at(i, j);
                    answer.at(j, i) = complMtrx3d.at(j, i);
                }
            }
        }

    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceTo2dPlateLayerComplMtrx : stiffMtrx size mismatch");
    }
}


void
StructuralMaterial :: reduceTo3dShellLayerComplMtrx(FloatMatrix &answer,
                                                    MatResponseForm form,
                                                    GaussPoint *gp,
                                                    FloatMatrix &complMtrx3d) const
//
// return material compliance matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}

{
    this->reduceTo2dPlateLayerComplMtrx(answer, form, gp, complMtrx3d);
}



void
StructuralMaterial :: reduceTo2dBeamLayerComplMtrx(FloatMatrix &answer,
                                                   MatResponseForm form,
                                                   GaussPoint *gp,
                                                   FloatMatrix &complMtrx3d) const
//
// return material compliance matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
{
    MaterialMode mode = gp->giveMaterialMode();

    if ( mode != _2dBeamLayer ) {
        OOFEM_ERROR("StructuralMaterial :: reduceTo2dBeamLayerComplMtrx : unsupported mode");
    }

    if ( ( complMtrx3d.isSquare() ) && ( complMtrx3d.giveNumberOfRows() == 6 ) ) {
        if ( form == ReducedForm ) {
            answer.resize(2, 2);

            answer.at(1, 1) = complMtrx3d.at(1, 1);
            answer.at(1, 2) = complMtrx3d.at(1, 5);
            answer.at(2, 1) = complMtrx3d.at(5, 1);
            answer.at(2, 2) = complMtrx3d.at(5, 5);
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = complMtrx3d.at(1, 1);
            answer.at(1, 5) = complMtrx3d.at(1, 5);
            answer.at(5, 1) = complMtrx3d.at(5, 1);
            answer.at(5, 5) = complMtrx3d.at(5, 5);
        }

    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceTo2dBeamLayerStiffMtrx: stiffMtrx3d size mismatch");
    }
}


void
StructuralMaterial :: reduceTo1dFiberComplMtrx(FloatMatrix &answer,
                                               MatResponseForm form,
                                               GaussPoint *gp,
                                               FloatMatrix &complMtrx3d) const
//
// return material compliance matrix for derived types of stressStreinState
// assumption sigma_z = 0.
//
// General strain vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
{
    MaterialMode mode = gp->giveMaterialMode();

    if ( mode != _1dFiber ) {
        OOFEM_ERROR("StructuralMaterial :: reduceTo1dFiberComplMtrx : unsupported mode");
    }

    if ( ( complMtrx3d.isSquare() ) && ( complMtrx3d.giveNumberOfRows() == 6 ) ) {
        if ( form == ReducedForm ) {
            answer.resize(3, 3);

            answer.at(1, 1) = complMtrx3d.at(1, 1);
            answer.at(1, 2) = complMtrx3d.at(1, 5);
            answer.at(1, 3) = complMtrx3d.at(1, 6);
            answer.at(2, 1) = complMtrx3d.at(5, 1);
            answer.at(2, 2) = complMtrx3d.at(5, 5);
            answer.at(2, 3) = complMtrx3d.at(5, 6);
            answer.at(3, 1) = complMtrx3d.at(6, 1);
            answer.at(3, 2) = complMtrx3d.at(6, 5);
            answer.at(3, 3) = complMtrx3d.at(6, 6);
        } else {
            answer.resize(6, 6);
            answer.zero();

            answer.at(1, 1) = complMtrx3d.at(1, 1);
            answer.at(1, 5) = complMtrx3d.at(1, 5);
            answer.at(1, 6) = complMtrx3d.at(1, 6);
            answer.at(5, 1) = complMtrx3d.at(5, 1);
            answer.at(5, 5) = complMtrx3d.at(5, 5);
            answer.at(5, 6) = complMtrx3d.at(5, 6);
            answer.at(6, 1) = complMtrx3d.at(6, 1);
            answer.at(6, 5) = complMtrx3d.at(6, 5);
            answer.at(6, 6) = complMtrx3d.at(6, 6);
        }

    } else {
        OOFEM_ERROR("StructuralMaterial :: reduceTo1dFiberComplMtrx: stiffMtrx3d size mismatch");
    }
}
#endif


void
StructuralMaterial :: givePlaneStressStiffMtrx(FloatMatrix &answer,
                                               MatResponseForm form, MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *atTime)
//
// returns Mat stiffness for PlaneStress
//
{
    FloatMatrix m3d;

    this->give3dMaterialStiffnessMatrix(m3d, mode, gp, atTime);
    this->reduceToPlaneStressStiffMtrx(answer, form, gp, m3d);
}

void
StructuralMaterial :: givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                               MatResponseForm form, MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *atTime)
//
// return material stiffness matrix for PlaneStrain mode
//
{
    FloatMatrix m3d;

    this->give3dMaterialStiffnessMatrix(m3d, mode, gp, atTime);
    this->reduceToPlaneStrainStiffMtrx(answer, form, gp, m3d);
}

void
StructuralMaterial :: give1dStressStiffMtrx(FloatMatrix &answer,
                                            MatResponseForm form, MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *atTime)
//
// return material stiffness matrix for 1d stress strain mode
//
{
    FloatMatrix m3d;

    this->give3dMaterialStiffnessMatrix(m3d, mode, gp, atTime);
    this->reduceTo1dStressStiffMtrx(answer, form, gp, m3d);
}


void
StructuralMaterial :: give2dBeamLayerStiffMtrx(FloatMatrix &answer,
                                               MatResponseForm form, MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *atTime)
//
// return material stiffness matrix for2dBeamLayer mode
//
{
    FloatMatrix m3d;

    this->give3dMaterialStiffnessMatrix(m3d, mode, gp, atTime);
    this->reduceTo2dBeamLayerStiffMtrx(answer, form, gp, m3d);
}


void
StructuralMaterial :: give2dPlateLayerStiffMtrx(FloatMatrix &answer,
                                                MatResponseForm form, MatResponseMode mode,
                                                GaussPoint *gp,
                                                TimeStep *atTime)
//
// return material stiffness matrix for 2dPlateLayer
//
{
    FloatMatrix m3d;

    this->give3dMaterialStiffnessMatrix(m3d, mode, gp, atTime);
    this->reduceTo2dPlateLayerStiffMtrx(answer, form, gp, m3d);
}

void
StructuralMaterial :: give1dFiberStiffMtrx(FloatMatrix &answer,
                                           MatResponseForm form, MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime)
//
// return material stiffness matrix for 2dPlateLayer
//
{
    FloatMatrix m3d;

    this->give3dMaterialStiffnessMatrix(m3d, mode, gp, atTime);
    this->reduceTo1dFiberStiffMtrx(answer, form, gp, m3d);
}


void
StructuralMaterial :: give3dShellLayerStiffMtrx(FloatMatrix &answer,
                                                MatResponseForm form, MatResponseMode mode,
                                                GaussPoint *gp,
                                                TimeStep *atTime)
//
// returns material stiffness matrix for 3dShellLayer
//
{
    FloatMatrix m3d;

    this->give3dMaterialStiffnessMatrix(m3d, mode, gp, atTime);
    this->reduceTo3dShellLayerStiffMtrx(answer, form, gp, m3d);
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
                                                        bool transpose) const
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
                                                        bool transpose) const
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
                                                             bool transpose) const
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
                                              const FloatArray &strainVector, bool transpose) const
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

    this->giveStrainVectorTranformationMtrx(tt, base, transpose);
    answer.beProductOf(tt, strainVector);
}


void
StructuralMaterial :: transformStressVectorTo(FloatArray &answer, const FloatMatrix &base,
                                              const FloatArray &stressVector, bool transpose) const
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

    this->giveStressVectorTranformationMtrx(tt, base, transpose);
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
StructuralMaterial :: setIPValue(const FloatArray &value, GaussPoint *aGaussPoint, InternalStateType type)
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(aGaussPoint) );
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
StructuralMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(aGaussPoint) );
    if ( type == IST_StressTensor ) {
        answer = status->giveStressVector();
        return 1;
    } else if ( type == IST_vonMisesStress ) {
        answer.resize(1);
        answer.at(1) = this->computeVonMisesStress( & status->giveStressVector() );
        return 1;
    } else if ( type == IST_StrainTensor ) {
        answer = status->giveStrainVector();
        return 1;
    } else if ( type == IST_StressTensorTemp ) {
        answer = status->giveTempStressVector();
        return 1;
    } else if ( type == IST_StrainTensorTemp ) {
        answer = status->giveTempStrainVector();
        return 1;
    } else if ( type == IST_PrincipalStressTensor || type == IST_PrincipalStressTempTensor ) {
        //         int indx;
        //         FloatArray st(6);
        FloatArray s;

        if ( type == IST_PrincipalStressTensor ) {
            s = status->giveStressVector();
        } else {
            s = status->giveTempStressVector();
        }

        //StructuralMaterial :: giveFullSymVectorForm(st, s, aGaussPoint->giveMaterialMode());
        this->computePrincipalValues(answer, s, principal_stress);
        return 1;
    } else if ( type == IST_PrincipalStrainTensor || type == IST_PrincipalStrainTempTensor ) {
        FloatArray st(6), s;

        if ( type == IST_PrincipalStrainTensor ) {
            s = status->giveStrainVector();
        } else {
            s = status->giveTempStrainVector();
        }

        StructuralMaterial :: giveFullSymVectorForm(st, s, aGaussPoint->giveMaterialMode());

        this->computePrincipalValues(answer, st, principal_strain);
        return 1;
    } else if ( type == IST_Temperature ) {
        /* add external source, if provided, such as staggered analysis */
        FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();
    FM_FieldPtr tf;
        int err;
        if ( ( tf = fm->giveField(FT_Temperature) ) ) {
            // temperature field registered
            FloatArray gcoords, et2;
            static_cast< StructuralElement * >( aGaussPoint->giveElement() )->computeGlobalCoordinates( gcoords, * aGaussPoint->giveCoordinates() );
            if ( ( err = tf->evaluateAt(answer, gcoords, VM_Total, atTime) ) ) {
                OOFEM_ERROR3("StructuralMaterial :: giveIPValue: tf->evaluateAt failed, element %d, error code %d", aGaussPoint->giveElement()->giveNumber(), err);
            }
        } else {
            answer.resize(1);
            answer.zero();
        }

        return 1;
    } else if ( type == IST_CylindricalStressTensor || type == IST_CylindricalStrainTensor ) {
        FloatArray gc, val = status->giveStressVector();
        FloatMatrix base(3, 3);
        static_cast< StructuralElement * >( aGaussPoint->giveElement() )->computeGlobalCoordinates( gc, * aGaussPoint->giveCoordinates() );
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
                transformStressVectorTo(answer, base, val, 0);
            } else {
                transformStrainVectorTo(answer, base, val, 0);
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
        return Material :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}


InternalStateValueType
StructuralMaterial :: giveIPValueType(InternalStateType type)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StressTensorTemp ) ||
        ( type == IST_CylindricalStressTensor ) ) {
        return ISVT_TENSOR_S3;
    }
    // strains components packed in engineering notation
    else if ( ( type == IST_StrainTensor ) || ( type == IST_StrainTensorTemp ) || ( type == IST_CylindricalStrainTensor ) ) {
        return ISVT_TENSOR_S3E;
    } else if ( ( type == IST_PrincipalStressTensor ) || ( type == IST_PrincipalStrainTensor ) ||
               ( type == IST_PrincipalStressTempTensor ) || ( type == IST_PrincipalStrainTempTensor ) ) {
        return ISVT_VECTOR;
    } else if ( ( type == IST_Temperature ) || ( type == IST_vonMisesStress ) ) {
        return ISVT_SCALAR;
    } else {
        return Material :: giveIPValueType(type);
    }
}


int
StructuralMaterial :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) ||
        ( type == IST_StressTensorTemp ) || ( type == IST_StrainTensorTemp ) ||
        ( type == IST_CylindricalStressTensor ) || ( type == IST_CylindricalStrainTensor ) ||
        ( type == IST_ShellForceMomentumTensor ) ) {
        this->giveStressStrainMask(answer, FullForm, mmode);
        return 1;
    } else if ( ( type == IST_PrincipalStressTensor ) || ( type == IST_PrincipalStrainTensor ) ||
               ( type == IST_PrincipalStressTempTensor ) || ( type == IST_PrincipalStrainTempTensor ) ) {
        if ( mmode == _3dMat || mmode == _PlaneStress || mmode == _PlaneStrain ) {
            answer.setValues(3, 1, 2, 3);
        } else if ( mmode == _1dMat ) {
            answer.setValues(3, 1, 0, 0);
        } else {
            return 0;
        }
        return 1;
    } else if ( type == IST_Temperature ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else {
        return Material :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}


int
StructuralMaterial :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) ||
        ( type == IST_StressTensorTemp ) || ( type == IST_StrainTensorTemp ) ||
        ( type == IST_CylindricalStressTensor ) || ( type == IST_CylindricalStrainTensor ) ||
        ( type == IST_ShellForceMomentumTensor ) ) {
        return StructuralMaterial :: giveSizeOfVoigtSymVector( aGaussPoint->giveMaterialMode() );
    } else if ( ( type == IST_PrincipalStressTensor ) || ( type == IST_PrincipalStrainTensor ) || ( type == IST_PrincipalPlasticStrainTensor ) ||
               ( type == IST_PrincipalStressTempTensor ) || ( type == IST_PrincipalStrainTempTensor ) ) {
        MaterialMode m = aGaussPoint->giveMaterialMode();
        if ( m == _3dMat || m == _PlaneStrain ) {
            return 3;
        } else if ( m == _PlaneStress ) {
            return 2;
        } else if ( m == _1dMat ) {
            return 1;
        } else {
            return 0;
        }
    } else if ( ( type == IST_Temperature ) || ( type == IST_vonMisesStress ) ) {
        return 1;
    } else if ( ( type == IST_DeformationGradientTensor ) || ( type == IST_FirstPKStressTensor ) ) {
        return this->giveSizeOfVoigtVector( aGaussPoint->giveMaterialMode() );
    } else {
        return Material :: giveIPValueSize(type, aGaussPoint);
    }
}


void
StructuralMaterial :: computeStressIndependentStrainVector(FloatArray &answer,
                                                           GaussPoint *gp, TimeStep *stepN, ValueModeType mode)
{
    FloatArray fullAnswer, et, e0, eigenstrain, answerTemper, answerEigenstrain;
    FloatMatrix GCS;
    MaterialMode matmode = gp->giveMaterialMode();
    StructuralCrossSection *crossSection =  dynamic_cast< StructuralCrossSection * >( gp->giveCrossSection() );
    Element *elem = gp->giveElement();
    StructuralElement *selem = dynamic_cast< StructuralElement * >( gp->giveElement() );


    answer.resize(0);
    answerTemper.resize(0);
    answerEigenstrain.resize(0);

    if ( stepN->giveIntrinsicTime() < this->castingTime ) {
        answer.zero();
        return;
    }

    //sum up all prescribed temperatures over an element
    //elem->computeResultingIPTemperatureAt(et, stepN, gp, mode);
    if ( selem ) {
        selem->computeResultingIPTemperatureAt(et, stepN, gp, mode);        // HUHU
    }

    //sum up all prescribed eigenstrain over an element
    if ( selem ) {
        selem->computeResultingIPEigenstrainAt(eigenstrain, stepN, gp, mode);
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
        if ( ( err = tf->evaluateAt(et2, gcoords, mode, stepN) ) ) {
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
        double thick, width;

        this->giveThermalDilatationVector(e0, gp, stepN);

        switch ( matmode ) {
        case _2dBeam:
            thick = crossSection->give(CS_Thickness);
            answerTemper.resize(3);
            answerTemper.zero();
            answerTemper.at(1) = e0.at(1) * ( et.at(1) - this->giveReferenceTemperature() );
            if ( et.giveSize() > 1 ) {
                answerTemper.at(2) = e0.at(1) * et.at(2) / thick;     // kappa_x
            }

            break;
        case _3dBeam:
            thick = crossSection->give(CS_Thickness);
            width = crossSection->give(CS_Width);
            answerTemper.resize(6);
            answerTemper.zero();
            answerTemper.at(1) = e0.at(1) * ( et.at(1) - this->giveReferenceTemperature() );
            if ( et.giveSize() > 1 ) {
                answerTemper.at(5) = e0.at(1) * et.at(2) / thick;     // kappa_y
                if ( et.giveSize() > 2 ) {
                    answerTemper.at(6) = e0.at(1) * et.at(3) / width;     // kappa_z
                }
            }

            break;
        case _2dPlate:
            thick = crossSection->give(CS_Thickness);
            if ( et.giveSize() > 1 ) {
                answerTemper.resize(5);
                answerTemper.zero();

                if ( et.giveSize() > 1 ) {
                    answerTemper.at(1) = e0.at(1) * et.at(2) / thick;     // kappa_x
                    answerTemper.at(2) = e0.at(2) * et.at(2) / thick;     // kappa_y
                }
            }

            break;
        case _3dShell:
            thick = crossSection->give(CS_Thickness);
            answerTemper.resize(8);
            answerTemper.zero();

            answerTemper.at(1) = e0.at(1) * ( et.at(1) - this->giveReferenceTemperature() );
            answerTemper.at(2) = e0.at(2) * ( et.at(1) - this->giveReferenceTemperature() );
            if ( et.giveSize() > 1 ) {
                answerTemper.at(4) = e0.at(1) * et.at(2) / thick;     // kappa_x
                answerTemper.at(5) = e0.at(2) * et.at(2) / thick;     // kappa_y
            }

            break;
        default:
            if ( e0.giveSize() ) {
                fullAnswer = e0;
                if ( mode == VM_Total ) {
                    fullAnswer.times(et.at(1) - this->referenceTemperature);
                } else {
                    fullAnswer.times( et.at(1) );
                }

                StructuralMaterial :: giveReducedSymVectorForm(answer, fullAnswer, gp->giveMaterialMode());
            }
        }
    }

    if ( eigenstrain.giveSize() ) { //found prescribed eigenstrain
        switch ( matmode ) {
        case _1dMat:
        case _3dMat:
        case _PlaneStress:
        case _PlaneStrain:
        case _3dRotContinuum:
        case _3dMicroplane:
            fullAnswer = eigenstrain;
            break;

        default:
            OOFEM_ERROR2( "StructuralMaterial :: Material mode %s for eigenstrains not supported", __MaterialModeToString(matmode) );
        }

        answerEigenstrain = fullAnswer;
        //this->giveReducedCharacteristicVector(answerEigenstrain, gp, fullAnswer);
    }

    //join temperature and eigenstrain vectors, compare vector sizes
    if ( answerTemper.giveSize() ) {
        answer = answerTemper;
        if ( answerEigenstrain.giveSize() ) {
            if ( answerTemper.giveSize() != answerEigenstrain.giveSize() ) {
                OOFEM_ERROR4( "StructuralMaterial :: Vector of temperature strains has the size %d which is different with the size of eigenstrain vector %d, element %d", answerTemper.giveSize(), answerEigenstrain.giveSize(), elem->giveNumber() );
            }

            answer.add(answerEigenstrain);
        }
    } else {
        if ( answerEigenstrain.giveSize() ) {
            answer = answerEigenstrain;
        }
    }
}


void 
StructuralMaterial :: giveFullSymVectorForm(FloatArray &answer, const FloatArray &vec, MaterialMode matMode)
{
    IntArray indx;
    answer.resize( StructuralMaterial :: giveVoigtSymVectorMask(indx, matMode) );
    answer.zero();
    answer.assemble(vec, indx);
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
StructuralMaterial :: giveFullVectorFormF(FloatArray &answer,const FloatArray &vec, MaterialMode matMode)
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
    answer.resize( indx.giveSize() );
    for ( int i = 1; i <= indx.giveSize(); i++ ) {
        answer.at(i) = vec.at( indx.at(i) );
    }
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
