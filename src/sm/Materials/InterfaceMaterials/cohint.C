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

#include "cohint.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {

REGISTER_Material(CohesiveInterfaceMaterial);

CohesiveInterfaceMaterial :: CohesiveInterfaceMaterial(int n, Domain *d) : StructuralInterfaceMaterial(n, d) { }

FloatArrayF<3>
CohesiveInterfaceMaterial :: giveEngTraction_3d(const FloatArrayF<3> &jump, GaussPoint *gp, TimeStep *tStep) const
{
    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus(gp) );

    double x = jump.at(1) + transitionOpening;

    FloatArrayF<3> answer;

    if ( stiffCoeffKn == 1. ){ //tension stiffness = compression stiffness
        answer.at(1) = kn*x;
    } else {
        //composed function from two atan's to have smooth intersection between tension and compression
        answer.at(1) = (M_PI/2. + atan(smoothMag*x))/M_PI*kn*stiffCoeffKn*x + (M_PI/2.-atan(smoothMag*x))/M_PI*kn*x;
    }

    // shear part of elastic stress-strain law
    answer.at(2) = ks * jump.at(2);
    answer.at(3) = ks * jump.at(3);

    // update gp
    status->letTempJumpBe(jump);
    status->letTempTractionBe(answer);
    return answer;
}


FloatMatrixF<3,3>
CohesiveInterfaceMaterial :: give3dStiffnessMatrix_Eng(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus(gp) );

//     double normalJump = status->giveTempJump().at(1);
//     if (normalJump > 0.) {
//         if(normalJump<transitionOpening){ // reduce traction in tension
//             double stiffTmp = kn*stiffCoeffKn + (kn - kn*stiffCoeffKn) * (1. - normalJump/transitionOpening);
//             answer.at(1, 1) = stiffTmp;
//         } else {
//             answer.at(1, 1) = kn * stiffCoeffKn;
//         }
//     } else {
//         // standard part of elastic stress-strain law
//         answer.at(1, 1) = kn;
//     }
//     
//     if ( rMode == SecantStiffness || rMode == TangentStiffness ) {
//         if ( normalJump + transitionOpening <= 0. ) { //local CS
//             answer.at(1, 1) = kn; //compression
//         } else {
//             answer.at(1, 1) = 0*kn*stiffCoeffKn; //tension
//         }
//     } else {
//         answer.at(1, 1) = kn; 
//     }

    double x = status->giveTempJump().at(1) + transitionOpening;

    double normal;
    if ( stiffCoeffKn == 1. ) { //tension stiffness = compression stiffness
        normal = kn;
    } else {
        //TangentStiffness by derivating traction with regards to x (=relative displacement)
        normal = (M_PI/2. + atan(smoothMag*x))/M_PI*kn*stiffCoeffKn + (M_PI/2.-atan(smoothMag*x))/M_PI*kn + smoothMag*kn*stiffCoeffKn*x/M_PI/(smoothMag*smoothMag*x*x+1) - smoothMag*kn*x/M_PI/(smoothMag*smoothMag*x*x+1);
    }

    return diag<3>({normal, ks, ks});
}


void
CohesiveInterfaceMaterial :: initializeFrom(InputRecord &ir)
{
    StructuralInterfaceMaterial :: initializeFrom(ir);

    // elastic parameters
    IR_GIVE_FIELD(ir, kn, _IFT_CohesiveInterfaceMaterial_kn);
    IR_GIVE_FIELD(ir, ks, _IFT_CohesiveInterfaceMaterial_ks);
    // compliance when in tension
    stiffCoeffKn = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, stiffCoeffKn, _IFT_CohesiveInterfaceMaterial_stiffCoeffKn);
    transitionOpening = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, transitionOpening, _IFT_CohesiveInterfaceMaterial_transitionopening);//Better term DisplacementOffset.

    smoothMag = 1.e+8;
    IR_GIVE_OPTIONAL_FIELD(ir, smoothMag, _IFT_CohesiveInterfaceMaterial_smoothMag);
}

void
CohesiveInterfaceMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralInterfaceMaterial :: giveInputRecord(input);
    input.setField(this->kn, _IFT_CohesiveInterfaceMaterial_kn);
    input.setField(this->ks, _IFT_CohesiveInterfaceMaterial_ks);
    input.setField(this->stiffCoeffKn, _IFT_CohesiveInterfaceMaterial_stiffCoeffKn);
    input.setField(this->stiffCoeffKn, _IFT_CohesiveInterfaceMaterial_transitionopening);
}

} // namespace oofem
