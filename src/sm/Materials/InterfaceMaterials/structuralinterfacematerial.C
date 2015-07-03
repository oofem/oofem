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

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"
#include "dynamicinputrecord.h"
#include "gausspoint.h"

namespace oofem {
StructuralInterfaceMaterial :: StructuralInterfaceMaterial(int n, Domain *d) : Material(n, d)
{
    this->useNumericalTangent = false;
}


int
StructuralInterfaceMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_InterfaceJump ) {
        answer = status->giveJump();
        answer.resizeWithValues(3); // In case some model is not storing all components.
        return 1;
    } else if ( type == IST_InterfaceTraction ) {
        answer = status->giveTraction();
        answer.resizeWithValues(3);
        return 1;
    } else if ( type == IST_InterfaceFirstPKTraction ) {
        answer = status->giveFirstPKTraction();
        answer = status->giveTempFirstPKTraction();
        answer.resizeWithValues(3);
        return 1;
    } else if ( type == IST_DeformationGradientTensor ) {
        answer.beVectorForm( status->giveF() );
        return 1;
    } else {
        return Material :: giveIPValue(answer, gp, type, tStep);
    }
}


IRResultType
StructuralInterfaceMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                    // Required by IR_GIVE_FIELD macro

    IR_GIVE_OPTIONAL_FIELD(ir, this->useNumericalTangent, _IFT_StructuralInterfaceMaterial_useNumericalTangent);

    return IRRT_OK;
}


void
StructuralInterfaceMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    Material :: giveInputRecord(input);
}


void
StructuralInterfaceMaterial :: give1dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix answer3D;
    give3dStiffnessMatrix_Eng(answer3D, mode, gp, tStep);
    answer.resize(1, 1);
    answer.at(1, 1) = answer3D.at(3, 3);

}


void
StructuralInterfaceMaterial :: give2dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix answer3D;
    give3dStiffnessMatrix_Eng(answer3D, mode, gp, tStep);
    IntArray mask = {1, 3};
    answer.beSubMatrixOf(answer3D, mask, mask);

#if 0
    answer3D.printYourself("analytical tangent");

    FloatMatrix answerNum;
    give2dStiffnessMatrix_Eng_Num(answerNum, mode, gp, tStep);
    answerNum.printYourself("numerical tangent");

    FloatMatrix comp;
    comp = answer3D;
    comp.subtract(answerNum);
    comp.printYourself("difference in numerical tangent to mat method");
#endif
}

void
StructuralInterfaceMaterial :: give3dStiffnessMatrix_Eng(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    // Default implementation. Will use large deformation dT/dj as stiffness
    give3dStiffnessMatrix_dTdj(answer, mode, gp, tStep);
}


void
StructuralInterfaceMaterial :: give1dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->give3dStiffnessMatrix_dTdj_Num(answer, gp, tStep);
    answer.resize(1, 1);
    answer.at(1, 1) = answer.at(3, 3);
}


void
StructuralInterfaceMaterial :: give2dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->give3dStiffnessMatrix_dTdj_Num(answer, gp, tStep);
    answer.resize(2, 2);
    answer.at(1, 1) = answer.at(1, 1);
    answer.at(1, 2) = answer.at(1, 3);
    answer.at(2, 1) = answer.at(3, 1);
    answer.at(2, 2) = answer.at(3, 3);
}


void
StructuralInterfaceMaterial :: give3dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->give3dStiffnessMatrix_dTdj_Num(answer, gp, tStep);
    OOFEM_WARNING("Using numerical tangent");
}


void
StructuralInterfaceMaterial :: giveFirstPKTraction_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump,
                                         const FloatMatrix &reducedF, TimeStep *tStep)
{
    FloatArray traction3D;
    FloatMatrix F(3, 3);
    F.at(1, 1) = 1.;
    F.at(2, 2) = 1.;
    F.at(3, 3) = reducedF.at(1, 1);
    
    this->giveFirstPKTraction_3d(traction3D, gp, { 0., 0., jump.at(1) }, F, tStep);
    answer = FloatArray{ traction3D.at(3) };
}

void
StructuralInterfaceMaterial :: giveFirstPKTraction_2d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump,
                                         const FloatMatrix &reducedF, TimeStep *tStep)
{
    FloatArray traction3D;
    FloatMatrix F(3, 3);
    F.at(1, 1) = reducedF.at(1, 1);
    F.at(1, 3) = reducedF.at(1, 2);
    F.at(2, 2) = 1.;
    F.at(3, 1) = reducedF.at(2, 1);
    F.at(3, 3) = reducedF.at(2, 2);

    this->giveFirstPKTraction_3d(traction3D, gp, { jump.at(1), 0., jump.at(2) }, F, tStep);
    answer = { traction3D.at(1), traction3D.at(3) };
}

void
StructuralInterfaceMaterial :: giveEngTraction_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
{
    FloatArray traction3D;
    this->giveEngTraction_3d(traction3D, gp, { 0., 0., jump.at(1) }, tStep);
    answer = FloatArray{ traction3D.at(3) };
}

void
StructuralInterfaceMaterial :: giveEngTraction_2d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
{
    FloatArray traction3D;
    this->giveEngTraction_3d(traction3D, gp, { jump.at(1), 0.0, jump.at(2) }, tStep);
    answer = { traction3D.at(1), traction3D.at(3) };
}

void
StructuralInterfaceMaterial :: giveEngTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, TimeStep *tStep)
{
    // Default implementation calls first PK version with F = I
    FloatMatrix F(3, 3);
    F.beUnitMatrix();
    giveFirstPKTraction_3d(answer, gp, jump, F, tStep);
}



// Numerical tangents
void
StructuralInterfaceMaterial :: give1dStiffnessMatrix_dTdj_Num(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Default implementation for computation of the numerical tangent
    // Computes the material stiffness using a central difference method

    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus( gp ) );
    double eps = 1.0e-9;
    const FloatMatrix &F = status->giveTempF();
    FloatArray T, TPlus, TMinus;
    FloatArray jumpPlus, jumpMinus, Kcolumn;
    FloatArray jump =  FloatArray{status->giveTempJump().at(3)};
    int dim = jump.giveSize();
    answer.resize( dim, dim );

    for ( int i = 1; i <= dim; i++ ) {
        jumpPlus = jumpMinus = jump;
        jumpPlus.at( i ) += eps;
        jumpMinus.at( i ) -= eps;
        this->giveFirstPKTraction_1d( TPlus, gp, jumpPlus, F, tStep );
        this->giveFirstPKTraction_1d( TMinus, gp, jumpMinus, F, tStep );

        Kcolumn.beDifferenceOf(TPlus, TMinus);
        answer.setColumn( Kcolumn, i );
    }
    answer.times( 1.0 / ( 2 * eps ) );
    this->giveFirstPKTraction_1d(T, gp, jump, F, tStep); // reset temp values by recomputing the stress
}

void
StructuralInterfaceMaterial :: give2dStiffnessMatrix_dTdj_Num(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Default implementation for computation of the numerical tangent
    // Computes the material stiffness using a central difference method

    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus( gp ) );
    double eps = 1.0e-9;
    const FloatMatrix &F = status->giveTempF();
    FloatArray T, TPlus, TMinus;
    FloatArray jumpPlus, jumpMinus, Kcolumn;
    FloatArray jump =  {status->giveTempJump().at(1), status->giveTempJump().at(3)};
    int dim = jump.giveSize();
    answer.resize( dim, dim );

    for ( int i = 1; i <= dim; i++ ) {
        jumpPlus = jumpMinus = jump;
        jumpPlus.at( i ) += eps;
        jumpMinus.at( i ) -= eps;
        this->giveFirstPKTraction_2d(TPlus, gp, jumpPlus, F, tStep);
        this->giveFirstPKTraction_2d(TMinus, gp, jumpMinus, F, tStep);

        Kcolumn.beDifferenceOf(TPlus, TMinus);
        answer.setColumn(Kcolumn, i);
    }
    answer.times( 1.0 / ( 2 * eps ) );
    this->giveFirstPKTraction_2d(T, gp, jump, F, tStep); // reset temp values by recomputing the stress
}

void
StructuralInterfaceMaterial :: give3dStiffnessMatrix_dTdj_Num(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Default implementation for computation of the numerical tangent
    // Computes the material stiffness using a central difference method

    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus( gp ) );
    double eps = 1.0e-9;
    const FloatMatrix &F = status->giveTempF();
    FloatArray T, TPlus, TMinus;
    FloatArray jumpPlus, jumpMinus, Kcolumn;
    const FloatArray &jump = status->giveTempJump();
    int dim = jump.giveSize();
    answer.resize( dim, dim );

    for(int i = 1; i <= dim; i++) {
        jumpPlus = jumpMinus = jump;
        jumpPlus.at( i ) += eps;
        jumpMinus.at( i ) -= eps;
        this->giveFirstPKTraction_3d( TPlus, gp, jumpPlus, F, tStep );
        this->giveFirstPKTraction_3d( TMinus, gp, jumpMinus, F, tStep );

        Kcolumn.beDifferenceOf(TPlus, TMinus);
        answer.setColumn(Kcolumn, i);
    }
    answer.times( 1.0 / ( 2 * eps ) );
    this->giveFirstPKTraction_3d(T, gp, jump, F, tStep); // reset temp values by recomputing the stress
}


void
StructuralInterfaceMaterial :: give1dStiffnessMatrix_Eng_Num(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Default implementation for computation of the numerical tangent d(sig)/d(jump)
    // Computes the material stiffness using a central difference method

    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus( gp ) );
    double eps = 1.0e-9;
    FloatArray t, tPlus, tMinus;
    FloatArray jumpPlus, jumpMinus, Kcolumn;
    FloatArray jump = FloatArray{status->giveTempJump().at(3)};
    int dim = jump.giveSize();
    answer.resize( dim, dim );
    answer.zero();

    for(int i = 1; i <= dim; i++) {
        jumpPlus = jumpMinus = jump;
        jumpPlus.at( i ) += eps;
        jumpMinus.at( i ) -= eps;
        this->giveEngTraction_1d(tPlus, gp, jumpPlus, tStep);
        this->giveEngTraction_1d(tMinus, gp, jumpMinus, tStep);

        Kcolumn.beDifferenceOf(tPlus, tMinus);
        answer.setColumn(Kcolumn, i);
    }
    answer.times( 1.0 / ( 2 * eps ) );
    this->giveEngTraction_1d( t, gp, jump, tStep ); // reset temp values by recomputing the stress
}

void
StructuralInterfaceMaterial :: give2dStiffnessMatrix_Eng_Num(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Default implementation for computation of the numerical tangent d(sig)/d(jump)
    // Computes the material stiffness using a central difference method

    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus( gp ) );
    double eps = 1.0e-9;
    FloatArray t, tPlus, tMinus;
    FloatArray tempJump, jumpPlus, jumpMinus, Kcolumn;
    FloatArray jump = {status->giveTempJump().at(1), status->giveTempJump().at(3)};
    int dim = jump.giveSize();
    answer.resize(dim, dim);
    answer.zero();

    for(int i = 1; i <= dim; i++) {
        jumpPlus = jumpMinus = jump;
        jumpPlus.at( i ) += eps;
        jumpMinus.at( i ) -= eps;
        this->giveEngTraction_2d(tPlus, gp, jumpPlus, tStep);
        this->giveEngTraction_2d(tMinus, gp, jumpMinus, tStep);

        Kcolumn.beDifferenceOf(tPlus, tMinus);
        answer.setColumn(Kcolumn, i);
    }
    answer.times( 1.0 / ( 2 * eps ) );
    this->giveEngTraction_2d( t, gp, jump, tStep ); // reset temp values by recomputing the stress
}

void
StructuralInterfaceMaterial :: give3dStiffnessMatrix_Eng_Num(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep)
{
    // Default implementation for computation of the numerical tangent d(sig)/d(jump)
    // Computes the material stiffness using a central difference method

    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus( gp ) );
    double eps = 1.0e-9;
    FloatArray t, tPlus, tMinus;
    FloatArray jumpPlus, jumpMinus, Kcolumn;
    const FloatArray &jump = status->giveTempJump();
    int dim = jump.giveSize();
    answer.resize(dim, dim);
    answer.zero();

    for(int i = 1; i <= dim; i++) {
        jumpPlus = jumpMinus = jump;
        jumpPlus.at( i ) += eps;
        jumpMinus.at( i ) -= eps;
        this->giveEngTraction_3d(tPlus, gp, jumpPlus, tStep);
        this->giveEngTraction_3d(tMinus, gp, jumpMinus, tStep);

        Kcolumn.beDifferenceOf(tPlus, tMinus);
        answer.setColumn(Kcolumn, i);
    }
    answer.times( 1.0 / ( 2 * eps ) );
    this->giveEngTraction_3d( t, gp, jump, tStep ); // reset temp values by recomputing the stress
}

} // end namespace oofem
