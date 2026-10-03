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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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
        return 1;
    } else if ( type == IST_InterfaceTraction ) {
        answer = status->giveTraction();
        return 1;
    } else if ( type == IST_InterfaceFirstPKTraction ) {
        answer = status->giveFirstPKTraction();
        return 1;
    } else if ( type == IST_DeformationGradientTensor ) {
        answer = to_voigt_form_33( status->giveF() );
        return 1;
    } else if ( type == IST_InterfaceNormal ) {
        answer = status->giveNormal();
        return 1;
    } else if ( type == IST_DamageScalar ) {
        answer.resize(1);
        answer.at(1) = status->giveDamage();
        return 1;
    } else {
        return Material :: giveIPValue(answer, gp, type, tStep);
    }
}


void
StructuralInterfaceMaterial :: initializeFrom(const std::shared_ptr<InputRecord> &ir)
{
    IR_GIVE_OPTIONAL_FIELD(ir, this->useNumericalTangent, _IFT_StructuralInterfaceMaterial_useNumericalTangent);
}


void
StructuralInterfaceMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    Material :: giveInputRecord(input);
}


FloatMatrixF<1,1>
StructuralInterfaceMaterial :: giveStiffnessMatrix_Eng_n(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto d = giveStiffnessMatrix_Eng_ntt(mode, gp, tStep);
    return {d(0,0)};
}


FloatMatrixF<2,2>
StructuralInterfaceMaterial :: giveStiffnessMatrix_Eng_nt(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto d = giveStiffnessMatrix_Eng_ntt(mode, gp, tStep);
    return d({0,1}, {0,1});
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

FloatMatrixF<3,3>
StructuralInterfaceMaterial :: giveStiffnessMatrix_Eng_ntt(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    // Default implementation. Will use large deformation dT/dj as stiffness
    return giveStiffnessMatrix_dTdj_ntt(mode, gp, tStep);
}


FloatMatrixF<1,1>
StructuralInterfaceMaterial :: giveStiffnessMatrix_dTdj_n(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto d = this->giveStiffnessMatrix_dTdj_Num_ntt(gp, tStep);
    return {d(0,0)};
}


FloatMatrixF<2,2>
StructuralInterfaceMaterial :: giveStiffnessMatrix_dTdj_nt(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    auto d = this->giveStiffnessMatrix_dTdj_Num_ntt(gp, tStep);
    return d({0,1}, {0,1});
}

FloatMatrixF<3,3>
StructuralInterfaceMaterial :: giveStiffnessMatrix_dTdj_ntt(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_WARNING("Using numerical tangent");
    return this->giveStiffnessMatrix_dTdj_Num_ntt(gp, tStep);
}


double
StructuralInterfaceMaterial :: giveFirstPKTraction_n(double jump, double reducedF, GaussPoint *gp, TimeStep *tStep) const
{
    return this->giveFirstPKTraction_ntt({ jump,  0., 0. }, diag<3>({reducedF, 1., 1.}), gp, tStep).at(1);
}

FloatArrayF<2>
StructuralInterfaceMaterial :: giveFirstPKTraction_nt(const FloatArrayF<2> &jump, const FloatMatrixF<2,2> &reducedF, GaussPoint *gp, TimeStep *tStep) const
{
    FloatMatrixF<3,3> F = {
        reducedF(0,0), reducedF(1,0), 0.,
        reducedF(0,1), reducedF(1,1), 0.,
        0., 0., 1.
    };

    auto traction3D = this->giveFirstPKTraction_ntt({ jump.at(1), jump.at(2), 0. }, F, gp, tStep);
    return { traction3D.at(1), traction3D.at(2) };
}

double
StructuralInterfaceMaterial :: giveEngTraction_n(double jump, GaussPoint *gp, TimeStep *tStep) const
{
    auto traction3D = this->giveEngTraction_ntt({ jump, 0., 0.}, gp, tStep);
    return traction3D.at(1);
}

FloatArrayF<2>
StructuralInterfaceMaterial :: giveEngTraction_nt(const FloatArrayF<2> &jump, GaussPoint *gp, TimeStep *tStep) const
{
    auto traction3D = this->giveEngTraction_ntt({ jump.at(1), jump.at(2), 0. }, gp, tStep);
    return { traction3D.at(1), traction3D.at(2) };
}

FloatArrayF<3>
StructuralInterfaceMaterial :: giveEngTraction_ntt(const FloatArrayF<3> &jump, GaussPoint *gp, TimeStep *tStep) const
{
    // Default implementation calls first PK version with F = I
    return giveFirstPKTraction_ntt(jump, eye<3>(), gp, tStep);
}



// Numerical tangents
FloatMatrixF<1,1>
StructuralInterfaceMaterial :: giveStiffnessMatrix_dTdj_Num_n(GaussPoint *gp, TimeStep *tStep) const
{
    // Default implementation for computation of the numerical tangent
    // Computes the material stiffness using a central difference method

    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus( gp ) );
    double eps = 1.0e-9;
    double F = status->giveTempF().at(1,1);
    double jump = status->giveTempJump().at(3);

    double tPlus = this->giveFirstPKTraction_n( jump + eps, F, gp, tStep );
    double tMinus = this->giveFirstPKTraction_n( jump - eps, F, gp, tStep );

    this->giveFirstPKTraction_n(jump, F, gp, tStep); // reset temp values by recomputing the stress
    return (tPlus - tMinus) / ( 2 * eps );
}

FloatMatrixF<2,2>
StructuralInterfaceMaterial :: giveStiffnessMatrix_dTdj_Num_nt(GaussPoint *gp, TimeStep *tStep) const
{
    // Default implementation for computation of the numerical tangent
    // Computes the material stiffness using a central difference method

    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus( gp ) );
    double eps = 1.0e-9;
    const auto &F3 = status->giveTempF();
    FloatMatrixF<2,2> F = {F3(0,0), F3(1,0), F3(0,1), F3(1,1)};
    FloatArrayF<2> jump = { status->giveTempJump().at(1), status->giveTempJump().at(2) };

    FloatMatrixF<2,2> answer;
    for ( int i = 1; i <= 2; i++ ) {
        auto jumpPlus = jump;
        auto jumpMinus = jump;
        jumpPlus.at( i ) += eps;
        jumpMinus.at( i ) -= eps;
        auto TPlus = this->giveFirstPKTraction_nt(jumpPlus, F, gp, tStep);
        auto TMinus = this->giveFirstPKTraction_nt(jumpMinus, F, gp, tStep);

        auto Kcolumn = TPlus - TMinus;
        answer.setColumn(Kcolumn, i);
    }
    answer *=  1.0 / ( 2 * eps );
    this->giveFirstPKTraction_nt(jump, F, gp, tStep); // reset temp values by recomputing the stress
    return answer;
}

FloatMatrixF<3,3>
StructuralInterfaceMaterial :: giveStiffnessMatrix_dTdj_Num_ntt(GaussPoint *gp, TimeStep *tStep) const
{
    // Default implementation for computation of the numerical tangent
    // Computes the material stiffness using a central difference method
    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus( gp ) );
    double eps = 1.0e-9;
    const auto &F = status->giveTempF();
    const auto &jump = status->giveTempJump();

    FloatMatrixF<3,3> answer;
    for(int i = 1; i <= 3; i++) {
        auto jumpPlus = jump;
        auto jumpMinus = jump;
        jumpPlus.at( i ) += eps;
        jumpMinus.at( i ) -= eps;
        auto TPlus = this->giveFirstPKTraction_ntt(jumpPlus, F, gp, tStep );
        auto TMinus = this->giveFirstPKTraction_ntt(jumpMinus, F, gp, tStep );

        auto Kcolumn = TPlus - TMinus;
        answer.setColumn(Kcolumn, i);
    }
    answer *= 1.0 / ( 2 * eps );
    this->giveFirstPKTraction_ntt(jump, F, gp, tStep); // reset temp values by recomputing the stress
    return answer;
}

FloatMatrixF<3,3>
StructuralInterfaceMaterial :: giveStiffnessMatrix_dTdj_Num_tnn(GaussPoint *gp, TimeStep *tStep) const
{
    // Default implementation for computation of the numerical tangent
    // Computes the material stiffness using a central difference method
    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus( gp ) );
    double eps = 1.0e-9;
    const auto &F = status->giveTempF();
    const auto &jump = status->giveTempJump();

    FloatMatrixF<3,3> answer;
    for(int i = 1; i <= 3; i++) {
        auto jumpPlus = jump;
        auto jumpMinus = jump;
        jumpPlus.at( i ) += eps;
        jumpMinus.at( i ) -= eps;
        auto TPlus = this->giveFirstPKTraction_tnn(jumpPlus, F, gp, tStep );
        auto TMinus = this->giveFirstPKTraction_tnn(jumpMinus, F, gp, tStep );

        auto Kcolumn = TPlus - TMinus;
        answer.setColumn(Kcolumn, i);
    }
    answer *= 1.0 / ( 2 * eps );
    this->giveFirstPKTraction_tnn(jump, F, gp, tStep); // reset temp values by recomputing the stress
    return answer;
}

FloatMatrixF<1,1>
StructuralInterfaceMaterial :: giveStiffnessMatrix_Eng_Num_n(GaussPoint *gp, TimeStep *tStep) const
{
    // Default implementation for computation of the numerical tangent d(sig)/d(jump)
    // Computes the material stiffness using a central difference method
    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus( gp ) );
    double eps = 1.0e-9;
    double jump = status->giveTempJump().at(1);

    auto tPlus = this->giveEngTraction_n(jump + eps, gp, tStep);
    auto tMinus = this->giveEngTraction_n(jump - eps, gp, tStep);

    this->giveEngTraction_n( jump, gp, tStep ); // reset temp values by recomputing the stress
    return (tPlus - tMinus) / ( 2 * eps );
}

FloatMatrixF<2,2>
StructuralInterfaceMaterial :: giveStiffnessMatrix_Eng_Num_nt(GaussPoint *gp, TimeStep *tStep) const
{
    // Default implementation for computation of the numerical tangent d(sig)/d(jump)
    // Computes the material stiffness using a central difference method

    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus( gp ) );
    double eps = 1.0e-12;
    FloatArrayF<2> jump = {status->giveTempJump().at(1), status->giveTempJump().at(2)};

    FloatMatrixF<2,2> answer;
    for(int i = 1; i <= 2; i++) {
        auto jumpPlus = jump;
        auto jumpMinus = jump;
        jumpPlus.at( i ) += eps;
        jumpMinus.at( i ) -= eps;
        auto tPlus = this->giveEngTraction_nt(jumpPlus, gp, tStep);
        auto tMinus = this->giveEngTraction_nt(jumpMinus, gp, tStep);

        auto Kcolumn = tPlus - tMinus;
        answer.setColumn(Kcolumn, i);
    }
    answer *= 1.0 / ( 2 * eps );
    this->giveEngTraction_nt( jump, gp, tStep ); // reset temp values by recomputing the stress
    return answer;
}

FloatMatrixF<3,3>
StructuralInterfaceMaterial :: giveStiffnessMatrix_Eng_Num_ntt(GaussPoint *gp, TimeStep *tStep) const
{
    // Default implementation for computation of the numerical tangent d(sig)/d(jump)
    // Computes the material stiffness using a central difference method

    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus( gp ) );
    double eps = 1.0e-9;
    const auto &jump = status->giveTempJump();

    FloatMatrixF<3,3> answer;
    for(int i = 1; i <= 3; i++) {
        auto jumpPlus = jump;
        auto jumpMinus = jump;
        jumpPlus.at( i ) += eps;
        jumpMinus.at( i ) -= eps;
        auto tPlus = this->giveEngTraction_ntt(jumpPlus, gp, tStep);
        auto tMinus = this->giveEngTraction_ntt(jumpMinus, gp, tStep);

        auto Kcolumn = tPlus - tMinus;
        answer.setColumn(Kcolumn, i);
    }
    answer *= 1.0 / ( 2 * eps );
    this->giveEngTraction_ntt(jump, gp, tStep ); // reset temp values by recomputing the stress
    return answer;
}

FloatMatrixF<3,3>
StructuralInterfaceMaterial :: giveStiffnessMatrix_Eng_Num_tnn(GaussPoint *gp, TimeStep *tStep) const
{
    // Default implementation for computation of the numerical tangent d(sig)/d(jump)
    // Computes the material stiffness using a central difference method

    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( this->giveStatus( gp ) );
    double eps = 1.0e-9;
    const auto &jump = status->giveTempJump();

    FloatMatrixF<3,3> answer;
    for(int i = 1; i <= 3; i++) {
        auto jumpPlus = jump;
        auto jumpMinus = jump;
        jumpPlus.at( i ) += eps;
        jumpMinus.at( i ) -= eps;
        auto tPlus = this->giveEngTraction_tnn(jumpPlus, gp, tStep);
        auto tMinus = this->giveEngTraction_tnn(jumpMinus, gp, tStep);

        auto Kcolumn = tPlus - tMinus;
        answer.setColumn(Kcolumn, i);
    }
    answer *= 1.0 / ( 2 * eps );
    this->giveEngTraction_tnn(jump, gp, tStep ); // reset temp values by recomputing the stress
    return answer;
}

} // end namespace oofem
