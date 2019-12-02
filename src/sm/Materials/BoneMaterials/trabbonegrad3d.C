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


#include "trabbonegrad3d.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "error.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Material(TrabBoneGrad3D);

TrabBoneGrad3D :: TrabBoneGrad3D(int n, Domain *d) : TrabBone3D(n, d), GradientDamageMaterialExtensionInterface(d)
{
}

bool
TrabBoneGrad3D :: hasMaterialModeCapability(MaterialMode mode) const
{
    return mode == _3dMat;
}
void
TrabBoneGrad3D :: giveStiffnessMatrix(FloatMatrix &answer,
                                      MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    OOFEM_ERROR("Shouldn't be called");
}

void
TrabBoneGrad3D :: giveGradientDamageStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMat:
        answer = give3dMaterialStiffnessMatrix(mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR( "givePDGradMatrix_uu : unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

void
TrabBoneGrad3D :: giveGradientDamageStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMat:
        give3dKappaMatrix(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

void
TrabBoneGrad3D :: giveGradientDamageStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMat:
        give3dGprime(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}


void
TrabBoneGrad3D :: giveGradientDamageStiffnessMatrix_dd_BB(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMat:
        giveInternalLength(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR( "unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

void
TrabBoneGrad3D :: giveNonlocalInternalForces_N_factor(double &answer, double nlDamageDrivingVariable, GaussPoint *gp, TimeStep *tStep)
{
    answer = nlDamageDrivingVariable;
}

void
TrabBoneGrad3D :: giveNonlocalInternalForces_B_factor(FloatArray &answer, const FloatArray &nlDamageDrivingVariable_grad, GaussPoint *gp, TimeStep *tStep)
{
    answer = nlDamageDrivingVariable_grad;
    answer.times(internalLength * internalLength);
}


void
TrabBoneGrad3D :: computeLocalDamageDrivingVariable(double &answer, GaussPoint *gp, TimeStep *tStep)
{
    auto status = static_cast< TrabBoneGrad3DStatus * >( this->giveStatus(gp) );
    answer = status->giveTempKappa();
}


FloatMatrixF<6,6>
TrabBoneGrad3D :: give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< TrabBoneGrad3DStatus * >( this->giveStatus(gp) );
    FloatMatrixF<6,6> answer;
    if ( mode == ElasticStiffness ) {
        auto compliance = this->constructAnisoComplTensor();
        auto elasticity = inv(compliance);
        answer = elasticity;
    } else if ( mode == SecantStiffness ) {
        if ( printflag ) {
            printf("secant\n");
        }

        auto compliance = this->constructAnisoComplTensor();
        auto elasticity = inv(compliance);
        auto tempDam = status->giveTempDam();
        answer = elasticity * (1.0 - tempDam);
    } else if ( mode == TangentStiffness ) {
        double kappa = status->giveKappa();
        double tempKappa = status->giveTempKappa();
        double tempDam = status->giveTempDam();
        if ( tempKappa > kappa ) {
            // plastic loading
            // Imports
            auto &tempEffectiveStress = status->giveTempEffectiveStress();
            auto &plasFlowDirec = status->givePlasFlowDirec();
            auto &SSaTensor = status->giveSSaTensor();
            auto beta = status->giveBeta();
            // Construction of the dyadic product tensor
            auto prodTensor = Tdot(SSaTensor, plasFlowDirec);
            // Construction of the tangent stiffness second term
            auto tempTensor2 = dot(SSaTensor, plasFlowDirec);
            auto secondTerm = dyad(tempTensor2, prodTensor) * (-( 1.0 - tempDam ) / beta);
            // Construction of the tangent stiffness
            auto tangentMatrix = SSaTensor * (1.0 - tempDam) + secondTerm;
            if ( tempDam > status->giveDam() ) {
                double nlKappa =  status->giveNonlocalCumulatedStrain();
                double kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
                double gPrime = TrabBone3D :: computeDamageParamPrime(kappa);
                // Construction of the tangent stiffness third term
                auto thirdTerm = dyad(tempEffectiveStress, prodTensor) * gPrime * ( ( 1. - mParam ) / beta );
                tangentMatrix += thirdTerm;
            }

            answer = tangentMatrix;
        } else {
            // elastic behavior with damage
            // Construction of the secant stiffness
            auto compliance = this->constructAnisoComplTensor();
            auto elasticity = inv(compliance);
            answer = elasticity * (1.0 - tempDam);
        }
    }

    double g = status->giveDensG();
    if ( g <= 0 ) {
        double factor = gammaL0 * pow(rho, rL) + gammaP0 *pow(rho, rP) * ( tDens - 1 ) * pow(g, tDens - 2);
        answer += I6_I6 * factor;
    }
    return answer;
}




void
TrabBoneGrad3D :: give3dKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(1, 6);
    answer.zero();
    if ( mode == TangentStiffness ) {
        auto status = static_cast< TrabBoneGrad3DStatus * >( this->giveStatus(gp) );

        double kappa = status->giveKappa();
        double tempKappa = status->giveTempKappa();
        double dKappa = tempKappa - kappa;

        if ( dKappa > 0.0 ) {
            auto &plasFlowDirec = status->givePlasFlowDirec();
            auto &SSaTensor = status->giveSSaTensor();
            double beta = status->giveBeta();
            auto prodTensor = Tdot(SSaTensor, plasFlowDirec);
            for ( int i = 1; i <= 6; i++ ) {
                answer.at(1, i) = prodTensor.at(i);
            }

            answer.times(1. / beta);
        }
    }
}



void
TrabBoneGrad3D :: give3dGprime(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(6, 1);
    answer.zero();
    if ( mode == TangentStiffness ) {
        auto status = static_cast< TrabBoneGrad3DStatus * >( this->giveStatus(gp) );
        double tempKappa = status->giveTempKappa();
        double damage = status->giveDam();
        double nlKappa = status->giveNonlocalCumulatedStrain();
        double kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
        double tempDamage = TrabBone3D :: computeDamageParam(kappa);

        if ( ( tempDamage - damage ) > 0 ) {
            auto &tempEffStress = status->giveTempEffectiveStress();
            for ( int i = 1; i <= 6; i++ ) {
                answer.at(i, 1) = tempEffStress.at(i);
            }

            double kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
            double gPrime = TrabBone3D :: computeDamageParamPrime(kappa);
            answer.times(gPrime * mParam);
        }
    }
}

void
TrabBoneGrad3D :: giveInternalLength(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(1, 1);
    answer.at(1, 1) = L;
}

void
TrabBoneGrad3D :: giveRealStressVectorGradientDamage(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalCumulatedStrain, TimeStep *tStep)
{
    auto status = static_cast< TrabBoneGrad3DStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);

    TrabBone3D :: performPlasticityReturn(gp, totalStrain, tStep);

    double tempDamage = computeDamage(gp, tStep);
    auto tempEffStress = status->giveTempEffectiveStress();
    answer1 = (1 - tempDamage) * tempEffStress;
    answer2 = status->giveTempKappa();

    if ( densCrit != 0 ) {
        answer1.add(computeDensificationStress(gp, totalStrain, tStep));
    }

    status->letTempStrainVectorBe(totalStrain);
    status->setNonlocalCumulatedStrain(nonlocalCumulatedStrain);

    status->setTempDam(tempDamage);
    status->setTempEffectiveStress(tempEffStress);
    status->letTempStressVectorBe(answer1);
}


double
TrabBoneGrad3D :: computeCumPlastStrain(GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< TrabBoneGrad3DStatus * >( this->giveStatus(gp) );
    double localCumPlastStrain = status->giveTempKappa();
    double nlCumPlastStrain = status->giveNonlocalCumulatedStrain();
    return mParam * nlCumPlastStrain + ( 1 - mParam ) * localCumPlastStrain;
}


void
TrabBoneGrad3D :: initializeFrom(InputRecord &ir)
{
    TrabBone3D :: initializeFrom(ir);

    IR_GIVE_OPTIONAL_FIELD(ir, L, _IFT_TrabBoneGrad3D_L);
    if ( L < 0.0 ) {
        L = 0.0;
    }

    mParam = 2.;
    IR_GIVE_OPTIONAL_FIELD(ir, mParam, _IFT_TrabBoneGrad3D_m);
}


TrabBoneGrad3DStatus :: TrabBoneGrad3DStatus(GaussPoint *g) :
    TrabBone3DStatus(g)
{}


void
TrabBoneGrad3DStatus :: printOutputAt(FILE *file, TimeStep *tStep) const
{
    TrabBone3DStatus :: printOutputAt(file, tStep);
}


void
TrabBoneGrad3DStatus :: initTempStatus()
{
    TrabBone3DStatus :: initTempStatus();
}


void
TrabBoneGrad3DStatus :: updateYourself(TimeStep *tStep)
{
    StructuralMaterialStatus :: updateYourself(tStep);
    this->kappa = this->tempKappa;
    this->dam = this->tempDam;
    this->tsed = this->tempTSED;
    this->plasDef = this->tempPlasDef;
    //TrabBone3DStatus :: updateYourself(tStep);
}
} // end namespace oofem
