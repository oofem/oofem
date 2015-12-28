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

TrabBoneGrad3D :: TrabBoneGrad3D(int n, Domain *d) : TrabBone3D(n, d), GradDpMaterialExtensionInterface(d)
{
    L = 0.;
}

TrabBoneGrad3D :: ~TrabBoneGrad3D()
{ }

int
TrabBoneGrad3D :: hasMaterialModeCapability(MaterialMode mode)
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
TrabBoneGrad3D :: givePDGradMatrix_uu(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMat:
        give3dMaterialStiffnessMatrix(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR( "givePDGradMatrix_uu : unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

void
TrabBoneGrad3D :: givePDGradMatrix_ku(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMat:
        give3dKappaMatrix(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR("unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

void
TrabBoneGrad3D :: givePDGradMatrix_uk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMat:
        give3dGprime(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR("unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

void
TrabBoneGrad3D :: givePDGradMatrix_kk(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _3dMat:
        giveInternalLength(answer, mode, gp, tStep);
        break;
    default:
        OOFEM_ERROR("unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

void
TrabBoneGrad3D :: givePDGradMatrix_LD(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    OOFEM_ERROR("unknown mode (%s)", __MaterialModeToString(mMode) );
}


void
TrabBoneGrad3D :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    double tempDam, beta, tempKappa, kappa;
    FloatArray tempEffectiveStress, tempTensor2, prodTensor, plasFlowDirec;
    FloatMatrix elasticity, compliance, SSaTensor, secondTerm, thirdTerm, tangentMatrix;
    TrabBoneGrad3DStatus *status = static_cast< TrabBoneGrad3DStatus * >( this->giveStatus(gp) );

    if ( mode == ElasticStiffness ) {
        this->constructAnisoComplTensor(compliance);
        elasticity.beInverseOf(compliance);
        answer = elasticity;
    } else if ( mode == SecantStiffness || ( ( mode == TangentStiffness ) && ( status->giveNsubsteps() > 1 ) ) ) {
        if ( printflag ) {
            printf("secant\n");
        }

        this->constructAnisoComplTensor(compliance);
        elasticity.beInverseOf(compliance);
        tempDam = status->giveTempDam();
        answer = elasticity;
        answer.times(1.0 - tempDam);
    } else if ( mode == TangentStiffness ) {
        kappa = status->giveKappa();
        tempKappa = status->giveTempKappa();
        tempDam = status->giveTempDam();
        if ( tempKappa > kappa ) {
            // plastic loading
            // Imports
            tempEffectiveStress = status->giveTempEffectiveStress();
            plasFlowDirec = status->givePlasFlowDirec();
            SSaTensor = status->giveSSaTensor();
            beta = status->giveBeta();
            // Construction of the dyadic product tensor
            prodTensor.beTProductOf(SSaTensor, plasFlowDirec);
            // Construction of the tangent stiffness second term
            tempTensor2.beProductOf(SSaTensor, plasFlowDirec);
            secondTerm.beDyadicProductOf(tempTensor2, prodTensor);
            secondTerm.times(-( 1.0 - tempDam ) / beta);
            // Construction of the tangent stiffness
            tangentMatrix = SSaTensor;
            tangentMatrix.times(1.0 - tempDam);
            tangentMatrix.add(secondTerm);
            if ( tempDam > status->giveDam() ) {
                double nlKappa =  status->giveNonlocalCumulatedStrain();
                kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
                double gPrime = TrabBone3D :: computeDamageParamPrime(kappa);
                // Construction of the tangent stiffness third term
                thirdTerm.beDyadicProductOf(tempEffectiveStress, prodTensor);
                thirdTerm.times(gPrime);
                thirdTerm.times( ( 1. - mParam ) / beta );
                tangentMatrix.add(thirdTerm);
            }

            answer = tangentMatrix;
        } else {
            // elastic behavior with damage
            // Construction of the secant stiffness
            this->constructAnisoComplTensor(compliance);
            elasticity.beInverseOf(compliance);
            answer = elasticity;
            answer.times(1.0 - tempDam);
        }
    }

    double g = status->giveDensG();
    if ( g <= 0 ) {
        double factor = gammaL0 * pow(rho, rL) + gammaP0 *pow(rho, rP) * ( tDens - 1 ) * pow(g, tDens - 2);
        tangentMatrix.resize(6, 6);
        tangentMatrix.zero();
        tangentMatrix.at(1, 1) = tangentMatrix.at(1, 2) = tangentMatrix.at(1, 3) = 1;
        tangentMatrix.at(2, 1) = tangentMatrix.at(2, 2) = tangentMatrix.at(2, 3) = 1;
        tangentMatrix.at(3, 1) = tangentMatrix.at(3, 2) = tangentMatrix.at(3, 3) = 1;
        tangentMatrix.times(factor);
        answer.add(tangentMatrix);
    }

    status->setSmtrx(answer);
}




void
TrabBoneGrad3D :: give3dKappaMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(1, 6);
    answer.zero();
    if ( mode == TangentStiffness ) {
        TrabBoneGrad3DStatus *status = static_cast< TrabBoneGrad3DStatus * >( this->giveStatus(gp) );
        double beta;
        FloatArray plasFlowDirec, prodTensor;
        FloatMatrix SSaTensor;

        double kappa = status->giveKappa();
        double tempKappa = status->giveTempKappa();
        double dKappa = tempKappa - kappa;

        if ( dKappa > 0.0 ) {
            plasFlowDirec = status->givePlasFlowDirec();
            SSaTensor = status->giveSSaTensor();
            beta = status->giveBeta();
            prodTensor.beTProductOf(SSaTensor, plasFlowDirec);
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
        TrabBoneGrad3DStatus *status = static_cast< TrabBoneGrad3DStatus * >( this->giveStatus(gp) );
        double damage, tempDamage;
        double nlKappa, kappa;
        FloatArray tempEffStress;
        double gPrime;
        double tempKappa = status->giveTempKappa();
        damage = status->giveDam();
        nlKappa =  status->giveNonlocalCumulatedStrain();
        kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
        tempDamage = TrabBone3D :: computeDamageParam(kappa);


        if ( ( tempDamage - damage ) > 0 ) {
            tempEffStress =  status->giveTempEffectiveStress();
            for ( int i = 1; i <= 6; i++ ) {
                answer.at(i, 1) = tempEffStress.at(i);
            }

            kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
            gPrime = TrabBone3D :: computeDamageParamPrime(kappa);
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
TrabBoneGrad3D :: giveRealStressVectorGrad(FloatArray &answer1, double &answer2, GaussPoint *gp, const FloatArray &totalStrain, double nonlocalCumulatedStrain, TimeStep *tStep)
{
    TrabBoneGrad3DStatus *status = static_cast< TrabBoneGrad3DStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);

    TrabBone3D :: performPlasticityReturn(gp, totalStrain, tStep);

    double tempDamage = computeDamage(gp, tStep);
    FloatArray tempEffStress = status->giveTempEffectiveStress();
    answer1.beScaled(1 - tempDamage, tempEffStress);
    answer2 = status->giveTempKappa();


    if ( densCrit != 0 ) {
        FloatArray densStress;
        computeDensificationStress(densStress, gp, totalStrain, tStep);
        answer1.add(densStress);
    }


    status->letTempStrainVectorBe(totalStrain);
    status->setNonlocalCumulatedStrain(nonlocalCumulatedStrain);

    status->setTempDam(tempDamage);
    status->setTempEffectiveStress(tempEffStress);
    status->letTempStressVectorBe(answer1);
}


void
TrabBoneGrad3D :: computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *tStep)
{
    TrabBoneGrad3DStatus *status = static_cast< TrabBoneGrad3DStatus * >( this->giveStatus(gp) );
    double localCumPlastStrain = status->giveTempKappa();
    double nlCumPlastStrain = status->giveNonlocalCumulatedStrain();
    kappa = mParam * nlCumPlastStrain + ( 1 - mParam ) * localCumPlastStrain;
}


IRResultType
TrabBoneGrad3D :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                             // Required by IR_GIVE_FIELD macro

    IR_GIVE_OPTIONAL_FIELD(ir, L, _IFT_TrabBoneGrad3D_L);
    if ( L < 0.0 ) {
        L = 0.0;
    }

    mParam = 2.;
    IR_GIVE_OPTIONAL_FIELD(ir, mParam, _IFT_TrabBoneGrad3D_m);

    return TrabBone3D :: initializeFrom(ir);
}


TrabBoneGrad3DStatus :: TrabBoneGrad3DStatus(int n, Domain *d, GaussPoint *g) :
    TrabBone3DStatus(n, d, g)
{
    nlKappa = 0;
}


TrabBoneGrad3DStatus :: ~TrabBoneGrad3DStatus()
{ }


void
TrabBoneGrad3DStatus :: printOutputAt(FILE *file, TimeStep *tStep)
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
    nss = 1;
    //TrabBone3DStatus :: updateYourself(tStep);
}
} // end namespace oofem
