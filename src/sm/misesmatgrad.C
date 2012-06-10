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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "misesmatgrad.h"
#include "stressvector.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "mathfem.h"
#include "error.h"


namespace oofem {
/////////////////////////////////////////////////////////////////
//gradient regularization of Mises plasticity coupled with isotropic damage////////////////////
/////////////////////////////////////////////////////////////////
double sig(double number) {
    if ( number < 0 ) {
        return -1;
    } else if ( number > 0 ) {
        return 1;
    } else {
        return 0;
    }
}

MisesMatGrad :: MisesMatGrad(int n, Domain *d) : MisesMat(n, d)
{
    R = 0.;
}


MisesMatGrad :: ~MisesMatGrad()
{}


int
MisesMatGrad :: hasMaterialModeCapability(MaterialMode mode)
{
    if ( ( mode == _1dMatGrad ) || ( mode == _PlaneStrainGrad ) || ( mode == _3dMatGrad ) ) {
        return 1;
    }

    return 0;
}


void
MisesMatGrad :: giveCharacteristicMatrix(FloatMatrix &answer,
                                         MatResponseForm form, MatResponseMode rMode, GaussPoint *gp, TimeStep *atTime)
//
// Returns characteristic material stiffness matrix of the receiver
//
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch ( mMode ) {
    case _1dMatGrad:
        if ( form == PDGrad_uu ) {
            give1dStressStiffMtrx(answer, form, rMode, gp, atTime);
            break;
        } else if ( form == PDGrad_ku ) {
            give1dKappaMatrix(answer, form, rMode, gp, atTime);
            break;
        } else if ( form == PDGrad_uk ) {
            give1dGprime(answer, form, rMode, gp, atTime);
            break;
        } else if ( form == PDGrad_kk ) {
            giveInternalLength(answer, form, rMode, gp, atTime);
            break;
        }

    case _PlaneStrainGrad:
        if ( form == PDGrad_uu ) {
            givePlaneStrainStiffMtrx(answer, form, rMode, gp, atTime);
            break;
        } else if ( form == PDGrad_ku ) {
            givePlaneStrainKappaMatrix(answer, form, rMode, gp, atTime);
            break;
        } else if ( form == PDGrad_uk ) {
            givePlaneStrainGprime(answer, form, rMode, gp, atTime);
            break;
        } else if ( form == PDGrad_kk ) {
            giveInternalLength(answer, form, rMode, gp, atTime);
            break;
        }

    case _PlaneStressGrad:
        _error2( "giveCharacteristicMatrix : unknown mode (%s)", __MaterialModeToString(mMode) );
        break;

    case _3dMatGrad:
        if ( form == PDGrad_uu ) {
            give3dMaterialStiffnessMatrix(answer, form, rMode, gp, atTime);
            break;
        } else if ( form == PDGrad_ku ) {
            give3dKappaMatrix(answer, form, rMode, gp, atTime);
            break;
        } else if ( form == PDGrad_uk ) {
            give3dGprime(answer, form, rMode, gp, atTime);
            break;
        } else if ( form == PDGrad_kk ) {
            giveInternalLength(answer, form, rMode, gp, atTime);
            break;
        }





    default:
        _error2( "giveCharacteristicMatrix : unknown mode (%s)", __MaterialModeToString(mMode) );
        return;
    }
}


void
MisesMatGrad :: give1dStressStiffMtrx(FloatMatrix &answer, MatResponseForm form,
                                      MatResponseMode mode,
                                      GaussPoint *gp,
                                      TimeStep *atTime)
{
    answer.resize(1, 1);
    answer.zero();
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    double E = lmat->give('E', gp);
    answer.at(1, 1) = E;
    if ( mode != TangentStiffness ) {
        return;
    }

    FloatArray stressVector;
    MisesMatStatus *status = ( MisesMatStatus * ) this->giveStatus(gp);
    double tempKappa = status->giveTempCumulativePlasticStrain();
    // increment of cumulative plastic strain as an indicator of plastic loading
    double dKappa = ( tempKappa - status->giveCumulativePlasticStrain() );
    /********************************************************************/
    double tempDamage = status->giveTempDamage();
    double damage = status->giveDamage();
    /*********************************************************************/
    double nlKappa =  status->giveTempStrainVector().at(2);
    double kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
    if ( dKappa <= 0.0 ) {
        answer.at(1, 1) = ( 1. - tempDamage ) * E;
        return;
    }


    // === plastic loading ===
    status->giveTempEffectiveStress(stressVector);
    double stress = stressVector.at(1);

    answer.at(1, 1) = ( 1. - tempDamage ) * E * H / ( E + H );
    if ( ( tempDamage - damage ) > 0 ) {
        answer.at(1, 1) = answer.at(1, 1) - ( 1. - mParam ) * computeDamageParamPrime(kappa) * E / ( E + H ) * sig(stress) * stress;
    }
}
void
MisesMatGrad :: givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    if ( mode != TangentStiffness ) {
        return;
    }

    MisesMatStatus *status = ( MisesMatStatus * ) this->giveStatus(gp);
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double kappa = status->giveCumulativePlasticStrain();
    double dKappa = tempKappa - kappa;
    double tempDamage = status->giveTempDamage();
    double damage = status->giveDamage();

    if ( dKappa <= 0.0 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
        return;
    }

    // === plastic loading ===
    // yield stress at the beginning of the step
    double sigmaY = sig0 + H * kappa;
    // trial deviatoric stress and its norm
    StressVector trialStressDev(_PlaneStrain);
    status->giveTrialStressDev(trialStressDev);
    double trialS = trialStressDev.computeStressNorm();
    // volumetric stress
    double trialStressVol;
    status->giveTrialStressVol(trialStressVol);
    // one correction term
    FloatMatrix stiffnessCorrection(4, 4);
    stiffnessCorrection.beDyadicProductOf(trialStressDev, trialStressDev);
    double factor = -2. * sqrt(6.) * G * G / trialS;
    double factor1 = factor * sigmaY / ( ( H + 3. * G ) * trialS * trialS );
    stiffnessCorrection.times(factor1);
    answer.add(stiffnessCorrection);
    // another correction term
    stiffnessCorrection.zero();
    stiffnessCorrection.at(1, 1) = stiffnessCorrection.at(2, 2) = stiffnessCorrection.at(3, 3) = 2. / 3.;
    stiffnessCorrection.at(1, 2) = stiffnessCorrection.at(1, 3) = stiffnessCorrection.at(2, 1) = -1. / 3.;
    stiffnessCorrection.at(2, 3) = stiffnessCorrection.at(3, 1) = stiffnessCorrection.at(3, 2) = -1. / 3.;
    stiffnessCorrection.at(4, 4) = 0.5;
    double factor2 = factor * dKappa;
    stiffnessCorrection.times(factor2);
    answer.add(stiffnessCorrection);
    //influence of damage
    answer.times(1 - tempDamage);
    if ( tempDamage > damage ) {
        FloatArray effStress, totalStress;
        totalStress = status->giveTempStressVector();
        status->giveTempEffectiveStress(effStress);
        double nlKappa = status->giveTempStrainVector().at(5);
        double kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
        double omegaPrime = computeDamageParamPrime(kappa);
        double scalar = -omegaPrime *sqrt(6.) * G / ( 3. * G + H ) / trialS;
        stiffnessCorrection.beDyadicProductOf(effStress, trialStressDev);
        stiffnessCorrection.times( scalar * ( 1. - mParam ) );
        answer.add(stiffnessCorrection);
    }
}


void
MisesMatGrad :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    MisesMatGradStatus *status = ( MisesMatGradStatus * ) this->giveStatus(gp);
    this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    // start from the elastic stiffness
    if ( mode != TangentStiffness ) {
        return;
    }

    double tempKappa = status->giveTempCumulativePlasticStrain();
    double kappa = status->giveCumulativePlasticStrain();
    double dKappa = tempKappa - kappa;

    if ( dKappa > 0.0 ) {
        double tempDamage = status->giveTempDamage();
        double damage = status->giveDamage();
        double sigmaY = sig0 + H * kappa;
        // trial deviatoric stress and its norm
        StressVector trialStressDev(_3dMat);
        /*****************************************************/
        double trialStressVol;
        status->giveTrialStressVol(trialStressVol);
        /****************************************************/
        status->giveTrialStressDev(trialStressDev);
        double trialS = trialStressDev.computeStressNorm();
        // one correction term
        FloatMatrix stiffnessCorrection(6, 6);
        stiffnessCorrection.beDyadicProductOf(trialStressDev, trialStressDev);
        double factor = -2. * sqrt(6.) * G * G / trialS;
        double factor1 = factor * sigmaY / ( ( H + 3. * G ) * trialS * trialS );
        stiffnessCorrection.times(factor1);
        answer.add(stiffnessCorrection);
        // another correction term
        stiffnessCorrection.bePinvID();
        double factor2 = factor * dKappa;
        stiffnessCorrection.times(factor2);
        answer.add(stiffnessCorrection);
        //influence of damage
        answer.times(1. - tempDamage);
        if ( tempDamage > damage ) {
            FloatArray effStress, totalStress;
            totalStress = status->giveTempStressVector();
            status->giveTempEffectiveStress(effStress);
            double nlKappa =  status->giveTempStrainVector().at(7);
            double kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
            double omegaPrime = computeDamageParamPrime(kappa);
            double scalar = -omegaPrime *sqrt(6.) * G / ( 3. * G + H ) / trialS;
            stiffnessCorrection.beDyadicProductOf(effStress, trialStressDev);
            stiffnessCorrection.times( scalar * ( 1. - mParam ) );
            answer.add(stiffnessCorrection);
        }
    }
}


void
MisesMatGrad :: give1dKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    answer.resize(1, 1);
    answer.zero();
    LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
    double E = lmat->give('E', gp);
    MisesMatGradStatus *status = ( MisesMatGradStatus * ) this->giveStatus(gp);
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double dKappa = tempKappa - status->giveCumulativePlasticStrain();
    FloatArray effStress(6);
    status->giveTempEffectiveStress(effStress);
    double stress = effStress.at(1);
    if ( dKappa > 0 ) {
        double trialS = sig(stress);
        double factor = trialS * E / ( E + H );
        answer.at(1, 1) = factor;
    }
}


void
MisesMatGrad :: givePlaneStrainKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    MisesMatGradStatus *status = ( MisesMatGradStatus * ) this->giveStatus(gp);
    answer.resize(1, 4);
    answer.zero();
    StressVector trialStressDev(_PlaneStrain);
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double dKappa = tempKappa - status->giveCumulativePlasticStrain();
    status->giveTrialStressDev(trialStressDev);
    if ( dKappa > 0 ) {
        double trialS = trialStressDev.computeStressNorm();
        answer.at(1, 1) = trialStressDev.at(1);
        answer.at(1, 2) = trialStressDev.at(2);
        answer.at(1, 3) = trialStressDev.at(3);
        answer.at(1, 4) = trialStressDev.at(4);
        double factor = sqrt(6.) * G / ( 3. * G + H ) / trialS;
        answer.times(factor);
    }
}


void
MisesMatGrad :: give3dKappaMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    MisesMatGradStatus *status = ( MisesMatGradStatus * ) this->giveStatus(gp);
    answer.resize(1, 6);
    StressVector trialStressDev(_3dMat);
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double dKappa = tempKappa - status->giveCumulativePlasticStrain();
    status->giveTrialStressDev(trialStressDev);
    if ( dKappa > 0 ) {
        double trialS = trialStressDev.computeStressNorm();
        for ( int i = 1; i <= 6; i++ ) {
            answer.at(1, i) = trialStressDev.at(i);
        }

        double factor = sqrt(6.) * G / ( 3. * G + H ) / trialS;
        answer.times(factor);
    }
}



void
MisesMatGrad :: give1dGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    MisesMatGradStatus *status = ( MisesMatGradStatus * ) this->giveStatus(gp);
    double damage, tempDamage;
    double nlKappa, kappa;
    double tempKappa = status->giveTempCumulativePlasticStrain();
    FloatArray tempEffStress, strain;
    double gPrime;
    answer.resize(1, 1);
    damage = status->giveDamage();
    tempDamage = status->giveTempDamage();
    nlKappa =  status->giveTempStrainVector().at(2);
    kappa = mParam * nlKappa + ( 1 - mParam ) * tempKappa;
    status->giveTempEffectiveStress(tempEffStress);
    if ( ( tempDamage - damage ) > 0 ) {
        answer.at(1, 1) = tempEffStress.at(1);
        gPrime = computeDamageParamPrime(kappa);
        answer.times(gPrime * mParam);
    } else {
        answer.zero();
    }
}
void
MisesMatGrad :: givePlaneStrainGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    MisesMatGradStatus *status = ( MisesMatGradStatus * ) this->giveStatus(gp);
    double damage, tempDamage;
    double nlKappa, kappa;
    answer.resize(4, 1);
    answer.zero();
    double tempKappa = status->giveTempCumulativePlasticStrain();
    FloatArray tempEffStress;
    double gPrime;
    damage = status->giveDamage();
    tempDamage = status->giveTempDamage();
    nlKappa =  status->giveTempStrainVector().at(5);
    kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
    status->giveTempEffectiveStress(tempEffStress);
    if ( ( tempDamage - damage ) > 0 ) {
        answer.at(1, 1) = tempEffStress.at(1);
        answer.at(2, 1) = tempEffStress.at(2);
        answer.at(3, 1) = tempEffStress.at(3);
        answer.at(4, 1) = tempEffStress.at(4);
        gPrime = computeDamageParamPrime(kappa);
        answer.times(gPrime * mParam);
    } else {
        answer.zero();
    }
}
void
MisesMatGrad :: give3dGprime(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    MisesMatGradStatus *status = ( MisesMatGradStatus * ) this->giveStatus(gp);
    answer.resize(6, 1);
    answer.zero();
    double damage, tempDamage;
    double nlKappa, kappa;
    FloatArray tempEffStress;
    double gPrime;
    double tempKappa = status->giveTempCumulativePlasticStrain();
    damage = status->giveDamage();
    tempDamage = status->giveTempDamage();
    nlKappa =  status->giveTempStrainVector().at(7);
    status->giveTempEffectiveStress(tempEffStress);
    kappa = mParam * nlKappa + ( 1. - mParam ) * tempKappa;
    if ( ( tempDamage - damage ) > 0 ) {
        for ( int i = 1; i <= 6; i++ ) {
            answer.at(i, 1) = tempEffStress.at(i);
        }

        gPrime = computeDamageParamPrime(kappa);
        answer.times(gPrime * mParam);
    } else {
        answer.zero();
    }
}

void
MisesMatGrad :: giveInternalLength(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    answer.resize(1, 1);
    answer.at(1, 1) = R;
}

void
MisesMatGrad :: giveRealStressVector(FloatArray &answer, MatResponseForm form, GaussPoint *gp,
                                     const FloatArray &totalStrain, TimeStep *atTime)
{
    MisesMatGradStatus *status = ( MisesMatGradStatus * ) this->giveStatus(gp);
    this->initGpForNewStep(gp);
    this->initTempStatus(gp);
    MaterialMode mode = gp->giveMaterialMode();
    MaterialMode plReturnMode;
    if ( mode == _1dMatGrad ) {
        plReturnMode = _1dMat;
    } else if ( mode == _PlaneStrainGrad ) {
        plReturnMode = _PlaneStrain;
    } else if ( mode == _PlaneStressGrad ) {
        plReturnMode = _PlaneStress;
    } else if ( mode == _3dMatGrad ) {
        plReturnMode = _3dMat;
    }

    double tempDam;
    FloatArray tempEffStress, totalStress, locTotalStrain;

    int size = totalStrain.giveSize();
    locTotalStrain = totalStrain;
    locTotalStrain.resize(size - 1);
    MisesMat ::  performPlasticityReturn(gp, locTotalStrain, plReturnMode);
    status->letTempStrainVectorBe(totalStrain);
    double localCumPlastStrain = status->giveTempCumulativePlasticStrain();
    tempDam = computeDamage(gp, atTime);
    status->giveTempEffectiveStress(tempEffStress);
    answer.beScaled( 1.0 - tempDam, tempEffStress);
    size = tempEffStress.giveSize();
    answer.resize(size + 1);
    answer.at(size + 1) = localCumPlastStrain;

    status->setTempDamage(tempDam);
    status->letTempEffectiveStressBe(tempEffStress);
    status->letTempStressVectorBe(answer);
}


void
MisesMatGrad :: computeCumPlastStrain(double &kappa, GaussPoint *gp, TimeStep *atTime)
{
    double nlCumPlastStrain;
    MisesMatGradStatus *status = ( MisesMatGradStatus * ) this->giveStatus(gp);
    double localCumPlastStrain = status->giveTempCumulativePlasticStrain();
    FloatArray strain;
    strain = status->giveTempStrainVector();
    int size = strain.giveSize();
    nlCumPlastStrain = strain.at(size);
    kappa = mParam * nlCumPlastStrain + ( 1 - mParam ) * localCumPlastStrain;
}


IRResultType
MisesMatGrad :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                             // Required by IR_GIVE_FIELD macro

    MisesMat :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, R, IFT_MisesMatGrad_r, "r"); // Macro
    if ( R < 0.0 ) {
        R = 0.0;
    }

    mParam = 2.;
    IR_GIVE_OPTIONAL_FIELD(ir, mParam, IFT_MisesMatGrad_m, "m"); // Macro

    return IRRT_OK;
}


MisesMatGradStatus :: MisesMatGradStatus(int n, Domain *d, GaussPoint *g) :
    MisesMatStatus(n, d, g)
{}


MisesMatGradStatus :: ~MisesMatGradStatus()
{}


void
MisesMatGradStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    StructuralMaterialStatus :: printOutputAt(file, tStep);
    fprintf(file, "status {");
    fprintf(file, "kappa %f, damage %f ", this->kappa, this->damage);
    fprintf(file, "}\n");
}


void
MisesMatGradStatus :: initTempStatus()
{
    // MisesMatStatus::initTempStatus();
    StructuralMaterialStatus :: initTempStatus();

    if ( plasticStrain.giveSize() == 0 ) {
        if ( gp->giveMaterialMode() == _1dMatGrad ) {
            plasticStrain.resize( ( ( StructuralMaterial * ) gp->giveMaterial() )->giveSizeOfReducedStressStrainVector(_1dMat) );
        } else if ( gp->giveMaterialMode() == _PlaneStrainGrad ) {
            plasticStrain.resize( ( ( StructuralMaterial * ) gp->giveMaterial() )->giveSizeOfReducedStressStrainVector(_PlaneStrain) );
        } else if ( gp->giveMaterialMode() == _PlaneStressGrad ) {
            plasticStrain.resize( ( ( StructuralMaterial * ) gp->giveMaterial() )->giveSizeOfReducedStressStrainVector(_PlaneStress) );
        } else if ( gp->giveMaterialMode() == _3dMatGrad ) {
            plasticStrain.resize( ( ( StructuralMaterial * ) gp->giveMaterial() )->giveSizeOfReducedStressStrainVector(_3dMat) );
        }

        plasticStrain.zero();
    }

    tempDamage = damage;
    tempPlasticStrain = plasticStrain;
    tempKappa = kappa;
    tempDefGrad = defGrad;
    tempLeftCauchyGreen = leftCauchyGreen;
    trialStressD.resize(0);
}


void
MisesMatGradStatus :: updateYourself(TimeStep *atTime)
{
    MisesMatStatus :: updateYourself(atTime);
}
} // end namespace oofem
