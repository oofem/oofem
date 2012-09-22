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

#include "druckerpragercutmat.h"
#include "isolinearelasticmaterial.h"
#include "gausspnt.h"
#include "intarray.h"
#include "stressvector.h"
#include "strainvector.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"

namespace oofem {
// constructor
DruckerPragerCutMat :: DruckerPragerCutMat(int n, Domain *d) : StructuralMaterial(n, d)
{
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
    sigT = 0.;
    H = 0.;
    omegaCrit = 0.;
    a=0.;
    yieldTol = 0.;
    newtonIter = 30;
    G = 0.;
    K = 0.;
}

// destructor
DruckerPragerCutMat :: ~DruckerPragerCutMat()
{ }

// specifies whether a given material mode is supported by this model
int
DruckerPragerCutMat :: hasMaterialModeCapability(MaterialMode mode)
{
    if ( ( mode == _3dMat ) || ( mode == _1dMat ) || ( mode == _PlaneStrain ) ) {
        return 1;
    }

    return 0;
}

// reads the model parameters from the input file
IRResultType
DruckerPragerCutMat :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // required by IR_GIVE_FIELD macro
    IRResultType result;                 // required by IR_GIVE_FIELD macro

    StructuralMaterial :: initializeFrom(ir);
    linearElasticMaterial->initializeFrom(ir); // takes care of elastic constants
    G = ( ( IsotropicLinearElasticMaterial * ) linearElasticMaterial )->giveShearModulus();
    K = ( ( IsotropicLinearElasticMaterial * ) linearElasticMaterial )->giveBulkModulus();

    IR_GIVE_FIELD(ir, tau0, IFT_DruckerPragerCutMat_tau0, "tau0"); // initial yield stress under pure shear
    IR_GIVE_FIELD(ir, alpha, IFT_DruckerPragerCutMat_alpha, "alpha"); // friction coefficient
    alphaPsi=alpha;
    IR_GIVE_OPTIONAL_FIELD(ir, alphaPsi, IFT_DruckerPragerCutMat_alphapsi, "alphapsi"); //dilatancy coefficient
    IR_GIVE_OPTIONAL_FIELD(ir, H, IFT_DruckerPragerCutMat_h, "h"); // hardening modulus
    IR_GIVE_OPTIONAL_FIELD(ir, sigT, IFT_DruckerPragerCutMat_sigT, "sigt"); // uniaxial tensile strength for cut-off
    IR_GIVE_OPTIONAL_FIELD(ir, omegaCrit, IFT_DruckerPragerCutMat_omegaCrit, "omega_crit"); // critical damage
    IR_GIVE_OPTIONAL_FIELD(ir, a, IFT_DruckerPragerCutMat_a, "a"); // exponent in damage law
    IR_GIVE_OPTIONAL_FIELD(ir, yieldTol, IFT_DruckerPragerCutMat_yieldTol, "yieldtol"); //tolerance of the error in the yield criterion
    IR_GIVE_OPTIONAL_FIELD(ir, newtonIter, IFT_DruckerPragerCutMat_newtonIter, "yieldtol"); //tolerance of the error in the yield criterion
    
    
    return IRRT_OK;
}

// creates a new material status corresponding to this class
MaterialStatus *
DruckerPragerCutMat :: CreateStatus(GaussPoint *gp) const
{
    DruckerPragerCutMatStatus *status;
    status = new DruckerPragerCutMatStatus(1, this->giveDomain(), gp);
    return status;
}


// returns the stress vector in 3d stress space
void
DruckerPragerCutMat :: giveRealStressVector(FloatArray &answer,
                                 MatResponseForm form,
                                 GaussPoint *gp,
                                 const FloatArray &totalStrain,
                                 TimeStep *atTime)
{
    MaterialMode mode = gp->giveMaterialMode();


    if ( ( mode == _3dMat ) || ( mode == _1dMat ) || ( mode == _PlaneStrain ) ) {
        giveRealStressVectorComputedFromStrain(answer, form, gp, totalStrain, atTime);
    } else {
        OOFEM_ERROR("MisesMat::giveRealStressVector : unknown material response mode");
    }
}

// returns the stress vector in 3d stress space
// computed from the previous plastic strain and current total strain
void
DruckerPragerCutMat :: giveRealStressVectorComputedFromStrain(FloatArray &answer,
                                                   MatResponseForm form,
                                                   GaussPoint *gp,
                                                   const FloatArray &totalStrain,
                                                   TimeStep *atTime)
{
    FloatArray reducedTotalStrainVector;
    
    DruckerPragerCutMatStatus *status = ( DruckerPragerCutMatStatus * ) this->giveStatus(gp);
    MaterialMode mode = gp->giveMaterialMode();
    this->initTempStatus(gp);
    this->initGpForNewStep(gp);
    // subtract stress independent part (temperature etc.)
    this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, atTime, VM_Total);
    this->performPlasticityReturn(gp, reducedTotalStrainVector, mode);
    double omega = computeDamage(gp, atTime);
    StressVector effStress(mode), totalStress(mode);
    status->giveTempEffectiveStress(effStress);
    totalStress = effStress;
    totalStress.times(1 - omega);
    answer = totalStress;
    status->setTempDamage(omega);
    status->letTempStrainVectorBe(totalStrain);
    status->letTempStressVectorBe(answer);
}

void DruckerPragerCutMat :: performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain, MaterialMode mode)
{
    double kappa, yieldValue, dKappa;
    FloatArray reducedStress;
    FloatArray strain, plStrain;
    //  MaterialMode mode = gp->giveMaterialMode();
    DruckerPragerCutMatStatus *status = ( DruckerPragerCutMatStatus * ) this->giveStatus(gp);
    StressVector fullStress(mode);
    this->initTempStatus(gp);
    this->initGpForNewStep(gp);
    // get the initial plastic strain and initial kappa from the status
    status->givePlasticStrain(plStrain);
    kappa = status->giveCumulativePlasticStrain();

    // === radial return algorithm ===
    if ( mode == _1dMat ) {
        LinearElasticMaterial *lmat = this->giveLinearElasticMaterial();
        double E = lmat->give('E', gp);
        /*trial stress*/
        fullStress.at(1) = E * ( totalStrain.at(1) - plStrain.at(1) );
        double tauY = tau0 + H * kappa;
        double trialS = fullStress.at(1);
        if(trialS > sigT) {//tension cut-off active
            status->setActiveSurface(1);
        } else if ( alpha*trialS + fabs(trialS)/sqrt(3.) - tauY > 0 ){//Drucker-Prager model active
            status->setActiveSurface(2);
            
        } else {//elastic
            status->setActiveSurface(0);
        }
        
        
        
        trialS = fabs(trialS);
        /*yield function*/
        yieldValue = trialS - tauY;
        // === radial return algorithm ===
        if ( yieldValue > 0 ) {
            dKappa = yieldValue / ( H + E );
            kappa += dKappa;
            fullStress.at(1) = fullStress.at(1) - dKappa *E *signum( fullStress.at(1) );
            plStrain.at(1) = plStrain.at(1) + dKappa *signum( fullStress.at(1) );
        }
    } else if ( ( mode == _PlaneStrain ) || ( mode == _3dMat ) ) {
        // elastic predictor
        StrainVector elStrain(totalStrain, mode);
        elStrain.subtract(plStrain);
        StrainVector elStrainDev(mode);
        double elStrainVol;
        elStrain.computeDeviatoricVolumetricSplit(elStrainDev, elStrainVol);
        StressVector trialStressDev(mode);
        elStrainDev.applyDeviatoricElasticStiffness(trialStressDev, G);
        /**************************************************************/
        double trialStressVol;
        trialStressVol = 3 * K * elStrainVol;
        /**************************************************************/

        // store the deviatoric and trial stress (reused by algorithmic stiffness)
        status->letTrialStressDevBe(trialStressDev);
        status->setTrialStressVol(trialStressVol);
        // check the yield condition at the trial state
        double trialS = trialStressDev.computeStressNorm();
        double tauY = tau0 + H * kappa;
        yieldValue = sqrt(3. / 2.) * trialS - tauY;
        if ( yieldValue > 0. ) {
            // increment of cumulative plastic strain
            dKappa = yieldValue / ( H + 3. * G );
            kappa += dKappa;
            StrainVector dPlStrain(mode);
            // the following line is equivalent to multiplication by scaling matrix P
            trialStressDev.applyDeviatoricElasticCompliance(dPlStrain, 0.5);
            // increment of plastic strain
            dPlStrain.times(sqrt(3. / 2.) * dKappa / trialS);
            plStrain.add(dPlStrain);
            // scaling of deviatoric trial stress
            trialStressDev.times(1. - sqrt(6.) * G * dKappa / trialS);
        }

        // assemble the stress from the elastically computed volumetric part
        // and scaled deviatoric part

        double stressVol = 3. * K * elStrainVol;
        trialStressDev.computeDeviatoricVolumetricSum(fullStress, stressVol);
    }

    // store the effective stress in status
    status->letTempEffectiveStressBe(fullStress);
    // store the plastic strain and cumulative plastic strain
    status->letTempPlasticStrainBe(plStrain);
    status->setTempCumulativePlasticStrain(kappa);
}


double
DruckerPragerCutMat :: computeDamageParam(double tempKappa)
{
    double tempDam;
    if ( tempKappa > 0. ) {
        tempDam = omegaCrit * ( 1.0 - exp(-a * tempKappa) );
    } else {
        tempDam = 0.;
    }

    return tempDam;
}

double
DruckerPragerCutMat :: computeDamageParamPrime(double tempKappa)
{
    double tempDam;
    if ( tempKappa >= 0. ) {
        tempDam = omegaCrit * a * exp(-a * tempKappa);
    } else {
        tempDam = 0.;
    }

    return tempDam;
}



double
DruckerPragerCutMat :: computeDamage(GaussPoint *gp,  TimeStep *atTime)
{
    double tempKappa, dam;
    DruckerPragerCutMatStatus *status = ( DruckerPragerCutMatStatus * ) this->giveStatus(gp);
    dam = status->giveDamage();
    computeCumPlastStrain(tempKappa, gp, atTime);
    double tempDam = computeDamageParam(tempKappa);
    if ( dam > tempDam ) {
        tempDam = dam;
    }

    return tempDam;
}


void DruckerPragerCutMat :: computeCumPlastStrain(double &tempKappa, GaussPoint *gp, TimeStep *atTime)
{
    DruckerPragerCutMatStatus *status = ( DruckerPragerCutMatStatus * ) this->giveStatus(gp);
    tempKappa = status->giveTempCumulativePlasticStrain();
}

void
DruckerPragerCutMat :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    MaterialMode mMode = gp->giveMaterialMode();

    if ( mMode == _3dMat ) {
        give3dSSMaterialStiffnessMatrix(answer, form, mode, gp, atTime);
    } else {
        OOFEM_ERROR("DruckerPragerCutMat::give3dMaterialStiffnessMatrix : unknown material response mode");
    }
}


// returns the consistent (algorithmic) tangent stiffness matrix
void
DruckerPragerCutMat :: give3dSSMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form,
                                            MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *atTime)
{
    // start from the elastic stiffness
    this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    if ( mode != TangentStiffness ) {
        return;
    }

    DruckerPragerCutMatStatus *status = ( DruckerPragerCutMatStatus * ) this->giveStatus(gp);
    double kappa = status->giveCumulativePlasticStrain();
    double tempKappa = status->giveTempCumulativePlasticStrain();
    // increment of cumulative plastic strain as an indicator of plastic loading
    double dKappa = tempKappa - kappa;

    if ( dKappa <= 0.0 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
        return;
    }

    // === plastic loading ===

    // yield stress at the beginning of the step
    double tauY = tau0 + H * kappa;

    // trial deviatoric stress and its norm
    StressVector trialStressDev(_3dMat);
    double trialStressVol;
    status->giveTrialStressVol(trialStressVol);
    status->giveTrialStressDev(trialStressDev);
    double trialS = trialStressDev.computeStressNorm();

    // one correction term
    FloatMatrix stiffnessCorrection(6, 6);
    stiffnessCorrection.beDyadicProductOf(trialStressDev, trialStressDev);
    double factor = -2. * sqrt(6.) * G * G / trialS;
    double factor1 = factor * tauY / ( ( H + 3. * G ) * trialS * trialS );
    stiffnessCorrection.times(factor1);
    answer.add(stiffnessCorrection);

    // another correction term
    stiffnessCorrection.bePinvID();
    double factor2 = factor * dKappa;
    stiffnessCorrection.times(factor2);
    answer.add(stiffnessCorrection);

    //influence of damage
    //    double omega = computeDamageParam(tempKappa);
    double omega = status->giveTempDamage();
    answer.times(1. - omega);
    FloatArray effStress;
    status->giveTempEffectiveStress(effStress);
    double omegaPrime = computeDamageParamPrime(tempKappa);
    double scalar = -omegaPrime *sqrt(6.) * G / ( 3. * G + H ) / trialS;
    stiffnessCorrection.beDyadicProductOf(effStress, trialStressDev);
    stiffnessCorrection.times(scalar);
    answer.add(stiffnessCorrection);
}

void
DruckerPragerCutMat :: give1dStressStiffMtrx(FloatMatrix &answer, MatResponseForm form,
                                  MatResponseMode mode,
                                  GaussPoint *gp,
                                  TimeStep *atTime)
{
    this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    FloatArray stressVector;
    DruckerPragerCutMatStatus *status = ( DruckerPragerCutMatStatus * ) this->giveStatus(gp);
    double kappa = status->giveCumulativePlasticStrain();
    // increment of cumulative plastic strain as an indicator of plastic loading
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double omega = status->giveTempDamage();
    double E = answer.at(1, 1);
    if ( mode != TangentStiffness ) {
        return;
    }


    if ( ( tempKappa - kappa ) <= 0.0 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
        answer.times(1 - omega);
        return;
    }


    // === plastic loading ===
    status->giveTempEffectiveStress(stressVector);
    double stress = stressVector.at(1);
    answer.resize(1, 1);
    answer.at(1, 1) = ( 1 - omega ) * E * H / ( E + H ) - computeDamageParamPrime(tempKappa) * E / ( E + H ) * stress * signum(stress);
}

void
DruckerPragerCutMat :: givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseForm form,
                                     MatResponseMode mode,
                                     GaussPoint *gp,
                                     TimeStep *atTime)
{
    this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    if ( mode != TangentStiffness ) {
        return;
    }

    DruckerPragerCutMatStatus *status = ( DruckerPragerCutMatStatus * ) this->giveStatus(gp);
    double tempKappa = status->giveTempCumulativePlasticStrain();
    double kappa = status->giveCumulativePlasticStrain();
    double dKappa = tempKappa - kappa;

    if ( dKappa <= 0.0 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
        return;
    }

    // === plastic loading ===
    // yield stress at the beginning of the step
    double tauY = tau0 + H * kappa;

    // trial deviatoric stress and its norm
    StressVector trialStressDev(_PlaneStrain);
    double trialStressVol;
    status->giveTrialStressVol(trialStressVol);
    status->giveTrialStressDev(trialStressDev);
    double trialS = trialStressDev.computeStressNorm();

    // one correction term
    FloatMatrix stiffnessCorrection(4, 4);
    stiffnessCorrection.beDyadicProductOf(trialStressDev, trialStressDev);
    double factor = -2. * sqrt(6.) * G * G / trialS;
    double factor1 = factor * tauY / ( ( H + 3. * G ) * trialS * trialS );
    stiffnessCorrection.times(factor1);
    answer.add(stiffnessCorrection);

    // another correction term
    stiffnessCorrection.resize(4, 4);
    stiffnessCorrection.zero();
    stiffnessCorrection.at(1, 1) = stiffnessCorrection.at(2, 2) = stiffnessCorrection.at(3, 3) = 2. / 3.;
    stiffnessCorrection.at(1, 2) = stiffnessCorrection.at(1, 3) = stiffnessCorrection.at(2, 1) = -1. / 3.;
    stiffnessCorrection.at(2, 3) = stiffnessCorrection.at(3, 1) = stiffnessCorrection.at(3, 2) = -1. / 3.;
    stiffnessCorrection.at(4, 4) = 0.5;
    double factor2 = factor * dKappa;
    stiffnessCorrection.times(factor2);
    answer.add(stiffnessCorrection);

    //influence of damage
    double omega = status->giveTempDamage();
    answer.times(1. - omega);
    FloatArray effStress;
    status->giveTempEffectiveStress(effStress);
    double omegaPrime = computeDamageParamPrime(tempKappa);
    double scalar = -omegaPrime *sqrt(6.) * G / ( 3. * G + H ) / trialS;
    stiffnessCorrection.beDyadicProductOf(effStress, trialStressDev);
    stiffnessCorrection.times(scalar);
    answer.add(stiffnessCorrection);
}

#ifdef __OOFEG
#endif

int
DruckerPragerCutMat :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    DruckerPragerCutMatStatus *status = ( DruckerPragerCutMatStatus * ) this->giveStatus(aGaussPoint);
    if ( type == IST_PlasticStrainTensor ) {
        answer  = * status->givePlasDef();
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = status->giveCumulativePlasticStrain();
        return 1;
    } else if ( ( type == IST_DamageScalar ) || (type == IST_DamageTensor ) ) {
        answer.resize(1);
        answer.at(1) = status->giveDamage();
        return 1;
    } else {
        return StructuralMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}




InternalStateValueType
DruckerPragerCutMat :: giveIPValueType(InternalStateType type)
{
    if ( type == IST_PlasticStrainTensor ) {
        return ISVT_TENSOR_S3;
    } else if ( ( type == IST_MaxEquivalentStrainLevel ) || ( type == IST_DamageScalar ) ) {
        return ISVT_SCALAR;
    } else {
        return StructuralMaterial :: giveIPValueType(type);
    }
}


int
DruckerPragerCutMat :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( type == IST_PlasticStrainTensor ) {
        if ( ( mmode == _3dMat ) || ( mmode == _3dMat_F ) ) {
            answer.resize(6);
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(3) = 3;
            answer.at(4) = 4;
            answer.at(5) = 5;
            answer.at(6) = 6;
            return 1;
        } else if ( mmode == _PlaneStrain ) {
            answer.resize(6);
            answer.at(1) = 1;
            answer.at(2) = 2;
            answer.at(3) = 3;
            answer.at(6) = 4;
            return 1;
        } else if ( mmode == _1dMat ) {
            answer.resize(1);
            answer.at(1) = 1;
            return 1;
        }
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else if ( ( type == IST_DamageScalar ) || ( type == IST_DamageTensor) ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    }

    return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
}


int
DruckerPragerCutMat :: giveIPValueSize(InternalStateType type, GaussPoint *gp)
{
    if ( type == IST_PlasticStrainTensor ) {
        MaterialMode mode = gp->giveMaterialMode();
        if (mode == _3dMat || mode == _3dMat_F)
            return 6;
        else if (mode == _PlaneStrain)
            return 4;
        else if (mode == _PlaneStress)
            return 3;
        else if (mode == _1dMat)
            return 1;
        else
            return 0;
    } else if ( type == IST_MaxEquivalentStrainLevel ) {
        return 1;
    } else if ( type == IST_DamageScalar || type == IST_DamageTensor) {
        return 1;
    } else {
        return StructuralMaterial :: giveIPValueSize(type, gp);
    }
}






//=============================================================================

DruckerPragerCutMatStatus :: DruckerPragerCutMatStatus (int n, Domain *d, GaussPoint *g) :
    StructuralMaterialStatus(n, d, g), plasticStrain(), tempPlasticStrain(), trialStressD()
{
    damage = tempDamage = 0.;
    kappa = tempKappa = 0.;
    effStress.resize(6);
    tempEffStress.resize(6);
}

DruckerPragerCutMatStatus :: ~DruckerPragerCutMatStatus ()
{ }

void
DruckerPragerCutMatStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    //int i, n;

    StructuralMaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, "status { ");
    /*
     * // this would not be correct, since the status is already updated and kappa is now the "final" value
     * if ( tempKappa > kappa ) {
     * fprintf(file, " Yielding, ");
     * } else {
     * fprintf(file, " Unloading, ");
     * }
     */

    /*********************************************************************************/

    fprintf(file, "damage %.4e", tempDamage);

    /********************************************************************************/
    /*
     * // print the plastic strain
     *  n = plasticStrain.giveSize();
     *  fprintf(file, " plastic_strains ");
     *  for ( i = 1; i <= n; i++ ) {
     *      fprintf( file, " % .4e", plasticStrain.at(i) );
     *  }
     */
    // print the cumulative plastic strain
    fprintf(file, ", kappa ");
    fprintf(file, " % .4e", kappa);

    fprintf(file, "}\n");
 
    fprintf(file, "}\n");
}


// initializes temporary variables based on their values at the previous equlibrium state
void DruckerPragerCutMatStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();

    if ( plasticStrain.giveSize() == 0 ) {
        plasticStrain.resize( ( ( StructuralMaterial * ) gp->giveMaterial() )->giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() ) );
        plasticStrain.zero();
    }

    tempDamage = damage;
    tempPlasticStrain = plasticStrain;
    tempKappa = kappa;
    trialStressD.resize(0); // to indicate that it is not defined yet
}


// updates internal variables when equilibrium is reached
void
DruckerPragerCutMatStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);

    plasticStrain = tempPlasticStrain;
    kappa = tempKappa;
    damage = tempDamage;
    trialStressD.resize(0); // to indicate that it is not defined any more

}


// saves full information stored in this status
// temporary variables are NOT stored
contextIOResultType
DruckerPragerCutMatStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data

    // write plastic strain (vector)
    if ( ( iores = plasticStrain.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write cumulative plastic strain (scalar)
    if ( !stream->write(& kappa, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // write damage (scalar)
    if ( !stream->write(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}



contextIOResultType
DruckerPragerCutMatStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full information stored in stream to this Status
//
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read plastic strain (vector)
    if ( ( iores = plasticStrain.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read cumulative plastic strain (scalar)
    if ( !stream->read(& kappa, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // read damage (scalar)
    if ( !stream->read(& damage, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK; // return succes
}
} // end namespace oofem
