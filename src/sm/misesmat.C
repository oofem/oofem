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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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
// file : misesmat.C


#include "misesmat.h"
#include "isolinearelasticmaterial.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "stressvector.h"
#include "strainvector.h"

#include "structuralcrosssection.h"
#include "mathfem.h"
#include "contextioerr.h"
#include "datastream.h"

namespace oofem {
// constructor
MisesMat :: MisesMat(int n, Domain *d) : StructuralMaterial(n, d)
{
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
    H = 0.;
    sig0 = 0.;
    G = 0.;
    K = 0.;
}

// destructor
MisesMat :: ~MisesMat()
{ }

// specifies whether a given material mode is supported by this model
int
MisesMat :: hasMaterialModeCapability(MaterialMode mode)
{
    if ( ( mode == _3dMat ) || ( mode == _3dMat_F ) ) {
        return 1;
    }

    return 0;
}

// reads the model parameters from the input file
IRResultType
MisesMat :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // required by IR_GIVE_FIELD macro
    IRResultType result;                 // required by IR_GIVE_FIELD macro

    StructuralMaterial :: initializeFrom(ir);
    linearElasticMaterial->initializeFrom(ir); // takes care of elastic constants
    G = ( ( IsotropicLinearElasticMaterial * ) linearElasticMaterial )->giveShearModulus();
    K = ( ( IsotropicLinearElasticMaterial * ) linearElasticMaterial )->giveBulkModulus();

    IR_GIVE_FIELD(ir, sig0, IFT_MisesMat_sig0, "sig0"); // uniaxial yield stress

    H = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, H, IFT_MisesMat_h, "h"); // hardening modulus


    return IRRT_OK;
}

// creates a new material status  corresponding to this class
MaterialStatus *
MisesMat :: CreateStatus(GaussPoint *gp) const
{
    MisesMatStatus *status;
    status = new MisesMatStatus(1, this->giveDomain(), gp);
    return status;
}

// returns the stress vector in 3d stress space
void
MisesMat :: giveRealStressVector(FloatArray &answer,
                                 MatResponseForm form,
                                 GaussPoint *gp,
                                 const FloatArray &totalStrain,
                                 TimeStep *atTime)
{
    MaterialMode mode = gp->giveMaterialMode();


    if ( mode == _3dMat ) {
        giveRealStressVectorComputedFromStrain(answer, form, gp, totalStrain, atTime);
    } else if ( mode == _3dMat_F ) {
        giveRealStressVectorComputedFromDefGrad(answer, form, gp, totalStrain, atTime);
    } else {
        OOFEM_ERROR("MisesMat::giveRealStressVector : unknown material response mode");
    }

    return;
}

// returns the stress vector in 3d stress space
// computed from the previous plastic strain and current deformation gradient
void
MisesMat :: giveRealStressVectorComputedFromDefGrad(FloatArray &answer,
                                                    MatResponseForm form,
                                                    GaussPoint *gp,
                                                    const FloatArray &totalDefGrad,
                                                    TimeStep *atTime)
{
    MisesMatStatus *status = ( MisesMatStatus * ) this->giveStatus(gp);

    this->initTempStatus(gp);
    this->initGpForNewStep(gp);

    double kappa, dKappa, yieldValue, mi;
    FloatMatrix F, oldF, invOldF;
    FloatArray totalStrain;
    FloatArray s(6);

    //converts deformation gradient from array into 3x3 matrix
    F.resize(3, 3);
    F.at(1, 1) = totalDefGrad(0);
    F.at(2, 1) = totalDefGrad(1);
    F.at(3, 1) = totalDefGrad(2);
    F.at(1, 2) = totalDefGrad(3);
    F.at(2, 2) = totalDefGrad(4);
    F.at(3, 2) = totalDefGrad(5);
    F.at(1, 3) = totalDefGrad(6);
    F.at(2, 3) = totalDefGrad(7);
    F.at(3, 3) = totalDefGrad(8);

    kappa = status->giveCumulativePlasticStrain();
    status->giveTempDefGrad(oldF);
    invOldF.beInverseOf(oldF);
    //relative deformation radient
    FloatMatrix f, f_T;
    f.resize(3, 3);
    f.beProductOf(F, invOldF);

    //compute elastic predictor
    FloatMatrix oldLeftCauchyGreen, trialLeftCauchyGreen, help;
    oldLeftCauchyGreen.resize(3, 3);
    trialLeftCauchyGreen.resize(3, 3);
    help.resize(3, 3);

    f.times( pow(f.giveDeterminant(), -1. / 3.) );
    f_T.beTranspositionOf(f);
    status->giveTempLeftCauchyGreen(oldLeftCauchyGreen);
    help.beProductOf(f, oldLeftCauchyGreen);
    trialLeftCauchyGreen.beProductOf(help, f_T);
    FloatMatrix def, F_T;
    F_T.beTranspositionOf(F);
    def.beProductOf(F, F_T);
    def.at(1, 1) = def.at(1, 1) - 1;
    def.at(2, 2) = def.at(2, 2) - 1;
    def.at(3, 3) = def.at(3, 3) - 1;
    def.times(1. / 2.);
    FloatArray vDef(6);
    vDef.at(1) = def.at(1, 1);
    vDef.at(2) = def.at(2, 2);
    vDef.at(3) = def.at(3, 3);
    vDef.at(4) = 2. * def.at(2, 3);
    vDef.at(5) = 2. * def.at(1, 3);
    vDef.at(6) = 2. * def.at(1, 2);

    StrainVector leftCauchyGreen(_3dMat);
    StrainVector leftCauchyGreenDev(_3dMat);
    double leftCauchyGreenVol;

    leftCauchyGreen(0) = trialLeftCauchyGreen.at(1, 1);
    leftCauchyGreen(1) = trialLeftCauchyGreen.at(2, 2);
    leftCauchyGreen(2) = trialLeftCauchyGreen.at(3, 3);
    leftCauchyGreen(3) = 2. * trialLeftCauchyGreen.at(2, 3);
    leftCauchyGreen(4) = 2. * trialLeftCauchyGreen.at(1, 3);
    leftCauchyGreen(5) = 2. * trialLeftCauchyGreen.at(1, 2);


    leftCauchyGreen.computeDeviatoricVolumetricSplit(leftCauchyGreenDev, leftCauchyGreenVol);
    StressVector trialStressDev(_3dMat);
    leftCauchyGreenDev.applyDeviatoricElasticStiffness(trialStressDev, G / 2.);
    s = trialStressDev;

    //check for plastic loading
    double trialS = trialStressDev.computeStressNorm();
    double sigmaY = sig0 + H * kappa;
    yieldValue = trialS - sqrt(2. / 3.) * sigmaY;


    //store deviatoric trial stress(reused by algorithmic stiffness)
    status->letTrialStressDevBe(trialStressDev);
    //the return-mapping algorithm
    double J = F.giveDeterminant();
    mi = leftCauchyGreenVol * G;
    if ( yieldValue > 0 ) {
        dKappa = ( yieldValue / ( 2 * mi ) ) / ( 1 + H / ( 3 * mi ) );
        FloatArray n = trialStressDev;
        n.times(2 * mi * dKappa / trialS);
        ////return map
        s = trialStressDev - n;
        kappa += sqrt(2. / 3.) * dKappa;


        //update of intermediate configuration
        trialLeftCauchyGreen.at(1, 1) = s(0) / G + leftCauchyGreenVol;
        trialLeftCauchyGreen.at(2, 2) = s(1) / G + leftCauchyGreenVol;
        trialLeftCauchyGreen.at(3, 3) = s(2) / G + leftCauchyGreenVol;
        trialLeftCauchyGreen.at(2, 3) = s(3) / G;
        trialLeftCauchyGreen.at(1, 3) = s(4) / G;
        trialLeftCauchyGreen.at(1, 2) = s(5) / G;
        trialLeftCauchyGreen.at(3, 2) = trialLeftCauchyGreen.at(2, 3);
        trialLeftCauchyGreen.at(3, 1) = trialLeftCauchyGreen.at(1, 3);
        trialLeftCauchyGreen.at(2, 1) = trialLeftCauchyGreen.at(1, 2);
        trialLeftCauchyGreen.times(J * J);
    }

    //addition of the elastic mean stress


    FloatMatrix kirchhoffStress;
    kirchhoffStress.resize(3, 3);
    kirchhoffStress.at(1, 1) = s(0) + 1. / 2. *  K * ( J * J - 1 );
    kirchhoffStress.at(2, 2) = s(1) + 1. / 2. *  K * ( J * J - 1 );
    kirchhoffStress.at(3, 3) = s(2) + 1. / 2. *  K * ( J * J - 1 );
    kirchhoffStress.at(2, 3) = s(3);
    kirchhoffStress.at(1, 3) = s(4);
    kirchhoffStress.at(1, 2) = s(5);
    kirchhoffStress.at(3, 2) = kirchhoffStress.at(2, 3);
    kirchhoffStress.at(3, 1) = kirchhoffStress.at(1, 3);
    kirchhoffStress.at(2, 1) = kirchhoffStress.at(1, 2);

    //transform Kirchhoff stress into Second Piola - Kirchhoff stress
    FloatMatrix iF(3, 3), iF_T(3, 3);
    iF.beInverseOf(F);
    iF_T.beTranspositionOf(iF);
    FloatMatrix S;
    help.resize(3, 3);
    S.resize(3, 3);

    help.beProductOf(iF, kirchhoffStress);
    S.beProductOf(help, iF_T);

    FloatMatrix Ep(3, 3);
    FloatMatrix E(3, 3);
    FloatArray ep(6);
    FloatArray e(6);
    this->computeGLPlasticStrain(F, Ep, trialLeftCauchyGreen, J);
    this->convertDefGradToGLStrain(F, E);

    ep(0) = Ep.at(1, 1);
    ep(1) = Ep.at(2, 2);
    ep(2) = Ep.at(3, 3);
    ep(3) = 2 * Ep.at(2, 3);
    ep(4) = 2 * Ep.at(1, 3);
    ep(5) = 2 * Ep.at(1, 2);

    e(0) = E.at(1, 1);
    e(1) = E.at(2, 2);
    e(2) = E.at(3, 3);
    e(3) = 2 * E.at(2, 3);
    e(4) = 2 * E.at(1, 3);
    e(5) = 2 * E.at(1, 2);

    answer.resize(6);
    answer(0) =  S.at(1, 1);
    answer(1) =  S.at(2, 2);
    answer(2) =  S.at(3, 3);
    answer(3) =  S.at(2, 3);
    answer(4) =  S.at(1, 3);
    answer(5) =  S.at(1, 2);


    status->setTrialStressVol(mi);
    status->letTempDefGradBe(F);
    status->letTempLeftCauchyGreenBe(trialLeftCauchyGreen);
    status->setTempCumulativePlasticStrain(kappa);
    status->letTempStressVectorBe(answer);
    status->letTempStrainVectorBe(e);
    status->letTempPlasticStrainBe(ep);
    return;
}

// converts the deformation gradient stored by columns in FloatArray F
// into the Green-Lagrange strain with components E11,E22,E33,E23,E31,E12 stored in FloatArray E
void
MisesMat :: convertDefGradToGLStrain(const FloatMatrix &F, FloatMatrix &E)
{
    E.resize(3, 3);
    FloatMatrix F_T;
    F_T.beTranspositionOf(F);
    E.beProductOf(F_T, F);
    E.at(1, 1) = E.at(1, 1) - 1;
    E.at(2, 2) = E.at(2, 2) - 1;
    E.at(3, 3) = E.at(3, 3) - 1;
    E.times(1. / 2.);
    return;
}
void
MisesMat :: computeGLPlasticStrain(const FloatMatrix &F, FloatMatrix &Ep, FloatMatrix b, double J)
{
    FloatMatrix I, F_T, invB, help;
    I.resize(3, 3);
    help.resize(3, 3);
    I.beUnitMatrix();
    I.times(-1);
    F_T.beTranspositionOf(F);
    invB.beInverseOf(b);
    help.beProductOf(F_T, invB);
    Ep.beProductOf(help, F);
    Ep.times( pow(J, -2. / 3.) );
    Ep.plus(I);
    Ep.times(1. / 2.);

    return;
}
// returns the stress vector in 3d stress space
// computed from the previous plastic strain and current total strain
void
MisesMat :: giveRealStressVectorComputedFromStrain(FloatArray &answer,
                                                   MatResponseForm form,
                                                   GaussPoint *gp,
                                                   const FloatArray &totalStrain,
                                                   TimeStep *atTime)
{
    double kappa, yieldValue, dKappa;
    FloatArray reducedStress;
    FloatArray strain, plStrain;
    StructuralCrossSection *crossSection = ( StructuralCrossSection * )
                                           ( gp->giveElement()->giveCrossSection() );
    MisesMatStatus *status = ( MisesMatStatus * ) this->giveStatus(gp);

    this->initTempStatus(gp);
    this->initGpForNewStep(gp);

    // subtract the stress-independent part of strain (e.g. thermal strain)
    this->giveStressDependentPartOfStrainVector(strain, gp, totalStrain,
                                                atTime, VM_Total);

    // get the initial plastic strain and initial kappa from the status
    status->givePlasticStrain(plStrain);
    kappa = status->giveCumulativePlasticStrain();
    // === radial return algorithm ===

    // elastic predictor
    StrainVector elStrain(totalStrain, _3dMat);
    elStrain.substract(plStrain);
    StrainVector elStrainDev(_3dMat);
    double elStrainVol;
    elStrain.computeDeviatoricVolumetricSplit(elStrainDev, elStrainVol);
    StressVector trialStressDev(_3dMat);
    elStrainDev.applyDeviatoricElasticStiffness(trialStressDev, G);
    /**************************************************************/
    double trialStressVol;
    trialStressVol = 3 * K * elStrainVol;
    /**************************************************************/

    // store the deviatoric trial stress (reused by algorithmic stiffness)
    status->letTrialStressDevBe(trialStressDev);
    status->setTrialStressVol(trialStressVol);


    // check the yield condition at the trial state
    double trialS = trialStressDev.computeStressNorm();
    double sigmaY = sig0 + H * kappa;
    yieldValue = sqrt(3. / 2.) * trialS - sigmaY;
    if ( yieldValue > 0. ) {
        // increment of cumulative plastic strain
        dKappa = yieldValue / ( H + 3. * G );
        kappa += dKappa;
        StrainVector dPlStrain(_3dMat);
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
    StressVector fullStress(_3dMat);
    double stressVol = 3. * K * elStrainVol;
    trialStressDev.computeDeviatoricVolumetricSum(fullStress, stressVol);

    // store the total strain in status
    status->letTempStrainVectorBe(totalStrain);

    // reduce the stress vector and store it in status
    crossSection->giveReducedCharacteristicVector(reducedStress, gp, fullStress);
    status->letTempStressVectorBe(reducedStress);

    // store the plastic strain and cumulative plastic strain
    status->letTempPlasticStrainBe(plStrain);
    status->setTempCumulativePlasticStrain(kappa);
    // reduce the stress vector (if required) and return it
    if ( form == FullForm ) {
        answer = fullStress;
    } else {
        answer = reducedStress;
    }

    return;
}





void
MisesMat :: give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    MaterialMode mMode = gp->giveMaterialMode();

    if ( mMode == _3dMat ) {
        give3dSSMaterialStiffnessMatrix(answer, form, mode, gp, atTime);
    } else if ( mMode == _3dMat_F ) {
        give3dLSMaterialStiffnessMatrix(answer, form, mode, gp, atTime);
    } else {
        OOFEM_ERROR("MisesMat::give3dMaterialStiffnessMatrix : unknown material response mode");
    }

    return;
}


// returns the consistent (algorithmic) tangent stiffness matrix
void
MisesMat :: give3dSSMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form,
                                            MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *atTime)
{
    // start from the elastic stiffness
    this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
    if ( mode != TangentStiffness ) {
        return;
    }

    MisesMatStatus *status = ( MisesMatStatus * ) this->giveStatus(gp);
    double kappa = status->giveCumulativePlasticStrain();
    // increment of cumulative plastic strain as an indicator of plastic loading
    double dKappa = status->giveTempCumulativePlasticStrain() - kappa;

    if ( dKappa <= 0.0 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
        return;
    }

    // === plastic loading ===

    // yield stress at the beginning of the step
    double sigmaY = sig0 + H * kappa;

    // trial deviatoric stress and its norm
    StressVector trialStressDev(_3dMat);

    status->giveTrialStressDev(trialStressDev);
    double trialS = trialStressDev.computeStressNorm();

    // one correction term
    FloatMatrix stiffnessCorrection(6, 6);
    stiffnessCorrection.beDyadicProductOf(trialStressDev, trialStressDev);
    double factor = -2. * sqrt(6.) * G * G / trialS;
    double factor1 = factor * sigmaY / ( ( H + 3. * G ) * trialS * trialS );
    stiffnessCorrection.times(factor1);
    answer.plus(stiffnessCorrection);

    // another correction term
    stiffnessCorrection.bePinvID();
    double factor2 = factor * dKappa;
    stiffnessCorrection.times(factor2);
    answer.plus(stiffnessCorrection);

    return;
}



void
MisesMat :: give3dLSMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    MisesMatStatus *status = ( MisesMatStatus * ) this->giveStatus(gp);
    // start from the elastic stiffness

    FloatMatrix I(6, 6);
    I.at(1, 1) = I.at(2, 2) = I.at(3, 3) = 1;
    I.at(4, 4) = I.at(5, 5) = I.at(6, 6) = 0.5;
    FloatArray delta(6);
    delta.at(1) = delta.at(2) = delta.at(3) = 1;

    FloatMatrix F, F_Tr;
    status->giveTempDefGrad(F);
    double J;
    J = F.giveDeterminant();

    StressVector trialStressDev(_3dMat);
    double trialStressVol;
    status->giveTrialStressVol(trialStressVol);
    status->giveTrialStressDev(trialStressDev);
    double trialS = trialStressDev.computeStressNorm();
    FloatArray n(6);
    n = trialStressDev;
    if ( trialS == 0 ) {
        n.resize(6);
    } else                 {
        n.times(1 / trialS);
    }


    FloatMatrix Cdev(6, 6);
    FloatMatrix C(6, 6);
    FloatMatrix help(6, 6);
    help.beDyadicProductOf(delta, delta);
    C = help;
    help.times(-1. / 3.);
    FloatMatrix C1 = I;
    C1.plus(help);
    C1.times(2 * trialStressVol);

    FloatMatrix n1(6, 6), n2(6, 6);
    n1.beDyadicProductOf(n, delta);
    n2.beDyadicProductOf(delta, n);
    help = n1;
    help.plus(n2);
    help.times(-2. / 3. * trialS);
    C1.plus(help);
    Cdev = C1;
    C.times(K * J * J);

    help = I;
    help.times( -K * ( J * J - 1 ) );
    C.plus(help);
    FloatMatrix Cvol = C;
    C.plus(C1);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    FloatMatrix invF(3, 3);
    FloatMatrix T(6, 6), tT(6, 6);

    invF.beInverseOf(F);
    //////////////////////////////////////////////////
    //first row of pull back transformation matrix
    T.at(1, 1) = invF.at(1, 1) * invF.at(1, 1);
    T.at(1, 2) = invF.at(1, 2) * invF.at(1, 2);
    T.at(1, 3) = invF.at(1, 3) * invF.at(1, 3);
    T.at(1, 4) = 2. * invF.at(1, 2) * invF.at(1, 3);
    T.at(1, 5) = 2. * invF.at(1, 1) * invF.at(1, 3);
    T.at(1, 6) = 2. * invF.at(1, 1) * invF.at(1, 2);
    //second row of pull back transformation matrix
    T.at(2, 1) = invF.at(2, 1) * invF.at(2, 1);
    T.at(2, 2) = invF.at(2, 2) * invF.at(2, 2);
    T.at(2, 3) = invF.at(2, 3) * invF.at(2, 3);
    T.at(2, 4) = 2. * invF.at(2, 2) * invF.at(2, 3);
    T.at(2, 5) = 2. * invF.at(2, 1) * invF.at(2, 3);
    T.at(2, 6) = 2. * invF.at(2, 1) * invF.at(2, 2);
    //third row of pull back transformation matrix
    T.at(3, 1) = invF.at(3, 1) * invF.at(3, 1);
    T.at(3, 2) = invF.at(3, 2) * invF.at(3, 2);
    T.at(3, 3) = invF.at(3, 3) * invF.at(3, 3);
    T.at(3, 4) = 2. * invF.at(3, 2) * invF.at(3, 3);
    T.at(3, 5) = 2. * invF.at(3, 1) * invF.at(3, 3);
    T.at(3, 6) = 2. * invF.at(3, 1) * invF.at(3, 2);
    //fourth row of pull back transformation matrix
    T.at(4, 1) = invF.at(2, 1) * invF.at(3, 1);
    T.at(4, 2) = invF.at(2, 2) * invF.at(3, 2);
    T.at(4, 3) = invF.at(2, 3) * invF.at(3, 3);
    T.at(4, 4) = ( invF.at(2, 2) * invF.at(3, 3) + invF.at(2, 3) * invF.at(3, 2) );
    T.at(4, 5) = ( invF.at(2, 1) * invF.at(3, 3) + invF.at(2, 3) * invF.at(3, 1) );
    T.at(4, 6) = ( invF.at(2, 1) * invF.at(3, 2) + invF.at(2, 2) * invF.at(3, 1) );
    //fifth row of pull back transformation matrix
    T.at(5, 1) = invF.at(1, 1) * invF.at(3, 1);
    T.at(5, 2) = invF.at(1, 2) * invF.at(3, 2);
    T.at(5, 3) = invF.at(1, 3) * invF.at(3, 3);
    T.at(5, 4) = ( invF.at(1, 2) * invF.at(3, 3) + invF.at(1, 3) * invF.at(3, 2) );
    T.at(5, 5) = ( invF.at(1, 1) * invF.at(3, 3) + invF.at(1, 3) * invF.at(3, 1) );
    T.at(5, 6) = ( invF.at(1, 1) * invF.at(3, 2) + invF.at(1, 2) * invF.at(3, 1) );
    //sixth row of pull back transformation matrix
    T.at(6, 1) = invF.at(1, 1) * invF.at(2, 1);
    T.at(6, 2) = invF.at(1, 2) * invF.at(2, 2);
    T.at(6, 3) = invF.at(1, 3) * invF.at(2, 3);
    T.at(6, 4) = ( invF.at(1, 2) * invF.at(2, 3) + invF.at(1, 3) * invF.at(2, 2) );
    T.at(6, 5) = ( invF.at(1, 1) * invF.at(2, 3) + invF.at(1, 3) * invF.at(2, 1) );
    T.at(6, 6) = ( invF.at(1, 1) * invF.at(2, 2) + invF.at(1, 2) * invF.at(2, 1) );
    ///////////////////////////////////////////
    tT.beTranspositionOf(T);
    answer.resize(6, 6);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    if ( mode != TangentStiffness ) {
        help.beProductOf(C, tT);
        answer.beProductOf(T, help);
        return;
    }

    double kappa = status->giveCumulativePlasticStrain();
    // increment of cumulative plastic strain as an indicator of plastic loading
    double dKappa = sqrt(3. / 2.) * ( status->giveTempCumulativePlasticStrain() - kappa );
    if ( dKappa <= 0.0 ) { // elastic loading - elastic stiffness plays the role of tangent stiffness
        help.beProductOf(C, tT);
        answer.beProductOf(T, help);
        return;
    }

    // === plastic loading ===
    // trial deviatoric stress and its norm


    double beta0, beta1, beta2, beta3, beta4;
    if ( trialS == 0 ) {
        beta1 = 0;
    } else {
        beta1 = 2 * trialStressVol * dKappa / trialS;
    }

    if ( trialStressVol == 0 ) {
        beta0 = 0;
        beta2 = 0;
        beta3 = beta1;
        beta4 = 0;
    } else   {
        beta0 = 1 + H / 3 / trialStressVol;
        beta2 = ( 1 - 1 / beta0 ) * 2. / 3. * trialS * dKappa / trialStressVol;
        beta3 = 1 / beta0 - beta1 + beta2;
        beta4 = ( 1 / beta0 - beta1 ) * trialS / trialStressVol;
    }

    FloatMatrix N;
    N.beDyadicProductOf(n, n);
    N.times(-2 * trialStressVol * beta3);
    answer.resize(6, 6);

    C1.times(-beta1);
    FloatMatrix mN(3, 3);
    mN.at(1, 1) = n.at(1);
    mN.at(1, 2) = n.at(6);
    mN.at(1, 3) = n.at(5);
    mN.at(2, 1) = n.at(6);
    mN.at(2, 2) = n.at(2);
    mN.at(2, 3) = n.at(4);
    mN.at(3, 1) = n.at(5);
    mN.at(3, 2) = n.at(4);
    mN.at(3, 3) = n.at(3);
    FloatMatrix mN2(3, 3);
    mN2.beProductOf(mN, mN);

    double volN2 = 1. / 3. * ( mN2.at(1, 1) + mN2.at(2, 2) + mN2.at(3, 3) );
    FloatArray devN2(6);
    devN2.at(1) = mN2.at(1, 1) - volN2;
    devN2.at(2) = mN2.at(2, 2) - volN2;

    devN2.at(3) = mN2.at(3, 3) - volN2;
    devN2.at(4) = mN2.at(2, 3);
    devN2.at(5) = mN2.at(1, 3);
    devN2.at(6) = mN2.at(1, 2);
    FloatMatrix nonSymPart, symP;
    nonSymPart.beDyadicProductOf(n, devN2);
    nonSymPart.times(-2 * trialStressVol * beta4);

    C.plus(C1);
    C.plus(N);
    C.plus(nonSymPart);
    help.beProductOf(C, tT);
    answer.beProductOf(T, help);
}


#ifdef __OOFEG
#endif


const FloatArray *
MisesMatStatus :: givePlasDef()
{
    return & plasticStrain;
}

int
MisesMat :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    MisesMatStatus *status = ( MisesMatStatus * ) this->giveStatus(aGaussPoint);
    if ( type == IST_PlasticStrainTensor ) {
        answer.resize(6);
        answer  = * status->givePlasDef();
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel )      {
        answer.resize(1);
        answer.at(1) = status->giveCumulativePlasticStrain();
        return 1;
    } else       {
        return StructuralMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}




InternalStateValueType
MisesMat :: giveIPValueType(InternalStateType type)
{
    if ( type == IST_PlasticStrainTensor ) {
        return ISVT_TENSOR_S3;
    } else if ( type == IST_MaxEquivalentStrainLevel )   {
        return ISVT_SCALAR;
    } else                                                                        {
        return StructuralMaterial :: giveIPValueType(type);
    }
}


int
MisesMat :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode)
{
    if ( type == IST_PlasticStrainTensor ) {
        answer.resize(6);
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 3;
        answer.at(4) = 4;
        answer.at(5) = 5;
        answer.at(6) = 6;
        return 1;
    } else if ( type == IST_MaxEquivalentStrainLevel )     {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else   {
        return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);
    }
}


int
MisesMat :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    if ( type == IST_PlasticStrainTensor ) {
        return 6;
    } else if ( type == IST_MaxEquivalentStrainLevel )  {
        return 1;
    } else                                                               {
        return StructuralMaterial :: giveIPValueSize(type, aGaussPoint);
    }
}

//=============================================================================

MisesMatStatus :: MisesMatStatus(int n, Domain *d, GaussPoint *g) :
    StructuralMaterialStatus(n, d, g), plasticStrain(), tempPlasticStrain(), trialStressD()
{
    kappa = tempKappa = 0.;
    defGrad.resize(3, 3);
    defGrad.at(1, 1) = defGrad.at(2, 2) = defGrad.at(3, 3) = 1;
    tempDefGrad.resize(3, 3);
    tempDefGrad.at(1, 1) = tempDefGrad.at(2, 2) = tempDefGrad.at(3, 3) = 1;
    leftCauchyGreen.resize(3, 3);
    leftCauchyGreen.at(1, 1) = leftCauchyGreen.at(2, 2) = leftCauchyGreen.at(3, 3) = 1;
}

MisesMatStatus :: ~MisesMatStatus()
{ }

void
MisesMatStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    int i, n;

    StructuralMaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, "status { ");
    // print the plastic strain
    n = plasticStrain.giveSize();
    fprintf(file, " plastic strains ");
    for ( i = 1; i <= n; i++ ) {
        fprintf( file, " % .4e", plasticStrain.at(i) );
    }

    // print the cumulative plastic strain
    fprintf(file, ", kappa ");
    fprintf(file, " % .4e", kappa);

    fprintf(file, "}\n");
    //print Left Cauchy - Green deformation tensor
    fprintf(file, " left Cauchy Green");
    fprintf( file, " % .4e", tempLeftCauchyGreen.at(1, 1) );
    fprintf( file, " % .4e", tempLeftCauchyGreen.at(2, 2) );
    fprintf( file, " % .4e", tempLeftCauchyGreen.at(3, 3) );
    fprintf( file, " % .4e", tempLeftCauchyGreen.at(2, 3) );
    fprintf( file, " % .4e", tempLeftCauchyGreen.at(1, 3) );
    fprintf( file, " % .4e", tempLeftCauchyGreen.at(1, 2) );

    //print deformation gradient
    fprintf(file, " Deformation Gradient");
    fprintf( file, " % .4e", tempDefGrad.at(1, 1) );
    fprintf( file, " % .4e", tempDefGrad.at(1, 2) );
    fprintf( file, " % .4e", tempDefGrad.at(1, 3) );
    fprintf( file, " % .4e", tempDefGrad.at(2, 1) );
    fprintf( file, " % .4e", tempDefGrad.at(2, 2) );
    fprintf( file, " % .4e", tempDefGrad.at(2, 3) );
    fprintf( file, " % .4e", tempDefGrad.at(3, 1) );
    fprintf( file, " % .4e", tempDefGrad.at(3, 2) );
    fprintf( file, " % .4e", tempDefGrad.at(3, 3) );


    fprintf(file, "}\n");
}


// initializes temporary variables based on their values at the previous equlibrium state
void MisesMatStatus :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();

    if ( plasticStrain.giveSize() == 0 ) {
        plasticStrain.resize( ( ( StructuralMaterial * ) gp->giveMaterial() )->
                             giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() ) );
        plasticStrain.zero();
    }

    tempPlasticStrain = plasticStrain;
    tempKappa = kappa;
    tempDefGrad = defGrad;
    tempLeftCauchyGreen = leftCauchyGreen;
    trialStressD.resize(0); // to indicate that it is not defined yet
}


// updates internal variables when equilibrium is reached
void
MisesMatStatus :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);

    plasticStrain = tempPlasticStrain;
    kappa = tempKappa;
    trialStressD.resize(0); // to indicate that it is not defined any more
    defGrad = tempDefGrad;
    leftCauchyGreen = tempLeftCauchyGreen;
}


// saves full information stored in this status
// temporary variables are NOT stored
contextIOResultType
MisesMatStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
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

    return CIO_OK;
}



contextIOResultType
MisesMatStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
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

    return CIO_OK; // return succes
}
} // end namespace oofem
