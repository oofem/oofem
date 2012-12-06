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

#include "concretedpm2.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "structuralms.h"
#include "gausspnt.h"
#include "intarray.h"
#include "datastream.h"
#include "contextioerr.h"
#include "timestep.h"
#include "structuralmaterial.h"
#include "isolinearelasticmaterial.h"
#include "structuralcrosssection.h"
#include "mathfem.h"

namespace oofem {
ConcreteDPM2Status :: ConcreteDPM2Status(int n, Domain *d, GaussPoint *gp) :
    StructuralMaterialStatus(n, d, gp),
    plasticStrain( gp->giveMaterialMode() ),
    tempPlasticStrain( gp->giveMaterialMode() )
{
    kappaP = tempKappaP = 0.;

    equivStrain = tempEquivStrain = 0.;
    equivStrainTension = tempEquivStrainTension = 0.;
    equivStrainCompression = tempEquivStrainCompression = 0;

    kappaDTension = tempKappaDTension = 0.;
    kappaDCompression = tempKappaDCompression = 0.;

    kappaDTensionOne = tempKappaDTensionOne = 0.;
    kappaDCompressionOne = tempKappaDCompressionOne = 0.;
    kappaDTensionTwo = tempKappaDTensionTwo = 0.;
    kappaDCompressionTwo = tempKappaDCompressionTwo = 0.;

    alpha = tempAlpha = 0.;

    damage = tempDamage = 0.;
    damageTension = tempDamageTension = 0.;
    damageCompression = tempDamageCompression = 0.;

    deltaLambda = 0.;
    state_flag = temp_state_flag = ConcreteDPM2Status :: ConcreteDPM2_Elastic;
    rateFactor = 1.;
    rateStrain = tempRateStrain = 0.;
}

ConcreteDPM2Status :: ~ConcreteDPM2Status()
{}

void
ConcreteDPM2Status :: initTempStatus()
{
    // Call the function of the parent class to initialize the variables defined there.
    StructuralMaterialStatus :: initTempStatus();

    //Plasticity part
    tempPlasticStrain = plasticStrain;
    tempKappaP = kappaP;

    tempAlpha = alpha;

    //Damage part
    tempEquivStrain = equivStrain;
    tempEquivStrainTension = equivStrainTension;
    tempEquivStrainCompression = equivStrainCompression;

    tempKappaDTension = kappaDTension;
    tempKappaDCompression = kappaDCompression;

    tempKappaDTensionOne = kappaDTensionOne;
    tempKappaDCompressionOne = kappaDCompressionOne;
    tempKappaDTensionTwo = kappaDTensionTwo;
    tempKappaDCompressionTwo = kappaDCompressionTwo;

    tempDamage = damage;
    tempDamageTension = damageTension;
    tempDamageCompression = damageCompression;

    temp_state_flag = state_flag;
    tempRateFactor = rateFactor;
    tempRateStrain = rateStrain;
}

void
ConcreteDPM2Status :: updateYourself(TimeStep *atTime)
{
    // Call corresponding function of the parent class to update
    // variables defined there.
    StructuralMaterialStatus :: updateYourself(atTime);

    // update variables defined in ConcreteDPM2Status

    alpha = tempAlpha;

    //Plasticity part
    plasticStrain = tempPlasticStrain;
    kappaP = tempKappaP;

    //Damage part
    equivStrain = tempEquivStrain;
    equivStrainTension = tempEquivStrainTension;
    equivStrainCompression = tempEquivStrainCompression;

    kappaDTension = tempKappaDTension;
    kappaDCompression = tempKappaDCompression;

    kappaDTensionOne = tempKappaDTensionOne;
    kappaDCompressionOne = tempKappaDCompressionOne;

    kappaDTensionTwo = tempKappaDTensionTwo;
    kappaDCompressionTwo = tempKappaDCompressionTwo;

    damage = tempDamage;
    damageTension = tempDamageTension;
    damageCompression = tempDamageCompression;

    state_flag = temp_state_flag;

    rateFactor = tempRateFactor;

    rateStrain = tempRateStrain;
}

void
ConcreteDPM2Status :: printOutputAt(FILE *file, TimeStep *tStep)
{
    // Call corresponding function of the parent class to print
    // variables defined there.
    StructuralMaterialStatus :: printOutputAt(file, tStep);

    fprintf(file, "\tstatus { ");

    // print status flag
    switch ( state_flag ) {
    case ConcreteDPM2Status :: ConcreteDPM2_Elastic:
        fprintf(file, "Elastic, ");
        break;
    case ConcreteDPM2Status :: ConcreteDPM2_Unloading:
        fprintf(file, "Unloading, ");
        break;
    case ConcreteDPM2Status :: ConcreteDPM2_Plastic:
        fprintf(file, "Plastic, ");
        break;
    case ConcreteDPM2Status :: ConcreteDPM2_Damage:
        fprintf(file, "Damage, ");
        break;
    case ConcreteDPM2Status :: ConcreteDPM2_PlasticDamage:
        fprintf(file, "PlasticDamage, ");
        break;
    }

    // print plastic strain vector and inelastic strain vector
    StrainVector plasticStrainVector( gp->giveMaterialMode() );
    giveFullPlasticStrainVector(plasticStrainVector);

    //    StructuralCrossSection *crossSection = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();

    StrainVector inelasticStrainVector( gp->giveMaterialMode() );
    ( ( StructuralCrossSection * )
      gp->giveCrossSection() )->giveFullCharacteristicVector(inelasticStrainVector, gp, strainVector);

    inelasticStrainVector.subtract(plasticStrainVector);
    inelasticStrainVector.times(damage);
    inelasticStrainVector.add(plasticStrainVector);
    inelasticStrainVector.times(le);

    fprintf(file, "plastic ");
    int n = plasticStrainVector.giveSize();
    for ( int i = 1; i <= n; i++ ) {
        fprintf( file, " % .4e", plasticStrainVector.at(i) );
    }

    fprintf(file, " inelastic ");
    n = inelasticStrainVector.giveSize();
    for ( int i = 1; i <= n; i++ ) {
        fprintf( file, " % .4e", inelasticStrainVector.at(i) );
    }

    fprintf(file, "}\n");


    fprintf(file, "\t\tkappaDTension ");
    // print hardening parameter
    fprintf(file, " % .4e\n", kappaDTension);

    fprintf(file, "\t\tkappaDCompression ");
    // print hardening parameter
    fprintf(file, " % .4e\n", kappaDCompression);

    fprintf(file, "\t\tkappaP ");
    // print hardening parameter
    fprintf(file, " % .4e\n", kappaP);

    fprintf(file, "\t\tkappaDTensionOne ");
    // print hardening parameter
    fprintf(file, " % .4e\n", kappaDTensionOne);

    fprintf(file, "\t\tkappaDCompressionOne ");
    // print hardening parameter
    fprintf(file, " % .4e\n", kappaDCompressionOne);

    fprintf(file, "\t\tkappaDTensionTwo ");
    // print hardening parameter
    fprintf(file, " % .4e\n", kappaDTensionTwo);

    fprintf(file, "\t\tkappaDCompressionTwo ");
    // print hardening parameter
    fprintf(file, " % .4e\n", kappaDCompressionTwo);

    fprintf(file, "\t\tdamage ");
    // print hardening parameter
    fprintf(file, " % .4e\n", damage);

    fprintf(file, "\t\tdamageTension ");
    // print hardening parameter
    fprintf(file, " % .4e\n", damageTension);

    fprintf(file, "\t\tdamageCompression ");
    // print hardening parameter
    fprintf(file, " % .4e\n", damageCompression);
}

contextIOResultType
ConcreteDPM2Status :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data

    if ( ( iores = plasticStrain.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream->write(& kappaP, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& equivStrain, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& equivStrainTension, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& equivStrainCompression, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& kappaDTensionOne, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& kappaDCompressionOne, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& kappaDTensionTwo, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& kappaDCompressionTwo, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& kappaDTension, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& kappaDCompression, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& damageTension, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& damageCompression, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& state_flag, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& le, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}


contextIOResultType
ConcreteDPM2Status :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( ( iores = plasticStrain.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream->read(& kappaP, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& equivStrain, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& equivStrainTension, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& equivStrainCompression, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& kappaDTensionOne, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& kappaDCompressionOne, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& kappaDTensionTwo, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& kappaDCompressionTwo, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& kappaDTension, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& kappaDCompression, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& damageTension, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& damageCompression, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& state_flag, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& le, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}



//   ********************************
//   *** CLASS DYNAMIC CONCRETE   ***
//   ********************************

#define IDM_ITERATION_LIMIT 1.e-8

ConcreteDPM2 :: ConcreteDPM2(int n, Domain *d) :
    StructuralMaterial(n, d)
{
    yieldTol = 0.;
    newtonIter = 0;
    matMode = _Unknown;

    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
}

ConcreteDPM2 :: ~ConcreteDPM2()
{
    delete linearElasticMaterial;
}

IRResultType
ConcreteDPM2 :: initializeFrom(InputRecord *ir)
{
    // Required by IR_GIVE_FIELD macro
      const char *__proc = "initializeFrom";
    IRResultType result;

    // call the corresponding service for the linear elastic material
    StructuralMaterial :: initializeFrom(ir);
    linearElasticMaterial->initializeFrom(ir);

    //isotropic flag
    isotropicFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, isotropicFlag, IFT_ConcreteDPM2_isoflag, "isoflag");

    double value;
    // elastic parameters
    IR_GIVE_FIELD(ir, eM, IFT_IsotropicLinearElasticMaterial_e, "e"); // Macro
    IR_GIVE_FIELD(ir, nu, IFT_IsotropicLinearElasticMaterial_n, "n"); // Macro
    propertyDictionary->add('E', eM);
    propertyDictionary->add('n', nu);

    IR_GIVE_FIELD(ir, value, IFT_IsotropicLinearElasticMaterial_talpha, "talpha"); // Macro
    propertyDictionary->add(tAlpha, value);

    gM = eM / ( 2. * ( 1. + nu ) );
    kM = eM / ( 3. * ( 1. - 2. * nu ) );

    IR_GIVE_FIELD(ir, fc, IFT_ConcreteDPM2_fc, "fc");
    IR_GIVE_FIELD(ir, ft, IFT_ConcreteDPM2_ft, "ft");

    this->e0 = this->ft / this->eM;

    // default parameters
    this->ecc = 0.525;
    IR_GIVE_OPTIONAL_FIELD(ir, ecc, IFT_ConcreteDPM2_ecc, "ecc");
    yieldHardInitial = 0.3;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldHardInitial, IFT_ConcreteDPM2_kinit, "kinit");
    //Inclination at transition point
    yieldHardPrimePeak = 0.5;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldHardPrimePeak, IFT_ConcreteDPM2_hp, "hp");

    if ( yieldHardPrimePeak < 0 ) {
        yieldHardPrimePeak = 0.;
        _warning("kPrimePeak cannot be less than zero\n");
    } else if ( yieldHardPrimePeak > ( 1. - yieldHardInitial ) )           {
        yieldHardPrimePeak = 1. - yieldHardInitial;
        _warning("kPrimePeak cannot be greater than 1.-kinit\n");
    }

    AHard = 8.e-2;
    IR_GIVE_OPTIONAL_FIELD(ir, AHard, IFT_ConcreteDPM2_ahard, "ahard");
    BHard = 3.e-3;
    IR_GIVE_OPTIONAL_FIELD(ir, BHard, IFT_ConcreteDPM2_bhard, "bhard");
    CHard = 2.;
    IR_GIVE_OPTIONAL_FIELD(ir, CHard, IFT_ConcreteDPM2_chard, "chard");
    DHard = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, DHard, IFT_ConcreteDPM2_dhard, "dhard");

    dilationConst = 0.85;
    IR_GIVE_OPTIONAL_FIELD(ir, dilationConst, IFT_ConcreteDPM2_dilation, "dilation");

    softeningType = 0; //0-Linear softening; 1-Bilinear softening
    IR_GIVE_OPTIONAL_FIELD(ir, softeningType, IFT_ConcreteDPM2_softeningType, "stype"); // Macro

    if ( softeningType > 1 ) {
        _error("softening type not implemented\n");
    }

    IR_GIVE_FIELD(ir, this->wf, IFT_ConcreteDPM2_wf, "wf");

    if ( softeningType == 1 ) {
        this->ftOne = 0.3 * this->ft;
        IR_GIVE_OPTIONAL_FIELD(ir, ftOne, IFT_ConcreteDPM2_ftOne, "ft1");
        this->wfOne = 0.15 * this->wf;
        IR_GIVE_OPTIONAL_FIELD(ir, this->wfOne, IFT_ConcreteDPM2_wfOne, "wf1");
    }

    this->ASoft = 15;
    IR_GIVE_OPTIONAL_FIELD(ir, ASoft, IFT_ConcreteDPM2_asoft, "asoft");

    this->BSoft = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, BSoft, IFT_ConcreteDPM2_bsoft, "bsoft");

    helem = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, helem, IFT_ConcreteDPM2_helem, "helem");


    //Compute m
    m = 3. * ( pow(fc, 2.) - pow(ft, 2.) ) / ( fc * ft ) * ecc / ( ecc + 1. );

    //Compute default value of dilationConst
    yieldTol = 1.e-10;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldTol, IFT_ConcreteDPM2_yieldtol, "yieldtol");
    newtonIter = 100;
    IR_GIVE_OPTIONAL_FIELD(ir, newtonIter, IFT_ConcreteDPM2_newtoniter, "newtoniter");

    //parameters for rate dependence
    strainRateFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, strainRateFlag, IFT_ConcreteDPM2_rateFlag, "rateflag");

    timeFactor = 1.;
    if ( strainRateFlag == 1 ) {
        IR_GIVE_FIELD(ir, fcZero, IFT_ConcreteDPM2_fcZero, "fczero");
        IR_GIVE_OPTIONAL_FIELD(ir, timeFactor, IFT_ConcreteDPM2_newtoniter, "timefactor");
    }

    return IRRT_OK;
}

int
ConcreteDPM2 :: hasMaterialModeCapability(MaterialMode mMode)
{
    if ( ( mMode == _3dMat ) ||
         ( mMode == _PlaneStrain ) ) {
        return 1;
    } else {
        return 0;
    }
}

void
ConcreteDPM2 :: giveRealStressVector(FloatArray &answer,
                                     MatResponseForm form,
                                     GaussPoint *gp,
                                     const FloatArray &strainVector,
                                     TimeStep *atTime)
{
    if ( matMode == _Unknown ) {
        matMode = gp->giveMaterialMode();
    }

    ConcreteDPM2Status *status = ( ConcreteDPM2Status * ) this->giveStatus(gp);

    // Initialize temp variables for this gauss point
    status->initTempStatus();

    status->letTempStrainVectorBe(strainVector);


    FloatArray testStrain;

    //  ConcreteDPM2Status *status = giveStatus (gp) ;

    StructuralCrossSection *crossSection = ( StructuralCrossSection * ) gp->giveElement()->giveCrossSection();

    StrainVector strain(strainVector, matMode);
    //Calculate the strain increment!
    StrainVector deltaStrain(matMode);

    //Calculate strain rate
    //Time step
    double deltaTime = 1.;
    if ( strainRateFlag == 1 ) {
        if ( atTime->giveTimeIncrement() == 0 ) { //Problem with the first step. For some reason the time increment is zero
            deltaTime = this->timeFactor;
        } else   {
            deltaTime = atTime->giveTimeIncrement() * this->timeFactor;
        }
    } else   {
        if ( atTime->giveTimeIncrement() == 0 ) {
            deltaTime = 1.;
        } else   {
            deltaTime = atTime->giveTimeIncrement();
        }
    }

    // perform plasticity return
    performPlasticityReturn(gp, strain);

    // compute the nominal stress
    StressVector stress(matMode);
    StressVector effectiveStress(matMode);

    // compute elastic strains and trial stress
    StrainVector elasticStrain = strain;
    StrainVector tempPlasticStrain(matMode);
    status->giveTempPlasticStrain(tempPlasticStrain);
    elasticStrain.subtract(tempPlasticStrain);
    elasticStrain.applyElasticStiffness(effectiveStress, eM, nu);


    StressVector effectiveStressTension(matMode);
    StressVector effectiveStressCompression(matMode);
    double alpha = computeAlpha(effectiveStressTension, effectiveStressCompression, effectiveStress);
    status->letTempAlphaBe(alpha);


    FloatArray damages(2);
    damages.zero();
    computeDamage(damages, strain, deltaTime, gp, atTime, alpha);

    //Split damage in a tension and compression part

    if ( isotropicFlag == 0 ) { //Default
        effectiveStressTension.times( 1. - damages.at(1) );
        effectiveStressCompression.times( 1. - damages.at(2) );
        stress = effectiveStressTension;
        stress.add(effectiveStressCompression);
    } else   { //Consider only tensile damage. Reduction to a fully isotropic model
        stress = effectiveStress;
        stress.times( 1. - damages.at(1) );
    }

    status->letTempStressVectorBe(stress);

    assignStateFlag(gp);

    if ( form == ReducedForm ) {
        answer = stress;
    } else {
        crossSection->giveFullCharacteristicVector(answer, gp, stress);
    }
}


void
ConcreteDPM2 :: computeDamage(FloatArray &answer,
                              const StrainVector &strain,
                              const double deltaTime,
                              GaussPoint *gp,
                              TimeStep *atTime,
                              const double tempAlpha)
{
    ConcreteDPM2Status *status = ( ConcreteDPM2Status * ) this->giveStatus(gp);

    double tempEquivStrain;
    double deltaPlasticStrainNorm;
    double tempDamageTension;
    double tempDamageCompression;

    double tempKappaDTension, tempKappaDCompression;
    double tempKappaDTensionOne, tempKappaDTensionTwo;
    double tempKappaDCompressionOne, tempKappaDCompressionTwo;

    double rateFactor;

    computeEquivalentStrain(tempEquivStrain, strain, gp, atTime);

    if ( ( status->giveDamageTension() == 0. ) && ( status->giveDamageCompression() == 0. ) ) {
        rateFactor = computeRateFactor(tempAlpha, deltaTime, gp, atTime);
    } else   {
        rateFactor = status->giveRateFactor();
    }

    //Compute equivalent strains for  tension and compression
    double tempEquivStrainTension = status->giveEquivStrainTension() + ( tempEquivStrain - status->giveEquivStrain() ) / rateFactor;
    double tempEquivStrainCompression = status->giveEquivStrainCompression() +
                                        ( tempAlpha * tempEquivStrain - status->giveAlpha() * status->giveEquivStrain() ) / rateFactor;

    double fTension = tempEquivStrainTension - status->giveKappaDTension();
    double fCompression = tempEquivStrainCompression - status->giveKappaDCompression();

    //If damage threshold is exceeded determine the rate factor from the previous step
    if ( ( tempEquivStrainTension > e0 || tempEquivStrainCompression > e0 ) &&
         ( ( status->giveDamageTension() == 0. ) && ( status->giveDamageCompression() == 0. ) ) ) {
        //Rate factor from last step
        rateFactor = status->giveRateFactor();

        //Recalculate the f's
        tempEquivStrainTension = status->giveEquivStrainTension() + ( tempEquivStrain - status->giveEquivStrain() ) / rateFactor;
        tempEquivStrainCompression = status->giveEquivStrainCompression() + ( tempAlpha * tempEquivStrain - status->giveAlpha() * status->giveEquivStrain() ) / rateFactor;

        fTension = tempEquivStrainTension - status->giveKappaDTension();
        fCompression = tempEquivStrainCompression - status->giveKappaDCompression();
    }


    status->letTempRateFactorBe(rateFactor);

    //Normalize the fs
    fTension = fTension / e0;
    fCompression = fCompression / e0;

    double ductilityMeasure = computeDuctilityMeasureDamage(strain, gp);
    double deltaPlasticStrainNormTension, deltaPlasticStrainNormCompression;

    if ( fTension <= 0. && fCompression <= 0. ) {
        //Neither tension nor compression is active

        tempKappaDTension = status->giveKappaDTension();
        tempKappaDTensionOne = status->giveKappaDTensionOne();
        tempKappaDTensionTwo = status->giveKappaDTensionTwo();

        tempKappaDCompression = status->giveKappaDCompression();
        tempKappaDCompressionOne = status->giveKappaDCompressionOne();
        tempKappaDCompressionTwo = status->giveKappaDCompressionTwo();

        tempDamageTension = status->giveDamageTension();
        tempDamageCompression = status->giveDamageCompression();
    } else if ( fTension > 0. && fCompression <= 0. )      { //Only tension is active
        //Update tension history variables
        tempKappaDTension = tempEquivStrainTension;
        deltaPlasticStrainNorm = computeDeltaPlasticStrainNormTension(tempKappaDTension, status->giveKappaDTension(), gp);
        tempKappaDTensionOne = status->giveKappaDTensionOne() + deltaPlasticStrainNorm / ductilityMeasure / rateFactor;
        tempKappaDTensionTwo = status->giveKappaDTensionTwo() + ( tempKappaDTension - status->giveKappaDTension() ) / ductilityMeasure;

        //Nothing changes for compression history variables
        tempKappaDCompression = status->giveKappaDCompression();
        tempKappaDCompressionOne = status->giveKappaDCompressionOne();
        tempKappaDCompressionTwo = status->giveKappaDCompressionTwo();

        //Initialise damage with tensile history variable
        this->initDamaged(tempKappaDTension, strain, gp);

        tempDamageTension = computeDamageParam(tempKappaDTension, tempKappaDTensionOne, tempKappaDTensionTwo, gp);
        tempDamageCompression = status->giveDamageCompression();
    } else if ( fTension <= 0. && fCompression > 0. )      {
        //Only compression is active

        //Nothing changes for the history variables in tension
        tempKappaDTension = status->giveKappaDTension();
        tempKappaDTensionOne = status->giveKappaDTensionOne();
        tempKappaDTensionTwo = status->giveKappaDTensionTwo();

        //Update compression history variables
        tempKappaDCompression = tempEquivStrainCompression;
        deltaPlasticStrainNormCompression = computeDeltaPlasticStrainNormCompression(tempAlpha, tempKappaDCompression, status->giveKappaDCompression(), gp);
        tempKappaDCompressionOne = status->giveKappaDCompressionOne() + deltaPlasticStrainNormCompression / ( ductilityMeasure * rateFactor );
        tempKappaDCompressionTwo = status->giveKappaDCompressionTwo() + ( tempKappaDCompression - status->giveKappaDCompression() ) / ductilityMeasure;

        //Determine damage parameters
        tempDamageTension = status->giveDamageTension();
        tempDamageCompression = computeDamageParam(tempKappaDCompression, tempKappaDCompressionOne, tempKappaDCompressionTwo, gp);
    } else if ( fTension > 0. && fCompression > 0. )      {
        //Both tension and compression is active

        //Update tension history variables
        tempKappaDTension = tempEquivStrainTension;
        deltaPlasticStrainNormTension = computeDeltaPlasticStrainNormTension(tempKappaDTension, status->giveKappaDTension(), gp);
        tempKappaDTensionOne = status->giveKappaDTensionOne() + deltaPlasticStrainNormTension / ( ductilityMeasure * rateFactor );
        tempKappaDTensionTwo = status->giveKappaDTensionTwo() + ( tempKappaDTension - status->giveKappaDTension() ) / ductilityMeasure;

        //Update the compression history variables
        tempKappaDCompression = tempEquivStrainCompression;
        deltaPlasticStrainNormCompression =
            computeDeltaPlasticStrainNormCompression(tempAlpha, tempKappaDTension, status->giveKappaDTension(), gp);
        tempKappaDCompressionOne = status->giveKappaDCompressionOne() + deltaPlasticStrainNormCompression / ( ductilityMeasure * rateFactor );
        tempKappaDCompressionTwo = status->giveKappaDCompressionTwo() +
                                   ( tempKappaDCompression - status->giveKappaDCompression() ) / ductilityMeasure;

        //Determine the damage parameters
        this->initDamaged(tempKappaDTension, strain, gp);

        tempDamageTension = computeDamageParam(tempKappaDTension, tempKappaDTensionOne, tempKappaDTensionTwo, gp);
        tempDamageCompression = computeDamageParam(tempKappaDCompression, tempKappaDCompressionOne, tempKappaDCompressionTwo, gp);
    }

    //Write all temp history variables to the status

    status->letTempEquivStrainBe(tempEquivStrain);

    //Tension
    status->letTempEquivStrainTensionBe(tempEquivStrainTension);
    status->letTempKappaDTensionBe(tempKappaDTension);
    status->letTempKappaDTensionOneBe(tempKappaDTensionOne);
    status->letTempKappaDTensionTwoBe(tempKappaDTensionTwo);
    status->letTempDamageTensionBe(tempDamageTension);

    //Compression
    status->letTempEquivStrainCompressionBe(tempEquivStrainCompression);
    status->letTempKappaDCompressionBe(tempKappaDCompression);
    status->letTempKappaDCompressionOneBe(tempKappaDCompressionOne);
    status->letTempKappaDCompressionTwoBe(tempKappaDCompressionTwo);
    status->letTempDamageCompressionBe(tempDamageCompression);

    answer.at(1) = tempDamageTension;
    answer.at(2) = tempDamageCompression;
}







double
ConcreteDPM2 :: computeRateFactor(const double alpha,
                                  const double deltaTime,
                                  GaussPoint *gp,
                                  TimeStep *atTime)
{
    if ( strainRateFlag == 0 ) {
        return 1;
    }

    ConcreteDPM2Status *status = ( ConcreteDPM2Status * ) this->giveStatus(gp);

    //Access old and new strain
    StrainVector oldStrain(status->giveStrainVector(), matMode);
    StrainVector strain(status->giveTempStrainVector(), matMode);

    //Determine the principal values of the strain
    FloatArray principalStrain;
    strain.computePrincipalValues(principalStrain);
    //Determine max and min value;
    double maxStrain = -1.e20, minStrain = 1.e20;
    for ( int k = 1; k <= principalStrain.giveSize(); k++ ) {
        //maximum
        if ( principalStrain.at(k) > maxStrain ) {
            maxStrain = principalStrain.at(k);
        }

        //minimum
        if ( principalStrain.at(k) < minStrain ) {
            minStrain = principalStrain.at(k);
        }
    }

    //Evaluate the equivalent strains
    double strainRate;
    double oldRateStrain =   status->giveRateStrain();
    if ( 1. - alpha > DYNCON_TOL ) { //Tension
        strainRate = ( maxStrain - oldRateStrain ) / deltaTime;
        status->letTempRateStrainBe(maxStrain);
    } else   { //Compression
        strainRate = ( minStrain - oldRateStrain ) / deltaTime;
        status->letTempRateStrainBe(minStrain);
    }

    printf("strainRate = %e\n", strainRate);

    //For tension
    double deltaS = 1. / ( 1. + 8. * this->fc / this->fcZero );
    double betaS = exp( ( 6 * deltaS - 2. ) * log(10.) );


    //This factor is different from the Model code.
    //The reason is that we use the strain rate of the equivalent strain.
    //1) this strain rate is always positive.

    double strainRateZeroTension = 1.e-6;
    double strainRateZeroCompression = -30.e-6;

    //For compression
    // Fip-Model code 1990 expressions modified to take into account that we have an equivalent strain!
    double alphaS = 1. / ( 5. + 9 * this->fc / this->fcZero );
    double gammaS = exp( ( 6.156 * alphaS - 2. ) * log(10.0) );

    double strainRateRatioTension, strainRateRatioCompression;

    double rateFactorTension = 1.;
    double rateFactorCompression = 1.;

    strainRateRatioTension = strainRate / strainRateZeroTension;

    //Tension
    if ( strainRate < 30.e-6 ) {
        rateFactorTension = 1.;
    } else if ( 30.e-6 < strainRate && strainRate < 1 )          {
        rateFactorTension = pow(strainRateRatioTension, deltaS);
    } else   {
        rateFactorTension =  betaS * pow(strainRateRatioTension, 0.333);
    }

    //Compression
    strainRateRatioCompression = strainRate / strainRateZeroCompression;
    if ( strainRate > -30.e-6 ) {
        rateFactorCompression = 1.;
    } else if ( -30.e-6 > strainRate && strainRate > -30 )      {
        rateFactorCompression = pow(strainRateRatioCompression, 1.026 * alphaS);
    } else   {
        rateFactorCompression =  gammaS * pow(strainRateRatioCompression, 0.333);
    }

    double rateFactor = ( 1. - alpha ) * rateFactorTension + alpha * rateFactorCompression;

    return rateFactor;
}



double
ConcreteDPM2 :: computeDeltaPlasticStrainNormTension(double tempKappaD, double kappaD, GaussPoint *gp)
{
    if ( matMode == _Unknown ) {
        matMode = gp->giveMaterialMode();
    }

    ConcreteDPM2Status *status = ( ConcreteDPM2Status * ) this->giveStatus(gp);

    StrainVector tempPlasticStrain(matMode);
    StrainVector plasticStrain(matMode);
    status->giveTempPlasticStrain(tempPlasticStrain);
    status->givePlasticStrain(plasticStrain);

    StrainVector deltaPlasticStrain(matMode);
    deltaPlasticStrain = tempPlasticStrain;
    deltaPlasticStrain.subtract(plasticStrain);

    double deltaPlasticStrainNorm = 0;

    //Distinguish pre-peak, peak and post-peak
    double factor = 0.;
    if ( tempKappaD < e0 ) {
        deltaPlasticStrainNorm = 0.;
    } else if ( tempKappaD > e0 && kappaD < e0 )       {
        factor = ( 1. - ( e0 - kappaD ) / ( tempKappaD - kappaD ) );
        deltaPlasticStrain.times(factor);
        deltaPlasticStrainNorm = deltaPlasticStrain.computeNorm();
    } else   {
        deltaPlasticStrainNorm = deltaPlasticStrain.computeNorm();
    }

    return deltaPlasticStrainNorm;
}



double
ConcreteDPM2 :: computeDeltaPlasticStrainNormCompression(const double tempAlpha, double tempKappaD, double kappaD, GaussPoint *gp)
{
    if ( matMode == _Unknown ) {
        matMode = gp->giveMaterialMode();
    }

    ConcreteDPM2Status *status = ( ConcreteDPM2Status * ) this->giveStatus(gp);

    StrainVector tempPlasticStrain(matMode);
    StrainVector plasticStrain(matMode);
    status->giveTempPlasticStrain(tempPlasticStrain);
    status->givePlasticStrain(plasticStrain);

    plasticStrain.times( status->giveAlpha() );
    tempPlasticStrain.times(tempAlpha);
    StrainVector deltaPlasticStrain(matMode);
    deltaPlasticStrain = tempPlasticStrain;
    deltaPlasticStrain.subtract(plasticStrain);

    double deltaPlasticStrainNorm = 0;

    //Distinguish pre-peak, peak and post-peak
    double factor = 0.;
    if ( tempKappaD < e0 ) {
        deltaPlasticStrainNorm = 0.;
    } else if ( tempKappaD > e0 && kappaD < e0 )       {
        factor = ( 1. - ( e0 - kappaD ) / ( tempKappaD - kappaD ) );
        deltaPlasticStrain.times(factor);
        deltaPlasticStrainNorm = deltaPlasticStrain.computeNorm();
    } else   {
        deltaPlasticStrainNorm = deltaPlasticStrain.computeNorm();
    }

    double tempKappaP = status->giveTempKappaP();
    const double yieldHardTwo = computeHardeningTwo(tempKappaP);
    double extraFactor = this->ft * yieldHardTwo * sqrt(2. / 3.) / this->rho / sqrt( 1. + 2. * pow(this->dilationConst, 2.) );

    return deltaPlasticStrainNorm * extraFactor;
}


void
ConcreteDPM2 :: computeEquivalentStrain(double &tempEquivStrain,
                                        const StrainVector &strain,
                                        GaussPoint *gp,
                                        TimeStep *atTime)
{
    //compute the elastic strain based one
    ConcreteDPM2Status *status = ( ConcreteDPM2Status * ) this->giveStatus(gp);

    //get temp plastic strain and tempKappa
    StrainVector tempPlasticStrain(matMode);
    status->giveTempPlasticStrain(tempPlasticStrain);

    StressVector elasticStress(matMode);
    // compute elastic strains and trial stress
    StrainVector elasticStrain = strain;
    elasticStrain.subtract(tempPlasticStrain);
    elasticStrain.applyElasticStiffness(elasticStress, eM, nu);

    //Compute trial coordinates
    double sigElastic, rhoElastic, thetaElastic;
    computeTrialCoordinates(elasticStress, sigElastic, rhoElastic, thetaElastic);

    double rFunction = ( 4. * ( 1. - pow(this->ecc, 2.) ) * pow(cos(thetaElastic), 2.) + pow(2. * this->ecc - 1., 2.) ) / ( 2. * ( 1. - pow(this->ecc, 2.) ) * cos(thetaElastic) + ( 2. * this->ecc - 1. ) * sqrt(4. * ( 1. - pow(this->ecc, 2.) ) * pow(cos(thetaElastic), 2.) + 5. * pow(this->ecc, 2.) - 4. * this->ecc) );

    double pHelp = -this->m * ( rhoElastic * rFunction / ( sqrt(6.) * fc ) + sigElastic / fc );

    double qHelp = -3. / 2. * pow(rhoElastic, 2.) / pow(this->fc, 2.);

    double help = -0.5 * pHelp + sqrt(pow(pHelp, 2.) / 4. - qHelp);

    //negative help values are not of interest and create problems since we might compute the square root of something
    if ( help > 0 ) {
        tempEquivStrain = help * e0;
    } else   {
        tempEquivStrain = 0.;
    }
}



double
ConcreteDPM2 :: computeDamageParam(double equivStrain, double kappaOne, double kappaTwo, GaussPoint *gp)
{
    ConcreteDPM2Status *status = ( ConcreteDPM2Status * ) this->giveStatus(gp);
    double omega;
    double le = status->giveLe();
    double help;

    if ( equivStrain > this->ft / this->eM ) {
        if ( softeningType == 0 ) { //linear
            omega = ( this->eM * equivStrain * this->wf - this->ft * this->wf + this->ft * kappaOne * le ) /
                    ( this->eM * equivStrain * this->wf - this->ft * le * kappaTwo );
        } else   { //bilinear: Calculate damage parameter for both parts of bilinear curve  and check which fulfils limits.
            omega = ( this->eM * equivStrain * this->wfOne - this->ft * this->wfOne - ( this->ftOne - this->ft ) * kappaOne * le ) /
                    ( this->eM * equivStrain * this->wfOne + ( this->ftOne - this->ft ) * le * kappaTwo );
            help = le * kappaOne + le * omega * kappaTwo;

            if ( help >= 0. && help < this->wfOne ) {
                return omega;
            }

            omega = ( this->eM * equivStrain * ( this->wf - this->wfOne ) - this->ftOne * ( this->wf - this->wfOne ) +
                      this->ftOne * kappaOne * le  - this->ftOne * this->wfOne ) /
                    ( this->eM * equivStrain * ( this->wf - this->wfOne )  - this->ftOne * le * kappaTwo );
            help = le * kappaOne + le * omega * kappaTwo;

            if ( help > this->wfOne && help < this->wf ) {
                return omega;
            }
        }
    } else   {
        omega = 0.;
    }

    if ( omega > 1 ) {
        omega = 1.;
        return omega;
    } else if ( omega < 0. )       {
        OOFEM_ERROR("ConcreteDPM2 :: computeDamageParam - omega is smaller than zero. Not possible\n");
    }

    return omega;
}


void
ConcreteDPM2 :: initDamaged(double kappaD,
                            const StrainVector &strain,
                            GaussPoint *gp)
{
    int i, indx = 1;
    double le;
    FloatArray principalStrains, crackPlaneNormal(3);
    FloatMatrix principalDir(3, 3);
    ConcreteDPM2Status *status = ( ConcreteDPM2Status * ) this->giveStatus(gp);

    if ( kappaD <= 0 ) {
        return;
    }

    if ( helem > 0. ) {
        status->setLe(this->helem);
    } else if ( ( status->giveDamageTension() == 0. ) && ( status->giveDamageCompression() == 0. ) )         {
        strain.computePrincipalValDir(principalStrains, principalDir);
        // finfd index of max positive principal strain
        for ( i = 2; i <= 3; i++ ) {
            if ( principalStrains.at(i) > principalStrains.at(indx) ) {
                indx = i;
            }
        }

        for ( i = 1; i <= 3; i++ ) {
            crackPlaneNormal.at(i) = principalDir.at(i, indx);
        }

        //This gives the element size
        le = gp->giveElement()->giveLenghtInDir(crackPlaneNormal);

        // remember le in cooresponding status
        status->setLe(le);
    }
}


double
ConcreteDPM2 :: computeDuctilityMeasureDamage(const StrainVector &strain, GaussPoint *gp)
{
    ConcreteDPM2Status *status = ( ConcreteDPM2Status * ) this->giveStatus(gp);
    StrainVector plasticStrain(matMode);
    StrainVector tempPlasticStrain(matMode);
    status->giveTempPlasticStrain(tempPlasticStrain);
    status->givePlasticStrain(plasticStrain);
    tempPlasticStrain.subtract(plasticStrain);
    StrainVector principalStrain(matMode);
    double ductilityMeasure;

    //Angle in uniaxial compression is atan(1./sqrt(6.))=0.387597
    double alphaZero = 0.40824829;

    double Rs = 0;
    if ( this->sig < 0. ) {
        if ( this->rho > 1.e-16 ) {
            Rs = -this->sig / ( alphaZero * this->rho );
        } else   { //Need to set a dummy valye
            Rs = -this->sig * 1.e16 / alphaZero;
        }
    } else   {
        Rs = 0;
    }

    //  printf("Rs=%e\n", Rs);

    //  if(Rs <1.0){
    ductilityMeasure = 1. + ( ASoft - 1. ) * pow(Rs, this->BSoft);

    return ductilityMeasure;
}


void
ConcreteDPM2 :: performPlasticityReturn(GaussPoint *gp,
                                        StrainVector &strain)
{
    ConcreteDPM2Status *status = ( ConcreteDPM2Status * ) this->giveStatus(gp);
    if ( matMode == _Unknown ) {
        matMode = gp->giveMaterialMode();
    }

    //get temp plastic strain and tempKappa
    StrainVector tempPlasticStrain(matMode);
    status->giveTempPlasticStrain(tempPlasticStrain);
    double tempKappaP = status->giveTempKappaP();

    // compute elastic strains and trial stress
    StressVector effectiveStress(matMode);
    StrainVector elasticStrain = strain;
    elasticStrain.subtract(tempPlasticStrain);
    elasticStrain.applyElasticStiffness(effectiveStress, eM, nu);

    //Compute trial coordinates
    computeTrialCoordinates(effectiveStress, this->sig, this->rho, this->thetaTrial);

    double yieldValue = computeYieldValue(this->sig, this->rho, this->thetaTrial, tempKappaP);

    StrainVector convergedStrain(matMode);

    //What happens here to thermal strains that are subtracted first?
    StrainVector oldStrain(status->giveStrainVector(), matMode);
    StrainVector tempStrain(matMode);
    StrainVector deltaStrain(matMode);

    // introduce a strange subincrementation flag
    int subIncrementFlag = 0;

    double apexStress = 0.;
    // choose correct stress return and update state flag
    if ( yieldValue  > yieldTol ) {
        apexStress = 0.;
        checkForVertexCase(apexStress, sig, tempKappaP);

        //Make the appropriate return
        if ( returnType == RT_Tension || returnType == RT_Compression ) {
            performVertexReturn(effectiveStress, apexStress, gp);
            if ( returnType == RT_Regular ) {
                //This was no real vertex case
                //get the original tempKappaP and stress
                tempKappaP = status->giveKappaP();
                elasticStrain.applyElasticStiffness(effectiveStress, eM, nu);
            }
        }

        if ( returnType == RT_Regular ) {
            //Attempt to implement subincrementation
            // initialize variables
            subIncrementFlag = 0;
            convergedStrain = oldStrain;
            tempStrain = strain;
            deltaStrain = oldStrain;
            deltaStrain.negated();
            deltaStrain.add(strain);
            //To get into the loop
            returnResult = RR_NotConverged;
            while ( returnResult == RR_NotConverged || subIncrementFlag == 1 ) {
                elasticStrain = tempStrain;
                elasticStrain.subtract(tempPlasticStrain);
                elasticStrain.applyElasticStiffness(effectiveStress, eM, nu);

                tempKappaP = performRegularReturn(effectiveStress, gp);

                if ( returnResult == RR_NotConverged ) {
                    OOFEM_LOG_INFO("Subincrementation required");
                    subIncrementFlag = 1;
                    deltaStrain.times(0.5);
                    tempStrain = convergedStrain;
                    tempStrain.add(deltaStrain);
                } else if ( returnResult == RR_Converged && subIncrementFlag == 1 )      {
                    effectiveStress.applyElasticCompliance(elasticStrain, eM, nu);
                    tempPlasticStrain = tempStrain;
                    tempPlasticStrain.subtract(elasticStrain);
                    status->letTempPlasticStrainBe(tempPlasticStrain);
                    status->letTempKappaPBe(tempKappaP);

                    subIncrementFlag = 0;
                    returnResult = RR_NotConverged;
                    convergedStrain = tempStrain;
                    deltaStrain = convergedStrain;
                    deltaStrain.negated();
                    deltaStrain.add(strain);
                    tempStrain = strain;
                } else   {
                    status->letTempKappaPBe(tempKappaP);
                }
            }
        }
    }

    // compute the plastic strains
    effectiveStress.applyElasticCompliance(elasticStrain, eM, nu);
    tempPlasticStrain = strain;
    tempPlasticStrain.subtract(elasticStrain);

    // update plastic strain
    status->letTempPlasticStrainBe(tempPlasticStrain);
}


bool
ConcreteDPM2 :: checkForVertexCase(double answer,
                                   const double sig,
                                   const double tempKappa)
{
    //Compute sigZero and compare with actual sig
    int l = 0;
    double FZero = 1.;
    double dFZeroDSigZero = 0.;
    double sigZero;
    const double yieldHardOne = computeHardeningOne(tempKappa);
    double apexCompression = 0.;

    //compressive apex
    if ( ( tempKappa < 1. ) && ( sig < 0. ) ) {
        sigZero = -15. * fc;
        while ( fabs(FZero) > yieldTol && l <= newtonIter ) {
            l++;
            FZero = pow( ( 1. - yieldHardOne ), 2. ) * pow( ( sigZero / fc ), 4. ) +
                    pow(yieldHardOne, 2.) * m * ( sigZero / fc ) - pow(yieldHardOne, 2.);

            dFZeroDSigZero = pow( ( 1. - yieldHardOne ), 2. ) * 4. * pow( ( sigZero / fc ), 3. ) / fc +
                             pow(yieldHardOne, 2.) * m / fc;

            sigZero = sigZero - FZero / dFZeroDSigZero;
        }

        if ( l < 15 && sigZero < 0. ) {
            apexCompression = sigZero;
        } else   {
            apexCompression = -15. * fc;
        }
    }

    if ( ( sig > 0. && tempKappa < 1. ) || ( sig > fc / m && tempKappa >= 1. ) ) {
        returnType = RT_Tension;
        answer = 0.;
    } else if ( sig < apexCompression )      {
        returnType = RT_Compression;
        answer = apexCompression;
    } else   {
        returnType = RT_Regular;
    }

    return false;
}



void
ConcreteDPM2 :: performVertexReturn(StressVector &effectiveStress,
                                    double apexStress,
                                    GaussPoint *gp)
{
    ConcreteDPM2Status *status = ( ConcreteDPM2Status * ) this->giveStatus(gp);
    StressVector deviatoricStressTrial(matMode);
    double sigTrial;
    StressVector stressTemp(matMode);
    StressVector deviatoricStress(matMode);
    double yieldValue = 0.;
    double yieldValueMid = 0.;
    double sig2 = 0.;
    double dSig;
    double sigMid;
    double sigAnswer;
    double ratioPotential;

    effectiveStress.computeDeviatoricVolumetricSplit(deviatoricStressTrial, sigTrial);
    const double rhoTrial = deviatoricStressTrial.computeSecondCoordinate();

    StrainVector deltaPlasticStrain(matMode);

    double tempKappaP = status->giveTempKappaP();
    const double kappaInitial = tempKappaP;

    if ( returnType == RT_Tension ) {
        sig2 = -0.1 * ft;
    } else if ( returnType == RT_Compression )      {
        sig2 = apexStress;
    }

    tempKappaP =
        computeTempKappa(kappaInitial, sigTrial, rhoTrial, sigTrial);

    yieldValue =
        computeYieldValue(sigTrial, 0., 0., tempKappaP);

    tempKappaP =
        computeTempKappa(kappaInitial, sigTrial, rhoTrial, sig2);

    yieldValueMid =
        computeYieldValue(sig2, 0., 0., tempKappaP);

    if ( yieldValue * yieldValueMid >= 0. ) {
        returnType = RT_Regular;
        return;
    }

    if ( yieldValue < 0.0 ) {
        dSig = sig2 - sigTrial;
        sigAnswer = sig2;
    } else   {
        dSig = sigTrial - sig2;
        sigAnswer = sig2;
    }

    for ( int j = 0; j < 250; j++ ) {
        dSig = 0.5 * dSig;

        sigMid = sigAnswer + dSig;


        tempKappaP =
            computeTempKappa(kappaInitial, sigTrial, rhoTrial, sigMid);

        yieldValueMid = computeYieldValue(sigMid, 0., 0., tempKappaP);

        if ( yieldValueMid <= 0. ) {
            sigAnswer = sigMid;
        }

        if ( fabs(yieldValueMid) < yieldTol && yieldValueMid <= 0. ) {
            for ( int i = 0; i < 3; i++ ) {
                effectiveStress(i) = sigAnswer;
            }

            for ( int i = 3; i < effectiveStress.giveSize(); i++ ) {
                effectiveStress(i) = 0.;
            }

            //Put tempKappaP in the status
            status->letTempKappaPBe(tempKappaP);

            ratioPotential =
                computeRatioPotential(sigAnswer, tempKappaP);


            double ratioTrial = rhoTrial / ( sigTrial - sigAnswer );

            if ( ( ( ( ratioPotential >= ratioTrial ) && returnType == RT_Tension ) ) ||
                 ( ( ratioPotential <= ratioTrial ) && returnType == RT_Compression ) ) {
                return;
            } else   {
                returnType = RT_Regular;
                return;
            }
        }
    }

    returnType = RT_Regular;
}


double
ConcreteDPM2 :: computeTempKappa(const double kappaInitial,
                                 const double sigTrial,
                                 const double rhoTrial,
                                 const double sig)
{
    //This function is called, if stress state is in vertex case
    double equivalentDeltaPlasticStrain;
    FloatArray deltaPlasticStrainPrincipal(3);
    rho = 0.;
    equivalentDeltaPlasticStrain = sqrt( 1. / 9. * pow( ( sigTrial - sig ) / ( kM ), 2. ) +
                                         pow(rhoTrial / ( 2. * gM ), 2.) );

    double thetaVertex = 0.;
    double ductilityMeasure = computeDuctilityMeasure(sig, rho, thetaVertex);

    return kappaInitial + equivalentDeltaPlasticStrain / ductilityMeasure;
}


double
ConcreteDPM2 :: computeDuctilityMeasure(const double sig,
                                        const double rho,
					const double theta)
{
    double thetaConst = pow(2. * cos(theta), 2.);
    double ductilityMeasure;
    double x = -( sig + fc / 3 ) / fc;
    if ( x < 0. ) {
        /*Introduce exponential help function which results in a smooth
         * transition. */
        double EHard = BHard - DHard;
        double FHard = ( BHard - DHard ) * CHard / ( AHard - BHard );
        ductilityMeasure = ( EHard * exp(x / FHard) + DHard ) / thetaConst;
    } else   {
        ductilityMeasure = ( AHard + ( BHard - AHard ) * exp( -x / ( CHard ) ) ) / thetaConst;
    }

    return ductilityMeasure;
}


double
ConcreteDPM2 :: computeRatioPotential(const double sig,
                                      const double tempKappa)
{
    //compute yieldHard and yieldSoft
    const double yieldHardOne = computeHardeningOne(tempKappa);
    const double yieldHardTwo = computeHardeningTwo(tempKappa);

    //Compute dilation parameter
    double AGParam = this->ft * yieldHardTwo * 3 / this->fc + m / 2;
    double BGParam =
        yieldHardTwo / 3. * ( 1. + this->ft / this->fc ) /
        ( log(AGParam) + log(this->dilationConst + 1.) - log(2 * this->dilationConst - 1.) - log(3. * yieldHardTwo + this->m / 2) );

    const double R = ( sig - ft / 3. * yieldHardTwo ) / fc / BGParam;
    const double mQ = AGParam * exp(R);

    const double Bl = sig / fc + rho / ( fc * sqrt(6.) );

    const double Al = ( 1. - yieldHardOne ) * pow(Bl, 2.) + sqrt(3. / 2.) * rho / fc;

    const double dgdsig = 4. * ( 1. - yieldHardOne ) / fc * Al * Bl + yieldHardOne * mQ / fc;

    const double dgdrho = Al / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. ) +
                          m * yieldHardOne / ( sqrt(6.) * fc );

    return dgdrho / dgdsig * 3. * ( 1. - 2. * nu ) / ( 1. + nu );
}


double
ConcreteDPM2 :: performRegularReturn(StressVector &effectiveStress,
                                     GaussPoint *gp)
{
    ConcreteDPM2Status *status = ( ConcreteDPM2Status * ) this->giveStatus(gp);

    MaterialMode matMode = gp->giveMaterialMode();
    StressVector trialStress(matMode), tempStress(matMode), deviatoricTrialStress(matMode);

    double yieldValue;

    FloatArray residuals(4), residualsNorm(4), dGDInv(2);
    FloatMatrix jacobian(4, 4), inverseOfJacobian(4, 4);
    FloatArray rVector(3), helpVector(3), helpVector2(3);
    double dKappaDDeltaLambda;
    FloatMatrix aMatrix(3, 3), aMatrixInv(3, 3);
    FloatArray deltaIncrement(3);
    deltaIncrement.zero();
    FloatArray deltaIncrementTemp(3);
    deltaIncrementTemp.zero();
    StrainVector tempPlasticStrain(matMode);
    StrainVector plasticStrain(matMode);

    double deltaLambda = 0.;
    double normOfResiduals = 0.;
    int iterationCount = 0;

    FloatArray jacobianTimesAnswerIncrement;
    FloatArray unknownsTrial;
    FloatArray residualsTrial;
    FloatArray deltaUnknowns;
    FloatArray jacobianTimesDeltaUnknowns;
    FloatArray unknowns(4);

    //Define stressVariables
    double trialSig, trialRho;

    trialStress = effectiveStress;

    //compute the principal directions of the stress
    FloatArray helpStress;
    FloatMatrix stressPrincipalDir;
    trialStress.computePrincipalValDir(helpStress, stressPrincipalDir);

    //compute invariants from stress state
    trialStress.computeDeviatoricVolumetricSplit(deviatoricTrialStress, trialSig);
    trialRho = deviatoricTrialStress.computeSecondCoordinate();

    sig = trialSig;
    rho = trialRho;

    // Write the old plastic strain to the temp ones
    status->giveTempPlasticStrain(plasticStrain);
    tempPlasticStrain = plasticStrain;

    // Do the same for kappa
    double kappaP = status->giveTempKappaP();
    double tempKappaP = kappaP;

    //initialise unknowns
    unknowns.at(1) = trialSig;
    unknowns.at(2) = trialRho;
    unknowns.at(3) = tempKappaP;
    unknowns.at(4) = 0.;

    // Look at the magnitudes of the residuals. You have to scale the yieldValue down.
    yieldValue = computeYieldValue(sig, rho, thetaTrial, tempKappaP);

    //initiate residuals
    residuals.zero();
    residuals.at(4) = yieldValue;

    normOfResiduals  = 1.; //just to get into the loop

    /* N.R. iteration for finding the correct plastic return which is found when the norm of the residuals are equal to zero*/

    while ( normOfResiduals > yieldTol ) {
        iterationCount++;
        if ( iterationCount == newtonIter ) {
            returnResult = RR_NotConverged;
            return 0.;
        }

        //Normalize residuals. Think about it more.
        residualsNorm.at(1) = residuals.at(1) / this->kM;
        residualsNorm.at(2) = residuals.at(2) / ( 2. * this->gM );
        residualsNorm.at(3) = residuals.at(3);
        residualsNorm.at(4) = residuals.at(4);

        normOfResiduals = residualsNorm.computeNorm();

        if ( isnan(normOfResiduals) ) {
            returnResult = RR_NotConverged;
            return 0.;
        }

        if ( normOfResiduals > yieldTol ) {
            // Test to run newton iteration using inverse of Jacobian
            computeJacobian(jacobian, sig, rho, tempKappaP, deltaLambda, gp);

            if ( computeInverseOfJacobian(inverseOfJacobian, jacobian) ) {
                returnResult = RR_NotConverged;
                return 0.;
            }

            deltaIncrement.beProductOf(inverseOfJacobian, residuals);
            deltaIncrement.negated();

            //compute trial values
            unknownsTrial = unknowns;
            residualsTrial = residuals;

            //compute Unknowns
            for ( int i = 0; i < 4; i++ ) {
                unknowns(i) = unknownsTrial(i) + deltaIncrement(i);
            }

            if ( unknowns.at(4) <= 0. ) { //Keep deltaLambda greater than zero!
                unknowns.at(4) = 0.;
            }

            if ( unknowns.at(2) <= 0. ) { //Keep rho greater than zero!
                unknowns.at(2) = 0.;
            }

            if ( unknowns.at(3) - kappaP <= 0. ) { //Keep deltaKappa greater than zero!
                unknowns.at(3) = kappaP;
            }

            //compute residuals
            sig = unknowns.at(1);
            rho  = unknowns.at(2);

            tempKappaP = unknowns.at(3);
            deltaLambda = unknowns.at(4);

            /* Compute the mVector holding the derivatives of the g function and the hardening function*/
            computeDGDInv(dGDInv, sig, rho, tempKappaP);
            dKappaDDeltaLambda = computeDKappaDDeltaLambda(sig, rho, tempKappaP);

            residuals.at(1) = sig - trialSig + this->kM *deltaLambda *dGDInv.at(1);
            residuals.at(2) = rho - trialRho + ( 2. * this->gM ) * deltaLambda * dGDInv.at(2);
            residuals.at(3) = -tempKappaP + kappaP + deltaLambda * dKappaDDeltaLambda;
            residuals.at(4) = computeYieldValue(sig, rho, thetaTrial, tempKappaP);
        }
    }

    returnResult = RR_Converged;

    StressVector stressPrincipal(_3dMat);
    StressVector stressTemp(_3dMat);
    stressPrincipal.zero();

    stressPrincipal(0) = sig + sqrt(2. / 3.) * rho * cos(thetaTrial);
    stressPrincipal(1) = sig + sqrt(2. / 3.) * rho * cos(thetaTrial - 2. * M_PI / 3.);
    stressPrincipal(2) = sig + sqrt(2. / 3.) * rho * cos(thetaTrial + 2. * M_PI / 3.);

    transformStressVectorTo(stressTemp, stressPrincipalDir, stressPrincipal, 1);

    if ( matMode == _PlaneStrain ) {
        effectiveStress(0) = stressTemp(0);
        effectiveStress(1) = stressTemp(1);
        effectiveStress(2) = stressTemp(2);
        effectiveStress(3) = stressTemp(5);
    } else   {
        effectiveStress = stressTemp;
    }

    return tempKappaP;
}



void
ConcreteDPM2 :: computeJacobian(FloatMatrix &answer,
                                const double sig,
                                const double rho,
                                const double kappa,
                                const double deltaLambda,
                                GaussPoint *gp)
{
    //Variables
    FloatArray dFDInv(2);
    computeDFDInv(dFDInv, sig, rho, kappa);

    FloatArray dGDInv(2);
    computeDGDInv(dGDInv, sig, rho, kappa);

    FloatMatrix dDGDDInv(2, 2);
    computeDDGDDInv(dDGDDInv, sig, rho, kappa);

    double dKappaDDeltaLambda;
    dKappaDDeltaLambda = computeDKappaDDeltaLambda(sig, rho, kappa);

    double dFDKappa;
    dFDKappa = computeDFDKappa(sig, rho, kappa);

    FloatArray dDGDInvDKappa(2);
    computeDDGDInvDKappa(dDGDInvDKappa, sig, rho, kappa);


    double dDKappaDDeltaLambdaDKappa;
    dDKappaDDeltaLambdaDKappa = computeDDKappaDDeltaLambdaDKappa(sig, rho, kappa);


    FloatArray dDKappaDDeltaLambdaDInv(2);
    computeDDKappaDDeltaLambdaDInv(dDKappaDDeltaLambdaDInv, sig, rho, kappa);


    /* Compute matrix*/
    answer.at(1, 1) = 1. + this->kM *deltaLambda *dDGDDInv.at(1, 1);
    answer.at(1, 2) = this->kM * deltaLambda * dDGDDInv.at(1, 2);
    answer.at(1, 3) = this->kM * deltaLambda * dDGDInvDKappa.at(1);
    answer.at(1, 4) = this->kM * dGDInv.at(1);
    /**/
    answer.at(2, 1) = 2. *this->gM *deltaLambda *dDGDDInv.at(2, 1);
    answer.at(2, 2) = 1. + 2. *this->gM *deltaLambda *dDGDDInv.at(2, 2);
    answer.at(2, 3) = 2. *this->gM *deltaLambda *dDGDInvDKappa.at(2);
    answer.at(2, 4) = 2. *this->gM *dGDInv.at(2);
    /**/
    answer.at(3, 1) = deltaLambda * dDKappaDDeltaLambdaDInv.at(1);
    answer.at(3, 2) = deltaLambda * dDKappaDDeltaLambdaDInv.at(2);
    answer.at(3, 3) = deltaLambda * dDKappaDDeltaLambdaDKappa - 1.;
    answer.at(3, 4) = dKappaDDeltaLambda;
    /**/
    answer.at(4, 1) = dFDInv.at(1);
    answer.at(4, 2) = dFDInv.at(2);
    answer.at(4, 3) = dFDKappa;
    answer.at(4, 4) = 0.;
}




double
ConcreteDPM2 :: computeYieldValue(const double sig,
                                  const double rho,
                                  const double theta,
                                  const double tempKappa) const
{
    //compute yieldHard
    const double yieldHardOne = computeHardeningOne(tempKappa);
    const double yieldHardTwo = computeHardeningTwo(tempKappa);


    //  compute elliptic function r
    const double rFunction = ( 4. * ( 1. - pow(ecc, 2.) ) * pow(cos(theta), 2.) +
                               pow( ( 2. * ecc - 1. ), 2. ) ) /
                             ( 2. * ( 1. - pow(ecc, 2.) ) * cos(theta) +
                               ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - pow(ecc, 2.) ) * pow(cos(theta), 2.)
                                                        + 5. * pow(ecc, 2.) - 4. * ecc) );

    //compute help function Al
    const double Al = ( 1. - yieldHardOne ) * pow( ( sig / fc + rho / ( sqrt(6.) * fc ) ), 2. ) +
                      sqrt(3. / 2.) * rho / fc;

    //Compute yield equation
    return pow(Al, 2.) +
           yieldHardOne * yieldHardTwo * m * ( sig / fc + rho * rFunction / ( sqrt(6.) * fc ) ) -
           pow(yieldHardOne, 2.) * pow(yieldHardTwo, 2.);
}


double
ConcreteDPM2 :: computeDFDKappa(const double sig,
                                const double rho,
                                const double tempKappa)
{
    const double theta = thetaTrial;
    //compute yieldHard and yieldSoft
    const double yieldHardOne = computeHardeningOne(tempKappa);
    const double yieldHardTwo = computeHardeningTwo(tempKappa);
    // compute the derivative of the hardening and softening laws
    const double dYieldHardOneDKappa = computeHardeningOnePrime(tempKappa);
    const double dYieldHardTwoDKappa = computeHardeningTwoPrime(tempKappa);

    //compute elliptic function r
    const double rFunction =
        ( 4. * ( 1. - ecc * ecc ) * cos(theta) * cos(theta) + ( 2. * ecc - 1. ) * ( 2. * ecc - 1. ) ) /
        ( 2 * ( 1. - ecc * ecc ) * cos(theta) + ( 2. * ecc - 1. ) *
          sqrt(4. * ( 1. - ecc * ecc ) * cos(theta) * cos(theta) + 5. * ecc * ecc - 4. * ecc) );

    //compute help functions Al, Bl
    const double Al = ( 1. - yieldHardOne ) * pow( ( sig / fc + rho / ( sqrt(6.) * fc ) ), 2. ) + sqrt(3. / 2.) * rho / fc;

    const double Bl = sig / fc + rho / ( fc * sqrt(6.) );

    const double dFDYieldHardOne = -2. *Al *pow(Bl, 2.)
                                   + yieldHardTwo * m * ( sig / fc + rho * rFunction / ( sqrt(6.) * fc ) ) - 2. *yieldHardOne *pow(yieldHardTwo, 2.);

    const double dFDYieldHardTwo = yieldHardOne * m * ( sig / fc + rho * rFunction / ( sqrt(6.) * fc ) ) - 2. *yieldHardTwo *pow(yieldHardOne, 2.);

    // compute dFDKappa
    double dFDKappa =  dFDYieldHardOne * dYieldHardOneDKappa +
                      dFDYieldHardTwo * dYieldHardTwoDKappa;



    /*
     * set dFDKappa to zero, if it becomes greater than zero.
     * dFDKappa can only be negative or zero in the converged state for
     * the case of hardenig and perfect plasticity. For trial stresses, however,
     * it might be psoitive, which may lead to convergency problems. Therefore,
     * it is set to zero in this cases.
     */
    if ( dFDKappa > 0. ) {
        dFDKappa = 0.;
    }

    return dFDKappa;
}

void
ConcreteDPM2 :: computeDFDInv(FloatArray &answer,
                              const double sig,
                              const double rho,
                              const double tempKappa) const
{
    const double theta = thetaTrial;

    //compute yieldHard
    const double yieldHardOne = computeHardeningOne(tempKappa);
    const double yieldHardTwo = computeHardeningTwo(tempKappa);

    //compute elliptic function r
    const double rFunction = ( 4. * ( 1. - ecc * ecc ) * cos(theta) * cos(theta) + ( 2. * ecc - 1. ) * ( 2. * ecc - 1. ) ) /
                             ( 2. * ( 1. - ecc * ecc ) * cos(theta) + ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - ecc * ecc ) * cos(theta) * cos(theta)
                                                                                               + 5. * ecc * ecc - 4. * ecc) );

    //compute help functions AL, BL
    const double AL = ( 1. - yieldHardOne ) * pow( ( sig / fc + rho / ( sqrt(6.) * fc ) ), 2. ) + sqrt(3. / 2.) * rho / fc;
    const double BL = sig / fc + rho / ( fc * sqrt(6.) );

    //compute dfdsig
    const double dfdsig = 4. * ( 1. - yieldHardOne ) / fc * AL * BL + yieldHardTwo * yieldHardOne * m / fc;

    //compute dfdrho
    const double dfdrho = AL / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * BL + 6. ) + rFunction * m * yieldHardTwo * yieldHardOne / ( sqrt(6.) * fc );

    answer(0) = dfdsig;
    answer(1) = dfdrho;
}


double
ConcreteDPM2 :: computeDKappaDDeltaLambda(const double sig,
                                          const double rho,
                                          const double tempKappa)
{
    //Variables
    double equivalentDGDStress;
    FloatArray dGDInv(2);
    computeDGDInv(dGDInv, sig, rho, tempKappa);
    FloatArray dGDStressPrincipal(3);

    equivalentDGDStress = sqrt( 1. / 3. * pow(dGDInv(0), 2.) +
                                pow(dGDInv(1), 2.) );

    double ductilityMeasure = computeDuctilityMeasure(sig, rho, this->thetaTrial);

    double dKappaDDeltaLambda = equivalentDGDStress / ductilityMeasure;
    return dKappaDDeltaLambda;
}



void
ConcreteDPM2 :: computeDDKappaDDeltaLambdaDInv(FloatArray &answer,
                                               const double sig,
                                               const double rho,
                                               const double tempKappa)
{
    //Variables
    double equivalentDGDStress;
    FloatArray dGDInv(2);
    FloatMatrix dDGDDInv(2, 2);
    FloatArray dGDStressPrincipal(3);
    FloatArray dEquivalentDGDStressDInv(2);
    FloatArray helpA(3);

    //Compute first and second derivative of plastic potential
    computeDGDInv(dGDInv, sig, rho, tempKappa);
    computeDDGDDInv(dDGDDInv, sig, rho, tempKappa);

    //Compute equivalentDGDStress
    equivalentDGDStress = sqrt( 1. / 3. * pow(dGDInv(0), 2.) +
                                pow(dGDInv(1), 2.) );

    //computeDuctilityMeasure
    double ductilityMeasure = computeDuctilityMeasure(sig, rho, this->thetaTrial);

    //Compute dEquivalentDGDStressDInv
    dEquivalentDGDStressDInv(0) =
        ( 2. / 3. * dGDInv(0) * dDGDDInv(0, 0) + 2. * dGDInv(1) * dDGDDInv(1, 0) ) / ( 2. * equivalentDGDStress );
    dEquivalentDGDStressDInv(1) =
        ( 2. / 3. * dGDInv(0) * dDGDDInv(0, 1) + 2. * dGDInv(1) * dDGDDInv(1, 1) ) / ( 2. * equivalentDGDStress );

    answer.zero();

    // compute the derivative of
    FloatArray dDuctilityMeasureDInv(2);
    computeDDuctilityMeasureDInv(dDuctilityMeasureDInv, sig, rho, tempKappa);

    answer(0) = ( dEquivalentDGDStressDInv(0) * ductilityMeasure - equivalentDGDStress * dDuctilityMeasureDInv(0) ) / pow(ductilityMeasure, 2.);

    answer(1) = ( dEquivalentDGDStressDInv(1) * ductilityMeasure - equivalentDGDStress * dDuctilityMeasureDInv(1) ) / pow(ductilityMeasure, 2.);
}



double
ConcreteDPM2 :: computeDDKappaDDeltaLambdaDKappa(const double sig,
                                                 const double rho,
                                                 const double tempKappa)
{
    //Variables
    double equivalentDGDStress;
    FloatArray dGDInv(2);
    FloatArray dDGDInvDKappa(2);
    FloatArray dGDStressPrincipal(3);
    FloatArray helpA(3);
    double dEquivalentDGDStressDKappa;

    //Compute first and second derivative of plastic potential
    computeDGDInv(dGDInv, sig, rho, tempKappa);
    computeDDGDInvDKappa(dDGDInvDKappa, sig, rho, tempKappa);

    equivalentDGDStress = sqrt( 1. / 3. * pow(dGDInv(0), 2.) +
                                pow(dGDInv(1), 2.) );

    //computeDuctilityMeasure
    double ductilityMeasure = computeDuctilityMeasure(sig, rho, this->thetaTrial);

    //Compute dEquivalentDGDStressDKappa
    dEquivalentDGDStressDKappa =
        ( 2. / 3. * dGDInv(0) * dDGDInvDKappa(0) + 2. * dGDInv(1) * dDGDInvDKappa(1) ) / ( 2. * equivalentDGDStress );

    // compute the derivative of
    double dDuctilityMeasureDKappa = 0.;

    double dDKappaDDeltaLambdaDKappa =
        ( dEquivalentDGDStressDKappa * ductilityMeasure -
          equivalentDGDStress * dDuctilityMeasureDKappa ) / pow(ductilityMeasure, 2.);

    return dDKappaDDeltaLambdaDKappa;
}


void
ConcreteDPM2 :: computeDDuctilityMeasureDInv(FloatArray &answer,
                                             const double sig,
                                             const double rho,
                                             const double tempKappa)
{

    double thetaConst = pow(2. * cos(thetaTrial), 2.);
    double x = ( -( sig + fc / 3 ) ) / fc;

    if ( x < 0. ) {
        double dXDSig = -1. / fc;
        /* Introduce exponential help function which results in a
         * smooth transition. */
        double EHard = BHard - DHard;
        double FHard = ( BHard - DHard ) * CHard / ( AHard - BHard );

        double dDuctilityMeasureDX = EHard / FHard *exp(x / FHard) / thetaConst;
        answer(0) = dDuctilityMeasureDX * dXDSig;
        answer(1) = 0.;
    } else   {
        double dXDSig = -1. / fc;
        double dDuctilityMeasureDX = -( BHard - AHard ) / ( CHard ) / thetaConst *exp( -x / ( CHard ) );
        answer(0) = dDuctilityMeasureDX * dXDSig;
        answer(1) = 0.;
    }
}



void
ConcreteDPM2 :: computeDGDInv(FloatArray &answer,
                              const double sig,
                              const double rho,
                              const double tempKappa)
{
    //compute yieldHard and yieldSoft
    const double yieldHardOne = computeHardeningOne(tempKappa);
    const double yieldHardTwo = computeHardeningTwo(tempKappa);

    //Compute dilation parameter
    double AGParam = this->ft * yieldHardTwo * 3 / this->fc + m / 2;
    double BGParam =
        yieldHardTwo / 3. * ( 1. + this->ft / this->fc ) /
        ( log(AGParam) + log(this->dilationConst + 1.) - log(2 * this->dilationConst - 1.) - log(3. * yieldHardTwo + this->m / 2) );

    const double R = ( sig - ft / 3. * yieldHardTwo ) / fc / BGParam;
    const double mQ = AGParam * exp(R);

    const double Bl = sig / fc + rho / ( fc * sqrt(6.) );

    const double Al = ( 1. - yieldHardOne ) * pow(Bl, 2.) + sqrt(3. / 2.) * rho / fc;

    const double dgdsig = 4. * ( 1. - yieldHardOne ) / fc * Al * Bl + yieldHardOne * mQ / fc;

    const double dgdrho = Al / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. ) +
                          m * yieldHardOne / ( sqrt(6.) * fc );

    answer(0) = dgdsig;
    answer(1) = dgdrho;
}



void
ConcreteDPM2 :: computeDDGDInvDKappa(FloatArray &answer,
                                     const double sig,
                                     const double rho,
                                     const double tempKappa)
{
    //Compute dilation parameter


    //compute yieldHard and yieldSoft
    const double yieldHardOne = computeHardeningOne(tempKappa);
    const double yieldHardTwo = computeHardeningTwo(tempKappa);

    const double dYieldHardOneDKappa = computeHardeningOnePrime(tempKappa);
    const double dYieldHardTwoDKappa = computeHardeningTwoPrime(tempKappa);

    //Compute dilation parameter
    double AGParam = this->ft * yieldHardTwo * 3 / this->fc + m / 2;
    double BGParam =
        yieldHardTwo / 3. * ( 1. + this->ft / this->fc ) /
        ( log(AGParam) + log(this->dilationConst + 1.) - log(2 * this->dilationConst - 1.) - log(3. * yieldHardTwo + this->m / 2) );

    const double R = ( sig - ft / 3. * yieldHardTwo ) / ( fc * BGParam );
    const double mQ = AGParam * exp(R);

    //Compute the derivative of mQ with respect to kappa

    //Derivative of AGParam
    double dAGParamDKappa = dYieldHardTwoDKappa * 3. * this->ft / this->fc;

    //Derivative of BGParam
    double BGParamTop = yieldHardTwo / 3. * ( 1. + this->ft / this->fc );
    double BGParamBottom = ( log(AGParam) + log(this->dilationConst + 1.) - log(2 * this->dilationConst - 1.) - log(3. * yieldHardTwo + this->m / 2) );
    double dBGParamTopDKappa = dYieldHardTwoDKappa / 3.;
    double dBGParamBottomDKappa = -3. * dYieldHardTwoDKappa / ( 3 * yieldHardTwo + m / 2. );
    double dBGParamDKappa = ( dBGParamTopDKappa * BGParamBottom - BGParamTop * dBGParamBottomDKappa ) / pow(BGParamBottom, 2.);

    //Derivative of R
    double RTop = ( sig - ft / 3. * yieldHardTwo );
    double RBottom = fc * BGParam;
    double dRTopDKappa = -this->ft / 3. * dYieldHardTwoDKappa;
    double dRBottomDKappa = this->fc * dBGParamDKappa;
    double dRDKappa = ( dRTopDKappa * RBottom - RTop * dRBottomDKappa ) / pow(RBottom, 2.);

    double dMQDKappa = dAGParamDKappa * exp(R) + AGParam *dRDKappa *exp(R);

    const double Bl = sig / fc + rho / ( fc * sqrt(6.) );

    const double Al = ( 1. - yieldHardOne ) * pow(Bl, 2.) + sqrt(3. / 2.) * rho / fc;

    const double dAlDYieldHard = -pow(Bl, 2.);

    const double dDGDSigDKappa =
        ( -4. * Al * Bl / fc + 4. * ( 1 - yieldHardOne ) / fc * dAlDYieldHard * Bl ) * dYieldHardOneDKappa +
        dYieldHardOneDKappa * mQ / fc + yieldHardOne * dMQDKappa / fc;

    const double dDGDRhoDKappa =
        ( dAlDYieldHard / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. ) -
          4. * Al / ( sqrt(6.) * fc ) * Bl + m / ( sqrt(6.) * fc ) ) * dYieldHardOneDKappa;

    answer(0) = dDGDSigDKappa;
    answer(1) = dDGDRhoDKappa;
}


void
ConcreteDPM2 :: computeDDGDDInv(FloatMatrix &answer,
                                const double sig,
                                const double rho,
                                const double tempKappa)
{
    //compute yieldHardOne and yieldSoft
    const double yieldHardOne = computeHardeningOne(tempKappa);
    const double yieldHardTwo = computeHardeningTwo(tempKappa);

    //CoQpute dilation parameter
    double AGParam = this->ft * yieldHardTwo * 3 / this->fc + m / 2;
    double BGParam =
        yieldHardTwo / 3. * ( 1. + this->ft / this->fc ) /
        ( log(AGParam) + log(this->dilationConst + 1.) - log(2 * this->dilationConst - 1.) - log(3. * yieldHardTwo + this->m / 2) );

    const double R = ( sig - ft / 3. * yieldHardTwo ) / fc / BGParam;

    const double dMQDSig = AGParam / ( BGParam * fc ) * exp(R);

    //compute help parameter Al and Bl and the corresponding derivatives
    const double Bl = sig / fc + rho / ( fc * sqrt(6.) );

    const double Al = ( 1. - yieldHardOne ) * pow(Bl, 2.) +
                      sqrt(3. / 2.) * rho / fc;

    const double dAlDSig = 2. * ( 1. - yieldHardOne ) * Bl / fc;
    const double dBlDSig = 1. / fc;

    const double dAlDRho = 2. * ( 1. - yieldHardOne ) * Bl / ( fc * sqrt(6.) ) + sqrt(3. / 2.) / fc;
    const double dBlDRho = 1. / ( fc * sqrt(6.) );

    //compute second derivatives of g
    const double ddgddSig = 4. * ( 1. - yieldHardOne ) / fc * ( dAlDSig * Bl + Al * dBlDSig ) +
                            yieldHardOne * dMQDSig / fc;

    const double ddgddRho = dAlDRho / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. ) +
                            Al * dBlDRho * 4. * ( 1. - yieldHardOne ) / ( sqrt(6.) * fc );

    const double ddgdSigdRho = 4. * ( 1. - yieldHardOne ) / fc * ( dAlDRho * Bl + Al * dBlDRho );

    const double ddgdRhodSig = dAlDSig / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. )
                               + Al / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * dBlDSig );

    answer(0, 0) = ddgddSig;
    answer(0, 1) = ddgdSigdRho;
    answer(1, 0) = ddgdRhodSig;
    answer(1, 1) = ddgddRho;
}



double
ConcreteDPM2 :: computeAlpha(StressVector &effectiveStressTension,
                             StressVector &effectiveStressCompression,
                             StressVector &effectiveStress)
{
    FloatMatrix stressPrincipalDir(3, 3);
    FloatArray help;
    StressVector principalStress(_3dMat);
    StressVector effectiveStressTest(_3dMat);
    effectiveStress.computePrincipalValDir(principalStress, stressPrincipalDir);

    //Split the principal values in a tension and a compression part
    StressVector principalStressTension(_3dMat);
    StressVector principalStressCompression(_3dMat);

    for ( int i = 1; i <= principalStress.giveSize(); i++ ) {
        if ( principalStress.at(i) >= 0 ) {
            principalStressTension.at(i) = principalStress.at(i);
        } else   {
            principalStressCompression.at(i) = principalStress.at(i);
        }
    }

    //Transform the tension and compression principal stresses back to the original coordinate system

    StressVector stressTemp(_3dMat);

    //Take care of type of stress state for tension
    transformStressVectorTo(stressTemp, stressPrincipalDir, principalStressTension, 1);
    if ( matMode == _PlaneStrain ) {
        effectiveStressTension(0) = stressTemp(0);
        effectiveStressTension(1) = stressTemp(1);
        effectiveStressTension(2) = stressTemp(2);
        effectiveStressTension(3) = stressTemp(5);
    } else   {
        effectiveStressTension = stressTemp;
    }

    //Take care of type of stress state for compression
    transformStressVectorTo(stressTemp, stressPrincipalDir, principalStressCompression, 1);
    if ( matMode == _PlaneStrain ) {
        effectiveStressCompression(0) = stressTemp(0);
        effectiveStressCompression(1) = stressTemp(1);
        effectiveStressCompression(2) = stressTemp(2);
        effectiveStressCompression(3) = stressTemp(5);
    } else   {
        effectiveStressCompression = stressTemp;
    }

    //Determine the two factors from the stress

    double squareNormOfPrincipalStress = 0.;
    for ( int i = 1; i <= principalStress.giveSize(); i++ ) {
        squareNormOfPrincipalStress += pow(principalStress.at(i), 2.);
    }

    double alphaTension = 0.;

    if ( squareNormOfPrincipalStress > 0 ) {
        for ( int i = 1; i <= principalStress.giveSize(); i++ ) {
            alphaTension += principalStressTension.at(i) *
                            ( principalStressTension.at(i) + principalStressCompression.at(i) ) / squareNormOfPrincipalStress;
        }
    }

    return 1. - alphaTension;
}




int
ConcreteDPM2 :: computeInverseOfJacobian(FloatMatrix &answer, const FloatMatrix &src)
// Receiver becomes inverse of given parameter src. If necessary, size is adjusted.
{
    int nRows = src.giveNumberOfRows();

    //gaussian elimination - slow but safe
    //
    int i, j, k;
    double piv, linkomb;
    FloatMatrix tmp = src;
    answer.zero();
    // initialize answer to be unity matrix;
    for ( i = 1; i <= nRows; i++ ) {
        answer.at(i, i) = 1.0;
    }

    // lower triangle elimination by columns
    for ( i = 1; i < nRows; i++ ) {
        piv = tmp.at(i, i);
        if ( fabs(piv) < 1.e-20 ) {
            //indication that return does not converge. Output flag for subincrementing
            return 1;
        }

        for ( j = i + 1; j <= nRows; j++ ) {
            linkomb = tmp.at(j, i) / tmp.at(i, i);
            for ( k = i; k <= nRows; k++ ) {
                tmp.at(j, k) -= tmp.at(i, k) * linkomb;
            }

            for ( k = 1; k <= nRows; k++ ) {
                answer.at(j, k) -= answer.at(i, k) * linkomb;
            }
        }
    }

    // upper triangle elimination by columns
    for ( i = nRows; i > 1; i-- ) {
        piv = tmp.at(i, i);
        for ( j = i - 1; j > 0; j-- ) {
            linkomb = tmp.at(j, i) / tmp.at(i, i);
            for ( k = i; k > 0; k-- ) {
                tmp.at(j, k) -= tmp.at(i, k) * linkomb;
            }

            for ( k = nRows; k > 0; k-- ) {
                answer.at(j, k) -= answer.at(i, k) * linkomb;
            }
        }
    }

    // diagonal scaling
    for ( i = 1; i <= nRows; i++ ) {
        for ( j = 1; j <= nRows; j++ ) {
            answer.at(i, j) /= tmp.at(i, i);
        }
    }

    return 0;
}

double
ConcreteDPM2 :: computeHardeningOne(const double kappa) const
{
    if ( kappa <= 0. ) {
        return yieldHardInitial;
    } else if ( kappa > 0. && kappa < 1. )      {
        return
            ( 1. - yieldHardInitial - yieldHardPrimePeak ) * pow(kappa, 3.)
            - ( 3. * ( 1. - yieldHardInitial ) - 3. * yieldHardPrimePeak ) * pow(kappa, 2.)
            + ( 3. * ( 1. - yieldHardInitial ) - 2. * yieldHardPrimePeak ) * kappa
            + yieldHardInitial;
    } else   {
        return 1.;
    }
}


double
ConcreteDPM2 :: computeHardeningOnePrime(const double kappa) const
{
    if ( kappa <= 0. ) {
        return 3. * ( 1 - yieldHardInitial ) - 2. * yieldHardPrimePeak;
    } else if ( kappa >= 0. && kappa < 1. )      {
        return
            3. * ( 1. - yieldHardInitial - yieldHardPrimePeak ) * pow(kappa, 2.)
            - 2. * ( 3. * ( 1. - yieldHardInitial ) - 3. * yieldHardPrimePeak ) * kappa
            + ( 3. * ( 1. - yieldHardInitial ) - 2. * yieldHardPrimePeak );
    } else   {
        return 0.;
    }
}


double
ConcreteDPM2 :: computeHardeningTwo(const double kappa) const
{
    if ( kappa <= 0. ) {
        return 1.;
    } else if ( kappa > 0. && kappa < 1. )      {
        return 1.;
    } else   {
        return 1. + ( kappa - 1. ) * yieldHardPrimePeak;
    }
}


double
ConcreteDPM2 :: computeHardeningTwoPrime(const double kappa) const
{
    if ( kappa <= 0. ) {
        return 0.;
    } else if ( kappa >= 0. && kappa < 1. ) {
        return 0.;
    } else   {
        return yieldHardPrimePeak;
    }
}



void
ConcreteDPM2 :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                              MatResponseForm form,
                                              MatResponseMode mode,
                                              GaussPoint *gp,
                                              TimeStep *atTime)
{
    if ( gp->giveMaterialMode() == _3dMat ||
         gp->giveMaterialMode() ==  _PlaneStrain ||
         gp->giveMaterialMode() == _3dRotContinuum  ) {
        if ( mode == ElasticStiffness ) {
            this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);
        } else if ( mode == SecantStiffness )      {
            computeSecantStiffness(answer, form, mode, gp, atTime);
        } else if ( mode == TangentStiffness )      {
            _error("Tangent stiffness not implemented. Use either elastic or secant stiffness.\n");
        }
    }
}

void
ConcreteDPM2 :: computeSecantStiffness(FloatMatrix &answer,
                                       MatResponseForm form,
                                       MatResponseMode mode,
                                       GaussPoint *gp,
                                       TimeStep *atTime)
{
    this->matMode = gp->giveMaterialMode();
    ConcreteDPM2Status *status = ( ConcreteDPM2Status * ) giveStatus(gp);

    double omegaTension;

    //  Damage parameters
    omegaTension = status->giveTempDamageTension();
    if ( omegaTension > 0.999999 ) {
        omegaTension = 0.999999;
    }

    this->giveLinearElasticMaterial()->giveCharacteristicMatrix(answer, form, mode, gp, atTime);


    if ( isotropicFlag == 1 ) {
        answer.times(1. - omegaTension);
        return;
    }


    StrainVector strain(status->giveTempStrainVector(), matMode);
    StrainVector tempPlasticStrain(matMode);
    status->giveTempPlasticStrain(tempPlasticStrain);
    StressVector effectiveStress(matMode);
    StrainVector elasticStrain = strain;
    elasticStrain.subtract(tempPlasticStrain);
    elasticStrain.applyElasticStiffness(effectiveStress, eM, nu);

    //Calculate the principal values of the effective stress
    StressVector fullEffectiveStress(_3dMat);
    effectiveStress.convertToFullForm(fullEffectiveStress);
    FloatMatrix stressPrincipalDir(3, 3);
    StressVector principalStress(_3dMat);
    fullEffectiveStress.computePrincipalValDir(principalStress, stressPrincipalDir);


    //exclude two special cases.
    if ( principalStress.containsOnlyZeroes() ) {
        return;
    }

    for ( int i = 1; i <= 3; i++ ) {
        if ( principalStress.at(i) < -DYNCON_TOL ) {
            return;
        }
    }


    answer.times(1. - omegaTension);
}



void
ConcreteDPM2 :: computeTrialCoordinates(const StressVector &stress, double &sigNew, double &rhoNew, double &thetaNew)
{
    StressVector deviatoricStress(matMode);
    stress.computeDeviatoricVolumetricSplit(deviatoricStress, sigNew);
    rhoNew = deviatoricStress.computeSecondCoordinate();
    thetaNew = deviatoricStress.computeThirdCoordinate();
}


void
ConcreteDPM2 :: assignStateFlag(GaussPoint *gp)
{
    ConcreteDPM2Status *status = ( ConcreteDPM2Status * ) this->giveStatus(gp);

    //Get kappaD from status to define state later on
    double damageTension = status->giveDamageTension();
    double damageCompression = status->giveDamageCompression();
    double tempDamageTension = status->giveTempDamageTension();
    double tempDamageCompression = status->giveTempDamageCompression();
    double kappaP = status->giveKappaP();
    double tempKappaP = status->giveTempKappaP();

    if ( tempKappaP > kappaP ) {
        if ( tempDamageTension > damageTension ||  tempDamageTension == 1. ||
             tempDamageCompression > damageCompression || tempDamageCompression == 1. ) {
            status->letTempStateFlagBe(ConcreteDPM2Status :: ConcreteDPM2_PlasticDamage);
        } else   {
            status->letTempStateFlagBe(ConcreteDPM2Status :: ConcreteDPM2_Plastic);
        }
    } else   {
        const int state_flag = status->giveStateFlag();
        if ( tempDamageTension > damageTension || tempDamageTension == 1. ||
             tempDamageCompression > damageCompression || tempDamageCompression == 1. ) {
            status->letTempStateFlagBe(ConcreteDPM2Status :: ConcreteDPM2_Damage);
        } else   {
            if ( state_flag == ConcreteDPM2Status :: ConcreteDPM2_Elastic ) {
                status->letTempStateFlagBe(ConcreteDPM2Status :: ConcreteDPM2_Elastic);
            } else   {
                status->letTempStateFlagBe(ConcreteDPM2Status :: ConcreteDPM2_Unloading);
            }
        }
    }
}

void
ConcreteDPM2 :: computeDRhoDStress(FloatArray &answer,
                                   const StressVector &stress) const
{
    int size = 6;
    //compute volumetric deviatoric split
    StressVector deviatoricStress(_3dMat);
    double volumetricStress;
    stress.computeDeviatoricVolumetricSplit(deviatoricStress, volumetricStress);
    double rho = deviatoricStress.computeSecondCoordinate();

    //compute the derivative of J2 with respect to the stress
    FloatArray dJ2DStress;
    dJ2DStress = deviatoricStress;
    for ( int i = 3; i < size; i++ ) {
        dJ2DStress(i) = deviatoricStress(i) * 2.0;
    }

    //compute the derivative of rho with respect to stress
    FloatArray dRhoDStress;
    dRhoDStress = dJ2DStress;
    dRhoDStress.times(1. / rho);

    answer = dRhoDStress;
}

void
ConcreteDPM2 :: computeDSigDStress(FloatArray &answer) const
{
    int size = 6;
    for ( int i = 0; i < 3; i++ ) {
        answer(i) = 1. / 3.;
    }

    for ( int i = 3; i < size; i++ ) {
        answer(i) = 0.;
    }
}


void
ConcreteDPM2 :: computeDDRhoDDStress(FloatMatrix &answer,
                                     const StressVector &stress) const

{
    int size = 6;

    //compute volumetric deviatoric split
    StressVector deviatoricStress(_3dMat);
    double volumetricStress;
    stress.computeDeviatoricVolumetricSplit(deviatoricStress, volumetricStress);
    double rho = deviatoricStress.computeSecondCoordinate();


    //compute first dericative of J2
    FloatArray dJ2dstress;
    dJ2dstress = deviatoricStress;
    for ( int i = 3; i < deviatoricStress.giveSize(); i++ ) {
        dJ2dstress(i) = deviatoricStress(i) * 2.;
    }

    //compute second derivative of J2
    FloatMatrix ddJ2ddstress(size, size);
    ddJ2ddstress.zero();
    for ( int i = 0; i < size; i++ ) {
        if ( i < 3 ) {
            ddJ2ddstress(i, i) = 2. / 3.;
        }

        if ( i > 2 ) {
            ddJ2ddstress(i, i) = 2.;
        }
    }

    ddJ2ddstress(0, 1) = -1. / 3.;
    ddJ2ddstress(0, 2) = -1. / 3.;
    ddJ2ddstress(1, 0) = -1. / 3.;
    ddJ2ddstress(1, 2) = -1. / 3.;
    ddJ2ddstress(2, 0) = -1. / 3.;
    ddJ2ddstress(2, 1) = -1. / 3.;

    //compute square of the first derivative of J2
    FloatMatrix dJ2DJ2(size, size);
    for ( int v = 0; v < size; v++ ) {
        for ( int w = 0; w < size; ++w ) {
            dJ2DJ2(v, w) = dJ2dstress(v) * dJ2dstress(w);
        }
    }

    //compute the second derivative of rho
    FloatMatrix ddRhoddStress;
    ddRhoddStress = ddJ2ddstress;
    ddRhoddStress.times(1. / rho);
    FloatMatrix help1;
    help1 = dJ2DJ2;
    help1.times( -1. / ( rho * rho * rho ) );
    ddRhoddStress.add(help1);
    answer = ddRhoddStress;
}

int
ConcreteDPM2 :: giveIPValue(FloatArray &answer,
                            GaussPoint *gp,
                            InternalStateType type,
                            TimeStep *atTime)
{
    ConcreteDPM2Status *status = ( ConcreteDPM2Status * ) this->giveStatus(gp);
    //  const ConcreteDPM2Status* status = giveStatus(gp) ;
    StrainVector plasticStrainVector(_3dMat);

    switch ( type ) {
    case IST_PlasticStrainTensor:
        answer.resize(6);
        status->giveFullPlasticStrainVector(plasticStrainVector);
        answer = plasticStrainVector;
        return 1;

    case IST_DamageTensor:
        answer.resize(2);
        answer.zero();
        answer.at(1) = status->giveDamageTension();
        answer.at(2) = status->giveDamageCompression();
        return 1;

    default:
        return StructuralMaterial :: giveIPValue(answer, gp, type, atTime);

    }
    return 0;
}

int
ConcreteDPM2 :: giveIPValueSize(InternalStateType type,
                                GaussPoint *gp)
{
    switch ( type ) {
    case IST_PlasticStrainTensor:
        return 6;

    case IST_DamageTensor:
        return 2;

    default:
        return StructuralMaterial :: giveIPValueSize(type, gp);

    }
}

int
ConcreteDPM2 :: giveIntVarCompFullIndx(IntArray &answer,
                                       InternalStateType type,
                                       MaterialMode mmode)
{
    switch ( type ) {
    case IST_PlasticStrainTensor:
        answer.resize(9);
        answer.zero();
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 3;
        answer.at(4) = 4;
        answer.at(5) = 5;
        answer.at(6) = 6;
        return 1;

    case IST_DamageTensor:
        answer.resize(9);
        answer.at(1) = 1;
        answer.at(2) = 2;
        return 1;

    default:
        return StructuralMaterial :: giveIntVarCompFullIndx(answer, type, mmode);

    }
}

InternalStateValueType
ConcreteDPM2 :: giveIPValueType(InternalStateType type)
{
    switch ( type ) {
    case IST_PlasticStrainTensor: // plastic strain tensor
        return ISVT_TENSOR_S3E;

    case IST_DamageTensor: // damage tensor used for internal variables
        return ISVT_TENSOR_S3;

    default:
        return StructuralMaterial :: giveIPValueType(type);

    }
}

MaterialStatus *
ConcreteDPM2 :: CreateStatus(GaussPoint *gp) const
{
    ConcreteDPM2Status *status =
        new  ConcreteDPM2Status(1, StructuralMaterial :: giveDomain(), gp);
    return status;
}
} //end of namespace
