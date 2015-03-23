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

#include "concretedpm2.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "../sm/Materials/structuralms.h"
#include "gausspoint.h"
#include "intarray.h"
#include "datastream.h"
#include "contextioerr.h"
#include "timestep.h"
#include "../sm/Materials/structuralmaterial.h"
#include "Materials/isolinearelasticmaterial.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "mathfem.h"
#include "classfactory.h"
#include <limits>

namespace oofem {
REGISTER_Material(ConcreteDPM2);

ConcreteDPM2Status :: ConcreteDPM2Status(int n, Domain *d, GaussPoint *gp) :
    StructuralMaterialStatus(n, d, gp),
    plasticStrain( gp->giveMaterialMode() ),
    tempPlasticStrain( gp->giveMaterialMode() )
{
    tempPlasticStrain.zero();
    plasticStrain.zero();
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

    damageTension = tempDamageTension = 0.;
    damageCompression = tempDamageCompression = 0.;

    deltaLambda = 0.;
    state_flag = temp_state_flag = ConcreteDPM2Status :: ConcreteDPM2_Elastic;
    rateFactor = 1.;
    rateStrain = tempRateStrain = 0.;

#ifdef keep_track_of_dissipated_energy
    stressWork = tempStressWork = 0.0;
    dissWork = tempDissWork = 0.0;
#endif
}

ConcreteDPM2Status :: ~ConcreteDPM2Status()
{ }

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

    tempDamageTension = damageTension;
    tempDamageCompression = damageCompression;

    temp_state_flag = state_flag;
    tempRateFactor = rateFactor;
    tempRateStrain = rateStrain;
#ifdef keep_track_of_dissipated_energy
    tempStressWork = stressWork;
    tempDissWork = dissWork;
#endif
}

void
ConcreteDPM2Status :: updateYourself(TimeStep *tStep)
{
    // Call corresponding function of the parent class to update
    // variables defined there.
    StructuralMaterialStatus :: updateYourself(tStep);

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

    damageTension = tempDamageTension;
    damageCompression = tempDamageCompression;

    state_flag = temp_state_flag;

    rateFactor = tempRateFactor;

    rateStrain = tempRateStrain;
#ifdef keep_track_of_dissipated_energy
    stressWork = tempStressWork;
    dissWork = tempDissWork;
#endif
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

    StrainVector inelasticStrainVector( gp->giveMaterialMode() ); ///@todo Is this material mode really correct? Shouldn't it be the full material mode?
    StructuralMaterial :: giveFullSymVectorForm( inelasticStrainVector, strainVector, gp->giveMaterialMode() );

    inelasticStrainVector.subtract(plasticStrainVector);
    inelasticStrainVector.times(damageTension);
    inelasticStrainVector.add(plasticStrainVector);
    inelasticStrainVector.times(le);

    fprintf(file, " plastic ");
    for ( auto &val : plasticStrainVector ) {
        fprintf( file, " %.10e", val );
    }

    fprintf(file, " inelastic ");
    for ( auto &val : inelasticStrainVector ) {
        fprintf( file, " %.10e", val );
    }

    fprintf(file, " equivStrain %.10e,", equivStrain);

    fprintf(file, " kappaDTension %.10e,", kappaDTension);

    fprintf(file, " kappaDCompression %.10e,", kappaDCompression);

    fprintf(file, " kappaP %.10e,", kappaP);

    fprintf(file, " kappaDTensionOne %.10e,", kappaDTensionOne);

    fprintf(file, " kappaDCompressionOne %.10e,", kappaDCompressionOne);

    fprintf(file, " kappaDTensionTwo %.10e,", kappaDTensionTwo);

    fprintf(file, " kappaDCompressionTwo %.10e,", kappaDCompressionTwo);

    fprintf(file, " damageTension %.10e,", damageTension);

    fprintf(file, " damageCompression %.10e,", damageCompression);

    fprintf(file, " alpha %.10e,", this->alpha);

#ifdef keep_track_of_dissipated_energy
    fprintf(file, " dissW %g, freeE %g, stressW %g ", this->dissWork, ( this->stressWork ) - ( this->dissWork ), this->stressWork);
#endif
    fprintf(file, "}\n");
}

contextIOResultType
ConcreteDPM2Status :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // save parent class status
    if ( ( iores = StructuralMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write raw data

    if ( ( iores = plasticStrain.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.write(kappaP) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(deltaLambda) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(dFDKappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(le) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(alpha) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(equivStrainTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(equivStrainCompression) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappaDTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappaDCompression) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappaDCompressionOne) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappaDTensionTwo) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(kappaDCompressionTwo) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    if ( !stream.write(damageTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(damageCompression) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(deltaEquivStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(rateFactor) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(rateStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    if ( !stream.write(state_flag) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#ifdef keep_track_of_dissipated_energy
    if ( !stream.write(stressWork) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(dissWork) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#endif
    return CIO_OK;
}


contextIOResultType
ConcreteDPM2Status :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    // read parent class status
    if ( ( iores = StructuralMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read raw data
    if ( ( iores = plasticStrain.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.read(kappaP) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(deltaLambda) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(dFDKappa) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(le) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(alpha) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(equivStrainTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(equivStrainCompression) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(kappaDTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(kappaDCompression) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(kappaDCompressionOne) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(kappaDTensionTwo) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(kappaDCompressionTwo) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    if ( !stream.read(damageTension) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(damageCompression) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(deltaEquivStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(rateFactor) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(rateStrain) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    if ( !stream.read(state_flag) ) {
        THROW_CIOERR(CIO_IOERR);
    }


#ifdef keep_track_of_dissipated_energy
    if ( !stream.read(stressWork) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(dissWork) ) {
        THROW_CIOERR(CIO_IOERR);
    }

#endif
    return CIO_OK;
}

#ifdef keep_track_of_dissipated_energy
void
ConcreteDPM2Status :: computeWork(GaussPoint *gp, double gf)
{
    MaterialMode matMode = gp->giveMaterialMode();
    StrainVector tempTotalstrain(tempStrainVector, matMode);
    StrainVector totalstrain(strainVector, matMode);

    //Calculate increase or decrease of total strain tensor during iteration/step
    StrainVector deltaTotalStrain(tempTotalstrain);
    deltaTotalStrain.subtract(totalstrain);

    //Ask for stress tensor at step
    StressVector stress(tempStressVector, matMode);
    int n = stress.giveSize();
    //Calculate increase/decrease in total work
    double dSW = ( tempStressVector.dotProduct(deltaTotalStrain, n) + stressVector.dotProduct(deltaTotalStrain, n) ) / 2.;

    double tempStressWork = this->giveTempStressWork() + dSW;

    //Calculate temporary elastic strain
    StrainVector tempElasticStrain(matMode);
    tempElasticStrain = tempTotalstrain;
    tempElasticStrain.subtract(tempPlasticStrain);

    //Calculate elastically stored energy density
    double We = tempStressVector.dotProduct(tempElasticStrain, n) / 2.;

    // dissipative work density
    tempDissWork = tempStressWork - We;

    // to avoid extremely small negative dissipation due to round-off error
    // (note: gf is the dissipation density at complete failure, per unit volume)
    //Rough estimation of gf which is ok for this purpose

    if ( fabs(tempDissWork) < 1.e-12 * gf ) {
        tempDissWork = 0.;
    }

    this->setTempStressWork(tempStressWork);
    this->setTempDissWork(tempDissWork);
}
#endif

//   ********************************
//   *** CLASS DYNAMIC CONCRETE   ***
//   ********************************

#define IDM_ITERATION_LIMIT 1.e-8

ConcreteDPM2 :: ConcreteDPM2(int n, Domain *d) :
    StructuralMaterial(n, d)
{
    yieldTol = 0.;
    yieldTolDamage = 0.;
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
    IRResultType result;

    // call the corresponding service for the linear elastic material
    StructuralMaterial :: initializeFrom(ir);
    linearElasticMaterial->initializeFrom(ir);

    //isotropic flag
    isotropicFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, isotropicFlag, _IFT_ConcreteDPM2_isoflag);

    double value;
    // elastic parameters
    IR_GIVE_FIELD(ir, eM, _IFT_IsotropicLinearElasticMaterial_e)
    IR_GIVE_FIELD(ir, nu, _IFT_IsotropicLinearElasticMaterial_n);
    propertyDictionary.add('E', eM);
    propertyDictionary.add('n', nu);

    IR_GIVE_FIELD(ir, value, _IFT_IsotropicLinearElasticMaterial_talpha);
    propertyDictionary.add(tAlpha, value);

    gM = eM / ( 2. * ( 1. + nu ) );
    kM = eM / ( 3. * ( 1. - 2. * nu ) );

    IR_GIVE_FIELD(ir, fc, _IFT_ConcreteDPM2_fc);
    IR_GIVE_FIELD(ir, ft, _IFT_ConcreteDPM2_ft);

    this->e0 = this->ft / this->eM;

    // default parameters
    this->ecc = 0.525;
    IR_GIVE_OPTIONAL_FIELD(ir, ecc, _IFT_ConcreteDPM2_ecc);
    yieldHardInitial = 0.3;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldHardInitial, _IFT_ConcreteDPM2_kinit);
    //Inclination at transition point
    yieldHardPrimePeak = 0.5;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldHardPrimePeak, _IFT_ConcreteDPM2_hp);

    if ( yieldHardPrimePeak < 0 ) {
        yieldHardPrimePeak = 0.;
        OOFEM_WARNING("kPrimePeak cannot be less than zero");
    } else if ( yieldHardPrimePeak > ( 1. - yieldHardInitial ) ) {
        yieldHardPrimePeak = 1. - yieldHardInitial;
        OOFEM_WARNING("kPrimePeak cannot be greater than 1.-kinit");
    }

    AHard = 8.e-2;
    IR_GIVE_OPTIONAL_FIELD(ir, AHard, _IFT_ConcreteDPM2_ahard);
    BHard = 3.e-3;
    IR_GIVE_OPTIONAL_FIELD(ir, BHard, _IFT_ConcreteDPM2_bhard);
    CHard = 2.;
    IR_GIVE_OPTIONAL_FIELD(ir, CHard, _IFT_ConcreteDPM2_chard);
    DHard = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, DHard, _IFT_ConcreteDPM2_dhard);
    dilationConst = 0.85;
    IR_GIVE_OPTIONAL_FIELD(ir, dilationConst, _IFT_ConcreteDPM2_dilation);

    softeningType = 1; //0-Linear softening; 1-Bilinear softening; 2-Exponential
    IR_GIVE_OPTIONAL_FIELD(ir, softeningType, _IFT_ConcreteDPM2_softeningType);

    if ( softeningType > 2 ) {
        OOFEM_ERROR("softening type not implemented");
    }

    IR_GIVE_FIELD(ir, this->wf, _IFT_ConcreteDPM2_wf);

    if ( softeningType == 1 ) {
        this->ftOne = 0.3 * this->ft;
        IR_GIVE_OPTIONAL_FIELD(ir, ftOne, _IFT_ConcreteDPM2_ftOne);
        this->wfOne = 0.15 * this->wf;
        IR_GIVE_OPTIONAL_FIELD(ir, this->wfOne, _IFT_ConcreteDPM2_wfOne);
    }

    this->efCompression = 100.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, this->efCompression, _IFT_ConcreteDPM2_efc);

    this->ASoft = 15;
    IR_GIVE_OPTIONAL_FIELD(ir, ASoft, _IFT_ConcreteDPM2_asoft);

    helem = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, helem, _IFT_ConcreteDPM2_helem);


    //Compute m
    m = 3. * ( pow(fc, 2.) - pow(ft, 2.) ) / ( fc * ft ) * ecc / ( ecc + 1. );

    //Compute default value of dilationConst
    yieldTol = 1.e-6;
    IR_GIVE_OPTIONAL_FIELD(ir, yieldTol, _IFT_ConcreteDPM2_yieldtol);

    yieldTolDamage = yieldTol * 10.;

    newtonIter = 100;
    IR_GIVE_OPTIONAL_FIELD(ir, newtonIter, _IFT_ConcreteDPM2_newtoniter);

    //parameters for rate dependence
    strainRateFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, strainRateFlag, _IFT_ConcreteDPM2_rateFlag);

    timeFactor = 1.;
    if ( strainRateFlag == 1 ) {
        IR_GIVE_FIELD(ir, fcZero, _IFT_ConcreteDPM2_fcZero);
        IR_GIVE_OPTIONAL_FIELD(ir, timeFactor, _IFT_ConcreteDPM2_timeFactor);
    }

    return IRRT_OK;
}

int
ConcreteDPM2 :: hasMaterialModeCapability(MaterialMode mMode)
{
    return mMode == _3dMat || mMode == _PlaneStrain || mMode == _1dMat;
}

void
ConcreteDPM2 :: giveRealStressVector_1d(FloatArray &answer,
                                        GaussPoint *gp,
                                        const FloatArray &strainVector,
                                        TimeStep *atTime)
{
    if ( matMode == _Unknown ) {
        matMode = gp->giveMaterialMode();
    }

    ConcreteDPM2Status *status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    // Initialize temp variables for this gauss point
    //    status->initTempStatus();
    //this->initGpForNewStep(gp);
    this->initTempStatus(gp);
    status->letTempStrainVectorBe(strainVector);

    StrainVector strain(strainVector, matMode);

    //Calculate strain rate
    //Time step
    double deltaTime = 1.;
    if ( strainRateFlag == 1 ) {
        if ( atTime->giveTimeIncrement() == 0 ) { //Problem with the first step. For some reason the time increment is zero
            deltaTime = this->timeFactor;
        } else {
            deltaTime = atTime->giveTimeIncrement() * this->timeFactor;
        }
    } else {
        if ( atTime->giveTimeIncrement() == 0 ) {
            deltaTime = 1.;
        } else {
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
    elasticStrain.subtract( status->giveTempPlasticStrain() );
    elasticStrain.applyElasticStiffness(effectiveStress, eM, nu);


    StressVector effectiveStressTension(matMode);
    StressVector effectiveStressCompression(matMode);
    double alpha = 0.;
    if ( effectiveStress.at(1) >= 0 ) { //1D tensile stress state
        alpha = 0.;
        effectiveStressTension = effectiveStress;
        effectiveStressCompression.zero();
    } else if ( effectiveStress.at(1) < 0 ) {      //1D compressive stress state
        alpha = 1.;
        effectiveStressTension.zero();
        effectiveStressCompression = effectiveStress;
    }

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
    } else { //Consider only tensile damage. Reduction to a fully isotropic model
        stress = effectiveStress;
        stress.times( 1. - damages.at(1) );
    }

    status->letTempStressVectorBe(stress);
#ifdef keep_track_of_dissipated_energy
    double gf = pow(ft, 2) / this->eM; //rough estimation only for this purpose
    status->computeWork(gp, gf);
#endif
    assignStateFlag(gp);
    answer = stress;
}

void
ConcreteDPM2 :: giveRealStressVector(FloatArray &answer,
                                     GaussPoint *gp,
                                     const FloatArray &strainVector,
                                     TimeStep *tStep)
{
    if ( matMode == _Unknown ) {
        matMode = gp->giveMaterialMode();
    }

    ConcreteDPM2Status *status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    // Initialize temp variables for this gauss point
    status->initTempStatus();

    status->letTempStrainVectorBe(strainVector);

    //  ConcreteDPM2Status *status = giveStatus (gp) ;

    StrainVector strain(strainVector, matMode);

    //Calculate strain rate
    //Time step
    double deltaTime = 1.;
    if ( strainRateFlag == 1 ) {
        if ( tStep->giveTimeIncrement() == 0 ) { //Problem with the first step. For some reason the time increment is zero
            deltaTime = this->timeFactor;
        } else {
            deltaTime = tStep->giveTimeIncrement() * this->timeFactor;
        }
    } else {
        if ( tStep->giveTimeIncrement() == 0 ) {
            deltaTime = 1.;
        } else {
            deltaTime = tStep->giveTimeIncrement();
        }
    }

    // perform plasticity return
    performPlasticityReturn(gp, strain);

    // compute the nominal stress
    StressVector stress(matMode);
    StressVector effectiveStress(matMode);

    // compute elastic strains and trial stress
    StrainVector elasticStrain = strain;
    elasticStrain.subtract( status->giveTempPlasticStrain() );
    elasticStrain.applyElasticStiffness(effectiveStress, eM, nu);


    StressVector effectiveStressTension(matMode);
    StressVector effectiveStressCompression(matMode);
    double alpha;

    alpha = computeAlpha(effectiveStressTension, effectiveStressCompression, effectiveStress);

    status->letTempAlphaBe(alpha);


    FloatArray damages(2);
    damages.zero();
    computeDamage(damages, strain, deltaTime, gp, tStep, alpha);

    //Split damage in a tension and compression part

    if ( isotropicFlag == 0 ) { //Default
        effectiveStressTension.times( 1. - damages.at(1) );
        effectiveStressCompression.times( 1. - damages.at(2) );
        stress = effectiveStressTension;
        stress.add(effectiveStressCompression);
    } else { //Consider only tensile damage. Reduction to a fully isotropic model
        stress = effectiveStress;
        stress.times( 1. - damages.at(1) );
    }

    status->letTempStressVectorBe(stress);
#ifdef keep_track_of_dissipated_energy
    double gf = pow(ft, 2) / this->eM; //rough estimation only for this purpose
    status->computeWork(gp, gf);
#endif
    assignStateFlag(gp);

    answer = stress;
}


void
ConcreteDPM2 :: computeDamage(FloatArray &answer,
                              const StrainVector &strain,
                              const double deltaTime,
                              GaussPoint *gp,
                              TimeStep *tStep,
                              const double tempAlpha)
{
    ConcreteDPM2Status *status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    double tempEquivStrain;
    double deltaPlasticStrainNorm;
    double tempDamageTension = 0.0;
    double tempDamageCompression = 0.0;

    double tempKappaDTension = 0.0, tempKappaDCompression = 0.0;
    double tempKappaDTensionOne = 0.0, tempKappaDTensionTwo = 0.0;
    double tempKappaDCompressionOne = 0.0, tempKappaDCompressionTwo = 0.0;

    double rateFactor;

    int unAndReloadingFlag = 0;
    double minEquivStrain = 0.;

    unAndReloadingFlag = checkForUnAndReloading(tempEquivStrain, minEquivStrain, gp);

    if ( ( status->giveDamageTension() == 0. ) && ( status->giveDamageCompression() == 0. ) ) {
        rateFactor = computeRateFactor(tempAlpha, deltaTime, gp, tStep);
    } else {
        rateFactor = status->giveRateFactor();
    }

    //Compute equivalent strains for  tension and compression
    double tempEquivStrainTension = 0.;
    double tempEquivStrainCompression = 0.;

    tempEquivStrainTension = status->giveEquivStrainTension() + ( tempEquivStrain - status->giveEquivStrain() ) / rateFactor;

    if ( unAndReloadingFlag == 0 ) { //Standard way
        tempEquivStrainCompression = status->giveEquivStrainCompression() + ( tempAlpha * ( tempEquivStrain - status->giveEquivStrain() ) ) / rateFactor;
    } else {
        tempEquivStrainCompression = status->giveEquivStrainCompression() + status->giveAlpha() * ( minEquivStrain - status->giveEquivStrain() ) / rateFactor + ( tempAlpha * ( tempEquivStrain - minEquivStrain ) ) / rateFactor;
    }


    //If damage threshold is exceeded determine the rate factor from the previous step
    if ( ( tempEquivStrainTension > e0 || tempEquivStrainCompression > e0 ) &&
         ( ( status->giveDamageTension() == 0. ) && ( status->giveDamageCompression() == 0. ) ) && !tStep->isTheFirstStep() ) {
        //Rate factor from last step
        rateFactor = status->giveRateFactor();

        tempEquivStrainTension = status->giveEquivStrainTension() + ( tempEquivStrain - status->giveEquivStrain() ) / rateFactor;
        if ( unAndReloadingFlag == 0 ) { //Standard way
            tempEquivStrainCompression = status->giveEquivStrainCompression() + ( tempAlpha * ( tempEquivStrain - status->giveEquivStrain() ) ) / rateFactor;
        } else {
            tempEquivStrainCompression = status->giveEquivStrainCompression() + status->giveAlpha() * ( minEquivStrain - status->giveEquivStrain() ) / rateFactor + ( tempAlpha * ( tempEquivStrain - minEquivStrain ) ) / rateFactor;
        }
    }

    status->letTempRateFactorBe(rateFactor);

    double fTension = tempEquivStrainTension - status->giveKappaDTension();
    double fCompression = tempEquivStrainCompression - status->giveKappaDCompression();

    //Normalize the fs
    fTension = fTension / e0;
    fCompression = fCompression / e0;

    double ductilityMeasure = computeDuctilityMeasureDamage(strain, gp);
    double deltaPlasticStrainNormTension, deltaPlasticStrainNormCompression;

    if ( fTension < -yieldTolDamage && fCompression < -yieldTolDamage ) {
        //Neither tension nor compression is active

        tempKappaDTension = status->giveKappaDTension();
        tempKappaDTensionOne = status->giveKappaDTensionOne();
        tempKappaDTensionTwo = status->giveKappaDTensionTwo();

        tempKappaDCompression = status->giveKappaDCompression();
        tempKappaDCompressionOne = status->giveKappaDCompressionOne();
        tempKappaDCompressionTwo = status->giveKappaDCompressionTwo();

        tempDamageTension = status->giveDamageTension();
        tempDamageCompression = status->giveDamageCompression();
    } else if ( fTension >= -yieldTolDamage && fCompression < -yieldTolDamage ) {   //Only tension is active
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

        tempDamageTension = computeDamageParamTension( tempKappaDTension, tempKappaDTensionOne, tempKappaDTensionTwo, status->giveLe(), status->giveDamageTension() );

        tempDamageCompression = status->giveDamageCompression();
    } else if ( fTension < -yieldTolDamage && fCompression >= -yieldTolDamage ) {
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
        tempDamageCompression = computeDamageParamCompression( tempKappaDCompression, tempKappaDCompressionOne, tempKappaDCompressionTwo, status->giveDamageCompression() );
    } else if ( fTension >= -yieldTolDamage && fCompression >= -yieldTolDamage ) {
        //Both tension and compression is active

        //Update tension history variables
        tempKappaDTension = tempEquivStrainTension;
        deltaPlasticStrainNormTension = computeDeltaPlasticStrainNormTension(tempKappaDTension, status->giveKappaDTension(), gp);
        tempKappaDTensionOne = status->giveKappaDTensionOne() + deltaPlasticStrainNormTension / ( ductilityMeasure * rateFactor );
        tempKappaDTensionTwo = status->giveKappaDTensionTwo() + ( tempKappaDTension - status->giveKappaDTension() ) / ductilityMeasure;

        //Update the compression history variables
        tempKappaDCompression = tempEquivStrainCompression;
        deltaPlasticStrainNormCompression =
            computeDeltaPlasticStrainNormCompression(tempAlpha, tempKappaDCompression, status->giveKappaDCompression(), gp);
        tempKappaDCompressionOne = status->giveKappaDCompressionOne() + deltaPlasticStrainNormCompression / ( ductilityMeasure * rateFactor );
        tempKappaDCompressionTwo = status->giveKappaDCompressionTwo() +
                                   ( tempKappaDCompression - status->giveKappaDCompression() ) / ductilityMeasure;

        //Determine the damage parameters
        this->initDamaged(tempKappaDTension, strain, gp);

        tempDamageTension = computeDamageParamTension( tempKappaDTension, tempKappaDTensionOne, tempKappaDTensionTwo, status->giveLe(), status->giveDamageTension() );

        tempDamageCompression = computeDamageParamCompression( tempKappaDCompression, tempKappaDCompressionOne, tempKappaDCompressionTwo, status->giveDamageCompression() );
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

int
ConcreteDPM2 :: checkForUnAndReloading(double &tempEquivStrain,
                                       double &minEquivStrain,
                                       GaussPoint *gp)
{
    ConcreteDPM2Status *status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    double sigEffective, rhoEffective, thetaEffective;

    //Access old and new strains
    StrainVector oldStrain(status->giveStrainVector(), matMode);
    StrainVector strain(status->giveTempStrainVector(), matMode);

    //Compute the temp equivalent strain
    StressVector tempEffectiveStress(matMode);
    StrainVector tempElasticStrain = strain;
    tempElasticStrain.subtract( status->giveTempPlasticStrain() );
    tempElasticStrain.applyElasticStiffness(tempEffectiveStress, eM, nu);
    computeTrialCoordinates(tempEffectiveStress, sigEffective, rhoEffective, thetaEffective);
    tempEquivStrain = computeEquivalentStrain(sigEffective, rhoEffective, thetaEffective);
    //Get the equivalent strain from the status
    double equivStrain = status->giveEquivStrain();

    //Compute the increment of effective stress
    StressVector effectiveStress(matMode);
    StrainVector elasticStrain = oldStrain;
    elasticStrain.subtract( status->givePlasticStrain() );
    elasticStrain.applyElasticStiffness(effectiveStress, eM, nu);
    StressVector deltaEffectiveStress(matMode);
    deltaEffectiveStress = tempEffectiveStress;
    deltaEffectiveStress.subtract(effectiveStress);

    //Compute equivalent strains for stress state slightly greater than the effective stress and smaller than the temp effective stress
    StressVector intermediateEffectiveStress(matMode);
    intermediateEffectiveStress = effectiveStress;
    //For slightly more than effective stress
    intermediateEffectiveStress.add(0.01 * deltaEffectiveStress);
    computeTrialCoordinates(intermediateEffectiveStress, sigEffective, rhoEffective, thetaEffective);
    double equivStrainPlus = computeEquivalentStrain(sigEffective, rhoEffective, thetaEffective);
    //For slightly less than temp effective stress
    intermediateEffectiveStress = effectiveStress;
    intermediateEffectiveStress.add(0.99 * deltaEffectiveStress);
    computeTrialCoordinates(intermediateEffectiveStress, sigEffective, rhoEffective, thetaEffective);
    double tempEquivStrainMinus = computeEquivalentStrain(sigEffective, rhoEffective, thetaEffective);

    //Check for unloading and reloading in the same step
    int unloadingFlag = 0;
    minEquivStrain = equivStrain;
    double midEquivStrain = 0.;

    unloadingFlag = 0;
    if ( ( equivStrain > equivStrainPlus && tempEquivStrain > tempEquivStrainMinus ) && ( fabs(equivStrainPlus - equivStrain) > yieldTolDamage / 100. && fabs(tempEquivStrainMinus - tempEquivStrain) > yieldTolDamage / 100. ) ) {
        unloadingFlag = 1;
        //Unloading and reloading takes place. Find the minimum equivalent strain by subincrementing the effective stress increment
        for ( double k = 1.0; k <= 100.0; k = k + 1.0 ) {
            intermediateEffectiveStress = effectiveStress;
            intermediateEffectiveStress.add(k / 100. * deltaEffectiveStress);
            computeTrialCoordinates(intermediateEffectiveStress, sigEffective, rhoEffective, thetaEffective);
            midEquivStrain = computeEquivalentStrain(sigEffective, rhoEffective, thetaEffective);

            if ( midEquivStrain <= minEquivStrain ) {
                minEquivStrain = midEquivStrain;
            } else {
                return unloadingFlag;
            }
        }
    }
    return unloadingFlag;
}


double
ConcreteDPM2 :: computeRateFactor(const double alpha,
                                  const double deltaTime,
                                  GaussPoint *gp,
                                  TimeStep *tStep)
{
    if ( strainRateFlag == 0 ) {
        return 1;
    }

    ConcreteDPM2Status *status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    //Access old and new strain
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
    } else { //Compression
        strainRate = ( minStrain - oldRateStrain ) / deltaTime;
        status->letTempRateStrainBe(minStrain);
    }

    //For tension
    double deltaS = 1. / ( 1. + 8. * this->fc / this->fcZero );
    double betaS = exp( ( 6 * deltaS - 2. ) * log(10.) );


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
    } else if ( 30.e-6 < strainRate && strainRate < 1 ) {
        rateFactorTension = pow(strainRateRatioTension, deltaS);
    } else {
        rateFactorTension =  betaS * pow(strainRateRatioTension, 0.333);
    }

    //Compression
    strainRateRatioCompression = strainRate / strainRateZeroCompression;
    if ( strainRate > -30.e-6 ) {
        rateFactorCompression = 1.;
    } else if ( -30.e-6 > strainRate && strainRate > -30 ) {
        rateFactorCompression = pow(strainRateRatioCompression, 1.026 * alphaS);
    } else {
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

    ConcreteDPM2Status *status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    const StrainVector &tempPlasticStrain = status->giveTempPlasticStrain();
    const StrainVector &plasticStrain = status->givePlasticStrain();

    StrainVector deltaPlasticStrain(matMode);
    deltaPlasticStrain = tempPlasticStrain;
    deltaPlasticStrain.subtract(plasticStrain);

    double deltaPlasticStrainNorm = 0;

    //Distinguish pre-peak, peak and post-peak

    double factor = 0.;
    if ( tempKappaD < e0 * ( 1. - yieldTolDamage ) ) {
        deltaPlasticStrainNorm = 0.;
    } else if ( tempKappaD > e0 * ( 1. - yieldTolDamage ) && kappaD < e0  * ( 1. - yieldTolDamage ) ) {
        factor = ( 1. - ( e0 - kappaD ) / ( tempKappaD - kappaD ) );
        deltaPlasticStrain.times(factor);
        deltaPlasticStrainNorm = deltaPlasticStrain.computeNorm();
    } else {
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

    ConcreteDPM2Status *status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    const StrainVector &tempPlasticStrain = status->giveTempPlasticStrain();
    const StrainVector &plasticStrain = status->givePlasticStrain();

    StrainVector deltaPlasticStrain(matMode);
    deltaPlasticStrain.add(tempAlpha, tempPlasticStrain);
    deltaPlasticStrain.add(-tempAlpha, plasticStrain);

    double deltaPlasticStrainNorm = 0;

    //Distinguish pre-peak, peak and post-peak
    double factor = 0.;
    if ( tempKappaD < e0 * ( 1. - yieldTolDamage ) ) {
        deltaPlasticStrainNorm = 0.;
    } else if ( tempKappaD > e0 * ( 1. - yieldTolDamage ) && kappaD < e0 * ( 1. - yieldTolDamage ) ) {
        factor = ( 1. - ( e0 - kappaD ) / ( tempKappaD - kappaD ) );
        deltaPlasticStrain.times(factor);
        deltaPlasticStrainNorm = deltaPlasticStrain.computeNorm();
    } else {
        deltaPlasticStrainNorm = deltaPlasticStrain.computeNorm();
    }

    double tempKappaP = status->giveTempKappaP();
    const double yieldHardTwo = computeHardeningTwo(tempKappaP);
    double extraFactor;
    if ( rho < 1.e-16 ) {
        extraFactor = this->ft * yieldHardTwo * sqrt(2. / 3.) / 1.e-16 / sqrt( 1. + 2. * pow(this->dilationConst, 2.) );
    } else {
        extraFactor = this->ft * yieldHardTwo * sqrt(2. / 3.) / this->rho / sqrt( 1. + 2. * pow(this->dilationConst, 2.) );
    }

    return deltaPlasticStrainNorm * extraFactor;
}


double
ConcreteDPM2 :: computeEquivalentStrain(const double sig,
                                        const double rho,
                                        const double theta)
{
    double rFunction = ( 4. * ( 1. - pow(this->ecc, 2.) ) * pow(cos(theta), 2.) + pow(2. * this->ecc - 1., 2.) ) / ( 2. * ( 1. - pow(this->ecc, 2.) ) * cos(theta) + ( 2. * this->ecc - 1. ) * sqrt(4. * ( 1. - pow(this->ecc, 2.) ) * pow(cos(theta), 2.) + 5. * pow(this->ecc, 2.) - 4. * this->ecc) );

    double pHelp = -this->m * ( rho * rFunction / ( sqrt(6.) * fc ) + sig / fc );

    double qHelp = -3. / 2. * pow(rho, 2.) / pow(this->fc, 2.);

    double help = -0.5 * pHelp + sqrt(pow(pHelp, 2.) / 4. - qHelp);

    double tempEquivStrain = 0.;
    if ( help > 0 ) {
        tempEquivStrain = help * e0;
    }
    return tempEquivStrain;
}


double
ConcreteDPM2 :: computeDamageParamTension(double equivStrain, double kappaOne, double kappaTwo, double le, double omegaOld)
{
    double omega = 0.;

    double help;
    if ( equivStrain > e0 * ( 1. - yieldTolDamage ) ) {
        if ( softeningType == 0 ) { //linear
            omega = ( this->eM * equivStrain * this->wf - this->ft * this->wf + this->ft * kappaOne * le ) /
                    ( this->eM * equivStrain * this->wf - this->ft * le * kappaTwo );
        } else if ( softeningType == 1 ) { //bilinear: Calculate damage parameter for both parts of bilinear curve  and check which fulfils limits.
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
        } else if ( softeningType == 2 ) { //exponential: Iterative solution
            omega = 1.; //Initial guess
            double residual = 0.;
            double dResidualDOmega = 0.;
            int nite = 0;

            do {
                nite++;


                residual  = ( 1 - omega ) * eM * equivStrain - ft *exp(-le * ( omega * kappaTwo + kappaOne ) / wf);
                dResidualDOmega = -eM * equivStrain + ft * le * kappaTwo / wf *exp(-le * ( omega * kappaTwo + kappaOne ) / wf);

                omega -= residual / dResidualDOmega;
                if ( nite > newtonIter ) {
                    OOFEM_ERROR("computeDamageParam: algorithm not converging");
                }
            } while ( fabs(residual / ft) >= 1.e-8 );
        }
    } else {
        omega = 0.;
    }


    if ( omega > 1. ) {
        omega = 1.;
    }

    if ( omega < 0. || omega < omegaOld ) {
        omega = omegaOld;
    }


    return omega;
}

double
ConcreteDPM2 :: computeDamageParamCompression(double equivStrain, double kappaOne, double kappaTwo, double omegaOld)
{
    if ( isotropicFlag == 1 ) {
        return 0.;
    }


    double omega = 1.;
    int nite = 0;
    double residual = 0.;
    double dResidualDOmega = 0.;
    if ( equivStrain > e0 * ( 1. - yieldTolDamage ) ) {
        do {
            nite++;

            residual = ( 1. - omega ) * this->eM * equivStrain - this->ft *exp(-( kappaOne + omega * kappaTwo ) / this->efCompression);
            dResidualDOmega = -this->eM * equivStrain + this->ft * kappaTwo / this->efCompression *exp(-( kappaOne + omega * kappaTwo ) / this->efCompression);

            omega -= residual / dResidualDOmega;
            if ( nite > newtonIter ) {
                OOFEM_ERROR("computeDamageParam: algorithm not converging");
            }
        } while ( fabs(residual / ft) >= 1.e-8 );
    } else {
        omega = 0.;
    }

    if ( omega > 1. ) {
        omega = 1.;
    }
    if ( omega < omegaOld || omega < 0. ) {
        omega = omegaOld;
    }

    return omega;
}



void
ConcreteDPM2 :: initDamaged(double kappaD,
                            const StrainVector &strain,
                            GaussPoint *gp)
{
    if ( kappaD <= 0. ) {
        return;
    }

    int i, indx = 1;
    double le;
    FloatArray principalStrains, crackPlaneNormal(3);
    FloatMatrix principalDir(3, 3);
    ConcreteDPM2Status *status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    if ( helem > 0. ) {
        status->setLe(helem);
    } else if ( matMode == _1dMat ) {
        le = gp->giveElement()->computeLength();
        status->setLe(le);
    } else if ( ( status->giveDamageTension() == 0. ) && ( status->giveDamageCompression() == 0. ) ) {
        strain.computePrincipalValDir(principalStrains, principalDir);
        // find index of max positive principal strain
        for ( i = 2; i <= 3; i++ ) {
            if ( principalStrains.at(i) > principalStrains.at(indx) ) {
                indx = i;
            }
        }

        for ( i = 1; i <= 3; i++ ) {
            crackPlaneNormal.at(i) = principalDir.at(i, indx);
        }

        // evaluate the projected element size
        le = gp->giveElement()->giveCharacteristicLength(crackPlaneNormal);
        if ( le == 0. ) {
            le = gp->giveElement()->computeMeanSize();
        }

        // store le in the corresponding status
        status->setLe(le);
    } else if ( status->giveLe() == 0. ) {
        // this happens if the status is initialized from a file
        // with nonzero damage
        // le determined as square root of element area or cube root of el. volume
        le = gp->giveElement()->computeMeanSize();
        status->setLe(le);
    }
}


double
ConcreteDPM2 :: computeDuctilityMeasureDamage(const StrainVector &strain, GaussPoint *gp)
{
    double ductilityMeasure;

    //Angle in uniaxial compression is atan(1./sqrt(6.))=0.387597
    double alphaZero = 0.40824829;

    double Rs = 0;
    if ( this->sig < 0. ) {
        if ( this->rho > 1.e-16 ) {
            Rs = -this->sig / ( alphaZero * this->rho );
        } else { //Need to set a dummy valye
            Rs = -this->sig * 1.e16 / alphaZero;
        }
    } else {
        Rs = 0;
    }

    ductilityMeasure = 1. + ( ASoft - 1. ) * Rs;

    return ductilityMeasure;
}


void
ConcreteDPM2 :: performPlasticityReturn(GaussPoint *gp,
                                        StrainVector &strain)
{
    ConcreteDPM2Status *status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );
    if ( matMode == _Unknown ) {
        matMode = gp->giveMaterialMode();
    }

    //get temp plastic strain and tempKappa
    StrainVector tempPlasticStrain = status->giveTempPlasticStrain();
    double tempKappaP = status->giveTempKappaP();

    // compute elastic strains and trial stress
    StressVector effectiveStress(matMode);
    StrainVector elasticStrain = strain;
    elasticStrain.subtract(tempPlasticStrain);
    elasticStrain.applyElasticStiffness(effectiveStress, eM, nu);

    //Compute trial coordinates
    computeTrialCoordinates(effectiveStress, this->sig, this->rho, this->thetaTrial);

    double yieldValue;

    StrainVector convergedStrain(matMode);

    //What happens here to thermal strains that are subtracted first?
    StrainVector oldStrain(status->giveStrainVector(), matMode);
    StrainVector tempStrain(matMode);
    StrainVector deltaStrain(matMode);

    // introduce a strange subincrementation flag
    int subIncrementFlag = 0;

    double apexStress = 0.;
    int subincrementcounter = 0;
    //Attempt to implement subincrementation
    // initialize variables
    subIncrementFlag = 0;
    convergedStrain = oldStrain;
    tempStrain = strain;
    deltaStrain.beDifferenceOf(strain, oldStrain);
    //To get into the loop
    returnResult = RR_NotConverged;
    while ( returnResult == RR_NotConverged || subIncrementFlag == 1 ) {
        elasticStrain = tempStrain;
        elasticStrain.subtract(tempPlasticStrain);
        elasticStrain.applyElasticStiffness(effectiveStress, eM, nu);
        double theta;
        computeTrialCoordinates(effectiveStress, this->sig, this->rho, theta);
        yieldValue = computeYieldValue(this->sig, this->rho, this->thetaTrial, tempKappaP);

        apexStress = 0.;

        if ( yieldValue > 0. ) {
            checkForVertexCase(apexStress, sig, tempKappaP);
            if ( returnType == RT_Tension || returnType == RT_Compression ) {
                tempKappaP = performVertexReturn(effectiveStress, apexStress, tempKappaP, gp);
            }
            if ( returnType == RT_Regular ) {
                tempKappaP = performRegularReturn(effectiveStress, gp);
            }
        } else {
            returnResult = RR_Converged;
            StrainVector plasticStrain(status->givePlasticStrain(), matMode);
            tempPlasticStrain = plasticStrain;
            status->letTempPlasticStrainBe(tempPlasticStrain);
            status->letTempKappaPBe(tempKappaP);
            break;
        }

        if ( returnResult == RR_NotConverged ) {
            subincrementcounter++;
            if ( subincrementcounter > 10 ) {
                OOFEM_LOG_INFO( "Unstable element %d \n", gp->giveElement()->giveGlobalNumber() );
                OOFEM_LOG_INFO( "Old strain vector %g %g %g %g %g %g  \n", oldStrain.at(1), oldStrain.at(2), oldStrain.at(3), oldStrain.at(4), oldStrain.at(5), oldStrain.at(6) );
                StrainVector help(status->giveTempPlasticStrain(), matMode);
                StressVector help1(matMode);
                double sig1, rho1, theta1;
                oldStrain.subtract(help);
                oldStrain.applyElasticStiffness(help1, eM, nu);
                OOFEM_LOG_INFO( "Old plastic strain vector %g %g %g %g %g %g  \n", help.at(1), help.at(2), help.at(3), help.at(4), help.at(5), help.at(6) );
                OOFEM_LOG_INFO( "New strain vector %g %g %g %g %g %g  \n", strain.at(1), strain.at(2), strain.at(3), strain.at(4), strain.at(5), strain.at(6) );
                computeTrialCoordinates(effectiveStress, this->sig, this->rho, theta);
                computeTrialCoordinates(help1, sig1, rho1, theta1);
                yieldValue = computeYieldValue(this->sig, this->rho, this->thetaTrial, tempKappaP);
                OOFEM_LOG_INFO("OLD Sig %g rho %g theta %g  \n", sig1, rho1, theta1);
                OOFEM_LOG_INFO("NEW Sig %g rho %g theta %g  \n", sig, rho, theta);
                if ( returnType == RT_Tension || returnType == RT_Compression ) {
                    OOFEM_LOG_INFO("Vertex case apexstress %g\n", apexStress);
                } else {
                    OOFEM_LOG_INFO("Regular case %g \n", 15.18);
                }
                OOFEM_LOG_INFO("KappaP old %g new %g yieldfun %g\n", status->giveTempKappaP(), tempKappaP, yieldValue);
                OOFEM_ERROR("ConcreteDPM2 :: performPlasticityReturn - Could not reach convergence with small deltaStrain, giving up.");
            } else if ( subincrementcounter > 9 && tempKappaP < 1. ) {
                tempKappaP = 1.;
                status->letTempKappaPBe(tempKappaP);
            }

            subIncrementFlag = 1;
            deltaStrain.times(0.5);
            tempStrain = convergedStrain;
            tempStrain.add(deltaStrain);
        } else if ( returnResult == RR_Converged && subIncrementFlag == 0 ) {
            effectiveStress.applyElasticCompliance(elasticStrain, eM, nu);
            tempPlasticStrain = strain;
            tempPlasticStrain.subtract(elasticStrain);
            status->letTempPlasticStrainBe(tempPlasticStrain);
            status->letTempKappaPBe(tempKappaP);
        } else if ( returnResult == RR_Converged && subIncrementFlag == 1 ) {
            OOFEM_LOG_INFO("Subincrementation %d required\n", subincrementcounter);
            subincrementcounter = 0;
            effectiveStress.applyElasticCompliance(elasticStrain, eM, nu);
            tempPlasticStrain = tempStrain;
            tempPlasticStrain.subtract(elasticStrain);
            status->letTempPlasticStrainBe(tempPlasticStrain);
            status->letTempKappaPBe(tempKappaP);

            subIncrementFlag = 0;
            returnResult = RR_NotConverged;
            convergedStrain = tempStrain;
            deltaStrain.beDifferenceOf(strain, convergedStrain);
            tempStrain = strain;
        }
    }
}


bool
ConcreteDPM2 :: checkForVertexCase(double &answer,
                                   const double sig,
                                   const double tempKappa)
{
    if ( matMode == _1dMat ) {
        returnType = RT_Regular;
        return false;
    }

    answer = 0.;
    if ( sig > 0. ) {
        returnType = RT_Tension;
    } else if ( sig < 0. &&  tempKappa < 1. ) {
        returnType = RT_Compression;
    } else {
        returnType = RT_Regular;
    }
    return false;
}



double
ConcreteDPM2 :: performVertexReturn(StressVector &effectiveStress,
                                    double apexStress, double tempKappaP,
                                    GaussPoint *gp)
{
    //    ConcreteDPM2Status *status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );
    StressVector deviatoricStressTrial(matMode);
    double sigTrial;
    double yieldValue = 0.;
    double yieldValueMid = 0.;
    double sig2 = 0.;
    double dSig;
    double sigMid;
    double sigAnswer;
    double ratioPotential;

    effectiveStress.computeDeviatoricVolumetricSplit(deviatoricStressTrial, sigTrial);
    const double rhoTrial = deviatoricStressTrial.computeSecondCoordinate();

    //    double tempKappaP = status->giveTempKappaP();
    const double kappaInitial = tempKappaP;

    sig2 = apexStress;

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
        returnResult = RR_NotConverged;
        return kappaInitial;
    }

    if ( yieldValue < 0.0 ) {
        dSig = sig2 - sigTrial;
        sigAnswer = sig2;
    } else {
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
            // //Put tempKappaP in the status
            // status->letTempKappaPBe(tempKappaP);

            ratioPotential =
                computeRatioPotential(sigAnswer, tempKappaP);


            double ratioTrial = rhoTrial / ( sigTrial - sigAnswer );

            if ( ( ( ( ratioPotential >= ratioTrial ) && returnType == RT_Tension ) ) ||
                 ( ( ratioPotential <= ratioTrial ) && returnType == RT_Compression ) ) {
                for ( int i = 0; i < 3; i++ ) {
                    effectiveStress(i) = sigAnswer;
                }

                for ( int i = 3; i < effectiveStress.giveSize(); i++ ) {
                    effectiveStress(i) = 0.;
                }
                returnResult = RR_Converged;
                return tempKappaP;
            } else {
                returnType = RT_Regular;
                returnResult = RR_NotConverged;
                return kappaInitial;
            }
        }
    }

    for ( int i = 0; i < 3; i++ ) {
        effectiveStress(i) = sigAnswer;
    }

    for ( int i = 3; i < effectiveStress.giveSize(); i++ ) {
        effectiveStress(i) = 0.;
    }
    returnResult = RR_Converged;
    //// Warning solution NOT CONVERGED!!!!
    return tempKappaP;
}


double
ConcreteDPM2 :: computeTempKappa(const double kappaInitial,
                                 const double sigTrial,
                                 const double rhoTrial,
                                 const double sig)
{
    //This function is called, if stress state is in vertex case
    double equivalentDeltaPlasticStrain;
    rho = 0.;
    equivalentDeltaPlasticStrain = sqrt( 1. / 9. * pow( ( sigTrial - sig ) / ( kM ), 2. ) +
                                         pow(rhoTrial / ( 2. * gM ), 2.) );

    double thetaVertex = M_PI / 3.;
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
    } else {
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

    const double dgdsig = 4. * ( 1. - yieldHardOne ) / fc * Al * Bl + pow(yieldHardOne, 2.) * mQ / fc;

    const double dgdrho = Al / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. ) +

      m *pow(yieldHardOne, 2.) / ( sqrt(6.) * fc );

    return dgdrho / dgdsig * 3. * ( 1. - 2. * nu ) / ( 1. + nu );
}


double
ConcreteDPM2 :: performRegularReturn(StressVector &effectiveStress,
                                     GaussPoint *gp)
{
    ConcreteDPM2Status *status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    MaterialMode matMode = gp->giveMaterialMode();
    StressVector trialStress(matMode), deviatoricTrialStress(matMode);

    double yieldValue;

    FloatArray residuals(4), residualsNorm(4), dGDInv(2);
    FloatMatrix jacobian(4, 4), inverseOfJacobian(4, 4);
    double dKappaDDeltaLambda;
    FloatArray deltaIncrement(3);
    FloatArray deltaIncrementTemp(3);

    double deltaLambda = 0.;
    double normOfResiduals = 0.;
    int iterationCount = 0;

    FloatArray unknownsTrial;
    FloatArray residualsTrial;
    FloatArray unknowns(4);

    if ( matMode == _1dMat ) { //Resize Matrices and Arrays for 1D
        residuals.resize(3);
        residualsNorm.resize(3);
        jacobian.resize(3, 3);
        inverseOfJacobian.resize(3, 3);
        unknowns.resize(3);
        deltaIncrement.resize(2); // I think this should be 3 and for the 3D case 4 but it does not matter
        deltaIncrementTemp.resize(2);
        dGDInv.resize(1);
    }

    deltaIncrement.zero();
    deltaIncrementTemp.zero();

    //Define stressVariables
    double trialSig, trialRho;

    trialStress = effectiveStress;

    //compute the principal directions of the stress
    FloatArray helpStress;
    FloatMatrix stressPrincipalDir;

    //compute invariants from stress state
    if ( matMode != _1dMat ) {
        trialStress.computePrincipalValDir(helpStress, stressPrincipalDir);
        trialStress.computeDeviatoricVolumetricSplit(deviatoricTrialStress, trialSig);
        trialRho = deviatoricTrialStress.computeSecondCoordinate();
    } else {  //1d case
        double angle; // this variable is used only as an input to the function computeTrialCoordinates and is not important and it is already calculated
        computeTrialCoordinates(trialStress, trialSig, trialRho, angle);
    }

    sig = trialSig;
    rho = trialRho;


    // Do the same for kappa
    double kappaP = status->giveTempKappaP();
    double tempKappaP = kappaP;

    //initialise unknowns
    if ( matMode != _1dMat ) {
        unknowns.at(1) = trialSig;
        unknowns.at(2) = trialRho;
        unknowns.at(3) = tempKappaP;
        unknowns.at(4) = 0.;
    } else {  //1D case
        unknowns.at(1) = trialSig * 3.; // It is calculated as the volumetric stress in this case sigma/3
        unknowns.at(2) = tempKappaP;
        unknowns.at(3) = 0.;
    }

    // Look at the magnitudes of the residuals. You have to scale the yieldValue down.
    yieldValue = computeYieldValue(sig, rho, thetaTrial, tempKappaP);

    //initiate residuals
    residuals.zero();
    residuals.at( residuals.giveSize() ) = yieldValue; //store in the last element of the array

    normOfResiduals  = 1.; //just to get into the loop

    /* N.R. iteration for finding the correct plastic return which is found when the norm of the residuals are equal to zero*/

    while ( normOfResiduals > yieldTol   ) {
        iterationCount++;
        if ( iterationCount == newtonIter ) {
            returnResult = RR_NotConverged;
            return status->giveTempKappaP();
        }

        if ( matMode != _1dMat ) {
            //Normalize residuals. Think about it more.
            residualsNorm.at(1) = residuals.at(1) / this->kM;
            residualsNorm.at(2) = residuals.at(2) / ( 2. * this->gM );
            residualsNorm.at(3) = residuals.at(3);
            residualsNorm.at(4) = residuals.at(4);
        } else {  //1D case
            residualsNorm.at(1) = residuals.at(1) / this->eM;
            residualsNorm.at(2) = residuals.at(2);
            residualsNorm.at(3) = residuals.at(3);
        }

        normOfResiduals = residualsNorm.computeNorm();

        if ( isnan(normOfResiduals) ) {
            returnResult = RR_NotConverged;
            return status->giveTempKappaP();
        }

        if ( normOfResiduals > yieldTol  ) {
            // Test to run newton iteration using inverse of Jacobian
            if ( matMode != _1dMat ) {
                computeJacobian(jacobian, sig, rho, tempKappaP, deltaLambda, gp);

                if ( computeInverseOfJacobian(inverseOfJacobian, jacobian) ) {
                    returnResult = RR_NotConverged;
                    return status->giveTempKappaP();
                }

                deltaIncrement.beProductOf(inverseOfJacobian, residuals);
                deltaIncrement.negated();

                //compute trial values
                unknownsTrial = unknowns;
                residualsTrial = residuals;

                //compute Unknowns
                for ( int i = 0; i < unknowns.giveSize(); i++ ) {
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
            } else {
                compute1dJacobian(jacobian, 3. * sig, tempKappaP, deltaLambda, gp);

                if ( computeInverseOfJacobian(inverseOfJacobian, jacobian) ) {
                    returnResult = RR_NotConverged;
                    return status->giveTempKappaP();
                }

                deltaIncrement.beProductOf(inverseOfJacobian, residuals);
                deltaIncrement.negated();

                //compute trial values
                unknownsTrial = unknowns;
                residualsTrial = residuals;

                //compute Unknowns
                for ( int i = 1; i <= unknownsTrial.giveSize(); i++ ) {
                    unknowns.at(i) = unknownsTrial.at(i) + deltaIncrement.at(i);
                }

                if ( unknowns.at(3) <= 0. ) { //Keep deltaLambda greater equal than zero!
                    unknowns.at(3) = 0.;
                }

                if ( unknowns.at(2) - kappaP <= 0. ) { //Keep deltaKappa greater equal than zero!
                    unknowns.at(2) = kappaP;
                }

                //compute residuals
                sig = unknowns.at(1) / 3.;
                rho = unknowns.at(1) * sqrt(2. / 3.); //for the 1d case

                tempKappaP = unknowns.at(2);
                deltaLambda = unknowns.at(3);

                /* Compute the mVector holding the derivatives of the g function and the hardening function*/
                computeDGDInv1d(dGDInv.at(1), unknowns.at(1), tempKappaP);
                dKappaDDeltaLambda = computeDKappaDDeltaLambda(unknowns.at(1), 1., tempKappaP);

                residuals.at(1) = 3. * ( sig - trialSig ) + this->eM *deltaLambda *dGDInv.at(1);

                residuals.at(2) = -tempKappaP + kappaP + deltaLambda * dKappaDDeltaLambda;
                residuals.at(3) = computeYieldValue(sig, rho, thetaTrial, tempKappaP);
            }
        }
    }

    returnResult = RR_Converged;
    if ( matMode == _1dMat ) {
        effectiveStress.at(1) = sig * 3;
    } else {
        StressVector stressPrincipal(_3dMat);
        StressVector stressTemp(_3dMat);
        stressPrincipal.zero();

        stressPrincipal(0) = sig + sqrt(2. / 3.) * rho * cos(thetaTrial);
        stressPrincipal(1) = sig + sqrt(2. / 3.) * rho * cos(thetaTrial - 2. * M_PI / 3.);
        stressPrincipal(2) = sig + sqrt(2. / 3.) * rho * cos(thetaTrial + 2. * M_PI / 3.);
        if ( matMode == _PlaneStrain ) {
            transformStressVectorTo(stressTemp, stressPrincipalDir, stressPrincipal, 1);
            effectiveStress(0) = stressTemp(0);
            effectiveStress(1) = stressTemp(1);
            effectiveStress(2) = stressTemp(2);
            effectiveStress(3) = stressTemp(5);
        } else {
            transformStressVectorTo(stressTemp, stressPrincipalDir, stressPrincipal, 1);
            effectiveStress = stressTemp;
        }
    }

    return tempKappaP;
}

void
ConcreteDPM2 :: compute1dJacobian(FloatMatrix &answer,
                                  const double totalsigma,
                                  const double kappa,
                                  const double deltaLambda,
                                  GaussPoint *gp)
{
    //Variables
    double dFDInv;
    computeDFDInv1d(dFDInv, totalsigma, kappa);

    double dGDInv;
    computeDGDInv1d(dGDInv, totalsigma, kappa);

    double dDGDDInv;
    computeDDGDDInv1d(dDGDDInv, totalsigma, kappa);

    double dKappaDDeltaLambda;
    dKappaDDeltaLambda = computeDKappaDDeltaLambda(totalsigma, 1, kappa);

    double dFDKappa;
    dFDKappa = computeDFDKappa(totalsigma, 1, kappa);

    double dDGDInvDKappa;
    computeDDGDInvDKappa1d(dDGDInvDKappa, totalsigma, kappa);

    double dDKappaDDeltaLambdaDKappa;
    dDKappaDDeltaLambdaDKappa = computeDDKappaDDeltaLambdaDKappa(totalsigma, 1, kappa);

    double dDKappaDDeltaLambdaDInv;
    computeDDKappaDDeltaLambdaDInv1d(dDKappaDDeltaLambdaDInv, totalsigma, kappa);

    /* Compute matrix*/
    /* 1st row */
    answer.at(1, 1) = 1. + this->eM * deltaLambda * dDGDDInv;
    answer.at(1, 2) = this->eM * deltaLambda * dDGDInvDKappa;
    answer.at(1, 3) = this->eM * dGDInv;
    /* 2nd row */
    answer.at(2, 1) = deltaLambda * dDKappaDDeltaLambdaDInv;
    answer.at(2, 2) = deltaLambda * dDKappaDDeltaLambdaDKappa - 1.;
    answer.at(2, 3) = dKappaDDeltaLambda;
    /* 3rd row */
    answer.at(3, 1) = dFDInv;
    answer.at(3, 2) = dFDKappa;
    answer.at(3, 3) = 0.;
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
           pow(yieldHardOne, 2.) * yieldHardTwo * m * ( sig / fc + rho * rFunction / ( sqrt(6.) * fc ) ) -
           pow(yieldHardOne, 2.) * pow(yieldHardTwo, 2.);
}


double
ConcreteDPM2 :: computeDFDKappa(const double sig,
                                const double rho,
                                const double tempKappa)
{
    double dFDKappa;
    const double theta = thetaTrial;
    //compute yieldHard and yieldSoft
    const double yieldHardOne = computeHardeningOne(tempKappa);
    const double yieldHardTwo = computeHardeningTwo(tempKappa);
    // compute the derivative of the hardening and softening laws
    const double dYieldHardOneDKappa = computeHardeningOnePrime(tempKappa);
    const double dYieldHardTwoDKappa = computeHardeningTwoPrime(tempKappa);
    //compute elliptic function r
    const double rFunction =
        ( 4. * ( 1. - pow(ecc, 2) ) * pow(cos(theta), 2) + pow( ( 2. * ecc - 1. ), 2 ) ) /
        ( 2 * ( 1. - pow(ecc, 2) ) * cos(theta) + ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - pow(ecc, 2) ) * pow(cos(theta), 2) + 5. * pow(ecc, 2) - 4. * ecc) );


    if ( matMode != _1dMat ) {
        //compute help functions Al, Bl
        const double Al = ( 1. - yieldHardOne ) * pow( ( sig / fc + rho / ( sqrt(6.) * fc ) ), 2. ) + sqrt(3. / 2.) * rho / fc;

        const double Bl = sig / fc + rho / ( fc * sqrt(6.) );

        const double dFDYieldHardOne = -2. *Al *pow(Bl, 2.)
                                       + 2. * yieldHardOne * yieldHardTwo * m * ( sig / fc + rho * rFunction / ( sqrt(6.) * fc ) ) - 2. *yieldHardOne *pow(yieldHardTwo, 2.);

        const double dFDYieldHardTwo = pow(yieldHardOne, 2.) * m * ( sig / fc + rho * rFunction / ( sqrt(6.) * fc ) ) - 2. *yieldHardTwo *pow(yieldHardOne, 2.);

        // compute dFDKappa
        dFDKappa =  dFDYieldHardOne * dYieldHardOneDKappa +
                   dFDYieldHardTwo * dYieldHardTwoDKappa;
    } else {  //1d case
        dFDKappa = -2 * pow(2 * sig / 3 / fc, 2) * ( sig / fc + pow(2 / 3 * sig / fc, 2) * ( 1 - yieldHardOne ) ) * dYieldHardOneDKappa +
                   ( 1 + rFunction ) * m * sig / 3 / fc * ( dYieldHardOneDKappa * 2. * yieldHardOne * yieldHardTwo + dYieldHardTwoDKappa * yieldHardOne ) -
                   2 * ( yieldHardOne * pow(yieldHardTwo, 2) * dYieldHardOneDKappa + yieldHardTwo * pow(yieldHardOne, 2) * dYieldHardTwoDKappa );
    }

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
ConcreteDPM2 :: computeDFDInv1d(double &answer,
                                const double sigma,
                                const double tempKappa) const
{
    const double theta = thetaTrial;
    //compute yieldHard
    const double yieldHardOne = computeHardeningOne(tempKappa);
    const double yieldHardTwo = computeHardeningTwo(tempKappa);
    //compute elliptic function r
    const double rFunction =  ( 4. * ( 1. - pow(ecc, 2) ) * pow(cos(theta), 2) + pow( ( 2. * ecc - 1. ), 2 ) ) /
                             ( 2. * ( 1. - pow(ecc, 2) ) * cos(theta) + ( 2. * ecc - 1. ) * sqrt(4. * ( 1. - pow(ecc, 2) ) * pow(cos(theta), 2) + 5. * pow(ecc, 2) - 4. * ecc) );
    answer = 2 * ( 1 / fc + 8 * sigma / pow(3 * fc, 2) * ( 1 - yieldHardOne ) ) * ( sigma / fc + pow(2 * sigma / 3 / fc, 2) * ( 1 - yieldHardOne ) ) + ( 1 + rFunction ) * m / ( 3 * fc ) * pow(yieldHardOne, 2) * yieldHardTwo;
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
    const double dfdsig = 4. * ( 1. - yieldHardOne ) / fc * AL * BL + yieldHardTwo *pow(yieldHardOne, 2.) * m / fc;

    //compute dfdrho
    const double dfdrho = AL / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * BL + 6. ) + rFunction *m *yieldHardTwo *pow(yieldHardOne, 2.) / ( sqrt(6.) * fc );

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
    double ductilityMeasure;

    if ( matMode != _1dMat ) {
        FloatArray dGDInv(2);
        computeDGDInv(dGDInv, sig, rho, tempKappa);

        equivalentDGDStress = sqrt( 1. / 3. * pow(dGDInv(0), 2.) +
                                    pow(dGDInv(1), 2.) );

        ductilityMeasure = computeDuctilityMeasure(sig, rho, this->thetaTrial);
    } else {  //1D case
        computeDGDInv1d(equivalentDGDStress, sig, tempKappa);
        equivalentDGDStress = fabs(equivalentDGDStress); // In 1D is the absolute value
        ductilityMeasure = computeDuctilityMeasure(sig / 3., sqrt(2. / 3.) * sig, this->thetaTrial);
    }

    double dKappaDDeltaLambda = equivalentDGDStress / ductilityMeasure;
    return dKappaDDeltaLambda;
}

void
ConcreteDPM2 :: computeDDKappaDDeltaLambdaDInv1d(double &answer,
                                                 const double sigma,
                                                 const double tempKappa)
{
    //Variables
    double dGDInv;
    double dDGDDInv;

    //Compute first and second derivative of plastic potential
    computeDGDInv1d(dGDInv, sigma, tempKappa);
    computeDDGDDInv1d(dDGDDInv, sigma, tempKappa);


    //computeDuctilityMeasure
    double ductilityMeasure = computeDuctilityMeasure(sigma / 3., sigma * sqrt(2. / 3.), this->thetaTrial);

    // compute the derivative of
    double dDuctilityMeasureDInv;
    computeDDuctilityMeasureDInv1d(dDuctilityMeasureDInv, sigma, tempKappa);
    if ( dGDInv < 0 ) {
        dDGDDInv = -dDGDDInv;
        dGDInv = -dGDInv;
    }

    answer = dDGDDInv / ductilityMeasure - dGDInv * dDuctilityMeasureDInv / pow(ductilityMeasure, 2);
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
    FloatArray dEquivalentDGDStressDInv(2);

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
    double dEquivalentDGDStressDKappa;
    double dDKappaDDeltaLambdaDKappa;
    double ductilityMeasure;
    if ( matMode != _1dMat ) {
        FloatArray dGDInv(2);
        FloatArray dDGDInvDKappa(2);

        //Compute first and second derivative of plastic potential
        computeDGDInv(dGDInv, sig, rho, tempKappa);
        computeDDGDInvDKappa(dDGDInvDKappa, sig, rho, tempKappa);

        equivalentDGDStress = sqrt( 1. / 3. * pow(dGDInv(0), 2.) +
                                    pow(dGDInv(1), 2.) );

        ductilityMeasure = computeDuctilityMeasure(sig, rho, this->thetaTrial);  //computeDuctilityMeasure
        //Compute dEquivalentDGDStressDKappa
        dEquivalentDGDStressDKappa =
            ( 2. / 3. * dGDInv(0) * dDGDInvDKappa(0) + 2. * dGDInv(1) * dDGDInvDKappa(1) ) / ( 2. * equivalentDGDStress );

        // compute the derivative of
        double dDuctilityMeasureDKappa = 0.;

        ///@todo Is this right? This is *NEVER* used.
        dDKappaDDeltaLambdaDKappa =
            ( dEquivalentDGDStressDKappa * ductilityMeasure -
              equivalentDGDStress * dDuctilityMeasureDKappa ) / pow(ductilityMeasure, 2.);
    } else {
        computeDGDInv1d(equivalentDGDStress, sig, tempKappa);
        computeDDGDInvDKappa1d(dEquivalentDGDStressDKappa, sig, tempKappa);
        if ( equivalentDGDStress < 0 ) {
            dEquivalentDGDStressDKappa = ( -1 ) * dEquivalentDGDStressDKappa;                        //We are differentiating the absolute value of the first derivative of G with respect to stress
        }

        ductilityMeasure = computeDuctilityMeasure(sig / 3., sig * sqrt(2. / 3.), this->thetaTrial);
    }

    dDKappaDDeltaLambdaDKappa = dEquivalentDGDStressDKappa / ductilityMeasure;
    return dDKappaDDeltaLambdaDKappa;
}

void
ConcreteDPM2 :: computeDDuctilityMeasureDInv1d(double &answer,
                                               const double sigma,
                                               const double tempKappa)
{
    double thetaConst = pow(2. * cos(thetaTrial), 2.);
    double x =  -( sigma + fc ) / ( 3 * fc ); //R hardening variable
    double dXDSigma = -1. / ( 3. * fc );
    if ( x < 0. ) {
        /* Introduce exponential help function which results in a
         * smooth transition. */
        double EHard = BHard - DHard;
        double FHard = ( BHard - DHard ) * CHard / ( AHard - BHard );

        double dDuctilityMeasureDX = EHard / FHard *exp(x / FHard) / thetaConst;
        answer = dDuctilityMeasureDX * dXDSigma;
    } else {
        double dDuctilityMeasureDX = ( AHard - BHard ) /  CHard  / thetaConst *exp(-x /  CHard);
        answer = dDuctilityMeasureDX * dXDSigma;
    }
}

void
ConcreteDPM2 :: computeDDuctilityMeasureDInv(FloatArray &answer,
                                             const double sig,
                                             const double rho,
                                             const double tempKappa)
{
    double thetaConst = pow(2. * cos(thetaTrial), 2.);
    double x = ( -( sig + fc / 3. ) ) / fc;

    if ( x < 0. ) {
        double dXDSig = -1. / fc;
        /* Introduce exponential help function which results in a
         * smooth transition. */
        double EHard = BHard - DHard;
        double FHard = ( BHard - DHard ) * CHard / ( AHard - BHard );

        double dDuctilityMeasureDX = EHard / FHard *exp(x / FHard) / thetaConst;
        answer(0) = dDuctilityMeasureDX * dXDSig;
        answer(1) = 0.;
    } else {
        double dXDSig = -1. / fc;
        double dDuctilityMeasureDX = -( BHard - AHard ) / ( CHard ) / thetaConst *exp( -x / ( CHard ) );
        answer(0) = dDuctilityMeasureDX * dXDSig;
        answer(1) = 0.;
    }
}


void
ConcreteDPM2 :: computeDGDInv1d(double &answer,
                                const double sigma,
                                const double tempKappa)
{
    //compute yieldHard and yieldSoft
    const double yieldHardOne = computeHardeningOne(tempKappa);
    const double yieldHardTwo = computeHardeningTwo(tempKappa);
    double AGParam = this->ft * yieldHardTwo * 3 / this->fc + m / 2;
    double BGParam =  yieldHardTwo / 3. * ( 1. + this->ft / this->fc ) / ( log(AGParam) + log(this->dilationConst + 1.) - log(2 * this->dilationConst - 1.) - log(3. * yieldHardTwo + this->m / 2) );
    const double R = ( sigma - yieldHardTwo * ft ) / ( 3 * fc * BGParam );
    const double mQ = AGParam * exp(R) / 3;
    answer = 2 * ( 1 / fc + 8 * sigma / pow(3 * fc, 2) * ( 1 - yieldHardOne ) ) * ( sigma / fc + pow(2 * sigma / ( 3 * fc ), 2) * ( 1 - yieldHardOne ) ) + pow(yieldHardOne, 2) / fc * ( m / 3 + mQ );
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

    const double dgdsig = 4. * ( 1. - yieldHardOne ) / fc * Al * Bl + pow(yieldHardOne, 2.) * mQ / fc;

    const double dgdrho = Al / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. ) +
                          m *pow(yieldHardOne, 2.) / ( sqrt(6.) * fc );

    answer(0) = dgdsig;
    answer(1) = dgdrho;
}

void
ConcreteDPM2 :: computeDDGDInvDKappa1d(double &answer,
                                       const double sigma,
                                       const double tempKappa)
{
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

    const double R = ( sigma - ft * yieldHardTwo ) / ( 3 * fc * BGParam );
    const double mQ = AGParam * exp(R) / 3;
    //Compute the derivative of mQ with respect to kappa

    //Derivative of AGParam
    double dAGParamDKappa = dYieldHardTwoDKappa * 3. * this->ft / this->fc;

    //Derivative of BGParam
    double BGParamTop = yieldHardTwo / 3. * ( 1. + this->ft / this->fc );
    double BGParamBottom = ( log(AGParam) + log(this->dilationConst + 1.) - log(2 * this->dilationConst - 1.) - log(3. * yieldHardTwo + this->m / 2) );
    double dBGParamTopDKappa1 = dYieldHardTwoDKappa * ( 1 + ft / fc ) / 3;
    double dBGParamBottomDKappa1 = BGParamBottom;
    double dBGParamTopDKappa2 = BGParamTop * ( dAGParamDKappa / AGParam - 3 * dYieldHardTwoDKappa / ( m / 2 + 3 * yieldHardTwo ) );
    double dBGParamBottomDKappa2 = pow(BGParamBottom, 2);

    double dBGParamDKappa = dBGParamTopDKappa1 / dBGParamBottomDKappa1 - dBGParamTopDKappa2 / dBGParamBottomDKappa2;
    double dMQDKappa = 1 / 3 * exp(R) * ( dAGParamDKappa - AGParam * ( ( sigma - ft * yieldHardTwo ) * dBGParamDKappa / ( 3 * fc * pow(BGParam, 2) ) + ft * dYieldHardTwoDKappa / 3 / fc / BGParam ) );

    answer = -8 / 9 * pow(sigma / fc, 2) * ( 1 / fc + 8 * sigma / pow(3 * fc, 2) * ( 1 - yieldHardOne ) ) * dYieldHardOneDKappa - sigma *pow(4 / 3 / fc, 2) * ( sigma / fc + pow(2 / 3 * sigma / fc, 2) * ( 1 - yieldHardOne ) ) * dYieldHardOneDKappa + 2 * dYieldHardOneDKappa * yieldHardOne / fc * ( this->m / 3 + mQ ) + yieldHardOne / fc * dMQDKappa;
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
        dYieldHardOneDKappa * 2 * yieldHardOne * mQ / fc + yieldHardOne * dMQDKappa / fc;

    const double dDGDRhoDKappa =
        ( dAlDYieldHard / ( sqrt(6.) * fc ) * ( 4. * ( 1. - yieldHardOne ) * Bl + 6. ) -
          4. * Al / ( sqrt(6.) * fc ) * Bl + m / ( sqrt(6.) * fc ) ) * 2 * yieldHardOne * dYieldHardOneDKappa;

    answer(0) = dDGDSigDKappa;
    answer(1) = dDGDRhoDKappa;
}

void
ConcreteDPM2 :: computeDDGDDInv1d(double &answer,
                                  const double sigma,
                                  const double tempKappa)
{
    //compute yieldHardOne and yieldSoft
    const double yieldHardOne = computeHardeningOne(tempKappa);
    const double yieldHardTwo = computeHardeningTwo(tempKappa);
    double AGParam = this->ft * yieldHardTwo * 3 / this->fc + m / 2;
    double BGParam =
        yieldHardTwo / 3. * ( 1. + this->ft / this->fc ) /
        ( log(AGParam) + log(this->dilationConst + 1.) - log(2 * this->dilationConst - 1.) - log(3. * yieldHardTwo + this->m / 2) );
    const double R = ( sigma - ft  * yieldHardTwo ) / ( 3 * fc * BGParam );
    const double dMQDSigma = AGParam / ( 9 * BGParam * fc ) * exp(R);
    answer = 2 * pow(1 / fc + 8 * sigma / pow(3 * fc, 2) * ( 1 - yieldHardOne ), 2) + pow(4 / 3 / fc, 2) * ( sigma / fc + pow(2 / 3 * sigma / fc, 2) * ( 1 - yieldHardOne ) ) * ( 1 - yieldHardOne ) + pow(yieldHardOne, 2) / fc * dMQDSigma;
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
                            pow(yieldHardOne, 2.) * dMQDSig / fc;

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
    StressVector principalStress(_3dMat);
    effectiveStress.computePrincipalValDir(principalStress, stressPrincipalDir);

    //Split the principal values in a tension and a compression part
    StressVector principalStressTension(_3dMat);
    StressVector principalStressCompression(_3dMat);

    for ( int i = 1; i <= principalStress.giveSize(); i++ ) {
        if ( principalStress.at(i) >= 0 ) {
            principalStressTension.at(i) = principalStress.at(i);
            principalStressCompression.at(i) = 0.;
        } else {
            principalStressTension.at(i) = 0.;
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
    } else {
        effectiveStressTension = stressTemp;
    }

    //Take care of type of stress state for compression
    transformStressVectorTo(stressTemp, stressPrincipalDir, principalStressCompression, 1);
    if ( matMode == _PlaneStrain ) {
        effectiveStressCompression(0) = stressTemp(0);
        effectiveStressCompression(1) = stressTemp(1);
        effectiveStressCompression(2) = stressTemp(2);
        effectiveStressCompression(3) = stressTemp(5);
    } else {
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
            linkomb = tmp.at(j, i) / piv;
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
    } else if ( kappa > 0. && kappa < 1. ) {
        return
            ( 1. - yieldHardInitial - yieldHardPrimePeak ) * pow(kappa, 3.)
            - ( 3. * ( 1. - yieldHardInitial ) - 3. * yieldHardPrimePeak ) * pow(kappa, 2.)
            + ( 3. * ( 1. - yieldHardInitial ) - 2. * yieldHardPrimePeak ) * kappa
            + yieldHardInitial;
    } else {
        return 1.;
    }
}


double
ConcreteDPM2 :: computeHardeningOnePrime(const double kappa) const
{
    if ( kappa <= 0. ) {
        return 3. * ( 1 - yieldHardInitial ) - 2. * yieldHardPrimePeak;
    } else if ( kappa >= 0. && kappa < 1. ) {
        return
            3. * ( 1. - yieldHardInitial - yieldHardPrimePeak ) * pow(kappa, 2.)
            - 2. * ( 3. * ( 1. - yieldHardInitial ) - 3. * yieldHardPrimePeak ) * kappa
            + ( 3. * ( 1. - yieldHardInitial ) - 2. * yieldHardPrimePeak );
    } else {
        return 0.;
    }
}


double
ConcreteDPM2 :: computeHardeningTwo(const double kappa) const
{
    if ( kappa <= 0. ) {
        return 1.;
    } else if ( kappa > 0. && kappa < 1. ) {
        return 1.;
    } else {
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
    } else {
        return yieldHardPrimePeak;
    }
}

void
ConcreteDPM2 :: give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    ConcreteDPM2Status *status = ( ConcreteDPM2Status * ) this->giveStatus(gp);
    double om;
    if ( mode == ElasticStiffness ) {
        answer.resize(1, 1);
        answer.at(1, 1) = eM;
    } else if ( mode == SecantStiffness ) {
        StrainVector strain(status->giveTempStrainVector(), _1dMat);
        if ( strain.at(1) >= 0 || isotropicFlag == 1 ) {
            om = status->giveTempDamageTension();
        } else {
            om = status->giveTempDamageCompression();
        }

        if ( om > 0.999999 ) {
            om = 0.999999;
        }

        answer.resize(1, 1);
        answer.at(1, 1) = eM;
        answer.times(1.0 - om);
    } else {
        OOFEM_WARNING("ConcreteDPM2 :: give1dStressStiffMtrx ... unknown type of stiffness (tangent stiffness not implemented for 1d). Elastic stiffness used!\n");
        answer.resize(1, 1);
        answer.at(1, 1) = eM;
    }


    return;
}

void
ConcreteDPM2 :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                              MatResponseMode mode,
                                              GaussPoint *gp,
                                              TimeStep *tStep)
{
    if ( gp->giveMaterialMode() == _3dMat ||
         gp->giveMaterialMode() ==  _PlaneStrain  ) {
        if ( mode == ElasticStiffness ) {
            FloatMatrix stiff;
            this->giveLinearElasticMaterial()->giveStiffnessMatrix(stiff, mode, gp, tStep);
            this->giveFullSymMatrixForm( answer, stiff, gp->giveMaterialMode() );
        } else if ( mode == SecantStiffness ) {
            computeSecantStiffness(answer, mode, gp, tStep);
        } else if ( mode == TangentStiffness ) {
            OOFEM_WARNING("unknown type of stiffness (tangent stiffness not implemented for 1d). Elastic stiffness used!");
            FloatMatrix stiff;
            this->giveLinearElasticMaterial()->giveStiffnessMatrix(stiff, mode, gp, tStep);
            this->giveFullSymMatrixForm( answer, stiff, gp->giveMaterialMode() );
        }
    }
}

void
ConcreteDPM2 :: computeSecantStiffness(FloatMatrix &answer,
                                       MatResponseMode mode,
                                       GaussPoint *gp,
                                       TimeStep *tStep)
{
    this->matMode = gp->giveMaterialMode();
    ConcreteDPM2Status *status = static_cast< ConcreteDPM2Status * >( giveStatus(gp) );

    double omegaTension;

    //  Damage parameters
    omegaTension = status->giveTempDamageTension();
    if ( omegaTension > 0.999999 ) {
        omegaTension = 0.999999;
    }

    FloatMatrix stiff;
    this->giveLinearElasticMaterial()->giveStiffnessMatrix(stiff, mode, gp, tStep);
    this->giveFullSymMatrixForm( answer, stiff, gp->giveMaterialMode() );

    if ( isotropicFlag == 1 ) {
        answer.times(1. - omegaTension);
        return;
    }


    StrainVector strain(status->giveTempStrainVector(), matMode);
    StressVector effectiveStress(matMode);
    StrainVector elasticStrain = strain;
    elasticStrain.subtract( status->giveTempPlasticStrain() );
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
    if ( matMode == _1dMat ) { //1d case
        sigNew = stress.computeFirstInvariant() / 3.;
        rhoNew = stress.computeFirstInvariant() * sqrt(2. / 3.);
        if ( sigNew >= 0 ) {
            thetaNew = 0.;
        } else {
            thetaNew = M_PI / 6;
        }
    } else {
        StressVector deviatoricStress(matMode);
        stress.computeDeviatoricVolumetricSplit(deviatoricStress, sigNew);
        rhoNew = deviatoricStress.computeSecondCoordinate();
        thetaNew = deviatoricStress.computeThirdCoordinate();
    }
}


void
ConcreteDPM2 :: assignStateFlag(GaussPoint *gp)
{
    ConcreteDPM2Status *status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

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
        } else {
            status->letTempStateFlagBe(ConcreteDPM2Status :: ConcreteDPM2_Plastic);
        }
    } else {
        const int state_flag = status->giveStateFlag();
        if ( tempDamageTension > damageTension || tempDamageTension == 1. ||
             tempDamageCompression > damageCompression || tempDamageCompression == 1. ) {
            status->letTempStateFlagBe(ConcreteDPM2Status :: ConcreteDPM2_Damage);
        } else {
            if ( state_flag == ConcreteDPM2Status :: ConcreteDPM2_Elastic ) {
                status->letTempStateFlagBe(ConcreteDPM2Status :: ConcreteDPM2_Elastic);
            } else {
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
                            TimeStep *tStep)
{
    ConcreteDPM2Status *status = static_cast< ConcreteDPM2Status * >( this->giveStatus(gp) );

    switch ( type ) {
    case IST_PlasticStrainTensor:
    {
        StrainVector e(_Unknown);
        status->giveFullPlasticStrainVector(e);
        answer = e;
    }
        return 1;

    case IST_DamageTensor:
        answer.resize(6);
        answer.zero();
        answer.at(1) = status->giveDamageTension();
        answer.at(2) = status->giveDamageCompression();
        answer.at(3) = status->giveDissWork();
        return 1;

#ifdef keep_track_of_dissipated_energy
    case IST_StressWorkDensity:
        answer.resize(1);
        answer.at(1) = status->giveStressWork();
        return 1;

    case IST_DissWorkDensity:
        answer.resize(1);
        answer.at(1) = status->giveDissWork();
        return 1;

    case IST_FreeEnergyDensity:
        answer.resize(1);
        answer.at(1) = status->giveStressWork() - status->giveDissWork();
        return 1;

#endif
    default:
        return StructuralMaterial :: giveIPValue(answer, gp, type, tStep);
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
