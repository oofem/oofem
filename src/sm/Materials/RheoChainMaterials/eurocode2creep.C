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

#include "eurocode2creep.h"
#include "mathfem.h"
#include "gausspoint.h"
#include "timestep.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"
#include "engngm.h"

namespace oofem {
REGISTER_Material(Eurocode2CreepMaterial);

IRResultType
Eurocode2CreepMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, fcm28, _IFT_Eurocode2CreepMaterial_fcm28);

    /* scaling factor transforming PREDICTED stiffness Ecm
     *  e.g. if the stiffness should be in MPa, then stiffnessFactor = 1.e6
     *  e.g. if the stiffness should be in Pa, then stiffnessFactor = 1. */
    this->stiffnessFactor = 1.e6;
    IR_GIVE_OPTIONAL_FIELD(ir, stiffnessFactor, _IFT_Eurocode2CreepMaterial_stiffnessFactor);

    IR_GIVE_FIELD(ir, t0, _IFT_Eurocode2CreepMaterial_t0); // duration of curing (in days)
    IR_GIVE_FIELD(ir, h0, _IFT_Eurocode2CreepMaterial_h0); // equivalent thickness, in mm !!!

    this->retardationSpectrumApproximation = false;
    if ( ir->hasField(_IFT_Eurocode2CreepMaterial_spectrum) ) {
        this->retardationSpectrumApproximation = true;
    }

    this->temperatureDependent = false;
    if ( ir->hasField(_IFT_Eurocode2CreepMaterial_temperatureDependent) ) {
        this->temperatureDependent = true;
    }


    int cemType;
    IR_GIVE_FIELD(ir, cemType, _IFT_Eurocode2CreepMaterial_cemType);

    double henv = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, henv, _IFT_Eurocode2CreepMaterial_henv);

    // read shrinkage type
    int sht = 0;
    IR_GIVE_FIELD(ir, sht, _IFT_Eurocode2CreepMaterial_shType);
    if ( sht > 3 ) {
        OOFEM_ERROR("unsupported shrinkage type");
    }
    this->shType = ( ec2ShrinkageType ) sht;

    this->computeElasticityStrengthParams(cemType);
    this->computeShrinkageParams(cemType, henv);
    this->computeCreepParams(cemType, henv);

    // the initialization of the upper class has to be done here
    // because RheoChainMat calls at initialization computeCharTimes
    // but it does not know which method should be chosen (retardation
    // spectrum switch is initialized here)
    // initialize Kelvin chain aging material
    KelvinChainMaterial :: initializeFrom(ir);

    return IRRT_OK;
}



void
Eurocode2CreepMaterial :: computeElasticityStrengthParams(int cemType)
{
    // Ecm28 in GPa (Eurocode formula given in Table 3.1
    this->Ecm28 = 22. *  pow( ( this->fcm28 * ( this->stiffnessFactor / 1.e6 ) / 10. ), 0.3 );
    this->Ecm28 *= 1.e9 / this->stiffnessFactor; // converted to stiffness units of the analysis

    // "s" is given in 3.2 in section 3.2.1
    if ( cemType == 1 ) { // class R
        this->s = 0.2;
    } else if ( cemType == 2 ) { // class N
        this->s = 0.25;
    } else if ( cemType == 3 ) { // class S
        this->s = 0.38;
    } else {
        OOFEM_ERROR("unsupported value of cement type");
    }
}


void
Eurocode2CreepMaterial :: computeShrinkageParams(int cemType, double henv)
{
    // drying shrinkage (used in B.11)
    double alpha_ds1, alpha_ds2;

    if ( cemType == 1 ) { // class R
        alpha_ds1 = 6;
        alpha_ds2 = 0.11;
    } else if ( cemType == 2 ) { // class N
        alpha_ds1 = 4;
        alpha_ds2 = 0.12;
    } else if ( cemType == 3 ) { // class S
        alpha_ds1 = 3;
        alpha_ds2 = 0.13;
    } else {
        OOFEM_ERROR("unsupported value of cement type");
    }


    if ( ( this->shType == EC2_TotalShrinkage ) || ( this->shType == EC2_DryingShrinkage ) ) {
        // kh from table 3.3
        if ( this->h0 >= 500 ) {
            this->kh = 0.7;
        } else if ( this->h0 >= 300 ) {
            this->kh = 0.75 - 0.05 * ( h0 - 300. ) / 200.;
        } else if ( this->h0 >= 200 ) {
            this->kh = 0.85 - 0.10 * ( h0 - 200. ) / 100.;
        } else if ( this->h0 >= 100 ) {
            this->kh = 1.00 - 0.15 * ( h0 - 100. ) / 100.;
        } else {
            this->kh = 1.;
        }


        double beta_RH = 1.55 * ( 1. - henv * henv * henv );
        this->eps_cd_0 = -0.85e-6 * ( ( 220. + 110. * alpha_ds1 ) * exp(-alpha_ds2 * this->fcm28 * ( this->stiffnessFactor / 1.e6 ) / 10.) * beta_RH );
    }

    if ( ( this->shType == EC2_TotalShrinkage ) || ( this->shType == EC2_AutogenousShrinkage ) ) {
        // characteristic strength in MPa
        double fck = this->fcm28 * ( this->stiffnessFactor / 1.e6 ) - 8.;
        // below 10 MPa no autogenous shrinkage
        fck = max(fck, 10.);

        this->eps_ca_infty = -2.5e-6 * ( fck - 10. ); // 3.12
    }

    return;
}

void
Eurocode2CreepMaterial :: computeCreepParams(int cemType, double henv)
{
    // "alpha_T_cement" influences the equivalent age (B.9)

    if ( cemType == 1 ) { // class R
        this->alpha_T_cement = 1.;
    } else if ( cemType == 2 ) { // class N
        this->alpha_T_cement = 0.;
    } else if ( cemType == 3 ) { // class S
        this->alpha_T_cement = -1.;
    } else {
        OOFEM_ERROR("unsupported value of cement type");
    }

    if ( this->fcm28 * ( this->stiffnessFactor / 1.e6 )  > 35. ) {
        double alpha_1, alpha_2, alpha_3;

        // B.8c
        alpha_1 = pow( ( 35. / ( this->fcm28 * ( this->stiffnessFactor / 1.e6 ) ) ), 0.7 );
        alpha_2 = pow( ( 35. / ( this->fcm28 * ( this->stiffnessFactor / 1.e6 ) ) ), 0.2 );
        alpha_3 = pow( ( 35. / ( this->fcm28 * ( this->stiffnessFactor / 1.e6 ) ) ), 0.5 );

        // influence of relative humidity on drying creep
        // B.3b
        this->phi_RH = ( 1. + alpha_1 * ( 1. - henv ) / ( 0.1 * pow(this->h0, 1. / 3.) ) ) * alpha_2;

        // B.8b
        // in EC there is 0.012 * RH, it is equivalent to 1.2 * henv
        this->beta_H = 1.5 * this->h0 * ( 1. + pow( ( 1.2 * henv ), 18. ) ) + 250. * alpha_3;
        this->beta_H = min(this->beta_H, 1500. * alpha_3);

        // convert to proper time units
        //this->beta_H *=  this->timeFactor;
        // this can't be done here because timeFactor has not been initialized yet! (done in RheoChM)
        // beta_H remains in days
    } else {
        // influence of relative humidity on drying creep
        // B.3a
        this->phi_RH = 1. + ( 1. - henv ) / ( 0.1 * pow(this->h0, 1. / 3.) );

        // B.8b
        // in EC there is 0.012 * RH, it is equivalent to 1.2 * henv
        this->beta_H = 1.5 * this->h0 * ( 1. + pow( ( 1.2 * henv ), 18. ) ) + 250.;
        this->beta_H = min(this->beta_H, 1500.);

        // convert to proper time units
        // this->beta_H *=  this->timeFactor;
        // this can't be done here because timeFactor has not been initialized yet! (done in RheoChM)
    }

    // B.4
    this->beta_fcm = 16.8 / sqrt( this->fcm28 * ( this->stiffnessFactor / 1.e6 ) );

    return;
}


double
Eurocode2CreepMaterial :: computeConcreteStrengthAtAge(double age)
{
    double beta = exp( this->s * ( 1. - sqrt( 28. / ( age / this->timeFactor ) ) ) );

    return beta * this->fcm28;
}

double
Eurocode2CreepMaterial :: computeMeanElasticModulusAtAge(double age)
{
    double fcm_at_age;
    double Ecm_at_age;

    fcm_at_age = this->computeConcreteStrengthAtAge(age);
    Ecm_at_age = pow( ( fcm_at_age / this->fcm28 ), 0.3 ) * this->Ecm28;

    return Ecm_at_age;
}

double
Eurocode2CreepMaterial :: computeCreepFunction(double t, double t_prime, GaussPoint *gp, TimeStep *tStep)
// compute the value of the creep function of concrete
// corresponding to the given load duration (in days)
// t-t_prime = duration of loading

{
    double J;     //return value


    // elastic modulus isn not temperature adjusted too
    J = 1. / this->computeMeanElasticModulusAtAge(t_prime);

    J += this->computeCreepCoefficient(t, t_prime, gp, tStep) / ( 1.05 * this->Ecm28 );

    return J;
}


double
Eurocode2CreepMaterial :: computeCreepCoefficient(double t, double t_prime, GaussPoint *gp, TimeStep *tStep)
// compute the value of the creep coefficient
// corresponding to the given load duration (in days)
// t-t_prime = duration of loading

{
    double phi;     //return value

    double beta_t0, beta_c;
    double t_prime_equiv;

    if ( temperatureDependent ) {
        t_prime_equiv = this->computeEquivalentAge(gp, tStep);
    } else {
        t_prime_equiv = t_prime;
    }

    // B.5
    beta_t0 = 1. / ( 0.1 + pow(t_prime_equiv / this->timeFactor, 0.2) );

    // B.7
    beta_c = pow( ( ( t - t_prime ) / ( this->beta_H * this->timeFactor + t - t_prime ) ), 0.3 );

    // t-t'_prime should not be temperature-adjusted (so says the EC)

    // B.1 + B.2
    phi = this->phi_RH * this->beta_fcm * beta_t0 * beta_c;


    return phi;
}



// implements B.10
double
Eurocode2CreepMaterial :: computeEquivalentMaturity(GaussPoint *gp, TimeStep *tStep) {
    Eurocode2CreepMaterialStatus *status = static_cast< Eurocode2CreepMaterialStatus * >( this->giveStatus(gp) );

    FloatArray et;

    double initMaturity, incrMaturity;
    double averageTemperature;

    //  Element *elem = gp->giveElement();
    StructuralElement *selem = dynamic_cast< StructuralElement * >( gp->giveElement() );

    selem->computeResultingIPTemperatureAt(et, tStep, gp, VM_Total);

    if ( tStep->isTheFirstStep() || ( tStep->giveIntrinsicTime() - tStep->giveTimeIncrement() - this->castingTime < 0. ) ) {
        initMaturity = this->relMatAge;
        averageTemperature = et.at(1);
    } else {
        initMaturity  = status->giveConcreteMaturity();
        averageTemperature = ( et.at(1) + status->giveTemperature() ) / 2.;
    }

    averageTemperature = min(averageTemperature, 80.);
    averageTemperature = max(averageTemperature, 0.);

    incrMaturity = exp( -1. * ( 4000. / ( 273. + averageTemperature ) - 13.65 ) ) * tStep->giveTimeIncrement();

    status->setTempConcreteMaturity(initMaturity + incrMaturity); // full increment - for updating
    status->setTempTemperature( et.at(1) );

    // middle of the step
    return initMaturity + incrMaturity / 2.;
}


// implements B.9
double
Eurocode2CreepMaterial :: computeEquivalentAge(GaussPoint *gp, TimeStep *tStep) {
    double maturity;
    double age;

    maturity = this->computeEquivalentMaturity(gp, tStep);

    age = maturity * pow( ( 9. / ( 2. + pow(maturity, 1.2) ) + 1. ), this->alpha_T_cement );
    age = max(age, 0.5);

    return age;
}



double
Eurocode2CreepMaterial :: giveEModulus(GaussPoint *gp, TimeStep *tStep)
{
    /*
     * This function returns the incremental modulus for the given time increment.
     * The modulus may also depend on the specimen geometry (gp - dependence).
     *
     * It is stored as "Einc" for further expected requests from other gaussPoints that correspond to the same material.
     *
     * Note: time -1 refers to the previous time.
     */

    // !!! chartime exponents are assumed to be equal to 1 !!!

    double chainStiffness = 0.;

    if ( !this->isActivated(tStep) ) {
        return 1.; // stresses are cancelled in giveRealStressVector;
    }

    // updateEparModuli is called in KelvinChainMaterial
    chainStiffness = KelvinChainMaterial :: giveEModulus(gp, tStep);

    if ( retardationSpectrumApproximation  ) { //retardation spectrum used
        double sum;
        double t_halfstep;

        sum = 1. / chainStiffness;     //  convert stiffness into compliance

        sum += 1 / this->EspringVal; // add zeroth unit

        t_halfstep = this->relMatAge + ( tStep->giveTargetTime() - 0.5 * tStep->giveTimeIncrement() );

        if ( t_halfstep <= 0. ) {
            OOFEM_ERROR("attempt to evaluate material stiffness at negative age");
        }

        sum += 1. / this->computeMeanElasticModulusAtAge(t_halfstep); // add initial compliance

        // convert to stiffness
        chainStiffness = 1. / sum;
    }

    return chainStiffness;
}


void
Eurocode2CreepMaterial :: computeCharTimes()
{
    /*
     * This function generates discrete characteristic times
     * according to rules which guarantee a good approximation
     * of the relaxation or creep function by the Dirichlet series
     */

    int j;

    if ( this->begOfTimeOfInterest == -1. ) {
        this->begOfTimeOfInterest = 0.001 * this->timeFactor;
    }

    if ( this->endOfTimeOfInterest == -1. ) {
        this->endOfTimeOfInterest = 10000. * this->timeFactor;
    }


    // perhaps different values should be used for moduli computed by the least-squares method
    // this gives really big mistakes for times near the interest boundary

    if ( retardationSpectrumApproximation  ) { //retardation spectrum used
        if ( this->begOfTimeOfInterest > 1.0 * this->timeFactor ) {
            OOFEM_WARNING("begOfTimeOfInterest was chosen bigger than 1 days, reseting its value to 1 days (could have lead to big errors in the numerical integration of the stiffness of the zeroth Kelvin unit (the retardation spectrum is very steep)");
            this->begOfTimeOfInterest = 1.0 * this->timeFactor;
        }

        this->tau1 = this->begOfTimeOfInterest;
    } else {
        this->tau1 = this->begOfTimeOfInterest / 10.;
    }

    //first retardation can be treated as equal to begOfTimeOfInterest



    if ( retardationSpectrumApproximation  ) { //retardation spectrum used
        //the last retardation time has to be bigger than approx sqrt(10) * endOfTimeOfInterest
        // e.g. if endoftimeofinterest is 10^2, then the last retardation time has to be at least 10^2.5
        if ( this->endOfTimeOfInterest > pow(10., 5.5) * timeFactor ) {
            OOFEM_WARNING("endOfTimeOfInterest was chosen bigger than 10.000 days, reseting to 10.000 days (the retardation spectrum is almost zero afterwards)");
            this->endOfTimeOfInterest = pow(10., 5.5) * timeFactor;
        }
    }

    j = 1;
    while ( this->endOfTimeOfInterest >=  this->tau1 * pow(10.0, ( double ) j - 1. - 0.5) ) {
        j++;
    }

    this->nUnits = j;

    this->charTimes.resize(this->nUnits);
    this->charTimes.zero();

    for ( int mu = 1; mu <= this->nUnits; mu++ ) {
        this->charTimes.at(mu) = this->tau1 * pow(10., mu - 1);

        if ( retardationSpectrumApproximation ) { // apply correction factors
            this->charTimes.at(mu) *= this->computeRetardationTimeCorrection(mu);
        }
    }
}



double
Eurocode2CreepMaterial :: computeRetardationTimeCorrection(int mu)
{
    return 1. + 0.555 * exp( -4. * pow( ( this->tau1 * pow( 10., double( mu - 1 ) ) / ( this->beta_H * this->timeFactor ) ), 2. ) );


    //return  1.;
}


double
Eurocode2CreepMaterial :: evaluateSpectrumAt(double tau)
{
    // second-order approximation of the retardation spectrum
    double L;

    double t = 2 * tau;
    double n = 0.3; // exponent in the model

    double f, df, ddf; // function "t/(beta+t)" and its derivatives
    double ddPhi; // second derivative of (f(t))^n

    f = t / ( this->beta_H * this->timeFactor + t );

    df = this->beta_H * this->timeFactor / ( ( this->beta_H * this->timeFactor + t ) * ( this->beta_H * this->timeFactor + t ) );

    ddf = -2. * this->beta_H * this->timeFactor / ( ( this->beta_H * this->timeFactor + t ) * ( this->beta_H * this->timeFactor + t ) * ( this->beta_H * this->timeFactor + t ) );

    ddPhi = n * ( n - 1. ) * pow( f, ( n - 2. ) ) * df * df + n *pow( f, ( n - 1. ) ) * ddf;

    L = -4. * tau * tau * ddPhi;

    return L;
}


void
Eurocode2CreepMaterial :: computeCharCoefficients(FloatArray &answer, double atTime, GaussPoint *gp, TimeStep *tStep)
{
    /*
     * If retardationSpectrumApproximation == true then analysis of continuous retardation spectrum is used for
     * computing characteristic coefficients (moduli) of Kelvin chain
     * Else least-squares method is used
     */
    if ( retardationSpectrumApproximation ) {
        // all moduli must be multiplied by g(t') / c - see equation (46) in Jirasek's retardation spectrum paper

        double coefficient;
        double equivalentAge;

        if ( temperatureDependent ) {
            equivalentAge = this->computeEquivalentAge(gp, tStep);
        } else {
            equivalentAge = atTime;
        }

        coefficient = 1.05 * this->Ecm28 * ( 0.1 + pow(equivalentAge / this->timeFactor, 0.2) ) / ( this->phi_RH * this->beta_fcm );


        // evaluate stiffness of the zero-th unit of the Kelvin chain
        // (aging elastic spring with retardation time = 0)
        // this is done employing Simpson's rule. the begOfTimeOfInterest cannot exceed 0.1 day
        // E0 = EspringVal = int( L, 0, tau1/sqrt(10) )

        double tau0 = this->tau1 / sqrt(10.0); // upper bound of the integral

        this->EspringVal = 1. / ( ( log(10.) / 3. ) * (
                                      this->evaluateSpectrumAt(tau0 * 1.e-8) + 4. * this->evaluateSpectrumAt(tau0 * 1.e-7) +
                                      2. * this->evaluateSpectrumAt(tau0 * 1.e-6) + 4. * this->evaluateSpectrumAt(tau0 * 1.e-5) +
                                      2. * this->evaluateSpectrumAt(tau0 * 1.e-4) + 4. * this->evaluateSpectrumAt(tau0 * 1.e-3) +
                                      2. * this->evaluateSpectrumAt(tau0 * 1.e-2) + 4. * this->evaluateSpectrumAt(tau0 * 1.e-1) +
                                      this->evaluateSpectrumAt(tau0) ) );

        this->EspringVal *= coefficient;

        // process remaining units
        answer.resize(nUnits);
        answer.zero();

        double tauMu;

        for ( int mu = 1; mu <= this->nUnits; mu++ ) {
            tauMu = this->tau1 * pow( 10., double( mu - 1 ) );

            answer.at(mu) =  2. / ( log(10.) * ( this->evaluateSpectrumAt( tauMu * pow(10., -sqrt(3.) / 6.) ) + this->evaluateSpectrumAt( tauMu * pow(10., sqrt(3.) / 6.) ) ) );
        }

        answer.times(coefficient);
    } else {   // moduli computed using the least-squares method
        KelvinChainMaterial :: computeCharCoefficients(answer, atTime, gp, tStep);
    }
}


void
Eurocode2CreepMaterial :: giveShrinkageStrainVector(FloatArray &answer,
                                                    GaussPoint *gp,
                                                    TimeStep *tStep,
                                                    ValueModeType mode)
{
    answer.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
    answer.zero();

    if ( ( this->shType == EC2_NoShrinkage ) || ( !this->isActivated(tStep) ) ) {
        return;
    }

    if ( ( mode != VM_Total ) && ( mode != VM_Incremental ) ) {
        OOFEM_ERROR("unsupported mode");
    }

    FloatArray eps_sh;

    // the relative matereial age can be negative to properly reflect the casting time
    double dryingTimeNow = ( tStep->giveTargetTime() + this->relMatAge - this->t0 ) / timeFactor;
    double autoShrTimeNow = ( this->relMatAge + tStep->giveTargetTime() ) / this->timeFactor;

    if ( mode == VM_Incremental ) {
        double dt = tStep->giveTimeIncrement() / this->timeFactor;

        if ( dryingTimeNow > 0. ) {
            if ( ( this->shType == EC2_TotalShrinkage ) || ( this->shType == EC2_DryingShrinkage ) ) {
                this->computeIncrementOfDryingShrinkageVector(eps_sh, gp, dryingTimeNow, dryingTimeNow - dt);
                answer.add(eps_sh);
            }
        }

        if ( autoShrTimeNow > ( this->castingTime + this->relMatAge ) ) {
            if ( ( this->shType == EC2_TotalShrinkage ) || ( this->shType == EC2_AutogenousShrinkage ) ) {
                this->computeIncrementOfAutogenousShrinkageVector(eps_sh, gp, autoShrTimeNow, autoShrTimeNow - dt);
                answer.add(eps_sh);
            }
        }
    } else { // total
      if ( dryingTimeNow > 0. ) {
        if ( ( this->shType == EC2_TotalShrinkage ) || ( this->shType == EC2_DryingShrinkage ) ) {
	  this->computeIncrementOfDryingShrinkageVector( eps_sh, gp, dryingTimeNow, max(0., ( this->relMatAge - this->t0 ) / this->timeFactor) );
	  answer.add(eps_sh);
        }
      }
	
      if ( autoShrTimeNow > ( this->castingTime + this->relMatAge ) ) {
	if ( ( this->shType == EC2_TotalShrinkage ) || ( this->shType == EC2_AutogenousShrinkage ) ) {
	  if (this->castingTime > 0.) {
	    this->computeIncrementOfAutogenousShrinkageVector(eps_sh, gp, autoShrTimeNow, (this->relMatAge+this->castingTime) / this->timeFactor);
	  } else {
	    this->computeIncrementOfAutogenousShrinkageVector(eps_sh, gp, autoShrTimeNow, (this->relMatAge) / this->timeFactor);
	  }
	  answer.add(eps_sh);
	}
      }
    }

    
    if ( answer.at(1) != answer.at(1) ) {
      OOFEM_ERROR("shrinkage is NaN: %f", answer.at(1));
    }

    return;
}

void
Eurocode2CreepMaterial :: computeIncrementOfDryingShrinkageVector(FloatArray &answer, GaussPoint *gp, double tNow, double tThen)
{
    int size;
    FloatArray fullAnswer;
    MaterialMode mode = gp->giveMaterialMode();

    if ( ( mode == _3dShell ) || ( mode ==  _3dBeam ) || ( mode == _2dPlate ) || ( mode == _2dBeam ) ) {
        size = 12;
    } else {
        size = 6;
    }

    fullAnswer.resize(size);
    fullAnswer.zero();

    if ( tNow > tThen ) {
        double dEpsSh;

        // B.10
        double beta_ds_now, beta_ds_then;

        beta_ds_now = tNow / ( tNow + 0.04 * pow(this->h0, 3. / 2.) );
        beta_ds_then = tThen / ( tThen + 0.04 * pow(this->h0, 3. / 2.) );

        dEpsSh = ( beta_ds_now - beta_ds_then ) * this->kh * this->eps_cd_0;

        fullAnswer.at(1) = fullAnswer.at(2) = fullAnswer.at(3) = dEpsSh;
    }

    StructuralMaterial :: giveReducedSymVectorForm( answer, fullAnswer, gp->giveMaterialMode() );
}


void
Eurocode2CreepMaterial :: computeIncrementOfAutogenousShrinkageVector(FloatArray &answer, GaussPoint *gp, double tNow, double tThen)
{
    int size;
    FloatArray fullAnswer;
    MaterialMode mode = gp->giveMaterialMode();

    if ( ( mode == _3dShell ) || ( mode ==  _3dBeam ) || ( mode == _2dPlate ) || ( mode == _2dBeam ) ) {
        size = 12;
    } else {
        size = 6;
    }

    fullAnswer.resize(size);
    fullAnswer.zero();


    if ( tNow > tThen ) {
        double dEpsAu;

        // B.13
        double beta_as_now, beta_as_then;

        beta_as_now = 1. - exp( -0.2 * sqrt(tNow) );
        beta_as_then = 1. - exp( -0.2 * sqrt(tThen) );

        dEpsAu = ( beta_as_now - beta_as_then ) * this->eps_ca_infty;

        fullAnswer.at(1) = fullAnswer.at(2) = fullAnswer.at(3) = dEpsAu;
    }

    StructuralMaterial :: giveReducedSymVectorForm( answer, fullAnswer, gp->giveMaterialMode() );
}


MaterialStatus *
Eurocode2CreepMaterial :: CreateStatus(GaussPoint *gp) const
/*
 * creates a new material status corresponding to this class
 */
{
    return new Eurocode2CreepMaterialStatus(1, this->giveDomain(), gp, nUnits);
}


void
Eurocode2CreepMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    KelvinChainMaterial :: giveRealStressVector(answer, gp, reducedStrain, tStep);

    if ( !this->isActivated(tStep) ) {
        return;
    }
}


/****************************************************************************************/
/**********     Eurocode2CreepMaterialStatus - HUMIDITY ****************************************/


Eurocode2CreepMaterialStatus :: Eurocode2CreepMaterialStatus(int n, Domain *d, GaussPoint *g, int nunits) :
    KelvinChainMaterialStatus(n, d, g, nunits) {
    tempMaturity = 0.;
    maturity = tempMaturity;

    tempTemperature = 0.;
    temperature = tempTemperature;
}

void
Eurocode2CreepMaterialStatus :: updateYourself(TimeStep *tStep)
{
    maturity = tempMaturity;
    tempMaturity = 0.;

    temperature = tempTemperature;
    tempTemperature = 0.;

    KelvinChainMaterialStatus :: updateYourself(tStep);
}


contextIOResultType
Eurocode2CreepMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
//
{
    contextIOResultType iores;

    if ( ( iores = KelvinChainMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write maturity value
    if ( !stream.write(maturity) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // write temperature value
    if ( !stream.write(temperature) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

contextIOResultType
Eurocode2CreepMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restore the state variables from a stream
//
{
    contextIOResultType iores;
    if ( ( iores = KelvinChainMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // read maturity value
    if ( !stream.read(maturity) ) {
        return CIO_IOERR;
    }

    // read temperature value
    if ( !stream.read(temperature) ) {
        return CIO_IOERR;
    }


    return CIO_OK;
}
} // end namespace oofem
