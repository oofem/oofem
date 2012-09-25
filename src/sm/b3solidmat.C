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

#include "b3solidmat.h"
#include "mathfem.h"
#include "gausspnt.h"
#include "structuralcrosssection.h"
#include "timestep.h"
#include "contextioerr.h"
#include "datastream.h"

namespace oofem {
IRResultType
B3SolidMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    //
    // NOTE
    //
    // this material model is unit-dependent!
    // units must be in MPa, m???, MN.
    //

    double fc, c, wc, ac, alpha1 = 1.0, alpha2 = 1.0;
    double initHum; //MPS shrinkage parameter - initial value of humidity (range 0.2-0.98)
    double finalHum; //MPS shrinkage parameter - final value of humidity (range 0.2-0.98)

    // EmoduliMode has to be set first because characteristic times are computed from RheoChM initialization
    this->EmoduliMode = 0;
    // retardation spectrum or least square method is used. Retardation spectrum is default (EmoduliMode==0)
    IR_GIVE_OPTIONAL_FIELD(ir, EmoduliMode, IFT_B3Material_emodulimode, "emodulimode");

    KelvinChainMaterial :: initializeFrom(ir);

    int mode = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, mode, IFT_B3Material_mode, "mode");

    if ( mode == 0 ) { // default - estimate model parameters q1,..,q5 from composition
        IR_GIVE_FIELD(ir, fc, IFT_B3Material_fc, "fc"); // 28-day standard cylinder compression strength [MPa]
        IR_GIVE_FIELD(ir,  c, IFT_B3Material_cc, "cc"); // cement content of concrete [kg/m^3]
        IR_GIVE_FIELD(ir, wc, IFT_B3Material_wc, "w/c"); // ratio (by weight) of water to cementitious material
        IR_GIVE_FIELD(ir, ac, IFT_B3Material_ac, "a/c"); // ratio (by weight) of aggregate to cement
        IR_GIVE_FIELD(ir, t0, IFT_B3Material_t0, "t0"); // age when drying begins [days]
    } else { // read model parameters for creep
        IR_GIVE_FIELD(ir, q1, IFT_B3Material_q1, "q1");
        IR_GIVE_FIELD(ir, q2, IFT_B3Material_q2, "q2");
        IR_GIVE_FIELD(ir, q3, IFT_B3Material_q3, "q3");
        IR_GIVE_FIELD(ir, q4, IFT_B3Material_q4, "q4");
    }

    // default value = 0; it can be used for basic creep without any external fields
    this->MicroPrestress = 0; //if = 1 computation exploiting Microprestress solidification theory is done;
    IR_GIVE_OPTIONAL_FIELD(ir, MicroPrestress, IFT_B3Material_microprestress, "microprestress");
    // is MPS theory used?
    if ( this->MicroPrestress == 1 ) {
        //?????????????????????????? co z nize uvedenych parametru dat volitelne a jake hodnoty?
        // microprestress-sol-theory: read data for microprestress evaluation
        // constant c0 [MPa^-1 * day^-1]
        IR_GIVE_FIELD(ir, c0, IFT_B3Material_c0, "c0");
        // constant c1 (=C1*R*T/M)
        IR_GIVE_FIELD(ir, c1, IFT_B3Material_c1, "c1");
        // tS0- necessary for the initial value of microprestress = S0 (age when drying begins)
        IR_GIVE_FIELD(ir, tS0, IFT_B3Material_ts0, "ts0");

        // microprestress-sol-theory: read data for inverse desorption isotherm
        IR_GIVE_FIELD(ir, w_h, IFT_B3Material_wh, "w_h");
        IR_GIVE_FIELD(ir, n, IFT_B3Material_ncoeff, "ncoeff");
        IR_GIVE_FIELD(ir, a, IFT_B3Material_a, "a");
        //?????????????????????????? co z vyse uvedenych parametru dat volitelne a jake hodnoty?
    }

    // read shrinkage mode
    int shm = 0;
    IR_GIVE_FIELD(ir, shm, IFT_B3Material_shmode, "shmode");
    this->shMode = ( b3ShModeType ) shm;

    if ( this->shMode == B3_PointShrinkage ) {   // #2 in enumerator
        // read additional shrinkage parameters
        IR_GIVE_FIELD(ir, es0, IFT_B3Material_es0, "es0");
        IR_GIVE_FIELD(ir, r, IFT_B3Material_r, "r");
        IR_GIVE_FIELD(ir, rprime, IFT_B3Material_rprime, "rprime");
        IR_GIVE_FIELD(ir, at, IFT_B3Material_at, "at");
        // read sorption isotherm data
        IR_GIVE_FIELD(ir, w_h, IFT_B3Material_wh, "w_h");
        IR_GIVE_FIELD(ir, n, IFT_B3Material_ncoeff, "ncoeff");
        IR_GIVE_FIELD(ir, a, IFT_B3Material_a, "a");
    } else if ( this->shMode == B3_PointShrinkageMPS ) {   // #3 in enumerator
        // microprestress-sol-theory: read data for shrinkage evaluation
        if ( this->MicroPrestress == 0 ) {
            _error("to use B3_PointShrinkageMPS - MicroPrestress must be = 1");       //or else no external fiels would be found
        }

        kSh = -1;
        initHum = -1;
        finalHum = -1;

        IR_GIVE_OPTIONAL_FIELD(ir, kSh, IFT_B3Material_ksh, "ksh");
        IR_GIVE_OPTIONAL_FIELD(ir, initHum, IFT_B3Material_initialhumidity, "inithum");
        IR_GIVE_OPTIONAL_FIELD(ir, finalHum, IFT_B3Material_finalhumidity, "finalhum");
        // either kSh or initHum and finalHum must be given in input record
        if ( !( ( this->kSh != -1 ) || ( ( initHum != -1 ) && ( finalHum != -1 ) ) ) ) {
            _error("either kSh or initHum and finalHum must be given in input record");
        }

        /*
         * if ( ( ( initHum < 0.2 ) || ( initHum > 0.98 ) || ( finalHum < 0.2 ) || ( finalHum > 0.98 ) ) && ( this->kSh == -1 ) ) {
         *  _error("initital humidity or final humidity out of range (0.2 - 0.98)");
         * }
         */

        if ( this->kSh == -1 ) {
            IR_GIVE_OPTIONAL_FIELD(ir, alpha1, IFT_B3Material_alpha1, "alpha1");       // influence of cement type
            IR_GIVE_OPTIONAL_FIELD(ir, alpha2, IFT_B3Material_alpha2, "alpha2");       // influence of curing type
            this->kSh = alpha1 * alpha2 * ( 1 - pow(finalHum, 3.) ) * ( 0.019 * pow(wc * c, 2.1) * pow(fc, -0.28) + 270 ) * 1.e-6 / fabs(initHum - finalHum);
        }
    } else if ( this->shMode == B3_AverageShrinkage ) {
        IR_GIVE_FIELD(ir, ks, IFT_B3Material_ks, "ks");   // cross-section shape factor
        /*
         * use ks = 1.0 for an infinite slab
         *        = 1.15 for an infinite cylinder
         *        = 1.25 for an infinite square prism
         *        = 1.30 for a sphere
         *        = 1.55 for a cube
         */
        IR_GIVE_FIELD(ir, hum, IFT_B3Material_hum, "hum");         // relative humidity of the environment
        IR_GIVE_FIELD(ir, vs, IFT_B3Material_vs, "vs");         // volume-to-surface ratio (in m???)
        if ( mode == 0 ) {       // default mode - estimate model parameters kt, EpsSinf, q5 from composition
            IR_GIVE_OPTIONAL_FIELD(ir, alpha1, IFT_B3Material_alpha1, "alpha1");             // influence of cement type
            IR_GIVE_OPTIONAL_FIELD(ir, alpha2, IFT_B3Material_alpha2, "alpha2");             // influence of curing type
        } else {         // read model parameters
            IR_GIVE_FIELD(ir, t0, IFT_B3Material_t0, "t0"); // age when drying begins [days]
            IR_GIVE_FIELD(ir, kt, IFT_B3Material_kt, "kt");
            IR_GIVE_FIELD(ir, EpsSinf, IFT_B3Material_EpsSinf, "epssinf");
            IR_GIVE_FIELD(ir, q5, IFT_B3Material_q5, "q5");
        }
    }

    IR_GIVE_FIELD(ir, talpha, IFT_B3Material_talpha, "talpha"); // Macro

    // evaluate the total water content [kg/m^3]
    w = wc * c;
    // estimate the conventional modulus ??? what if fc is not given ???
    E28 = 4733. * sqrt(fc);
    // estimate parameters from composition
    if ( mode == 0 ) {
        this->predictParametersFrom(fc, c, wc, ac, t0, alpha1, alpha2);
    }

    return IRRT_OK;
}

void
B3SolidMaterial :: predictParametersFrom(double fc, double c, double wc, double ac,
                                         double t0, double alpha1, double alpha2)
{
    /*
     * Prediction of model parameters - estimation from concrete composition
     * and strength
     *
     * fc   - 28-day standard cylinder compression strength in MPa
     * c    - cement content of concrete  in kg/m^3
     * wc   - ratio (by weight) of water to cementitious material
     * ac   - ratio (by weight) of aggregate to cement
     * t0   - age when drying begins (in days)
     * alpha1, alpha2   - influence of cement type and curing
     *
     * The prediction of the material parameters of the present model
     * is restricted to Portland cement concretes with the following
     * parameter ranges:
     *
     * 2500 <= fc <= 10000 [psi]   (psi = 6895 Pa)
     * 10   <= c  <= 45    [lb ft-3] (lb ft-3 = 16.03 kg m-3)
     * 0.3  <= wc <= 0.85
     * 2.5  <= ac <= 13.5
     *
     *
     * alpha1 = 1.0  for type I cement;
     *        = 0.85 for type II cement;
     *        = 1.1  for type III cement;
     *
     * alpha2 = 0.75 for steam-cured specimens;
     *        = 1.2  for specimens sealed during curing;
     *        = 1.0  for specimens cured in water or 100% relative humidity.
     *
     */

    // Basic creep parameters

    // q1  = 0.6e6 / E28; // replaced by the formula dependent on fc
    q1  = 126.74271 / sqrt(fc);
    q2  = 185.4 * pow(c, 0.5) * pow(fc, -0.9); // [1/TPa]
    q3  = 0.29 * pow(wc, 4.) * q2;
    q4  = 20.3 * pow(ac, -0.7);

    // Shrinkage parameters

    if ( this->shMode == B3_AverageShrinkage ) {
        kt = 85000 * pow(t0, -0.08) * pow(fc, -0.25); // 85000-> result in [days/m^2] or 8.5-> result in [days/cm^2]
        EpsSinf = alpha1 * alpha2 * ( 1.9e-2 * pow(w, 2.1) * pow(fc, -0.28) + 270. ); // exponent corrected: -0.25 -> -0.28

        // Drying creep parameter
        q5 = 7.57e5 * ( 1. / fc ) * pow(EpsSinf, -0.6);
    }

    char buff [ 1024 ];
    sprintf(buff, "q1=%lf q2=%lf q3=%lf q4=%lf q5=%lf kt=%lf EpsSinf=%lf", q1, q2, q3, q4, q5, kt, EpsSinf);
    OOFEM_LOG_DEBUG("B3mat[%d]: estimated params: %s\n", this->number, buff);
}

void
B3SolidMaterial :: giveThermalDilatationVector(FloatArray &answer,
                                               GaussPoint *gp,  TimeStep *tStep)
//
// returns a FloatArray(6) of initial strain vector
// eps_0 = {exx_0, eyy_0, ezz_0, gyz_0, gxz_0, gxy_0}^T
// caused by unit temperature in direction of
// gp (element) local axes
//
{
    answer.resize(6);
    answer.zero();
    answer.at(1) = ( talpha );
    answer.at(2) = ( talpha );
    answer.at(3) = ( talpha );
}

double
B3SolidMaterial :: giveEModulus(GaussPoint *gp, TimeStep *atTime)
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

    double v, eta;
    double sum = 0.0;
    double t_halfstep;

    t_halfstep = relMatAge + ( atTime->giveTargetTime() - 0.5 * atTime->giveTimeIncrement() ) / timeFactor;
    v = computeSolidifiedVolume(gp, atTime);
    eta = this->computeFlowTermViscosity(gp, atTime);     //evaluated in the middle of the time-step

    if ( EparVal.giveSize() == 0 ) {
        this->updateEparModuli(gp, t_halfstep);
    }

    if ( this->EmoduliMode == 0 ) { //retardation spectrum used
        // first Kelvin units of the Kelvin chain will be computed
        sum = KelvinChainMaterial :: giveEModulus(gp, atTime);
        // the first aging elastic spring must be added
        sum += 1 / EspringVal;
    } else {   //least-squares method used
        sum = KelvinChainMaterial :: giveEModulus(gp, atTime);
    }

    return 1. / ( q1 * 1.e-6 + sum / v  + 0.5 * ( atTime->giveTimeIncrement() / timeFactor ) / eta );
}

void
B3SolidMaterial :: computeCharTimes()
{
    /*
     * This function generates discrete characteristic times
     * according to rules which guarantee a good approximation
     * of the relaxation or creep function by the Dirichlet series
     */

    int mu;
    double Tau1;
    int j;

    if ( this->begOfTimeOfInterest == -1 ) {
        this->begOfTimeOfInterest = 0.1;         //default value
    }

    // perhaps different values should be used for moduli computed by the least-squares method
    // this gives really big mistakes for times near the interest boundary

    //first retardation time found in given by formula 0.3 * begOfTimeOfInterest
    Tau1 = 0.3 * this->begOfTimeOfInterest;

    //last retardation time has to be bigger than 0.5 * endOfTimeOfInterest
    this->endOfTimeOfInterest = RheoChainMaterial :: giveEndOfTimeOfInterest();

    j = 1;
    //while ( 0.5 * this->endOfTimeOfInterest >= Tau1 * pow10(j-1) ) {
    while ( 0.5 * this->endOfTimeOfInterest >= Tau1 * pow( 10.0, ( double )(j - 1) ) ) {
        j++;
    }


    this->nUnits = j;

    this->charTimes.resize(this->nUnits);

    for ( mu = 1; mu <= this->nUnits; mu++ ) {
        charTimes.at(mu) = Tau1 * pow(10., mu - 1);
    }
}


void
B3SolidMaterial :: computeCharCoefficients(FloatArray &answer, GaussPoint *gp, double atTime)
{
    /*
     * If EmoduliMode == 0 then analysis of continuous retardation spectrum is used for
     * computing characteristic coefficients (moduli) of Kelvin chain
     * Else least-squares method is used
     */
    if ( this->EmoduliMode == 0 ) {
        int mu;
        double tau0, tauMu;
        // constant "lambda0" is assumed to be equal to 1 day (typical value)
        // constant "n" is assumed to be equal to 0.1 (typical value)

        // modulus of elasticity of the first unit of Kelvin chain.
        // (aging elastic spring with retardation time = 0)
        tau0 = pow(2 * this->giveCharTime(1) / sqrt(10.0), 0.1);
        EspringVal = 1.e6 / ( q2 * log(1.0 + tau0) - q2 * tau0 / ( 10.0 + 10.0 * tau0 ) );

        // evaluation of moduli of elasticity for the remaining units
        // (aging kelvin units with retardation time tauMu)
        answer.resize(nUnits);
        for ( mu = 1; mu <= this->nUnits; mu++ ) {
            tauMu = pow(2 * this->giveCharTime(mu), 0.1);
            answer.at(mu) = 10.e6 * pow(1 + tauMu, 2) / ( log(10.0) * q2 * tauMu * ( 0.9 + tauMu ) );
            this->charTimes.at(mu) *= 1.35;
        }

        answer.at(nUnits) /= 1.2;   //last unit moduli reduction
    } else {   // moduli computed using the least-squares method
        int i, j, r, rSize;
        double taui, tauj, sum, tti, ttj, sumRhs;
        FloatArray rhs(this->nUnits);
        FloatMatrix A(this->nUnits, this->nUnits);

        const FloatArray &rTimes = this->giveDiscreteTimes();
        rSize = rTimes.giveSize();
        FloatArray discreteComplianceFunctionVal(rSize);

        // compute values of the compliance function at specified times rTimes
        // (can be done directly, since the compliance function is available)

        for ( i = 1; i <= rSize; i++ ) {
            discreteComplianceFunctionVal.at(i) = this->computeNonAgingCreepFunction( gp, rTimes.at(i) );
        }

        // assemble the matrix of the set of linear equations
        // for computing the optimal compliances
        // !!! chartime exponents are assumed to be equal to 1 !!!
        for ( i = 1; i <= this->nUnits; i++ ) {
            taui = this->giveCharTime(i);
            for ( j = 1; j <= this->nUnits; j++ ) {
                tauj = this->giveCharTime(j);
                for ( sum = 0., r = 1; r <= rSize; r++ ) {
                    tti = rTimes.at(r) / taui;
                    ttj = rTimes.at(r) / tauj;
                    sum += ( 1. - exp(-tti) ) * ( 1. - exp(-ttj) );
                }

                A.at(i, j) = sum;
            }

            // assemble rhs
            // !!! chartime exponents are assumed to be equal to 1 !!!
            for ( sumRhs = 0., r = 1; r <= rSize; r++ ) {
                tti = rTimes.at(r) / taui;
                sumRhs += ( 1. - exp(-tti) ) * discreteComplianceFunctionVal.at(r);
            }

            rhs.at(i) = sumRhs;
        }

        // solve the linear system
        A.solveForRhs(rhs, answer);

        // convert compliances into moduli
        for ( i = 1; i <= this->nUnits; i++ ) {
            answer.at(i) = 1.e6 / answer.at(i);
        }
    }
}

double
B3SolidMaterial :: computeNonAgingCreepFunction(GaussPoint *gp, double loadDuration)
// compute the value of the creep function of the non-aging solidifying constituent
// corresponding to the given load duration (in days)

{
    double Phi;     //return value
    double lambda0 = 1.0;     // standard value [day]
    double n = 0.1;     // standard value

    Phi = q2 * log( 1 + pow(loadDuration / lambda0, n) );

    return Phi;
}

//original function from B3Material - useless here
double
B3SolidMaterial :: computeCreepFunction(GaussPoint *gp, double atTime, double ofAge)
{
    // computes the value of creep function at time atTime
    // when load is acting from ofAge
    // WARNING: area returned by crossSection is assumed to be in [m^2].

    double Qf, Z, r, Q, C0, TauSh, St1, St2, H1, H2, Cd;
    double n, m;

    // standard values of exponents
    m = 0.5;
    n = 0.1;

    // basic creep - approximation of the exact B3 model by closed-form functions
    Qf = 1. / ( 0.086 * pow(ofAge, 2. / 9.) + 1.21 * pow(ofAge, 4. / 9.) );
    Z  = pow(ofAge, -m) * log( 1. + pow(atTime - ofAge, n) );
    r  = 1.7 * pow(ofAge, 0.12) + 8.0;
    Q  = Qf * pow( ( 1. + pow( ( Qf / Z ), r ) ), -1. / r );

    C0 = q2 * Q + q3 *log( 1. + pow(atTime - ofAge, n) ) + q4 *log(atTime / ofAge);


    Cd = 0.0;
    if ( this->shMode == B3_AverageShrinkage ) {
        // additional creep due to drying

        TauSh = kt * pow(ks * 2.0 * vs, 2.);
        if ( ( atTime - t0 ) >= 0 ) {
            St1  = tanh( pow( ( atTime - t0 ) / TauSh, 1. / 2. ) );
        } else {
            St1 = 0.0;
        }

        if ( ( ofAge - t0 ) >= 0 ) {
            St2  = tanh( pow( ( ofAge - t0 ) / TauSh, 1. / 2. ) );
        } else {
            St2 = 0.0;
        }

        H1  = 1. - ( 1. - hum ) * St1;
        H2  = 1. - ( 1. - hum ) * St2;
        Cd = q5 * pow( ( exp(-8.0 * H1) - exp(-8.0 * H2) ), 0.5 );
    }

    return 1.e-6 * ( q1 + C0 + Cd );
}


void
B3SolidMaterial :: giveShrinkageStrainVector(FloatArray &answer,
                                             MatResponseForm form,
                                             GaussPoint *gp,
                                             TimeStep *atTime,
                                             ValueModeType mode)
{
    FloatArray prevAnswer;

    if ( this->shMode == B3_NoShrinkage ) {
        if ( form == FullForm ) {
            answer.resize(6);
            answer.zero();
        } else {
            answer.resize( this->giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() ) );
            answer.zero();
        }

        return;
    }

    if ( ( mode != VM_Total ) && ( mode != VM_Incremental ) ) {
        _error("giveShrinkageStrainVector: unsupported mode");
    }

    if ( this->shMode == B3_AverageShrinkage ) {
        this->computeTotalAverageShrinkageStrainVector(answer, form, gp, atTime);

        if ( ( mode == VM_Incremental ) && ( !atTime->isTheFirstStep() ) ) {
            this->computeTotalAverageShrinkageStrainVector( prevAnswer, form, gp, atTime->givePreviousStep() );
            answer.subtract(prevAnswer);
        }
    } else if ( this->shMode == B3_PointShrinkageMPS ) {
        this->computePointShrinkageStrainVectorMPS(answer, form, gp, atTime);
    } else {
        this->computeShrinkageStrainVector(answer, form, gp, atTime, mode);
    }
}

void
B3SolidMaterial :: computeTotalAverageShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
                                                            GaussPoint *gp, TimeStep *atTime)
{
    /*
     * returns average shrinkage strain vector of cross section at drying
     */

    double TauSh, St, kh, help, E607, Et0Tau, EpsShInf, EpsSh;
    double time = relMatAge + atTime->giveTargetTime() / timeFactor;
    int size = 6;
    FloatArray fullAnswer;
    MaterialMode mode = gp->giveMaterialMode();

    if ( ( mode == _3dShell ) || ( mode ==  _3dBeam ) || ( mode == _2dPlate ) || ( mode == _2dBeam ) ) {
        size = 12;
    } else {
        size = 6;
    }

    fullAnswer.resize(size);
    fullAnswer.zero();

    // size dependence
    TauSh = kt * pow(ks * 2.0 * vs, 2.);
    // time curve
    St  = tanh( pow( ( time - t0 ) / TauSh, 1. / 2. ) );
    // humidity dependence
    if ( hum <= 0.98 ) {
        kh = 1. - pow(hum, 3);
    } else if ( hum == 1 ) {
        kh = -0.2;              // swelling in water
    } else {
        // linear interpolation for 0.98 <= h <= 1.
        help = 1. - pow(0.98, 3);
        kh = help + ( -0.2 - help ) / ( 1. - 0.98 ) * ( hum - 0.98 );
    }

    // time dependence of ultimate shrinkage
    E607 = E28 * pow(607 / ( 4. + 0.85 * 607 ), 0.5);
    Et0Tau = E28 * pow( ( t0 + TauSh ) / ( 4. + 0.85 * ( t0 + TauSh ) ), 0.5 );

    EpsShInf = EpsSinf * E607 / Et0Tau;
    // mean shrinkage in the cross section:
    EpsSh = -EpsShInf * kh * St;

    fullAnswer.at(1) = fullAnswer.at(2) = fullAnswer.at(3) = EpsSh * 1.e-6;

    if ( form == FullForm ) {
        answer = fullAnswer;
        return;
    }

    ( ( StructuralCrossSection * ) gp->giveCrossSection() )->
    giveReducedCharacteristicVector(answer, gp, fullAnswer);
}


void
B3SolidMaterial :: computePointShrinkageStrainVectorMPS(FloatArray &answer, MatResponseForm form,
                                                        GaussPoint *gp, TimeStep *atTime)
{
    /* dEpsSh/dt = kSh * dh/dt   (h = humidity)
     * ->> EpsSh = kSh * h_difference
     */
    double humDiff, EpsSh;
    int size;
    FloatArray fullAnswer;
    MaterialMode mode = gp->giveMaterialMode();

    if ( ( mode == _3dShell ) || ( mode ==  _3dBeam ) || ( mode == _2dPlate ) || ( mode == _2dBeam ) ) {
        size = 12;
    } else {
        size = 6;
    }

    humDiff = this->giveHumidityIncrement(gp, atTime);
    EpsSh = humDiff * kSh;

    fullAnswer.resize(size);
    fullAnswer.zero();
    fullAnswer.at(1) = fullAnswer.at(2) = fullAnswer.at(3) = EpsSh;

    if ( form == FullForm ) {
        answer = fullAnswer;
        return;
    }

    ( ( StructuralCrossSection * ) gp->giveCrossSection() )->
    giveReducedCharacteristicVector(answer, gp, fullAnswer);
}


double
B3SolidMaterial :: computeSolidifiedVolume(GaussPoint *gp, TimeStep *atTime)
// compute the relative volume of the solidified material at given age (in days)
{
    double m, lambda0, alpha;
    double v;     //return value
    double atAge;     // (equivalent age)

    // standard values of exponents - empirical constants
    m = 0.5;
    lambda0 = 1;     //[day]
    alpha = q3 / q2;

    atAge = relMatAge + ( atTime->giveTargetTime() - 0.5 * atTime->giveTimeIncrement() ) / timeFactor;
    v = 1 / ( alpha + pow(lambda0 / atAge, m) );

    return v;
}


double
B3SolidMaterial :: computeFlowTermViscosity(GaussPoint *gp, TimeStep *atTime)
//used for evaluation of the flow term viscosity (term containing q4)
{
    double eta, S, tHalfStep;

    if ( this->MicroPrestress == 1 ) {
        S = this->computeMicroPrestress(gp, atTime, 0); //microprestress in the middle of the time-step
        eta = 1.e6 / ( q4 * c0 * S );
    } else if ( this->MicroPrestress == 0 ) {
        tHalfStep = relMatAge + ( atTime->giveTargetTime() - 0.5 * atTime->giveTimeIncrement() ) / timeFactor;
        eta = 1.e6 * tHalfStep / q4;
    } else {
        _error("computeFlowTermViscosity - mode is not supported");
    }

    return eta;
}


void
B3SolidMaterial :: giveEigenStrainVector(FloatArray &answer, MatResponseForm form,
                                         GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
//
// computes the strain due to creep at constant stress during the increment
// (in fact, the INCREMENT of creep strain is computed for mode == VM_Incremental)
//
{
    double v, eta;
    FloatArray help, reducedAnswer, sigma;
    FloatMatrix C;
    KelvinChainMaterialStatus *status = ( KelvinChainMaterialStatus * ) this->giveStatus(gp);

    // !!! chartime exponents are assumed to be equal to 1 !!!
    if ( mode == VM_Incremental ) {
        v = computeSolidifiedVolume(gp, atTime);
        eta = this->computeFlowTermViscosity(gp, atTime);         //evaluated too in the middle of the time-step

        sigma = status->giveStressVector();         //stress vector at the beginning of time-step
        this->giveUnitComplianceMatrix(C, ReducedForm, gp, atTime);
        reducedAnswer.resize( C.giveNumberOfRows() );

        reducedAnswer.beProductOf(C, sigma);
        reducedAnswer.times( atTime->giveTimeIncrement() / ( timeFactor * eta ) );

        //computes non-aging creep component of the Kelvin Chain
        KelvinChainMaterial :: giveEigenStrainVector(help, ReducedForm, gp, atTime, mode);

        help.times(1 / v);
        reducedAnswer.add(help);

        if ( form == ReducedForm ) {
            answer =  reducedAnswer;
            return;
        }

        // expand the strain to full form if requested
        ( ( StructuralCrossSection * ) gp->giveCrossSection() )->
        giveFullCharacteristicVector(answer, gp, reducedAnswer);
    } else {
        /* error - total mode not implemented yet */
        _error("giveEigenStrainVector - mode is not supported");
    }
}


void
B3SolidMaterial :: computeShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
                                                GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
{
    // free shrinkage at material point, requires staggered analaysis
    // additional material parameters required:
    //  es0     - final shrinkage at material point
    //  r       - coefficient
    //  rprime  - coefficient
    //  at      - coeff relating stress-induced thermal strain and shrinkage
    double sv, sn, et0, et, wrate = 0.0, trate = 0.0, h1;
    double time = relMatAge + atTime->giveTargetTime() / timeFactor;
    int i, err, tflag = 0, wflag = 0;
    KelvinChainMaterialStatus *status = ( KelvinChainMaterialStatus * ) this->giveStatus(gp);
    int size = 6;
    FloatArray fullAnswer;
    MaterialMode mmode = gp->giveMaterialMode();

    if ( ( mmode == _3dShell ) || ( mmode ==  _3dBeam ) || ( mmode == _2dPlate ) || ( mmode == _2dBeam ) ) {
        size = 12;
    } else {
        size = 6;
    }

    fullAnswer.resize(size);
    fullAnswer.zero();



    /* ask for humidity and temperature from external sources, if provided */
    FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();
    Field *tf;
    FloatArray gcoords, et2, ei2, stressVector, fullStressVector;

    if ( ( tf = fm->giveField(FT_Temperature) ) ) {
        // temperature field registered
        gp->giveElement()->computeGlobalCoordinates( gcoords, * gp->giveCoordinates() );
        if ( ( err = tf->evaluateAt(et2, gcoords, VM_Incremental, atTime) ) ) {
            _error2("computeShrinkageStrainVector: tf->evaluateAt failed, error value %d", err);
        }

        trate = et2.at(1);
        tflag = 1;
    }

    if ( ( tf = fm->giveField(FT_HumidityConcentration) ) ) {
        // temperature field registered
        gp->giveElement()->computeGlobalCoordinates( gcoords, * gp->giveCoordinates() );
        if ( ( err = tf->evaluateAt(et2, gcoords, VM_Total, atTime) ) ) {
            _error2("computeShrinkageStrainVector: tf->evaluateAt failed, error value %d", err);
        }

        if ( ( err = tf->evaluateAt(ei2, gcoords, VM_Incremental, atTime) ) ) {
            _error2("computeShrinkageStrainVector: tf->evaluateAt failed, error value %d", err);
        }

        // convert water mass to relative humidity
        wrate = this->inverse_sorption_isotherm( et2.at(1) ) - this->inverse_sorption_isotherm( et2.at(1) - ei2.at(1) );
        wflag = 1;
    }

    if ( ( tflag == 0 ) || ( wflag == 0 ) ) {
        _error("computeTotalShrinkageStrainVector: external fields not found");
    }

    if ( status->giveStressVector().giveSize() ) {
        stressVector      = status->giveStressVector();
    } else {
        stressVector.resize( this->giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() ) );
        stressVector.zero();
    }

    giveFullCharacteristicVector(fullStressVector, gp, stressVector);
    // compute volumetric stress
    sv = 0.0;
    for ( i = 1; i <= 3; i++ ) {
        sv += stressVector.at(i);
    }

    et = 1. / this->computeCreepFunction(gp, time + 0.01, time);
    et0 = 1. / this->computeCreepFunction(gp, t0 + 0.01, t0);

    h1 = es0 * ( et0 / et );
    sn = sgn(wrate + at * trate);
    // compute increment of shrinkage strain
    fullAnswer.at(1) = h1 * ( 1.0 + sn * ( r * fullStressVector.at(1) + rprime * sv ) ) * ( wrate + at * trate );
    fullAnswer.at(2) = h1 * ( 1.0 + sn * ( r * fullStressVector.at(2) + rprime * sv ) ) * ( wrate + at * trate );
    fullAnswer.at(3) = h1 * ( 1.0 + sn * ( r * fullStressVector.at(3) + rprime * sv ) ) * ( wrate + at * trate );

    if ( mode == VM_Incremental ) {
        if ( form == FullForm ) {
            answer = fullAnswer;
            return;
        }

        ( ( StructuralCrossSection * ) gp->giveCrossSection() )->
        giveReducedCharacteristicVector(answer, gp, fullAnswer);
        return;
    } else { // total values required
        FloatArray ssv, fssv;
        if ( status->giveShrinkageStrainVector()->giveSize() == 0 ) {
            ssv.resize( this->giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() ) );
            ssv.zero();
        } else {
            ssv = * status->giveShrinkageStrainVector();
        }

        this->giveFullCharacteristicVector(fssv, gp, ssv);
        // add increment to total values
        fullAnswer.add(fssv);

        if ( form == FullForm ) {
            answer = fullAnswer;
            return;
        }

        ( ( StructuralCrossSection * ) gp->giveCrossSection() )->
        giveReducedCharacteristicVector(answer, gp, fullAnswer);
        return;
    }
}


double
B3SolidMaterial :: inverse_sorption_isotherm(double w)
// Function calculates relative humidity from water content (inverse relation form sorption isotherm).
// Relative humidity (phi) is from range 0.2 - 0.98 !!!
// sorption isotherm by C. R. Pedersen (1990), Combined heat and moisture transfer in building constructions,
//                      PhD-thesis, Technical University of Denmark, Lyngby.
// w (kg/kg) ... water content
// phi ... relative humidity
// w_h, n, a ... constants obtained from experiments
{
    double phi;

    // relative humidity
    phi = exp( a * ( 1.0 - pow( ( w_h / w ), ( n ) ) ) );

    /*
     * if ( ( phi < 0.2 ) || ( phi > 0.98 ) ) {
     *  _error3("inverse_sorption_isotherm : Relative humidity h = %e (w=%e) is out of range", phi, w);
     * }
     */

    return phi;
}


double
B3SolidMaterial :: computeMicroPrestress(GaussPoint *gp, TimeStep *atTime, int option)
{
    //return 1.;//!!!!!!!!!!!!!!!!!!!!!!!!!!
    /* numerically solves non-linear differential equation in form:
     * dS/dt + c0*S^2 = -c1/h * dh/dt
     * where "S" is the microprestress, "c0" and "c1" constants and "h" humidity
     *
     * If option == 0, microprestress is evaluated in the middle of the time step (used for stiffnesses)
     * If option == 1, MPS is evaluated at the end of the time step. (Used for updating)
     */
    //UNCOMMENT THIS SECTION TO USE GENERALIZED TRAPEZOIDAL RULE
    /* new value of S - "Stemp" is computed using the generalized trapezoidal rule
     * with weights "1-alpha" and "alpha".
     */

    /*double S, Stemp; // old and new microprestress
     * double humOld, humNew; // previous and new value of humidity
     * double a, b, c, D; // auxiliary parameters used for solving quadratic equation
     *
     * double alpha = 0.5; // 0 for forward Euler method, 1 for backward, 0.5 for trapezoidal rule
     * double deltaT;
     * deltaT = atTime->giveTimeIncrement()/timeFactor;
     *
     * if(atTime->isTheFirstStep()){
     *      // if the time step is the first one, ask for an initial microprestress
     *      S = this->giveInitMicroPrestress();
     * } else {
     *      B3SolidMaterialStatus *status = ( B3SolidMaterialStatus * ) this->giveStatus(gp);
     *      S = status->giveMPS();
     * }
     *
     * humOld = this->giveHumidity(gp, atTime) - this->giveHumidityIncrement(gp, atTime);
     *
     * if (option == 1) {
     *      humNew = this->giveHumidity(gp, atTime);
     * } else if (option == 0) {
     *      deltaT /= 2;
     *      humNew = humOld + 0.5*(this->giveHumidityIncrement(gp, atTime)); //linearly approximated
     * } else {
     *      _error("computeMicroPrestress: invalid option parameter");
     * }
     *
     * a = c0 * alpha * deltaT;
     * b = 1.;
     * c = -S + c1* (1-alpha)*(humNew/humOld -1) + c1*alpha*(1- humOld/humNew) + deltaT*pow(S, 2.)*(1-alpha)*c0;
     *
     * D = pow(b, 2.) - 4 * a * c; // discriminant
     * Stemp = ( -b + sqrt(D) )/ (2 * a); // only positive root of quadr. equation is solved
     *
     * return Stemp;*/

    // THIS SECTION USES PARTICULAR SOLUTION OF THE GOVERNING DIFFERENTIAL EQUATION
    double S, Stemp;     // old and new microprestress
    double humOld, humNew;     // previous and new value of humidity
    double A;     // auxiliary variable
    double RHS;     // constant RHS of the diff equation
    double dHdt;     // first time difference of the humidity
    double deltaT;     // length of the time step

    deltaT = atTime->giveTimeIncrement() / timeFactor;
    if ( atTime->isTheFirstStep() ) {
        // if the time step is the first one, ask for an initial microprestress
        S = this->giveInitMicroPrestress();
    } else {
        B3SolidMaterialStatus *status =
            ( B3SolidMaterialStatus * ) this->giveStatus(gp);
        S = status->giveMPS();
    }

    humOld = this->giveHumidity(gp, atTime) - this->giveHumidityIncrement(gp,
                                                                          atTime);
    if ( option == 1 ) {  // hum in the end of the time step
        humNew = this->giveHumidity(gp, atTime);
    } else if ( option == 0 ) {  // hum in the middle of the time step
        deltaT /= 2;
        humNew = humOld + 0.5 * ( this->giveHumidityIncrement(gp, atTime) );       //linearly approximated
    } else {
        _error("computeMicroPrestress: invalid option parameter");
    }

    //following section is used if humidity remains constant
    if ( ( humNew - humOld ) == 0. ) {
        Stemp = 1 / ( 1 / S + c0 * deltaT );
    } else {     // the following section is used only if there is a change in humidity
        dHdt = ( humNew - humOld ) / deltaT;
        RHS = fabs(c1 * dHdt / humNew);
        A = sqrt(RHS / c0);
        Stemp = A * ( 1 - 2 * ( A - S ) / ( ( A - S ) + ( A + S ) * exp( 2 * deltaT * sqrt(RHS * c0) ) ) );
    }

    return Stemp;
}


double
B3SolidMaterial :: giveInitMicroPrestress(void)
{
    double S0;
    S0 = 1 / ( c0 * tS0 );
    return S0;
}


double
B3SolidMaterial :: giveHumidity(GaussPoint *gp, TimeStep *atTime) //computes humidity at given TimeStep
{
    double humidity;
    int err, wflag = 0;

    /* ask for humidity from external sources, if provided */
    FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();
    Field *tf;
    FloatArray gcoords, et2;

    if ( ( tf = fm->giveField(FT_HumidityConcentration) ) ) {
        // humidity field registered
        gp->giveElement()->computeGlobalCoordinates( gcoords, * gp->giveCoordinates() );
        if ( ( err = tf->evaluateAt(et2, gcoords, VM_Total, atTime) ) ) {
            _error2("giveHumidity: tf->evaluateAt failed, error value %d", err);
        }

        // convert water mass to relative humidity
        humidity = this->inverse_sorption_isotherm( et2.at(1) );
        wflag = 1;
    }

    if ( wflag == 0 ) {
        _error("giveHumidity: external fields not found");
    }

    return humidity;
}


double
B3SolidMaterial :: giveHumidityIncrement(GaussPoint *gp, TimeStep *atTime) //computes humidity increment at given TimeStep
{
    double humIncrement;
    int err, wflag = 0;

    /* ask for humidity from external sources, if provided */
    FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();
    Field *tf;
    FloatArray gcoords, et2, ei2;

    if ( ( tf = fm->giveField(FT_HumidityConcentration) ) ) {
        // humidity field registered
        gp->giveElement()->computeGlobalCoordinates( gcoords, * gp->giveCoordinates() );
        if ( ( err = tf->evaluateAt(et2, gcoords, VM_Total, atTime) ) ) {
            _error2("giveHumidityIncrement: tf->evaluateAt failed, error value %d", err);
        }

        if ( ( err = tf->evaluateAt(ei2, gcoords, VM_Incremental, atTime) ) ) {
            _error2("giveHumidityIncrement: tf->evaluateAt failed, error value %d", err);
        }

        // convert water mass to relative humidity
        humIncrement = this->inverse_sorption_isotherm( et2.at(1) ) - this->inverse_sorption_isotherm( et2.at(1) - ei2.at(1) );
        wflag = 1;
    }

    if ( wflag == 0 ) {
        _error("giveHumidityIncrement: external fields not found");
    }

    return humIncrement;
}


MaterialStatus *
B3SolidMaterial :: CreateStatus(GaussPoint *gp) const
/*
 * creates a new material status corresponding to this class
 */
{
    return new B3SolidMaterialStatus(1, this->giveDomain(), gp, nUnits);
}


void
B3SolidMaterial :: updateYourself(GaussPoint *gp, TimeStep *tNow)
{
    B3SolidMaterialStatus *status = ( B3SolidMaterialStatus * ) this->giveStatus(
        gp);
    if ( this->MicroPrestress == 1 ) {
        status->setMPS( this->computeMicroPrestress(gp, tNow, 1) );
    }

    // now we call Kelvin to update itself
    // at the end of updating it will call material status to be updated.
    KelvinChainMaterial :: updateYourself(gp, tNow);
}


/****************************************************************************************/
/**********     B3SolidMaterialStatus - HUMIDITY ****************************************/


B3SolidMaterialStatus :: B3SolidMaterialStatus(int n, Domain *d, GaussPoint *g, int nunits) :
    KelvinChainMaterialStatus(n, d, g, nunits) {}

void
B3SolidMaterialStatus :: updateYourself(TimeStep *tStep)
{
    //sem se vlozi updatovani microprestressu
    microprestress_old = microprestress_new;
    microprestress_new = 0.;

    KelvinChainMaterialStatus :: updateYourself(tStep);
}


contextIOResultType
B3SolidMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
//
{
    contextIOResultType iores;

    if ( stream == NULL ) {
        _error("saveContext : can't write into NULL stream");
    }

    if ( ( iores = KelvinChainMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write microprestress value
    if ( !stream->write(& microprestress_old, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

contextIOResultType
B3SolidMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restore the state variables from a stream
//
{
    contextIOResultType iores;
    if ( ( iores = KelvinChainMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // write microprestress value
    if ( !stream->read(& microprestress_old, 1) ) {
        return CIO_IOERR;
    }

    return CIO_OK;
}
} // end namespace oofem
