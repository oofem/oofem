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


#ifndef __MAKEDEPEND
 #include <math.h>
#endif
#include "mps.h"
#include "mathfem.h"
#include "gausspnt.h"
#include "structuralcrosssection.h"
#include "timestep.h"
#include "contextioerr.h"
#include "datastream.h"

#include "rankinemat.h"

namespace oofem {
/****************************************************************************************/
/**************     MPSMaterialStatus     ***********************************************/


MPSMaterialStatus :: MPSMaterialStatus(int n, Domain *d, GaussPoint *g, int nunits) :
    KelvinChainSolidMaterialStatus(n, d, g, nunits)
{
    hum = -1.;
    hum_increment = -1.;
    T = -1.;
    T_increment = -1.;
}

void
MPSMaterialStatus :: updateYourself(TimeStep *tStep)
{
    equivalentTime = equivalentTimeTemp;
    equivalentTimeTemp = -1.;

    flowTermViscosity = flowTermViscosityTemp;
    flowTermViscosityTemp = -1.;

    KelvinChainSolidMaterialStatus :: updateYourself(tStep);
}

contextIOResultType
MPSMaterialStatus :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
//
{
    contextIOResultType iores;

    if ( stream == NULL ) {
        _error("saveContext : can't write into NULL stream");
    }

    if ( ( iores = KelvinChainSolidMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    if ( !stream->write(& equivalentTime, 1) ) { // write equivalent time
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& flowTermViscosity, 1) ) { // write viscosity of the aging dashpot
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

contextIOResultType
MPSMaterialStatus :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restore the state variables from a stream
//
{
    contextIOResultType iores;
    if ( ( iores = KelvinChainSolidMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream->read(& equivalentTime, 1) ) { // restore equivalentTime
        return CIO_IOERR;
    }

    if ( !stream->read(& flowTermViscosity, 1) ) {  // restore equivalentTime
        return CIO_IOERR;
    }

    return CIO_OK;
}

//   **********************************************************************************************
//   ***                         CLASS MPS Material                                             ***
//   ***  (concrete creep and shrinkage described by the microprestress solidification theory)  ***
//   **********************************************************************************************


IRResultType
MPSMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    double fc, c, wc, ac;

    /* scaling factor transforming PREDICTED stiffnesses q1, q2, q3 and q4 to desired
     *  default value is 1.0 = no change
     *  e.g. if the stiffness should be in MPa, then stiffnessFactor = 1.e6 */
    double stiffnessFactor;

    KelvinChainSolidMaterial :: initializeFrom(ir);

    // checking timefactor - for MPS material must be equal to 1.
    double tf;
    IR_GIVE_FIELD(ir, tf, IFT_RheoChainMaterial_timefactor, "timefactor");
    if ( tf != 1. ) {
        _error("initializeFrom: for MPS material timefactor must be equal to 1.");
    }

    int mode = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, mode, IFT_MPSMaterial_mode, "mode");

    if ( mode == 0 ) { // default - estimate model parameters q1,..,q4 from composition
        IR_GIVE_FIELD(ir, fc, IFT_MPSMaterial_fc, "fc"); // 28-day standard cylinder compression strength [MPa]
        IR_GIVE_FIELD(ir,  c, IFT_MPSMaterial_cc, "cc"); // cement content of concrete [kg/m^3]
        IR_GIVE_FIELD(ir, wc, IFT_MPSMaterial_wc, "w/c"); // ratio (by weight) of water to cementitious material
        IR_GIVE_FIELD(ir, ac, IFT_MPSMaterial_ac, "a/c"); // ratio (by weight) of aggregate to cement

        stiffnessFactor = 1.;
        IR_GIVE_OPTIONAL_FIELD(ir, stiffnessFactor, IFT_MPSMaterial_stiffnessfactor, "stiffnessfactor"); // ratio (by weight) of aggregate to cement

        this->predictParametersFrom(fc, c, wc, ac, stiffnessFactor);

    } else { // read model parameters for creep
        IR_GIVE_FIELD(ir, q1, IFT_MPSMaterial_q1, "q1");
        IR_GIVE_FIELD(ir, q2, IFT_MPSMaterial_q2, "q2");
        IR_GIVE_FIELD(ir, q3, IFT_MPSMaterial_q3, "q3");
        IR_GIVE_FIELD(ir, q4, IFT_MPSMaterial_q4, "q4");
    }

    IR_GIVE_FIELD(ir, talpha, IFT_MPSMaterial_talpha, "talpha"); // Macro
    IR_GIVE_FIELD(ir, lambda0, IFT_MPSMaterial_lambda0, "lambda0"); // Macro

    int type = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, type, IFT_MPSMaterial_coupledanalysistype, "coupledanalysistype");
    if ( type >= 2 ) {
        _error("initializeFrom: CoupledAnalysisType must be equal to 0 or 1");
    }

    this->CoupledAnalysis = ( coupledAnalysisType ) type;
    this->kSh = 0.;

    if ( this->CoupledAnalysis == MPS ) {
        // muS- the only parameter necessary for evaluation of the flow-term viscosity
        // muS = c0 * c1 * q4 [1/(Pa*s)]
        IR_GIVE_FIELD(ir, muS, IFT_MPSMaterial_mus, "mus");
        IR_GIVE_FIELD(ir, kappaT, IFT_MPSMaterial_kappat, "kappat"); //[-] replaces ln(h) in diff equation for MPS
        this->ct = 0.;
        IR_GIVE_OPTIONAL_FIELD(ir, ct, IFT_MPSMaterial_ct, "ct");

        // age when drying or thermal changes begin [days] - necessary to determine initial viscosity
        IR_GIVE_FIELD(ir, t0, IFT_MPSMaterial_t0, "t0");
        IR_GIVE_OPTIONAL_FIELD(ir, kSh, IFT_MPSMaterial_ksh, "ksh");

        // Parameters for desorption isotherm
        IR_GIVE_FIELD(ir, w_h, IFT_MPSMaterial_wh, "w_h");
        IR_GIVE_FIELD(ir, n, IFT_MPSMaterial_ncoeff, "ncoeff");
        IR_GIVE_FIELD(ir, a, IFT_MPSMaterial_a, "a");

        this->roomTemperature = 298.; // reference room temperature for MPS algorithm [K]

        // ACTIVATION ENERGIES
        this->QEtoR = 2700.; // [K] value according to Bazant 1995
        IR_GIVE_OPTIONAL_FIELD(ir, QEtoR, IFT_MPSMaterial_qetor, "qetor");
        this->QRtoR = 5000.; // [K] value according to Bazant 1995
        IR_GIVE_OPTIONAL_FIELD(ir, QRtoR, IFT_MPSMaterial_qrtor, "qrtor");
        this->QStoR = 3000.; // [K] value according to Bazant 1995
        IR_GIVE_OPTIONAL_FIELD(ir, QStoR, IFT_MPSMaterial_qstor, "qstor");

        this->alphaE = 10.; // according to Bazant 1995
        IR_GIVE_OPTIONAL_FIELD(ir, alphaE, IFT_MPSMaterial_alphae, "alphae");
        this->alphaR = 0.1; // according to Bazant 1995
        IR_GIVE_OPTIONAL_FIELD(ir, alphaR, IFT_MPSMaterial_alphar, "alphar");
        this->alphaS = 0.1; // according to Bazant 1995
        IR_GIVE_OPTIONAL_FIELD(ir, alphaS, IFT_MPSMaterial_alphas, "alphas");
    }

    return IRRT_OK;
}

void
MPSMaterial :: updateYourself(GaussPoint *gp, TimeStep *atTime)
{
    MPSMaterialStatus *status = ( MPSMaterialStatus * ) this->giveStatus(gp);

    if ( this->CoupledAnalysis == MPS ) {
        status->setEquivalentTime( this->computeEquivalentTime(gp, atTime, 1) );
        // set humidity and temperature to zero
        status->setHum(-1.);
        status->setHumIncrement(-1.);
        status->setT(-1.);
        status->setTIncrement(-1.);
    }

    // now we call Solidifying Kelvin Chain to update itself
    // at the end of updating it will call material status to be updated.
    KelvinChainSolidMaterial :: updateYourself(gp, atTime);
}

void
MPSMaterial :: giveThermalDilatationVector(FloatArray &answer,
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

void
MPSMaterial :: giveShrinkageStrainVector(FloatArray &answer,
                                         MatResponseForm form,
                                         GaussPoint *gp,
                                         TimeStep *atTime,
                                         ValueModeType mode)
{
    if ( ( mode != VM_Total ) && ( mode != VM_Incremental ) ) {
        _error("giveShrinkageStrainVector: unsupported mode");
    }


    if ( this->kSh == 0. ) {
        if ( form == FullForm ) {
            answer.resize(6);
            answer.zero();
        } else {
            answer.resize( this->giveSizeOfReducedStressStrainVector( gp->giveMaterialMode() ) );
            answer.zero();
        }
    } else {
        this->computePointShrinkageStrainVector(answer, form, gp, atTime);
    }
}

MaterialStatus *
MPSMaterial :: CreateStatus(GaussPoint *gp) const
/*
 * creates a new material status corresponding to this class
 */
{
    return new MPSMaterialStatus(1, this->giveDomain(), gp, nUnits);
}


void
MPSMaterial :: predictParametersFrom(double fc, double c, double wc, double ac, double stiffnessFactor)
{
    /*
     * Prediction of model parameters - estimation from concrete composition
     * and strength
     *
     * fc   - 28-day standard cylinder compression strength in MPa
     * c    - cement content of concrete  in kg/m^3
     * wc   - ratio (by weight) of water to cementitious material
     * ac   - ratio (by weight) of aggregate to cement
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
     */

    // Basic creep parameters
    q1  = 1.e-12 * stiffnessFactor * 126.74271 / ( sqrt(fc) );
    q2  = 1.e-12 * stiffnessFactor * 185.4 * __OOFEM_POW(c, 0.5) * __OOFEM_POW(fc, -0.9);
    q3  = 1.e-6 *  stiffnessFactor * 0.29 * __OOFEM_POW(wc, 4.) * q2;
    q4  = 1.e-12 * stiffnessFactor * 20.3 * __OOFEM_POW(ac, -0.7);

    char buff [ 1024 ];
    sprintf(buff, "q1=%lf q2=%lf q3=%lf q4=%lf", q1, q2, q3, q4);
    OOFEM_LOG_DEBUG("MPS[%d]: estimated params: %s\n", this->number, buff);
}

void
MPSMaterial :: computeCharTimes()
{
    /*
     * This function generates discrete characteristic times
     * according to rules which guarantee a good approximation
     * of the relaxation or creep function by the Dirichlet series
     */

    int mu;
    double Tau1;
    int j;


    if ( this->begOfTimeOfInterest == -1. ) {
        this->begOfTimeOfInterest = 0.01 * lambda0;
    }

    if ( this->endOfTimeOfInterest == -1. ) {
        this->endOfTimeOfInterest = 10000. * lambda0;
    }

    // perhaps different values should be used for moduli computed by the least-squares method
    // this gives really big mistakes for times near the interest boundary

    //first retardation time found in given by formula 0.3 * begOfTimeOfInterest
    Tau1 = 0.3 * this->begOfTimeOfInterest;

    //last retardation time has to be bigger than 0.5 * endOfTimeOfInterest
    this->endOfTimeOfInterest = RheoChainMaterial :: giveEndOfTimeOfInterest();

    j = 1;
    while ( 0.5 * this->endOfTimeOfInterest >= Tau1 * pow( 10.0, ( double )(j - 1) ) ) {
        j++;
    }

    this->nUnits = j;

    this->charTimes.resize(this->nUnits);

    for ( mu = 1; mu <= this->nUnits; mu++ ) {
        charTimes.at(mu) = Tau1 * __OOFEM_POW(10., mu - 1);
    }
}


void
MPSMaterial :: computeCharCoefficients(FloatArray &answer, GaussPoint *gp, double atTime)
{
    int mu;
    double tau0, tauMu;

    // constant "n" is assumed to be equal to 0.1 (typical value)

    // modulus of elasticity of the first unit of Kelvin chain.
    // (aging elastic spring with retardation time = 0)
    double lambda0ToPowN = __OOFEM_POW(lambda0, 0.1);
    tau0 = __OOFEM_POW(2 * this->giveCharTime(1) / sqrt(10.0), 0.1);
    EspringVal = 1. / ( q2 * log(1.0 + tau0 / lambda0ToPowN) - q2 * tau0 / ( 10.0 * lambda0ToPowN + 10.0 * tau0) );

    // evaluation of moduli of elasticity for the remaining units
    // (Solidifying kelvin units with retardation times tauMu)
    answer.resize(nUnits);
    for ( mu = 1; mu <= this->nUnits; mu++ ) {
        tauMu = __OOFEM_POW(2 * this->giveCharTime(mu), 0.1);
        answer.at(mu) = 10. * __OOFEM_POW(1 + tauMu / lambda0ToPowN, 2) / ( log(10.0) * q2 * ( tauMu / lambda0ToPowN ) * ( 0.9 + tauMu / lambda0ToPowN) );
        this->charTimes.at(mu) *= 1.35;
    }

    answer.at(nUnits) /= 1.2;   // modulus of the last unit is reduced
}


double
MPSMaterial :: giveEModulus(GaussPoint *gp, TimeStep *atTime)
{
    /*
     * This function returns the incremental modulus for the given time increment.
     * The modulus may also depend on the specimen geometry (gp - dependence).
     *
     * It is stored as "Einc" for further expected requests from other gaussPoints that correspond to the same material.
     *
     * Note: time -1 refers to the previous time.
     */

    double sum = 0.0;
    double v;

    double Cf; // incremental viscous flow compliance
    double eta, dt;
    double dEtaR, etaR, L;


    if ( EparVal.giveSize() == 0 ) {
        this->updateEparModuli(gp, 0.); // stiffnesses are time independent (evaluated at time t = 0.)
    }

    // contribution of the solidifying Kelving chain
    sum = KelvinChainSolidMaterial :: giveEModulus(gp, atTime);

    v = computeSolidifiedVolume(gp, atTime);

    dt = atTime->giveTimeIncrement();
    eta = this->computeFlowTermViscosity(gp, atTime);

    //incremental viscous flow compliance

    if ( this->CoupledAnalysis == Basic ) {
        Cf = 0.5 * ( dt ) / eta;
    } else if ( this->CoupledAnalysis == MPS ) {
        MPSMaterialStatus *status = ( MPSMaterialStatus * ) this->giveStatus(gp);

        // TRAPEZOIDAL INTEGRATION RULE
        if ( atTime->isTheFirstStep() ) {
            etaR = this->giveInitViscosity(atTime) /  this->computePsiR(gp, atTime, 0.);
        } else {
            etaR = status->giveFlowTermViscosity() /  this->computePsiR(gp, atTime, 0.);
        }

        dEtaR =  eta /  this->computePsiR(gp, atTime, 1.) - etaR;
        if (  fabs(dEtaR) > 1.e-4 * etaR ) {
            L = log(1 + dEtaR / etaR);
            Cf = dt * ( 1. - etaR * L / dEtaR ) / dEtaR;
        } else {
            Cf = dt * ( 0.5 - dEtaR / ( 3 * etaR ) ) / etaR;
        }

        // TRAPEZOIDAL INTEGRATION RULE

        // MIDPOINT INTEGRATION RULE
        // if ( atTime->isTheFirstStep() ) {
        //  Cf = dt * this->computePsiR(gp, atTime, 2.) / (eta + this->giveInitViscosity(atTime) );
        // } else {
        //  Cf = dt * this->computePsiR(gp, atTime, 2.) / (eta + status->giveFlowTermViscosity() );
        // }
        // TRAPEZOIDAL INTEGRATION RULE
    } else {
        _error("giveEModulus - mode is not supported");
    }

    Einc = 1 / ( q1 + 1. / ( EspringVal * v ) + sum  +  Cf );

    return Einc;
}


double
MPSMaterial :: computeSolidifiedVolume(GaussPoint *gp, TimeStep *atTime)
// compute the relative volume of the solidified material at given age (in days)
{
    double m, alpha;
    double v;     //return value
    double atAge;     // (equivalent age)

    // standard values of exponents - empirical constants
    m = 0.5;

    alpha = q3 / q2;

    if ( this->CoupledAnalysis == Basic ) {
        atAge = relMatAge + ( atTime->giveTargetTime() - 0.5 * atTime->giveTimeIncrement() );
    } else {
        atAge = computeEquivalentTime(gp, atTime, 0);
    }

    v = 1 / ( alpha + __OOFEM_POW(lambda0 / atAge, m) );

    return v;
}


double
MPSMaterial :: computeBetaMu(GaussPoint *gp, TimeStep *atTime, double Mu)
{
    double betaMu;
    double deltaT;
    double tauMu;


    deltaT = atTime->giveTimeIncrement();

    if ( this->CoupledAnalysis == MPS ) {
        deltaT *=  0.5 * ( this->computePsiR(gp, atTime, 0.) + this->computePsiR(gp, atTime, 1.) );
    }

    tauMu = this->giveCharTime(Mu);

    if ( deltaT / tauMu > 30 ) {
        betaMu = 0;
    } else {
        betaMu = exp(-( deltaT ) / tauMu);
    }

    return betaMu;
}

double
MPSMaterial :: computeLambdaMu(GaussPoint *gp, TimeStep *atTime, double Mu)
{
    double lambdaMu;
    double deltaT;
    double tauMu;

    deltaT = atTime->giveTimeIncrement();

    if ( this->CoupledAnalysis == MPS ) {
        deltaT *=  0.5 * ( this->computePsiR(gp, atTime, 0) + this->computePsiR(gp, atTime, 1) );
    }

    tauMu = this->giveCharTime(Mu);

    if ( deltaT / tauMu < 1.e-5 ) {
        lambdaMu = 1 - 0.5 * ( deltaT / tauMu ) + 1 / 6 * ( __OOFEM_POW(deltaT / tauMu, 2) ) - 1 / 24 * ( __OOFEM_POW(deltaT / tauMu, 3) );
    } else if ( deltaT / tauMu > 30 ) {
        lambdaMu = tauMu / deltaT;
    } else {
        lambdaMu = ( 1.0 -  exp(-( deltaT ) / tauMu) ) * tauMu / deltaT;
    }

    return lambdaMu;
}

double
MPSMaterial :: computeFlowTermViscosity(GaussPoint *gp, TimeStep *atTime)
{
    double eta, tHalfStep;

    double prevEta, PsiS, A, B, e, dt;
    double T_new, T_old, H_new, H_old;
    double reductFactor;

    if ( this->CoupledAnalysis == Basic ) {
        tHalfStep = relMatAge + ( atTime->giveTargetTime() - 0.5 * atTime->giveTimeIncrement() );
        eta = tHalfStep / q4;
    } else if ( this->CoupledAnalysis == MPS ) {
        MPSMaterialStatus *status = ( MPSMaterialStatus * ) this->giveStatus(gp);

        // check whether this viscosity has been already computed
        if ( status->giveFlowTermViscosityTemp() != -1.  &&  !( atTime->isTheFirstStep() ) ) {
            return status->giveFlowTermViscosityTemp();

            // no, the viscosity needs to be evaluated now
        } else {
            if ( atTime->isTheFirstStep() ) { // if the time step is the first one, ask for an initial value of viscosity
                prevEta = this->giveInitViscosity(atTime);
            } else {
                // asks for the value of viscosity from the end of the last time-step
                prevEta = status->giveFlowTermViscosity();
            }

            dt = atTime->giveTimeIncrement();
            // evaluate auxiliary factors A and B
            T_new = this->giveTemperature(gp, atTime, 1);
            T_old = this->giveTemperature(gp, atTime, 0);
            H_new = this->giveHumidity(gp, atTime, 1);
            H_old = this->giveHumidity(gp, atTime, 0);

            PsiS = this->computePsiS(gp, atTime); // evaluated in the middle of the time step

            // original version
            //A = sqrt( muS * fabs( T_new * log(H_new) - T_old * log(H_old) ) / ( dt * this->roomTemperature ) );

            if ( this->ct == 0. ) {
                reductFactor = 1.;
            } else if ( ( status->giveTmax() - T_new < 0. ) || atTime->isTheFirstStep() ) {
                status->setTmax(T_new);
                reductFactor = 1.;
            } else {
                reductFactor = exp( - this->ct * fabs (T_new - status->giveTmax() ) );
            }

            A = sqrt( muS * ( kappaT * reductFactor * fabs(T_new - T_old) + 0.5*(T_new + T_old) * fabs(log(H_new) - log(H_old)) ) / ( dt * this->roomTemperature ) );
            B = sqrt(PsiS / this->q4);

            if ( ( A * B * dt ) > 1.e-6 ) {
                e = exp(-2 * A * B * dt);
                eta = ( B / A ) * ( B * ( 1. - e ) + A * prevEta * ( 1. + e ) ) / ( B * ( 1. + e ) + A * prevEta * ( 1. - e ) );
            } else {
                eta = ( prevEta + B * B * dt ) / ( 1. + A * A * prevEta * dt );
            }


            // and now we store into the material status new value of viscosity
            status->setFlowTermViscosityTemp(eta);
        }
    } else {
        _error("computeFlowTermViscosity - mode is not supported");
    }


    return eta;
}

// returns initial value of the flow term viscosity
double
MPSMaterial :: giveInitViscosity(TimeStep *atTime)
{
    if ( ( t0 - atTime->giveTimeIncrement() ) < 0 ) {
        _error("giveInitViscosity - length of the first time step must be bigger than t0");
    }

    return ( t0 - atTime->giveTimeIncrement() ) / q4;
}


void
MPSMaterial :: giveEigenStrainVector(FloatArray &answer, MatResponseForm form,
                                     GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
//
// computes the strain due to creep at constant stress during the increment
// (in fact, the INCREMENT of creep strain is computed for mode == VM_Incremental)
//
{
    FloatArray KelvinEigenStrain, reducedAnswer, sigma;
    FloatMatrix C;

    double eta, dt;
    double dEtaR, etaR, L;

    MPSMaterialStatus *status = ( MPSMaterialStatus * ) this->giveStatus(gp);

    if ( mode == VM_Incremental ) {
        sigma = status->giveStressVector();       //stress vector at the beginning of time-step
        this->giveUnitComplianceMatrix(C, ReducedForm, gp, atTime);
        reducedAnswer.resize( C.giveNumberOfRows() );
        reducedAnswer.beProductOf(C, sigma);

        // flow strain increment at constant stress
        dt = atTime->giveTimeIncrement();
        eta = this->computeFlowTermViscosity(gp, atTime);

        if ( this->CoupledAnalysis == Basic ) {
            reducedAnswer.times(dt / eta);
        } else if ( this->CoupledAnalysis == MPS ) {
            // TRAPEZOIDAL INTEGRATION RULE
            if ( atTime->isTheFirstStep() ) {
                etaR = this->giveInitViscosity(atTime) /  this->computePsiR(gp, atTime, 0.);
            } else {
                etaR = status->giveFlowTermViscosity() /  this->computePsiR(gp, atTime, 0.);
            }

            dEtaR =  eta /  this->computePsiR(gp, atTime, 1.) - etaR;

            if (  fabs(dEtaR) > 1.e-4 * etaR ) {
                L = log(1 + dEtaR / etaR);
                reducedAnswer.times(dt * L / dEtaR);
            } else {
                reducedAnswer.times(dt * ( 1 - 0.5 * dEtaR / etaR +  dEtaR * dEtaR / ( 3 * etaR * etaR ) )  / etaR);
            }

            // TRAPEZOIDAL INTEGRATION RULE

            // MIDPOINT INTEGRATION RULE
            // if ( atTime->isTheFirstStep() ) {
            //   reducedAnswer.times( dt * 2. * this->computePsiR(gp, atTime, 2.) / ( eta + this->giveInitViscosity(atTime) ) );
            // } else {
            //   reducedAnswer.times( dt * 2. * this->computePsiR(gp, atTime, 2.) / ( eta + status->giveFlowTermViscosity() ) );
            // }
            // MIDPOINT INTEGRATION RULE
        } else {
            _error("giveEigenStrainVector - mode is not supported")
        }

        //computes creep component of the Kelvin Chain
        KelvinChainSolidMaterial :: giveEigenStrainVector(KelvinEigenStrain, ReducedForm, gp, atTime, mode);
        reducedAnswer.add(KelvinEigenStrain);

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
MPSMaterial :: computePointShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
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

    humDiff = this->giveHumidity(gp, atTime, 3);
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
MPSMaterial :: inverse_sorption_isotherm(double w)
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

    /*if ( ( phi < 0.2 ) || ( phi > 0.98 ) ) {
     * _error3("inverse_sorption_isotherm : Relative humidity h = %e (w=%e) is out of range", phi, w);
     * }*/
    //if ( phi < 0.20 ){ phi = 0.2;}
    //if ( phi > 0.98 ){ phi = 0.98;}

    return phi;
}

double
MPSMaterial :: giveHumidity(GaussPoint *gp, TimeStep *atTime, int option)
{
    double H_tot, H_inc;

    MPSMaterialStatus *status = ( MPSMaterialStatus * ) this->giveStatus(gp);

    // compute humidity and its increment if the step is first or humidity has not been yet computed
    if ( ( status->giveHum() == -1. ) && ( status->giveHumIncrement() == -1. ) ) {
        FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();
        Field *tf;
        int err, wflag;
        FloatArray gcoords;
        FloatArray et2, ei2; // total and incremental values of water mass

        if ( ( tf = fm->giveField(FT_HumidityConcentration) ) ) {
            gp->giveElement()->computeGlobalCoordinates( gcoords, * gp->giveCoordinates() );
            if ( ( err = tf->evaluateAt(et2, gcoords, VM_Total, atTime) ) ) {
                _error2("giveHumidity: tf->evaluateAt failed, error value %d", err);
            }

            if ( ( err = tf->evaluateAt(ei2, gcoords, VM_Incremental, atTime) ) ) {
                _error2("giveHumidity: tf->evaluateAt failed, error value %d", err);
            }

            // convert water mass to relative humidity
            H_tot = this->inverse_sorption_isotherm( et2.at(1) );
            H_inc = this->inverse_sorption_isotherm( et2.at(1) ) - this->inverse_sorption_isotherm( et2.at(1) - ei2.at(1) );
            wflag = 1;
        }

        if ( wflag == 0 ) {
            _error("giveHumidity: external fields not found");
        }

        // write field values to status
        status->setHum(H_tot);
        status->setHumIncrement(H_inc);
    } else {
        H_tot =  status->giveHum();
        H_inc =  status->giveHumIncrement();
    }

    switch ( option ) {
    case 0: return H_tot - H_inc; // returns relative humidity on the BEGINNING of the time step

    case 1: return H_tot; // returns relative humidity in the END of the time step

    case 2: return H_tot - 0.5 * H_inc; // returns relative humidity in the middle of the time step = AVERAGE

    case 3: return H_inc; // returns relative humidity INCREMENT

    default: _error2("giveHumidity: option  %d not supported", option);
    }
    return 0.; // happy compiler
}

double
MPSMaterial :: giveTemperature(GaussPoint *gp, TimeStep *atTime, int option)
{
    double T_tot, T_inc;
    MPSMaterialStatus *status = ( MPSMaterialStatus * ) this->giveStatus(gp);

    // compute humidity and its increment if the step is first or humidity has not been yet computed
    if ( ( status->giveT() == -1. ) && ( status->giveTIncrement() == -1. ) ) {
        FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();
        Field *tf;
        int err, tflag;
        FloatArray gcoords;
        FloatArray et1, ei1; // total and incremental values of temperature

        if ( ( tf = fm->giveField(FT_Temperature) ) ) {
            gp->giveElement()->computeGlobalCoordinates( gcoords, * gp->giveCoordinates() );
            if ( ( err = tf->evaluateAt(et1, gcoords, VM_Total, atTime) ) ) {
                _error2("giveTemperature: tf->evaluateAt failed, error value %d", err);
            }

            if ( ( err = tf->evaluateAt(ei1, gcoords, VM_Incremental, atTime) ) ) {
                _error2("giveTemperature: tf->evaluateAt failed, error value %d", err);
            }

            T_tot = et1.at(1);
            T_inc = ei1.at(1);

            tflag = 1;
        }

        if ( tflag == 0 ) {
            _error("giveTemperature: external fields not found");
        }

        // write field values to status
        status->setT(T_tot);
        status->setTIncrement(T_inc);
    } else {
        T_tot =  status->giveT();
        T_inc =  status->giveTIncrement();
    }

    switch ( option ) {
    case 0: return T_tot - T_inc; // returns temperature on the BEGINNING of the time step

    case 1: return T_tot; // returns temperature in the END of the time step

    case 2: return T_tot - 0.5 * T_inc; // returns temperature in the middle of the time step = AVERAGE

    case 3: return T_inc; // returns temperature INCREMENT

    default: _error2("giveTemperature: option  %d not supported", option);
    }
    return 0.; // happy compiler
}

double
MPSMaterial :: computePsiR(GaussPoint *gp, TimeStep *atTime, double option)
{
    double T, H;
    T = this->giveTemperature(gp, atTime, option);
    H = this->giveHumidity(gp, atTime, option);

    double betaRH = alphaR + ( 1. - alphaR ) * H * H;
    double betaRT = exp( QRtoR * ( 1. / this->roomTemperature - 1. / T ) );
    return betaRH * betaRT;
}

double
MPSMaterial :: computePsiS(GaussPoint *gp, TimeStep *atTime)
{
    double AverageTemp, AverageHum;

    AverageTemp = this->giveTemperature(gp, atTime, 2);
    AverageHum = this->giveHumidity(gp, atTime, 2);

    double betaSH = alphaS + ( 1. - alphaS ) * AverageHum * AverageHum;
    double betaRT = exp( QStoR * ( 1. / this->roomTemperature - 1. /  AverageTemp ) );

    return betaSH * betaRT;
}

double
MPSMaterial :: computePsiE(GaussPoint *gp, TimeStep *atTime)
{
    double AverageTemp, AverageHum;

    AverageTemp = this->giveTemperature(gp, atTime, 2);
    AverageHum = this->giveHumidity(gp, atTime, 2);

    double betaEH = 1. / ( 1. +  pow( ( alphaE * ( 1. - AverageHum ) ), 4. ) );
    double betaET = exp( QEtoR * ( 1. /  this->roomTemperature - 1. / AverageTemp ) );

    return betaEH * betaET;
}

double
MPSMaterial :: computeEquivalentTime(GaussPoint *gp, TimeStep *atTime, int option)
{
    double tEquiv = 0.;
    double PsiE;

    PsiE = computePsiE(gp, atTime);

    if ( atTime->isTheFirstStep() ) {
        if ( option == 0 ) { // gives time in the middle of the timestep
            return relMatAge - atTime->giveTimeIncrement() + PsiE * ( 0.5 * atTime->giveTimeIncrement() );
        } else if ( option == 1 ) { // gives time in the middle of the timestep - for UPDATING
            return relMatAge - atTime->giveTimeIncrement() + PsiE *atTime->giveTimeIncrement();
        } else {
            _error("computeEquivalentTime - mode is not supported")
        }
    } else {
        MPSMaterialStatus *status = ( MPSMaterialStatus * ) this->giveStatus(gp);
        tEquiv = status->giveEquivalentTime();

        if ( option == 0 ) { // gives time in the middle of the timestep
            tEquiv = tEquiv + PsiE *  0.5 * atTime->giveTimeIncrement();
        } else if ( option == 1 ) { // gives time in the middle of the timestep - for UPDATING
            tEquiv = tEquiv + PsiE *atTime->giveTimeIncrement();
        } else {
            _error("computeEquivalentTime - mode is not supported")
        }

        return tEquiv;
    }
    return tEquiv; // happy compiler
}
} // end namespace oofem
