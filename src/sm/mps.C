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

#include "mps.h"
#include "mathfem.h"
#include "gausspoint.h"
#include "timestep.h"
#include "engngm.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(MPSMaterial);

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
        OOFEM_ERROR("can't write into NULL stream");
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
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    double fc, c, wc, ac;

    /* scaling factor transforming PREDICTED stiffnesses q1, q2, q3 and q4 to desired
     *  default value is 1.0 = no change
     *  e.g. if the stiffness should be in MPa, then stiffnessFactor = 1.e6 */
    double stiffnessFactor;

    KelvinChainSolidMaterial :: initializeFrom(ir);

    // checking timefactor - for MPS material must be equal to 1.
    double tf;
    IR_GIVE_FIELD(ir, tf, _IFT_RheoChainMaterial_timefactor);
    if ( tf != 1. ) {
        OOFEM_ERROR("for MPS material timefactor must be equal to 1.");
    }

    int mode = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, mode, _IFT_MPSMaterial_mode);

    if ( mode == 0 ) { // default - estimate model parameters q1,..,q4 from composition
        IR_GIVE_FIELD(ir, fc, _IFT_MPSMaterial_fc); // 28-day standard cylinder compression strength [MPa]
        IR_GIVE_FIELD(ir,  c, _IFT_MPSMaterial_cc); // cement content of concrete [kg/m^3]
        IR_GIVE_FIELD(ir, wc, _IFT_MPSMaterial_wc); // ratio (by weight) of water to cementitious material
        IR_GIVE_FIELD(ir, ac, _IFT_MPSMaterial_ac); // ratio (by weight) of aggregate to cement

        stiffnessFactor = 1.;
        IR_GIVE_OPTIONAL_FIELD(ir, stiffnessFactor, _IFT_MPSMaterial_stiffnessfactor); // ratio (by weight) of aggregate to cement

        this->predictParametersFrom(fc, c, wc, ac, stiffnessFactor);
    } else { // read model parameters for creep
        IR_GIVE_FIELD(ir, q1, _IFT_MPSMaterial_q1);
        IR_GIVE_FIELD(ir, q2, _IFT_MPSMaterial_q2);
        IR_GIVE_FIELD(ir, q3, _IFT_MPSMaterial_q3);
        IR_GIVE_FIELD(ir, q4, _IFT_MPSMaterial_q4);
    }

    IR_GIVE_FIELD(ir, talpha, _IFT_MPSMaterial_talpha);
    IR_GIVE_FIELD(ir, lambda0, _IFT_MPSMaterial_lambda0);

    int type = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, type, _IFT_MPSMaterial_coupledanalysistype);
    if ( type >= 2 ) {
        OOFEM_ERROR("CoupledAnalysisType must be equal to 0 or 1");
    }

    this->CoupledAnalysis = ( coupledAnalysisType ) type;
    this->kSh = 0.;

    if ( this->CoupledAnalysis == MPS ) {
        // muS- the only parameter necessary for evaluation of the flow-term viscosity
        // muS = c0 * c1 * q4 [1/(Pa*s)]
        IR_GIVE_FIELD(ir, muS, _IFT_MPSMaterial_mus);
        IR_GIVE_FIELD(ir, kappaT, _IFT_MPSMaterial_kappat); //[-] replaces ln(h) in diff equation for MPS
        this->ct = 0.;
        IR_GIVE_OPTIONAL_FIELD(ir, ct, _IFT_MPSMaterial_ct);

        // age when drying or thermal changes begin [days] - necessary to determine initial viscosity
        IR_GIVE_FIELD(ir, t0, _IFT_MPSMaterial_t0);
        IR_GIVE_OPTIONAL_FIELD(ir, kSh, _IFT_MPSMaterial_ksh);

        // Parameters for desorption isotherm
        IR_GIVE_FIELD(ir, w_h, _IFT_MPSMaterial_wh);
        IR_GIVE_FIELD(ir, n, _IFT_MPSMaterial_ncoeff);
        IR_GIVE_FIELD(ir, a, _IFT_MPSMaterial_a);

        this->roomTemperature = 298.; // reference room temperature for MPS algorithm [K]

        // ACTIVATION ENERGIES
        this->QEtoR = 2700.; // [K] value according to Bazant 1995
        IR_GIVE_OPTIONAL_FIELD(ir, QEtoR, _IFT_MPSMaterial_qetor);
        this->QRtoR = 5000.; // [K] value according to Bazant 1995
        IR_GIVE_OPTIONAL_FIELD(ir, QRtoR, _IFT_MPSMaterial_qrtor);
        this->QStoR = 3000.; // [K] value according to Bazant 1995
        IR_GIVE_OPTIONAL_FIELD(ir, QStoR, _IFT_MPSMaterial_qstor);

        this->alphaE = 10.; // according to Bazant 1995
        IR_GIVE_OPTIONAL_FIELD(ir, alphaE, _IFT_MPSMaterial_alphae);
        this->alphaR = 0.1; // according to Bazant 1995
        IR_GIVE_OPTIONAL_FIELD(ir, alphaR, _IFT_MPSMaterial_alphar);
        this->alphaS = 0.1; // according to Bazant 1995
        IR_GIVE_OPTIONAL_FIELD(ir, alphaS, _IFT_MPSMaterial_alphas);
    }

    return IRRT_OK;
}


void
MPSMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    KelvinChainSolidMaterial :: giveRealStressVector(answer, gp, reducedStrain, tStep);

    MPSMaterialStatus *status = static_cast< MPSMaterialStatus * >( this->giveStatus(gp) );

    if ( this->CoupledAnalysis == MPS ) {
        status->setEquivalentTime( this->computeEquivalentTime(gp, tStep, 1) );
        // set humidity and temperature to zero
        status->setHum(-1.);
        status->setHumIncrement(-1.);
        status->setT(-1.);
        status->setTIncrement(-1.);
    }
}

#if 0
void
MPSMaterial :: updateYourself(GaussPoint *gp, TimeStep *tStep)
{
    MPSMaterialStatus *status = static_cast< MPSMaterialStatus * >( this->giveStatus(gp) );

    if ( this->CoupledAnalysis == MPS ) {
        status->setEquivalentTime( this->computeEquivalentTime(gp, tStep, 1) );
        // set humidity and temperature to zero
        status->setHum(-1.);
        status->setHumIncrement(-1.);
        status->setT(-1.);
        status->setTIncrement(-1.);
    }

    // now we call Solidifying Kelvin Chain to update itself
    // at the end of updating it will call material status to be updated.
    KelvinChainSolidMaterial :: updateYourself(gp, tStep);
}
#endif


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
                                         GaussPoint *gp,
                                         TimeStep *tStep,
                                         ValueModeType mode)
{
    if ( ( mode != VM_Total ) && ( mode != VM_Incremental ) ) {
        OOFEM_ERROR("unsupported mode");
    }


    if ( this->kSh == 0. ) {
        answer.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
        answer.zero();
    } else {
        this->computePointShrinkageStrainVector(answer, gp, tStep);
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
    q2  = 1.e-12 * stiffnessFactor * 185.4 * pow(c, 0.5) * pow(fc, -0.9);
    q3 = 0.29 * pow(wc, 4.) * q2;
    q4  = 1.e-12 * stiffnessFactor * 20.3 * pow(ac, -0.7);

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
    while ( 0.5 * this->endOfTimeOfInterest >= Tau1 * pow( 10.0, ( double ) ( j - 1 ) ) ) {
        j++;
    }

    this->nUnits = j;

    this->charTimes.resize(this->nUnits);

    for ( mu = 1; mu <= this->nUnits; mu++ ) {
        charTimes.at(mu) = Tau1 * pow(10., mu - 1);
    }
}


void
MPSMaterial :: computeCharCoefficients(FloatArray &answer, double tStep)
{
    int mu;
    double tau0, tauMu;

    // constant "n" is assumed to be equal to 0.1 (typical value)

    // modulus of elasticity of the first unit of Kelvin chain.
    // (aging elastic spring with retardation time = 0)
    double lambda0ToPowN = pow(lambda0, 0.1);
    tau0 = pow(2 * this->giveCharTime(1) / sqrt(10.0), 0.1);
    EspringVal = 1. / ( q2 * log(1.0 + tau0 / lambda0ToPowN) - q2 * tau0 / ( 10.0 * lambda0ToPowN + 10.0 * tau0 ) );

    // evaluation of moduli of elasticity for the remaining units
    // (Solidifying kelvin units with retardation times tauMu)
    answer.resize(nUnits);
    for ( mu = 1; mu <= this->nUnits; mu++ ) {
        tauMu = pow(2 * this->giveCharTime(mu), 0.1);
        answer.at(mu) = 10. * pow(1 + tauMu / lambda0ToPowN, 2) / ( log(10.0) * q2 * ( tauMu / lambda0ToPowN ) * ( 0.9 + tauMu / lambda0ToPowN ) );
        this->charTimes.at(mu) *= 1.35;
    }

    answer.at(nUnits) /= 1.2;   // modulus of the last unit is reduced
}


double
MPSMaterial :: giveEModulus(GaussPoint *gp, TimeStep *tStep)
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

    double Cf = 0.0; // incremental viscous flow compliance
    double eta, dt;
    double dEtaR, etaR, L;

    // contribution of the solidifying Kelving chain
    sum = KelvinChainSolidMaterial :: giveEModulus(gp, tStep);

    v = computeSolidifiedVolume(gp, tStep);

    dt = tStep->giveTimeIncrement();
    eta = this->computeFlowTermViscosity(gp, tStep);

    //incremental viscous flow compliance

    if ( this->CoupledAnalysis == Basic ) {
        Cf = 0.5 * ( dt ) / eta;
    } else if ( this->CoupledAnalysis == MPS ) {
        MPSMaterialStatus *status = static_cast< MPSMaterialStatus * >( this->giveStatus(gp) );

        // TRAPEZOIDAL INTEGRATION RULE
        if ( tStep->isTheFirstStep() ) {
            etaR = this->giveInitViscosity(tStep) /  this->computePsiR(gp, tStep, 0);
        } else {
            etaR = status->giveFlowTermViscosity() /  this->computePsiR(gp, tStep, 0);
        }

        dEtaR =  eta /  this->computePsiR(gp, tStep, 1) - etaR;
        if (  fabs(dEtaR) > 1.e-4 * etaR ) {
            L = log(1 + dEtaR / etaR);
            Cf = dt * ( 1. - etaR * L / dEtaR ) / dEtaR;
        } else {
            Cf = dt * ( 0.5 - dEtaR / ( 3 * etaR ) ) / etaR;
        }

        // TRAPEZOIDAL INTEGRATION RULE

        // MIDPOINT INTEGRATION RULE
        // if ( tStep->isTheFirstStep() ) {
        //  Cf = dt * this->computePsiR(gp, tStep, 2.) / (eta + this->giveInitViscosity(tStep) );
        // } else {
        //  Cf = dt * this->computePsiR(gp, tStep, 2.) / (eta + status->giveFlowTermViscosity() );
        // }
        // TRAPEZOIDAL INTEGRATION RULE
    } else {
        OOFEM_ERROR("mode is not supported");
    }

    return 1. / ( q1 + 1. / ( EspringVal * v ) + sum  +  Cf );
}


double
MPSMaterial :: computeSolidifiedVolume(GaussPoint *gp, TimeStep *tStep)
// compute the relative volume of the solidified material at given age (in days)
{
    double m, alpha;
    double atAge;     // (equivalent age)

    // standard values of exponents - empirical constants
    m = 0.5;

    alpha = q3 / q2;

    if ( this->CoupledAnalysis == Basic ) {
        atAge = relMatAge + ( tStep->giveTargetTime() - 0.5 * tStep->giveTimeIncrement() );
    } else {
        atAge = computeEquivalentTime(gp, tStep, 0);
    }

    return 1. / ( alpha + pow(lambda0 / atAge, m) );
}


double
MPSMaterial :: computeBetaMu(GaussPoint *gp, TimeStep *tStep, int Mu)
{
    double betaMu;
    double deltaT;
    double tauMu;


    deltaT = tStep->giveTimeIncrement();

    if ( this->CoupledAnalysis == MPS ) {
        deltaT *=  0.5 * ( this->computePsiR(gp, tStep, 0) + this->computePsiR(gp, tStep, 1) );
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
MPSMaterial :: computeLambdaMu(GaussPoint *gp, TimeStep *tStep, int Mu)
{
    double lambdaMu;
    double deltaT;
    double tauMu;

    deltaT = tStep->giveTimeIncrement();

    if ( this->CoupledAnalysis == MPS ) {
        deltaT *=  0.5 * ( this->computePsiR(gp, tStep, 0) + this->computePsiR(gp, tStep, 1) );
    }

    tauMu = this->giveCharTime(Mu);

    if ( deltaT / tauMu < 1.e-5 ) {
        lambdaMu = 1 - 0.5 * ( deltaT / tauMu ) + 1 / 6 * ( pow(deltaT / tauMu, 2) ) - 1 / 24 * ( pow(deltaT / tauMu, 3) );
    } else if ( deltaT / tauMu > 30 ) {
        lambdaMu = tauMu / deltaT;
    } else {
        lambdaMu = ( 1.0 -  exp(-( deltaT ) / tauMu) ) * tauMu / deltaT;
    }

    return lambdaMu;
}

double
MPSMaterial :: computeFlowTermViscosity(GaussPoint *gp, TimeStep *tStep)
{
    double eta = 0.0, tHalfStep;

    double prevEta, PsiS, A, B, e, dt;
    double T_new, T_old, H_new, H_old;
    double reductFactor;

    if ( this->CoupledAnalysis == Basic ) {
        tHalfStep = relMatAge + ( tStep->giveTargetTime() - 0.5 * tStep->giveTimeIncrement() );
        eta = tHalfStep / q4;
    } else if ( this->CoupledAnalysis == MPS ) {
        MPSMaterialStatus *status = static_cast< MPSMaterialStatus * >( this->giveStatus(gp) );

        // check whether this viscosity has been already computed
        if ( status->giveFlowTermViscosityTemp() != -1.  &&  !( tStep->isTheFirstStep() ) ) {
            return status->giveFlowTermViscosityTemp();

            // no, the viscosity needs to be evaluated now
        } else {
            if ( tStep->isTheFirstStep() ) { // if the time step is the first one, ask for an initial value of viscosity
                prevEta = this->giveInitViscosity(tStep);
            } else {
                // asks for the value of viscosity from the end of the last time-step
                prevEta = status->giveFlowTermViscosity();
            }

            dt = tStep->giveTimeIncrement();
            // evaluate auxiliary factors A and B
            T_new = this->giveTemperature(gp, tStep, 1);
            T_old = this->giveTemperature(gp, tStep, 0);
            H_new = this->giveHumidity(gp, tStep, 1);
            H_old = this->giveHumidity(gp, tStep, 0);

            PsiS = this->computePsiS(gp, tStep); // evaluated in the middle of the time step

            // original version
            //A = sqrt( muS * fabs( T_new * log(H_new) - T_old * log(H_old) ) / ( dt * this->roomTemperature ) );

            if ( this->ct == 0. ) {
                reductFactor = 1.;
            } else if ( ( status->giveTmax() - T_new < 0. ) || tStep->isTheFirstStep() ) {
                status->setTmax(T_new);
                reductFactor = 1.;
            } else {
                reductFactor = exp( -this->ct * fabs( T_new - status->giveTmax() ) );
            }

            A = sqrt( muS * ( kappaT * reductFactor * fabs(T_new - T_old) + 0.5 * ( T_new + T_old ) * fabs( log(H_new) - log(H_old) ) ) / ( dt * this->roomTemperature ) );
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
        OOFEM_ERROR("mode is not supported");
    }


    return eta;
}

// returns initial value of the flow term viscosity
double
MPSMaterial :: giveInitViscosity(TimeStep *tStep)
{
    if ( ( t0 - tStep->giveTimeIncrement() ) < 0 ) {
        OOFEM_ERROR("length of the first time step must be bigger than t0");
    }

    return ( t0 - tStep->giveTimeIncrement() ) / q4;
}


void
MPSMaterial :: giveEigenStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode)
//
// computes the strain due to creep at constant stress during the increment
// (in fact, the INCREMENT of creep strain is computed for mode == VM_Incremental)
//
{
    FloatArray KelvinEigenStrain, reducedAnswer, sigma;
    FloatMatrix C;

    double eta, dt;
    double dEtaR, etaR, L;

    MPSMaterialStatus *status = static_cast< MPSMaterialStatus * >( this->giveStatus(gp) );

    if ( mode == VM_Incremental ) {
        sigma = status->giveStressVector();       //stress vector at the beginning of time-step
        this->giveUnitComplianceMatrix(C, gp, tStep);
        reducedAnswer.resize( C.giveNumberOfRows() );
        reducedAnswer.beProductOf(C, sigma);

        // flow strain increment at constant stress
        dt = tStep->giveTimeIncrement();
        eta = this->computeFlowTermViscosity(gp, tStep);

        if ( this->CoupledAnalysis == Basic ) {
            reducedAnswer.times(dt / eta);
        } else if ( this->CoupledAnalysis == MPS ) {
            // TRAPEZOIDAL INTEGRATION RULE
            if ( tStep->isTheFirstStep() ) {
                etaR = this->giveInitViscosity(tStep) /  this->computePsiR(gp, tStep, 0);
            } else {
                etaR = status->giveFlowTermViscosity() /  this->computePsiR(gp, tStep, 0);
            }

            dEtaR =  eta /  this->computePsiR(gp, tStep, 1) - etaR;

            if (  fabs(dEtaR) > 1.e-4 * etaR ) {
                L = log(1 + dEtaR / etaR);
                reducedAnswer.times(dt * L / dEtaR);
            } else {
                reducedAnswer.times(dt * ( 1 - 0.5 * dEtaR / etaR +  dEtaR * dEtaR / ( 3 * etaR * etaR ) )  / etaR);
            }

            // TRAPEZOIDAL INTEGRATION RULE

            // MIDPOINT INTEGRATION RULE
            // if ( tStep->isTheFirstStep() ) {
            //   reducedAnswer.times( dt * 2. * this->computePsiR(gp, tStep, 2) / ( eta + this->giveInitViscosity(tStep) ) );
            // } else {
            //   reducedAnswer.times( dt * 2. * this->computePsiR(gp, tStep, 2) / ( eta + status->giveFlowTermViscosity() ) );
            // }
            // MIDPOINT INTEGRATION RULE
        } else {
            OOFEM_ERROR("mode is not supported")
        }

        //computes creep component of the Kelvin Chain
        KelvinChainSolidMaterial :: giveEigenStrainVector(KelvinEigenStrain, gp, tStep, mode);
        reducedAnswer.add(KelvinEigenStrain);

        answer = reducedAnswer;
        return;
    } else {
        /* error - total mode not implemented yet */
        OOFEM_ERROR("mode is not supported");
    }
}


void
MPSMaterial :: computePointShrinkageStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
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

    humDiff = this->giveHumidity(gp, tStep, 3);
    EpsSh = humDiff * kSh;

    fullAnswer.resize(size);
    fullAnswer.zero();
    fullAnswer.at(1) = fullAnswer.at(2) = fullAnswer.at(3) = EpsSh;

    StructuralMaterial :: giveReducedSymVectorForm( answer, fullAnswer, gp->giveMaterialMode() );
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
    // relative humidity
    double phi = exp( a * ( 1.0 - pow( ( w_h / w ), ( n ) ) ) );

    /*if ( ( phi < 0.2 ) || ( phi > 0.98 ) ) {
     * OOFEM_ERROR("Relative humidity h = %e (w=%e) is out of range", phi, w);
     * }*/
    //if ( phi < 0.20 ){ phi = 0.2;}
    //if ( phi > 0.98 ){ phi = 0.98;}

    return phi;
}

double
MPSMaterial :: giveHumidity(GaussPoint *gp, TimeStep *tStep, int option)
{
    double H_tot = 0.0, H_inc = 0.0;

    MPSMaterialStatus *status = static_cast< MPSMaterialStatus * >( this->giveStatus(gp) );

    // compute humidity and its increment if the step is first or humidity has not been yet computed
    if ( ( status->giveHum() == -1. ) && ( status->giveHumIncrement() == -1. ) ) {
        FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();

        FM_FieldPtr tf;
        int err, wflag = 0;
        FloatArray gcoords;
        FloatArray et2, ei2; // total and incremental values of water mass

        if ( ( tf = fm->giveField(FT_HumidityConcentration) ) ) {
            gp->giveElement()->computeGlobalCoordinates( gcoords, * gp->giveNaturalCoordinates() );
            if ( ( err = tf->evaluateAt(et2, gcoords, VM_Total, tStep) ) ) {
                OOFEM_ERROR("tf->evaluateAt failed, error value %d", err);
            }

            if ( ( err = tf->evaluateAt(ei2, gcoords, VM_Incremental, tStep) ) ) {
                OOFEM_ERROR("tf->evaluateAt failed, error value %d", err);
            }

            // convert water mass to relative humidity
            H_tot = this->inverse_sorption_isotherm( et2.at(1) );
            H_inc = this->inverse_sorption_isotherm( et2.at(1) ) - this->inverse_sorption_isotherm( et2.at(1) - ei2.at(1) );
            wflag = 1;
        }

        if ( wflag == 0 ) {
            OOFEM_ERROR("external fields not found");
        }

        // write field values to status
        status->setHum(H_tot);
        status->setHumIncrement(H_inc);
    } else {
        H_tot = status->giveHum();
        H_inc = status->giveHumIncrement();
    }

    switch ( option ) {
    case 0: return H_tot - H_inc; // returns relative humidity on the BEGINNING of the time step

    case 1: return H_tot; // returns relative humidity in the END of the time step

    case 2: return H_tot - 0.5 * H_inc; // returns relative humidity in the middle of the time step = AVERAGE

    case 3: return H_inc; // returns relative humidity INCREMENT

    default: OOFEM_ERROR("option  %d not supported", option);
    }
    return 0.; // happy compiler
}

double
MPSMaterial :: giveTemperature(GaussPoint *gp, TimeStep *tStep, int option)
{
    double T_tot = 0.0, T_inc = 0.0;
    MPSMaterialStatus *status = static_cast< MPSMaterialStatus * >( this->giveStatus(gp) );

    // compute humidity and its increment if the step is first or humidity has not been yet computed
    if ( ( status->giveT() == -1. ) && ( status->giveTIncrement() == -1. ) ) {
        FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();

        FM_FieldPtr tf;
        int err, tflag = 0;
        FloatArray gcoords;
        FloatArray et1, ei1; // total and incremental values of temperature

        if ( ( tf = fm->giveField(FT_Temperature) ) ) {
            gp->giveElement()->computeGlobalCoordinates( gcoords, * gp->giveNaturalCoordinates() );
            if ( ( err = tf->evaluateAt(et1, gcoords, VM_Total, tStep) ) ) {
                OOFEM_ERROR("tf->evaluateAt failed, error value %d", err);
            }

            if ( ( err = tf->evaluateAt(ei1, gcoords, VM_Incremental, tStep) ) ) {
                OOFEM_ERROR("tf->evaluateAt failed, error value %d", err);
            }

            T_tot = et1.at(1);
            T_inc = ei1.at(1);

            tflag = 1;
        }

        if ( tflag == 0 ) {
            OOFEM_ERROR("external fields not found");
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

    default: OOFEM_ERROR("option  %d not supported", option);
    }
    return 0.; // happy compiler
}

double
MPSMaterial :: computePsiR(GaussPoint *gp, TimeStep *tStep, int option)
{
    double T, H;
    T = this->giveTemperature(gp, tStep, option);
    H = this->giveHumidity(gp, tStep, option);

    double betaRH = alphaR + ( 1. - alphaR ) * H * H;
    double betaRT = exp( QRtoR * ( 1. / this->roomTemperature - 1. / T ) );
    return betaRH * betaRT;
}

double
MPSMaterial :: computePsiS(GaussPoint *gp, TimeStep *tStep)
{
    double AverageTemp, AverageHum;

    AverageTemp = this->giveTemperature(gp, tStep, 2);
    AverageHum = this->giveHumidity(gp, tStep, 2);

    double betaSH = alphaS + ( 1. - alphaS ) * AverageHum * AverageHum;
    double betaRT = exp( QStoR * ( 1. / this->roomTemperature - 1. /  AverageTemp ) );

    return betaSH * betaRT;
}

double
MPSMaterial :: computePsiE(GaussPoint *gp, TimeStep *tStep)
{
    double AverageTemp, AverageHum;

    AverageTemp = this->giveTemperature(gp, tStep, 2);
    AverageHum = this->giveHumidity(gp, tStep, 2);

    double betaEH = 1. / ( 1. +  pow( ( alphaE * ( 1. - AverageHum ) ), 4. ) );
    double betaET = exp( QEtoR * ( 1. /  this->roomTemperature - 1. / AverageTemp ) );

    return betaEH * betaET;
}

double
MPSMaterial :: computeEquivalentTime(GaussPoint *gp, TimeStep *tStep, int option)
{
    double tEquiv = 0.;
    double PsiE;

    PsiE = computePsiE(gp, tStep);

    if ( tStep->isTheFirstStep() ) {
        if ( option == 0 ) { // gives time in the middle of the timestep
            return relMatAge - tStep->giveTimeIncrement() + PsiE * ( 0.5 * tStep->giveTimeIncrement() );
        } else if ( option == 1 ) { // gives time in the middle of the timestep - for UPDATING
            return relMatAge - tStep->giveTimeIncrement() + PsiE *tStep->giveTimeIncrement();
        } else {
            OOFEM_ERROR("mode is not supported")
        }
    } else {
        MPSMaterialStatus *status = static_cast< MPSMaterialStatus * >( this->giveStatus(gp) );
        tEquiv = status->giveEquivalentTime();

        if ( option == 0 ) { // gives time in the middle of the timestep
            tEquiv = tEquiv + PsiE *  0.5 * tStep->giveTimeIncrement();
        } else if ( option == 1 ) { // gives time in the middle of the timestep - for UPDATING
            tEquiv = tEquiv + PsiE *tStep->giveTimeIncrement();
        } else {
            OOFEM_ERROR("mode is not supported")
        }

        return tEquiv;
    }
    return tEquiv; // happy compiler
}
} // end namespace oofem
