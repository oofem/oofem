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

#include "b3solidmat.h"
#include "mathfem.h"
#include "gausspoint.h"
#include "timestep.h"
#include "contextioerr.h"
#include "datastream.h"
#include "classfactory.h"
#include "engngm.h"
#include "b3mat.h"

namespace oofem {
REGISTER_Material(B3SolidMaterial);

void
B3SolidMaterial :: initializeFrom(InputRecord &ir)
{
    // ph!!!
    //KelvinChainMaterial :: initializeFrom(ir);

    //
    // NOTE
    //
    // this material model is unit-dependent!
    // units must be in MPa, m???, MN.
    //

    double fc = -1.0, c = -1.0, wc = -1.0, ac = -1.0, alpha1 = 1.0, alpha2 = 1.0;
    double initHum; //MPS shrinkage parameter - initial value of humidity (range 0.2-0.98)
    double finalHum; //MPS shrinkage parameter - final value of humidity (range 0.2-0.98)

    // EmoduliMode has to be set first because characteristic times are computed from RheoChM initialization
    this->EmoduliMode = 0;
    // retardation spectrum or least square method is used. Retardation spectrum is default (EmoduliMode==0)
    IR_GIVE_OPTIONAL_FIELD(ir, EmoduliMode, _IFT_B3SolidMaterial_emodulimode);

    // characteristic time, usually 1 for analysis running in days
    this->lambda0 = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, lambda0, _IFT_B3SolidMaterial_lambda0);

   KelvinChainMaterial :: initializeFrom(ir);

    int mode = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, mode, _IFT_B3Material_mode);

    if ( mode == 0 ) { // default - estimate model parameters q1,..,q5 from composition
        IR_GIVE_FIELD(ir, fc, _IFT_B3Material_fc); // 28-day standard cylinder compression strength [MPa]
        IR_GIVE_FIELD(ir,  c, _IFT_B3Material_cc); // cement content of concrete [kg/m^3]
        IR_GIVE_FIELD(ir, wc, _IFT_B3Material_wc); // ratio (by weight) of water to cementitious material
        IR_GIVE_FIELD(ir, ac, _IFT_B3Material_ac); // ratio (by weight) of aggregate to cement
        IR_GIVE_FIELD(ir, t0, _IFT_B3Material_t0); // age when drying begins [days]
    } else { // read model parameters for creep
        IR_GIVE_FIELD(ir, q1, _IFT_B3Material_q1);
        IR_GIVE_FIELD(ir, q2, _IFT_B3Material_q2);
        IR_GIVE_FIELD(ir, q3, _IFT_B3Material_q3);
        IR_GIVE_FIELD(ir, q4, _IFT_B3Material_q4);
    }

    // default value = 0; it can be used for basic creep without any external fields
    this->MicroPrestress = 0; //if = 1 computation exploiting Microprestress solidification theory is done;
    IR_GIVE_OPTIONAL_FIELD(ir, MicroPrestress, _IFT_B3SolidMaterial_microprestress);

    // is MPS theory used?
    if ( this->MicroPrestress == 1 ) {
        // microprestress-sol-theory: read data for microprestress evaluation
        // constant c0 [MPa^-1 * day^-1]
        IR_GIVE_FIELD(ir, c0, _IFT_B3SolidMaterial_c0);
        // constant c1 (=C1*R*T/M)
        IR_GIVE_FIELD(ir, c1, _IFT_B3SolidMaterial_c1);
        // t0- necessary for the initial value of microprestress = S0 (age when drying begins)
	//        IR_GIVE_FIELD(ir, t0, _IFT_B3Material_t0);
	this->tS0 = this->relMatAge;

        // microprestress-sol-theory: read data for inverse desorption isotherm
        IR_GIVE_FIELD(ir, w_h, _IFT_B3Material_wh);
        IR_GIVE_FIELD(ir, n, _IFT_B3Material_ncoeff);
        IR_GIVE_FIELD(ir, a, _IFT_B3Material_a);
    }

    // read shrinkage mode
    int shm = 0;
    IR_GIVE_FIELD(ir, shm, _IFT_B3Material_shmode);
    this->shMode = ( b3ShModeType ) shm;

    if ( this->shMode == B3_PointShrinkage ) {   // #2 in enumerator

        if ( this->MicroPrestress == 0 ) {
            throw ValueInputException(ir, _IFT_B3SolidMaterial_microprestress, "to use B3_PointShrinkageMPS - MicroPrestress must be = 1"); //or else no external fiels would be found
        }

        kSh = -1;
        initHum = -1;
        finalHum = -1;

        IR_GIVE_OPTIONAL_FIELD(ir, kSh, _IFT_B3SolidMaterial_ksh);
        IR_GIVE_OPTIONAL_FIELD(ir, initHum, _IFT_B3SolidMaterial_initialhumidity);
        IR_GIVE_OPTIONAL_FIELD(ir, finalHum, _IFT_B3SolidMaterial_finalhumidity);

        // either kSh or initHum and finalHum must be given in input record
        if ( !( ( this->kSh != -1 ) || ( ( initHum != -1 ) && ( finalHum != -1 ) ) ) ) {
            throw ValueInputException(ir, "none", "either kSh or initHum and finalHum must be given in input record");
        }
        
        if ( ( initHum < 0.2 || initHum > 0.98 ) && this->kSh == -1 ) {
            throw ValueInputException(ir, _IFT_B3SolidMaterial_initialhumidity, "initital humidity out of range (0.2 - 0.98)");
        }
        if ( ( finalHum < 0.2 || finalHum > 0.98 ) && this->kSh == -1 ) {
            throw ValueInputException(ir, _IFT_B3SolidMaterial_finalhumidity, "final humidity out of range (0.2 - 0.98)");
        }
        

        if ( this->kSh == -1 ) { // predict kSh from composition
            IR_GIVE_OPTIONAL_FIELD(ir, alpha1, _IFT_B3Material_alpha1);       // influence of cement type
            IR_GIVE_OPTIONAL_FIELD(ir, alpha2, _IFT_B3Material_alpha2);       // influence of curing type
            this->kSh = alpha1 * alpha2 * ( 1 - pow(finalHum, 3.) ) * ( 0.019 * pow(wc * c, 2.1) * pow(fc, -0.28) + 270 ) * 1.e-6 / fabs(initHum - finalHum);
        }

    } else if ( this->shMode == B3_AverageShrinkage ) {
        IR_GIVE_FIELD(ir, ks, _IFT_B3Material_ks);   // cross-section shape factor
        /*
         * use ks = 1.0 for an infinite slab
         *        = 1.15 for an infinite cylinder
         *        = 1.25 for an infinite square prism
         *        = 1.30 for a sphere
         *        = 1.55 for a cube
         */
        IR_GIVE_FIELD(ir, hum, _IFT_B3Material_hum);         // relative humidity of the environment
        IR_GIVE_FIELD(ir, vs, _IFT_B3Material_vs);         // volume-to-surface ratio (in m???)
        if ( mode == 0 ) {       // default mode - estimate model parameters kt, EpsSinf, q5 from composition
            IR_GIVE_OPTIONAL_FIELD(ir, alpha1, _IFT_B3Material_alpha1);             // influence of cement type
            IR_GIVE_OPTIONAL_FIELD(ir, alpha2, _IFT_B3Material_alpha2);             // influence of curing type
        } else {         // read model parameters
            IR_GIVE_FIELD(ir, t0, _IFT_B3Material_t0); // age when drying begins [days]
            IR_GIVE_FIELD(ir, kt, _IFT_B3Material_kt);
            IR_GIVE_FIELD(ir, EpsSinf, _IFT_B3Material_EpsSinf);
            IR_GIVE_FIELD(ir, q5, _IFT_B3Material_q5);
        }
    }

    // evaluate the total water content [kg/m^3]
    w = wc * c;
    // estimate the conventional modulus ??? what if fc is not given ???
    E28 = 4733. * sqrt(fc);
    // estimate parameters from composition
    if ( mode == 0 ) {
        this->predictParametersFrom(fc, c, wc, ac, t0, alpha1, alpha2);
    }
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
    q1  = 126.74271 / sqrt(fc) * 1e-6;
    q2  = 185.4 * pow(c, 0.5) * pow(fc, -0.9) * 1.e-6; 
    q3  = 0.29 * pow(wc, 4.) * q2;
    q4  = 20.3 * pow(ac, -0.7) * 1.e-6;

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


double
B3SolidMaterial :: giveEModulus(GaussPoint *gp, TimeStep *tStep) const
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

    double v = computeSolidifiedVolume(tStep);
    double eta = this->computeFlowTermViscosity(gp, tStep);     //evaluated in the middle of the time-step

    ///@todo THREAD UNSAFE!
    double t_halfstep = this->relMatAge - this->castingTime + ( tStep->giveTargetTime() - 0.5 * tStep->giveTimeIncrement() ) / timeFactor;
    this->updateEparModuli( t_halfstep, gp, tStep );

    if ( this->EmoduliMode == 0 ) { //retardation spectrum used
        double sum;
        //  compute compliance of the Kelvin chain
        sum = 1. / KelvinChainMaterial :: giveEModulus(gp, tStep);
        // add compliance of non-aging spring
        sum += 1. / EspringVal;
        // convert to stiffness
        chainStiffness = 1. / sum;

    } else {   //least-squares method used
        chainStiffness = KelvinChainMaterial :: giveEModulus(gp, tStep);
    }

    return 1. / ( q1 +  1. / (chainStiffness * v)  + 0.5 * ( tStep->giveTimeIncrement() / timeFactor ) / eta );
}


void
B3SolidMaterial :: updateEparModuli(double tPrime, GaussPoint *gp, TimeStep *tStep) const
{
    /*
     * Since the elastic moduli are constant in time it is necessary to update them only once
     * on the beginning of the computation
     */
    if ( this->EparVal.isEmpty() ) {
        RheoChainMaterial :: updateEparModuli(tPrime, gp, tStep);
    }
}

void
B3SolidMaterial :: computeCharTimes()
{
    /*
     * This function generates discrete characteristic times
     * according to rules which guarantee a good approximation
     * of the relaxation or creep function by the Dirichlet series
     */

    double Tau1;
    int j;

    if ( this->begOfTimeOfInterest == -1. ) {
        this->begOfTimeOfInterest = 0.001 * lambda0;
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
    this->charTimes.zero();

    for ( int mu = 1; mu <= this->nUnits; mu++ ) {
        charTimes.at(mu) = Tau1 * pow(10., mu - 1);
    }
}


FloatArray
B3SolidMaterial :: computeCharCoefficients(double tPrime, GaussPoint *gp, TimeStep *tStep) const
{
    /*
     * If EmoduliMode == 0 then analysis of continuous retardation spectrum is used for
     * computing characteristic coefficients (moduli) of Kelvin chain
     * Else least-squares method is used
     */
    if ( this->EmoduliMode == 0 ) {
        // constant "n" is assumed to be equal to 0.1 (typical value)

        // modulus of elasticity of the first unit of Kelvin chain.
        // (aging elastic spring with retardation time = 0)
        double lambda0ToPowN = pow(lambda0, 0.1);
        double tau0 = pow(2 * this->giveCharTime(1) / sqrt(10.0), 0.1);
        this->EspringVal = 1. / ( q2 * log(1.0 + tau0 / lambda0ToPowN) - q2 * tau0 / ( 10.0 * lambda0ToPowN + 10.0 * tau0 ) );
        
        // evaluation of moduli of elasticity for the remaining units
        // (Solidifying kelvin units with retardation times tauMu)
        FloatArray answer(nUnits);
        for ( int mu = 1; mu <= this->nUnits; mu++ ) {
            double tauMu = pow(2 * this->giveCharTime(mu), 0.1);
            answer.at(mu) = 10. * pow(1 + tauMu / lambda0ToPowN, 2) / ( log(10.0) * q2 * ( tauMu / lambda0ToPowN ) * ( 0.9 + tauMu / lambda0ToPowN ) );
            this->charTimes.at(mu) *= 1.35;
        }
      
        answer.at(nUnits) /= 1.2;   // modulus of the last unit is reduced
        return answer;

    } else {   // moduli computed using the least-squares method
    
        return KelvinChainMaterial :: computeCharCoefficients(tPrime, gp, tStep);
   
    }
}

double
B3SolidMaterial :: computeCreepFunction(double t, double t_prime, GaussPoint *gp, TimeStep *tStep) const
// compute the value of the creep function of the non-aging solidifying constituent
// corresponding to the given load duration (in days)
// t-t_prime = duration of loading

{
    //double lambda0 = 1.0;     // standard value [day]
    double n = 0.1;     // standard value
    return q2 * log( 1 + pow( (t-t_prime) / lambda0, n) );
}


void
B3SolidMaterial :: giveShrinkageStrainVector(FloatArray &answer,
                                             GaussPoint *gp,
                                             TimeStep *tStep,
                                             ValueModeType mode) const
{
    FloatArray prevAnswer;

    if ( this->shMode == B3_NoShrinkage ) {
        answer.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
        answer.zero();
        return;
    }

    if ( ( mode != VM_Total ) && ( mode != VM_Incremental ) ) {
        OOFEM_ERROR("unsupported mode");
    }

    if ( this->shMode == B3_AverageShrinkage ) {
        this->computeTotalAverageShrinkageStrainVector(answer, gp, tStep);

        if ( ( mode == VM_Incremental ) && ( !tStep->isTheFirstStep() ) ) {
            this->computeTotalAverageShrinkageStrainVector( prevAnswer, gp, tStep->givePreviousStep() );
            answer.subtract(prevAnswer);
        }
    } else {
        this->computePointShrinkageStrainVector(answer, gp, tStep);
    }
}

void
B3SolidMaterial :: computeTotalAverageShrinkageStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) const
{
    /*
     * returns average shrinkage strain vector of cross section at drying
     */

    double TauSh, St, kh, help, E607, Et0Tau, EpsShInf, EpsSh;
    double time = this->relMatAge - this->castingTime + tStep->giveTargetTime() / timeFactor;
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

    StructuralMaterial :: giveReducedSymVectorForm( answer, fullAnswer, gp->giveMaterialMode() );
}


void
B3SolidMaterial :: computePointShrinkageStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) const
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

    humDiff = this->giveHumidityIncrement(gp, tStep);
    EpsSh = humDiff * kSh;

    fullAnswer.resize(size);
    fullAnswer.zero();
    fullAnswer.at(1) = fullAnswer.at(2) = fullAnswer.at(3) = EpsSh;

    StructuralMaterial :: giveReducedSymVectorForm( answer, fullAnswer, gp->giveMaterialMode() );
}


double
B3SolidMaterial :: computeSolidifiedVolume(TimeStep *tStep) const
// compute the relative volume of the solidified material at given age (in days)
{
    double m, alpha;
    double v;     //return value
    double atAge;     // (equivalent age)

    // standard values of exponents - empirical constants
    m = 0.5;
    //lambda0 = 1;     //[day]
    alpha = q3 / q2;

    atAge = this->relMatAge - this->castingTime + ( tStep->giveTargetTime() - 0.5 * tStep->giveTimeIncrement() ) / timeFactor;
    v = 1 / ( alpha + pow(lambda0 / atAge, m) );

    return v;
}


double
B3SolidMaterial :: computeFlowTermViscosity(GaussPoint *gp, TimeStep *tStep) const
//used for evaluation of the flow term viscosity (term containing q4)
{
    if ( this->MicroPrestress == 1 ) {
        double S = this->computeMicroPrestress(gp, tStep, 0); //microprestress in the middle of the time-step
        return 1. / ( q4 * c0 * S );
        //static_cast< B3SolidMaterialStatus* >( gp->giveMaterialStatus() )->setMPS(S);
    } else if ( this->MicroPrestress == 0 ) {
        double tHalfStep = this->relMatAge - this->castingTime + ( tStep->giveTargetTime() - 0.5 * tStep->giveTimeIncrement() ) / timeFactor;
        return 1. * tHalfStep / q4;
    } else {
        OOFEM_ERROR("mode is not supported");
    }
}


void
B3SolidMaterial :: giveEigenStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const
//
// computes the strain due to creep at constant stress during the increment
// (in fact, the INCREMENT of creep strain is computed for mode == VM_Incremental)
//
{
    auto status = static_cast< KelvinChainMaterialStatus * >( this->giveStatus(gp) );

    // !!! chartime exponents are assumed to be equal to 1 !!!
    if ( mode == VM_Incremental ) {
        double v = computeSolidifiedVolume(tStep);
        double eta = this->computeFlowTermViscosity(gp, tStep);         //evaluated too in the middle of the time-step

        FloatArray help, reducedAnswer;
        //        sigma = status->giveStressVector();         //stress vector at the beginning of time-step
        auto sigma = status->giveViscoelasticStressVector();         //stress vector at the beginning of time-step
        FloatMatrix C;
        this->giveUnitComplianceMatrix(C, gp, tStep);
        reducedAnswer.resize( C.giveNumberOfRows() );
        reducedAnswer.zero();

        reducedAnswer.beProductOf(C, sigma);
        reducedAnswer.times( tStep->giveTimeIncrement() / ( timeFactor * eta ) );

        //computes non-aging creep component of the Kelvin Chain
        KelvinChainMaterial :: giveEigenStrainVector(help, gp, tStep, mode);

        help.times(1 / v);
        reducedAnswer.add(help);

        answer = reducedAnswer;
    } else {
        /* error - total mode not implemented yet */
        OOFEM_ERROR("mode is not supported");
    }
}


double
B3SolidMaterial :: inverse_sorption_isotherm(double w) const
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

    /*
     * if ( ( phi < 0.2 ) || ( phi > 0.98 ) ) {
     *  OOFEM_ERROR("Relative humidity h = %e (w=%e) is out of range", phi, w);
     * }
     */

    return phi;
}


double
B3SolidMaterial :: computeMicroPrestress(GaussPoint *gp, TimeStep *tStep, int option) const
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
     * deltaT = tStep->giveTimeIncrement()/timeFactor;
     *
     * if(tStep->isTheFirstStep()){
     *      // if the time step is the first one, ask for an initial microprestress
     *      S = this->giveInitMicroPrestress();
     * } else {
     *      B3SolidMaterialStatus *status = static_cast< B3SolidMaterialStatus * >( this->giveStatus(gp) );
     *      S = status->giveMPS();
     * }
     *
     * humOld = this->giveHumidity(gp, tStep) - this->giveHumidityIncrement(gp, tStep);
     *
     * if (option == 1) {
     *      humNew = this->giveHumidity(gp, tStep);
     * } else if (option == 0) {
     *      deltaT /= 2;
     *      humNew = humOld + 0.5*(this->giveHumidityIncrement(gp, tStep)); //linearly approximated
     * } else {
     *      OOFEM_ERROR("invalid option parameter");
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

    deltaT = tStep->giveTimeIncrement() / timeFactor;
    if ( tStep->isTheFirstStep() ) {
        // if the time step is the first one, ask for an initial microprestress
        S = this->giveInitMicroPrestress();
    } else {
        B3SolidMaterialStatus *status = static_cast< B3SolidMaterialStatus * >( this->giveStatus(gp) );
        S = status->giveMPS();
    }

    humOld = this->giveHumidity(gp, tStep) - this->giveHumidityIncrement(gp, tStep);
    if ( option == 1 ) {  // hum in the end of the time step
        humNew = this->giveHumidity(gp, tStep);
    } else if ( option == 0 ) {  // hum in the middle of the time step
        deltaT /= 2;
        humNew = humOld + 0.5 * ( this->giveHumidityIncrement(gp, tStep) );       //linearly approximated
    } else {
        OOFEM_ERROR("invalid option parameter");
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
B3SolidMaterial :: giveInitMicroPrestress() const
{
    return 1 / ( c0 * tS0 );
}


double
B3SolidMaterial :: giveHumidity(GaussPoint *gp, TimeStep *tStep) const //computes humidity at given TimeStep
{
    /* ask for humidity from external sources, if provided */
    FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();
    FieldPtr tf;

    if ( ( tf = fm->giveField(FT_HumidityConcentration) ) ) {
        // humidity field registered
        FloatArray gcoords, et2;
        int err;
        gp->giveElement()->computeGlobalCoordinates( gcoords, gp->giveNaturalCoordinates() );
        if ( ( err = tf->evaluateAt(et2, gcoords, VM_Total, tStep) ) ) {
            OOFEM_ERROR("tf->evaluateAt failed, error value %d", err);
        }

        // convert water mass to relative humidity
        return this->inverse_sorption_isotherm( et2.at(1) );
    } else {
        OOFEM_ERROR("external fields not found");
        return 0.;
    }

}


double
B3SolidMaterial :: giveHumidityIncrement(GaussPoint *gp, TimeStep *tStep) const //computes humidity increment at given TimeStep
{
    /* ask for humidity from external sources, if provided */
    FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();
    FieldPtr tf;

    if ( ( tf = fm->giveField(FT_HumidityConcentration) ) ) {
        // humidity field registered
        FloatArray gcoords, et2, ei2;
        int err;
        gp->giveElement()->computeGlobalCoordinates( gcoords, gp->giveNaturalCoordinates() );
        if ( ( err = tf->evaluateAt(et2, gcoords, VM_Total, tStep) ) ) {
            OOFEM_ERROR("tf->evaluateAt failed, error value %d", err);
        }

        if ( ( err = tf->evaluateAt(ei2, gcoords, VM_Incremental, tStep) ) ) {
            OOFEM_ERROR("tf->evaluateAt failed, error value %d", err);
        }

        // convert water mass to relative humidity
        return this->inverse_sorption_isotherm( et2.at(1) ) - this->inverse_sorption_isotherm( et2.at(1) - ei2.at(1) );
    } else {
        OOFEM_ERROR("external fields not found");
        return 0.;
    }
}


MaterialStatus *
B3SolidMaterial :: CreateStatus(GaussPoint *gp) const
/*
 * creates a new material status corresponding to this class
 */
{
    return new B3SolidMaterialStatus(gp, nUnits);
}


void
B3SolidMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    KelvinChainMaterial :: giveRealStressVector(answer, gp, reducedStrain, tStep);

    B3SolidMaterialStatus *status = static_cast< B3SolidMaterialStatus * >( this->giveStatus(gp) );
    if ( this->MicroPrestress == 1 ) {
        status->setMPS( this->computeMicroPrestress(gp, tStep, 1) );
    }
}


/****************************************************************************************/
/**********     B3SolidMaterialStatus - HUMIDITY ****************************************/


B3SolidMaterialStatus :: B3SolidMaterialStatus(GaussPoint *g, int nunits) :
    KelvinChainMaterialStatus(g, nunits) { }

void
B3SolidMaterialStatus :: updateYourself(TimeStep *tStep)
{
    microprestress_old = microprestress_new;
    microprestress_new = 0.;

    KelvinChainMaterialStatus :: updateYourself(tStep);
}


void
B3SolidMaterialStatus :: saveContext(DataStream &stream, ContextMode mode)
{
    KelvinChainMaterialStatus :: saveContext(stream, mode);

    if ( !stream.write(microprestress_old) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}

void
B3SolidMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode)
{
    KelvinChainMaterialStatus :: restoreContext(stream, mode);

    if ( !stream.read(microprestress_old) ) {
        THROW_CIOERR(CIO_IOERR);
    }
}
} // end namespace oofem
