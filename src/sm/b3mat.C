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

#include "b3mat.h"
#include "mathfem.h"
#include "gausspnt.h"
#include "structuralcrosssection.h"
#include "timestep.h"

namespace oofem {
IRResultType
B3Material :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    //
    // NOTE
    //
    // this material model is unit dependent!
    // units must be in Mpa, m, MN.
    //
    //
    double fc, c, wc, ac, alpha1 = 0.0, alpha2 = 0.0;
    MaxwellChainMaterial :: initializeFrom(ir);

    int mode = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, mode, IFT_B3Material_mode, "mode");
    if ( mode == 0 ) { // default, determine raw params q1,..,q5 from composition
        IR_GIVE_FIELD(ir, fc, IFT_B3Material_fc, "fc"); // 28-day standard cylinder compression strength in MPa
        IR_GIVE_FIELD(ir, c, IFT_B3Material_cc, "cc"); // cement content of concrete  in kg m^-3.
        IR_GIVE_FIELD(ir, wc, IFT_B3Material_wc, "w/c"); // ratio (by weight) of water to cementitious material
        IR_GIVE_FIELD(ir, ac, IFT_B3Material_ac, "a/c"); // ratio (by weight) of agregate to cement
        IR_GIVE_FIELD(ir, t0, IFT_B3Material_t0, "t0"); // age when drying begins (in days)
    } else { // read raw Basic creep parameters
        IR_GIVE_FIELD(ir, q1, IFT_B3Material_q1, "q1");
        IR_GIVE_FIELD(ir, q2, IFT_B3Material_q2, "q2");
        IR_GIVE_FIELD(ir, q3, IFT_B3Material_q3, "q3");
        IR_GIVE_FIELD(ir, q4, IFT_B3Material_q4, "q4");
    }

    // read shrinkage mode
    int shm = 0;
    IR_GIVE_FIELD(ir, shm, IFT_B3Material_shmode, "shmode");
    this->shMode = ( b3ShModeType ) shm;

    if ( this->shMode == B3_PointShrinkage ) {
        // read additional shrinkage params required
        IR_GIVE_FIELD(ir, es0, IFT_B3Material_es0, "es0");
        IR_GIVE_FIELD(ir, r, IFT_B3Material_r, "r");
        IR_GIVE_FIELD(ir, rprime, IFT_B3Material_rprime, "rprime");
        IR_GIVE_FIELD(ir, at, IFT_B3Material_at, "at");
        // read sorption isotherm data
        IR_GIVE_FIELD(ir, w_h, IFT_B3Material_wh, "w_h"); // Macro
        IR_GIVE_FIELD(ir, n, IFT_B3Material_ncoeff, "ncoeff"); // Macro
        IR_GIVE_FIELD(ir, a, IFT_B3Material_a, "a"); // Macro
    } else if ( this->shMode == B3_AverageShrinkage ) {
        if ( mode == 0 ) { // default mode
            IR_GIVE_FIELD(ir, alpha1, IFT_B3Material_alpha1, "alpha1"); // shrinkage parameter
            IR_GIVE_FIELD(ir, alpha2, IFT_B3Material_alpha2, "alpha2"); // shrinkage parameter
            IR_GIVE_FIELD(ir, ks, IFT_B3Material_ks, "ks"); //cross section shape factor
            /*
             * use ks = 1.0 for an infinite slab
             * = 1.15 for an infinite cylinder
             * = 1.25 for an infinite square prism
             * = 1.30 for a sphere
             * = 1.55 for a cube
             */
            IR_GIVE_FIELD(ir, hum, IFT_B3Material_hum, "hum"); // relative humidity of the environment
            IR_GIVE_FIELD(ir, vs, IFT_B3Material_vs, "vs"); // volume to surface ratio (in m)
        } else { // read raw params
            IR_GIVE_FIELD(ir, kt, IFT_B3Material_kt, "kt"); // shrinkage parameter
            IR_GIVE_FIELD(ir, EpsSinf, IFT_B3Material_EpsSinf, "EpsSinf"); // shrinkage parameter
            IR_GIVE_FIELD(ir, q5, IFT_B3Material_q5, "q5"); // shrinkage parameter
            IR_GIVE_FIELD(ir, vs, IFT_B3Material_vs, "vs"); // volume to surface ratio (in m)
            IR_GIVE_FIELD(ir, ks, IFT_B3Material_ks, "ks"); //cross section shape factor
            IR_GIVE_FIELD(ir, hum, IFT_B3Material_hum, "hum"); // relative humidity of the environment
        }
    }

    IR_GIVE_FIELD(ir, talpha, IFT_B3Material_talpha, "talpha"); // Macro


    w = wc * c;
    E28 = 4734. * sqrt(fc); // or 4733. ?
    if ( mode == 0 ) {
        this->predictParametersFrom(fc, c, wc, ac, t0, alpha1, alpha2);
    }

    return IRRT_OK;
}

void
B3Material :: giveThermalDilatationVector(FloatArray &answer,
                                          GaussPoint *gp,  TimeStep *tStep)
//
// returns strain vector
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
B3Material :: predictParametersFrom(double fc, double c, double wc, double ac,
                                    double t0, double alpha1, double alpha2)
{
    /*
     * Prediction of model parameters - estimation from concrete composition
     * and strength
     *
     * fc   - 28-day standard cylinder compression strength in MPa
     * c    - cement content of concrete  in kg m^-3.
     * wc   - ratio (by weight) of water to cementitious material
     * ac   - ratio (by weight) of agregate to cement
     * t0   - age when drying begins (in days)
     * alpha1, alpha2   - shrinkage parameters
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
     * alpha1   = 1.0  for type I cement;
     *        = 0.85 for type II cement;
     *  = 1.1  for type III cement;
     *
     * alpha2   = 0.75 for steam-cured specimens;
     *        = 1.2  for specimens sealed during curing;
     *  = 1.0  for specimens cured in water or 100% relative humidity.
     *
     */

    // Basic creep parameters

    // q1  = 0.6e6 / E28;
    q1  = 127 * pow(fc, -0.5);
    q2  = 185.4 * pow(c, 0.5) * pow(fc, -0.9);
    q3  = 0.29 * pow(wc, 4.) * q2;
    q4  = 20.3 * pow(ac, -0.7);

    // Shrinkage

    if ( this->shMode == B3_AverageShrinkage ) {
        // the exact value converted from US units would be 85220
        // but this is the SI formula presented in Inelastic Analysis
        kt = 85000 * pow(t0, -0.08) * pow(fc, -0.25);
        EpsSinf = alpha1 * alpha2 * ( 1.9e-2 * pow(w, 2.1) * pow(fc, -0.28) + 270. );

        // Creep at drying

        q5 = 7.57e5 * ( 1. / fc ) * pow(EpsSinf, -0.6);
    }

    char buff [ 1024 ];
    sprintf(buff, "q1=%lf q2=%lf q3=%lf q4=%lf q5=%lf kt=%lf EpsSinf=%lf", q1, q2, q3, q4, q5, kt, EpsSinf);
    OOFEM_LOG_DEBUG("B3mat[%d]: estimated params: %s\n", this->number, buff);
}


double
B3Material :: computeCreepFunction(GaussPoint *gp, double atTime, double ofAge)
{
    // computes the value of creep function at time ofAge
    // when load is acting from atTime
    // WARNING: Area returned by crossSection is assumed to be in [m^2].


    double Qf, Z, r, Q, C0, TauSh, St1, St2, H1, H2, Cd;
    double n, m;


    m = 0.5;
    n = 0.1;

    // basic creep

    Qf = 1. / ( 0.086 * pow(ofAge, 2. / 9.) + 1.21 * pow(ofAge, 4. / 9.) );
    Z  = pow(ofAge, -m) * log( 1. + pow(atTime - ofAge, n) );
    r  = 1.7 * pow(ofAge, 0.12) + 8.0;
    Q  = Qf * pow( ( 1. + pow( ( Qf / Z ), r ) ), -1. / r );

    C0 = q2 * Q + q3 *log( 1. + pow(atTime - ofAge, n) ) + q4 *log(atTime / ofAge);


    if ( this->shMode == B3_AverageShrinkage ) {
        // Aditional creep due to drying

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
    } else {
        Cd = 0.0;
    }

    return 1.e-6 * ( q1 + C0 + Cd );
}



void
B3Material :: giveShrinkageStrainVector(FloatArray &answer,
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
    } else {
        this->computeShrinkageStrainVector(answer, form, gp, atTime, mode);
    }
}

void
B3Material :: computeTotalAverageShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
                                                       GaussPoint *gp, TimeStep *atTime)
{
    /*
     * returns average shrinkage strain vector of cross-section at drying
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
        help = 1. - pow(hum, 3);
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
B3Material :: computeShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
                                           GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
{
    // additional material parameters required:
    //  es0     - final shrinkage at material point
    //  r       - coefficient
    //  rprime  - coefficient
    //  at      - coeff relating stress-induced thermal strain and shrinkage
    double sv, sn, et0, et, wrate = 0.0, trate = 0.0, h1;
    double time = relMatAge + atTime->giveTargetTime() / timeFactor;
    int i, err, tflag = 0, wflag = 0;
    MaxwellChainMaterialStatus *status = ( MaxwellChainMaterialStatus * ) this->giveStatus(gp);
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
    // commpute volumetric stress
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

    //fullAnswer.at(4) = h1*(sn*(r* fullStressVector.at(4)))*(wrate+at*trate);
    //fullAnswer.at(5) = h1*(sn*(r* fullStressVector.at(5)))*(wrate+at*trate);
    //fullAnswer.at(6) = h1*(sn*(r* fullStressVector.at(6)))*(wrate+at*trate);

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
B3Material :: inverse_sorption_isotherm(double w)
{
    // phi ... relative humidity
    // w_h, n, a ... constants obtained from experiments
    double phi;

    // relative humidity
    phi = exp( a * ( 1.0 - pow( ( w_h / w ), ( n ) ) ) );

    if ( ( phi < 0.2 ) || ( phi > 0.98 ) ) {
        _error3("inverse_sorption_isotherm : Relative humidity h = %e (w=%e) is out of range", phi, w);
    }

    return phi;
}
} // end namespace oofem
