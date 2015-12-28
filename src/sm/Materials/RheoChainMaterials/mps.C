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
 *               Copyright (C) 1993 - 2015   Borek Patzak
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

/***************************************************************/
/**************     MPSMaterialStatus     **********************/
/***************************************************************/



MPSMaterialStatus :: MPSMaterialStatus(int n, Domain *d, GaussPoint *gp, int nunits) :
    KelvinChainSolidMaterialStatus(n, d, gp, nunits)
{
    hum = -1.;
    hum_increment = -1.;
    T = -1.;
    T_increment = -1.;
    storedEmodulusFlag = false;
    storedEmodulus = -1.;
    equivalentTime = 0.;
    equivalentTimeTemp = 0.;
    flowTermViscosityTemp = -1.;

#ifdef keep_track_of_strains
    dryingShrinkageStrain = 0.;
    tempDryingShrinkageStrain = 0.;
    autogenousShrinkageStrain = 0.;
    tempAutogenousShrinkageStrain = 0.;

    int rsize = StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() );
    creepStrain.resize(rsize);
    creepStrain.zero();
    creepStrainIncrement.resize(rsize);
    creepStrainIncrement.zero();
#endif
}

void
MPSMaterialStatus :: updateYourself(TimeStep *tStep)
{
    equivalentTime = equivalentTimeTemp;
    equivalentTimeTemp = 0.;

    flowTermViscosity = flowTermViscosityTemp;
    flowTermViscosityTemp = -1.;

    storedEmodulusFlag = false;
    storedEmodulus = -1.;

    hum = -1.;
    hum_increment = -1.;

    T = -1.;
    T_increment = -1.;

#ifdef keep_track_of_strains
    dryingShrinkageStrain = tempDryingShrinkageStrain;
    tempDryingShrinkageStrain = 0.;

    creepStrain.add(creepStrainIncrement);
    creepStrainIncrement.zero();

    autogenousShrinkageStrain = tempAutogenousShrinkageStrain;
    tempAutogenousShrinkageStrain = 0.;
#endif

    KelvinChainSolidMaterialStatus :: updateYourself(tStep);
}


void
MPSMaterialStatus :: initTempStatus()
{
    KelvinChainSolidMaterialStatus :: initTempStatus();

    equivalentTimeTemp = 0.;

    flowTermViscosityTemp = -1.;

    storedEmodulusFlag = false;
    storedEmodulus = -1.;

    hum = -1.;
    hum_increment = -1.;

    T = -1.;
    T_increment = -1.;

#ifdef keep_track_of_strains
    tempDryingShrinkageStrain = 0.;
    creepStrainIncrement.zero();
    tempAutogenousShrinkageStrain = 0.;
#endif

}

contextIOResultType
MPSMaterialStatus :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full information stored in this Status
//
{
    contextIOResultType iores;

    if ( ( iores = KelvinChainSolidMaterialStatus :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    if ( !stream.write(equivalentTime) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(flowTermViscosity) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

contextIOResultType
MPSMaterialStatus :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restore the state variables from a stream
//
{
    contextIOResultType iores;
    if ( ( iores = KelvinChainSolidMaterialStatus :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.read(equivalentTime) ) {
        return CIO_IOERR;
    }

    if ( !stream.read(flowTermViscosity) ) {
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
    this->stiffnessFactor = 1.e6;
    IR_GIVE_OPTIONAL_FIELD(ir, stiffnessFactor, _IFT_MPSMaterial_stiffnessfactor);

    result = KelvinChainSolidMaterial :: initializeFrom(ir);
    if ( result != IRRT_OK ) return result;

    // checking timefactor - for MPS material must be equal to 1.
    double tf;
    IR_GIVE_FIELD(ir, tf, _IFT_RheoChainMaterial_timefactor);
    if ( tf != 1. ) {
        OOFEM_WARNING("for MPS material timefactor must be equal to 1.");
        return IRRT_BAD_FORMAT;
    }

    // initialize exponent p or p_tilde of the governing equation
    if ( ir->hasField(_IFT_MPSMaterial_p_tilde) &&  ir->hasField(_IFT_MPSMaterial_p) ) {
        OOFEM_ERROR("for MPS p and p_tilde cannot be defined at the same time");
    }

    double p_tilde = 2.;
    IR_GIVE_OPTIONAL_FIELD(ir, p_tilde, _IFT_MPSMaterial_p_tilde);

    p = p_tilde / ( p_tilde - 1. );
    IR_GIVE_OPTIONAL_FIELD(ir, p, _IFT_MPSMaterial_p);

    if ( p > 100. ) {
        p = 100.;
    }

    int mode = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, mode, _IFT_MPSMaterial_mode);



    if ( mode == 0 ) { // default - estimate model parameters q1,..,q4 from composition
        IR_GIVE_FIELD(ir, fc, _IFT_MPSMaterial_fc); // 28-day standard cylinder compression strength [MPa]
        IR_GIVE_FIELD(ir,  c, _IFT_MPSMaterial_cc); // cement content of concrete [kg/m^3]
        IR_GIVE_FIELD(ir, wc, _IFT_MPSMaterial_wc); // ratio (by weight) of water to cementitious material
        IR_GIVE_FIELD(ir, ac, _IFT_MPSMaterial_ac); // ratio (by weight) of aggregate to cement

        this->predictParametersFrom(fc, c, wc, ac);
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
    if ( type >= 4 ) {
        OOFEM_WARNING("CoupledAnalysisType must be equal to 0, 1, 2 or 3");
        return IRRT_BAD_FORMAT;
    }

    this->CoupledAnalysis = ( coupledAnalysisType ) type;


    if ( ( this->CoupledAnalysis == MPS_full ) || ( this->CoupledAnalysis == MPS_humidity ) ||  ( this->CoupledAnalysis == MPS_temperature ) ) {
        // age when drying or thermal changes begin [days] - necessary to determine initial viscosity
        IR_GIVE_FIELD(ir, t0, _IFT_MPSMaterial_t0);

        if ( this->p < 100. ) {
            // muS- the only parameter necessary for evaluation of the flow-term viscosity
            // muS = c0 * c1^{p-1} * q4 * (p-1)^p [1/(Pa*s)]
            IR_GIVE_FIELD(ir, muS, _IFT_MPSMaterial_mus);
        } else { // p = 100
            // dimensionless constant for p = infty
            IR_GIVE_FIELD(ir, k3, _IFT_MPSMaterial_k3);
        }
    }

    if ( ( this->CoupledAnalysis == MPS_full ) || ( this->CoupledAnalysis == MPS_humidity ) ) {
        this->kSh = 0.;
        IR_GIVE_OPTIONAL_FIELD(ir, kSh, _IFT_MPSMaterial_ksh);

        this->sh_a = 1.;
        this->sh_hC = 1.;
        this->sh_n = 1.;

        IR_GIVE_OPTIONAL_FIELD(ir, sh_a, _IFT_MPSMaterial_sh_a);
        IR_GIVE_OPTIONAL_FIELD(ir, sh_hC, _IFT_MPSMaterial_sh_hC);
        IR_GIVE_OPTIONAL_FIELD(ir, sh_n, _IFT_MPSMaterial_sh_n);

        this->alphaE = 10.; // according to Bazant 1995
        IR_GIVE_OPTIONAL_FIELD(ir, alphaE, _IFT_MPSMaterial_alphae);
        this->alphaR = 0.1; // according to Bazant 1995
        IR_GIVE_OPTIONAL_FIELD(ir, alphaR, _IFT_MPSMaterial_alphar);
        this->alphaS = 0.1; // according to Bazant 1995
        IR_GIVE_OPTIONAL_FIELD(ir, alphaS, _IFT_MPSMaterial_alphas);


        /// if sortion parameters are provided then it is assumed that the transport analysis uses moisture content/ratio; if not then it exports relative humidity
        // Parameters for desorption isotherm
        //this->w_h =  this->n = this->a = 0.;
        //IR_GIVE_OPTIONAL_FIELD(ir, w_h, _IFT_MPSMaterial_wh);
        //IR_GIVE_OPTIONAL_FIELD(ir, n, _IFT_MPSMaterial_ncoeff);
        //IR_GIVE_OPTIONAL_FIELD(ir, a, _IFT_MPSMaterial_a);
    }

    if ( ( this->CoupledAnalysis == MPS_full ) || ( this->CoupledAnalysis == MPS_temperature ) ) {
        if ( this->CoupledAnalysis == MPS_temperature ) {
            IR_GIVE_FIELD(ir, kTm, _IFT_MPSMaterial_ktm); //[-] replaces ln(h) in diff equation for MPS
        } else {
            this->kTm = -1.;
            IR_GIVE_OPTIONAL_FIELD(ir, kTm, _IFT_MPSMaterial_ktm); //[-] replaces ln(h) in diff equation for MPS
        }

        this->kTc = this->kTm;
        IR_GIVE_OPTIONAL_FIELD(ir, kTc, _IFT_MPSMaterial_ktc);


        this->roomTemperature = 298.15; // reference room temperature for MPS algorithm [K] = 25 C

        /// flag - if true, external fields and reference temperature are in Celsius
        this->temperScaleDifference = 0.;
        if ( ir->hasField(_IFT_MPSMaterial_temperInCelsius) ) {
            this->temperScaleDifference = 273.15;
        }

        // ACTIVATION ENERGIES
        this->QEtoR = 2700.; // [K] value according to Bazant 1995
        IR_GIVE_OPTIONAL_FIELD(ir, QEtoR, _IFT_MPSMaterial_qetor);
        this->QRtoR = 5000.; // [K] value according to Bazant 1995
        IR_GIVE_OPTIONAL_FIELD(ir, QRtoR, _IFT_MPSMaterial_qrtor);
        this->QStoR = 3000.; // [K] value according to Bazant 1995
        IR_GIVE_OPTIONAL_FIELD(ir, QStoR, _IFT_MPSMaterial_qstor);
    }

    // autogenous shrinkage - final amplitude
    // according to NEVILLE: Properties of concrete, typical value is from -40e-6 to -100e-6 for ordinary concretes
    // however for low w/c ratios it can reach up to -700e-6

    // autogenous shrinkage according to fib
    this->eps_cas0 = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, eps_cas0, _IFT_MPSMaterial_eps_cas0);

    if ( ir->hasField(_IFT_MPSMaterial_alpha_as) &&  ir->hasField(_IFT_MPSMaterial_fc) ) {
        double alpha_as;
        // table 5.1-13 from fib2010
        // alpha_as = 800 for cem 32.5N
        // alpha_as = 700 for cem 32.5R and 42.5 N
        // alpha_as = 600 for cem 42.5R, 52.5 N and 52.5 R
        IR_GIVE_FIELD(ir, alpha_as, _IFT_MPSMaterial_alpha_as);
        // mean compressive strenght at the age of 28 days
        IR_GIVE_FIELD(ir, fc, _IFT_MPSMaterial_fc);

        // fib2010 equations 5.1-78
        eps_cas0 = -alpha_as *pow( ( fc / 10. ) / ( 6. + fc / 10. ), 2.5 ) * 1e-6;
    }

    // autogenous shrinkage according to B4
    b4_eps_au_infty = 0.;
    b4_tau_au = 0.;
    b4_alpha = 0.;
    b4_r_t = 0.;
  
    if ( ir->hasField(_IFT_MPSMaterial_B4_cem_type) ) {
        /// auxiliary parameters for autogenous shrinkage according to B4 model
        double b4_r_alpha = 0., b4_eps_au_cem = 0., b4_tau_au_cem = 0., b4_r_ea, b4_r_ew, b4_r_tw;
        int b4_cem_type;

        IR_GIVE_FIELD(ir, b4_cem_type, _IFT_MPSMaterial_B4_cem_type);
        // from Table 2
        if ( b4_cem_type == 0 ) { //R = regular cement
            b4_tau_au_cem = 1.;
            b4_r_tw = 3.;
            b4_r_t = -4.5;
            b4_r_alpha = 1.;
            b4_eps_au_cem = 210.e-6;
            b4_r_ea = -0.75;
            b4_r_ew = -3.5;
        } else if ( b4_cem_type == 1 ) { //RS = rapid hardening
            b4_tau_au_cem = 41.;
            b4_r_tw = 3.;
            b4_r_t = -4.5;
            b4_r_alpha = 1.4;
            b4_eps_au_cem = -84.e-6;
            b4_r_ea = -0.75;
            b4_r_ew = -3.5;
        } else if ( b4_cem_type == 2 ) { //SL = slow hardening
            b4_tau_au_cem = 1.;
            b4_r_tw = 3.;
            b4_r_t = -4.5;
            b4_r_alpha = 1.;
            b4_eps_au_cem = 0.e-6;
            b4_r_ea = -0.75;
            b4_r_ew = -3.5;
        } else {
            OOFEM_ERROR("Unknown cement type (b4_cem_type)");
        }

        IR_GIVE_FIELD(ir, wc, _IFT_MPSMaterial_wc);
        IR_GIVE_FIELD(ir, ac, _IFT_MPSMaterial_ac);

        if ( ir->hasField( _IFT_MPSMaterial_B4_eps_au_infty) ) { 
          IR_GIVE_FIELD(ir, b4_eps_au_infty, _IFT_MPSMaterial_B4_eps_au_infty);
        } else { 
          b4_eps_au_infty = -b4_eps_au_cem *pow(ac / 6., b4_r_ea) *  pow(wc / 0.38, b4_r_ew);
        }

        if ( ir->hasField( _IFT_MPSMaterial_B4_tau_au) ) { 
          IR_GIVE_FIELD(ir, b4_tau_au, _IFT_MPSMaterial_B4_tau_au); // must be in units of the analysis
        } else {
          b4_tau_au = b4_tau_au_cem * pow(wc / 0.38, b4_r_tw);
          b4_tau_au *= lambda0; // converted to desired time unit
        }

        if ( ir->hasField( _IFT_MPSMaterial_B4_alpha) ) {
          IR_GIVE_FIELD(ir, b4_alpha, _IFT_MPSMaterial_B4_alpha);
        } else {
          b4_alpha = b4_r_alpha * wc / 0.38;
        }

        // possibility to overwrite default value of exponent r_t
        IR_GIVE_OPTIONAL_FIELD(ir, b4_r_t, _IFT_MPSMaterial_B4_r_t);
    }

    if ( ( eps_cas0 != 0. ) && ( b4_eps_au_infty != 0. ) ) {
        OOFEM_ERROR("autogenous shrinkage cannot be described according to fib and B4 simultaneously");
    }

    return IRRT_OK;
}


void
MPSMaterial :: giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
{
    KelvinChainSolidMaterial :: giveRealStressVector(answer, gp, reducedStrain, tStep);

    if ( !this->isActivated(tStep) ) {
        return;
    }

    MPSMaterialStatus *status = static_cast< MPSMaterialStatus * >( this->giveStatus(gp) );

    if ( ( this->CoupledAnalysis == MPS_full ) ||  ( this->CoupledAnalysis == MPS_humidity ) || ( this->CoupledAnalysis == MPS_temperature ) ) {
        status->setEquivalentTimeTemp( this->computeEquivalentTime(gp, tStep, 1) );
    }
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
    MaterialMode mMode =  gp->giveMaterialMode();

    answer.resize(6);
    answer.zero();

    if ( ( mMode ==  _2dLattice ) || ( mMode ==  _3dLattice ) ) {
        answer.at(1) = ( talpha );
    } else {
        answer.at(1) = ( talpha );
        answer.at(2) = ( talpha );
        answer.at(3) = ( talpha );
    }
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

#ifndef keep_track_of_strains
    if ( mode == VM_Total ) {
        OOFEM_ERROR("for total formulation of shrinkage strains need to define keep_track_of_strains");
    }
#endif



    MaterialMode mMode = gp->giveMaterialMode();

    FloatArray dry_shr, eps_as;
    double dryShrIncr = 0.;
    double autoShrIncr = 0.;

    answer.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
    answer.zero();

    if ( ( this->CoupledAnalysis == Basic ) || ( this->CoupledAnalysis == MPS_temperature ) || ( !this->isActivated(tStep) ) ) {
        ;
    } else { // compute drying shrinkage
        this->computePointShrinkageStrainVector(dry_shr, gp, tStep);
    }

#ifdef keep_track_of_strains
    MPSMaterialStatus *status = static_cast< MPSMaterialStatus * >( this->giveStatus(gp) );
    if ( dry_shr.giveSize() >= 1 ) {
        dryShrIncr = dry_shr.at(1);
        status->setTempDryingShrinkageStrain(status->giveDryingShrinkageStrain() + dryShrIncr);
    }
#endif

    if ( ( this->eps_cas0 != 0. ) || ( this->b4_eps_au_infty != 0. ) ) {
        if ( this->eps_cas0 != 0. ) {
            this->computeFibAutogenousShrinkageStrainVector(eps_as, gp, tStep);
        } else if ( this->b4_eps_au_infty != 0. ) {
            this->computeB4AutogenousShrinkageStrainVector(eps_as, gp, tStep);
        }

#ifdef keep_track_of_strains
        if ( eps_as.giveSize() >= 1 ) {
            autoShrIncr = eps_as.at(1);
            status->setTempAutogenousShrinkageStrain(status->giveAutogenousShrinkageStrain() + autoShrIncr);
        }
#endif
    }

    if ( mode == VM_Incremental ) {
        if ( dry_shr.giveSize() >= 1 ) {
            answer.add(dry_shr);
        }

        if ( eps_as.giveSize() >= 1 ) {
            answer.add(eps_as);
        }
    } else { //( mode == VM_Total )
        FloatArray fullAnswer;
        int size;

        if ( ( mMode == _3dShell ) || ( mMode ==  _3dBeam ) || ( mMode == _2dPlate ) || ( mMode == _2dBeam ) ) {
            size = 12;
        } else {
            size = 6;
        }

        fullAnswer.resize(size);
        fullAnswer.zero();

        if ( ( mMode ==  _2dLattice ) || ( mMode ==  _3dLattice ) ) {
            fullAnswer.at(1) = status->giveAutogenousShrinkageStrain() + autoShrIncr + status->giveDryingShrinkageStrain() + dryShrIncr;
        } else {
            fullAnswer.at(1) = fullAnswer.at(2) = fullAnswer.at(3) = status->giveAutogenousShrinkageStrain() + autoShrIncr + status->giveDryingShrinkageStrain() + dryShrIncr;
        }

        StructuralMaterial :: giveReducedSymVectorForm( answer, fullAnswer, gp->giveMaterialMode() );
    }

    return;
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
MPSMaterial :: predictParametersFrom(double fc, double c, double wc, double ac)
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
    q1  = 1.e-12 * this->stiffnessFactor * 126.74271 / ( sqrt(fc) );
    q2  = 1.e-12 * this->stiffnessFactor * 185.4 * pow(c, 0.5) * pow(fc, -0.9);
    q3 = 0.29 * pow(wc, 4.) * q2;
    q4  = 1.e-12 * this->stiffnessFactor * 20.3 * pow(ac, -0.7);

    char buff [ 1024 ];
    sprintf(buff, "q1=%lf q2=%lf q3=%lf q4=%lf", q1, q2, q3, q4);
    OOFEM_LOG_DEBUG("MPS[%d]: estimated params: %s\n", this->number, buff);
}



double
MPSMaterial :: computeCreepFunction(double t, double t_prime)
{
    // computes the value of creep function at time t
    // when load is acting from time t_prime
    // t-t_prime = duration of loading

    double Qf, Z, r, Q, C0;
    double n, m;

    m = 0.5;
    n = 0.1;

    // basic creep

    Qf = 1. / ( 0.086 * pow(t_prime, 2. / 9.) + 1.21 * pow(t_prime, 4. / 9.) );
    Z  = pow(t_prime, -m) * log( 1. + pow(t - t_prime, n) );
    r  = 1.7 * pow(t_prime, 0.12) + 8.0;
    Q  = Qf * pow( ( 1. + pow( ( Qf / Z ), r ) ), -1. / r );

    C0 = q2 * Q + q3 *log( 1. + pow(t - t_prime, n) ) + q4 *log(t / t_prime);

    return  ( q1 + C0 );
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
    answer.zero();
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

    double Emodulus;
    MPSMaterialStatus *status = static_cast< MPSMaterialStatus * >( this->giveStatus(gp) );

    if ( !this->isActivated(tStep) ) {
        return 1.; // stresses are cancelled in giveRealStressVector;
    }

    if ( status->giveStoredEmodulusFlag() ) {
        Emodulus = status->giveStoredEmodulus();
    } else {
        if ( EparVal.giveSize() == 0 ) {
            this->updateEparModuli(0.); // stiffnesses are time independent (evaluated at time t = 0.)
        }

        // contribution of the solidifying Kelving chain
        sum = KelvinChainSolidMaterial :: giveEModulus(gp, tStep);

        v = computeSolidifiedVolume(gp, tStep);

        dt = tStep->giveTimeIncrement();
        eta = this->computeFlowTermViscosity(gp, tStep);

        //incremental viscous flow compliance

        if ( this->CoupledAnalysis == Basic ) {
            Cf = 0.5 * ( dt ) / eta;
        } else if ( ( this->CoupledAnalysis == MPS_full ) || ( this->CoupledAnalysis == MPS_humidity ) || ( this->CoupledAnalysis == MPS_temperature ) ) {
            MPSMaterialStatus *status = static_cast< MPSMaterialStatus * >( this->giveStatus(gp) );

            // TRAPEZOIDAL INTEGRATION RULE
            // is the first time step or the material has been just activated (i.e. the previous time was less than casting time)
            if ( tStep->isTheFirstStep() || ( tStep->giveIntrinsicTime() - tStep->giveTimeIncrement() - this->castingTime < 0. ) ) {
                //        if ( tStep->isTheFirstStep() ) {
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

        Emodulus = 1. / ( q1 + 1. / ( EspringVal * v ) + sum  +  Cf );
        status->storeEmodulus(Emodulus);
        status->setEmodulusFlag(true);
    }

    if ( Emodulus < 1. ) {
        OOFEM_ERROR("Incremental modulus is negative %f", Emodulus);
    }

    return Emodulus;
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

    if ( ( this->CoupledAnalysis == MPS_full ) || ( this->CoupledAnalysis == MPS_humidity ) || ( this->CoupledAnalysis == MPS_temperature ) ) {
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

    if ( ( this->CoupledAnalysis == MPS_full ) || ( this->CoupledAnalysis == MPS_humidity ) || ( this->CoupledAnalysis == MPS_temperature ) ) {
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
    double eta = 0., tHalfStep;

    double prevEta, PsiS, e, dt;
    double A = 0., B = 0.;
    double T_new = 0., T_old = 0., H_new = 0., H_old = 0.;
    double kT = 0.;


    /*  double eta, tHalfStep;
     *
     * double prevEta, PsiS, e, dt;
     * double A, B;
     * double T_new, T_old, H_new, H_old;
     * double kT;
     */

    if ( this->CoupledAnalysis == Basic ) {
        tHalfStep = relMatAge + ( tStep->giveTargetTime() - 0.5 * tStep->giveTimeIncrement() );
        eta = tHalfStep / q4;
    } else if ( ( this->CoupledAnalysis == MPS_full ) || ( this->CoupledAnalysis == MPS_humidity ) || ( this->CoupledAnalysis == MPS_temperature ) ) {
        MPSMaterialStatus *status = static_cast< MPSMaterialStatus * >( this->giveStatus(gp) );

        // check whether this viscosity has been already computed
        if ( status->giveFlowTermViscosityTemp() != -1. ) {
            return status->giveFlowTermViscosityTemp();

            // no, the viscosity needs to be evaluated now
        } else {
            // the first time step or the material has been just activated (i.e. the previous time was less than casting time)
            //  ask for an initial value of viscosity
            if ( tStep->isTheFirstStep() || ( tStep->giveIntrinsicTime() - tStep->giveTimeIncrement() - this->castingTime < 0. ) ) {
                prevEta = this->giveInitViscosity(tStep);
            } else {
                // asks for the value of viscosity from the end of the last time-step
                prevEta = status->giveFlowTermViscosity();
            }

            dt = tStep->giveTimeIncrement();
            PsiS = this->computePsiS(gp, tStep); // evaluated in the middle of the time step

            if ( ( this->CoupledAnalysis == MPS_full ) || ( this->CoupledAnalysis == MPS_humidity ) ) {
                H_new = this->giveHumidity(gp, tStep, 1);
                H_old = this->giveHumidity(gp, tStep, 0);
            }

            if ( ( this->CoupledAnalysis == MPS_full ) || ( this->CoupledAnalysis == MPS_temperature ) ) {
                T_new = this->giveTemperature(gp, tStep, 1);
                T_old = this->giveTemperature(gp, tStep, 0);

                if ( ( status->giveTmax() - T_new <= 0. ) || tStep->isTheFirstStep() ||  ( tStep->giveIntrinsicTime() - tStep->giveTimeIncrement() - this->castingTime < 0. ) ) {
                    status->setTmax(T_new);
                    kT = kTm;
                } else {
                    kT = kTc;
                }
            }


            if ( p == 2 ) { // ANALYTICAL SOLUTION FOR "p" = 2
                // evaluate auxiliary factors A and B

                if ( this->CoupledAnalysis == MPS_full ) {
                    if ( kTm == -1. ) {
                        A = sqrt( muS * fabs( ( T_new + T_old ) * ( H_new - H_old ) / ( H_new + H_old ) + log( ( H_new + H_old ) / 2. ) * ( T_new - T_old ) ) / ( dt * this->roomTemperature ) );
                    } else {
                        A = sqrt( muS * fabs( ( T_new + T_old ) * ( H_new - H_old ) / ( H_new + H_old ) - kT * ( T_new - T_old ) ) / ( dt * this->roomTemperature ) );
                    }
                    B = sqrt(PsiS / this->q4);
                } else if ( this->CoupledAnalysis == MPS_humidity ) {
                    // original version of the RHS for constant && room temperature
                    A = sqrt(muS * ( fabs( log(H_new) - log(H_old) ) ) / dt);
                    B = sqrt(PsiS / this->q4);
                } else if ( this->CoupledAnalysis == MPS_temperature ) {
                    A = sqrt( muS * fabs( kT * ( T_new - T_old ) ) / ( dt * this->roomTemperature ) );
                    B = sqrt(PsiS / this->q4);
                } else {
                    OOFEM_ERROR("mode is not supported");
                }

                if ( ( A * B * dt ) > 1.e-6 ) {
                    e = exp(-2 * A * B * dt);
                    eta = ( B / A ) * ( B * ( 1. - e ) + A * prevEta * ( 1. + e ) ) / ( B * ( 1. + e ) + A * prevEta * ( 1. - e ) );
                } else {
                    eta = ( prevEta + B * B * dt ) / ( 1. + A * A * prevEta * dt );
                }
            } else if ( ( p < 100 ) && ( p > 2 ) ) { // SOLUTION FOR VARIOUS "p" EXPLOITING NEXTON'S METHOD
                double iterTol = 1.e-12;
                double deltaEta;
                double deltaDeltaEta;
                double f, df;
                int k;

                if ( this->CoupledAnalysis == MPS_full ) {
                    if ( kTm == -1. ) {
                        A =  pow( muS, 1. / ( p - 1. ) ) * fabs( ( T_new + T_old ) * ( H_new - H_old ) / ( H_new + H_old ) + log( ( H_new + H_old ) / 2. ) * ( T_new - T_old ) ) / ( dt * this->roomTemperature );
                    } else {
                        A =  pow( muS, 1. / ( p - 1. ) ) * fabs( ( T_new + T_old ) * ( H_new - H_old ) / ( H_new + H_old ) - kT * ( T_new - T_old ) ) / ( dt * this->roomTemperature );
                    }
                    B = PsiS / this->q4;
                } else if ( this->CoupledAnalysis == MPS_humidity ) {
                    A =  pow( muS, 1. / ( p - 1. ) ) * ( fabs( log(H_new) - log(H_old) ) ) / dt;
                    B = PsiS / this->q4;
                } else if ( this->CoupledAnalysis == MPS_temperature ) {
                    A =  pow( muS, 1. / ( p - 1. ) ) * fabs( kT * ( T_new - T_old ) ) / ( dt * this->roomTemperature );
                    B = PsiS / this->q4;
                } else {
                    OOFEM_ERROR("mode is not supported");
                }

                deltaEta = 0.;
                deltaDeltaEta = 0.;
                k = 1;
                double relError = 1.;

                while ( relError  > iterTol ) {
                    f = deltaEta / dt + A *pow( prevEta + 0.5 *deltaEta, p / ( p - 1. ) ) - B;
                    df = 1 / dt + A *p *pow( prevEta + 0.5 *deltaEta, 1. / ( p - 1. ) ) / ( 2. * ( p - 1. ) );

                    deltaDeltaEta = -f / df;
                    deltaEta += deltaDeltaEta;

                    if ( k++ > 50 ) {
                        OOFEM_ERROR("iterative Newton method not converging");
                    }

                    relError = fabs(deltaDeltaEta / deltaEta);
                }

                eta = prevEta + deltaEta;
            } else if  ( p < 0. ) { // SOLUTION FOR NEGATIVE "p"
                //            } else if  (p < 0.)  { // SOLUTION FOR NEGATIVE "p"
                double iterTol = 1.e-6;
                double deltaEta;
                double deltaDeltaEta;
                double f, df;
                int k;

                if ( this->CoupledAnalysis == MPS_full ) {
                    if ( kTm == -1. ) {
                        A =  pow( muS, 1. / ( p - 1. ) ) * fabs( T_new * log(H_new) - T_old * log(H_old) ) / ( dt * this->roomTemperature );
                    } else {
                        A =  ( kT * fabs(T_new - T_old) + 0.5 * ( T_new + T_old ) * fabs( log(H_new) - log(H_old) ) ) * pow( muS, 1. / ( p - 1. ) ) /  ( dt * this->roomTemperature );
                    }

                    B = PsiS / this->q4;
                } else if ( this->CoupledAnalysis == MPS_humidity ) {
                    A =  ( fabs( log(H_new) - log(H_old) ) ) * pow( muS, 1. / ( p - 1. ) ) /  dt;
                    B = PsiS / this->q4;
                } else if ( this->CoupledAnalysis == MPS_temperature ) {
                    A =  kT * fabs(T_new - T_old) * pow( muS, 1. / ( p - 1. ) ) /  ( dt * this->roomTemperature );
                    B = PsiS / this->q4;
                } else {
                    OOFEM_ERROR("mode is not supported");
                }

                deltaEta = 0.;
                deltaDeltaEta = 0.;
                k = 1;
                double relError = 1.;
                double minEta = 1.e-6;

                while ( relError  > iterTol ) {
                    f = deltaEta / dt + A *pow( prevEta + deltaEta, p / ( p - 1. ) )  - B;
                    df = 1 / dt + A *pow( prevEta + deltaEta, 1. / ( p - 1. ) ) * p / ( p - 1. );

                    deltaDeltaEta = -f / df;
                    deltaEta += deltaDeltaEta;

                    if ( ( prevEta + deltaEta ) < 0 ) {
                        deltaEta = minEta - prevEta;
                    }

                    if ( k++ > 50 ) {
                        OOFEM_ERROR("iterative Newton method not converging");
                    }


                    relError = fabs(deltaDeltaEta / deltaEta);
                }

                eta = prevEta + deltaEta;
            } else { // ANALYTICAL SOLUTION FOR "p" = INFINITY
                // version for infinite value of parameter p
                // evaluate auxiliary factors A and B
                if ( this->CoupledAnalysis == MPS_full ) {
                    // original version of the RHS for variable/elevated temperature
                    if ( kTm == -1. ) {
                        A =  k3 * fabs( T_new * log(H_new) - T_old * log(H_old) ) / ( dt * this->roomTemperature );
                    } else {
                        A =  k3 * ( kT * fabs(T_new - T_old) + 0.5 * ( T_new + T_old ) * fabs( log(H_new) - log(H_old) ) ) / ( dt * this->roomTemperature );
                    }

                    B = PsiS / this->q4;
                } else if ( this->CoupledAnalysis == MPS_humidity ) {
                    // original version of the RHS for constant && room temperature
                    // compiler warning:
                    A =  k3 * ( fabs( log(H_new) - log(H_old) ) ) / dt;
                    B = PsiS / this->q4;
                } else if ( this->CoupledAnalysis == MPS_temperature ) {
                    A =  k3 * kT * fabs(T_new - T_old) / ( dt * this->roomTemperature );
                    B = PsiS / this->q4;
                } else {
                    OOFEM_ERROR("mode is not supported");
                }

                //if ( A < 1.e-14 ) {
                if ( A < 1.e-10 ) {
                    eta = prevEta + B * dt;
                } else {
                    if ( ( A * dt ) > 30. ) {
                        eta = prevEta;
                    } else {
                        eta = ( prevEta - B / A ) * exp(-A * dt) + B / A;
                    }
                }
            }

            // and now we store into the material status new value of viscosity
            status->setFlowTermViscosityTemp(eta);
        }
    } else {
        OOFEM_ERROR("mode is not supported");
    }

    if ( eta <= 0. ) {
        OOFEM_ERROR("trying to return negative viscosity %f", eta);
    }

    if ( eta != eta ) {
        OOFEM_ERROR("eta is NaN: %f", eta);
    }


    return eta;
}

// returns initial value of the flow term viscosity
double
MPSMaterial :: giveInitViscosity(TimeStep *tStep)
{
    if ( ( t0 - tStep->giveTimeIncrement() ) <= 0. ) {
        OOFEM_ERROR("length of the first time step increment %e must be smaller than t0 %e", tStep->giveTimeIncrement(), t0);
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
        //        sigma = status->giveStressVector();       //stress vector at the beginning of time-step
        sigma = status->giveViscoelasticStressVector();
        this->giveUnitComplianceMatrix(C, gp, tStep);

        reducedAnswer.beProductOf(C, sigma);

        // flow strain increment at constant stress
        dt = tStep->giveTimeIncrement();
        eta = this->computeFlowTermViscosity(gp, tStep);

        if ( this->CoupledAnalysis == Basic ) {
            reducedAnswer.times(dt / eta);
        } else if ( ( this->CoupledAnalysis == MPS_full ) || ( this->CoupledAnalysis == MPS_humidity ) || ( this->CoupledAnalysis == MPS_temperature ) ) {
            // TRAPEZOIDAL INTEGRATION RULE
            //            if ( tStep->isTheFirstStep() ) {
            // is the first time step or the material has been just activated (i.e. the previous time was less than casting time)
            if ( tStep->isTheFirstStep() || ( tStep->giveIntrinsicTime() - tStep->giveTimeIncrement() - this->castingTime < 0. ) ) {
                etaR = this->giveInitViscosity(tStep) / this->computePsiR(gp, tStep, 0);
            } else {
                etaR = status->giveFlowTermViscosity() / this->computePsiR(gp, tStep, 0);
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

#ifdef keep_track_of_strains
        status->setCreepStrainIncrement(answer);
#endif

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
    MaterialMode mMode = gp->giveMaterialMode();
    double h, kShFactor;

    if ( ( mMode == _3dShell ) || ( mMode ==  _3dBeam ) || ( mMode == _2dPlate ) || ( mMode == _2dBeam ) ) {
        size = 12;
    } else {
        size = 6;
    }

    if ( this->sh_a != 1. ) {
        h = this->giveHumidity(gp, tStep, 2);
        kShFactor =  sh_a + ( 1. - sh_a ) / ( 1. + pow( ( 1. - h ) / ( 1. - sh_hC ), sh_n ) );
    } else {
        kShFactor = 1.;
    }

    humDiff = this->giveHumidity(gp, tStep, 3);

    EpsSh = humDiff * kSh * kShFactor;

    fullAnswer.resize(size);
    fullAnswer.zero();

    if ( ( mMode ==  _2dLattice ) || ( mMode ==  _3dLattice ) ) {
        fullAnswer.at(1) = EpsSh;
    } else {
        fullAnswer.at(1) = fullAnswer.at(2) = fullAnswer.at(3) = EpsSh;
    }

    StructuralMaterial :: giveReducedSymVectorForm( answer, fullAnswer, gp->giveMaterialMode() );
}

void
MPSMaterial :: computeFibAutogenousShrinkageStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    /* from fib MC 2010: equations 5.1-76, 5.1-79
     * eps_cas(t) = eps_cas0 * beta_as(t)
     * beta_as(t) = 1 - exp( -0.2 * sqrt(t) )
     *
     * ->> d_eps_cas(t) = eps_cas0 * [ -exp(-0.2 * sqrt(t_e_n+1) ) + exp(-0.2 * sqrt(t_e_n) ) ]
     */

    double eps_cas;

    int size;
    FloatArray fullAnswer;
    MaterialMode mMode = gp->giveMaterialMode();

    double t_equiv_beg, t_equiv_end;


    if ( ( mMode == _3dShell ) || ( mMode ==  _3dBeam ) || ( mMode == _2dPlate ) || ( mMode == _2dBeam ) ) {
        size = 12;
    } else {
        size = 6;
    }


    //    if ( tStep->isTheFirstStep() ) {
    // is the first time step or the material has been just activated (i.e. the previous time was less than casting time)
    if ( tStep->isTheFirstStep() || ( tStep->giveIntrinsicTime() - tStep->giveTimeIncrement() - this->castingTime < 0. ) ) {
        t_equiv_beg = relMatAge - tStep->giveTimeIncrement();
    } else {
        MPSMaterialStatus *status = static_cast< MPSMaterialStatus * >( this->giveStatus(gp) );
        t_equiv_beg = status->giveEquivalentTime();
    }
    // time must be converted into days
    t_equiv_beg /= this->lambda0;

    t_equiv_end = this->computeEquivalentTime(gp, tStep, 1);
    // time must be converted into days
    t_equiv_end /= this->lambda0;

    eps_cas = eps_cas0 * ( -exp( -0.2 * sqrt(t_equiv_end) ) + exp( -0.2 * sqrt(t_equiv_beg) ) );

    fullAnswer.resize(size);
    fullAnswer.zero();

    if ( ( mMode ==  _2dLattice ) || ( mMode ==  _3dLattice ) ) {
        fullAnswer.at(1) = eps_cas;
    } else {
        fullAnswer.at(1) = fullAnswer.at(2) = fullAnswer.at(3) = eps_cas;
    }


    StructuralMaterial :: giveReducedSymVectorForm( answer, fullAnswer, gp->giveMaterialMode() );
}


void
MPSMaterial :: computeB4AutogenousShrinkageStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    /* from Bazant's B4 model (ver. Dec 2014): equations 23-25 and Table 2
     * eps_au(t) = eps_au_infty * [ 1 + (tau_au / t_eq)^alpha ]^r_t
     * alpha = r_alpha * (wc/0.38)
     * eps_au_infty = -eps_au_cem * (ac/6)^r_ea * (wc/0.38)^r_ew
     * tau_au = tau_au_cem * (wc/0.38)^r_tw
     */

    double eps_au;

    int size;
    FloatArray fullAnswer;
    MaterialMode mMode = gp->giveMaterialMode();

    double t_equiv_beg, t_equiv_end;

    if ( ( mMode == _3dShell ) || ( mMode ==  _3dBeam ) || ( mMode == _2dPlate ) || ( mMode == _2dBeam ) ) {
        size = 12;
    } else {
        size = 6;
    }

    //    if ( tStep->isTheFirstStep() ) {

    // is the first time step or the material has been just activated (i.e. the previous time was less than casting time)
    if ( tStep->isTheFirstStep() || ( tStep->giveIntrinsicTime() - tStep->giveTimeIncrement() - this->castingTime < 0. ) ) {
        t_equiv_beg = relMatAge - tStep->giveTimeIncrement();
    } else {
        MPSMaterialStatus *status = static_cast< MPSMaterialStatus * >( this->giveStatus(gp) );
        t_equiv_beg = status->giveEquivalentTime();
    }

    t_equiv_end = this->computeEquivalentTime(gp, tStep, 1);
    //  t_equiv_beg, t_equiv_end and b4_tau_au are in time-units of the analysis
    eps_au = b4_eps_au_infty * ( pow(1. + pow(b4_tau_au / t_equiv_end, b4_alpha), b4_r_t) -  pow(1. + pow(b4_tau_au / t_equiv_beg, b4_alpha), b4_r_t) );

    fullAnswer.resize(size);
    fullAnswer.zero();

    if ( ( mMode ==  _2dLattice ) || ( mMode ==  _3dLattice ) ) {
        fullAnswer.at(1) = eps_au;
    } else {
        fullAnswer.at(1) = fullAnswer.at(2) = fullAnswer.at(3) = eps_au;
    }


    StructuralMaterial :: giveReducedSymVectorForm( answer, fullAnswer, gp->giveMaterialMode() );
}


// double
// MPSMaterial :: inverse_sorption_isotherm(double w)
// Function calculates relative humidity from water content (inverse relation form sorption isotherm).
// Relative humidity (phi) is from range 0.2 - 0.98 !!!
// sorption isotherm by C. R. Pedersen (1990), Combined heat and moisture transfer in building constructions,
//                      PhD-thesis, Technical University of Denmark, Lyngby.
// w (kg/kg) ... water content
// phi ... relative humidity
// w_h, n, a ... constants obtained from experiments
//{
// relative humidity
//   double phi = exp( a * ( 1.0 - pow( ( w_h / w ), ( n ) ) ) );
//
//   /*if ( ( phi < 0.2 ) || ( phi > 0.98 ) ) {
//     * OOFEM_ERROR("Relative humidity h = %e (w=%e) is out of range", phi, w);
//     * }*/
//   //if ( phi < 0.20 ){ phi = 0.2;}
//    //if ( phi > 0.98 ){ phi = 0.98;}
//
//   return phi;
//}

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
            gp->giveElement()->computeGlobalCoordinates( gcoords, gp->giveNaturalCoordinates() );
            if ( ( err = tf->evaluateAt(et2, gcoords, VM_Total, tStep) ) ) {
                OOFEM_ERROR("tf->evaluateAt failed, error value %d", err);
            }

            if ( ( err = tf->evaluateAt(ei2, gcoords, VM_Incremental, tStep) ) ) {
                OOFEM_ERROR("tf->evaluateAt failed, error value %d", err);
            }

            H_tot = et2.at(1);
            H_inc = ei2.at(1);

            // convert water mass to relative humidity
            // H_tot = this->inverse_sorption_isotherm( et2.at(1) );
            // H_inc = this->inverse_sorption_isotherm( et2.at(1) ) - this->inverse_sorption_isotherm( et2.at(1) - ei2.at(1) );
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

    // compute temperature and its increment if the step is first or temperature has not been yet computed
    if ( ( status->giveT() == -1. ) && ( status->giveTIncrement() == -1. ) ) {
        FieldManager *fm = domain->giveEngngModel()->giveContext()->giveFieldManager();

        FM_FieldPtr tf;
        int err, tflag = 0;
        FloatArray gcoords;
        FloatArray et1, ei1; // total and incremental values of temperature

        if ( ( tf = fm->giveField(FT_Temperature) ) ) {
            gp->giveElement()->computeGlobalCoordinates( gcoords, gp->giveNaturalCoordinates() );
            if ( ( err = tf->evaluateAt(et1, gcoords, VM_Total, tStep) ) ) {
                OOFEM_ERROR("tf->evaluateAt failed, error value %d", err);
            }

            if ( ( err = tf->evaluateAt(ei1, gcoords, VM_Incremental, tStep) ) ) {
                OOFEM_ERROR("tf->evaluateAt failed, error value %d", err);
            }

            T_tot = et1.at(1) + temperScaleDifference;
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
    double H, T;
    double betaRH, betaRT;

    if ( this->CoupledAnalysis == MPS_humidity ) {
        H = this->giveHumidity(gp, tStep, option);
        betaRH = alphaR + ( 1. - alphaR ) * H * H;
        return betaRH;
    } else if ( this->CoupledAnalysis == MPS_temperature ) {
        T = this->giveTemperature(gp, tStep, option);
        betaRT = exp( QRtoR * ( 1. / this->roomTemperature - 1. / T ) );
        return betaRT;
    } else if ( this->CoupledAnalysis == MPS_full ) {
        H = this->giveHumidity(gp, tStep, option);
        T = this->giveTemperature(gp, tStep, option);
        betaRH = alphaR + ( 1. - alphaR ) * H * H;
        betaRT = exp( QRtoR * ( 1. / this->roomTemperature - 1. / T ) );
        return betaRH * betaRT;
    } else {
        OOFEM_ERROR("mode is not supported");
        return 0.;
    }
}

double
MPSMaterial :: computePsiS(GaussPoint *gp, TimeStep *tStep)
{
    double AverageHum, AverageTemp;
    double betaSH, betaST;

    if ( this->CoupledAnalysis == MPS_humidity ) {
        AverageHum = this->giveHumidity(gp, tStep, 2);
        betaSH = alphaS + ( 1. - alphaS ) * AverageHum * AverageHum;
        return betaSH;
    } else if ( this->CoupledAnalysis == MPS_temperature ) {
        AverageTemp = this->giveTemperature(gp, tStep, 2);
        betaST = exp( QStoR * ( 1. / this->roomTemperature - 1. /  AverageTemp ) );
        return betaST;
    } else if ( this->CoupledAnalysis == MPS_full ) {
        AverageHum = this->giveHumidity(gp, tStep, 2);
        AverageTemp = this->giveTemperature(gp, tStep, 2);
        betaSH = alphaS + ( 1. - alphaS ) * AverageHum * AverageHum;
        betaST = exp( QStoR * ( 1. / this->roomTemperature - 1. /  AverageTemp ) );
        return betaSH * betaST;
    } else {
        OOFEM_ERROR("mode is not supported");
        return 0.;
    }
}

double
MPSMaterial :: computePsiE(GaussPoint *gp, TimeStep *tStep)
{
    double AverageHum, AverageTemp;
    double betaEH, betaET;

    if ( this->CoupledAnalysis == MPS_humidity ) {
        AverageHum = this->giveHumidity(gp, tStep, 2);
        betaEH = 1. / ( 1. +  pow( ( alphaE * ( 1. - AverageHum ) ), 4. ) );
        return betaEH;
    } else if ( this->CoupledAnalysis == MPS_temperature ) {
        AverageTemp = this->giveTemperature(gp, tStep, 2);
        betaET = exp( QEtoR * ( 1. /  this->roomTemperature - 1. / AverageTemp ) );
        return betaET;
    } else if ( this->CoupledAnalysis == MPS_full ) {
        AverageHum = this->giveHumidity(gp, tStep, 2);
        AverageTemp = this->giveTemperature(gp, tStep, 2);
        betaEH = 1. / ( 1. +  pow( ( alphaE * ( 1. - AverageHum ) ), 4. ) );
        betaET = exp( QEtoR * ( 1. /  this->roomTemperature - 1. / AverageTemp ) );
        return betaEH * betaET;
    } else {
        OOFEM_ERROR(" mode is not supported");
        return 0.;
    }
}


double
MPSMaterial :: computeEquivalentTime(GaussPoint *gp, TimeStep *tStep, int option)
{
    double tEquiv = 0.;
    double PsiE;

    PsiE = computePsiE(gp, tStep);

    // is the first time step or the material has been just activated (i.e. the previous time was less than casting time)
    if ( tStep->isTheFirstStep() || ( tStep->giveIntrinsicTime() - tStep->giveTimeIncrement() - this->castingTime < 0. ) ) {
        if ( option == 0 ) { // gives time in the middle of the timestep
            return relMatAge - tStep->giveTimeIncrement() + PsiE * ( 0.5 * tStep->giveTimeIncrement() );
        } else if ( option == 1 ) { // gives time in the end of the timestep - for UPDATING
            return relMatAge - tStep->giveTimeIncrement() + PsiE *tStep->giveTimeIncrement();
        } else {
            OOFEM_ERROR("mode is not supported")
        }
    } else {
        MPSMaterialStatus *status = static_cast< MPSMaterialStatus * >( this->giveStatus(gp) );
        tEquiv = status->giveEquivalentTime();

        if ( option == 0 ) { // gives time in the middle of the timestep
            tEquiv = tEquiv + PsiE *  0.5 * tStep->giveTimeIncrement();
        } else if ( option == 1 ) { // gives time on the end of the timestep - for UPDATING
            tEquiv = tEquiv + PsiE *tStep->giveTimeIncrement();
        } else {
            OOFEM_ERROR("mode is not supported")
        }

        return tEquiv;
    }
    return tEquiv; // happy compiler
}


int
MPSMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    MPSMaterialStatus *status = static_cast< MPSMaterialStatus * >( this->giveStatus(gp) );
    if ( type == IST_DryingShrinkageTensor ) {
        answer.resize(6);
        answer.zero();
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveDryingShrinkageStrain();
        return 1;
    } else if ( type == IST_AutogenousShrinkageTensor ) {
        answer.resize(6);
        answer.zero();
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveAutogenousShrinkageStrain();
        return 1;
    } else if ( type == IST_TotalShrinkageTensor ) {
        answer.resize(6);
        answer.zero();
        answer.at(1) = answer.at(2) = answer.at(3) = status->giveAutogenousShrinkageStrain() + status->giveDryingShrinkageStrain();
        return 1;
    } else if ( type == IST_CreepStrainTensor ) {
        StructuralMaterial :: giveFullSymVectorForm( answer, status->giveCreepStrain(), gp->giveMaterialMode() );

        return 1;
    } else {
        return RheoChainMaterial :: giveIPValue(answer, gp, type, tStep);
    }

    return 1; // to make the compiler happy
}
} // end namespace oofem
