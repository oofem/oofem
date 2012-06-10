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

#include "hemotkmat.h"
#include "flotmtrx.h"
#include "gausspnt.h"
#include "mathfem.h"

namespace oofem {
int
HeMoTKMaterial :: hasMaterialModeCapability(MaterialMode mode)
{
    if ( ( mode == _2dHeMo ) || ( mode == _3dHeMo ) ) {
        return 1;
    }

    return 0;
}


IRResultType
HeMoTKMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    // double value ;

    this->Material :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, a_0, IFT_HeMoTKMaterial_a_0, "a_0"); // Macro
    IR_GIVE_FIELD(ir, nn, IFT_HeMoTKMaterial_nn, "nn"); // Macro
    IR_GIVE_FIELD(ir, phi_c, IFT_HeMoTKMaterial_phi_c, "phi_c"); // Macro
    IR_GIVE_FIELD(ir, delta_wet, IFT_HeMoTKMaterial_delta_wet, "delta_wet"); // Macro

    IR_GIVE_FIELD(ir, w_h, IFT_HeMoTKMaterial_w_h, "w_h"); // Macro
    IR_GIVE_FIELD(ir, n, IFT_HeMoTKMaterial_n, "n"); // Macro
    IR_GIVE_FIELD(ir, a, IFT_HeMoTKMaterial_a, "a"); // Macro

    IR_GIVE_FIELD(ir, latent, IFT_HeMoTKMaterial_latent, "latent"); // Macro
    IR_GIVE_FIELD(ir, c, IFT_HeMoTKMaterial_c, "c"); // Macro
    IR_GIVE_FIELD(ir, rho, IFT_HeMoTKMaterial_rho, "rho"); // Macro
    IR_GIVE_FIELD(ir, chi_eff, IFT_HeMoTKMaterial_chi_eff, "chi_eff"); // Macro

    IR_GIVE_FIELD(ir, por, IFT_HeMoTKMaterial_por, "por"); // Macro
    IR_GIVE_FIELD(ir, rho_gws, IFT_HeMoTKMaterial_rho_gws, "rho_gws"); // Macro

    return IRRT_OK;
}


double
HeMoTKMaterial :: give(int aProperty, GaussPoint *gp)
//
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
//
{
    return this->Material :: give(aProperty, gp);
}



void
HeMoTKMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime)
{
    /*
     * returns constitutive matrix of receiver
     */
    if ( ( mode == Conductivity_ww ) || ( mode == Conductivity_hh ) || ( mode == Conductivity_hw ) || ( mode == Conductivity_wh ) ) {
        this->computeConductivityMtrx(answer, mode, gp, atTime);
    } else {
        _error2( "giveCharacteristicMatrix : unknown mode (%s)", __MatResponseModeToString(mode) );
    }
}


double
HeMoTKMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *atTime)
{
    return this->computeCapacityCoeff(mode, gp, atTime);
}


void HeMoTKMaterial :: computeConductivityMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    MaterialMode mmode = gp->giveMaterialMode();
    switch ( mmode ) {
    case _2dHeMo:
        this->matcond2d(answer, gp, mode, atTime);
        return;

    case _3dHeMo:
        this->matcond3d(answer, gp, mode, atTime);
        return;

    default:
        _error("Unsupported MaterialMode");
    }
}


void
HeMoTKMaterial :: matcond1d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime)
//  function creates conductivity matrix of the
//  isotropic heat material for 1D problems
//
//  d - conductivity matrix of the material
//  25.9.2001
{
    double k = 0.0, w = 0.0, t = 0.0;
    TransportMaterialStatus *status = ( TransportMaterialStatus * ) this->giveStatus(gp);
    FloatArray s;


    //  w = Tm->ip[ipp].av[0];
    //  t = Tm->ip[ipp].av[1];
    s = status->giveTempStateVector();
    if ( s.isEmpty() ) {
        _error("matcond1d: undefined state vector");
    }

    w = s.at(2);
    t = s.at(1);

    if ( mode == Conductivity_ww ) {
        k = perm_ww(w, t);
    } else if ( mode == Conductivity_wh ) {
        k = perm_wt(w, t);
    } else if ( mode == Conductivity_hw ) {
        k = perm_ww(w, t) * get_latent(w, t);
    } else if ( mode == Conductivity_hh ) {
        k = get_chi(w, t) + get_latent(w, t) * perm_wt(w, t);
    } else {
        _error("Unknown MatResponseMode");
    }

    d.resize(1, 1);
    d.at(1, 1) = k;
}

void
HeMoTKMaterial :: matcond2d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime)
//  function creates conductivity matrix of the
//  isotropic heat material for 2D problems
//
//  d - conductivity matrix of the material
//  25.9.2001
{
    double k = 0.0, w = 0.0, t = 0.0;
    TransportMaterialStatus *status = ( TransportMaterialStatus * ) this->giveStatus(gp);
    FloatArray s;


    //  w = Tm->ip[ipp].av[0];
    //  t = Tm->ip[ipp].av[1];
    s = status->giveTempStateVector();
    if ( s.isEmpty() ) {
        _error("matcond2d: undefined state vector");
    }

    w = s.at(2);
    t = s.at(1);

    if ( mode == Conductivity_ww ) {
        k = perm_ww(w, t);
    } else if ( mode == Conductivity_wh ) {
        k = perm_wt(w, t);
    } else if ( mode == Conductivity_hw ) {
        k = perm_ww(w, t) * get_latent(w, t);
    } else if ( mode == Conductivity_hh ) {
        k = get_chi(w, t) + get_latent(w, t) * perm_wt(w, t);
    } else {
        _error("Unknown MatResponseMode");
    }

    d.resize(2, 2);
    d.at(1, 1) = k;
    d.at(1, 2) = 0.0;
    d.at(2, 1) = 0.0;
    d.at(2, 2) = k;
}

void
HeMoTKMaterial :: matcond3d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime)
//  function creates conductivity matrix of the
//  isotropic heat material for 3D problems
//
//  d - conductivity matrix of the material
//  25.9.2001
{
    double k = 0.0, w = 0.0, t = 0.0;
    TransportMaterialStatus *status = ( TransportMaterialStatus * ) this->giveStatus(gp);
    FloatArray s;


    //  w = Tm->ip[ipp].av[0];
    //  t = Tm->ip[ipp].av[1];
    s = status->giveTempStateVector();
    if ( s.isEmpty() ) {
        _error("matcond3d: undefined state vector");
    }

    w = s.at(2);
    t = s.at(1);

    if ( mode == Conductivity_ww ) {
        k = perm_ww(w, t);
    } else if ( mode == Conductivity_wh ) {
        k = perm_wt(w, t);
    } else if ( mode == Conductivity_hw ) {
        k = perm_ww(w, t) * get_latent(w, t);
    } else if ( mode == Conductivity_hh ) {
        k = get_chi(w, t) + get_latent(w, t) * perm_wt(w, t);
    } else {
        _error("Unknown MatResponseMode");
    }

    d.resize(3, 3);
    d.at(1, 1) = k;
    d.at(1, 2) = 0.0;
    d.at(1, 3) = 0.0;
    d.at(2, 1) = 0.0;
    d.at(2, 2) = k;
    d.at(2, 3) = 0.0;
    d.at(3, 1) = 0.0;
    d.at(3, 2) = 0.0;
    d.at(3, 3) = k;
}


double HeMoTKMaterial :: computeCapacityCoeff(MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    if ( mode == Capacity_ww ) {
        return 1.0 * rho;
    } else if ( mode == Capacity_wh ) {
        return 0.0;
    } else if ( mode == Capacity_hw ) {
        TransportMaterialStatus *status = ( TransportMaterialStatus * ) this->giveStatus(gp);
        FloatArray s;
        double w, t;

        s = status->giveTempStateVector();
        if ( s.isEmpty() ) {
            _error("computeCapacityCoeff: undefined state vector");
        }

        w = s.at(2);
        t = s.at(1);
        return get_b(w, t) * get_latent(w, t);
    } else if ( mode == Capacity_hh ) {
        TransportMaterialStatus *status = ( TransportMaterialStatus * ) this->giveStatus(gp);
        FloatArray s;
        double w, t;

        s = status->giveTempStateVector();
        if ( s.isEmpty() ) {
            _error("computeCapacityCoeff: undefined state vector");
        }

        w = s.at(2);
        t = s.at(1);
        return get_ceff(w, t);
    } else {
        _error("Unknown MatResponseMode");
    }

    return 0.0; // to make compiler happy
}


double
HeMoTKMaterial :: giveHumidity(GaussPoint *gp, ValueModeType mode)
{
    FloatArray tempState = ( ( TransportMaterialStatus * ) giveStatus(gp) )->giveTempStateVector();
    if ( tempState.giveSize() < 2 ) {
        _error("giveHumidity: undefined moisture status!");
    }

    FloatArray state = ( ( TransportMaterialStatus * ) giveStatus(gp) )->giveStateVector();

    if ( mode == VM_Total ) {
        return inverse_sorption_isotherm( tempState.at(2) );
    } else if ( mode == VM_Incremental ) {
        return inverse_sorption_isotherm( tempState.at(2) ) - inverse_sorption_isotherm( state.at(2) );
    } else if ( mode == VM_Velocity ) { // VM_Previous
        return inverse_sorption_isotherm( state.at(2) );
    }

    return 1.;
}




double
HeMoTKMaterial :: perm_ww(double w, double t)
{
    // Function calculates permability water content-water content (k_ww)
    // phi ... relative humidity
    // delta_gw ... water vapor permeability
    // dphi_dw ... differentiation of relative with respect to water content
    // p_gws ... saturation water vapor pressure
    double k_ww, phi, delta_gw, dphi_dw, p_gws;

    phi = inverse_sorption_isotherm(w);
    delta_gw = give_delta_gw(phi);
    dphi_dw = give_dphi_dw(w);
    p_gws = give_p_gws(t);

    k_ww = delta_gw * p_gws * dphi_dw;

    return ( k_ww );
}

double
HeMoTKMaterial :: perm_wt(double w, double t)
{
    // Function calculates permability water content-temperature (k_wt)
    // delta_gw ... water vapor permeability
    // d_pgw_d_t ... differentiation of water vapor presssure with respect to temperature
    double k_wt, phi, delta_gw, dpgw_dt;

    phi = inverse_sorption_isotherm(w);
    delta_gw = give_delta_gw(phi);
    dpgw_dt = give_dpgw_dt(t, phi);

    k_wt = delta_gw * dpgw_dt;

    return ( k_wt );
}


double
HeMoTKMaterial :: give_delta_gw(double phi)
// Function calculates the water vapor permeability delta_gw. Relative humidity (phi) is from range 0.2 - 0.98 !!!
// delta_gw ... by Z. P. Bazant and L. J. Najjar (1972), Nonlinear water diffusion in nonsaturated concrete,
//              MATERIAUX ET CONSTRUCTIONS, Vol. 5, No. 25, pp. 3 -20.
// phi ... relative humidity
// a_0, nn, phi_c, delta_wet ... constants obtained from experiments
{
    double delta_gw;

    if ( ( phi < 0.2 ) || ( phi > 0.98 ) ) {
        _error("give_delta_gw : Relative humidity is out of range");
    }

    // water vapor permeability
    delta_gw = delta_wet * ( a_0 + ( 1.0 - a_0 ) / ( 1.0 + pow( ( 1.0 - phi ) / ( 1.0 - phi_c ), nn ) ) );

    return ( delta_gw );
}



double
HeMoTKMaterial :: sorption_isotherm(double phi)
// Function calculates water content in medium from relative humidity.  Relative humidity (phi) is from range 0.2 - 0.98 !!!
// sorption isotherm by C. R. Pedersen (1990), Combined heat and moisture transfer in building constructions,
//                      PhD-thesis, Technical University of Denmark, Lingby.
// w (kg/kg) ... water content
// phi ... relative humidity
// w_h, n, a ... constants obtained from experiments
{
    double w;

    if ( ( phi < 0.2 ) || ( phi > 0.98 ) ) {
        OOFEM_ERROR2("HeMoTKMaterial :: sorption_isotherm : Relative humidity %.3f is out of range", phi);
    }

    // water content
    w = w_h * pow( ( 1.0 - log(phi) / a ), ( -1.0 / n ) );

    return ( w );
}



double
HeMoTKMaterial :: inverse_sorption_isotherm(double w)
// Function calculates relative humidity from water content (inverse relation form sorption isotherm).
// Relative humidity (phi) is from range 0.2 - 0.98 !!!
// sorption isotherm by C. R. Pedersen (1990), Combined heat and moisture transfer in building constructions,
//                      PhD-thesis, Technical University of Denmark, Lingby.
// w (kg/kg) ... water content
// phi ... relative humidity
// w_h, n, a ... constants obtained from experiments
{
    double phi;

    // relative humidity
    phi = exp( a * ( 1.0 - pow( ( w_h / w ), ( n ) ) ) );

    if ( ( phi < 0.2 ) || ( phi > 0.98 ) ) {
        OOFEM_ERROR2("HeMoTKMaterial :: inverse_sorption_isotherm : Relative humidity %.3f is out of range", phi);
    }

    return ( phi );
}



double
HeMoTKMaterial :: give_dphi_dw(double w)
// Function calculates differentiation of relative humidity with respect to water content (inverse relation form sorption isotherm).
// Relative humidity (phi) is from range 0.2 - 0.98 !!!
// sorption isotherm by C. R. Pedersen (1990), Combined heat and moisture transfer in building constructions,
//                      PhD-thesis, Technical University of Denmark, Lingby.
// w (kg/kg) ... water content
// phi ... relative humidity
// w_h, n, a ... constants obtained from experiments
{
    double dphi_dw;

    // diferentiation of realative humidity with respect to water content
    dphi_dw = exp( a * ( 1.0 - pow( ( w_h / w ), n ) ) ) * a * n * pow(w_h, n) * pow( w, ( -1.0 - n ) );

    return ( dphi_dw );
}

double
HeMoTKMaterial :: give_dpgw_dt(double t, double phi)
// Function calculates differentiation of water vapor pressure with respect to temperature
// saturation water vapor pressure by C. R. Pedersen (1990), Combined heat and moisture transfer in building constructions,
//                                 PhD-thesis, Technical University of Denmark, Lingby.
// t ... temperature [K]
// phi ... relative humidity
{
    double dp_gw_dt, dp_gws_dt;

    //differentiation of saturation water vapor pressure (analytical expression) with respect to temperature
    dp_gws_dt = exp( 23.5771 - 4042.9 / ( t - 37.58 ) ) * 4042.9 / ( t - 37.58 ) / ( t - 37.58 );

    dp_gw_dt = phi * dp_gws_dt;

    return ( dp_gw_dt );
}


double
HeMoTKMaterial :: give_p_gws(double t)
// Function calculates saturation water vapor pressure
// saturation water vapor pressure by C. R. Pedersen (1990), Combined heat and moisture transfer in building constructions,
//                                 PhD-thesis, Technical University of Denmark, Lingby.
// t ... temperature [K]
// phi ... relative humidity
{
    double p_gws;

    //saturation water vapor pressure ... analytical expression
    p_gws = exp( 23.5771 - 4042.9 / ( t - 37.58 ) );

    return ( p_gws );
}

double
HeMoTKMaterial :: get_latent(double w, double t)
{
    // Function calculates latent heat of evaporation

    double l;

    l = latent; //zatim!!!!

    return ( l );
}

double
HeMoTKMaterial :: get_b(double w, double t)
{
    // Function calculates coefficient b
    // sat ... degree of saturation

    double b, sat, phi, dphi_dw;

    phi = inverse_sorption_isotherm(w);
    dphi_dw = give_dphi_dw(w);
    sat = get_sat(w, t);

    b = por * rho_gws * ( phi + ( 1 - sat ) * dphi_dw );

    return ( b );
}


double
HeMoTKMaterial :: get_sat(double w, double t)
{
    // Function calculates degree of saturation
    // sat ... degree of saturation

    double sat;

    sat = 1.0; //zatim!!!!

    return ( sat );
}

double
HeMoTKMaterial :: get_ceff(double w, double t)
{
    // Function calculates effective thermal capacity

    double ceff;

    ceff = c * rho; //zatim!!!!

    return ( ceff );
}


double
HeMoTKMaterial :: get_chi(double w, double t)
{
    // Function calculates effective thermal conductivity

    double chi;

    chi = chi_eff; //zatim!!!!

    return ( chi );
}

bool
HeMoTKMaterial :: isCharacteristicMtrxSymmetric(MatResponseMode mode)
{
    if ( ( mode == Conductivity_ww ) || ( mode == Conductivity_hh ) || ( mode == Conductivity_hw ) || ( mode == Conductivity_wh ) ) {
        return false;
    } else {
        _error2( "isCharacteristicMtrxSymmetric : unknown mode (%s)", __MatResponseModeToString(mode) );
    }

    return false; // to make compiler happy
}

int
HeMoTKMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
// IST_Humidity overriden to use inverse_sorption_isotherm
{
    if ( type == IST_Humidity ) {
        answer.resize(1);
        answer.at(1) = giveHumidity(aGaussPoint, VM_Velocity); // VM_Previous = equilibrated value of humidity
        return 1;
    } else {
        return TransportMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}
} // end namespace oofem
