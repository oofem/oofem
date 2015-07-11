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

#include "hemotkmat.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(HeMoTKMaterial);

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
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, a_0, _IFT_HeMoTKMaterial_a_0);
    IR_GIVE_FIELD(ir, nn, _IFT_HeMoTKMaterial_nn);
    IR_GIVE_FIELD(ir, phi_c, _IFT_HeMoTKMaterial_phi_c);
    IR_GIVE_FIELD(ir, delta_wet, _IFT_HeMoTKMaterial_delta_wet);

    IR_GIVE_FIELD(ir, w_h, _IFT_HeMoTKMaterial_w_h);
    IR_GIVE_FIELD(ir, n, _IFT_HeMoTKMaterial_n);
    IR_GIVE_FIELD(ir, a, _IFT_HeMoTKMaterial_a);

    IR_GIVE_FIELD(ir, latent, _IFT_HeMoTKMaterial_latent);
    IR_GIVE_FIELD(ir, c, _IFT_HeMoTKMaterial_c);
    IR_GIVE_FIELD(ir, rho, _IFT_HeMoTKMaterial_rho);
    IR_GIVE_FIELD(ir, chi_eff, _IFT_HeMoTKMaterial_chi_eff);

    IR_GIVE_FIELD(ir, por, _IFT_HeMoTKMaterial_por);
    IR_GIVE_FIELD(ir, rho_gws, _IFT_HeMoTKMaterial_rho_gws);

    return Material :: initializeFrom(ir);
}


double
HeMoTKMaterial :: give(int aProperty, GaussPoint *gp)
//
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
//
{
    return Material :: give(aProperty, gp);
}


void
HeMoTKMaterial :: giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep)
{
    TransportMaterialStatus *ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );

    double w = field.at(2);
    double t = field.at(1);

    FloatArray ans_w, ans_t;
    FloatArray grad_w, grad_t;
    int size = grad.giveSize() / 2;
    for ( int i = 1; i <= size; ++i ) {
        grad_w.at(i) = grad.at(i);
    }
    for ( int i = size + 1; i <= size * 2; ++i ) {
        grad_t.at(i) = grad.at(i);
    }
    ans_w.beScaled(perm_ww(w, t), grad_w);
    ans_w.beScaled(perm_wt(w, t), grad_t);
    ans_t.beScaled(perm_ww(w, t) * get_latent(w, t), grad_w);
    ans_t.beScaled(get_chi(w, t) + get_latent(w, t) * perm_wt(w, t), grad_t);

    answer.resize(size * 2);
    answer.zero();
    answer.addSubVector(ans_w, 1);
    answer.addSubVector(ans_t, size + 1);

    ms->setTempField(field);
    ms->setTempGradient(grad);
    ms->setTempFlux(answer);
}


void
HeMoTKMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *tStep)
{
    /*
     * returns constitutive matrix of receiver
     */
    if ( ( mode == Conductivity_ww ) || ( mode == Conductivity_hh ) || ( mode == Conductivity_hw ) || ( mode == Conductivity_wh ) ) {
        this->computeConductivityMtrx(answer, mode, gp, tStep);
    } else {
        OOFEM_ERROR("unknown mode (%s)", __MatResponseModeToString(mode) );
    }
}


double
HeMoTKMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *tStep)
{
    return this->computeCapacityCoeff(mode, gp, tStep);
}


void HeMoTKMaterial :: computeConductivityMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mmode = gp->giveMaterialMode();
    switch ( mmode ) {
    case _2dHeMo:
        this->matcond2d(answer, gp, mode, tStep);
        return;

    case _3dHeMo:
        this->matcond3d(answer, gp, mode, tStep);
        return;

    default:
        OOFEM_ERROR("Unsupported MaterialMode");
    }
}


void
HeMoTKMaterial :: matcond1d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *tStep)
//  function creates conductivity matrix of the
//  isotropic heat material for 1D problems
//
//  d - conductivity matrix of the material
//  25.9.2001
{
    double k = 0.0, w = 0.0, t = 0.0;
    TransportMaterialStatus *status = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    FloatArray s;


    //  w = Tm->ip[ipp].av[0];
    //  t = Tm->ip[ipp].av[1];
    s = status->giveTempField();
    if ( s.isEmpty() ) {
        OOFEM_ERROR("undefined state vector");
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
        OOFEM_ERROR("Unknown MatResponseMode");
    }

    d.resize(1, 1);
    d.at(1, 1) = k;
}

void
HeMoTKMaterial :: matcond2d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *tStep)
//  function creates conductivity matrix of the
//  isotropic heat material for 2D problems
//
//  d - conductivity matrix of the material
//  25.9.2001
{
    double k = 0.0, w = 0.0, t = 0.0;
    TransportMaterialStatus *status = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    FloatArray s;


    //  w = Tm->ip[ipp].av[0];
    //  t = Tm->ip[ipp].av[1];
    s = status->giveTempField();
    if ( s.isEmpty() ) {
        OOFEM_ERROR("undefined state vector");
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
        OOFEM_ERROR("Unknown MatResponseMode");
    }

    d.resize(2, 2);
    d.at(1, 1) = k;
    d.at(1, 2) = 0.0;
    d.at(2, 1) = 0.0;
    d.at(2, 2) = k;
}

void
HeMoTKMaterial :: matcond3d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *tStep)
//  function creates conductivity matrix of the
//  isotropic heat material for 3D problems
//
//  d - conductivity matrix of the material
//  25.9.2001
{
    double k = 0.0, w = 0.0, t = 0.0;
    TransportMaterialStatus *status = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    FloatArray s;


    //  w = Tm->ip[ipp].av[0];
    //  t = Tm->ip[ipp].av[1];
    s = status->giveTempField();
    if ( s.isEmpty() ) {
        OOFEM_ERROR("undefined state vector");
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
        OOFEM_ERROR("Unknown MatResponseMode");
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


double HeMoTKMaterial :: computeCapacityCoeff(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    if ( mode == Capacity_ww ) {
        return 1.0 * rho;
    } else if ( mode == Capacity_wh ) {
        return 0.0;
    } else if ( mode == Capacity_hw ) {
        TransportMaterialStatus *status = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
        FloatArray s;
        double w, t;

        s = status->giveTempField();
        if ( s.isEmpty() ) {
            OOFEM_ERROR("undefined state vector");
        }

        w = s.at(2);
        t = s.at(1);
        return get_b(w, t) * get_latent(w, t);
    } else if ( mode == Capacity_hh ) {
        TransportMaterialStatus *status = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
        FloatArray s;
        double w, t;

        s = status->giveTempField();
        if ( s.isEmpty() ) {
            OOFEM_ERROR("undefined state vector");
        }

        w = s.at(2);
        t = s.at(1);
        return get_ceff(w, t);
    } else {
        OOFEM_ERROR("Unknown MatResponseMode");
    }

    return 0.0; // to make compiler happy
}


double
HeMoTKMaterial :: giveHumidity(GaussPoint *gp, ValueModeType mode)
{
    TransportMaterialStatus *ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    const FloatArray &tempState = ms->giveTempField();
    if ( tempState.giveSize() < 2 ) {
        OOFEM_ERROR("undefined moisture status!");
    }

    const FloatArray &state = ms->giveField();

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
        OOFEM_ERROR("Relative humidity is out of range");
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
        OOFEM_ERROR("Relative humidity %.3f is out of range", phi);
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
        OOFEM_ERROR("Relative humidity %.3f is out of range", phi);
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
        OOFEM_ERROR("unknown mode (%s)", __MatResponseModeToString(mode) );
    }

    return false; // to make compiler happy
}

int
HeMoTKMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
// IST_Humidity overriden to use inverse_sorption_isotherm
{
    if ( type == IST_Humidity ) {
        answer.resize(1);
        answer.at(1) = giveHumidity(gp, VM_Velocity); // VM_Previous = equilibrated value of humidity
        return 1;
    } else {
        return TransportMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}
} // end namespace oofem
