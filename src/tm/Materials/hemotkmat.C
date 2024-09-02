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

#include "tm/Materials/hemotkmat.h"
#include "floatmatrixf.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(HeMoTKMaterial);

bool
HeMoTKMaterial :: hasMaterialModeCapability(MaterialMode mode) const
{
    return mode == _2dHeMo || mode == _3dHeMo;
}


void
HeMoTKMaterial :: initializeFrom(InputRecord &ir)
{
    Material :: initializeFrom(ir);

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
}


double
HeMoTKMaterial :: give(int aProperty, GaussPoint *gp) const
{
    return Material :: give(aProperty, gp);
}


std::pair<FloatArrayF<3>, FloatArrayF<3>>
HeMoTKMaterial :: computeHeMoFlux3D(const FloatArrayF<3> &grad_t, const FloatArrayF<3> &grad_w, double t, double w, GaussPoint *gp, TimeStep *tStep) const
{
    auto ms = static_cast< HeMoTransportMaterialStatus * >( this->giveStatus(gp) );

    auto ans_w = perm_ww(w, t) * grad_w + perm_wt(w, t) * grad_t;
    auto ans_t = perm_ww(w, t) * get_latent(w, t) * grad_w + (get_chi(w, t) + get_latent(w, t) * perm_wt(w, t)) * grad_t;

    ms->setTempTemperature(t);
    ms->setTempTemperatureGradient(grad_t);
    ms->setTempHeatFlux(ans_t);

    ms->setTempHumidity(w);
    ms->setTempHumidityGradient(grad_w);
    ms->setTempHumidityFlux(ans_w);

    return {ans_t, ans_w};
}

FloatMatrixF<3,3>
HeMoTKMaterial :: computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< HeMoTransportMaterialStatus * >( this->giveStatus(gp) );

    double t = status->giveTempTemperature();
    double w = status->giveTempHumidity();

    double k = 0.0;
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

    return k * eye<3>();
}


double
HeMoTKMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *tStep) const
{
    return this->computeCapacityCoeff(mode, gp, tStep);
}


double HeMoTKMaterial :: computeCapacityCoeff(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    if ( mode == Capacity_ww ) {
        return 1.0 * rho;
    } else if ( mode == Capacity_wh ) {
        return 0.0;
    } else if ( mode == Capacity_hw ) {
        auto status = static_cast< HeMoTransportMaterialStatus * >( this->giveStatus(gp) );

        double t = status->giveTempTemperature();
        double w = status->giveTempHumidity();
        return get_b(w, t) * get_latent(w, t);
    } else if ( mode == Capacity_hh ) {
        auto status = static_cast< HeMoTransportMaterialStatus * >( this->giveStatus(gp) );

        double t = status->giveTempTemperature();
        double w = status->giveTempHumidity();
        return get_ceff(w, t);
    } else {
        OOFEM_ERROR("Unknown MatResponseMode");
    }

    // return 0.0; // to make compiler happy
}


double
HeMoTKMaterial :: giveHumidity(GaussPoint *gp, ValueModeType mode) const
{
    auto status = static_cast< HeMoTransportMaterialStatus * >( this->giveStatus(gp) );

    double w = status->giveHumidity();
    double tempw = status->giveTempHumidity();

    if ( mode == VM_Total ) {
        return inverse_sorption_isotherm( tempw );
    } else if ( mode == VM_Incremental ) {
        return inverse_sorption_isotherm( tempw ) - inverse_sorption_isotherm( w );
    } else if ( mode == VM_Velocity ) { // VM_Previous
        return inverse_sorption_isotherm( w );
    }

    return 1.;
}


double
HeMoTKMaterial :: perm_ww(double w, double t) const
{
    // Function calculates permability water content-water content (k_ww)
    // phi ... relative humidity
    // delta_gw ... water vapor permeability
    // dphi_dw ... differentiation of relative with respect to water content
    // p_gws ... saturation water vapor pressure
    double phi = inverse_sorption_isotherm(w);
    double delta_gw = give_delta_gw(phi);
    double dphi_dw = give_dphi_dw(w);
    double p_gws = give_p_gws(t);

    return delta_gw * p_gws * dphi_dw;
}

double
HeMoTKMaterial :: perm_wt(double w, double t) const
{
    // Function calculates permability water content-temperature (k_wt)
    // delta_gw ... water vapor permeability
    // d_pgw_d_t ... differentiation of water vapor presssure with respect to temperature
    double phi = inverse_sorption_isotherm(w);
    double delta_gw = give_delta_gw(phi);
    double dpgw_dt = give_dpgw_dt(t, phi);

    return delta_gw * dpgw_dt;
}


double
HeMoTKMaterial :: give_delta_gw(double phi) const
// Function calculates the water vapor permeability delta_gw. Relative humidity (phi) is from range 0.2 - 0.98 !!!
// delta_gw ... by Z. P. Bazant and L. J. Najjar (1972), Nonlinear water diffusion in nonsaturated concrete,
//              MATERIAUX ET CONSTRUCTIONS, Vol. 5, No. 25, pp. 3 -20.
// phi ... relative humidity
// a_0, nn, phi_c, delta_wet ... constants obtained from experiments
{
    if ( phi < 0.2 || phi > 0.98 ) {
        OOFEM_ERROR("Relative humidity is out of range");
    }

    // water vapor permeability
    return delta_wet * ( a_0 + ( 1.0 - a_0 ) / ( 1.0 + pow( ( 1.0 - phi ) / ( 1.0 - phi_c ), nn ) ) );
}


double
HeMoTKMaterial :: sorption_isotherm(double phi) const
// Function calculates water content in medium from relative humidity.  Relative humidity (phi) is from range 0.2 - 0.98 !!!
// sorption isotherm by C. R. Pedersen (1990), Combined heat and moisture transfer in building constructions,
//                      PhD-thesis, Technical University of Denmark, Lingby.
// w (kg/kg) ... water content
// phi ... relative humidity
// w_h, n, a ... constants obtained from experiments
{
    if ( phi < 0.2 || phi > 0.98 ) {
        OOFEM_ERROR("Relative humidity %.3f is out of range", phi);
    }

    // water content
    return w_h * pow( ( 1.0 - log(phi) / a ), ( -1.0 / n ) );
}


double
HeMoTKMaterial :: inverse_sorption_isotherm(double w) const
// Function calculates relative humidity from water content (inverse relation form sorption isotherm).
// Relative humidity (phi) is from range 0.2 - 0.98 !!!
// sorption isotherm by C. R. Pedersen (1990), Combined heat and moisture transfer in building constructions,
//                      PhD-thesis, Technical University of Denmark, Lingby.
// w (kg/kg) ... water content
// phi ... relative humidity
// w_h, n, a ... constants obtained from experiments
{
    // relative humidity
    double phi = exp( a * ( 1.0 - pow( ( w_h / w ), ( n ) ) ) );

    if ( phi < 0.2 || phi > 0.98 ) {
        OOFEM_ERROR("Relative humidity %.3f is out of range", phi);
    }

    return phi;
}


double
HeMoTKMaterial :: give_dphi_dw(double w) const
// Function calculates differentiation of relative humidity with respect to water content (inverse relation form sorption isotherm).
// Relative humidity (phi) is from range 0.2 - 0.98 !!!
// sorption isotherm by C. R. Pedersen (1990), Combined heat and moisture transfer in building constructions,
//                      PhD-thesis, Technical University of Denmark, Lingby.
// w (kg/kg) ... water content
// phi ... relative humidity
// w_h, n, a ... constants obtained from experiments
{
    // diferentiation of realative humidity with respect to water content
    return exp( a * ( 1.0 - pow( ( w_h / w ), n ) ) ) * a * n * pow(w_h, n) * pow( w, ( -1.0 - n ) );
}

double
HeMoTKMaterial :: give_dpgw_dt(double t, double phi) const
// Function calculates differentiation of water vapor pressure with respect to temperature
// saturation water vapor pressure by C. R. Pedersen (1990), Combined heat and moisture transfer in building constructions,
//                                 PhD-thesis, Technical University of Denmark, Lingby.
// t ... temperature [K]
// phi ... relative humidity
{
    //differentiation of saturation water vapor pressure (analytical expression) with respect to temperature
    double dp_gws_dt = exp( 23.5771 - 4042.9 / ( t - 37.58 ) ) * 4042.9 / ( t - 37.58 ) / ( t - 37.58 );
    double dp_gw_dt = phi * dp_gws_dt;
    return dp_gw_dt;
}


double
HeMoTKMaterial :: give_p_gws(double t) const
// Function calculates saturation water vapor pressure
// saturation water vapor pressure by C. R. Pedersen (1990), Combined heat and moisture transfer in building constructions,
//                                 PhD-thesis, Technical University of Denmark, Lingby.
// t ... temperature [K]
// phi ... relative humidity
{
    //saturation water vapor pressure ... analytical expression
    return exp( 23.5771 - 4042.9 / ( t - 37.58 ) );
}

double
HeMoTKMaterial :: get_latent(double w, double t) const
{
    // Function calculates latent heat of evaporation
    return latent; //zatim!!!!
}

double
HeMoTKMaterial :: get_b(double w, double t) const
{
    // Function calculates coefficient b
    // sat ... degree of saturation
    double phi = inverse_sorption_isotherm(w);
    double dphi_dw = give_dphi_dw(w);
    double sat = get_sat(w, t);

    return por * rho_gws * ( phi + ( 1 - sat ) * dphi_dw );
}


double
HeMoTKMaterial :: get_sat(double w, double t) const
{
    // Function calculates degree of saturation
    // sat ... degree of saturation
    return 1.0; //zatim!!!!
}

double
HeMoTKMaterial :: get_ceff(double w, double t) const
{
    // Function calculates effective thermal capacity
    return c * rho; //zatim!!!!
}


double
HeMoTKMaterial :: get_chi(double w, double t) const
{
    // Function calculates effective thermal conductivity
    return chi_eff; //zatim!!!!
}

bool
HeMoTKMaterial :: isCharacteristicMtrxSymmetric(MatResponseMode mode) const
{
    if ( mode == Conductivity_ww || mode == Conductivity_hh || mode == Conductivity_hw || mode == Conductivity_wh ) {
        return false;
    } else {
        OOFEM_ERROR("unknown mode (%s)", __MatResponseModeToString(mode) );
    }

    // return false; // to make compiler happy
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
