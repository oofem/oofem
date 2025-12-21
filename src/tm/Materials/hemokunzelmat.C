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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "tm/Materials/hemokunzelmat.h"
#include "floatmatrixf.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(HeMoKunzelMaterial);

bool
HeMoKunzelMaterial :: hasMaterialModeCapability(MaterialMode mode) const
{
    return mode == _2dHeMo || mode == _3dHeMo;
}


void
HeMoKunzelMaterial :: initializeFrom(InputRecord &ir)
{
    Material :: initializeFrom(ir);

    int type = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, type, _IFT_HeMoKunzelMaterial_iso_type);
    this->Isotherm = ( isothermType ) type;

    if ( this->Isotherm == Hansen ) {
        IR_GIVE_FIELD(ir, iso_wh, _IFT_HeMoKunzelMaterial_iso_wh);
        IR_GIVE_FIELD(ir, iso_n, _IFT_HeMoKunzelMaterial_iso_n);
        IR_GIVE_FIELD(ir, iso_a, _IFT_HeMoKunzelMaterial_iso_a);
    } else if ( this->Isotherm == Kunzeliso ) {
        IR_GIVE_FIELD(ir, iso_wh, _IFT_HeMoKunzelMaterial_iso_wh);
        IR_GIVE_FIELD(ir, iso_b, _IFT_HeMoKunzelMaterial_iso_b);
    } else {
        OOFEM_ERROR("Unknown Isotherm type");
    }


    IR_GIVE_FIELD(ir, type, _IFT_HeMoKunzelMaterial_permeability_type);

    if ( type >= 3 ) {
        OOFEM_ERROR("initializeFrom: isothermType must be equal to 0, 1, 2");
    }

    this->Permeability = ( permeabilityType ) type;

    if ( this->Permeability == Multilin_h ) {
        IR_GIVE_FIELD(ir, perm_h, _IFT_HeMoKunzelMaterial_perm_h);
        IR_GIVE_FIELD(ir, perm_Dwh, _IFT_HeMoKunzelMaterial_perm_dwh);

        if ( !( perm_h.giveSize() == perm_Dwh.giveSize() ) ) {
            OOFEM_ERROR("initializeFrom: the size of 'perm_h' and 'perm_dw(h)' must be the same");
        }

        for ( int i = 1; i < perm_h.giveSize(); i++ ) {
            if ( ( perm_h.at(i) < 0. ) || ( perm_h.at(i) > 1. ) ) {
                OOFEM_ERROR("initializeFrom: perm_h must be in the range <0; 1>");
            }
        }
    } else if ( this->Permeability == Multilin_wV ) {
        IR_GIVE_FIELD(ir, perm_wV, _IFT_HeMoKunzelMaterial_perm_wv);
        IR_GIVE_FIELD(ir, perm_DwwV, _IFT_HeMoKunzelMaterial_perm_dwwv);

        if ( !( perm_wV.giveSize() == perm_DwwV.giveSize() ) ) {
            OOFEM_ERROR("initializeFrom: the size of 'perm_h' and 'perm_dw(h)' must be the same");
        }

        for ( int i = 1; i < perm_wV.giveSize(); i++ ) {
            if ( ( perm_wV.at(i) < 0. ) || ( perm_wV.at(i) > 1. ) ) {
                OOFEM_ERROR("initializeFrom: perm_wv must be in the range <0; 1>");
            }
        }
    } else if ( this->Permeability == Kunzelperm ) {
        IR_GIVE_FIELD(ir, A, _IFT_HeMoKunzelMaterial_a);
    } else {
        OOFEM_ERROR("initializeFrom: permeabilityType must be equal to 0, 1 or 2");
    }

    IR_GIVE_FIELD(ir, mu, _IFT_HeMoKunzelMaterial_mu);

    PL = 101325.;
    IR_GIVE_OPTIONAL_FIELD(ir, PL, _IFT_HeMoKunzelMaterial_pl);

    rhoH2O = 1000.;
    IR_GIVE_OPTIONAL_FIELD(ir, rhoH2O, _IFT_HeMoKunzelMaterial_rhoh2o);

    IR_GIVE_FIELD(ir, lambda0, _IFT_HeMoKunzelMaterial_lambda0);
    IR_GIVE_FIELD(ir, b, _IFT_HeMoKunzelMaterial_b);
    IR_GIVE_FIELD(ir, cs, _IFT_HeMoKunzelMaterial_cs);

    //     cw = 4.1813e3;
    cw = 4183.;
    IR_GIVE_OPTIONAL_FIELD(ir, cw, _IFT_HeMoKunzelMaterial_cw);

    hv = 2.5e6;
    IR_GIVE_OPTIONAL_FIELD(ir, hv, _IFT_HeMoKunzelMaterial_hv);
}


double
HeMoKunzelMaterial :: give(int aProperty, GaussPoint *gp) const
{
    return Material :: give(aProperty, gp);
}


std::pair<FloatArrayF<3>, FloatArrayF<3>>
HeMoKunzelMaterial :: computeHeMoFlux3D(const FloatArrayF<3> &grad_t, const FloatArrayF<3> &grad_w, double t, double h, GaussPoint *gp, TimeStep *tStep) const
{
    auto ms = static_cast< HeMoTransportMaterialStatus * >( this->giveStatus(gp) );

    ms->setTempTemperature(t);
    ms->setTempHumidity(h);

    ms->setTempTemperatureGradient(grad_t);
    ms->setTempHumidityGradient(grad_w);

    auto ans_t = -perm_hm(h, t) * grad_w - perm_hh(h, t) * grad_t;
    auto ans_w = -perm_mm(h, t) * grad_w - perm_mh(h, t) * grad_t;
    
    ms->setTempHeatFlux(ans_t);
    ms->setTempHumidityFlux(ans_w);

    return {ans_t, ans_w};
}


FloatMatrixF<3,3>
HeMoKunzelMaterial :: computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto status = static_cast< HeMoTransportMaterialStatus * >( this->giveStatus(gp) );

    double t = status->giveTempTemperature();
    double h = status->giveTempHumidity();

    double k = 0.0;
    if ( mode == Conductivity_ww ) {
        k = perm_mm(h, t);
    } else if ( mode == Conductivity_wh ) {
        k = perm_mh(h, t);
    } else if ( mode == Conductivity_hw ) {
        k = perm_hm(h, t);
    } else if ( mode == Conductivity_hh ) {
        k = perm_hh(h, t);
    } else {
        OOFEM_ERROR("Unknown MatResponseMode");
    }

    return k * eye<3>();
}


double
HeMoKunzelMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                              GaussPoint *gp,
                                              TimeStep *tStep) const
{
    return this->computeCapacityCoeff(mode, gp, tStep);
}


double HeMoKunzelMaterial :: computeCapacityCoeff(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    //     if  (gp->giveElement()->giveNumber() == 4)
    //       double bzzz = 20;

    if ( mode == Capacity_ww ) {
        auto status = static_cast< HeMoTransportMaterialStatus * >( this->giveStatus(gp) );

        double h = status->giveTempHumidity();
        return this->giveMoistureContentDerivative(h);
    } else if ( mode == Capacity_wh ) {
        return 0.0;
    } else if ( mode == Capacity_hw ) {
        return 0.0;
    } else if ( mode == Capacity_hh ) {
        auto status = static_cast< HeMoTransportMaterialStatus * >( this->giveStatus(gp) );
 
        double h = status->giveTempHumidity();
        double w = this->giveMoistureContent(h);

        double dHs_dT = cs * give('d', nullptr);
        double dHw_dT = cw * w;

        return dHs_dT + dHw_dT;

        // CONSTANT	return 1.7e6;
    } else {
        OOFEM_ERROR("Unknown MatResponseMode");
    }

    // return 0.0; // to make compiler happy
}


double
HeMoKunzelMaterial :: giveMoistureContent(double h) const
{
    if ( h < 0.0 || h > 1.00 ) {
        OOFEM_ERROR("HeMoKunzelMaterial :: giveMoistureContent : Relative humidity %.3f is out of range", h);
    }

    if ( this->Isotherm == Hansen ) {
        return iso_wh * pow( ( 1.0 - log(h) / iso_a ), ( -1.0 / iso_n ) );
    } else if ( this->Isotherm == Kunzeliso ) {
        return iso_wh * ( iso_b - 1. ) * h / ( iso_b - h );
    } else {
        OOFEM_ERROR("Unknown Isotherm type");
    }
}

double
HeMoKunzelMaterial :: giveMoistureContentDerivative(double h) const
{
    if ( h < 0.0 || h > 1.00 ) {
        OOFEM_ERROR("HeMoKunzelMaterial :: giveMoistureContentDerivative : Relative humidity %.3f is out of range", h);
    }

    if ( this->Isotherm == Hansen ) {
        return iso_wh / ( iso_n * iso_a * h ) * pow( ( 1.0 - log(h) / iso_a ), ( -( 1.0 + iso_n ) / iso_n ) );
    } else if ( this->Isotherm == Kunzeliso ) {
        return iso_wh * ( iso_b - 1. ) * iso_b / ( ( iso_b - h ) * ( iso_b - h ) );
    } else {
        OOFEM_ERROR("Unknown Isotherm type");
    }
}


double
HeMoKunzelMaterial :: computeWaterVaporPerm(double T) const
{
    // vapor diffusion coefficient in air [kg m^-1 s^-1 Pa^-1]
    double delta = 2.0 * 1.e-7 * pow(T, 0.81) / PL;
    return delta / mu;
}

double
HeMoKunzelMaterial :: computeSatVaporPressure(double T) const
{
    double T_C = T - 273.15;

    double T0, a;
    if ( T_C < 0. ) {
        a = 22.44;
        T0 = 272.44;
    } else {
        a = 17.08;
        T0 = 234.18;
    }

    return 611. * exp( a * T_C / ( T0 + T_C ) );
}

double
HeMoKunzelMaterial :: computeSatVaporPressureDerivative(double T) const
{
    double T_C = T - 273.15;

    double T0, a;
    if ( T_C < 0. ) {
        a = 22.44;
        T0 = 272.44;
    } else {
        a = 17.08;
        T0 = 234.18;
    }

    return 611. *a *T0 *exp( a *T_C / ( T0 + T_C ) ) / ( ( T0 + T_C ) * ( T0 + T_C ) );
}


double
HeMoKunzelMaterial :: computeDw(double h) const
{
    if ( this->Permeability == Multilin_h ) {
        double tol = 1.e-10;
        for ( int i = 1; i <= perm_h.giveSize(); i++ ) {
            if ( ( h - perm_h.at(i) ) < tol ) {
                return perm_Dwh.at(i - 1) + ( perm_Dwh.at(i) - perm_Dwh.at(i - 1) ) * ( h - perm_h.at(i - 1) ) / ( perm_h.at(i) - perm_h.at(i - 1) );
            }
        }
    } else if ( this->Permeability == Multilin_wV ) {
        double wV = this->giveMoistureContent(h) / rhoH2O;
        double tol = 1.e-10;
        for ( int i = 1; i <= perm_wV.giveSize(); i++ ) {
            if ( ( wV - perm_wV.at(i) ) < tol ) {
                return perm_DwwV.at(i - 1) + ( perm_DwwV.at(i) - perm_DwwV.at(i - 1) ) * ( wV - perm_wV.at(i - 1) ) / ( perm_wV.at(i) - perm_wV.at(i - 1) );
            }
        }
    } else if ( this->Permeability == Kunzelperm ) {
        double w;
        w = this->giveMoistureContent(h);
        return 3.8 * ( A / iso_wh ) * ( A / iso_wh ) * pow(1000., w / iso_wh - 1.);
    } else {
        OOFEM_ERROR("initializeFrom: permeabilityType must be equal to 0, 1 or 2");
    }

    return 0.;
}


double
HeMoKunzelMaterial :: perm_mm(double h, double T) const
{
    // Function calculates permability relative humidity - relative humidity (k_mm)
    // h     ... relative humidity
    // Dw    ... capillary transport coefficient [m^2/s]
    // dw_dh ... derivative of moisture storage function [kg/m^3]
    // p_sat ... saturation water vapor pressure

    double dw_dh = this->giveMoistureContentDerivative(h);

    double Dw = this->computeDw(h);

    double deltap = this->computeWaterVaporPerm(T);
    double p_sat = this->computeSatVaporPressure(T);

    return Dw * dw_dh + deltap * p_sat;
}

double
HeMoKunzelMaterial :: perm_mh(double h, double T) const
{
    // Function calculates permeability relative humidity-temperature (k_mh)
    // deltap   ... water vapor permeability
    // dpsat_dt ... differentiation of saturation water vapor presssure with respect to temperature
    double delta_p = computeWaterVaporPerm(T);
    double dpsat_dT = computeSatVaporPressureDerivative(T);

    return delta_p * h * dpsat_dT;
}


double
HeMoKunzelMaterial :: perm_hm(double h, double T) const
{
    // Function calculates permability temperature-relative humidity (k_hm)
    double deltap = this->computeWaterVaporPerm(T);
    double p_sat = this->computeSatVaporPressure(T);
    return hv * deltap * p_sat;
}

double
HeMoKunzelMaterial :: perm_hh(double h, double T) const
{
    // Function calculates permability water temperature-temperature (k_hh)

    double w = this->giveMoistureContent(h);
    double lambda = lambda0 * ( 1. + b * w / give('d', nullptr) );
    double deltap = this->computeWaterVaporPerm(T);
    double dpsat_dT = computeSatVaporPressureDerivative(T);

    return lambda + hv * deltap * h * dpsat_dT;
}

bool
HeMoKunzelMaterial :: isCharacteristicMtrxSymmetric(MatResponseMode mode) const
{
    if ( mode == Conductivity_ww || mode == Conductivity_hh || mode == Conductivity_hw || mode == Conductivity_wh ) {
        return false;
    } else {
        OOFEM_ERROR( "isCharacteristicMtrxSymmetric : unknown mode (%s)", __MatResponseModeToString(mode) );
    }

    // return false; // to make compiler happy
}

int
HeMoKunzelMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
// IST_Humidity overriden to use inverse_sorption_isotherm
{
    double humidity = static_cast< HeMoTransportMaterialStatus * >( this->giveStatus(gp) )->giveHumidity();
    if ( type == IST_Humidity ) {
        answer.resize(1);
        answer.at(1) = humidity;
        return 1;
    } else if ( type == IST_MoistureContent ) {
        answer.resize(1);
        answer.at(1) = giveMoistureContent(humidity);
        return 1;
    } else {
        return TransportMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}

} // end namespace oofem
