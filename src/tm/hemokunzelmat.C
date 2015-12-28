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

#include "hemokunzelmat.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(HeMoKunzelMaterial);

int
HeMoKunzelMaterial :: hasMaterialModeCapability(MaterialMode mode)
{
    if ( ( mode == _2dHeMo ) || ( mode == _3dHeMo ) ) {
        return 1;
    }

    return 0;
}


IRResultType
HeMoKunzelMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

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


    return Material :: initializeFrom(ir);
}


double
HeMoKunzelMaterial :: give(int aProperty, GaussPoint *gp)
//
// Returns the value of the property aProperty (e.g. the Young's modulus
// 'E') of the receiver.
//
{
    return Material :: give(aProperty, gp);
}


void
HeMoKunzelMaterial :: giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep)
{
    TransportMaterialStatus *ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );

    ms->setTempField(field);
    ms->setTempGradient(grad);

    double h = field.at(2);
    double t = field.at(1);

    int size = grad.giveSize() / 2;
    FloatArray ans_h, ans_m;
    FloatArray grad_h(size), grad_m(size);
    for ( int i = 1; i <= size; ++i ) {
        grad_h.at(i) = grad.at(i);
    }
    for ( int i = 1; i <= 2; ++i ) {
        grad_m.at(i) = grad.at(i+size);
    }

    ans_m.add(-perm_mm(h, t), grad_m);
    ans_m.add(-perm_mh(h, t), grad_h);
    ans_h.add(-perm_hm(h, t), grad_m);
    ans_h.add(-perm_hh(h, t), grad_h);

    answer.resize(size * 2);
    answer.zero();
    answer.addSubVector(ans_h, 1);
    answer.addSubVector(ans_m, size + 1);

    ms->setTempFlux(answer);
}


void
HeMoKunzelMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
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
        OOFEM_ERROR( "giveCharacteristicMatrix : unknown mode (%s)", __MatResponseModeToString(mode) );
    }
}


double
HeMoKunzelMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                              GaussPoint *gp,
                                              TimeStep *atTime)
{
    return this->computeCapacityCoeff(mode, gp, atTime);
}


void HeMoKunzelMaterial :: computeConductivityMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
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
        OOFEM_ERROR("Unsupported MaterialMode");
    }
}


void
HeMoKunzelMaterial :: matcond1d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime)
{
    double k = 0.0, h = 0.0, t = 0.0;
    TransportMaterialStatus *status = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    FloatArray s;

    //     s = status->giveTempStateVector();
    s = status->giveTempField();
    if ( s.isEmpty() ) {
        OOFEM_ERROR("matcond1d: undefined state vector");
    }

    h = s.at(2);
    t = s.at(1);

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

    d.resize(1, 1);
    d.at(1, 1) = k;
}

void
HeMoKunzelMaterial :: matcond2d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime)
{
    double k = 0.0, h = 0.0, t = 0.0;
    TransportMaterialStatus *status = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    FloatArray s;

    //     s = status->giveTempStateVector();
    s = status->giveTempField();
    if ( s.isEmpty() ) {
        OOFEM_ERROR("matcond2d: undefined state vector");
    }

    h = s.at(2);
    t = s.at(1);

    //     if  (gp->giveElement()->giveNumber() == 4)
    //       double bzzz = 20;

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

    d.resize(2, 2);
    d.at(1, 1) = k;
    d.at(1, 2) = 0.0;
    d.at(2, 1) = 0.0;
    d.at(2, 2) = k;
}

void
HeMoKunzelMaterial :: matcond3d(FloatMatrix &d, GaussPoint *gp, MatResponseMode mode, TimeStep *atTime)
{
    double k = 0.0, h = 0.0, t = 0.0;
    TransportMaterialStatus *status = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    FloatArray s;


    //     s = status->giveTempStateVector();
    s = status->giveTempField();
    if ( s.isEmpty() ) {
        OOFEM_ERROR("matcond3d: undefined state vector");
    }

    h = s.at(2);
    t = s.at(1);

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


double HeMoKunzelMaterial :: computeCapacityCoeff(MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    //     if  (gp->giveElement()->giveNumber() == 4)
    //       double bzzz = 20;


    if ( mode == Capacity_ww ) {
        TransportMaterialStatus *status = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
        FloatArray s;
        double h;
        double dw_dh;

        //       s = status->giveTempStateVector();
        s = status->giveTempField();
        if ( s.isEmpty() ) {
            OOFEM_ERROR("computeCapacityCoeff: undefined state vector");
        }

        h = s.at(2);
        dw_dh = this->sorptionIsothermDerivative(h);

        return dw_dh;

        // CONSTANT
        //return 10.;
    } else if ( mode == Capacity_wh ) {
        return 0.0;
    } else if ( mode == Capacity_hw ) {
        return 0.0;
    } else if ( mode == Capacity_hh ) {
        TransportMaterialStatus *status = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
        FloatArray s;
        double h, w;
        double dHs_dT, dHw_dT;

        //s = status->giveTempStateVector();
        s = status->giveTempField();
        if ( s.isEmpty() ) {
            OOFEM_ERROR("computeCapacityCoeff: undefined state vector");
        }

        h = s.at(2);

        w = this->sorptionIsotherm(h);

        dHs_dT = cs * give('d', NULL);
        dHw_dT = cw * w;

        return ( dHs_dT + dHw_dT );

        // CONSTANT	return 1.7e6;
    } else {
        OOFEM_ERROR("Unknown MatResponseMode");
    }

    return 0.0; // to make compiler happy
}


double
HeMoKunzelMaterial :: sorptionIsotherm(double h)
{
    double w = 0.;

    if ( ( h < 0.0 ) || ( h > 1.00 ) ) {
        OOFEM_ERROR("HeMoKunzelMaterial :: sorptionIsotherm : Relative humidity %.3f is out of range", h);
    }

    if ( this->Isotherm == Hansen ) {
        w = iso_wh * pow( ( 1.0 - log(h) / iso_a ), ( -1.0 / iso_n ) );
    } else if ( this->Isotherm == Kunzeliso ) {
        w = iso_wh * ( iso_b - 1. ) * h / ( iso_b - h );
    } else {
        OOFEM_ERROR("Unknown Isotherm type");
    }

    return ( w );
}

double
HeMoKunzelMaterial :: sorptionIsothermDerivative(double h)
{
    double dw_dh = 0.;

    if ( ( h < 0.0 ) || ( h > 1.00 ) ) {
        OOFEM_ERROR("HeMoKunzelMaterial :: sorptionIsothermDerivative : Relative humidity %.3f is out of range", h);
    }

    if ( this->Isotherm == Hansen ) {
        dw_dh = iso_wh / ( iso_n * iso_a * h ) * pow( ( 1.0 - log(h) / iso_a ), ( -( 1.0 + iso_n ) / iso_n ) );
    } else if ( this->Isotherm == Kunzeliso ) {
        dw_dh = iso_wh * ( iso_b - 1. ) * iso_b / ( ( iso_b - h ) * ( iso_b - h ) );
    } else {
        OOFEM_ERROR("Unknown Isotherm type");
    }

    return ( dw_dh );
}


double
HeMoKunzelMaterial :: computeWaterVaporPerm(double T)
{
    /// vapor diffusion coefficient in air [kg m^-1 s^-1 Pa^-1]
    double delta;
    double deltap;

    delta = 2.0 * 1.e-7 * pow(T, 0.81) / PL;
    deltap = delta / mu;

    return ( deltap );
}

double
HeMoKunzelMaterial :: computeSatVaporPressure(double T)
{
    double p_sat;
    double T0, a;
    double T_C = T - 273.15;

    if ( T_C < 0. ) {
        a = 22.44;
        T0 = 272.44;
    } else {
        a = 17.08;
        T0 = 234.18;
    }

    p_sat = 611. * exp( a * T_C / ( T0 + T_C ) );

    return p_sat;
}

double
HeMoKunzelMaterial :: computeSatVaporPressureDerivative(double T)
{
    double dpsat_dT;
    double T0, a;
    double T_C = T - 273.15;

    if ( T_C < 0. ) {
        a = 22.44;
        T0 = 272.44;
    } else {
        a = 17.08;
        T0 = 234.18;
    }

    dpsat_dT = 611. *a *T0 *exp( a *T_C / ( T0 + T_C ) ) / ( ( T0 + T_C ) * ( T0 + T_C ) );

    return dpsat_dT;
}


double
HeMoKunzelMaterial :: computeDw(double h)
{
    double Dw = 0.;

    if ( this->Permeability == Multilin_h ) {
        double tol = 1.e-10;
        for ( int i = 1; i <= perm_h.giveSize(); i++ ) {
            if ( ( h - perm_h.at(i) ) < tol ) {
                Dw = perm_Dwh.at(i - 1) + ( perm_Dwh.at(i) - perm_Dwh.at(i - 1) ) * ( h - perm_h.at(i - 1) ) / ( perm_h.at(i) - perm_h.at(i - 1) );
                break;
            }
        }
    } else if ( this->Permeability == Multilin_wV ) {
        double wV = this->sorptionIsotherm(h) / rhoH2O;
        double tol = 1.e-10;
        for ( int i = 1; i <= perm_wV.giveSize(); i++ ) {
            if ( ( wV - perm_wV.at(i) ) < tol ) {
                Dw = perm_DwwV.at(i - 1) + ( perm_DwwV.at(i) - perm_DwwV.at(i - 1) ) * ( wV - perm_wV.at(i - 1) ) / ( perm_wV.at(i) - perm_wV.at(i - 1) );
                break;
            }
        }
    } else if ( this->Permeability == Kunzelperm ) {
        double w;
        w = this->sorptionIsotherm(h);
        Dw = 3.8 * ( A / iso_wh ) * ( A / iso_wh ) * pow(1000., w / iso_wh - 1.);
    } else {
        OOFEM_ERROR("initializeFrom: permeabilityType must be equal to 0, 1 or 2");
    }

    return Dw;
}


double
HeMoKunzelMaterial :: perm_mm(double h, double T)
{
    // Function calculates permability relative humidity - relative humidity (k_mm)
    // h     ... relative humidity
    // Dw    ... capillary transport coefficient [m^2/s]
    // dw_dh ... derivative of moisture storage function [kg/m^3]
    // p_sat ... saturation water vapor pressure

    double k_mm;
    double dw_dh, Dw;
    double deltap, p_sat;

    dw_dh = this->sorptionIsothermDerivative(h);

    Dw = this->computeDw(h);

    deltap = this->computeWaterVaporPerm(T);
    p_sat = this->computeSatVaporPressure(T);

    k_mm = Dw * dw_dh + deltap * p_sat;

    return ( k_mm );

    //return 5.e-8;
}

double
HeMoKunzelMaterial :: perm_mh(double h, double T)
{
    // Function calculates permeability relative humidity-temperature (k_mh)
    // deltap   ... water vapor permeability
    // dpsat_dt ... differentiation of saturation water vapor presssure with respect to temperature

    double k_mh;
    double delta_p, dpsat_dT;

    delta_p = computeWaterVaporPerm(T);
    dpsat_dT = computeSatVaporPressureDerivative(T);

    k_mh = delta_p * h * dpsat_dT;

    //     return 2.e-7;

    return k_mh;
}


double
HeMoKunzelMaterial :: perm_hm(double h, double T)
{
    // Function calculates permability temperature-relative humidity (k_hm)

    double k_hm;
    double deltap, p_sat;

    deltap = this->computeWaterVaporPerm(T);
    p_sat = this->computeSatVaporPressure(T);

    k_hm = hv * deltap * p_sat;

    return ( k_hm );

    //   return 0.1;
    //    return 0.;
}

double
HeMoKunzelMaterial :: perm_hh(double h, double T)
{
    // Function calculates permability water temperature-temperature (k_hh)

    double k_hh;
    double lambda, deltap, dpsat_dT, w;

    w = this->sorptionIsotherm(h);
    lambda = lambda0 * ( 1. + b * w / give('d', NULL) );
    deltap = this->computeWaterVaporPerm(T);
    dpsat_dT = computeSatVaporPressureDerivative(T);

    k_hh = lambda + hv * deltap * h * dpsat_dT;

    return ( k_hh );

    // return 2.;
}

bool
HeMoKunzelMaterial :: isCharacteristicMtrxSymmetric(MatResponseMode mode)
{
    if ( ( mode == Conductivity_ww ) || ( mode == Conductivity_hh ) || ( mode == Conductivity_hw ) || ( mode == Conductivity_wh ) ) {
        return false;
    } else {
        OOFEM_ERROR( "isCharacteristicMtrxSymmetric : unknown mode (%s)", __MatResponseModeToString(mode) );
    }

    return false; // to make compiler happy
}

int
HeMoKunzelMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime)
// IST_Humidity overriden to use inverse_sorption_isotherm
{
    if ( type == IST_Humidity ) {
        answer.resize(1);
        answer.at(1) = giveHumidity(gp, VM_Velocity); // VM_Previous = equilibrated value of humidity
        return 1;
    } else {
        return TransportMaterial :: giveIPValue(answer, gp, type, atTime);
    }
}

double
HeMoKunzelMaterial :: giveHumidity(GaussPoint *gp, ValueModeType mode)
{
    TransportMaterialStatus *ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    const FloatArray &tempState = ms->giveTempField();
    if ( tempState.giveSize() < 2 ) {
        OOFEM_ERROR("Undefined moisture status");
    }

    const FloatArray &state = ms->giveField();

    if ( mode == VM_Total ) {
        return tempState.at(2);
    } else if ( mode == VM_Incremental ) {
        return tempState.at(2) - state.at(2);
    } else if ( mode == VM_Velocity ) { // VM_Previous
        return state.at(2);
    } else {
        OOFEM_ERROR("Undefined moisture mode");
    }

    return 1.;
}
} // end namespace oofem
