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

#include "tm/Materials/hemobaznajmat.h"
#include "floatmatrixf.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {

REGISTER_Material( HeMoBazNajMaterial );

bool
HeMoBazNajMaterial :: hasMaterialModeCapability(MaterialMode mode) const
{
    return mode == _2dHeMo || mode == _3dHeMo;
}


void
HeMoBazNajMaterial :: initializeFrom(InputRecord &ir)
{
    Material :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, C1, _IFT_HeMoBazNajMaterial_c1);
    IR_GIVE_FIELD(ir, n, _IFT_HeMoBazNajMaterial_n);
    IR_GIVE_FIELD(ir, alpha0, _IFT_HeMoBazNajMaterial_alpha0);
    IR_GIVE_FIELD(ir, hC, _IFT_HeMoBazNajMaterial_hc);

    this->moistureCapacity = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, moistureCapacity, _IFT_HeMoBazNajMaterial_capa);

    IR_GIVE_FIELD(ir, heatConductivity, _IFT_HeMoBazNajMaterial_k);
    IR_GIVE_FIELD(ir, heatCapacity, _IFT_HeMoBazNajMaterial_c);
}


double
HeMoBazNajMaterial :: give(int aProperty, GaussPoint *gp) const
{
    return this->Material :: give(aProperty, gp);
}


std::pair<FloatArrayF<3>, FloatArrayF<3>>
HeMoBazNajMaterial :: computeHeMoFlux3D(const FloatArrayF<3> &grad_t, const FloatArrayF<3> &grad_w, double t, double h, GaussPoint *gp, TimeStep *tStep) const
{
    auto ms = static_cast< HeMoTransportMaterialStatus * >( this->giveStatus(gp) );

    auto ans_w = perm_mm(h, t) * grad_w + perm_mh(h, t) * grad_t;
    auto ans_t = perm_hm(h, t) * grad_w + perm_hh(h, t) * grad_t;

    ms->setTempTemperature(t);
    ms->setTempTemperatureGradient(grad_t);
    ms->setTempHeatFlux(ans_t);
    ms->setTempHumidity(h);
    ms->setTempHumidityGradient(grad_w);
    ms->setTempHumidityFlux(ans_w);

    return {ans_t, ans_w};
}


FloatMatrixF<3,3>
HeMoBazNajMaterial :: computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
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
HeMoBazNajMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *atTime) const
{
    return this->computeCapacityCoeff(mode, gp, atTime);
}


double HeMoBazNajMaterial :: computeCapacityCoeff(MatResponseMode mode, GaussPoint *gp, TimeStep *atTime) const
{
    if ( mode == Capacity_ww ) {
        return this->moistureCapacity;

    } else if ( mode == Capacity_wh ) {
        return 0.0;

    } else if ( mode == Capacity_hw ) {
        return 0.0;

    } else if ( mode == Capacity_hh ) {
        return this->heatCapacity;

    } else {
        OOFEM_ERROR("Unknown MatResponseMode");
    }
}


double
HeMoBazNajMaterial :: perm_mm(double h, double T) const
{
    return C1 * ( alpha0 + ( 1. - alpha0 ) / ( 1. + pow( ( 1. - h ) / ( 1. - hC ), n ) ) );
}

double
HeMoBazNajMaterial :: perm_mh(double h, double T) const
{
    return 0.;
}

double
HeMoBazNajMaterial :: perm_hm(double h, double T) const
{
    return 0.;
}

double
HeMoBazNajMaterial :: perm_hh(double h, double T) const
{
    return this->heatConductivity;
}

bool
HeMoBazNajMaterial :: isCharacteristicMtrxSymmetric(MatResponseMode mode) const
{
    if ( mode == Conductivity_ww || mode == Conductivity_hh || mode == Conductivity_hw || mode == Conductivity_wh ) {
        return true;
    } else {
        OOFEM_ERROR( "isCharacteristicMtrxSymmetric : unknown mode (%s)", __MatResponseModeToString(mode) );
    }

    // return false; // to make compiler happy
}

int
HeMoBazNajMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime)
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
} // end namespace oofem
