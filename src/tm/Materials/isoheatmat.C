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

#include "tm/Materials/isoheatmat.h"
#include "floatmatrix.h"
#include "function.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "engngm.h"

namespace oofem {
REGISTER_Material(IsotropicHeatTransferMaterial);

IsotropicHeatTransferMaterial :: IsotropicHeatTransferMaterial(int n, Domain *d) : TransportMaterial(n, d) { }


void
IsotropicHeatTransferMaterial :: initializeFrom(InputRecord &ir)
{
    Material :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, conductivity, _IFT_IsotropicHeatTransferMaterial_k);
    IR_GIVE_FIELD(ir, capacity, _IFT_IsotropicHeatTransferMaterial_c);
    IR_GIVE_OPTIONAL_FIELD(ir, maturityT0, _IFT_IsotropicHeatTransferMaterial_maturityT0);
    IR_GIVE_OPTIONAL_FIELD(ir, density, _IFT_IsotropicHeatTransferMaterial_d);
}

double
IsotropicHeatTransferMaterial :: giveProperty(int aProperty, GaussPoint *gp, TimeStep *tStep) const
{
    if ( aProperty == 'k' ) { //thermal conductivity [W/m/K]   
        return conductivity.eval( { { "te", giveTemperature(gp) }, { "t", tStep->giveIntrinsicTime() } }, this->giveDomain(), gp, giveTemperature(gp) );
    } else if ( aProperty == 'c' ) { //mass-specific heat capacity [J/kg/K]
        return capacity.eval( { { "te", giveTemperature(gp) }, { "t", tStep->giveIntrinsicTime() } }, this->giveDomain(), gp, giveTemperature(gp) );
    } else if ( aProperty == 'd' && density.isDefined() ) { //density [kg/m3]
        return density.eval( { { "te", giveTemperature(gp) }, { "t", tStep->giveIntrinsicTime() } }, this->giveDomain(), gp, giveTemperature(gp) );
    } else if ( aProperty == HeatCapaCoeff ) { //volume-specific heat capacity [J/m3/K]
        return ( this->giveProperty('c', gp, tStep) * this->giveProperty('d', gp, tStep) );
    }

    return this->Material :: give(aProperty, gp);
}


FloatArrayF<3>
IsotropicHeatTransferMaterial :: computeFlux3D(const FloatArrayF<3> &grad, double field, GaussPoint *gp, TimeStep *tStep) const
{
    auto ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );

    ms->setTempField(field);
    ms->setTempGradient(grad);

    ///@todo Shouldn't the conductivity typically depend on the primary field and/or its gradient?
    auto answer = -this->giveIsotropicConductivity(gp, tStep) * grad;
    ms->setTempFlux(answer);
    return answer;
}


FloatMatrixF<3,3>
IsotropicHeatTransferMaterial :: computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    double cond = this->giveIsotropicConductivity(gp, tStep);
    return cond * eye<3>();
}


double
IsotropicHeatTransferMaterial :: giveIsotropicConductivity(GaussPoint *gp, TimeStep *tStep) const
{
    return giveProperty('k', gp, tStep);
}

double
IsotropicHeatTransferMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                                         GaussPoint *gp,
                                                         TimeStep *tStep) const
{
    if ( mode == Capacity ) {
        return ( this->giveProperty('c', gp, tStep) * this->giveProperty('d', gp, tStep) );
    } else {
        OOFEM_ERROR("unknown mode (%s)", __MatResponseModeToString(mode) );
    }

    // return 0.;
}


int
IsotropicHeatTransferMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if (  type == IST_HydrationDegree ) {
        answer.resize(1);
        answer.at(1) = 0.;
        return 1;
    }

    if ( type == IST_Temperature ) {
        answer = FloatArray{ this->giveTemperature(gp) };
        return 1;
    } else if ( type == IST_Density ) {
        answer = FloatArray{ this->giveProperty('d', gp, tStep) };
        return 1;
    } else if ( type == IST_HeatCapacity ) {
        answer = FloatArray{ this->giveProperty('c', gp, tStep) };
        return 1;
    } else if ( type == IST_ThermalConductivityIsotropic ) {
        answer = FloatArray{ this->giveProperty('k', gp, tStep) };
        return 1;
    } else if ( type == IST_EnergyMassCapacity ) {
        answer = FloatArray{ this->giveProperty('c', gp, tStep) * this->giveProperty('d', gp, tStep) * this->giveTemperature(gp) };
        return 1;
    }

    return TransportMaterial :: giveIPValue(answer, gp, type, tStep);
}

double IsotropicHeatTransferMaterial :: giveTemperature(GaussPoint *gp) const
{
    auto ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    return ms->giveTempField();
}

} // end namespace oofem
