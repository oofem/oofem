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

#include "isoheatmat.h"
#include "floatmatrix.h"
#include "function.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "engngm.h"

namespace oofem {
REGISTER_Material(IsotropicHeatTransferMaterial);

IsotropicHeatTransferMaterial :: IsotropicHeatTransferMaterial(int n, Domain *d) : TransportMaterial(n, d)
{
    // constructor
    maturityT0 = 0.;
}

IsotropicHeatTransferMaterial :: ~IsotropicHeatTransferMaterial() {
    // destructor
}

IRResultType
IsotropicHeatTransferMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, conductivity, _IFT_IsotropicHeatTransferMaterial_k);
    IR_GIVE_FIELD(ir, capacity, _IFT_IsotropicHeatTransferMaterial_c);
    IR_GIVE_OPTIONAL_FIELD(ir, maturityT0, _IFT_IsotropicHeatTransferMaterial_maturityT0);
    IR_GIVE_OPTIONAL_FIELD(ir, density, _IFT_IsotropicHeatTransferMaterial_d);

    return Material :: initializeFrom(ir);
}

double
IsotropicHeatTransferMaterial :: give(int aProperty, GaussPoint *gp, TimeStep *tStep)
//
// Returns the value of the property aProperty (e.g. 'k' the conductivity of the receiver).
//
{
    if ( aProperty == 'k' ) { //thermal conductivity [W/m/K]   
        return conductivity.eval( { { "te", giveTemperature(gp) }, { "t", tStep->giveIntrinsicTime() } }, this->giveDomain(), gp, giveTemperature(gp) );
    } else if ( aProperty == 'c' ) { //mass-specific heat capacity [J/kg/K]
        return capacity.eval( { { "te", giveTemperature(gp) }, { "t", tStep->giveIntrinsicTime() } }, this->giveDomain(), gp, giveTemperature(gp) );
    } else if ( aProperty == 'd' && density.isDefined() ) { //density [kg/m3]
        return density.eval( { { "te", giveTemperature(gp) }, { "t", tStep->giveIntrinsicTime() } }, this->giveDomain(), gp, giveTemperature(gp) );
    } else if ( aProperty == HeatCapaCoeff ) { //volume-specific heat capacity [J/m3/K]
        return ( this->give('c', gp, tStep) * this->give('d', gp, tStep) );
    }

    return this->Material :: give(aProperty, gp);
}


void
IsotropicHeatTransferMaterial :: giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep)
{
    TransportMaterialStatus *ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );

    ms->setTempField(field);
    ms->setTempGradient(grad);

    ///@todo Shouldn't the conductivity typically depend on the primary field and/or its gradient?
    answer.beScaled(-this->giveIsotropicConductivity(gp, tStep), grad);

    ms->setTempFlux(answer);
}


void
IsotropicHeatTransferMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                                          MatResponseMode mode,
                                                          GaussPoint *gp,
                                                          TimeStep *tStep)
{
    /*
     * returns constitutive (conductivity) matrix of receiver
     */
    MaterialMode mMode = gp->giveMaterialMode();
    double cond = this->giveIsotropicConductivity(gp, tStep);

    switch  ( mMode ) {
    case _1dHeat:
        answer.resize(1, 1);
        answer.at(1, 1) = cond;
    case _2dHeat:
        answer.resize(2, 2);
        answer.at(1, 1) = cond;
        answer.at(2, 2) = cond;
        return;

    case _3dHeat:
        answer.resize(3, 3);
        answer.at(1, 1) = cond;
        answer.at(2, 2) = cond;
        answer.at(3, 3) = cond;
        return;

    default:
        OOFEM_ERROR("unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}

double
IsotropicHeatTransferMaterial :: giveIsotropicConductivity(GaussPoint *gp, TimeStep *tStep) {
    return give('k', gp, tStep);
}

double
IsotropicHeatTransferMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                                         GaussPoint *gp,
                                                         TimeStep *tStep)
{
    if ( mode == Capacity ) {
        return ( this->give('c', gp, tStep) * this->give('d', gp, tStep) );
    } else {
        OOFEM_ERROR("unknown mode (%s)", __MatResponseModeToString(mode) );
    }

    return 0.;
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
        answer = FloatArray{ this->give('d', gp, tStep) };
        return 1;
    } else if ( type == IST_HeatCapacity ) {
        answer = FloatArray{ this->give('c', gp, tStep) };
        return 1;
    } else if ( type == IST_ThermalConductivityIsotropic ) {
        answer = FloatArray{ this->give('k', gp, tStep) };
        return 1;
    }

    return TransportMaterial :: giveIPValue(answer, gp, type, tStep);
}






MaterialStatus *
IsotropicHeatTransferMaterial :: CreateStatus(GaussPoint *gp) const
{
    return new IsotropicHeatTransferMaterialStatus(1, domain, gp);
}

IsotropicHeatTransferMaterialStatus :: IsotropicHeatTransferMaterialStatus(int n, Domain *d, GaussPoint *g) : TransportMaterialStatus(n, d, g)
{
    //constructor
}


IsotropicHeatTransferMaterialStatus :: ~IsotropicHeatTransferMaterialStatus()
{
    //destructor
}

void
IsotropicHeatTransferMaterialStatus :: updateYourself(TimeStep *tStep)
{
    TransportMaterialStatus :: updateYourself(tStep);
}

double IsotropicHeatTransferMaterial :: giveTemperature(GaussPoint *gp)
{
    IsotropicHeatTransferMaterialStatus *ms = static_cast< IsotropicHeatTransferMaterialStatus * >( this->giveStatus(gp) );
    return ms->giveTempField().at(1);
}

} // end namespace oofem
