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

#include "isomoisturemat.h"
#include "floatmatrix.h"
#include "gausspoint.h"

namespace oofem {
IRResultType
IsotropicMoistureTransferMaterial :: initializeFrom(InputRecord *ir)
{
    this->Material :: initializeFrom(ir);

    return IRRT_OK;
}

/*
 * double
 * IsotropicMoistureTransferMaterial :: give(int aProperty, GaussPoint *gp)
 * //
 * // Returns the value of the property aProperty.
 * //
 * {
 * if ( aProperty == 'c' ) { // moisture permeability [kg/(m s)]
 *    return permeability;
 * } else if ( aProperty == 'k' ) { // moisture capacity [kg / m^3]
 *    return moistureCapacity;
 * }
 *
 * return this->Material :: give(aProperty, gp);
 * }
 */

void
IsotropicMoistureTransferMaterial :: giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep)
{
    TransportMaterialStatus *ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );

    ///@todo Shouldn't the permeability typically depend on the primary field and/or its gradient?
    answer.beScaled(-this->givePermeability(gp, tStep), grad);

    ms->setTempField(field);
    ms->setTempGradient(grad);
    ms->setTempFlux(answer);
}

void
IsotropicMoistureTransferMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                                              MatResponseForm form,
                                                              MatResponseMode mode,
                                                              GaussPoint *gp,
                                                              TimeStep *atTime)
{
    /*
     * returns constitutive (conductivity) matrix of receiver
     */

    double permeability;
    permeability = this->givePermeability(gp, atTime);

    MaterialMode mMode = gp->giveMaterialMode();
    switch  ( mMode ) {
    case _1dHeat:
        answer.resize(1, 1);
        answer.at(1, 1) = permeability;
    case _2dHeat:
        answer.resize(2, 2);
        answer.at(1, 1) = permeability;
        answer.at(2, 2) = permeability;
        return;

    case _3dHeat:
        answer.resize(3, 3);
        answer.at(1, 1) = permeability;
        answer.at(2, 2) = permeability;
        answer.at(3, 3) = permeability;
        return;

    default:
        _error2( "giveCharacteristicMatrix : unknown mode (%s)", __MaterialModeToString(mMode) );
    }

    return;
}


double
IsotropicMoistureTransferMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                                             GaussPoint *gp,
                                                             TimeStep *atTime)
{
    if ( mode == Capacity ) {
        return ( this->giveMoistureCapacity(gp, atTime) );
    } else {
        _error2( "giveCharacteristicValue : unknown mode (%s)", __MatResponseModeToString(mode) );
    }

    return 0.;
}

int
IsotropicMoistureTransferMaterial :: giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint)
{
    if ( type == IST_HydrationDegree ) {
        return 1;
    } else if ( type == IST_Humidity ) {
        return 1;
    } else {
        return TransportMaterial :: giveIPValueSize(type, aGaussPoint);
    }
}


int
IsotropicMoistureTransferMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    if (  type == IST_HydrationDegree ) {
        answer.resize(1);
        answer.at(1) = 0.;
        return 1;
    }

    /* else if (  type == IST_Humidity ) {
     * FloatArray state = ( ( TransportMaterialStatus * ) giveStatus(aGaussPoint) )->giveStateVector();
     * if ( state.giveSize() < 1 ) {
     *  _error("computeWaterChange: undefined moisture status!");
     * }
     *
     * answer.resize(1);
     * answer.at(1) =  state.at(1);
     * return 1;
     * }*/

    return TransportMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
}
} // end namespace oofem
