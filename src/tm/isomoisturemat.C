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

#include "isomoisturemat.h"
#include "floatmatrix.h"
#include "gausspoint.h"

namespace oofem {
IRResultType
IsotropicMoistureTransferMaterial :: initializeFrom(InputRecord *ir)
{
    return Material :: initializeFrom(ir);
}

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
                                                              MatResponseMode mode,
                                                              GaussPoint *gp,
                                                              TimeStep *tStep)
{
    /*
     * returns constitutive (conductivity) matrix of receiver
     */

    double permeability;
    permeability = this->givePermeability(gp, tStep);

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
        OOFEM_ERROR("unknown mode (%s)", __MaterialModeToString(mMode) );
    }

    return;
}


double
IsotropicMoistureTransferMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                                             GaussPoint *gp,
                                                             TimeStep *tStep)
{
    if ( mode == Capacity ) {
        return ( this->giveMoistureCapacity(gp, tStep) );
    } else {
        OOFEM_ERROR("unknown mode (%s)", __MatResponseModeToString(mode) );
    }

    return 0.;
}


int
IsotropicMoistureTransferMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if (  type == IST_HydrationDegree ) {
        answer.resize(1);
        answer.at(1) = 0.;
        return 1;
    }

    /* else if (  type == IST_Humidity ) {
     * FloatArray state = static_cast< TransportMaterialStatus * >( giveStatus(gp) )->giveStateVector();
     * if ( state.giveSize() < 1 ) {
     *  OOFEM_ERROR("undefined moisture status!");
     * }
     *
     * answer.resize(1);
     * answer.at(1) =  state.at(1);
     * return 1;
     * }*/

    return TransportMaterial :: giveIPValue(answer, gp, type, tStep);
}
} // end namespace oofem
