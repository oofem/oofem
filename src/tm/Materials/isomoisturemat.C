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

#include "tm/Materials/isomoisturemat.h"
#include "floatmatrix.h"
#include "gausspoint.h"

namespace oofem {
void
IsotropicMoistureTransferMaterial :: initializeFrom(InputRecord &ir)
{
    Material :: initializeFrom(ir);
}


FloatArrayF<3>
IsotropicMoistureTransferMaterial :: computeFlux3D(const FloatArrayF<3> &grad, double field, GaussPoint *gp, TimeStep *tStep) const
{
    auto ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    ms->setTempGradient(grad);
    ms->setTempField(field);

    auto answer = -this->givePermeability(gp, tStep) * grad;
    ms->setTempFlux(answer);
    return answer;
}


FloatMatrixF<3,3>
IsotropicMoistureTransferMaterial :: computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    double permeability = this->givePermeability(gp, tStep);
    return permeability * eye<3>();
}


double
IsotropicMoistureTransferMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                                             GaussPoint *gp,
                                                             TimeStep *tStep) const
{
    if ( mode == Capacity ) {
        return this->giveMoistureCapacity(gp, tStep);
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
    } else if (  type == IST_MoistureContent ) {
        FloatArray humidity;
        answer.resize(1);
        this->giveIPValue(humidity, gp, IST_Humidity, tStep);
        answer.at(1) = giveMoistureContent(humidity.at(1));
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
