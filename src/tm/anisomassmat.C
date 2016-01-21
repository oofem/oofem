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

#include "anisomassmat.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "classfactory.h"

#include <cstdlib>

namespace oofem {
REGISTER_Material(AnisotropicMassTransferMaterial);

IRResultType
AnisotropicMassTransferMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                          // Required by IR_GIVE_FIELD macro

    ///@todo Why hardcode this for 2d ? Just take the whole matrix as the input instead and not worry about it.
#if 0
    IR_GIVE_FIELD(ir, k, _IFT_AnisotropicMassTransferMaterial_c);     // Read permeability matrix c from input file
#else
    FloatArray temp;
    IR_GIVE_FIELD(ir, temp, _IFT_AnisotropicMassTransferMaterial_c);     // Read permeability matrix c from input file
    k.resize(2, 2);
    k.at(1, 1) = temp.at(1);
    k.at(1, 2) = temp.at(2);
    k.at(2, 1) = temp.at(3);
    k.at(2, 2) = temp.at(4);
#endif

    return Material :: initializeFrom(ir);
}


void
AnisotropicMassTransferMaterial :: giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep)
{
    TransportMaterialStatus *ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );

    answer.beProductOf(k, grad);
    answer.negated();

    ms->setTempField(field);
    ms->setTempGradient(grad);
    ms->setTempFlux(answer);
}


void
AnisotropicMassTransferMaterial :: giveCharacteristicMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch  ( mMode ) {
    case _1dHeat:
    case _2dHeat:
    case _3dHeat:
        answer = k;
        return;

    default:
        OOFEM_ERROR("unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}


double
AnisotropicMassTransferMaterial :: giveCharacteristicValue(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep)
{
    OOFEM_ERROR("unknown mode (%s)", __MatResponseModeToString(mode) );

    return 0.;
}
} // end namespace oofem
