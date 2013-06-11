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

#include "anisomassmat.h"
#include "domain.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "classfactory.h"

#include <cstdlib>

namespace oofem {

REGISTER_Material( AnisotropicMassTransferMaterial );

IRResultType
AnisotropicMassTransferMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";       // Required by IR_GIVE_FIELD macro
    IRResultType result;                          // Required by IR_GIVE_FIELD macro

    this->Material :: initializeFrom(ir);

    FloatArray temp;

    ///@todo Why hardcode this for 2d ? Just take the whole matrix as the input instead and not worry about it.
    IR_GIVE_FIELD(ir, temp, _IFT_AnisotropicMassTransferMaterial_c);     // Read permeability matrix c from input file
    k.resize(2, 2);
    k.at(1, 1) = temp.at(1);
    k.at(1, 2) = temp.at(2);
    k.at(2, 1) = temp.at(3);
    k.at(2, 2) = temp.at(4);

    return IRRT_OK;
}


void
AnisotropicMassTransferMaterial :: giveFluxVector(FloatArray& answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep)
{
    TransportMaterialStatus *ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );

    answer.beProductOf(k, grad);
    answer.negated();

    ms->setTempField(field);
    ms->setTempGradient(grad);
    ms->setTempFlux(answer);
}


void
AnisotropicMassTransferMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                                            MatResponseForm form,
                                                            MatResponseMode mode,
                                                            GaussPoint *gp,
                                                            TimeStep *atTime)
{
    MaterialMode mMode = gp->giveMaterialMode();
    switch  ( mMode ) {
    case _1dHeat:
    case _2dHeat:
    case _3dHeat:
        answer = k;
        return;

    default:
        _error2( "giveCharacteristicMatrix : unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}


double
AnisotropicMassTransferMaterial :: giveCharacteristicValue(MatResponseMode mode, GaussPoint *gp, TimeStep *atTime)
{
    _error2( "giveCharacteristicValue : unknown mode (%s)", __MatResponseModeToString(mode) );

    return 0.;
}


int
AnisotropicMassTransferMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    TransportMaterialStatus *ms = static_cast< TransportMaterialStatus * >( this->giveStatus(aGaussPoint) );
    FloatMatrix temp;

    switch ( type ) {
    case IST_Velocity:
        answer = ms->giveFlux();
        break;
    case IST_PressureGradient:
        answer = ms->giveGradient();
        break;
    default:
      return TransportMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
    return 1;
}

} // end namespace oofem
