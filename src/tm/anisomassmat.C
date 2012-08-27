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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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
#include "flotmtrx.h"
#include "gausspnt.h"
#ifndef __MAKEDEPEND
 #include <stdlib.h>
#endif

namespace oofem {
AnisotropicMassTransferMaterialStatus :: AnisotropicMassTransferMaterialStatus(int n, Domain *d, GaussPoint *g) : TransportMaterialStatus(n, d, g)
{
    pressureGradient.resize(2);
    seepageVelocity.resize(2);
}

void
AnisotropicMassTransferMaterialStatus :: setPressureGradient(const FloatArray &gradP)
{
    pressureGradient = gradP;
}

void
AnisotropicMassTransferMaterialStatus :: setSeepageValocity(const FloatArray &w)
{
    seepageVelocity = w;
}

IRResultType
AnisotropicMassTransferMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";       // Required by IR_GIVE_FIELD macro
    IRResultType result;                          // Required by IR_GIVE_FIELD macro

    this->Material :: initializeFrom(ir);

    FloatArray temp;

    IR_GIVE_FIELD(ir, temp, IFT_AnisotropicMassTransferMaterial_c, "c");     // Read permeability matrix c from input file
    k.resize(2, 2);
    k.at(1, 1) = temp.at(1);
    k.at(1, 2) = temp.at(2);
    k.at(2, 1) = temp.at(3);
    k.at(2, 2) = temp.at(4);

    return IRRT_OK;
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
        answer.resize(1, 1);
        answer.at(1, 1) = k.at(1, 1);
    case _2dHeat:
        answer.resize(2, 2);
        answer.at(1, 1) = k.at(1, 1);
        answer.at(1, 2) = k.at(1, 2);
        answer.at(2, 1) = k.at(2, 1);
        answer.at(2, 2) = k.at(2, 2);
        return;

    case _3dHeat:
        _error2( "giveCharacteristicMatrix : 3D capabilitiet not yet implemented (%s)", __MaterialModeToString(mMode) );
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

void
AnisotropicMassTransferMaterial :: giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep)
{
    AnisotropicMassTransferMaterialStatus *thisMaterialStatus;
    thisMaterialStatus = ( ( AnisotropicMassTransferMaterialStatus * ) this->giveStatus(gp) );

    thisMaterialStatus->setPressureGradient(eps);

    answer.beProductOf(k, eps);
    answer.times(-1.0);

    thisMaterialStatus->setSeepageValocity(answer);
}

int
AnisotropicMassTransferMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    AnisotropicMassTransferMaterialStatus *thisMaterialStatus;
    thisMaterialStatus = ( ( AnisotropicMassTransferMaterialStatus * ) this->giveStatus(aGaussPoint) );
    FloatMatrix temp;

    switch ( type ) {
    case IST_Velocity:
        answer = thisMaterialStatus->giveSeepageVelocity();
        break;
    case IST_PressureGradient:
        answer = thisMaterialStatus->giveGradP();
        break;
    default:
      return TransportMaterial :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
    return 1;
}
} // end namespace oofem
