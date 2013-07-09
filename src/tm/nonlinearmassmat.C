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

#include "nonlinearmassmat.h"
#include "domain.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {

REGISTER_Material( NonlinearMassTransferMaterial );

IRResultType
NonlinearMassTransferMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    this->Material :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, C, _IFT_NonlinearMassTransferMaterial_c);
    IR_GIVE_FIELD(ir, alpha, _IFT_NonlinearMassTransferMaterial_alpha);

    return IRRT_OK;
}

void
NonlinearMassTransferMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                                          MatResponseMode mode,
                                                          GaussPoint *gp,
                                                          TimeStep *atTime)
{
    MaterialMode mMode = gp->giveMaterialMode();
    TransportMaterialStatus *status = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    FloatArray eps = status->giveTempGradient();
    double gradPNorm;
    FloatMatrix t1, t2;

    gradPNorm = eps.computeNorm();

    t1.beDyadicProductOf(eps, eps);
    if ( gradPNorm != 0.0 ) {
        t1.times( C * alpha * pow(gradPNorm, alpha - 2) );
    }

    switch  ( mMode ) {
    case _1dHeat:
        t2.resize(1, 1);
        t2.at(1, 1) = 1;
        break;
    case _2dHeat:
        t2.resize(2, 2);
        t2.at(1, 1) = t2.at(2, 2) = 1;
        break;
    case _3dHeat:
        t2.resize(3, 3);
        t2.at(1, 1) = t2.at(2, 2) = t2.at(3, 3) = 1;
        break;
    default:
        _error2( "giveCharacteristicMatrix : unknown mode (%s)", __MaterialModeToString(mMode) );
    }

    answer.beEmptyMtrx();
    answer.add(t1);
    answer.add(1 + C * pow(gradPNorm, alpha), t2);

}

double
NonlinearMassTransferMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                                         GaussPoint *gp,
                                                         TimeStep *atTime)
{
    return 0.;
}

void
NonlinearMassTransferMaterial :: giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep)
{
    TransportMaterialStatus *ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );

    double gradPNorm = grad.computeNorm();
    answer.beScaled( -(1. + C * pow(gradPNorm, alpha)), grad);

    ms->setTempGradient(grad);
    ms->setTempField(field);
    ms->setTempFlux(answer);
}

int
NonlinearMassTransferMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    TransportMaterialStatus *ms = static_cast< TransportMaterialStatus * >( this->giveStatus(aGaussPoint) );
 
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
