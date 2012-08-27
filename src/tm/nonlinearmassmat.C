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

#include "nonlinearmassmat.h"
#include "domain.h"
#include "flotmtrx.h"
#include "gausspnt.h"
#include <math.h>
#ifndef __MAKEDEPEND
 #include <stdlib.h>
#endif

namespace oofem {
IRResultType
NonlinearMassTransferMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                  // Required by IR_GIVE_FIELD macro

    this->Material :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, C, IFT_NonlinearMassTransferMaterial_c, "c");                     // Constant C
    IR_GIVE_FIELD(ir, alpha, IFT_NonlinearMassTransferMaterial_alpha, "alpha"); // Constant alpha

    return IRRT_OK;
}

void
NonlinearMassTransferMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                                          MatResponseForm form,
                                                          MatResponseMode mode,
                                                          GaussPoint *gp,
                                                          TimeStep *atTime)
{
    MaterialMode mMode = gp->giveMaterialMode();
    AnisotropicMassTransferMaterialStatus *status = ( ( AnisotropicMassTransferMaterialStatus * ) this->giveStatus(gp) );
    FloatArray eps = status->giveGradP();
    FloatArray eps2;
    double gradPNorm;
    FloatMatrix t1, t2, op;

    switch  ( mMode ) {
    case _1dHeat:
        answer.resize(1, 1);
    case _2dHeat:
        answer.resize(2, 2);

        gradPNorm = sqrt( eps.at(1) * eps.at(1) + eps.at(2) * eps.at(2) );

        eps2.resize(2);
        eps2.at(1) = eps.at(1);
        eps2.at(2) = eps.at(2);

        t1.beDyadicProductOf(eps2, eps2);

        if ( gradPNorm != 0.0 ) {
            t1.times( C * alpha * pow(gradPNorm, alpha - 2) );
        }

        t2.resize(2, 2);
        t2.at(1, 1) = t2.at(2, 2) = 1;
        t2.times( 1 + C * pow(gradPNorm, alpha) );

        answer.add(t1);
        answer.add(t2);
        answer.times(1.0);

        return;

    default:
        _error2( "giveCharacteristicMatrix : unknown mode (%s)", __MaterialModeToString(mMode) );
    }
}


double
NonlinearMassTransferMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                                         GaussPoint *gp,
                                                         TimeStep *atTime)
{
    return 0.;
}

void
NonlinearMassTransferMaterial :: giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep)
{
    AnisotropicMassTransferMaterialStatus *thisMaterialStatus;
    thisMaterialStatus = ( ( AnisotropicMassTransferMaterialStatus * ) this->giveStatus(gp) );

    thisMaterialStatus->setPressureGradient(eps);

    double gradPNorm = sqrt( eps.at(1) * eps.at(1) + eps.at(2) * eps.at(2) );

    answer.resize(2);
    answer = eps;
    answer.times( 1 + C * pow(gradPNorm, alpha) );
    answer.times(-1.0);

    thisMaterialStatus->setSeepageValocity(answer);
}

int
NonlinearMassTransferMaterial :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
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
