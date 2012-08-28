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
#include "mathfem.h"

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
NonlinearMassTransferMaterial :: giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep)
{
    AnisotropicMassTransferMaterialStatus *thisMaterialStatus;
    thisMaterialStatus = ( ( AnisotropicMassTransferMaterialStatus * ) this->giveStatus(gp) );

    thisMaterialStatus->setPressureGradient(eps);

    double gradPNorm = eps.computeNorm();

    answer = eps;
    answer.times( 1 + C * pow(gradPNorm, alpha) );
    answer.negated();

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
