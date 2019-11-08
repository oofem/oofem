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

#include "tm/Materials/nonlinearmassmat.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"
#include <floatmatrixf.h>

namespace oofem {
REGISTER_Material(NonlinearMassTransferMaterial);

void
NonlinearMassTransferMaterial :: initializeFrom(InputRecord &ir)
{
    Material :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, C, _IFT_NonlinearMassTransferMaterial_c);
    IR_GIVE_FIELD(ir, alpha, _IFT_NonlinearMassTransferMaterial_alpha);
}

double
NonlinearMassTransferMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                                         GaussPoint *gp,
                                                         TimeStep *tStep) const
{
    return 0.;
}

FloatArrayF<3>
NonlinearMassTransferMaterial :: computeFlux3D(const FloatArrayF<3> &grad, double field, GaussPoint *gp, TimeStep *tStep) const
{
    auto ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    ms->setTempGradient(grad);
    ms->setTempField(field);
    
    double gradPNorm = norm(grad);
    auto answer = -( 1. + C * pow(gradPNorm, alpha) ) * grad;
    ms->setTempFlux(answer);
    return answer;
}


FloatMatrixF<3,3>
NonlinearMassTransferMaterial :: computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    auto ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    const auto &eps = ms->giveTempGradient();
    double gradPNorm = norm(eps);

    auto scale = gradPNorm != 0. ? C * alpha * pow(gradPNorm, alpha - 2) : 1.;
    auto t1 = scale * dyad(eps, eps);
    return t1 + (1 + C * pow(gradPNorm, alpha)) * eye<3>();
}

int
NonlinearMassTransferMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    auto ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );

    switch ( type ) {
    case IST_Velocity:
        answer = ms->giveFlux();
        break;
    case IST_PressureGradient:
        answer = ms->giveGradient();
        break;
    default:
        return TransportMaterial :: giveIPValue(answer, gp, type, tStep);
    }
    return 1;
}
} // end namespace oofem
