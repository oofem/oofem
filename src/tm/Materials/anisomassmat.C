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

#include "tm/Materials/anisomassmat.h"
#include "floatmatrix.h"
#include "floatmatrixf.h"
#include "gausspoint.h"
#include "classfactory.h"

#include <cstdlib>

namespace oofem {
REGISTER_Material(AnisotropicMassTransferMaterial);

IRResultType
AnisotropicMassTransferMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                          // Required by IR_GIVE_FIELD macro

    FloatArray c;
    IR_GIVE_FIELD(ir, c, _IFT_AnisotropicMassTransferMaterial_c);     // Read permeability matrix c from input file
    if ( c.giveSize() == 4 ) {
        k = {
            c[0], c[1], 0.,
            c[2], c[3], 0.,
            0., 0., 0.,
        };
    } else if ( c.giveSize() == 9 ) {
        k = {
            c[0], c[1], c[2],
            c[3], c[4], c[5],
            c[6], c[7], c[8],
        };
    }
    return TransportMaterial :: initializeFrom(ir);
}


FloatArrayF<3>
AnisotropicMassTransferMaterial :: computeFlux3D(const FloatArrayF<3> &grad, double field, GaussPoint *gp, TimeStep *tStep) const
{
    TransportMaterialStatus *ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    auto answer = - dot(k, grad);
    ms->setTempField(field);
    ms->setTempGradient(grad);
    ms->setTempFlux(answer);
    return answer;
}


FloatMatrixF<3,3>
AnisotropicMassTransferMaterial :: computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    return k;
}


double
AnisotropicMassTransferMaterial :: giveCharacteristicValue(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("unknown mode (%s)", __MatResponseModeToString(mode) );

    return 0.;
}
} // end namespace oofem
