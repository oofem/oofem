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

void
AnisotropicMassTransferMaterial :: initializeFrom(InputRecord &ir)
{
    TransportMaterial :: initializeFrom(ir);

    FloatMatrix c;
    IR_GIVE_FIELD(ir, c, _IFT_AnisotropicMassTransferMaterial_c);     // Read permeability matrix c from input file
    if ( c.giveNumberOfColumns() != 3 && c.giveNumberOfRows() != 3 ) {
        throw ValueInputException(ir, _IFT_AnisotropicMassTransferMaterial_c, "c must be a 3x3");
    }
    k = c;
}


FloatArrayF<3>
AnisotropicMassTransferMaterial :: computeFlux3D(const FloatArrayF<3> &grad, double field, GaussPoint *gp, TimeStep *tStep) const
{
    auto ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
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
}
} // end namespace oofem
