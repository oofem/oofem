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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "tm/Materials/twophasemat.h"
#include "floatmatrix.h"
#include "function.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "engngm.h"

namespace oofem {
REGISTER_Material(TwoPhaseMaterial);

void
TwoPhaseMaterial :: initializeFrom(InputRecord &ir)
{
    IR_GIVE_FIELD(ir, this->slaveMaterial, _IFT_TwoPhaseMaterial_mat);
    if ( this->slaveMaterial.giveSize() != 2 ) {
        throw ValueInputException(ir, _IFT_TwoPhaseMaterial_mat, "must have two values");
    }
}


FloatArrayF<3>
TwoPhaseMaterial :: computeFlux3D(const FloatArrayF<3> &grad, double field, GaussPoint *gp, TimeStep *tStep) const
{
    auto ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    ms->setTempField(field);
    ms->setTempGradient(grad);

    double vof = this->giveVof(gp, tStep);
    auto v0 = this->giveMaterial(0)->computeFlux3D(grad, field, gp, tStep);
    auto v1 = this->giveMaterial(1)->computeFlux3D(grad, field, gp, tStep);
    auto answer = (1.0 - vof) * v0 + vof * v1;
    ms->setTempFlux(answer);
    return answer;
}


FloatMatrixF<3,3>
TwoPhaseMaterial :: computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    double vof = this->giveVof(gp, tStep);
    auto v0 = this->giveMaterial(0)->computeTangent3D(mode, gp, tStep);
    auto v1 = this->giveMaterial(1)->computeTangent3D(mode, gp, tStep);
    return (1.0 - vof) * v0 + vof * v1;
}

double
TwoPhaseMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *tStep) const
{
    double vof = this->giveVof(gp, tStep);
    auto v0 = this->giveMaterial(0)->giveCharacteristicValue(mode, gp, tStep);
    auto v1 = this->giveMaterial(1)->giveCharacteristicValue(mode, gp, tStep);
    return (1.0 - vof) * v0 + vof * v1;
}

TransportMaterial *
TwoPhaseMaterial :: giveMaterial(int i) const
{
    return static_cast< TransportMaterial * >( domain->giveMaterial( slaveMaterial[i] ) );
}

double 
TwoPhaseMaterial :: giveVof (GaussPoint* gp, TimeStep* tStep) const 
{
    return static_cast< TransportElement * >(gp->giveElement()) ->computeVof(tStep);
}


} // end namespace oofem
