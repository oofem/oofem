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

#include "tm/Materials/twophasemat.h"
#include "floatmatrix.h"
#include "function.h"
#include "gausspoint.h"
#include "classfactory.h"
#include "engngm.h"

namespace oofem {
REGISTER_Material(TwoPhaseMaterial);

TwoPhaseMaterial :: TwoPhaseMaterial(int n, Domain *d) : TransportMaterial(n, d)
{
}

TwoPhaseMaterial :: ~TwoPhaseMaterial() {
    // destructor
}


IRResultType
TwoPhaseMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->slaveMaterial, _IFT_TwoPhaseMaterial_mat);
    if ( this->slaveMaterial.giveSize() != 2 ) {
        OOFEM_WARNING("mat array should have two values");
        return IRRT_BAD_FORMAT;
    }

    return IRRT_OK;
}


void
TwoPhaseMaterial :: giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep)
{
    TransportMaterialStatus *ms = static_cast< TransportMaterialStatus * >( this->giveStatus(gp) );
    double vof = this->giveVof(gp, tStep);
    FloatArray v0, v1; 
    this->giveMaterial(0)->giveFluxVector(v0, gp, grad, field, tStep);
    this->giveMaterial(1)->giveFluxVector(v1, gp, grad, field, tStep);

    answer = (1.0 - vof) * v0 + vof * v1;

    ms->setTempFlux(answer);
}


void
TwoPhaseMaterial :: giveCharacteristicMatrix(FloatMatrix &answer,
                                             MatResponseMode mode,
                                             GaussPoint *gp,
                                             TimeStep *tStep)
{
    FloatMatrix v1;
    double vof = this->giveVof(gp, tStep);

    this->giveMaterial(0)->giveCharacteristicMatrix(answer, mode, gp, tStep);
    this->giveMaterial(1)->giveCharacteristicMatrix(v1, mode, gp, tStep);

    answer.times(1.0 - vof);
    v1.times(vof);
    answer.add(v1);
}

double
TwoPhaseMaterial :: giveCharacteristicValue(MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *tStep)
{
    double vof = this->giveVof(gp, tStep);
    return (1.0 - vof) * 
    this->giveMaterial(0)->giveCharacteristicValue(mode, gp, tStep) + vof * this->giveMaterial(1)->giveCharacteristicValue(mode, gp, tStep);
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
