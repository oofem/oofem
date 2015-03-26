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

#include "isolinmoisturemat.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(IsotropicLinMoistureTransferMaterial);

IRResultType
IsotropicLinMoistureTransferMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, permeability, _IFT_IsotropicLinMoistureTransferMaterial_perm);
    IR_GIVE_FIELD(ir, moistureCapacity, _IFT_IsotropicLinMoistureTransferMaterial_capa);

    return IsotropicMoistureTransferMaterial :: initializeFrom(ir);
}


double
IsotropicLinMoistureTransferMaterial :: giveMoistureCapacity(GaussPoint *gp, TimeStep *tStep)
{
    return this->moistureCapacity;
}

double
IsotropicLinMoistureTransferMaterial :: givePermeability(GaussPoint *gp, TimeStep *tStep)
{
    return this->permeability;
}
} // end namespace oofem
