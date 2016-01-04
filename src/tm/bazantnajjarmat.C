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

#include "bazantnajjarmat.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(BazantNajjarMoistureTransferMaterial);

IRResultType
BazantNajjarMoistureTransferMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, C1, _IFT_BazantNajjarMoistureTransferMaterial_c1);
    IR_GIVE_FIELD(ir, n, _IFT_BazantNajjarMoistureTransferMaterial_n);
    IR_GIVE_FIELD(ir, alpha0, _IFT_BazantNajjarMoistureTransferMaterial_alpha0);
    IR_GIVE_FIELD(ir, hC, _IFT_BazantNajjarMoistureTransferMaterial_hc);

    this->moistureCapacity = 1.;
    IR_GIVE_OPTIONAL_FIELD(ir, moistureCapacity, _IFT_BazantNajjarMoistureTransferMaterial_capa);

    return IsotropicMoistureTransferMaterial :: initializeFrom(ir);
}


double
BazantNajjarMoistureTransferMaterial :: giveMoistureCapacity(GaussPoint *gp, TimeStep *tStep)
{
    return this->moistureCapacity;
}

double
BazantNajjarMoistureTransferMaterial :: givePermeability(GaussPoint *gp, TimeStep *tStep)
{
    double permeability;
    double humidity = this->giveHumidity(gp, VM_Total);

    humidity = max(humidity, 0.0);
    humidity = min(humidity, 1.0);

    permeability = C1 * ( alpha0 + ( 1. - alpha0 ) / ( 1. + pow( ( 1. - humidity ) / ( 1. - hC ), n ) ) );
    return permeability;
}

double
BazantNajjarMoistureTransferMaterial :: giveHumidity(GaussPoint *gp, ValueModeType mode)
{
    FloatArray tempState = static_cast< TransportMaterialStatus * >( giveStatus(gp) )->giveTempField();
    if ( ( tempState.at(1) > 1.0 ) || ( tempState.at(1) < 0.0 ) ) {
         OOFEM_WARNING("Relative humidity %.5f is out of range 0.0 - 1.0", tempState.at(1) );
    }
    return tempState.at(1);
}
} // end namespace oofem
