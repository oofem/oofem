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

#include "usrdeftempfield.h"
#include "timestep.h"
#include "classfactory.h"

#include <sstream>

namespace oofem {
REGISTER_BoundaryCondition(UserDefinedTemperatureField);

void
UserDefinedTemperatureField :: computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode)
// Returns the value of the receiver at time and given position respecting the mode.
{
    int err;
    double result;

    if ( ( mode != VM_Incremental ) && ( mode != VM_Total ) ) {
        OOFEM_ERROR("unknown mode (%s)", __ValueModeTypeToString(mode) );
    }

    answer.resize(this->size);
    std :: ostringstream buff;
    for ( int i = 1; i <= size; i++ ) {
        buff << "x=" << coords.at(1) << ";y=" << coords.at(2) << ";z=" << coords.at(3) <<
        ";t=" << tStep->giveTargetTime() << ";" << ftExpression [ i - 1 ];
        result = myParser.eval(buff.str().c_str(), err);
        if ( err ) {
            OOFEM_ERROR("parser syntax error");
        }

        answer.at(i) = result;

        if ( ( mode == VM_Incremental ) && ( !tStep->isTheFirstStep() ) ) {
            buff << "x=" << coords.at(1) << ";y=" << coords.at(2) << ";z=" << coords.at(3) <<
            ";t=" << ( tStep->giveTargetTime() - tStep->giveTimeIncrement() ) << ";" << ftExpression [ i - 1 ];
            result = myParser.eval(buff.str().c_str(), err);
            if ( err ) {
                OOFEM_ERROR("parser syntax error");
            }

            answer.at(i) -= result;
        }
    }
}

IRResultType
UserDefinedTemperatureField :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, size, _IFT_UserDefinedTemperatureField_size);
    if ( size > 3 ) {
        size = 3;
    }

    if ( size > 0 ) {
        IR_GIVE_FIELD(ir, ftExpression [ 0 ], _IFT_UserDefinedTemperatureField_t1);
    }

    if ( size > 1 ) {
        IR_GIVE_FIELD(ir, ftExpression [ 1 ], _IFT_UserDefinedTemperatureField_t2);
    }

    if ( size > 2 ) {
        IR_GIVE_FIELD(ir, ftExpression [ 2 ], _IFT_UserDefinedTemperatureField_t3);
    }

    return IRRT_OK;
}
} // end namespace oofem
