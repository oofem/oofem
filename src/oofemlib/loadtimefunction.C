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

/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

#include "loadtimefunction.h"
#include "timestep.h"
#include "dynamicinputrecord.h"

namespace oofem {

LoadTimeFunction :: LoadTimeFunction ( int n, Domain* d ) :
    FEMComponent( n, d ),
    initialValue( 0.)
{
}

double
LoadTimeFunction :: evaluate(TimeStep *atTime, ValueModeType mode)
{
    if ( mode == VM_Total ) {
        return this->__at( atTime->giveIntrinsicTime() );
    } else if ( mode == VM_Velocity ) {
        return this->__derAt( atTime->giveIntrinsicTime() );
    } else if ( mode == VM_Acceleration ) {
        return this->__accelAt( atTime->giveIntrinsicTime() );
    } else if ( mode == VM_Incremental ) {
        //return this->__at( atTime->giveTime() ) - this->__at( atTime->giveTime() - atTime->giveTimeIncrement() );

        if ( atTime->isTheFirstStep() ) {
            return this->__at(atTime->giveIntrinsicTime() - this->initialValue);
        } else {
            return this->__at( atTime->giveIntrinsicTime() ) - this->__at( atTime->giveIntrinsicTime() - atTime->giveTimeIncrement() );
        }
    } else {
        _error2("LoadTimeFunction:: evaluate: unsupported mode(%d)", mode);
    }

    return 0.;
}


IRResultType
LoadTimeFunction :: initializeFrom(InputRecord *ir)
{
    //
    // instanciates receiver according to input record
    //
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;              // Required by IR_GIVE_FIELD macro


    this->initialValue = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, this->initialValue, _IFT_LoadTimeFunction_initialvalue);

    return IRRT_OK;
}


void
LoadTimeFunction :: giveInputRecord(DynamicInputRecord& input)
{
    FEMComponent :: giveInputRecord(input);
    input.setField(this->initialValue, _IFT_LoadTimeFunction_initialvalue);
}


} // end namespace oofem
