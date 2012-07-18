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

/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

#include "loadtime.h"
#include "timestep.h"

namespace oofem {
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


    IR_GIVE_OPTIONAL_FIELD(ir, initialValue, IFT_LoadTimeFunction_initialvalue, "initialvalue"); // Macro

    return IRRT_OK;
}


IRResultType
LoadTimeFunction :: initializeFrom(InputRecord *ir, DataReader *dr)
{
    return initializeFrom(ir);
}


int
LoadTimeFunction :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    sprintf(buff, " initialvalue %e", this->initialValue);
    str += buff;

    return 1;
}
} // end namespace oofem
