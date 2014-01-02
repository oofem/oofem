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

#include "function.h"
#include "timestep.h"
#include "dynamicinputrecord.h"

namespace oofem {
Function :: Function(int n, Domain *d) :
    FEMComponent(n, d),
    initialValue(0.)
{}

double
Function :: evaluate(TimeStep *tStep, ValueModeType mode)
{
    if ( mode == VM_Total ) {
        return this->evaluateAtTime( tStep->giveIntrinsicTime() );
    } else if ( mode == VM_Velocity ) {
        return this->evaluateVelocityAtTime( tStep->giveIntrinsicTime() );
    } else if ( mode == VM_Acceleration ) {
        return this->evaluateAccelerationAtTime( tStep->giveIntrinsicTime() );
    } else if ( mode == VM_Incremental ) {
        //return this->evaluateAtTime( tStep->giveTime() ) - this->evaluateAtTime( tStep->giveTime() - tStep->giveTimeIncrement() );

        if ( tStep->isTheFirstStep() ) {
            return this->evaluateAtTime(tStep->giveIntrinsicTime() - this->initialValue);
        } else {
            return this->evaluateAtTime( tStep->giveIntrinsicTime() ) - this->evaluateAtTime( tStep->giveIntrinsicTime() - tStep->giveTimeIncrement() );
        }
    } else {
        _error2("Function:: evaluate: unsupported mode(%d)", mode);
    }

    return 0.;
}


IRResultType
Function :: initializeFrom(InputRecord *ir)
{
    //
    // instanciates receiver according to input record
    //
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;              // Required by IR_GIVE_FIELD macro


    this->initialValue = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, this->initialValue, _IFT_Function_initialvalue);

    return IRRT_OK;
}


void
Function :: giveInputRecord(DynamicInputRecord &input)
{
    FEMComponent :: giveInputRecord(input);
    input.setField(this->initialValue, _IFT_Function_initialvalue);
}


double
Function :: evaluateAtTime(double t)
{
    //std::map< std::string, double >valDict {{t, "t"}};
    std::map< std::string, double >valDict;
    valDict["t"] = t;
    FloatArray v;
    this->evaluate(v, valDict);
    if ( v.giveSize() != 1 ) {
        OOFEM_ERROR2("%s :: evaluateAtTime - Function doesn't return scalar results.", this->giveClassName());
    }
    return v.at(1);
}

///@todo Move operator from C++11 would be nice here.
void
Function :: evaluate(FloatArray &answer, std::map< std::string, double > &valDict)
{
    answer.resize(1);
    answer.at(1) = this->evaluateAtTime(valDict["t"]);
}

} // end namespace oofem
