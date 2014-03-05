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

#include "function.h"
#include "timestep.h"
#include "dynamicinputrecord.h"

namespace oofem {
Function :: Function(int n, Domain *d) :
    FEMComponent(n, d),
    initialValue(0.)
{ }

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
        OOFEM_ERROR("unsupported mode(%d)", mode);
    }

    return 0.;
}


IRResultType
Function :: initializeFrom(InputRecord *ir)
{
    //
    // instanciates receiver according to input record
    //
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
    std :: map< std :: string, FunctionArgument >valDict;
    valDict.insert( std :: make_pair("t", t) );
    FloatArray v;
    this->evaluate(v, valDict);
    ///@todo This should be possible and nice to use if we have C++11
    //this->evaluate(v, {{t, "t"}});
    if ( v.giveSize() != 1 ) {
        OOFEM_ERROR("Function doesn't return scalar results.");
    }
    return v.at(1);
}


void
Function :: evaluate(FloatArray &answer, std :: map< std :: string, FunctionArgument > &valDict)
{
    std :: map< std :: string, FunctionArgument > :: iterator it = valDict.find("t");
#ifdef DEBUG
    if ( it == valDict.end() ) {
        OOFEM_ERROR("Missing necessary argument \"t\"");
    }
#endif
    answer.resize(1);
    answer.at(1) = this->evaluateAtTime(it->second.val0);
}
} // end namespace oofem
