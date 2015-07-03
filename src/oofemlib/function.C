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
#include "error.h"

namespace oofem {
Function :: Function(int n, Domain *d) :
    FEMComponent(n, d)
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
        return this->evaluateAtTime( tStep->giveTargetTime() ) - this->evaluateAtTime( tStep->giveTargetTime() - tStep->giveTimeIncrement() );
    } else if (mode == VM_Intermediate) {
      return this->evaluateAtTime( tStep->giveIntrinsicTime() );
    } else {
        OOFEM_ERROR("unsupported mode(%d)", mode);
    }

    return 0.;
}


double
Function :: evaluateAtTime(double t)
{
    FloatArray v;
    this->evaluate(v, {{"t",t}});
    if ( v.giveSize() != 1 ) {
        OOFEM_ERROR("Function doesn't return scalar results.");
    }
    return v.at(1);
}


void
Function :: evaluate(FloatArray &answer, const std :: map< std :: string, FunctionArgument > &valDict)
{
    auto it = valDict.find("t");
    if ( it == valDict.end() ) {
        OOFEM_ERROR("Missing necessary argument \"t\"");
    }
    answer = FloatArray{this->evaluateAtTime(it->second.val0)};
}

double
Function :: evaluate(const std :: map< std :: string, FunctionArgument > &valDict)
{
    FloatArray ans;
    this->evaluate(ans, valDict);
    if ( ans.giveSize() != 1 ) {
        OOFEM_ERROR("Function does not return scalar value");
    }
    return ans[0];
}

} // end namespace oofem
