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

#include "periodicpiecewiselinfunction.h"
#include "mathfem.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "domain.h"

namespace oofem {
REGISTER_Function(PeriodicPiecewiseLinFunction);

double PeriodicPiecewiseLinFunction :: evaluateAtTime(double time)
// Returns the value of the receiver at time 'time'. 'time' should be
// one of the dates of the receiver (currently there is no interpola-
// tion between two points).
{
    double add, last;

    if ( !this->dates.giveSize() ) {
        OOFEM_ERROR("Undefined dates and values!");
    }

    if ( addTF && !domain->giveFunction(addTF) ) {
        OOFEM_ERROR("Undefined time function to add!");
    }

    if ( addTF ) {
        add = domain->giveFunction(addTF)->evaluateAtTime(time);
    } else {
        add = 0.;
    }

    // periodicity
    last = dates.at( this->dates.giveSize() ); // time of last date
    if ( ( period >= 0.0 ) && ( time > last ) ) {
        double d = ( time - last ) / period; // periods after last
        time = last + ( d - floor(d) - 1. ) * period;
    }

    return add + PiecewiseLinFunction :: evaluateAtTime(time);
}


double PeriodicPiecewiseLinFunction :: evaluateVelocityAtTime(double time)
// Returns the value of the receiver at time 'time'. 'time' should be
// one of the dates of the receiver (currently there is no interpola-
// tion between two points).
{
    double add, last;

    if ( !this->dates.giveSize() ) {
        OOFEM_ERROR("Undefined dates and values!");
    }

    if ( addTF && !domain->giveFunction(addTF) ) {
        OOFEM_ERROR("Undefined time function to add!");
    }

    if ( addTF ) {
        add = domain->giveFunction(addTF)->evaluateVelocityAtTime(time);
    } else {
        add = 0.;
    }

    // periodicity
    last = dates.at( this->dates.giveSize() ); // time of last date
    if ( ( period >= 0.0 ) && ( time > last ) ) {
        double d = ( time - last ) / period; // periods after last
        time = last + ( d - floor(d) - 1. ) * period;
    }

    return add + PiecewiseLinFunction :: evaluateVelocityAtTime(time);
}

IRResultType
PeriodicPiecewiseLinFunction :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    period = -1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, period, _IFT_PeriodicPiecewiseLinFunction_period);
    addTF = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, addTF, _IFT_PeriodicPiecewiseLinFunction_addtf);

    return PiecewiseLinFunction :: initializeFrom(ir);
}


void PeriodicPiecewiseLinFunction :: giveInputRecord(DynamicInputRecord &input)
{
    PiecewiseLinFunction :: giveInputRecord(input);
    input.setField(this->period, _IFT_PeriodicPiecewiseLinFunction_period);
    input.setField(this->addTF, _IFT_PeriodicPiecewiseLinFunction_addtf);
}
} // end namespace oofem
