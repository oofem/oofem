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

#include "piecewisper.h"
#include "mathfem.h"

namespace oofem {
double PeriodicPiecewiseLinFunction :: __at(double time)
// Returns the value of the receiver at time 'time'. 'time' should be
// one of the dates of the receiver (currently there is no interpola-
// tion between two points).
{
    double add, last;

    if ( !numberOfPoints ) {
        _error("at: Undefined dates and values!");
    }

    if ( addTF && !domain->giveLoadTimeFunction(addTF) ) {
        _error("at: Undefined time function to add!");
    }

    if ( addTF ) {
        add = domain->giveLoadTimeFunction(addTF)->__at(time);
    } else {
        add = 0.;
    }

    // periodicity
    last = dates.at(numberOfPoints); // time of last date
    if ( ( period >= 0.0 ) && ( time > last ) ) {
        double d = ( time - last ) / period; // periods after last
        time = last + ( d - floor(d) - 1. ) * period;
    }

    return add + PiecewiseLinFunction :: __at(time);
}


double PeriodicPiecewiseLinFunction :: __derAt(double time)
// Returns the value of the receiver at time 'time'. 'time' should be
// one of the dates of the receiver (currently there is no interpola-
// tion between two points).
{
    double add, last;

    if ( !numberOfPoints ) {
        _error("derAt: Undefined dates and values!");
    }

    if ( addTF && !domain->giveLoadTimeFunction(addTF) ) {
        _error("derAt: Undefined time function to add!");
    }

    if ( addTF ) {
        add = domain->giveLoadTimeFunction(addTF)->__derAt(time);
    } else {
        add = 0.;
    }

    // periodicity
    last = dates.at(numberOfPoints); // time of last date
    if ( ( period >= 0.0 ) && ( time > last ) ) {
        double d = ( time - last ) / period; // periods after last
        time = last + ( d - floor(d) - 1. ) * period;
    }

    return add + PiecewiseLinFunction :: __derAt(time);
}

IRResultType
PeriodicPiecewiseLinFunction :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    PiecewiseLinFunction :: initializeFrom(ir);

    period = -1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, period, IFT_PeriodicPiecewiseLinFunction_period, "period"); // Macro
    addTF = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, addTF, IFT_PeriodicPiecewiseLinFunction_addtf, "addtf"); // Macroo

    return IRRT_OK;
}


int
PeriodicPiecewiseLinFunction :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    PiecewiseLinFunction :: giveInputRecordString(str, keyword);
    sprintf(buff, " period %e", this->period);
    str += buff;
    sprintf(buff, " addTF %d", this->addTF);
    str += buff;

    return 1;
}
} // end namespace oofem
