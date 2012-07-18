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

#include "piecewis.h"
#include "mathfem.h"

namespace oofem {
#define PiecewiseLinFunction_PRECISION 1.e-12

double PiecewiseLinFunction :: __at(double time)
// Returns the value of the receiver at time 'time'. 'time' should be
// one of the dates of the receiver (currently there is no interpola-
// tion between two points).
{
    const double precision = PiecewiseLinFunction_PRECISION;
    double xa, xb, ya, yb;
    int i;

    if ( !numberOfPoints ) {
        //      this -> getPoints() ;
        _error("at: Undefined dates and values");
    }

    for ( i = 1; i <= numberOfPoints; i++ ) {
        if ( fabs(dates.at(i) - time) < precision ) {
            return values.at(i);
        } else if ( dates.at(i) > time ) {
            if ( i == 1 ) {
                OOFEM_WARNING3("PiecewiseLinFunction :: __at: computational time %f is out of given time %f, extrapolating value(s)", time, dates.at(i) );
                return 0.;
            }

            xa = dates.at(i - 1);
            xb = dates.at(i);
            ya = values.at(i - 1);
            yb = values.at(i);

            return ya + ( time - xa ) * ( yb - ya ) / ( xb - xa );
        }
    }

    return 0.;
}

double PiecewiseLinFunction :: __derAt(double time)
// Returns the derivative of the receiver at time 'time'. 'time' should be
// one of the dates of the receiver (currently there is no interpola-
// tion between two points).
{
    const double precision = PiecewiseLinFunction_PRECISION;
    double xa, xb, ya, yb;
    int i;

    if ( !numberOfPoints ) {
        //      this -> getPoints() ;
        _error("at: Undefined dates and values");
    }

    for ( i = 1; i <= numberOfPoints; i++ ) {
        if ( fabs(dates.at(i) - time) < precision ) {
            if ( i < numberOfPoints ) {
                return ( values.at(i + 1) - values.at(i) ) / ( dates.at(i + 1) - dates.at(i) );
            } else {
                return ( values.at(i) - values.at(i - 1) ) / ( dates.at(i) - dates.at(i - 1) );
            }
        } else if ( dates.at(i) > time ) {
            if ( i == 1 ) {
                return 0.;
            }

            xa = dates.at(i - 1);
            xb = dates.at(i);
            ya = values.at(i - 1);
            yb = values.at(i);

            return ( yb - ya ) / ( xb - xa );
        }
    }

    return 0.;
}

IRResultType
PiecewiseLinFunction :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    LoadTimeFunction :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, numberOfPoints, IFT_PiecewiseLinFunction_npoints, "npoints"); // Macro

    IR_GIVE_FIELD(ir, dates, IFT_PiecewiseLinFunction_t, "t"); // Macro
    IR_GIVE_FIELD(ir, values, IFT_PiecewiseLinFunction_ft, "f(t)"); // Macro

    return IRRT_OK;
}


int
PiecewiseLinFunction :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];
    int i;

    LoadTimeFunction :: giveInputRecordString(str, keyword);
    sprintf(buff, " npoints %d", this->numberOfPoints);
    str += buff;
    sprintf( buff, " t %d", this->dates.giveSize() );
    str += buff;
    for ( i = 1; i <= this->dates.giveSize(); i++ ) {
        sprintf( buff, " %e", this->dates.at(i) );
        str += buff;
    }

    sprintf( buff, " f(t) %d", this->values.giveSize() );
    str += buff;
    for ( i = 1; i <= this->values.giveSize(); i++ ) {
        sprintf( buff, " %e", this->values.at(i) );
        str += buff;
    }

    return 1;
}
} // end namespace oofem
