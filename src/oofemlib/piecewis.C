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

#include <fstream>
#include <ios>
#include <sstream>

namespace oofem {
#define PiecewiseLinFunction_PRECISION 1.e-12

PiecewiseLinFunction::PiecewiseLinFunction ( int i, Domain* d ) : LoadTimeFunction ( i, d ), dates(), values()
{
}

double PiecewiseLinFunction :: __at(double time)
// Returns the value of the receiver at time 'time'. 'time' should be
// one of the dates of the receiver (currently there is no interpola-
// tion between two points).
{
    const double precision = PiecewiseLinFunction_PRECISION;
    double xa, xb, ya, yb;

    if ( this->dates.giveSize() == 0 ) {
        _error("at: Undefined dates and values");
    }

    for ( int i = 1; i <= this->dates.giveSize(); i++ ) {
        if ( fabs(this->dates.at(i) - time) < precision ) {
            return this->values.at(i);
        } else if ( this->dates.at(i) > time ) {
            if ( i == 1 ) {
                OOFEM_WARNING3("PiecewiseLinFunction :: __at: computational time %f is out of given time %f, extrapolating value(s)", time, dates.at(i) );
                return 0.;
            }

            xa = this->dates.at(i - 1);
            xb = this->dates.at(i);
            ya = this->values.at(i - 1);
            yb = this->values.at(i);

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

    if ( this->dates.giveSize() == 0 ) {
        _error("at: Undefined dates and values");
    }

    for ( int i = 1; i <= this->dates.giveSize(); i++ ) {
        if ( fabs(dates.at(i) - time) < precision ) {
            if ( i < this->dates.giveSize() ) {
                return ( this->values.at(i + 1) - this->values.at(i) ) / ( this->dates.at(i + 1) - this->dates.at(i) );
            } else {
                return ( this->values.at(i) - this->values.at(i - 1) ) / ( this->dates.at(i) - this->dates.at(i - 1) );
            }
        } else if ( dates.at(i) > time ) {
            if ( i == 1 ) {
                return 0.;
            }

            xa = this->dates.at(i - 1);
            xb = this->dates.at(i);
            ya = this->values.at(i - 1);
            yb = this->values.at(i);

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

    // Optional means, read data from external file (useful for very large sets of data)
    if ( ir->hasField( IFT_PiecewiseLinFunction_dataFile, "datafile") ) {
        std::list< double > t, ft;
        // Open the file;
        std::string fname;
        IR_GIVE_FIELD(ir, fname, IFT_PiecewiseLinFunction_dataFile, "datafile");
        std::ifstream file (fname.c_str(), std::ios::in);
        if ( !file.is_open() ) OOFEM_ERROR2("PieceWiseLinFunction :: initializeFrom - Failed to open data file: %s\n", fname.c_str());
        // Data should be stored in two columns (or just interleaved)
        double temp_t, temp_ft;
        std :: string sLine = "";
        while ( !file.eof() ){
            getline(file, sLine);
            if(sLine[0]=='#'){
                continue;
            }
            std :: stringstream ss1(sLine);
            ss1 >> temp_t >> temp_ft;
            t.push_back(temp_t);
            ft.push_back(temp_ft);
        }
        
        // Copy data over the float arrays
        dates.resize(t.size());
        values.resize(ft.size());
        std::list< double >::iterator it_t = t.begin(), it_ft = ft.begin();
        for ( int i = 1; i <= (int)t.size(); ++i, ++it_t, ++it_ft ) {
            dates.at(i) = *it_t;
            values.at(i) = *it_ft;
        }
    } else {
        int numberOfPoints;
        IR_GIVE_FIELD(ir, numberOfPoints, IFT_PiecewiseLinFunction_npoints, "npoints");
        IR_GIVE_FIELD(ir, dates, IFT_PiecewiseLinFunction_t, "t");
        IR_GIVE_FIELD(ir, values, IFT_PiecewiseLinFunction_ft, "f(t)");
    }

    return IRRT_OK;
}


int
PiecewiseLinFunction :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    LoadTimeFunction :: giveInputRecordString(str, keyword);
    sprintf(buff, " npoints %d", this->dates.giveSize());
    str += buff;
    sprintf( buff, " t %d", this->dates.giveSize() );
    str += buff;
    for ( int i = 1; i <= this->dates.giveSize(); i++ ) {
        sprintf( buff, " %e", this->dates.at(i) );
        str += buff;
    }

    sprintf( buff, " f(t) %d", this->values.giveSize() );
    str += buff;
    for ( int i = 1; i <= this->values.giveSize(); i++ ) {
        sprintf( buff, " %e", this->values.at(i) );
        str += buff;
    }

    return 1;
}

} // end namespace oofem
