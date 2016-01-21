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

#include "piecewiselinfunction.h"
#include "mathfem.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

#include <fstream>
#include <sstream>

namespace oofem {
#define PiecewiseLinFunction_PRECISION 1.e-12

REGISTER_Function(PiecewiseLinFunction);

PiecewiseLinFunction :: PiecewiseLinFunction(int i, Domain *d) : Function(i, d), dates(), values()
{ }

double PiecewiseLinFunction :: evaluateAtTime(double time)
// Returns the value of the receiver at time 'time'. 'time' should be
// one of the dates of the receiver (currently there is no interpola-
// tion between two points).
{
    const double precision = PiecewiseLinFunction_PRECISION;
    double xa, xb, ya, yb;

    if ( this->dates.giveSize() == 0 ) {
        OOFEM_ERROR("Undefined dates and values");
    }

    for ( int i = 1; i <= this->dates.giveSize(); i++ ) {
        if ( fabs(this->dates.at(i) - time) < precision ) {
            return this->values.at(i);
        } else if ( this->dates.at(i) > time ) {
            if ( i == 1 ) {
                OOFEM_WARNING("computational time %f is out of given time %f, using closest value", time, dates.at(i) );
                return this->dates.at(i);
            }

            xa = this->dates.at(i - 1);
            xb = this->dates.at(i);
            ya = this->values.at(i - 1);
            yb = this->values.at(i);

            return ya + ( time - xa ) * ( yb - ya ) / ( xb - xa );
        }
    }

    OOFEM_WARNING("computational time %f is out of given time, using closest value", time );
    return this->values.at(this->values.giveSize());
}

double PiecewiseLinFunction :: evaluateVelocityAtTime(double time)
// Returns the derivative of the receiver at time 'time'. 'time' should be
// one of the dates of the receiver (currently there is no interpola-
// tion between two points).
{
    const double precision = PiecewiseLinFunction_PRECISION;
    double xa, xb, ya, yb;

    if ( this->dates.giveSize() == 0 ) {
        OOFEM_ERROR("Undefined dates and values");
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
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // Optional means, read data from external file (useful for very large sets of data)
    if ( ir->hasField(_IFT_PiecewiseLinFunction_dataFile) ) {
        std :: list< double >t, ft;
        // Open the file;
        std :: string fname;
        IR_GIVE_FIELD(ir, fname, _IFT_PiecewiseLinFunction_dataFile);
        std :: ifstream file(fname.c_str(), std :: ios :: in);
        if ( !file.is_open() ) {
            OOFEM_ERROR("Failed to open data file: %s\n", fname.c_str());
        }
        // Data should be stored in two columns (or just interleaved)
        double temp_t, temp_ft;
        std :: string sLine = "";
        while ( !file.eof() ) {
            getline(file, sLine);
            if ( sLine [ 0 ] == '#' ) {
                continue;
            }
            std :: stringstream ss1(sLine);
            ss1 >> temp_t >> temp_ft;
            t.push_back(temp_t);
            ft.push_back(temp_ft);
        }

        // Copy data over the float arrays
        dates.resize( t.size() );
        values.resize( ft.size() );
        std :: list< double > :: iterator it_t = t.begin(), it_ft = ft.begin();
        for ( int i = 1; i <= ( int ) t.size(); ++i, ++it_t, ++it_ft ) {
            dates.at(i) = * it_t;
            values.at(i) = * it_ft;
        }
    } else {
        int numberOfPoints;
        IR_GIVE_OPTIONAL_FIELD(ir, numberOfPoints, "npoints");
        IR_GIVE_FIELD(ir, dates, _IFT_PiecewiseLinFunction_t);
        IR_GIVE_FIELD(ir, values, _IFT_PiecewiseLinFunction_ft);
    }

    return Function :: initializeFrom(ir);
}


void PiecewiseLinFunction :: giveInputRecord(DynamicInputRecord &input)
{
    Function :: giveInputRecord(input);
    input.setField(this->dates, _IFT_PiecewiseLinFunction_t);
    input.setField(this->values, _IFT_PiecewiseLinFunction_ft);
}



contextIOResultType
PiecewiseLinFunction :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    if ( mode & CM_Definition ) {
        dates.storeYourself(stream);
        values.storeYourself(stream);
    }

    return CIO_OK;
}


contextIOResultType
PiecewiseLinFunction :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    if ( mode & CM_Definition ) {
        dates.restoreYourself(stream);
        values.restoreYourself(stream);
    }

    return CIO_OK;
}
} // end namespace oofem
