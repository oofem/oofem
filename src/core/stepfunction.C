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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "stepfunction.h"
#include "mathfem.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

#include <fstream>
#include <sstream>

namespace oofem {
#define StepFunction_PRECISION 1.e-12

REGISTER_Function(StepFunction);

StepFunction :: StepFunction(int i, Domain *d) : Function(i, d), datevalues()
{ }

double StepFunction :: evaluateAtTime(double time)
/* The function is defined as:
 * @f[
 * f(t) = \begin{cases}
 * v_1, & t < t_2 \\
 * v_i, & t_i \leq t < t_{i+1}, i= 2, \ldots, n-1 \\
 * v_n, & t \geq t_n
 * \end{cases}
 * @f]
 * where @f$ (t_i, v_i), i=1,\ldots,n @f$ are the points stored in 'dates' and 'values'.
 */
{
    int size = (int) this->datevalues.size();
    if ( size == 0 ) {
        OOFEM_ERROR("Undefined dates and values");
    } else if ( size == 1) {
        return std::get<1>(this->datevalues.at(0));
    }
    // fast track 
    if (time < std::get<0>(this->datevalues.at(1))) {
        return std::get<1>(this->datevalues.at(0));
    } else if (time > std::get<0>(this->datevalues.at(size-1))) {
        return std::get<1>(this->datevalues.at(size-1));
    }

    // search
    auto it = std::upper_bound(
        datevalues.begin(), datevalues.end(), std::make_tuple(time, 0.0),
        [](const std::tuple<double, double>& a, const std::tuple<double, double>& b) {
            return std::get<0>(a) < std::get<0>(b);
        }
    );
    
    if (it != datevalues.end()) {
        return std::get<1>(*(it - 1));
    } else {
        return std::get<1>(this->datevalues.at(size-1));
    }

    /*
    const double precision = StepFunction_PRECISION;
    int size = (int) this->dates.size();

    if ( size == 0 ) {
        OOFEM_ERROR("Undefined dates and values");
    } else if ( size == 1) {
        return this->values.at(1);
    } else {
        for ( int i = 2; i <= size; i++ ) {
            double date = this->dates.at(i);
            if ( fabs(date - time) < precision ) {
                return this->values.at(i);
            } else if ( date > time ) {
                return this->values.at(i - 1);
            }
        }
        return this->values.at(size);
    }
    */
}

double StepFunction :: evaluateVelocityAtTime(double time)
// Returns the derivative of the receiver at time 'time'. 'time' should be
// one of the dates of the receiver (currently there is no interpola-
// tion between two points).
{
    return 0.;
}

void
StepFunction :: initializeFrom(InputRecord &ir)
{
    Function :: initializeFrom(ir);

    // Optional means, read data from external file (useful for very large sets of data)
    if ( ir.hasField(_IFT_StepFunction_dataFile) ) {
        datevalues.clear();
        // Open the file;
        std :: string fname;
        IR_GIVE_FIELD(ir, fname, _IFT_StepFunction_dataFile);
        if (fname != ""){ //allows empty filename when the function will not be used intentionally
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
                datevalues.push_back(std::tuple<double, double>(temp_t, temp_ft));
            }
            file.close();
        }
    } else {
        int numberOfPoints;
        IR_GIVE_OPTIONAL_FIELD(ir, numberOfPoints, "npoints");
        FloatArray dates, values;
        IR_GIVE_FIELD(ir, dates, _IFT_StepFunction_t);
        IR_GIVE_FIELD(ir, values, _IFT_StepFunction_ft);
        if ( dates.size() != values.size() ) {
            OOFEM_ERROR ("Size of t and f(t) arrays must be equal");
        }
        // copy data to internal structure
        datevalues.clear();
        for ( int i = 0; i < dates.size(); i++ ) {
            datevalues.push_back(std::tuple<double, double>(dates(i), values(i)));
        }
    }
}


void StepFunction :: giveInputRecord(DynamicInputRecord &input)
{
    Function :: giveInputRecord(input);

    int numberOfPoints = (int) this->datevalues.size();
    input.setField(numberOfPoints, "npoints");
    FloatArray dates(numberOfPoints), values(numberOfPoints);
    for ( int i = 0; i < numberOfPoints; i++ ) {
        dates(i) = std::get<0>(this->datevalues.at(i));
        values(i) = std::get<1>(this->datevalues.at(i));
    }
    input.setField(dates, _IFT_StepFunction_t);
    input.setField(values, _IFT_StepFunction_ft);
}



void
StepFunction :: saveContext(DataStream &stream, ContextMode mode)
{
    if ( mode & CM_Definition ) {
        int numberOfPoints = (int) this->datevalues.size();
        FloatArray dates(numberOfPoints), values(numberOfPoints);
        for ( int i = 0; i < numberOfPoints; i++ ) {
            dates(i) = std::get<0>(this->datevalues.at(i));
            values(i) = std::get<1>(this->datevalues.at(i));
        }
        dates.storeYourself(stream);
        values.storeYourself(stream);
    }
}


void
StepFunction :: restoreContext(DataStream &stream, ContextMode mode)
{
    if ( mode & CM_Definition ) {
        FloatArray dates, values;
        dates.restoreYourself(stream);
        values.restoreYourself(stream);
        if ( dates.size() != values.size() ) {
            OOFEM_ERROR ("Size of t and f(t) arrays must be equal");
        }
        // copy data to internal structure
        datevalues.clear();
        for ( int i = 0; i < dates.size(); i++ ) {
            datevalues.push_back(std::tuple<double, double>(dates(i), values(i)));
        }
    }
}
} // end namespace oofem
