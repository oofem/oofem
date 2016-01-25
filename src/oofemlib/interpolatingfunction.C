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

#include "interpolatingfunction.h"
#include "floatarray.h"
#include "classfactory.h"
#include "error.h"

#include <iostream>
#include <fstream>

namespace oofem {
REGISTER_Function(InterpolatingFuction);

InterpolatingFuction :: InterpolatingFuction(int num, Domain *d) : Function(num, d)
{ }

InterpolatingFuction :: ~InterpolatingFuction()
{ }


void
InterpolatingFuction :: evaluate(FloatArray &answer, const std :: map< std :: string, FunctionArgument > &valDict)
{
    auto it = valDict.find("x");
    if ( it == valDict.end() ) {
        OOFEM_ERROR("Coordinate needed for evaluating function");
    }
    const FloatArray &globalCoordinates = it->second.val1;
    // Map the corresponding random value from the field
    int countXDown = 0, countXUp = 0, countYDown = 0, countYUp = 0;
    double randomVariable = 0.;
    //Determine x count
    double helpX = 0., helpY = 0.;
    int i = 0, exitFlag = 0;

    //Check if gp lies in random field
    double randomXStart = field(0);
    double randomXEnd = field( 3 * ( numberReal(0) - 1 ) * numberReal(1) );
    double randomYStart = field(1);
    double randomYEnd = field(3 * ( numberReal(0) - 1 ) + 1);

    if ( !( globalCoordinates.at(1) > randomXStart && globalCoordinates.at(1) < randomXEnd &&
            globalCoordinates.at(2) > randomYStart && globalCoordinates.at(2) < randomYEnd ) ) {
        randomVariable = 1;
    } else {
        //Determine the corner points of the square over which we are going to interpolate
        while ( exitFlag == 0 ) {
            if ( field( 3 * i * numberReal(1) ) > globalCoordinates.at(1) ) {
                if ( i == 0 ) { //Check soundness
                    OOFEM_ERROR("i is zero");
                } else if ( i == numberReal(0) ) {
                    OOFEM_ERROR("i is equal to realNumber");
                }

                exitFlag = 1;
                countXDown = i - 1;
                countXUp = i;
                helpX = ( globalCoordinates.at(1) - field( 3 * countXDown * numberReal(1) ) )
                        / ( field( 3 * countXUp * numberReal(1) ) - field( 3 * countXDown * numberReal(1) ) );
            }

            i++;
        }

        //Determine y count
        i = 0, exitFlag = 0;
        while ( exitFlag == 0 ) {
            if ( field(3 * i + 1) > globalCoordinates.at(2) ) {
                if ( i == 0 ) {
                    OOFEM_ERROR("i is zero");
                } else if ( i == numberReal(0) ) {
                    OOFEM_ERROR("i is equal to realNumber");
                }

                exitFlag = 1;
                countYDown = i - 1;
                countYUp = i;
                helpY = ( globalCoordinates.at(2) - field(3 * countYDown + 1) )
                        / ( field(3 * countYUp + 1) - field(3 * countYDown + 1) );
            }

            i++;
        }

        //Do the interpolation
        if ( randomVariable == 0. ) {
            randomVariable =
                ( 1. - helpX ) * ( 1. - helpY ) * field(3 * ( countXDown * numberReal(1) + countYDown ) + 2) +
            helpX * ( 1. - helpY ) * field(3 * ( countXUp * numberReal(1) + countYDown ) + 2) +
            helpX *helpY *field(3 * ( countXUp *numberReal ( 1 ) + countYUp ) + 2) +
            ( 1. - helpX ) * helpY * field(3 * ( countXDown * numberReal(1) + countYUp ) + 2);
        }
    }

    if ( randomVariable <= 0. ) {
        randomVariable = 1.e-8;
    }

    answer = FloatArray{randomVariable};
}

double
InterpolatingFuction :: evaluateAtTime(double t)
{
    OOFEM_ERROR("InterpolatingFunction needs coordinates to evaluate.");
    return 0.;
}



IRResultType
InterpolatingFuction :: initializeFrom(InputRecord *ir)
{
    std :: string name;
    IRResultType result;              // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, name, _IFT_InterpolatingFuction_filename);

    std :: ifstream inputField( name.c_str() );

    if ( !inputField.is_open() ) {
        OOFEM_WARNING("Unable to open file %s", name.c_str());
        return IRRT_BAD_FORMAT;
    }

    double deltaX, deltaY;
    numberReal.resize(2);
    numberReal.zero();
    inputField >> numberReal(0) >> numberReal(1) >> deltaX >> deltaY;

    field.resize( 3 * numberReal(0) * numberReal(1) );
    field.zero();

    //Read in coordinates and field values
    for ( int i = 0; i < numberReal(0) * numberReal(1); i++ ) {
        inputField >> field(3 * i) >> field(3 * i + 1) >> field(3 * i + 2);
    }

    return IRRT_OK;
}
} // end namespace oofem
