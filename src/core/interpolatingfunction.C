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
InterpolatingFuction :: evaluate(FloatArray &answer, const std :: map< std :: string, FunctionArgument > &valDict, GaussPoint *gp, double param)
{
    double randomVariable = 0.;

    auto it = valDict.find("x");
    if ( it == valDict.end() ) {
        OOFEM_ERROR("Coordinate needed for evaluating function");
    }
    const FloatArray &globalCoordinates = it->second.val1;

    if(this->dimension == 2){//2D Case
    // Map the corresponding random value from the field
    int countXDown = 0, countXUp = 0, countYDown = 0, countYUp = 0;
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
    }
    else if(this->dimension == 3){//3D case
    // Map the corresponding random value from the field
    int countXDown = 0, countXUp = 0, countYDown = 0, countYUp = 0, countZDown = 0, countZUp = 0;
    //Determine x count
    double helpX = 0., helpY = 0., helpZ = 0.;
    int i = 0, exitFlag = 0;

    //Check if gp lies in random field
    double randomXStart = field.at(1);
    double randomXEnd = field(4 * ( numberReal.at(1)*numberReal.at(2)*numberReal.at(3) - 1 ));
    double randomYStart = field.at(2);
    double randomYEnd = field(4 * ( numberReal.at(1)*numberReal.at(2)*numberReal.at(3) - 1 ) + 1);
    double randomZStart = field.at(3);
    double randomZEnd = field(4 * ( numberReal.at(1)*numberReal.at(2)*numberReal.at(3) - 1 ) + 2);

    double fieldX0Y0Z0, fieldX1Y0Z0, fieldX0Y1Z0, fieldX0Y0Z1, fieldX1Y1Z0, fieldX1Y1Z1, fieldX0Y1Z1, fieldX1Y0Z1;
    double field00, field01, field10, field11;
    double field0,field1;

    if ( !( globalCoordinates.at(1) >= randomXStart && globalCoordinates.at(1) <= randomXEnd &&
            globalCoordinates.at(2) >= randomYStart && globalCoordinates.at(2) <= randomYEnd &&
	    globalCoordinates.at(3) >= randomZStart && globalCoordinates.at(3) <= randomZEnd ) ) {
        randomVariable = 1;
	printf("Warning: External field too small\n");
    } else {
        //Determine the corner points of the square over which we are going to interpolate
      //Determine xCount
        while ( exitFlag == 0 ) {
	  if ( field.at( 4 * i * numberReal.at(2)*numberReal.at(3) +1) > globalCoordinates.at(1) ) {
                if ( i == 0 ) { //Check soundness
                    OOFEM_ERROR("error in externalfieldgenerator.C: i is zero");
                } else if ( i > numberReal.at(1) ) {
		   OOFEM_ERROR("error in externalfieldgenerator.C: i is greater than realNumber1");
                }

                exitFlag = 1;
                countXDown = i - 1;
                countXUp = i;
                helpX = ( globalCoordinates.at(1) - field.at( 4 * countXDown * numberReal.at(2) * numberReal.at(3)+1 ) )
		  / ( field.at( 4 * countXUp * numberReal.at(2) * numberReal.at(3)+1 ) - field.at( 4 * countXDown * numberReal.at(2) * numberReal.at(3) +1) );
            }

            i++;
        }

        //Determine y count
        i = 0, exitFlag = 0;
        while ( exitFlag == 0 ) {
	  if ( field.at(4 *i*numberReal.at(3) + 2) > globalCoordinates.at(2) ) {
                if ( i == 0 ) {
                    OOFEM_ERROR("error in externalfieldgenerator.C: i is zero");
                } else if ( i > numberReal.at(2) ) {//Not sure that this is a problem
                    OOFEM_ERROR("error in externalfieldgenerator.C: i is greater than realNumber2");
                }

                exitFlag = 1;
                countYDown = i - 1;
                countYUp = i;
                helpY = ( globalCoordinates.at(2) - field.at(4 * countYDown*numberReal.at(3) + 2) )
		  / ( field.at(4 * countYUp*numberReal.at(3) + 2) - field.at(4 * countYDown*numberReal.at(3) + 2) );
            }

            i++;
        }


        //Determine z count
        i = 0, exitFlag = 0;
        while ( exitFlag == 0 ) {
            if ( field.at(4 * i + 3) > globalCoordinates.at(3) ) {
                if ( i == 0 ) {
                    OOFEM_ERROR("error in externalfieldgenerator.C: i is zero");
                } else if ( i > numberReal.at(3) ) {
                    OOFEM_ERROR("error in externalfieldgenerator.C: i is equal to realNumber");
                }

                exitFlag = 1;
                countZDown = i - 1;
                countZUp = i;
                helpZ = ( globalCoordinates.at(3) - field.at(4 * countZDown + 3) )
                        / ( field.at(4 * countZUp + 3) - field.at(4 * countZDown + 3) );
            }

            i++;
        }


        //Do the interpolation in 3D.
        if ( randomVariable == 0. ) {

	  fieldX0Y0Z0 = field.at(4 * ( countXDown * numberReal.at(2)*numberReal.at(3) + countYDown*numberReal.at(3) + countZDown ) + 4);
	  fieldX1Y0Z0 = field.at(4 * ( countXUp * numberReal.at(2)*numberReal.at(3) + countYDown*numberReal.at(3) + countZDown ) + 4);
	  fieldX0Y1Z0 = field.at(4 * ( countXDown * numberReal.at(2)*numberReal.at(3) + countYUp*numberReal.at(3) + countZDown ) + 4);
	  fieldX0Y0Z1 = field.at(4 * ( countXDown * numberReal.at(2)*numberReal.at(3) + countYDown*numberReal.at(3) + countZUp ) + 4);
	  fieldX1Y1Z0 = field.at(4 * ( countXUp * numberReal.at(2)*numberReal.at(3) + countYUp*numberReal.at(3) + countZDown ) + 4);
	  fieldX1Y1Z1 = field.at(4 * ( countXUp * numberReal.at(2)*numberReal.at(3) + countYUp*numberReal.at(3) + countZUp ) + 4);
	  fieldX0Y1Z1 = field.at(4 * ( countXDown * numberReal.at(2)*numberReal.at(3) + countYUp*numberReal.at(3) + countZUp ) + 4);
	  fieldX1Y0Z1 = field.at(4 * ( countXUp * numberReal.at(2)*numberReal.at(3) + countYDown*numberReal.at(3) + countZUp ) + 4);

	  field00 =  fieldX0Y0Z0 * (1. - helpX) +  fieldX1Y0Z0 * helpX;
	  field10 =  fieldX0Y1Z0 * (1. - helpX) + fieldX1Y1Z0 * helpX;
	  field01 =  fieldX0Y0Z1 * (1. - helpX) + fieldX1Y0Z1 * helpX;
	  field11 =  fieldX0Y1Z1 * (1. - helpX) + fieldX1Y1Z1 * helpX;

	  field0 = field00 * (1.-helpY) + field10*helpY;
	  field1 = field01 * (1.-helpY) + field11*helpY;

	  randomVariable = field0*(1. - helpZ) + field1*helpZ;
        }
    }
    }
    else{
      OOFEM_ERROR("Unknown dimension. Must be 2 or 3.");
    }

    answer = FloatArray{randomVariable};
}

double
InterpolatingFuction :: evaluateAtTime(double t)
{
    OOFEM_ERROR("InterpolatingFunction needs coordinates to evaluate.");
}



void
InterpolatingFuction :: initializeFrom(InputRecord &ir)
{
    std :: string name;

    this->dimension = 2;
    IR_GIVE_OPTIONAL_FIELD(ir, dimension, _IFT_InterpolatingFuction_dim);

    IR_GIVE_FIELD(ir, name, _IFT_InterpolatingFuction_filename);

    std :: ifstream inputField( name.c_str() );

    if ( !inputField.is_open() ) {
        throw ValueInputException(ir, _IFT_InterpolatingFuction_filename, "Unable to open file: " + name);
    }

    if(this->dimension == 2){//2D field
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
    }
    else{//3D field
      double deltaX, deltaY, deltaZ;
      numberReal.resize(3);
      numberReal.zero();
      inputField >> numberReal(0) >> numberReal(1)  >> numberReal(2) >> deltaX >> deltaY >> deltaZ;
      field.resize( 4 * numberReal.at(1) * numberReal.at(2) * numberReal.at(3) );
      field.zero();

      //Read in coordinates and field values
      for ( int i = 0; i < numberReal.at(1) * numberReal.at(2) * numberReal.at(3) ; i++ ) {
        inputField >> field.at(4 * i+1) >> field.at(4 * i + 2) >> field.at(4 * i + 3) >> field.at(4 * i + 4);
      }
    }
}
} // end namespace oofem
