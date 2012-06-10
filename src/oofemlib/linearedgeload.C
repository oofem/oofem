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

#include "linearedgeload.h"
#include "flotarry.h"
#include "mathfem.h"

namespace oofem {
IRResultType
LinearEdgeLoad :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    BoundaryLoad :: initializeFrom(ir);
    if ( componentArray.giveSize() != nDofs * 2 ) {
        _error("instanciateFrom: componentArray size mismatch");
    }

    int fType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, fType, IFT_LinearEdgeLoad_formulation, "formulation"); // Macro
    if ( fType == 1 ) {
        this->formulation = BL_GlobalFormulation;
        // read start and end coordinates
        IR_GIVE_FIELD(ir, startCoords, IFT_LinearEdgeLoad_startcoord, "sc"); // Macro
        IR_GIVE_FIELD(ir, endCoords, IFT_LinearEdgeLoad_endcoord, "ec"); // Macro
        if ( startCoords.isEmpty() || endCoords.isEmpty() ) {
            _error("instanciateFrom: coordinates not specified");
        }
    } else {
        this->formulation = BL_EntityFormulation;
    }

    return IRRT_OK;
}


int
LinearEdgeLoad :: giveInputRecordString(std :: string &str, bool keyword)
{
    int i;
    char buff [ 1024 ];

    BoundaryLoad :: giveInputRecordString(str, keyword);
    sprintf(buff, " formulation %d", ( int ) this->formulation);
    str += buff;
    if ( this->formulation == BL_GlobalFormulation ) {
        sprintf( buff, " sc %d", this->startCoords.giveSize() );
        str += buff;
        for ( i = 1; i <= this->startCoords.giveSize(); i++ ) {
            sprintf( buff, " %e", this->startCoords.at(i) );
            str += buff;
        }

        sprintf( buff, " ec %d", this->endCoords.giveSize() );
        str += buff;
        for ( i = 1; i <= this->endCoords.giveSize(); i++ ) {
            sprintf( buff, " %e", this->endCoords.at(i) );
            str += buff;
        }
    }

    return 1;
}




void
LinearEdgeLoad :: computeNArray(FloatArray &answer, FloatArray &coords) const
{
    // compute local isoparametric coordinates of given point
    double ksi;

    if ( formulation == BL_GlobalFormulation ) {
        int i;
        double length = endCoords.distance(startCoords);
        double dl     = coords.distance(startCoords);
        double eta = dl / length;
        ksi    = ( dl - 0.5 * length ) / ( 0.5 * length );
        FloatArray dir = endCoords;

        dir.subtract(startCoords);

        if ( ( ksi < -1.0 ) ||  ( ksi > 1.0 ) ) {
            _warning2("computeNArray: point out of receiver, skipped", 1);
            answer.resize(2);
            answer.zero();
        }

        for ( i = 1; i <= dir.giveSize(); i++ ) {
            if ( fabs( startCoords.at(i) + dir.at(i) * eta - coords.at(i) ) > 1.e-6 ) {
                _warning2("computeNArray: point out of receiver, skipped", 1);
                answer.resize(2);
                answer.zero();
            }
        }
    } else {
        ksi = coords.at(1);
    }

    double n1, n2;

    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;

    answer.resize(2);

    answer.at(1) = n1;
    answer.at(2) = n2;
}
} // end namespace oofem
