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

#include "linearedgeload.h"
#include "floatarray.h"
#include "mathfem.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_BoundaryCondition(LinearEdgeLoad);

void
LinearEdgeLoad :: initializeFrom(InputRecord &ir)
{
    BoundaryLoad :: initializeFrom(ir);

    int fType = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, fType, _IFT_LinearEdgeLoad_formulation);
    if ( fType == 1 ) {
        this->formulation = FT_Global;
        // read start and end coordinates
        IR_GIVE_FIELD(ir, startCoords, _IFT_LinearEdgeLoad_startcoord);
        IR_GIVE_FIELD(ir, endCoords, _IFT_LinearEdgeLoad_endcoord);
        if ( startCoords.isEmpty() || endCoords.isEmpty() ) {
            throw ValueInputException(ir, "sc/ec", "coordinates not specified");
        }
    } else {
        this->formulation = FT_Entity;
    }
}


void LinearEdgeLoad :: giveInputRecord(DynamicInputRecord &input)
{
    BoundaryLoad :: giveInputRecord(input);
    input.setField(this->formulation, _IFT_LinearEdgeLoad_formulation);
    if ( this->formulation == FT_Global ) {
        input.setField(this->startCoords, _IFT_LinearEdgeLoad_startcoord);
        input.setField(this->endCoords, _IFT_LinearEdgeLoad_endcoord);
    }
}


void
LinearEdgeLoad :: computeNArray(FloatArray &answer, const FloatArray &coords) const
{
    // compute local isoparametric coordinates of given point
    double ksi;

    if ( formulation == FT_Global ) {
        double length = distance(endCoords, startCoords);
        double dl     = distance(coords, startCoords);
        double eta = dl / length;
        ksi    = ( dl - 0.5 * length ) / ( 0.5 * length );
        FloatArray dir = endCoords;

        dir.subtract(startCoords);

        if ( ( ksi < -1.0 ) ||  ( ksi > 1.0 ) ) {
            OOFEM_WARNING("point out of receiver, skipped", 1);
            answer.resize(2);
            answer.zero();
        }

        for ( int i = 1; i <= dir.giveSize(); i++ ) {
            if ( fabs( startCoords.at(i) + dir.at(i) * eta - coords.at(i) ) > 1.e-6 ) {
                OOFEM_WARNING("point out of receiver, skipped", 1);
                answer.resize(2);
                answer.zero();
            }
        }
    } else {
        ksi = coords.at(1);
    }

    answer.resize(2);

    answer.at(1) = ( 1. - ksi ) * 0.5;
    answer.at(2) = ( 1. + ksi ) * 0.5;
}
} // end namespace oofem
