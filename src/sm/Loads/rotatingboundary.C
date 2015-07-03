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

#include "rotatingboundary.h"
#include "dofmanager.h"
#include "dof.h"
#include "mathfem.h"
#include "function.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"

namespace oofem {
REGISTER_BoundaryCondition(RotatingBoundary);

double RotatingBoundary :: give(Dof *dof, ValueModeType mode, double time)
{
    DofIDItem id = dof->giveDofID();
    FloatArray *coords = dof->giveDofManager()->giveCoordinates();
    FloatArray answer, newcoords;
    double theta = 0.;

    if ( mode == VM_Total ) {
        theta = this->giveTimeFunction()->evaluateAtTime(time);
    } else if ( mode == VM_Velocity ) {
        theta = this->giveTimeFunction()->evaluateVelocityAtTime(time);
    } else if ( mode == VM_Acceleration ) {
        theta = this->giveTimeFunction()->evaluateAccelerationAtTime(time);
    } else {
        OOFEM_ERROR("Should not be called for value mode type then total, velocity, or acceleration.");
    }

    if ( axis.giveSize() != 3 ) {
        OOFEM_ERROR("Size of rotation axis != 3.");
    }

    if ( center.giveSize() == 0 ) {
        center.resize( coords->giveSize() );
        center.zero();
    }

    if ( coords == NULL || coords->giveSize() != center.giveSize() ) {
        OOFEM_ERROR("Size of coordinate system different from center of rotation.");
    }

    double &nx = axis.at(1);
    double &ny = axis.at(2);
    double &nz = axis.at(3);

    if ( coords->giveSize() == 1 ) {
        R.resize(1, 1);
        R.at(1, 1) = cos(theta) + nx * nx * ( 1 - cos(theta) );
    }
    if ( coords->giveSize() == 2 ) {
        R.resize(2, 2);
        R.at(1, 1) = cos(theta) + nx * nx * ( 1 - cos(theta) );
        R.at(1, 2) = nx * ny * ( 1 - cos(theta) ) - nz *sin(theta);
        R.at(2, 1) = ny * nx * ( 1 - cos(theta) ) + nz *sin(theta);
        R.at(2, 2) = cos(theta) + ny * ny * ( 1 - cos(theta) );
    } else if ( coords->giveSize() == 3  ) {
        R.resize(3, 3);

        R.at(1, 1) = cos(theta) + nx * nx * ( 1 - cos(theta) );
        R.at(1, 2) = nx * ny * ( 1 - cos(theta) ) - nz *sin(theta);
        R.at(1, 3) = nx * nz * ( 1 - cos(theta) ) + ny *sin(theta);

        R.at(2, 1) = ny * nx * ( 1 - cos(theta) ) + nz *sin(theta);
        R.at(2, 2) = cos(theta) + ny * ny * ( 1 - cos(theta) );
        R.at(2, 3) = ny * nz * ( 1 - cos(theta) ) - nx *sin(theta);

        R.at(3, 1) = nz * nx * ( 1 - cos(theta) ) - ny *sin(theta);
        R.at(3, 2) = nz * ny * ( 1 - cos(theta) ) + nx *sin(theta);
        R.at(3, 3) = cos(theta) + nz * nz * ( 1 - cos(theta) );
    } else {
        OOFEM_ERROR("Size of coordinate system has to be 1, 2 or 3.");
    }

    newcoords.beDifferenceOf(center, *coords);
    answer.beProductOf(R, newcoords);
    answer.add(center);
    answer.subtract(* coords);

    switch ( id ) {
    case D_u:
        return answer.at(1);

    case D_v:
        return answer.at(2);

    case D_w:
        return answer.at(3);

    default:
        return 0.0;
    }
}

IRResultType
RotatingBoundary :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, axis, _IFT_RotatingBoundary_axis);
    axis.normalize();

    IR_GIVE_OPTIONAL_FIELD(ir, center, _IFT_RotatingBoundary_center);

    return GeneralBoundaryCondition :: initializeFrom(ir);
}

void
RotatingBoundary :: giveInputRecord(DynamicInputRecord &input)
{
    GeneralBoundaryCondition :: giveInputRecord(input);
    input.setField(this->axis, _IFT_RotatingBoundary_axis);
    input.setField(this->center, _IFT_RotatingBoundary_center);
}
} // end namespace oofem
