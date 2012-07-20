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

#include "rotatingboundary.h"
#include "engngm.h"
#include "dofmanager.h"
#include "mathfem.h"
#include "loadtime.h"

namespace oofem {
double RotatingBoundary :: give(Dof *dof, ValueModeType mode, TimeStep *stepN)
// Returns the value at stepN of the prescribed value of the kinematic
// unknown 'u'. Returns 0 if 'u' has no prescribed value.
{
    DofIDItem id = dof->giveDofID();
    FloatArray *coords = dof->giveDofManager()->giveCoordinates();
    FloatArray answer, newcoords;
    double theta;

    theta = this->giveLoadTimeFunction()->evaluate(stepN, mode);

    if ( axis.giveSize() != 3 ) {
        OOFEM_ERROR("RotatingBoundary :: give - Size of rotation axis != 3.");
    }

    if ( center.giveSize() == 0 ) {
        center.resize( coords->giveSize() );
        center.zero();
    }

    if ( coords == NULL || coords->giveSize() != center.giveSize() ) {
        OOFEM_ERROR("RotatingBoundary :: give - Size of coordinate system different from center of rotation.");
    }

    double &nx = axis.at(1);
    double &ny = axis.at(2);
    double &nz = axis.at(3);

    if ( coords->giveSize() == 1 ) {
        R.resize(1, 1, 0);
        R.at(1, 1) = cos( theta ) + nx * nx * ( 1 - cos( theta ) );
    }
    if ( coords->giveSize() == 2 ) {
        R.resize(2, 2, 0);
        R.at(1, 1) = cos( theta ) + nx * nx * ( 1 - cos( theta ) );
        R.at(1, 2) = nx * ny * ( 1 - cos( theta ) ) - nz * sin( theta );
        R.at(2, 1) = ny * nx * ( 1 - cos( theta ) ) + nz * sin( theta );
        R.at(2, 2) = cos( theta ) + ny * ny * ( 1 - cos( theta ) );
    } else if ( coords->giveSize() == 3  ) {
        R.resize(3, 3, 0);

        R.at(1, 1) = cos( theta ) + nx * nx * ( 1 - cos( theta ) );
        R.at(1, 2) = nx * ny * ( 1 - cos( theta ) ) - nz * sin( theta );
        R.at(1, 3) = nx * nz * ( 1 - cos( theta ) ) + ny * sin( theta );

        R.at(2, 1) = ny * nx * ( 1 - cos( theta ) ) + nz * sin( theta );
        R.at(2, 2) = cos( theta ) + ny * ny * ( 1 - cos( theta ) );
        R.at(2, 3) = ny * nz * ( 1 - cos( theta ) ) - nx * sin( theta );

        R.at(3, 1) = nz * nx * ( 1 - cos( theta ) ) - ny * sin( theta );
        R.at(3, 2) = nz * ny * ( 1 - cos( theta ) ) + nx * sin( theta );
        R.at(3, 3) = cos( theta ) + nz * nz * ( 1 - cos( theta ) );
    } else {
        OOFEM_ERROR("RotatingBoundary :: give - Size of coordinate system has to be 1, 2 or 3.");
    }

    newcoords.subtract( center );
    newcoords.add( *coords );
    answer.beProductOf(R, newcoords );
    answer.add( center );
    answer.subtract( *coords );

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
    return 0.0;
}

IRResultType
RotatingBoundary :: initializeFrom(InputRecord *ir)
// Sets up the dictionary where the receiver stores the conditions it
// imposes.
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    GeneralBoundaryCondition :: initializeFrom(ir);

    isImposedTimeFunction = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, isImposedTimeFunction, IFT_BoundaryCondition_IsImposedTimeFunct, "isimposedtimefunction"); // Macro

    IR_GIVE_FIELD(ir, axis, IFT_RotatingBoundary_axis, "axis"); // Macro
    axis.normalize();

    IR_GIVE_OPTIONAL_FIELD(ir, center, IFT_RotatingBoundary_center, "center"); // Macro

    return IRRT_OK;
}

int
RotatingBoundary :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    GeneralBoundaryCondition :: giveInputRecordString(str, keyword);

    sprintf(buff, " isimposedtimefunction %d ", this->isImposedTimeFunction);
    str += buff;

    return 1;
}
} // end namespace oofem
