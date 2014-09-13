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

#include "tf1.h"
#include "timestep.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_BoundaryCondition(TF1);

void
TF1 :: computeValueAt(FloatArray &answer, TimeStep *tStep, const FloatArray &coords, ValueModeType mode)
// Returns the value of the receiver at time and given position respecting the mode.
{
    FloatArray cd(3);
    double t, k, c, h, result;
    answer.resize(1);

    k = 2.0e-4;
    c = 1.0;
    h = 200.0;


    t = tStep->giveTargetTime();
    for ( int i = 1; i <= coords.giveSize(); i++ ) {
        cd.at(i) = coords.at(i);
    }

    result = -1.e-5 - k *macbra( atan ( c *t + cd.at ( 2 ) - h ) );

    if ( ( mode == VM_Incremental ) && ( !tStep->isTheFirstStep() ) ) {
        t = tStep->giveTargetTime() - tStep->giveTimeIncrement();
        result -= -1.e-5 - k *macbra( atan ( c *t + cd.at ( 2 ) - h ) );
    }

    answer.at(1) = result;
}

IRResultType
TF1 :: initializeFrom(InputRecord *ir)
{
    return IRRT_OK;
}
} // end namespace oofem
