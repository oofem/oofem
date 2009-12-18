/* $Header: /home/cvs/bp/oofem/sm/src/tf1.C,v 1.4 2003/04/06 14:08:31 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

#include "tf1.h"
#include "timestep.h"
#include "mathfem.h"

namespace oofem {

void
TF1 :: computeValueAt(FloatArray &answer, TimeStep *stepN, FloatArray &coords, ValueModeType mode)
// Returns the value of the receiver at time and given position respecting the mode.
{
    int i;
    FloatArray cd(3);
    double t, k, c, h, result;
    answer.resize(1);

    k = 2.0e-4;
    c = 1.0;
    h = 200.0;


    t = stepN->giveTime();
    for ( i = 1; i <= coords.giveSize(); i++ ) {
        cd.at(i) = coords.at(i);
    }

    result = -1.e-5 - k *macbra( atan(c * t + cd.at(2) - h) );

    if ( ( mode == VM_Incremental ) && ( !stepN->isTheFirstStep() ) ) {
        t = stepN->giveTime() - stepN->giveTimeIncrement();
        result -= -1.e-5 - k *macbra( atan(c * t + cd.at(2) - h) );
    }

    answer.at(1) = result;

    return;
}

IRResultType
TF1 :: initializeFrom(InputRecord *ir)
{
    return IRRT_OK;
}

} // end namespace oofem
