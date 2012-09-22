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

#include "constantedgeload.h"
#include "loadtime.h"
#include "flotarry.h"
#include "timestep.h"

namespace oofem {
IRResultType
ConstantEdgeLoad :: initializeFrom(InputRecord *ir)
{
    BoundaryLoad :: initializeFrom(ir);
    if ( componentArray.giveSize() != nDofs ) {
        _error("instanciateFrom: componentArray size mismatch");
    }

    return IRRT_OK;
}

void
ConstantEdgeLoad :: computeValueAt(FloatArray &answer, TimeStep *stepN, FloatArray &coords, ValueModeType mode)
{
    // we overload general implementation on the boundary load level due
    // to implementation efficiency

    double factor;

    if ( ( mode != VM_Total ) && ( mode != VM_Incremental ) ) {
        _error("computeValueAt: mode not supported");
    }

    // ask time distribution
    /*
     * factor = this -> giveLoadTimeFunction() -> at(stepN->giveTime()) ;
     * if ((mode==VM_Incremental) && (!stepN->isTheFirstStep()))
     * //factor -= this->giveLoadTimeFunction()->at(stepN->givePreviousStep()->giveTime()) ;
     * factor -= this->giveLoadTimeFunction()->at(stepN->giveTime()-stepN->giveTimeIncrement()) ;
     */
    factor = this->giveLoadTimeFunction()->evaluate(stepN, mode);
    answer = componentArray;
    answer.times(factor);
    
    if ( !isImposed(stepN) ){
        answer.zero();
    }
}
} // end namespace oofem
