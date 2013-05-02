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

#include "deadweight.h"
#include "timestep.h"
#include "loadtimefunction.h"
#include "classfactory.h"

namespace oofem {

REGISTER_BoundaryCondition( DeadWeight );

void DeadWeight :: computeValueAt(FloatArray& answer, TimeStep* atTime, FloatArray& coords, ValueModeType mode)
{
    computeComponentArrayAt(answer, atTime, mode);
}

void DeadWeight :: setDeadWeighComponents(const FloatArray& newComponents)
{
    this->componentArray.at(1) = newComponents.at(1);
    this->componentArray.at(2) = newComponents.at(2);
}


} // end namespace oofem
