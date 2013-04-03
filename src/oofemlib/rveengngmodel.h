/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
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


#ifndef rveengngmodel_h
#define rveengngmodel_h

#include "floatarray.h"
#include "timestep.h"

namespace oofem {

/**
 * The rveEngngModel is an abstract class introducing functions for setting boundary conditions on an RVE and computing characteristic data.
 * This class is used when creating an engineering model with RVE capabilities
 *
 * @author Carl Sandström
 */
class rveEngngModel
{
public:
    /**
     * Abstract method for setting boundary condition on RVE.
     * @param type Specifies the type of boundary condition. Type are defines in classes which inherit rveEngngModel
     * @param value Includes data for the boundary condition.
     */
    virtual void rveSetBoundaryConditions(int type, FloatArray value) = 0;
    /**
     * Abstract method for computing data on the RVE.
     * @param type Specifies the data requested. The type is defined in classes which inherit rveEngngModel
     * @param value The value pertinent to the boundary condition
     * @param answer The response produced by the specified boundary conditions
     * @param atTime Pertinent timestep
     */
    virtual void rveGiveCharacteristicData(int type, void *value, void *answer, TimeStep *atTime) = 0;
};

}

#endif // rveengngmodel_h
