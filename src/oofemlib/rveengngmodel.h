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


#ifndef rveengngmodel_h
#define rveengngmodel_h

#include "oofemcfg.h"
#include "floatarray.h"
#include "timestep.h"

namespace oofem {
/**
 * The rveEngngModel is an abstract class introducing functions for setting boundary conditions on an RVE and computing characteristic data.
 * This class is used when creating an engineering model with RVE capabilities
 *
 * @author Carl Sandström
 */
class OOFEM_EXPORT rveEngngModel
{
public:

    rveEngngModel() { }
    virtual ~rveEngngModel() { }

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
     * @param tStep Pertinent timestep
     */
    virtual void rveGiveCharacteristicData(int type, void *value, void *answer, TimeStep *tStep) = 0;
};
}

#endif // rveengngmodel_h
