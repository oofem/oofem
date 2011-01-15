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

// Initialization module reading data related to Gauss points from a specified file

//
// class gpInitModule
//

#ifndef gpinitmodule_h

 #ifndef __MAKEDEPEND
  #include <stdio.h>
 #endif
 #include "initmodule.h"
 #include "domain.h"
 #include "engngm.h"

namespace oofem {
/**
 * Represents GP (Gauss point) initialization module.
 * This module reads certain internal variables
 * of all Gauss points from a file.
 * In this way, one can specify e.g. initial damage
 * and initial stresses computed by another model.
 */
class GPInitModule : public InitModule
{
public:
    /// Constructor. Creates empty GPInitModule.
    GPInitModule(int n, EngngModel *e);

    /// Destructor
    ~GPInitModule();

    /// Initializes receiver acording to object description stored in input record
    virtual IRResultType initializeFrom(InputRecord *ir);
    /**
     * Reads the specific data (initial state of individual Gauss points)
     * from a file. Each Gauss point is represented by one line in that file.
     * The line should contain the following data:
     * element number
     * Gauss point number
     * coordinates (array, first number of coordinates, then individual coordinates)
     * number of (groups of) variables that will follow
     * for each group:
     *   variable identification number (according to internalstatetype.h)
     *   variable values (array, first number of variables, then their values)
     */
    void doInit();
    /// Returns class name of the receiver.
    virtual const char *giveClassName() const { return "GPInitModuleClass"; }
};

 #define gpinitmodule_h
#endif
} // namespace oofem

