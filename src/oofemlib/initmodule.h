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

#ifndef initmodule_h
#define initmodule_h

#include "inputrecord.h"

namespace oofem {
class EngngModel;
class TimeStep;

/**
 * Represents init module - a base class for all init modules. InitModule is an abstraction
 * for module performing some specific kind of initialization. The modules can declare necessary component
 * services using the interface concept. The basic class declares the basic services (the general
 * interface).
 * The initialization modules are maintained by InitModuleManager.
 * The initialization for is done only once, at simulation startup by one of above
 * described method.
 */
class InitModule
{
protected:
    /// Number.
    int number;
    /// Problem pointer.
    EngngModel *emodel;
    /// Initialization file.
    FILE *initStream;

public:
    /// Constructor. Creates empty Init Module with number n.
    InitModule(int n, EngngModel *e);
    /// Destructor
    virtual ~InitModule();
    /// Initializes receiver according to object description stored in input record.
    virtual IRResultType initializeFrom(InputRecord *ir) = 0;
    /// Reads the input. Abstract service.
    virtual void doInit() = 0;
    /// Returns class name of the receiver.
    virtual const char *giveClassName() const { return "InitModule"; }
};
} // end namespace oofem
#endif // initmodule_h



