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

#ifndef initmodule_h
#define initmodule_h

#include "oofemcfg.h"
#include "inputrecord.h"

#include <cstdio>

///@name Input fields for init module
//@{
#define _IFT_InitModule_initfilename "initfile"
//@}

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
class OOFEM_EXPORT InitModule
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
    InitModule(int n, EngngModel * e);
    /// Destructor
    virtual ~InitModule();
    /// Initializes receiver according to object description stored in input record.
    virtual IRResultType initializeFrom(InputRecord *ir) = 0;
    /// Reads the input. Abstract service.
    virtual void doInit() = 0;
    /// Returns class name of the receiver.
    virtual const char *giveClassName() const { return "InitModule"; }
    /// Error printing helper.
    std :: string errorInfo(const char *func) const { return std :: string(giveClassName()) + func; }
};
} // end namespace oofem
#endif // initmodule_h
