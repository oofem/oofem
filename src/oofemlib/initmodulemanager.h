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

#ifndef initmodulemanager_h
#define initmodulemanager_h

#include "modulemanager.h"
#include "initmodule.h"
#include "datareader.h"

namespace oofem {
class EngngModel;

/**
 * Class representing and implementing InitModuleManager. It is attribute of EngngModel.
 * It manages the init modules, which perform module - specific init oprations.
 */
class OOFEM_EXPORT InitModuleManager : public ModuleManager< InitModule >
{
public:
    InitModuleManager(EngngModel * emodel);
    virtual ~InitModuleManager();

    InitModule *CreateModule(const char *name, int n, EngngModel *emodel);

    /**
     * Performs the initialization of individual modules.
     * Loops over all modules and calls corresponding doInit module service.
     */
    void doInit();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual const char *giveClassName() const { return "InitModuleManager"; }
};
} // end namespace oofem
#endif // initmodulemanager_h
