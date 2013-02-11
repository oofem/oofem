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

#ifndef exportmodulemanager_h
#define exportmodulemanager_h

#include "modulemanager.h"
#include "exportmodule.h"

namespace oofem {
class EngngModel;

/**
 * Class representing and implementing ExportModuleManager. It is attribute of EngngModel.
 * It manages the export output modules, which perform module - specific output operations.
 */
class ExportModuleManager : public ModuleManager< ExportModule >
{
public:
    ExportModuleManager(EngngModel *emodel);
    virtual ~ExportModuleManager();

    virtual IRResultType initializeFrom(InputRecord *ir);
    ExportModule *CreateModuleOfType(const char *name, int num, EngngModel *emodel);

    /**
     * Writes the output. Loops over all modules and calls corresponding doOutput module service.
     * @param tStep Time step.
     */
    void doOutput(TimeStep *tStep);
    /**
     * Initializes output manager. The corresponding initialize module services are called.
     */
    void initialize();
    /**
     * Terminates the receiver, the corresponding terminate module services are called.
     */
    void terminate();
    virtual const char *giveClassName() const { return "ExportModuleManager"; }
};
} // end namespace oofem
#endif // exportmodulemanager_h
