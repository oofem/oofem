/* $Header: /home/cvs/bp/oofem/oofemlib/src/exportmodulemanager.h,v 1.4.4.1 2004/04/05 15:19:43 bp Exp $ */
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


//   *********************************
//   *** CLASS ExportModuleManager ***
//   *********************************


#ifndef exportmodulemanager_h
#define exportmodulemanager_h

#include "alist.h"
#include "modulemanager.h"
#include "exportmodule.h"
#include "datareader.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <time.h>
#endif

namespace oofem {
class EngngModel;

/**
 * Class representing and implementing ExportModuleManager. It is attribute of EngngModel.
 * It manages the export output modules, which perform module - specific output oprations.
 */
class ExportModuleManager : public ModuleManager< ExportModule >
{
private:
public:
    ExportModuleManager(EngngModel *emodel);
    ~ExportModuleManager();

    /**
     * Instanciates the receiver from input record. Called from instanciateYourself to initialize yourself
     * from corresponding record. Should be caled before instanciateYourself.
     */
    IRResultType initializeFrom(InputRecord *ir);

    /** Creates new instance of module of given name, with given number, belonging to given EngngModel */
    ExportModule *CreateModuleOfType(char *name, int num, EngngModel *emodel);

    /**
     * Writes the output. Loops over all modules and calls corresponding doOutput module service.
     * @param tStep time step.
     */
    void              doOutput(TimeStep *tStep);
    /**
     * Initializes output manager. The corresponding initialize module services are called.
     */
    void              initialize();
    /**
     * Terminates the receiver, the corresponding terminate module services are called.
     */
    void              terminate();
    const char *giveClassName() const { return "ExportModuleManager"; }

protected:
};
} // end namespace oofem
#endif // exportmodulemanager_h
