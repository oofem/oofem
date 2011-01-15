/* $Header: /home/cvs/bp/oofem/oofemlib/src/initmodulemanager.h,v 1.4.4.1 2004/04/05 15:19:43 bp Exp $ */
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
//   *** CLASS InitModuleManager ***
//   *********************************


#ifndef initmodulemanager_h
#define initmodulemanager_h

#include "alist.h"
#include "modulemanager.h"
#include "initmodule.h"
#include "datareader.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <time.h>
#endif

namespace oofem {
class EngngModel;

/**
 * Class representing and implementing InitModuleManager. It is attribute of EngngModel.
 * It manages the init modules, which perform module - specific init oprations.
 */
class InitModuleManager : public ModuleManager< InitModule >
{
private:
public:
    InitModuleManager(EngngModel *emodel);
    ~InitModuleManager();

    /**
     * Instanciates the receiver from input record. Called from instanciateYourself to initialize yourself
     * from corresponding record. Should be caled before instanciateYourself.
     */
    IRResultType initializeFrom(InputRecord *ir);

    /** Creates new instance of module of given name, with number n, belonging to given EngngModel */
    InitModule *CreateModuleOfType(char *name, int n, EngngModel *emodel);

    /**
     * Performs the initialization of individual modules. Loops over all modules and calls corresponding doInit module service.
     * @param tStep time step.
     */
    void              doInit();
    const char *giveClassName() const { return "InitModuleManager"; }

protected:
};
} // end namespace oofem
#endif // initmodulemanager_h
