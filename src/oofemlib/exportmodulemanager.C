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

#include "exportmodulemanager.h"
#include "modulemanager.h"
#include "engngm.h"
#include "exportmodule.h"

#include "usrdefsub.h"
#include "util.h"
#include "oofem_limits.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <time.h>
#endif

namespace oofem {
ExportModuleManager :: ExportModuleManager(EngngModel *emodel) : ModuleManager< ExportModule >(emodel)
{}

ExportModuleManager :: ~ExportModuleManager()
{}

IRResultType
ExportModuleManager :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;              // Required by IR_GIVE_FIELD macro

    this->numberOfModules = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfModules, IFT_ExportModuleManager_nmodules, "nmodules"); // Macro
    return IRRT_OK;
}



ExportModule *ExportModuleManager :: CreateModuleOfType(const char *name, int n, EngngModel *emodel) {
    return CreateUsrDefExportModuleOfType(name, n, emodel);
}

void
ExportModuleManager :: doOutput(TimeStep *tStep)
{
    int i;

    for ( i = 1; i <= numberOfModules; i++ ) {
        ( ( ExportModule * ) this->giveModule(i) )->doOutput(tStep);
    }
}

void
ExportModuleManager :: initialize()
{
    int i;

    for ( i = 1; i <= numberOfModules; i++ ) {
        ( ( ExportModule * ) this->giveModule(i) )->initialize();
    }
}


void
ExportModuleManager :: terminate()
{
    int i;

    for ( i = 1; i <= numberOfModules; i++ ) {
        ( ( ExportModule * ) this->giveModule(i) )->terminate();
    }
}
} // end namespace oofem
