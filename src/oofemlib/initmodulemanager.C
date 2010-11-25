/* $Header: /home/cvs/bp/oofem/oofemlib/src/initmodulemanager.C,v 1.8.4.1 2004/04/05 15:19:43 bp Exp $ */
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

#include "initmodulemanager.h"
#include "modulemanager.h"
#include "engngm.h"
#include "exportmodule.h"
#include "strreader.h"

#include "usrdefsub.h"
#include "util.h"
#include "oofem_limits.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <time.h>
#endif

namespace oofem {
InitModuleManager :: InitModuleManager(EngngModel *emodel) : ModuleManager< InitModule >(emodel)
{}

InitModuleManager :: ~InitModuleManager()
{}

IRResultType
InitModuleManager :: initializeFrom(InputRecord *ir) {
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;              // Required by IR_GIVE_FIELD macro

    this->numberOfModules = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfModules, IFT_InitModuleManager_nmodules, "ninitmodules"); // Macro
    return IRRT_OK;
}

InitModule *InitModuleManager :: CreateModuleOfType(char *name, EngngModel *emodel) {
    return CreateUsrDefInitModuleOfType(name, emodel);
}

void
InitModuleManager :: doInit()
{
    int i;

    for ( i = 1; i <= numberOfModules; i++ ) {
        ( ( InitModule * ) this->giveModule(i) )->doInit();
    }
}
} // end namespace oofem
