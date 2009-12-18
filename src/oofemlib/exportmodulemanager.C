/* $Header: /home/cvs/bp/oofem/oofemlib/src/exportmodulemanager.C,v 1.8.4.1 2004/04/05 15:19:43 bp Exp $ */
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

#include "exportmodulemanager.h"
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

ExportModuleManager :: ExportModuleManager(EngngModel *emodel)
{
    this->emodel = emodel;

    moduleList = new AList< ExportModule >(0);
    numberOfModules = 0;
}

ExportModuleManager :: ~ExportModuleManager()
{
    delete moduleList;
}

ExportModule *
ExportModuleManager :: giveExportModule(int num)
// Returns the n-th module.
{
    ExportModule *elem = NULL;

    if ( moduleList->includes(num) ) {
        elem = moduleList->at(num);
    } else {
        OOFEM_ERROR2("ExportModuleManager::giveOuputModule: No module no. %d defined", num);
    }

    return elem;
}

int
ExportModuleManager :: instanciateYourself(DataReader *dr, InputRecord *ir)
{
    const char *__proc = "instanciateYourself"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                     // Required by IR_GIVE_FIELD macro

    int i;
    char name [ MAX_NAME_LENGTH ];
    ExportModule *module;
    InputRecord *mir;

    //this->initializeFrom (ir);

    // read modules
    moduleList->growTo(numberOfModules);
    for ( i = 0; i < numberOfModules; i++ ) {
        mir = dr->giveInputRecord(DataReader :: IR_expModuleRec, i + 1);
        result = mir->giveRecordKeywordField(name, MAX_NAME_LENGTH);
        if ( result != IRRT_OK ) {
            IR_IOERR(giveClassName(), __proc, IFT_RecordIDField, "", mir, result);
        }

        // read type of module
        module = CreateUsrDefExportModuleOfType(name, emodel);
        if ( module == NULL ) {
            OOFEM_ERROR2("ExportModuleManager::instanciateYourself: unknown module (%s)", name);
        }

        module->initializeFrom(mir);
        moduleList->put(i + 1, module);
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Instanciated output modules ", nnode)
#  endif
    return 1;
}


IRResultType
ExportModuleManager ::  initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    this->numberOfModules = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfModules, IFT_ExportModuleManager_nmodules, "nmodules"); // Macro
    return IRRT_OK;
}


void
ExportModuleManager :: doOutput(TimeStep *tStep)
{
    int i;

    for ( i = 1; i <= numberOfModules; i++ ) {
        this->giveExportModule(i)->doOutput(tStep);
    }
}

void
ExportModuleManager :: initialize()
{
    int i;

    for ( i = 1; i <= numberOfModules; i++ ) {
        this->giveExportModule(i)->initialize();
    }
}


void
ExportModuleManager :: terminate()
{
    int i;

    for ( i = 1; i <= numberOfModules; i++ ) {
        this->giveExportModule(i)->terminate();
    }
}

} // end namespace oofem
