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

#include "exportmodulemanager.h"
#include "modulemanager.h"
#include "exportmodule.h"
#include "classfactory.h"

namespace oofem {
ExportModuleManager :: ExportModuleManager(EngngModel *emodel) : ModuleManager< ExportModule >(emodel)
{ }

ExportModuleManager :: ~ExportModuleManager()
{ }

IRResultType
ExportModuleManager :: initializeFrom(InputRecord *ir)
{
    IRResultType result;              // Required by IR_GIVE_FIELD macro

    this->numberOfModules = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfModules, _IFT_ModuleManager_nmodules);
    return IRRT_OK;
}

ExportModule *ExportModuleManager :: CreateModule(const char *name, int n, EngngModel *emodel)
{
    return classFactory.createExportModule(name, n, emodel);
}

void
ExportModuleManager :: doOutput(TimeStep *tStep, bool substepFlag)
{
    for ( auto &module: moduleList ) {
        if ( substepFlag ) {
            if ( module->testSubStepOutput() ) {
                module->doOutput(tStep);
            }
        } else {
            module->doOutput(tStep);
        }
    }
}

void
ExportModuleManager :: initialize()
{
    for ( auto &module: moduleList ) {
        module->initialize();
    }
}


void
ExportModuleManager :: terminate()
{
    for ( auto &module: moduleList ) {
        module->terminate();
    }
}
} // end namespace oofem
