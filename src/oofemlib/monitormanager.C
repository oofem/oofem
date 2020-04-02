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

#include "monitormanager.h"
#include "modulemanager.h"
#include "monitor.h"
#include "classfactory.h"

namespace oofem {
MonitorManager :: MonitorManager(EngngModel *emodel) : ModuleManager< Monitor >(emodel)
{ }

MonitorManager :: ~MonitorManager()
{ }

void
MonitorManager :: initializeFrom(InputRecord &ir)
{
    this->numberOfModules = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfModules, "nmonitors");
}

std::unique_ptr<Monitor> MonitorManager :: CreateModule(const char *name, int n, EngngModel *emodel)
{
    return classFactory.createMonitor(name, n);
}

void
MonitorManager :: update(TimeStep *tStep, Monitor::MonitorEvent event)
{
    for ( auto &module: moduleList ) {
        module->update(this->emodel, tStep, event);
    }
}

} // end namespace oofem
