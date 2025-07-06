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

#ifndef monitormanager_h
#define monitormanager_h

#include "modulemanager.h"
#include "monitor.h"

namespace oofem {
class EngngModel;

/**
 * Class representing and implementing MonitorManager. It is attribute of EngngModel.
 * It manages individual monitors, which perform monitor - specific output operations.
 * Monitors can be called at user specific events, such as time step termination, etc. 
 * The event type is preopagated to each monitor.
 */
class OOFEM_EXPORT MonitorManager : public ModuleManager< Monitor >
{
public:
    MonitorManager(EngngModel * emodel);
    virtual ~MonitorManager();

    void initializeFrom(InputRecord &ir) override;
    std::unique_ptr<Monitor> CreateModule(const char *name, int num, EngngModel *emodel) override;

    /**
     * Writes the output. Loops over all modules and calls corresponding doOutput module service.
     * @param tStep Time step.
     * @param substepFlag is set to true, only the modules with substepFlag set to true will be processed.
     */
    void update(TimeStep *tStep, Monitor::MonitorEvent event);
    /**
     * Initializes output manager. The corresponding initialize module services are called.
     */
    //void initialize();
    /**
     * Terminates the receiver, the corresponding terminate module services are called.
     */
    //void terminate();
    const char *giveClassName() const override { return "MonitorManager"; }
};
} // end namespace oofem
#endif // monitormanager_h
