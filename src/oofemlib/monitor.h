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

#ifndef monitor_h
#define monitor_h

#include "oofemcfg.h"

namespace oofem {
class EngngModel;
class TimeStep;
class InputRecord;

///@name Input fields for Monitor
//@{
#define _IFT_DummyMonitor_Name "dummymonitor"


/**
 * Class representing monitor, an abstract class inplementing solution monitor.
 * Monitors are managed by MonitorManager and are invoked at user controlled events, such as solution step termoination.
 * Derived classes can perform specific actions (updating model, plotting proogress, etc) 
 */
class OOFEM_EXPORT Monitor
{
protected:
    int number;
public:
    enum MonitorEvent {
        TimeStepTermination
    };
public:
    Monitor(int n) {
        this->number = n;
    }
    virtual ~Monitor() {}

    /// Initializes receiver according to object description stored in input record.
    virtual void initializeFrom(InputRecord &ir) = 0;
    /**
     * updates the monitor state. This can also mean updating received eModel state.
     * @param eModel instance of EngngModel
     * @param tStep time step
     * @param event event type 
     */
    virtual void update(EngngModel* eModel, TimeStep *tStep, MonitorEvent event) = 0;
 
    /// Returns class name of the receiver.
    virtual const char *giveClassName() const = 0;

};

/** Dummy monitor */
class OOFEM_EXPORT DummyMonitor : public Monitor
{
    public:
      DummyMonitor (int n) : Monitor (n) {}

        /// Initializes receiver according to object description stored in input record.
        void initializeFrom(InputRecord &ir) override;
        void update(EngngModel* eModel, TimeStep *tStep, MonitorEvent event) override;
        virtual const char *giveClassName() const override {return "DummuMonitor"; }
};

} // end namespace oofem
#endif // monitor_h
