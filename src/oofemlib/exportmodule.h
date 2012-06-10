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

#ifndef exportmodule_h
#define exportmodule_h

#include "intarray.h"
#include "dynalist.h"
#include "inputrecord.h"

namespace oofem {
class EngngModel;
class TimeStep;
class Range;

/**
 * Represents export output module - a base class for all output modules. ExportModule is an abstraction
 * for module performing some specific kind of output. The modules can declare necessary component
 * services using the interface concept. The basic class declares the basic services (the general
 * interface) and implements the services intended to filter output to certain time steps.
 * The output modules are maintained by ExportModuleManager.
 * The output for given time step is done only if this step is selected by one of above
 * described method.
 */
class ExportModule
{
protected:
    /// Component number.
    int number;
    /// Problem pointer.
    EngngModel *emodel;
    /// Indicates all steps selection.
    bool tstep_all_out_flag;
    /// User timeStep Output step. Indicates every tstep_step_out-th step selected.
    int tstep_step_out;
    /// List of user selected step numbers.
    dynaList< Range >tsteps_out;

    /// Indicates all domains.
    bool domain_all_flag;
    /// Domain selection mask.
    IntArray domainMask;

public:

    /// Constructor. Creates empty Output Manager with number n.
    ExportModule(int n, EngngModel *e);
    /// Destructor
    virtual ~ExportModule();
    /// Initializes receiver according to object description stored in input record.
    virtual IRResultType initializeFrom(InputRecord *ir);
    /**
     * Writes the output. Abstract service.
     * @param tStep time step.
     */
    virtual void doOutput(TimeStep *tStep) = 0;
    /**
     * Initializes receiver.
     * The init file messages should be printed.
     */
    virtual void initialize() { }
    /**
     * Terminates the receiver.
     * The terminating messages should be printed.
     * All the streams should be closed.
     */
    virtual void terminate() { }
    /// Returns class name of the receiver.
    virtual const char *giveClassName() const { return "ExportModule"; }

protected:
    /**
     * Gives the appropriate name (minus specific file extension).
     * @param tStep Active time step.
     */
    std::string giveOutputBaseFileName(TimeStep *tStep);
    /**
     * Tests if given time step output is required.
     * @param tStep Time step to check.
     * @return True if output required.
     */
    bool testTimeStepOutput(TimeStep *tStep);
    /**
     * Test if domain output is required.
     * @return True if required.
     */
    bool testDomainOutput(int n);

    /// Prints simple error message and exits.
    void error(const char *file, int line, const char *format, ...) const;
};
} // end namespace oofem
#endif // exportmodule_h



