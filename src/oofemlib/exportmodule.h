/* $Header: /home/cvs/bp/oofem/oofemlib/src/exportmodule.h,v 1.5.4.1 2004/04/05 15:19:43 bp Exp $ */
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

//
// class ExportModule
//

#ifndef exportmodule_h
#define exportmodule_h

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <string.h>
#endif
#include "intarray.h"
#include "dynalist.h"
#include "range.h"
#include "inputrecord.h"

namespace oofem {
class EngngModel;
class TimeStep;

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

    /// Component number
    int number;
    /// Problem pointer
    EngngModel *emodel;
    /// Indicates all steps selection
    int tstep_all_out_flag;
    /// User timeStep Output step. Indicates every tstep_step_out-th step selected.
    int tstep_step_out;
    /// List of user selected step numbers.
    dynaList< Range >tsteps_out;

    /// Indicates all domains
    int domain_all_flag;
    /// Domain selection mask
    IntArray domainMask;

public:

    /// Constructor. Creates empty Output Manager with number n. By default all components are selected.
    ExportModule(int n, EngngModel *e);
    /// Destructor
    virtual ~ExportModule();
    /// Initializes receiver acording to object description stored in input record.
    virtual IRResultType initializeFrom(InputRecord *ir);
    /**
     * Writes the output. Abstract service.
     * @param tStep time step.
     */
    virtual void                  doOutput(TimeStep *tStep) = 0;
    /**
     * Initializes receiver.
     * The init file messages should be printed.
     */
    virtual void                  initialize() { }
    /**
     * Terminates the receiver.
     * The terminating messages should be printed.
     * All the streams should be closed.
     */
    virtual void                  terminate() { }
    /// Returns class name of the receiver.
    virtual const char *giveClassName() const { return "ExportModule"; }




protected:
    /**
     * Tests if given time step output is required.
     * @return nonzero if output required.
     */
    int testTimeStepOutput(TimeStep *);
    /**
     * Test if domain output is required.
     * @return nonzero if required
     */
    int testDomainOutput(int n);

    /// prints simple error message and exits
    void error(const char *file, int line, const char *format, ...) const;

protected:
};
} // end namespace oofem
#endif // exportmodule_h



