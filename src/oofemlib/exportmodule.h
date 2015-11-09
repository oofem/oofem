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

#ifndef exportmodule_h
#define exportmodule_h

#include "oofemcfg.h"
#include "intarray.h"
#include "inputrecord.h"
#include "range.h"
#include "set.h"

#include <list>

///@name Input fields for export module
//@{
#define _IFT_ExportModule_tstepall "tstep_all"
#define _IFT_ExportModule_tstepstep "tstep_step"
#define _IFT_ExportModule_tstepsout "tsteps_out"
#define _IFT_ExportModule_subtstepsout "subtsteps_out"
#define _IFT_ExportModule_domainall "domain_all"
#define _IFT_ExportModule_domainmask "domain_mask"
#define _IFT_ExportModule_regionsets "regionsets"
#define _IFT_ExportModule_timescale "timescale"
//@}

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
class OOFEM_EXPORT ExportModule
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
    std :: list< Range >tsteps_out;
    /**
     * Flag turning output in solution step substeps/itarations. Allows to visualize the
     * varibles during equilibrium iterations, etc. Usefull for debugging.
     * In this case the export module output name will contain timeStep substep number.
     *  This tStep substep must be set/managed by corresponding engineering model.
     */
    bool tstep_substeps_out_flag;

    /// Indicates all domains.
    bool domain_all_flag;
    /// Domain selection mask.
    IntArray domainMask;

    /// regions represented by sets
    IntArray regionSets;

    /// Scaling time in output, e.g. conversion from seconds to hours
    double timeScale;

    /// Returns number of regions (aka regionSets)
    int giveNumberOfRegions();

    /// Default region set
    Set defaultElementSet;

    /// Returns element set
    Set *giveRegionSet(int i);

public:

    /// Constructor. Creates empty Output Manager with number n.
    ExportModule(int n, EngngModel * e);
    /// Destructor
    virtual ~ExportModule();
    /// Initializes receiver according to object description stored in input record.
    virtual IRResultType initializeFrom(InputRecord *ir);
    /**
     * Writes the output. Abstract service.
     * @param tStep Time step.
     * @param forcedOutput If true, no testTimeStepOutput should be done.
     */
    virtual void doOutput(TimeStep *tStep, bool forcedOutput = false) = 0;
    /**
     * Writes the output. Abstract service.
     * @param tStep time step.
     */
    void doForcedOutput(TimeStep *tStep) { doOutput(tStep, true); }
    /**
     * Initializes receiver.
     * The init file messages should be printed.
     */
    /**
     * Returns true if module is configured to export indvidual substep/iterations.
     */
    bool testSubStepOutput() { return this->tstep_substeps_out_flag; }

    virtual void initialize();

    /**
     * Fill regionSets with all elements if regionSets is initially empty.
     * Implementation depends on derived classes.
    */
    virtual void initializeElementSet();

    /**
     * Terminates the receiver.
     * The terminating messages should be printed.
     * All the streams should be closed.
     */
    virtual void terminate() { }
    /// Returns class name of the receiver.
    virtual const char *giveClassName() const = 0;

protected:
    /**
     * Gives the appropriate name (minus specific file extension).
     * @param tStep Active time step.
     */
    std :: string giveOutputBaseFileName(TimeStep *tStep);
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

    /// Returns string for prepending output (used by error reporting macros).
    std :: string errorInfo(const char *func) const;
};
} // end namespace oofem
#endif // exportmodule_h
