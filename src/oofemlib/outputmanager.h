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

#ifndef outputmanager_h
#define outputmanager_h

#include "domain.h"
#include "engngm.h"
#include "intarray.h"
#include "dynalist.h"
#include "range.h"

namespace oofem {
/**
 * Represents output manager. It controls and manages the time step output.
 * Allows to filter output to certain time steps, dof managers and elements.
 * For particular component (time step, dofmanager or element) it is possible
 * to filter its output using following modes:
 * - All component can be selected.
 * - Particular component can be selected.
 *   Particular time steps can be selected by specifying their list or by specifying outputStep - the each n-th
 *   step will be selected starting from firstStep (see engngModel).
 *   Particular dofmanagers and elements can be selected by specifying their list. The member of list can be
 *   component number or range of component numbers.
 * - Particular components can be excluded from selection.
 *
 * The output will be done for particular component if it is selected and is not excluded.
 * The output for given time step is done only if this step is selected by one of above
 * described method. The output for dofmanagers and elements in given step is done only if
 * particular time step is selected and if particular dofmanager or element is selected.
 */
class OutputManager
{
protected:
    /// Domain pointer.
    Domain *domain;
    /// Indicates all steps selection.
    int tstep_all_out_flag;
    /// User timeStep Output step. Indicates every tstep_step_out-th step selected.
    int tstep_step_out;
    /// List of user selected step numbers.
    dynaList< Range >tsteps_out;

    /// Indicates all dofmanagers are selected.
    int dofman_all_out_flag;
    /// List of dofmanager numbers or their ranges being selected.
    dynaList< Range >dofman_out;
    /// List of dofmanager numbers or their ranges being excluded.
    dynaList< Range >dofman_except;

    /// Indicates all elements are selected.
    int element_all_out_flag;
    /// List of element numbers or their ranges being selected.
    dynaList< Range >element_out;
    /// List of element numbers or their ranges being excluded.
    dynaList< Range >element_except;

public:

    /// Creates empty Output Manager. By default all components are selected.
    OutputManager(Domain *d);
    /// Initializes receiver according to object description stored in input record.
    IRResultType initializeFrom(InputRecord *ir);

    /**
     * Does the dofmanager output.
     * All selected dofmanagers are requested for doing their output using printOutputAt service.
     */
    void doDofManOutput(FILE *, TimeStep *);
    /**
     * Does the element output.
     * All selected elements are requested for doing their output using printOutputAt service.
     */
    void doElementOutput(FILE *, TimeStep *);

    /**
     * Tests if given dof manager is required to do its output for given time step.
     * @return nonzero if output required.
     */
    int testDofManOutput(int, TimeStep *);
    /**
     * Tests if given element is required to do its output for given time step.
     * @return nonzero if output required.
     */
    int testElementOutput(int, TimeStep *);
    /**
     * Tests if given dof manager is required to do its output.
     * The time step is not considered.
     * @return nonzero if output required.
     */
    int _testDofManOutput(int number);
    /**
     * Tests if given element  is required to do its output.
     * The time step is not considered.
     * @return nonzero if output required.
     */
    int _testElementOutput(int number);
    /**
     * Tests if given time step output is required.
     * @return nonzero if output required.
     */
    int testTimeStepOutput(TimeStep *);
    /// prints simple error message and exits
    const char *giveClassName() const { return "OutputManager"; }

    /**
     * Receiver becomes shallow copy of the argument. Shallow here means that only
     * tstep_all_out_flag, tstep_step_out, dofman_all_out_flag, and element_all_out_flag are copied.
     */
    void beCopyOf(OutputManager *om) {
        this->tstep_all_out_flag = om->tstep_all_out_flag;
        this->tstep_step_out = om->tstep_step_out;
        this->dofman_all_out_flag = om->dofman_all_out_flag;
        this->element_all_out_flag = om->element_all_out_flag;
    }
};
} // end namespace oofem
#endif // outputmanager_h
