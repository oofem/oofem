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

#include "outputmanager.h"
#include "timestep.h"
#include "engngm.h"
#include "element.h"
#include "dofmanager.h"

namespace oofem {
OutputManager :: OutputManager(Domain *d) : dofman_out(), dofman_except(), element_out(), element_except()
{
    domain = d;
    tstep_all_out_flag = 0;
    tstep_step_out = 0;
    dofman_all_out_flag = element_all_out_flag = 0;
}

IRResultType
OutputManager :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // test if record is of OutputManager
    if ( !ir->hasField(IFT_OutputManager_name, "outputmanager") ) {
        OOFEM_ERROR("OutputManager::instanciateFrom: bad record");
    }

    tstep_all_out_flag  = ir->hasField(IFT_OutputManager_tstepall, "tstep_all");
    tstep_step_out      = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, tstep_step_out, IFT_OutputManager_tstepstep, "tstep_step"); // Macro

    dofman_all_out_flag = ir->hasField(IFT_OutputManager_dofmanall, "dofman_all");
    element_all_out_flag = ir->hasField(IFT_OutputManager_elementall, "element_all");

    tsteps_out.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, tsteps_out, IFT_OutputManager_tstepsout, "tsteps_out"); // Macro
    dofman_out.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, dofman_out, IFT_OutputManager_dofmanoutput, "dofman_output"); // Macro
    dofman_except.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, dofman_except, IFT_OutputManager_dofmanexcept, "dofman_except"); // Macro
    element_out.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, element_out, IFT_OutputManager_elementoutput, "element_output"); // Macro
    element_except.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, element_except, IFT_OutputManager_elementexcept, "element_except"); // Macro

    return IRRT_OK;
}

void
OutputManager :: doDofManOutput(FILE *file, TimeStep *tStep)
{
    int i, ndofman = domain->giveNumberOfDofManagers();

    if ( !testTimeStepOutput(tStep) ) {
        return;
    }

    fprintf(file, "\n\nDofManager output:\n------------------\n");

    if ( dofman_all_out_flag   && dofman_except.isEmpty() ) {
        for ( i = 1; i <= ndofman; i++ ) {
#ifdef __PARALLEL_MODE
            // test for null dof in parallel mode
            if ( domain->giveDofManager(i)->giveParallelMode() == DofManager_null ) {
                continue;
            }

#endif
            domain->giveDofManager(i)->printOutputAt(file, tStep);
        }
    } else {
        for ( i = 1; i <= ndofman; i++ ) {
            if ( _testDofManOutput(i) ) {
                domain->giveDofManager(i)->printOutputAt(file, tStep);
            }
        }
    }

    fprintf(file, "\n\n");
}

void
OutputManager :: doElementOutput(FILE *file, TimeStep *tStep)
{
    int i, nelem = domain->giveNumberOfElements();
    if ( !testTimeStepOutput(tStep) ) {
        return;
    }

    fprintf(file, "\n\nElement output:\n---------------\n");

    if ( element_all_out_flag   && element_except.isEmpty() ) {
        for ( i = 1; i <= nelem; i++ ) {
#ifdef __PARALLEL_MODE
            // test for remote element in parallel mode
            if ( domain->giveElement(i)->giveParallelMode() == Element_remote ) {
                continue;
            }

#endif

            domain->giveElement(i)->printOutputAt(file, tStep);
        }
    } else {
        for ( i = 1; i <= nelem; i++ ) {
            if ( _testElementOutput(i) ) {
                domain->giveElement(i)->printOutputAt(file, tStep);
            }
        }
    }

    fprintf(file, "\n\n");
}

int
OutputManager :: _testDofManOutput(int number)
{
    int selected = 0;

#ifdef __PARALLEL_MODE
    // test for null dof in parallel mode
    if ( domain->giveDofManager(number)->giveParallelMode() == DofManager_null ) {
        return 0;
    }

#endif

    // test all_select flag on
    if ( dofman_all_out_flag ) {
        selected  = 1;
    } else {
        // test for particular dofman selection
        dynaList< Range > :: iterator dofmanOutIter;

        int _label = domain->giveDofManager(number)->giveLabel();
        for ( dofmanOutIter = dofman_out.begin(); dofmanOutIter != dofman_out.end(); ++dofmanOutIter ) {
            if ( ( * dofmanOutIter ).test(_label) ) {
                selected  = 1;
                break;
            }
        }
    }

    // if not selected exit
    if ( !selected ) {
        return 0;
    }

    // if selected check exclude list
    dynaList< Range > *list2 = & ( this->dofman_except );
    dynaList< Range > :: iterator dofmanExceptIter;
    int _label = domain->giveDofManager(number)->giveLabel();

    for ( dofmanExceptIter = list2->begin(); dofmanExceptIter != list2->end(); ++dofmanExceptIter ) {
        // test if excluded
        if ( ( * dofmanExceptIter ).test(_label) ) {
            return 0;
        }
    }

    // dofman not excluded;
    return selected;
}


int
OutputManager :: _testElementOutput(int number)
{
    int selected = 0;

#ifdef __PARALLEL_MODE
    // test for remote element in parallel mode
    if ( domain->giveElement(number)->giveParallelMode() == Element_remote ) {
        return 0;
    }

#endif


    // test all_select flag on
    if ( element_all_out_flag ) {
        selected  = 1;
    } else {
        // test for particular element selection
        dynaList< Range > :: iterator elemOutIter;
        int _label = domain->giveElement(number)->giveLabel();

        for ( elemOutIter = element_out.begin(); elemOutIter != element_out.end(); ++elemOutIter ) {
            if ( ( * elemOutIter ).test(_label) ) {
                selected  = 1;
                break;
            }
        }
    }

    // if not selected exit
    if ( !selected ) {
        return 0;
    }

    // if selected check exclude list
    dynaList< Range > :: iterator elemExceptIter;
    int _label = domain->giveElement(number)->giveLabel();

    for ( elemExceptIter = element_except.begin(); elemExceptIter != element_except.end(); ++elemExceptIter ) {
        // test if excluded
        if ( ( * elemExceptIter ).test(_label) ) {
            return 0;
        }
    }

    // element not excluded;
    return selected;
}

int
OutputManager :: testTimeStepOutput(TimeStep *tStep)
{
    if ( tstep_all_out_flag ) {
        return 1;
    }

    if ( tstep_step_out ) {
        if ( ( ( tStep->giveNumber() - domain->giveEngngModel()->giveNumberOfFirstStep() ) % tstep_step_out ) == 0 ) {
            return 1;
        }
    }

    dynaList< Range > :: iterator tstepsIter;
    for ( tstepsIter = tsteps_out.begin(); tstepsIter != tsteps_out.end(); ++tstepsIter ) {
        // test if INCLUDED
        if ( ( * tstepsIter ).test( tStep->giveNumber() ) ) {
            return 1;
        }
    }

    return 0;
}

int
OutputManager :: testDofManOutput(int num, TimeStep *tStep)
{
    if ( testTimeStepOutput(tStep) ) {
        if ( _testDofManOutput(num) ) {
            return 1;
        }
    }

    return 0;
}

int
OutputManager :: testElementOutput(int num, TimeStep *tStep)
{
    if ( testTimeStepOutput(tStep) ) {
        if ( _testElementOutput(num) ) {
            return 1;
        }
    }

    return 0;
}
} // end namespace oofem
