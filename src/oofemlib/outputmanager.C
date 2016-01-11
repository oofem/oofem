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

#include "outputmanager.h"
#include "timestep.h"
#include "engngm.h"
#include "element.h"
#include "dofmanager.h"
#include "range.h"

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
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    tstep_all_out_flag  = ir->hasField(_IFT_OutputManager_tstepall);
    tstep_step_out      = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, tstep_step_out, _IFT_OutputManager_tstepstep);

    dofman_all_out_flag = ir->hasField(_IFT_OutputManager_dofmanall);
    element_all_out_flag = ir->hasField(_IFT_OutputManager_elementall);

    tsteps_out.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, tsteps_out, _IFT_OutputManager_tstepsout);
    dofman_out.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, dofman_out, _IFT_OutputManager_dofmanoutput);
    dofman_except.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, dofman_except, _IFT_OutputManager_dofmanexcept);
    element_out.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, element_out, _IFT_OutputManager_elementoutput);
    element_except.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, element_except, _IFT_OutputManager_elementexcept);

    return IRRT_OK;
}

void
OutputManager :: doDofManOutput(FILE *file, TimeStep *tStep)
{
    int ndofman = domain->giveNumberOfDofManagers();

    if ( !testTimeStepOutput(tStep) ) {
        return;
    }

    fprintf(file, "\n\nDofManager output:\n------------------\n");

    if ( dofman_all_out_flag   && dofman_except.empty() ) {
        for ( int i = 1; i <= ndofman; i++ ) {
            // test for null dof in parallel mode
            if ( domain->giveDofManager(i)->giveParallelMode() == DofManager_null ) {
                continue;
            }

            domain->giveDofManager(i)->printOutputAt(file, tStep);
        }
    } else {
        for ( int i = 1; i <= ndofman; i++ ) {
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
    if ( !testTimeStepOutput(tStep) ) {
        return;
    }

    fprintf(file, "\n\nElement output:\n---------------\n");

    if ( element_all_out_flag   && element_except.empty() ) {
        for ( auto &elem : domain->giveElements() ) {
            // test for remote element in parallel mode
            if ( elem->giveParallelMode() == Element_remote ) {
                continue;
            }

            elem->printOutputAt(file, tStep);
        }
    } else {
        int nelem = domain->giveNumberOfElements();
        for ( int i = 1; i <= nelem; i++ ) {
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

    // test for null dof in parallel mode
    if ( domain->giveDofManager(number)->giveParallelMode() == DofManager_null ) {
        return 0;
    }

    // test all_select flag on
    if ( dofman_all_out_flag ) {
        selected  = 1;
    } else {
        // test for particular dofman selection
        int _label = domain->giveDofManager(number)->giveLabel();
        for ( Range &range: dofman_out ) {
            if ( range.test(_label) ) {
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
    int _label = domain->giveDofManager(number)->giveLabel();

    for ( Range &range: this->dofman_except ) {
        // test if excluded
        if ( range.test(_label) ) {
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

    // test for remote element in parallel mode
    if ( domain->giveElement(number)->giveParallelMode() == Element_remote ) {
        return 0;
    }

    // test all_select flag on
    if ( element_all_out_flag ) {
        selected  = 1;
    } else {
        // test for particular element selection
        int _label = domain->giveElement(number)->giveLabel();

        for ( Range &range: this->element_out ) {
            if ( range.test(_label) ) {
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
    int _label = domain->giveElement(number)->giveLabel();

    for ( Range &range: element_except ) {
        // test if excluded
        if ( range.test(_label) ) {
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

    for ( Range &range: this->tsteps_out ) {
        // test if INCLUDED
        if ( range.test( tStep->giveNumber() ) ) {
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

void
OutputManager :: beCopyOf(OutputManager *om)
{
    this->tstep_all_out_flag = om->tstep_all_out_flag;
    this->tstep_step_out = om->tstep_step_out;
    this->dofman_all_out_flag = om->dofman_all_out_flag;
    this->element_all_out_flag = om->element_all_out_flag;
}
} // end namespace oofem
