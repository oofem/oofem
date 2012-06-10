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

#include "exportmodule.h"
#include "timestep.h"
#include "engngm.h"
#include "oofem_limits.h"
#include "range.h"

#ifndef __MAKEDEPEND
 #include <cstdarg>
#endif

namespace oofem {
ExportModule :: ExportModule(int n, EngngModel *e) : tsteps_out(), domainMask()
{
    this->number = n;
    emodel = e;
}


ExportModule :: ~ExportModule()
{ }


IRResultType
ExportModule :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    tstep_all_out_flag = ir->hasField(IFT_ExportModule_tstepall, "tstep_all");

    tstep_step_out = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, tstep_step_out, IFT_ExportModule_tstepstep, "tstep_step"); // Macro

    IR_GIVE_OPTIONAL_FIELD(ir, tsteps_out, IFT_ExportModule_tstepsout, "tsteps_out"); // Macro

    domain_all_flag = ir->hasField(IFT_ExportModule_domainall, "domain_all");

    if ( !domain_all_flag ) {
        domainMask.resize(0);
        IR_GIVE_OPTIONAL_FIELD(ir, domainMask, IFT_ExportModule_domainmask, "domain_mask"); // Macro
    }

    return IRRT_OK;
}

std::string
ExportModule :: giveOutputBaseFileName(TimeStep *tStep)
{
    char fext[100];
#ifdef __PARALLEL_MODE
    if (this->emodel->giveNumberOfProcesses() > 1)
        sprintf( fext, "_%03d.m%d.%d", emodel->giveRank(), this->number, tStep->giveNumber() );
    else
#endif
    sprintf( fext, ".m%d.%d", this->number, tStep->giveNumber() );
    return this->emodel->giveOutputBaseFileName() + fext;
}

bool
ExportModule :: testTimeStepOutput(TimeStep *tStep)
{
    if ( tstep_all_out_flag ) {
        return true;
    }

    if ( tstep_step_out ) {
        //if (((tStep->giveNumber()-emodel->giveNumberOfFirstStep()) % tstep_step_out) == 0) return 1;
        if ( ( ( tStep->giveNumber() ) % tstep_step_out ) == 0 ) {
            return 1;
        }
    }

    dynaList< Range > :: iterator tstepsIter;
    for ( tstepsIter = tsteps_out.begin(); tstepsIter != tsteps_out.end(); ++tstepsIter ) {
        // test if INCLUDED
        if ( ( * tstepsIter ).test( tStep->giveNumber() ) ) {
            return true;
        }
    }

    return 0;
}

bool
ExportModule :: testDomainOutput(int n)
{
    if ( domain_all_flag ) {
        return true;
    }

    return domainMask.findFirstIndexOf(n);
}

void ExportModule :: error(const char *file, int line, const char *format, ...) const
{
    char buffer [ MAX_ERROR_MSG_LENGTH ];
    va_list args;

    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);

    __OOFEM_ERROR3(file, line, "Class: %s\n%s", this->giveClassName(), buffer);
}
} // end namespace oofem
