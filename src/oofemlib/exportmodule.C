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

#include "exportmodule.h"
#include "timestep.h"
#include "engngm.h"
#include "domain.h"

#include <cstdarg>


namespace oofem {
ExportModule :: ExportModule(int n, EngngModel *e) : tsteps_out(), domainMask(), regionSets(), defaultElementSet( 0, e->giveDomain(1) )
{
    this->number = n;
    emodel = e;
    regionSets.resize(0);
    timeScale = 1.;
}


ExportModule :: ~ExportModule()
{ }


IRResultType
ExportModule :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    tstep_all_out_flag = ir->hasField(_IFT_ExportModule_tstepall);

    tstep_step_out = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, tstep_step_out, _IFT_ExportModule_tstepstep);

    IR_GIVE_OPTIONAL_FIELD(ir, tsteps_out, _IFT_ExportModule_tstepsout);

    tstep_substeps_out_flag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, tstep_substeps_out_flag, _IFT_ExportModule_subtstepsout);

    domain_all_flag = ir->hasField(_IFT_ExportModule_domainall);

    if ( !domain_all_flag ) {
        domainMask.clear();
        IR_GIVE_OPTIONAL_FIELD(ir, domainMask, _IFT_ExportModule_domainmask);
    }

    IR_GIVE_OPTIONAL_FIELD(ir, regionSets, _IFT_ExportModule_regionsets); // Macro

    IR_GIVE_OPTIONAL_FIELD(ir, timeScale, _IFT_ExportModule_timescale); // Macro

    return IRRT_OK;
}

void ExportModule :: initialize(){
  initializeElementSet();
}

void ExportModule :: initializeElementSet(){
    if ( regionSets.isEmpty() ){
        // default: whole domain region
        defaultElementSet.clear();
        defaultElementSet.setDomain(emodel->giveDomain(1));
        defaultElementSet.addAllElements();
    }
}


int ExportModule :: giveNumberOfRegions()
{
    // Returns number of regions (aka sets)
  if (regionSets.isEmpty())
    return 1; // defaultElementSet
  else 
    return this->regionSets.giveSize();
}

Set *ExportModule :: giveRegionSet(int i)
{
  if (regionSets.isEmpty()) {
    return & this->defaultElementSet;
  }  else {
    int setid = regionSets.at(i);
    if ( setid > 0 ) {
        return emodel->giveDomain(1)->giveSet(setid);
    } else {
        OOFEM_ERROR ("Bad set index");
        return & this->defaultElementSet;
    }
  }
}

std :: string
ExportModule :: giveOutputBaseFileName(TimeStep *tStep)
{
    char fext [ 100 ];

    if ( this->testSubStepOutput() ) {
        // include tStep version in output file name
        if ( this->emodel->isParallel() && this->emodel->giveNumberOfProcesses() > 1 ) {
            sprintf( fext, "_%03d.m%d.%d.%d", emodel->giveRank(), this->number, tStep->giveNumber(), tStep->giveSubStepNumber() );
        } else {
            sprintf( fext, ".m%d.%d.%d", this->number, tStep->giveNumber(), tStep->giveSubStepNumber() );
        }
        return this->emodel->giveOutputBaseFileName() + fext;
    } else {
        if ( this->emodel->isParallel() && this->emodel->giveNumberOfProcesses() > 1 ) {
            sprintf( fext, "_%03d.m%d.%d", emodel->giveRank(), this->number, tStep->giveNumber() );
        } else {
            sprintf( fext, ".m%d.%d", this->number, tStep->giveNumber() );
        }
        return this->emodel->giveOutputBaseFileName() + fext;
    }
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

    for ( auto &step: tsteps_out ) {
        // test if INCLUDED
        if ( step.test( tStep->giveNumber() ) ) {
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

std :: string ExportModule :: errorInfo(const char *func) const
{
    return std :: string(this->giveClassName()) + "::" + func;
}
} // end namespace oofem
