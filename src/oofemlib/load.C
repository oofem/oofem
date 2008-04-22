/* $Header: /home/cvs/bp/oofem/oofemlib/src/load.C,v 1.8 2003/04/18 10:37:07 bp Exp $ */
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


//   file LOAD.C

#include "load.h"
#include "deadwght.h"
#include "nodload.h"
#include "boundary.h"
#include "initial.h"
//#include "temperatureload.h"
#include "debug.h"
#include "verbose.h"
#include "usrdefsub.h"
#include "timestep.h"
#ifndef __MAKEDEPEND
#include <stdlib.h>
#endif

Load :: Load(int i, Domain *aDomain) :
    GeneralBoundaryCondition(i, aDomain), componentArray(), dofExcludeMask()
    // Constructor. Creates a load with number i, belonging to aDomain.
{
    loadTimeFunction = 0;
}


FloatArray &Load :: giveComponentArray()
// Returns the array that contains the components of the receivert. If
// this array does not exist yet, forms it by reading its values in the
// data file.
{
    return componentArray;
}


void
Load :: computeComponentArrayAt(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
// Returns an array, the force induced at stepN by the receiver.
{
    // FloatArray* force ;
    double factor;

    /*
     * factor = this -> giveLoadTimeFunction() -> at(stepN->giveTime()) ;
     * if ((mode!=VM_Incremental)&&(mode!=VM_Total))
     * _error ("computeComponentArrayAt: unknown mode");
     * if ((mode==VM_Incremental) && (!stepN->isTheFirstStep()))
     * //   factor -= this->giveLoadTimeFunction()->at(stepN->givePreviousStep()->giveTime()) ;
     * factor -= this->giveLoadTimeFunction()->at(stepN->giveTime()-stepN->giveTimeIncrement());
     */
    factor = this->giveLoadTimeFunction()->evaluate(stepN, mode);
    answer  = this->giveComponentArray();
    answer.times(factor);
    return;
}


IRResultType
Load :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

#  ifdef VERBOSE
    // VERBOSE_PRINT1 ("Instanciating load ",number)
#  endif
    GeneralBoundaryCondition :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, componentArray, IFT_Load_components, "components"); // Macro

    int size = componentArray.giveSize();
    dofExcludeMask.resize(size);
    dofExcludeMask.zero();
    IR_GIVE_OPTIONAL_FIELD(ir, dofExcludeMask, IFT_Load_dofexcludemask, "dofexcludemask"); // Macro
    if ( dofExcludeMask.giveSize() != size ) {
        _error("initializeFrom: dofExcludeMask and componentArray size mismatch");
    } else {
        for ( int i = 1; i <= size; i++ ) {
            if ( dofExcludeMask.at(i) ) {
                componentArray.at(i) = 0.0;
            }
        }
    }

    return IRRT_OK;
}


int
Load :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    GeneralBoundaryCondition :: giveInputRecordString(str, keyword);
    sprintf( buff, " components %d", this->componentArray.giveSize() );
    str += buff;
    for ( int i = 1; i <= this->componentArray.giveSize(); i++ ) {
        sprintf( buff, " %e", this->componentArray.at(i) );
        str += buff;
    }

    return 1;
}


int
Load :: isDofExcluded(int indx)
{
    if ( ( indx > 0 ) && ( indx <= dofExcludeMask.giveSize() ) ) {
        return dofExcludeMask.at(indx);
    } else {
        _error("isDofExcluded: dof indx out of range");
    }

    return 0;
}



