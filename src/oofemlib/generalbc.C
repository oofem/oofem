/* $Header: /home/cvs/bp/oofem/oofemlib/src/generalbc.C,v 1.4.4.1 2004/04/05 15:19:43 bp Exp $ */
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

#include "generalbc.h"
#include "deadwght.h"
#include "nodload.h"
#include "boundary.h"
#include "initial.h"
#include "debug.h"
#include "verbose.h"
#include "usrdefsub.h"
#include "timestep.h"
#include "bcvaltype.h"

#ifndef __MAKEDEPEND
#include <stdlib.h>
#ifdef HAVE_STRINGS_H
#include <strings.h>
#endif
#endif

namespace oofem {

GeneralBoundaryCondition :: GeneralBoundaryCondition(int n, Domain *d) : FEMComponent(n, d)
{
    loadTimeFunction = 0;
}



LoadTimeFunction *GeneralBoundaryCondition :: giveLoadTimeFunction()
// Returns the load-time function of the receiver. Reads its number in the
// data file if has not been done yet.
{
    if ( !loadTimeFunction ) {
        _error("giveLoadTimeFunction: LoadTimeFunction is not defined");
    }

    return domain->giveLoadTimeFunction(loadTimeFunction);
}


GeneralBoundaryCondition *GeneralBoundaryCondition :: ofType(char *aClass)
// Returns a new load, which has the same number than the receiver,
// but belongs to aClass (NodalLoad, DeadWeight,..).
{
    GeneralBoundaryCondition *newBC;

    if ( !strncasecmp(aClass, "boundarycondition", 5) ) {
        newBC = new BoundaryCondition(number, domain);
    } else if ( !strncasecmp(aClass, "deadweight", 5) )   {
        newBC = new DeadWeight(number, domain);
    }
    //else if (! strncasecmp(aClass,"initialcondition",5))
    //   newBC = new InitialCondition(number,domain) ;
    else if ( !strncasecmp(aClass, "nodalload", 5) ) {
        newBC = new NodalLoad(number, domain);
    }
    // else if (!strncasecmp(aClass,"temperatureload",12))
    //  newBC = new TemperatureLoad(number,domain) ;
    else {
        // last resort - call aditional user defined subroutine
        newBC = CreateUsrDefBoundaryConditionOfType(aClass, number, domain);
        if ( newBC == NULL ) {
            _error2("ofType:  unknown bc type (%s)", aClass);
            exit(0);
        }
    }

    return newBC;
}

IRResultType
GeneralBoundaryCondition :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;           // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, loadTimeFunction, IFT_GeneralBoundaryCondition_LoadTimeFunct, "loadtimefunction"); // Macro
    if ( loadTimeFunction <= 0 ) {
        _error("initializeFrom: bad loadtimefunction id");
    }

    int val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_GeneralBoundaryCondition_valType, "valtype"); // Macro
    valType = ( bcValType ) val;


    return IRRT_OK;
}


int
GeneralBoundaryCondition :: giveInputRecordString(std :: string &str, bool keyword)
{
    char buff [ 1024 ];

    FEMComponent :: giveInputRecordString(str, keyword);
    sprintf(buff, " loadtimefunction %d", this->loadTimeFunction);
    str += buff;

    return 1;
}

} // end namespace oofem
