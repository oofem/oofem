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

#include "initmodule.h"
//#include "timestep.h"
//#include "engngm.h"
//#include "logger.h"
#include "oofem_limits.h"

#ifndef __MAKEDEPEND
 #include <stdarg.h>
#endif

namespace oofem {
InitModule :: InitModule(EngngModel *e)
{
    emodel = e;
}


InitModule :: ~InitModule()
{ }


IRResultType
InitModule :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    char initFileName [ MAX_FILENAME_LENGTH ];
    IR_GIVE_FIELD2(ir, initFileName, IFT_InitModule_initfilename, "initfile", MAX_FILENAME_LENGTH); // Macro
    if ( ( initStream = fopen(initFileName, "r") ) == NULL ) {
        OOFEM_ERROR2("InitModule::initializeFrom: failed to open file %s", initFileName);
    }

    return IRRT_OK;
}
} // end namespace oofem
