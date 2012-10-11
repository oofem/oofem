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
 *               Copyright (C) 1993 - 2010   Borek Patzak
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

#include "load.h"
#include "reinforcement.h"
#include "nodload.h"
#include "boundary.h"
#include "initial.h"
//#include "temperatureload.h"
#include "verbose.h"
#include "usrdefsub.h"
#include "timestep.h"
#ifndef __MAKEDEPEND
 #include <stdlib.h>
#endif

namespace oofem {




IRResultType
Reinforcement :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

#  ifdef VERBOSE
    // VERBOSE_PRINT1 ("Instanciating load ",number)
#  endif
    //Load :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, porosity, IFT_Reinfocement_porosity, "porosity"); // Macro
    IR_GIVE_FIELD(ir, shapefactor,  IFT_Reinfocement_permeability, "shapefactor"); // Macro
    IR_GIVE_FIELD(ir, permeability, IFT_Reinfocement_shapeFactor, "permeability"); // Macro
    
    return IRRT_OK;
}




} // end namespace oofem
