/* $Header: /home/cvs/bp/oofem/tm/src/transportelement.C,v 1.3.4.1 2004/04/05 15:19:53 bp Exp $ */
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

#include "fmelement.h"
#include "domain.h"
#include "timestep.h"
#include "node.h"
#include "dof.h"
#include "load.h"
#include "boundaryload.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "intarray.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "debug.h"
#include "verbose.h"

#include "elementside.h"
#include "mathfem.h"
#ifndef __MAKEDEPEND
#include <stdlib.h>
#include <stdio.h>
#endif

FMElement :: FMElement(int n, Domain *aDomain) :
    Element(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{ }


FMElement :: ~FMElement()
// Destructor.
{ }

IRResultType
FMElement :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                            // Required by IR_GIVE_FIELD macro

    Element :: initializeFrom(ir);
    IR_GIVE_OPTIONAL_FIELD(ir, boundarySides, IFT_SUPGElement_bsides, "bsides"); // Macro
    if ( !boundarySides.isEmpty() ) {
        IR_GIVE_FIELD(ir, boundaryCodes, IFT_SUPGElement_bcodes, "bcodes"); // Macro
    }

    //this -> computeGaussPoints();
    return IRRT_OK;
}

