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

#include "fmelement.h"

namespace oofem {
FMElement :: FMElement(int n, Domain *aDomain) :
    Element(n, aDomain)
{ }


FMElement :: ~FMElement()
{ }

IRResultType
FMElement :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom";
    IRResultType result;

    IR_GIVE_OPTIONAL_FIELD(ir, boundarySides, IFT_SUPGElement_bsides, "bsides"); // Macro
    // if ( !boundarySides.isEmpty() ) {
    //  IR_GIVE_FIELD(ir, boundaryCodes, IFT_SUPGElement_bcodes, "bcodes"); // Macro
    //}

    return Element :: initializeFrom(ir);
}
} // end namespace oofem
