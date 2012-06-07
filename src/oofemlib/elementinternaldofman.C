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

#include "elementinternaldofman.h"
#include "dof.h"
#include "intarray.h"
#include "verbose.h"

namespace oofem {

ElementDofManager :: ElementDofManager(int n, Domain *aDomain, Element* elem) :
    DofManager(n, aDomain)
{
    this->element = elem;
}


ElementDofManager :: ~ElementDofManager()
{ }


IRResultType ElementDofManager :: initializeFrom(InputRecord *ir)
// Gets from the source line from the data file all the data of the receiver.
{
#  ifdef VERBOSE
    // VERBOSE_PRINT1("Instanciating node ",number)
#  endif

    return DofManager :: initializeFrom(ir);
}


void ElementDofManager :: printYourself()
// Prints the receiver on screen.
{
    int i;

    printf("InternalElementDofManager %d \n", number);
    for ( i = 0; i < numberOfDofs; i++ ) {
        if ( dofArray [ i ] ) {
            dofArray [ i ]->printYourself();
        } else {
            printf("dof %d is nil \n", i + 1);
        }
    }

    loadArray.printYourself();
    printf("\n");
}

} // end namespace oofem
