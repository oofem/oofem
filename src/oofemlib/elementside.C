/* $Header: /home/cvs/bp/oofem/oofemlib/src/elementside.C,v 1.6 2003/04/06 14:08:24 bp Exp $ */
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


//   file ElementSide.C

#include "elementside.h"
#include "dof.h"
#include "nodload.h"
#include "timestep.h"

#include "flotarry.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "debug.h"
#include "verbose.h"
#ifndef __MAKEDEPEND
#include <math.h>
#include <stdlib.h>
#endif

namespace oofem {

ElementSide :: ElementSide(int n, Domain *aDomain) :
    DofManager(n, aDomain)
    // Constructor. Creates a node with number n, belonging to aDomain.
{ }




ElementSide :: ~ElementSide()
// Destructor.
{ }

IRResultType
ElementSide :: initializeFrom(InputRecord *ir)
// Gets from the source line from the data file all the data of the receiver.
{
    // int j, size ;
    // IntArray *dofIDArry;

#  ifdef VERBOSE
    // VERBOSE_PRINT1("Instanciating node ",number)
#  endif

    return DofManager :: initializeFrom(ir);
}



void ElementSide :: printYourself()
// Prints the receiver on screen.
{
    int i;
    // double x,y ;

    printf("Side %d \n", number);
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

void
ElementSide :: computeTransformation(FloatMatrix &answer, const IntArray *map)
{
    //
    // computes transfromation of receiver from global cs to nodal (user-defined) cs.
    // Note: implementation rely on D_u, D_v and D_w (R_u, R_v, R_w) order in cltypes.h
    // file. Do not change their order and do not insert any values between these values.
    //
    //

    int i, size;

    if ( map == NULL ) {
        size = numberOfDofs;
    } else {
        size = map->giveSize();
    }


    // response for all local dofs is computed

    answer.resize(size, size);
    answer.zero();

    for ( i = 1; i <= size; i++ ) {
        answer.at(i, i) = 1.0;
    }
}

} // end namespace oofem
