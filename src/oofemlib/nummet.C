/* $Header: /home/cvs/bp/oofem/oofemlib/src/nummet.C,v 1.4 2003/04/06 14:08:25 bp Exp $ */
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

//
// file: nummet.cc
//

#include "nummet.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "ldltfact.h"


/*
 * NumericalMethod* NumericalMethod :: ofType (char* aClass)
 * // Returns a new element, which has the same number than the receiver,
 * // but belongs to aClass (PlaneStrain, or Truss2D,..).
 * {
 * NumericalMethod* newNModel ;
 *
 * if (! strncasecmp(aClass,"ldltfactorization",16))
 *    newNModel = new LDLTFactorization(number,domain,engngModel) ;
 * //   else if (! strncasecmp(aClass,"t2",2))
 * //      newNModel = new Truss2D(number,domain) ;
 * else {
 *    printf ("%s : unknown method type \n",aClass) ;
 *    exit(0) ;}
 *
 * return newNModel ;
 * }
 *
 */
