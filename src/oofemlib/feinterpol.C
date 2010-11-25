/* $Header: /home/cvs/bp/oofem/oofemlib/src/feinterpol.h,v 1.1 2003/04/06 14:08:24 bp Exp $ */
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

#include "feinterpol.h"
#include "element.h"

namespace oofem {
/*
 * void
 * FEInterpolation :: nodes2coords(Domain *d, IntArray &nodes, const FloatArray **c, int n)
 * {
 *  int i, nnode = nodes.giveSize();
 *  if ( n < nnode ) {
 *      OOFEM_ERROR("FEInterpolation::nodes2coords: size mismatch");
 *  }
 *
 *  for ( i = 0; i < nnode; i++ ) {
 *      c [ i ] = d->giveNode( nodes(i) )->giveCoordinates();
 *  }
 * }
 */

int FEIElementGeometryWrapper :: giveNumberOfVertices() const { return elem->giveNumberOfNodes(); }
const FloatArray *FEIElementGeometryWrapper :: giveVertexCoordinates(int i) const { return elem->giveNode(i)->giveCoordinates(); }
} // end namespace oofem
