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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#include "traxisym1_ht.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "engngm.h"
#include "structuralms.h"
#include "load.h"
#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

namespace oofem {
TrAxisym1_ht :: TrAxisym1_ht(int n, Domain *aDomain, ElementMode em) :
    Tr1_ht(n, aDomain, em)
    // Constructor.
{ }

TrAxisym1_ht :: ~TrAxisym1_ht()
// Destructor
{ }

double
TrAxisym1_ht :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
    double area, weight;

    weight  = aGaussPoint->giveWeight();
    area    = this->giveArea();

    return 2.0 *area *weight *this->computeRadiusAt(aGaussPoint);
}

double
TrAxisym1_ht :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double dx, dy, length, radius;
    Node *nodeA, *nodeB;
    int aNode = 0, bNode = 0;
    FloatMatrix n;

    if ( iEdge == 1 ) { // edge between nodes 1 2
        aNode = 1;
        bNode = 2;
    } else if ( iEdge == 2 ) { // edge between nodes 2 3
        aNode = 2;
        bNode = 3;
    } else if ( iEdge == 3 ) { // edge between nodes 3 4
        aNode = 3;
        bNode = 1;
    } else {
        _error("computeEdgeVolumeAround: wrong egde number");
    }

    nodeA   = this->giveNode(aNode);
    nodeB   = this->giveNode(bNode);

    dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
    dy      = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
    length = sqrt(dx * dx + dy * dy);
    this->computeEgdeNMatrixAt(n, gp);
    radius = n.at(1, 1) * nodeA->giveCoordinate(1) + n.at(1, 2) * nodeB->giveCoordinate(1);
    return 0.5 *length *radius *gp->giveWeight();
}

double
TrAxisym1_ht :: computeRadiusAt(GaussPoint *gp)
{
    double r = 0.0;
    int i;
    FloatMatrix n;

    this->computeNmatrixAt( n, gp->giveCoordinates() );
    for ( i = 1; i <= 3; i++ ) {
        r += n.at(1, i) * this->giveNode(i)->giveCoordinate(1);
    }

    return r;
}
} // end namespace oofem
