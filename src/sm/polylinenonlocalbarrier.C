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

#include "polylinenonlocalbarrier.h"
#include "domain.h"
#include "node.h"
#include "intarray.h"
#include "flotarry.h"
#include "mathfem.h"

namespace oofem {
PolylineNonlocalBarrier :: PolylineNonlocalBarrier(int n, Domain *aDomain) :
    NonlocalBarrier(n, aDomain), vertexNodes()
    // Constructor. Creates an element with number n, belonging to aDomain.
{
    localXCoordIndx = 1;
    localYCoordIndx = 2;
}


PolylineNonlocalBarrier :: ~PolylineNonlocalBarrier()
// Destructor.
{ }


void
PolylineNonlocalBarrier :: applyConstraint(const FloatArray &c1, const FloatArray &c2, double &weight,
                                           bool &shieldFlag, NonlocalMaterialExtensionInterface *nei)
{
    if ( this->isActivated(c1, c2) ) {
        weight = 0.0;
        shieldFlag = true;
    } else {
        shieldFlag = false;
    }
}


bool
PolylineNonlocalBarrier :: isActivated(const FloatArray &c1, const FloatArray &c2)
{
    int size = vertexNodes.giveSize();
    int indx;
    double xc1, xc2, xa, xb, yc1, yc2, ya, yb;
    double a11, a12, a21, a22, b1, b2, det, t, s;
    Node *A, *B;

    int mci = max(localXCoordIndx, localYCoordIndx);
    if ( ( c1.giveSize() > mci ) || ( c2.giveSize() > mci ) ) {
        _error("PolylineNonlocalBarrier::isActivated: local coordinate index size violation");
    }

    xc1 = c1.at(localXCoordIndx);
    yc1 = c1.at(localYCoordIndx);
    xc2 = c2.at(localXCoordIndx);
    yc2 = c2.at(localYCoordIndx);

    for ( indx = 1; indx < size; indx++ ) {
        A = domain->giveNode( vertexNodes.at(indx) );
        B = domain->giveNode( vertexNodes.at(indx + 1) );

        xa = A->giveCoordinate(localXCoordIndx);
        ya = A->giveCoordinate(localYCoordIndx);
        xb = B->giveCoordinate(localXCoordIndx);
        yb = B->giveCoordinate(localYCoordIndx);

        a11 = xc2 - xc1;
        a12 = xa - xb;
        a21 = yc2 - yc1;
        a22 = ya - yb;
        b1  = xa - xc1;
        b2  = ya - yc1;
        det = a11 * a22 - a21 * a12;
        if ( det == 0 ) {
            continue;
        }

        t = ( b1 * a22 - b2 * a12 ) / det;
        if ( t < 0. || t > 1. ) {
            continue;
        }

        s = ( -b1 * a21 + b2 * a11 ) / det;
        if ( s >= 0. && s <= 1. ) {
            return true;
        }
    }

    return false;
}

IRResultType
PolylineNonlocalBarrier :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, vertexNodes, IFT_PolylineNonlocalBarrier_vertexnodes, "vertexnodes"); // Macro

    // default: polyline in xy plane
    localXCoordIndx = 1;
    localYCoordIndx = 2;

    IR_GIVE_OPTIONAL_FIELD(ir, localXCoordIndx, IFT_PolylineNonlocalBarrier_xcoordindx, "xcoordindx"); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, localYCoordIndx, IFT_PolylineNonlocalBarrier_ycoordindx, "ycoordindx"); // Macro

    return IRRT_OK;
}
} // end namespace oofem
