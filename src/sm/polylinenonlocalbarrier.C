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
 *               Copyright (C) 1993 - 2025   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "polylinenonlocalbarrier.h"
#include "domain.h"
#include "node.h"
#include "intarray.h"
#include "floatarray.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_NonlocalBarrier(PolylineNonlocalBarrier)

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
PolylineNonlocalBarrier :: applyConstraint(const double cl, const FloatArray &c1, const FloatArray &c2, double &weight,
                                           bool &shieldFlag, const NonlocalMaterialExtensionInterface &nei)
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

double
PolylineNonlocalBarrier :: calculateMinimumDistanceFromBoundary(const FloatArray &coords)
{
    double min = 1.e10;
    double tempDistance;
    //Loop over all linear sections forming the nonlocal boundary to find the minimum distance
    for ( int i = 1; i < vertexNodes.giveSize(); i++ ) {
        //Get the coordinates of the vertices
        FloatArray coordsA(coords);
        FloatArray coordsB(coords);
        for ( int j = 1; j <= coords.giveSize(); j++ ) {
            coordsA.at(j) =  domain->giveNode( vertexNodes.at(i) )->giveCoordinate(j);
            coordsB.at(j) = domain->giveNode( vertexNodes.at(i + 1) )->giveCoordinate(j);
        }

        //Get the distance from the line segment described by the vertices
        tempDistance = giveDistancePointLine(coordsA, coordsB, coords);
        if ( min > tempDistance ) { //Check if it is smaller than the minimum value
            min = tempDistance;
        }
    }

    return min;
}

double
PolylineNonlocalBarrier :: giveDistancePointLine(const FloatArray &coordsA, const FloatArray &coordsB, const FloatArray &coordsGP)
{
    FloatArray lineAGP(coordsGP); //Line start A, Line End Gauss Point
    lineAGP.subtract(coordsA);
    FloatArray lineAB(coordsB); //Line Start A, Line End B
    lineAB.subtract(coordsA);

    if ( lineAB.computeNorm() == 0. ) { //Check if A,B coincide
        return lineAGP.computeNorm();
    }

    //Since vectors AB and AP are collinear and have the same direction: AP=scaleFactor*AB
    double scaleFactor;
    scaleFactor = ( lineAGP.dotProduct(lineAB) ) / lineAB.computeSquaredNorm();
    if ( scaleFactor < 0. || scaleFactor > 1. ) { //Check if P is outside line segment AB
        FloatArray lineBGP(coordsGP); //Line start B, Line End Gauss Point
        lineBGP.subtract(coordsB);
        //Return minimum of A-Gauss Point and B-Gauss Point
        if ( lineAGP.computeNorm() < lineBGP.computeNorm() ) {
            return lineAGP.computeNorm();
        } else {
            return lineBGP.computeNorm();
        }
    } else {
        // Find coordinates of Point P = A + AB*scaleFactor
        lineAB.times(scaleFactor);
        FloatArray coordsP(coordsA);
        coordsP.add(lineAB);
        //Calculate distance from Point P to Gauss Point
        FloatArray linePGP(coordsP); //Line start P, Line End Gauss Point
        linePGP.subtract(coordsGP);
        return linePGP.computeNorm();
    }
}


void
PolylineNonlocalBarrier :: initializeFrom(InputRecord &ir)
{
    IR_GIVE_FIELD(ir, vertexNodes, _IFT_PolylineNonlocalBarrier_vertexnodes);

    // default: polyline in xy plane
    localXCoordIndx = 1;
    localYCoordIndx = 2;

    IR_GIVE_OPTIONAL_FIELD(ir, localXCoordIndx, _IFT_PolylineNonlocalBarrier_xcoordindx);
    IR_GIVE_OPTIONAL_FIELD(ir, localYCoordIndx, _IFT_PolylineNonlocalBarrier_ycoordindx);
}
} // end namespace oofem
