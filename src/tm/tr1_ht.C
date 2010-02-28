/* $Header: /home/cvs/bp/oofem/tm/src/tr1_ht.C,v 1.2 2003/04/23 14:22:15 bp Exp $ */
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

#include "tr1_ht.h"
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
 #include <math.h>
 #include <stdio.h>
#endif

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

namespace oofem {
Tr1_ht :: Tr1_ht(int n, Domain *aDomain, ElementMode em) :
    TransportElement(n, aDomain, em)
    // Constructor.
{
    numberOfDofMans  = 3;
    area = -1.0;
    numberOfGaussPoints = 1;
}

Tr1_ht :: ~Tr1_ht()
// Destructor
{ }


double
Tr1_ht :: giveArea()
// returns the area occupied by the receiver
{
    if ( area > 0 ) {
        return area;         // check if previously computed
    }

    Node *node1, *node2, *node3;
    double x1, x2, x3, y1, y2, y3;

    node1 = this->giveNode(1);
    node2 = this->giveNode(2);
    node3 = this->giveNode(3);

    x1 = node1->giveCoordinate(1);
    x2 = node2->giveCoordinate(1);
    x3 = node3->giveCoordinate(1);

    y1 = node1->giveCoordinate(2);
    y2 = node2->giveCoordinate(2);
    y3 = node3->giveCoordinate(2);

    return ( area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 ) );
}

void
Tr1_ht :: computeNSubMatrixAt(FloatMatrix &answer, FloatArray *coords)
// Returns the displacement interpolation matrix {N} of the receiver,
// evaluated at aGaussPoint.
{
    double l1, l2, l3;

    l1 = coords->at(1);
    l2 = coords->at(2);
    l3 = 1.0 - l1 - l2;

    answer.resize(1, 3);
    answer.zero();

    answer.at(1, 1) = l1;
    answer.at(1, 2) = l2;
    answer.at(1, 3) = l3;

    return;
}

void
Tr1_ht :: computeNmatrixAt(FloatMatrix &answer, FloatArray *coords)
{
    this->computeNSubMatrixAt(answer, coords);
}


void
Tr1_ht :: computeGradientMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint)
{
    Node *node1, *node2, *node3;
    double x1, x2, x3, y1, y2, y3, area;

    node1 = this->giveNode(1);
    node2 = this->giveNode(2);
    node3 = this->giveNode(3);

    x1 = node1->giveCoordinate(1);
    x2 = node2->giveCoordinate(1);
    x3 = node3->giveCoordinate(1);

    y1 = node1->giveCoordinate(2);
    y2 = node2->giveCoordinate(2);
    y3 = node3->giveCoordinate(2);

    area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );

    answer.resize(2, 3);

    answer.at(1, 1) = y2 - y3;
    answer.at(1, 2) = y3 - y1;
    answer.at(1, 3) = y1 - y2;

    answer.at(2, 1) = x3 - x2;
    answer.at(2, 2) = x1 - x3;
    answer.at(2, 3) = x2 - x1;

    answer.times( 1. / ( 2. * area ) );
    return;
}

void
Tr1_ht :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Triangle, numberOfGaussPoints, _2dHeat);
    }
}

void
Tr1_ht ::   giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const {
    // returns DofId mask array for inode element node.
    // DofId mask array determines the dof ordering requsted from node.
    // DofId mask array contains the DofID constants (defined in cltypes.h)
    // describing physical meaning of particular DOFs.
    answer.resize(1);
    answer.at(1) = T_f;
}


IRResultType
Tr1_ht :: initializeFrom(InputRecord *ir)
{
    //const char *__keyword, *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    //IRResultType result;                               // Required by IR_GIVE_FIELD macro

    this->Element :: initializeFrom(ir);

    numberOfGaussPoints = 1;
    //IR_GIVE_OPTIONAL_FIELD (ir, numberOfGaussPoints, "nip"); // Macro
    //if ( numberOfGaussPoints != 1) numberOfGaussPoints = 1;

    this->computeGaussPoints();
    return IRRT_OK;
}


double
Tr1_ht :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
    double area, weight;

    weight  = aGaussPoint->giveWeight();
    area    = this->giveArea();

    return 2.0 *area *weight *this-> giveCrossSection()->give('t');
}

void
Tr1_ht :: computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    /*
     *
     * computes interpolation matrix for element edge.
     * we assemble locally this matrix for only nonzero
     * shape functions.
     * (for example only two nonzero shape functions for 2 dofs are
     * necessary for linear plane stress tringle edge).
     * These nonzero shape functions are then mapped to
     * global element functions.
     *
     * Using mapping technique will allow to assemble shape functions
     * without regarding particular side
     */

    double ksi, n1, n2;
    answer.resize(1, 2);
    answer.zero();

    ksi = gp->giveCoordinate(1);
    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;

    answer.at(1, 1) = n1;
    answer.at(1, 2) = n2;

    return;
}


double
Tr1_ht :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double dx, dy, length, thick;
    Node *nodeA, *nodeB;
    int aNode = 0, bNode = 0;

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
    thick = this->giveCrossSection()->give('t');
    return 0.5 *length *thick *gp-> giveWeight();
}


void
Tr1_ht :: giveEdgeDofMapping(IntArray &answer, int iEdge)
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    answer.resize(2);
    if ( iEdge == 1 ) { // edge between nodes 1,2
        answer.at(1) = 1;
        answer.at(2) = 2;
    } else if ( iEdge == 2 ) { // edge between nodes 2 3
        answer.at(1) = 2;
        answer.at(2) = 3;
    } else if ( iEdge == 3 ) { // edge between nodes 3 4
        answer.at(1) = 3;
        answer.at(2) = 1;
    } else {
        _error("giveEdgeDofMapping: wrong edge number");
    }

    return;
}

void
Tr1_ht :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    double n1, n2, ksi = gp->giveCoordinate(1);
    Node *nodeA, *nodeB;
    int aNode = 0, bNode = 0;

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
        _error("computeEdgeIpGlobalCoords: wrong egde number");
    }

    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;

    nodeA   = this->giveNode(aNode);
    nodeB   = this->giveNode(bNode);

    answer.resize(2);
    answer.at(1) = n1 * nodeA->giveCoordinate(1) + n2 *nodeB-> giveCoordinate(1);
    answer.at(2) = n1 * nodeA->giveCoordinate(2) + n2 *nodeB-> giveCoordinate(2);
}

void
Tr1_ht :: computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *atTime, ValueModeType mode)
{
    this->computeInternalSourceRhsSubVectorAt(answer, atTime, mode, 1);
}

int
Tr1_ht :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    double l1, l2, l3;

    l1 = lcoords.at(1);
    l2 = lcoords.at(2);
    l3 = 1.0 - l1 - l2;

    answer.resize(2);
    answer.at(1) = l1 * this->giveNode(1)->giveCoordinate(1) + l2 *this-> giveNode(2)->giveCoordinate(1) +
                   l3 *this-> giveNode(3)->giveCoordinate(1);
    answer.at(2) = l1 * this->giveNode(1)->giveCoordinate(2) + l2 *this-> giveNode(2)->giveCoordinate(2) +
                   l3 *this-> giveNode(3)->giveCoordinate(2);

    return 1;
}

Interface *
Tr1_ht :: giveInterface(InterfaceType interface)
{
    if ( interface == SpatialLocalizerInterfaceType ) {
        return ( SpatialLocalizerInterface * ) this;
    } else if ( interface == EIPrimaryFieldInterfaceType ) {
        return ( EIPrimaryFieldInterface * ) this;
    } else if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return ( ZZNodalRecoveryModelInterface * ) this;
    }

    return NULL;
}

int
Tr1_ht :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( ( type == IST_TemperatureFlow ) ) {
        return 3;
    }

    return 0;
}

void
Tr1_ht :: ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint, InternalStateType type)
{
    int i;
    FloatMatrix n;
    this->computeNmatrixAt( n, aGaussPoint->giveCoordinates() );

    if ( this->giveIPValueSize(type, aGaussPoint) ) {
        answer.resize(1, 3);
    } else {
        return;
    }

    for ( i = 1; i <= 3; i++ ) {
        answer.at(1, i)  = n.at(1, i);
    }

    return;
}


int
Tr1_ht :: SpatialLocalizerI_containsPoint(const FloatArray &coords) {
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}

double
Tr1_ht :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray lcoords(3), gcoords;
    double dist;
    int size, gsize;

    lcoords.at(1) = lcoords.at(2) = lcoords.at(3) = 1. / 3.;
    this->computeGlobalCoordinates(gcoords, lcoords);

    if ( ( size = coords.giveSize() ) < ( gsize = gcoords.giveSize() ) ) {
        _error("SpatialLocalizerI_giveDistanceFromParametricCenter: coordinates size mismatch");
    }

    if ( size == gsize ) {
        dist = coords.distance(gcoords);
    } else {
        FloatArray helpCoords = coords;

        helpCoords.resize(gsize);
        dist = helpCoords.distance(gcoords);
    }

    return dist;
}


#define POINT_TOL 1.e-3

int
Tr1_ht :: computeLocalCoordinates(FloatArray &answer, const FloatArray &coords)
{
    Node *node1, *node2, *node3;
    double area, x1, x2, x3, y1, y2, y3;

    node1 = this->giveNode(1);
    node2 = this->giveNode(2);
    node3 = this->giveNode(3);

    x1 = node1->giveCoordinate(1);
    x2 = node2->giveCoordinate(1);
    x3 = node3->giveCoordinate(1);

    y1 = node1->giveCoordinate(2);
    y2 = node2->giveCoordinate(2);
    y3 = node3->giveCoordinate(2);

    area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );

    answer.resize(3);

    answer.at(1) = ( ( x2 * y3 - x3 * y2 ) + ( y2 - y3 ) * coords.at(1) + ( x3 - x2 ) * coords.at(2) ) / 2. / area;
    answer.at(2) = ( ( x3 * y1 - x1 * y3 ) + ( y3 - y1 ) * coords.at(1) + ( x1 - x3 ) * coords.at(2) ) / 2. / area;
    answer.at(3) = ( ( x1 * y2 - x2 * y1 ) + ( y1 - y2 ) * coords.at(1) + ( x2 - x1 ) * coords.at(2) ) / 2. / area;


    for ( int i = 1; i <= 3; i++ ) {
        if ( answer.at(i) < ( 0. - POINT_TOL ) ) {
            return 0;
        }

        if ( answer.at(i) > ( 1. + POINT_TOL ) ) {
            return 0;
        }
    }

    return 1;
}
} // end namespace oofem
