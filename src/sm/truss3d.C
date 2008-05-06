/* $Header: /home/cvs/bp/oofem/sm/src/truss3d.C,v 1.5 2003/04/06 14:08:32 bp Exp $ */
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

//   file Truss3d.C

#include "truss3d.h"
#include "domain.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"

#include "engngm.h"
#ifndef __MAKEDEPEND
#include <stdlib.h>
#include <math.h>
#endif

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#endif


//
// upravit teplota u geom. nelinearity
//

Truss3d :: Truss3d(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain)
    // Constructor.
{
    numberOfDofMans     = 2;
    rotationMatrix      = NULL;
    length              = 0.;
    //  referenceNode       - 0  ;
}

Interface *
Truss3d :: giveInterface(InterfaceType interface)
{
    if ( interface == DirectErrorIndicatorRCInterfaceType ) {
        return ( DirectErrorIndicatorRCInterface * ) this;
    }

    return NULL;
}

void
Truss3d :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
    double coeff, l, x1, x2, y1, y2, z1, z2;
    // FloatMatrix* answer;

    x1 = this->giveNode(1)->giveCoordinate(1);
    y1 = this->giveNode(1)->giveCoordinate(2);
    z1 = this->giveNode(1)->giveCoordinate(3);
    x2 = this->giveNode(2)->giveCoordinate(1);
    y2 = this->giveNode(2)->giveCoordinate(2);
    z2 = this->giveNode(2)->giveCoordinate(3);

    // answer = new FloatMatrix(1,4);
    answer.resize(1, 6);

    answer.at(1, 1) = x1 - x2;
    answer.at(1, 2) = y1 - y2;
    answer.at(1, 3) = z1 - z2;
    answer.at(1, 4) = x2 - x1;
    answer.at(1, 5) = y2 - y1;
    answer.at(1, 6) = z2 - z1;

    l = this->giveLength();
    coeff = 1.0 / l / l;
    answer.times(coeff);

    return;
}

void
Truss3d :: computeNLBMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint, int i)
//
// Returns nonlinear part of geometrical equations of the receiver at gp.
//
// Returns A matrix (see Bittnar & Sejnoha Num. Met. Mech. part II, chap 9)
{
    double coeff, l;

    l = this->giveLength();
    coeff = 1.0 / l / l;

    answer.resize(6, 6);
    answer.zero();

    answer.at(1, 1) = coeff;
    answer.at(1, 4) = coeff * ( -1 );
    answer.at(2, 2) = coeff;
    answer.at(2, 5) = coeff * ( -1 );
    answer.at(3, 3) = coeff;
    answer.at(3, 6) = coeff * ( -1 );
    answer.at(4, 1) = coeff * ( -1 );
    answer.at(4, 4) = coeff;
    answer.at(5, 2) = coeff * ( -1 );
    answer.at(5, 5) = coeff;
    answer.at(6, 3) = coeff * ( -1 );
    answer.at(6, 6) = coeff;

    return;
}



void Truss3d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    numberOfIntegrationRules = 1;
    integrationRulesArray = new IntegrationRule * [ 1 ];
    integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 2);
    integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Line, 1, _1dMat);
}



void
Truss3d :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver. This expression is
// valid in both local and global axes.
{
    Material *mat;
    double halfMass;

    mat        = this->giveMaterial();
    halfMass   = mat->give('d') * this->giveCrossSection()->give('A') * this->giveLength() / 2.;
    answer.resize(6, 6);
    answer.zero();
    answer.at(1, 1) = halfMass;
    answer.at(2, 2) = halfMass;
    answer.at(3, 3) = halfMass;
    answer.at(4, 4) = halfMass;
    answer.at(5, 5) = halfMass;
    answer.at(6, 6) = halfMass;

    if ( this->updateRotationMatrix() ) {
        answer.rotatedWith(* this->rotationMatrix);
    }

    return;
}


void
Truss3d :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
{
    double ksi, n1, n2;
    //FloatMatrix* answer ;

    ksi = aGaussPoint->giveCoordinate(1);
    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;
    //answer = new FloatMatrix(2,4) ;
    answer.resize(3, 6);
    answer.zero();

    answer.at(1, 1) = n1;
    answer.at(1, 4) = n2;
    answer.at(2, 2) = n1;
    answer.at(2, 5) = n2;
    answer.at(3, 3) = n1;
    answer.at(3, 6) = n2;

    return;
}

int
Truss3d :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    double ksi, n1, n2;

    ksi = lcoords.at(1);
    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;

    answer.resize(3);
    answer.at(1) = n1 * this->giveNode(1)->giveCoordinate(1) + n2 * this->giveNode(2)->giveCoordinate(1);
    answer.at(2) = n1 * this->giveNode(1)->giveCoordinate(2) + n2 * this->giveNode(2)->giveCoordinate(2);
    answer.at(3) = n1 * this->giveNode(1)->giveCoordinate(3) + n2 * this->giveNode(2)->giveCoordinate(3);

    return 1;
}


double Truss3d :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    double weight  = aGaussPoint->giveWeight();
    return 0.5 * this->giveLength() * weight * this->giveCrossSection()->give('A');
}


double Truss3d :: giveLength()
// Returns the length of the receiver.
{
    double dx, dy, dz;
    Node *nodeA, *nodeB;

    if ( length == 0. ) {
        nodeA   = this->giveNode(1);
        nodeB   = this->giveNode(2);
        dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
        dy      = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
        dz      = nodeB->giveCoordinate(3) - nodeA->giveCoordinate(3);
        length  = sqrt(dx * dx + dy * dy + dz * dz);
    }

    return length;
}


int
Truss3d :: giveLocalCoordinateSystem(FloatMatrix &answer)
//
// returns a unit vectors of local coordinate system at element
// stored rowwise (mainly used by some materials with ortho and anisotrophy)
//
{
    FloatArray lx(3), ly(3), lz(3), help(3);
    double length = this->giveLength();
    Node *nodeA, *nodeB;
    int i;

    // if (referenceNode == 0)
    // _error ("instanciateFrom: wrong reference node specified");

    answer.resize(3, 3);
    answer.zero();
    nodeA  = this->giveNode(1);
    nodeB  = this->giveNode(2);
    // refNode= this->giveDomain()->giveNode (this->referenceNode);

    for ( i = 1; i <= 3; i++ ) {
        lx.at(i) = ( nodeB->giveCoordinate(i) - nodeA->giveCoordinate(i) ) / length;
        // help.at(i)=(refNode->giveCoordinate(i)-nodeA->giveCoordinate(i));
    }

    int minIndx = 1;
    for ( i = 2; i <= 3; i++ ) {
        if ( lx.at(i) < lx.at(minIndx) ) {
            minIndx = i;
        }
    }

    help.zero();
    help.at(minIndx) = 1.0;

    lz.beVectorProductOf(lx, help);
    lz.normalize();
    ly.beVectorProductOf(lz, lx);
    ly.normalize();

    for ( i = 1; i <= 3; i++ ) {
        answer.at(1, i) = lx.at(i);
        answer.at(2, i) = ly.at(i);
        answer.at(3, i) = lz.at(i);
    }

    return 1;
}

IRResultType
Truss3d :: initializeFrom(InputRecord *ir)
{
    this->NLStructuralElement :: initializeFrom(ir);
    this->computeGaussPoints();
    return IRRT_OK;
}


void
Truss3d ::   giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const {
    // returns DofId mask array for inode element node.
    // DofId mask array determines the dof ordering requsted from node.
    // DofId mask array contains the DofID constants (defined in cltypes.h)
    // describing physical meaning of particular DOFs.
    //IntArray* answer = new IntArray (2);
    answer.resize(3);

    answer.at(1) = D_u;
    answer.at(2) = D_v;
    answer.at(3) = D_w;

    return;
}



void
Truss3d :: computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint)
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

    this->computeNmatrixAt(aGaussPoint, answer);
    return;
}


void
Truss3d :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    if ( iEdge != 1 ) {
        _error("giveEdgeDofMapping: wrong edge number");
    }


    answer.resize(6);
    answer.at(1) = 1;
    answer.at(2) = 2;
    answer.at(3) = 3;
    answer.at(4) = 4;
    answer.at(5) = 5;
    answer.at(6) = 6;

    return;
}

double
Truss3d ::   computeEdgeVolumeAround(GaussPoint *aGaussPoint, int iEdge)
{
    if ( iEdge != 1 ) { // edge between nodes 1 2
        _error("computeEdgeVolumeAround: wrong egde number");
    }

    double weight  = aGaussPoint->giveWeight();
    return 0.5 * this->giveLength() * weight;
}

int
Truss3d :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    // returns transformation matrix from
    // edge local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)
    //
    FloatMatrix lcs;
    this->giveLocalCoordinateSystem(lcs);
    answer.beTranspositionOf(lcs);

    return 1;
}


#ifdef __OOFEG
void Truss3d :: drawRawGeometry(oofegGraphicContext &gc)
{
    GraphicObj *go;
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* poin */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void Truss3d :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
{
    GraphicObj *go;
    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    double defScale = gc.getDefScale();
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* poin */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);

    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}
#endif



