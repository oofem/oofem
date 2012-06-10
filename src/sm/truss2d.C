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

#include "truss2d.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "mathfem.h"

#ifdef __OOFEG
 #include "engngm.h"
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {

Truss2d :: Truss2d(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain)
{
    numberOfDofMans     = 2;
    length              = 0.;
    pitch               = 10.;   // a dummy value
    cs_mode             = 0;
}


void
Truss2d :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
    double l, x1, x2, z1, z2;
    //determine in which plane the truss is defined
    int c1, c2;
    resolveCoordIndices(c1, c2);

    x1 = this->giveNode(1)->giveCoordinate(c1);
    z1 = this->giveNode(1)->giveCoordinate(c2);
    x2 = this->giveNode(2)->giveCoordinate(c1);
    z2 = this->giveNode(2)->giveCoordinate(c2);

    answer.resize(1, 4);

    answer.at(1, 1) = x1 - x2;
    answer.at(1, 2) = z1 - z2;
    answer.at(1, 3) = x2 - x1;
    answer.at(1, 4) = z2 - z1;

    l = this->giveLength();
    answer.times( 1.0 / l / l );
}

void
Truss2d :: computeNLBMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint, int i)
//
// Returns nonlinear part of geometrical equations of the receiver at gp.
//
// Returns A matrix (see Bittnar & Sejnoha Num. Met. Mech. part II, chap 9)
{
    double coeff, l;

    l = this->giveLength();
    coeff = 1.0 / l / l;

    answer.resize(4, 4);
    answer.zero();

    answer.at(1, 1) = coeff;
    answer.at(1, 3) = coeff * ( -1 );
    answer.at(2, 2) = coeff;
    answer.at(2, 4) = coeff * ( -1 );
    answer.at(3, 1) = coeff * ( -1 );
    answer.at(3, 3) = coeff;
    answer.at(4, 2) = coeff * ( -1 );
    answer.at(4, 4) = coeff;
}



void Truss2d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 2);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Line, 1, _1dMat);
    }
}


void
Truss2d :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver. This expression is
// valid in both local and global axes.
{
    Material *mat;
    double halfMass;
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    answer.resize(4, 4);
    answer.zero();
    if ( !isActivated(tStep) ) {
        return;
    }

    mat        = this->giveMaterial();
    halfMass   =  mat->give('d', gp) * this->giveCrossSection()->give(CS_Area) * this->giveLength() / 2.;
    answer.at(1, 1) = halfMass;
    answer.at(2, 2) = halfMass;
    answer.at(3, 3) = halfMass;
    answer.at(4, 4) = halfMass;
}


void
Truss2d :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
{
    double ksi, n1, n2;

    ksi = aGaussPoint->giveCoordinate(1);
    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;

    answer.resize(2, 4);
    answer.zero();
    answer.at(1, 1) = n1;
    answer.at(1, 3) = n2;
    answer.at(2, 2) = n1;
    answer.at(2, 4) = n2;
}


int
Truss2d :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    //determine in which plane the truss is defined
    int c1, c2;
    resolveCoordIndices(c1, c2);

    double ksi, n1, n2;

    ksi = lcoords.at(1);
    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;

    answer.resize(3);
    answer.at(c1) = n1 * this->giveNode(1)->giveCoordinate(c1) + n2 *this->giveNode(2)->giveCoordinate(c1);
    answer.at(c2) = n1 * this->giveNode(1)->giveCoordinate(c2) + n2 *this->giveNode(2)->giveCoordinate(c2);

    return 1;
}


double Truss2d :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    double weight  = aGaussPoint->giveWeight();
    return 0.5 * this->giveLength() * weight * this->giveCrossSection()->give(CS_Area);
}


double Truss2d :: giveLength()
// Returns the length of the receiver.
{
    //determine in which plane the truss is defined
    int c1, c2;
    resolveCoordIndices(c1, c2);

    double dx, dz;
    Node *nodeA, *nodeB;

    if ( length == 0. ) {
        nodeA   = this->giveNode(1);
        nodeB   = this->giveNode(2);
        dx      = nodeB->giveCoordinate(c1) - nodeA->giveCoordinate(c1);
        dz      = nodeB->giveCoordinate(c2) - nodeA->giveCoordinate(c2);
        length  = sqrt(dx * dx + dz * dz);
    }

    return length;
}


double Truss2d :: givePitch()
// Returns the pitch of the receiver.
{
    double xA, xB, zA, zB;
    Node *nodeA, *nodeB;
    //determine in which plane the truss is defined
    int c1, c2;
    resolveCoordIndices(c1, c2);

    if ( pitch == 10. ) {             // 10. : dummy initialization value
        nodeA  = this->giveNode(1);
        nodeB  = this->giveNode(2);
        xA     = nodeA->giveCoordinate(c1);
        xB     = nodeB->giveCoordinate(c1);
        zA     = nodeA->giveCoordinate(c2);
        zB     = nodeB->giveCoordinate(c2);
        pitch  = atan2(zB - zA, xB - xA);
    }

    return pitch;
}


int
Truss2d :: giveLocalCoordinateSystem(FloatMatrix &answer)
//
// returns a unit vectors of local coordinate system at element
// stored rowwise (mainly used by some materials with ortho and anisotrophy)
//
{
    double sine, cosine;

    answer.resize(3, 3);
    answer.zero();

    sine = sin( this->givePitch() );
    cosine = cos(pitch);

    answer.at(1, 1) = cosine;
    answer.at(1, 2) = sine;
    answer.at(2, 1) = -sine;
    answer.at(2, 2) = cosine;
    answer.at(3, 3) = 1.0;

    return 1;
}


void
Truss2d :: resolveCoordIndices(int &c1, int &c2)
{
    if ( cs_mode == 0 ) {
        //xz-plane
        c1 = 1;
        c2 = 3;
    } else if ( cs_mode == 1 ) {
        //xy-plane
        c1 = 1;
        c2 = 2;
    } else if ( cs_mode == 2 ) {
        //yz-plane
        c1 = 2;
        c2 = 3;
    } else {
        _error("resolveCoordIndices: Unknow cs_mode");
    }
}

IRResultType
Truss2d :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                            // Required by IR_GIVE_FIELD macro
    this->NLStructuralElement :: initializeFrom(ir);

    cs_mode = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, cs_mode, IFT_Truss2d_cs, "cs"); // Macro

    if ( cs_mode != 0 && cs_mode != 1 && cs_mode != 2 ) {
        _error("Unsupported value of cs_mode");
    }

    this->computeGaussPoints();
    return IRRT_OK;
}


void
Truss2d :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    if ( cs_mode == 0 ) {
        answer.setValues(2, D_u, D_w);
    } else if ( cs_mode == 1 ) {
        answer.setValues(2, D_u, D_v);
    } else if ( cs_mode == 2 ) {
        answer.setValues(2, D_v, D_w);
    }
}

void
Truss2d :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    computeGlobalCoordinates( answer, * ( gp->giveCoordinates() ) );
}

void
Truss2d :: computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint)
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
}


void
Truss2d :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    if ( iEdge != 1 ) {
        _error("giveEdgeDofMapping: wrong edge number");
    }


    answer.resize(4);
    answer.at(1) = 1;
    answer.at(2) = 2;
    answer.at(3) = 3;
    answer.at(4) = 4;
}

double
Truss2d ::   computeEdgeVolumeAround(GaussPoint *aGaussPoint, int iEdge)
{
    if ( iEdge != 1 ) { // edge between nodes 1 2
        _error("computeEdgeVolumeAround: wrong egde number");
    }

    double weight  = aGaussPoint->giveWeight();
    return 0.5 * this->giveLength() * weight;
}

int
Truss2d :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    // returns transformation matrix from
    // edge local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)
    //
    //double dx,dy, length ;
    double sine, cosine;

    answer.resize(2, 2);
    answer.zero();

    sine           = sin( this->givePitch() );
    cosine         = cos(pitch);

    answer.at(1, 1) = cosine;
    answer.at(1, 2) = -sine;
    answer.at(2, 1) = sine;
    answer.at(2, 2) = cosine;

    return 1;
}


#ifdef __OOFEG
void Truss2d :: drawRawGeometry(oofegGraphicContext &gc)
{
    int c1, c2;
    resolveCoordIndices(c1, c2);

    GraphicObj *go;
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* point */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    if ( cs_mode == 0 ) {
        p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(c1);
        p [ 0 ].y = 0.;
        p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(c2);
        p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(c1);
        p [ 1 ].y = 0.;
        p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(c2);
    } else if ( cs_mode == 1 ) {
        p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(c1);
        p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(c2);
        p [ 0 ].z = 0.;
        p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(c1);
        p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(c2);
        p [ 1 ].z = 0.;
    } else if ( cs_mode == 2 ) {
        p [ 0 ].x = 0.;
        p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(c1);
        p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(c2);
        p [ 1 ].x = 0.;
        p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(c1);
        p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(c2);
    }

    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void Truss2d :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
{
    int c1, c2;
    resolveCoordIndices(c1, c2);

    GraphicObj *go;
    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    double defScale = gc.getDefScale();
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* point */
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
} // end namespace oofem
