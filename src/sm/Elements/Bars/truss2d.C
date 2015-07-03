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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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

#include "Elements/Bars/truss2d.h"
#include "CrossSections/structuralcrosssection.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"
#include "fei1dlin.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(Truss2d);

FEI1dLin Truss2d :: interp(1);  

Truss2d :: Truss2d(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain)
{
    numberOfDofMans     = 2;
    length              = 0.;
    pitch               = 10.;   // a dummy value
    cs_mode             = 0;
}


void
Truss2d :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
    // eps = B*a = [-cos -sin cos sin]/l *[u1 v1 u2 v2]^t
    // cos = (x2 - x1) / l
    // sin = (z2 - z1) / l

    double l, x1, x2, z1, z2;
    //determine in which plane the truss is defined
    int c1 = 0, c2 = 0;
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

    l = this->computeLength();
    answer.times(1.0 / l / l);
}


void
Truss2d :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    // Will be the same as the regular B-matrix
    this->computeBmatrixAt(gp, answer);
}



void Truss2d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 2) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], 1, this);
    }
}


void
Truss2d :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver. This expression is
// valid in both local and global axes.
{
    answer.resize(4, 4);
    answer.zero();
    if ( !isActivated(tStep) ) {
        return;
    }

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    double density = this->giveStructuralCrossSection()->give('d', gp);
    double halfMass = density * this->giveCrossSection()->give(CS_Area, gp) * this->computeLength() * 0.5;
    answer.at(1, 1) = halfMass;
    answer.at(2, 2) = halfMass;
    answer.at(3, 3) = halfMass;
    answer.at(4, 4) = halfMass;
}


void
Truss2d :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp.
{
    double ksi, n1, n2;

    ksi = iLocCoord.at(1);
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
    int c1 = 0, c2 = 0;
    resolveCoordIndices(c1, c2);

    double ksi, n1, n2;

    ksi = lcoords.at(1);
    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;

    answer.resize(3);
    answer.zero();
    answer.at(c1) = n1 * this->giveNode(1)->giveCoordinate(c1) + n2 *this->giveNode(2)->giveCoordinate(c1);
    answer.at(c2) = n1 * this->giveNode(1)->giveCoordinate(c2) + n2 *this->giveNode(2)->giveCoordinate(c2);

    return 1;
}


double Truss2d :: computeVolumeAround(GaussPoint *gp)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    double weight  = gp->giveWeight();
    return 0.5 * this->computeLength() * weight * this->giveCrossSection()->give(CS_Area, gp);
}


double Truss2d :: computeLength()
// Returns the length of the receiver.
{
    //determine in which plane the truss is defined
    int c1 = 0, c2 = 0;
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
    int c1 = 0, c2 = 0;
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
        OOFEM_ERROR("Unknow cs_mode");
    }
}

IRResultType
Truss2d :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                            // Required by IR_GIVE_FIELD macro

    cs_mode = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, cs_mode, _IFT_Truss2d_cs);

    if ( cs_mode != 0 && cs_mode != 1 && cs_mode != 2 ) {
        OOFEM_WARNING("Unsupported value of cs_mode");
        return IRRT_BAD_FORMAT;
    }

    return NLStructuralElement :: initializeFrom(ir);
}


void
Truss2d :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( cs_mode == 0 ) {
        answer = {D_u, D_w};
    } else if ( cs_mode == 1 ) {
        answer = {D_u, D_v};
    } else if ( cs_mode == 2 ) {
        answer = {D_v, D_w};
    }
}


void
Truss2d :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveRealStress_1d(answer, gp, strain, tStep);
}


void
Truss2d :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveStiffnessMatrix_1d(answer, rMode, gp, tStep);
}


void
Truss2d :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    this->computeGlobalCoordinates( answer, gp->giveNaturalCoordinates() );
}

void
Truss2d :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
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

    this->computeNmatrixAt(gp->giveSubPatchCoordinates(), answer);
}


void
Truss2d :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    if ( iEdge != 1 ) {
        OOFEM_ERROR("wrong edge number");
    }


    answer.resize(4);
    answer.at(1) = 1;
    answer.at(2) = 2;
    answer.at(3) = 3;
    answer.at(4) = 4;
}

double
Truss2d ::   computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    if ( iEdge != 1 ) { // edge between nodes 1 2
        OOFEM_ERROR("wrong egde number");
    }

    double weight  = gp->giveWeight();
    return 0.5 * this->computeLength() * weight;
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


FEInterpolation *Truss2d :: giveInterpolation() const 
{
    // return interpolator
    return & interp; 
}


#ifdef __OOFEG
void Truss2d :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
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


void Truss2d :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    int c1, c2;
    resolveCoordIndices(c1, c2);

    GraphicObj *go;
    double defScale = gc.getDefScale();
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* point */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);

    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, defScale);

    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, defScale);

    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}
#endif
} // end namespace oofem
