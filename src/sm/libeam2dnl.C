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

#include "libeam2dnl.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "structuralms.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "flotarry.h"
#include "mathfem.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
LIBeam2dNL :: LIBeam2dNL(int n, Domain *aDomain) : NLStructuralElement(n, aDomain), LayeredCrossSectionInterface()
{
    numberOfDofMans     = 2;
    length              = 0.;
    pitch               = 10.;   // a dummy value
}

Interface *
LIBeam2dNL :: giveInterface(InterfaceType interface)
{
    if ( interface == LayeredCrossSectionInterfaceType ) {
        return ( LayeredCrossSectionInterface * ) this;
    }

    return NULL;
}


void
LIBeam2dNL :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the strain matrix of the receiver.
{
    double l, ksi, n1, n2, n1x, n2x, n3x;

    l    = this->giveLength();
    ksi  = aGaussPoint->giveCoordinate(1);

    n1    = 0.5 * ( 1.0 - ksi );
    n2    = 0.5 * ( 1.0 + ksi );
    n1x   = -1.0 / l;
    n2x   =  1.0 / l;
    n3x   = -4.0 * ksi / l;

    answer.resize(3, 6);
    answer.zero();

    // eps_x
    answer.at(1, 1) =  n1x;
    answer.at(1, 4) =  n2x;
    // kappa
    answer.at(2, 3) =  n1x;
    answer.at(2, 6) =  n2x;
    // gamma_xz
    answer.at(3, 2) =  n1x;
    answer.at(3, 3) =  n1 - n3x * l / 8.0;
    answer.at(3, 5) =  n2x;
    answer.at(3, 6) =  n2 + n3x * l / 8.0;
}


void
LIBeam2dNL :: computeNLBMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint, int i)
{
    //
    // Returns nonlinear part of geometrical equations of the receiver at gp.
    //
    // Returns A matrix (see Bittnar & Sejnoha Num. Met. Mech. part II, chap 9)
    double l, l8, ll88, ksi, n1x, n2x, n3x;

    l    = this->giveLength();
    ksi  = aGaussPoint->giveCoordinate(1);

    //n1    = 0.5*(1.0 - ksi);
    //n2    = 0.5*(1.0 + ksi);
    n1x   = -1.0 / l;
    n2x   =  1.0 / l;
    n3x   = -4.0 * ksi / l;
    l8    = l / 8.;
    ll88  = l8 * l8;

    answer.resize(6, 6);
    answer.zero();
    /*
     * if (i==1){
     * answer.at(1,1)= n1x*n1x;
     * answer.at(1,4)= n1x*n2x;
     * answer.at(2,2)= n1x*n1x;
     * answer.at(2,3)=-n1x*n3x*l8;
     * answer.at(2,5)= n1x*n2x;
     * answer.at(2,6)= n1x*n3x*l8;
     * answer.at(3,2)=-n3x*l8*n1x;
     * answer.at(3,3)= n3x*n3x*ll88;
     * answer.at(3,5)=-n3x*l8*n2x;
     * answer.at(3,6)=-n3x*n3x*ll88;
     *
     * answer.at(4,1)= n2x*n1x;
     * answer.at(4,4)= n2x*n2x;
     * answer.at(5,2)= n2x*n1x;
     * answer.at(5,3)=-n2x*n3x*l8;
     * answer.at(5,5)= n2x*n2x;
     * answer.at(5,6)= n2x*n3x*l8;
     * answer.at(6,2)= n3x*l8*n1x;
     * answer.at(6,3)=-n3x*n3x*ll88;
     * answer.at(6,5)= n3x*l8*n2x;
     * answer.at(6,2)= n3x*n3x*ll88;
     *
     * }
     *
     * if (i==2){
     * answer.at(1,2)=n1x*n1x;
     * answer.at(1,5)=n1x*n2x;
     * answer.at(2,1)=n1x*n1x;
     * answer.at(2,4)=n1x*n2x;
     * answer.at(4,2)=n2x*n1x;
     * answer.at(4,5)=n2x*n2x;
     * answer.at(5,1)=n2x*n1x;
     * answer.at(5,4)=n2x*n2x;
     *
     * }
     *
     * if (i==3){
     * answer.at(1,3)=n1x*n1;
     * answer.at(1,6)=n1x*n2;
     * answer.at(3,1)=n1*n1x;
     * answer.at(3,4)=n1*n2x;
     * answer.at(4,3)=n2x*n1;
     * answer.at(4,6)=n2x*n2;
     * answer.at(6,1)=n2*n1x;
     * answer.at(6,4)=n2*n2x;
     * }
     */

    /* Working
     *
     * // large deflection
     * // based on von Karman equation
     *
     * if (i==1) {
     * answer.at(2,2)= n1x*n1x;
     * answer.at(2,3)=-n1x*n3x*l8;
     * answer.at(2,5)= n1x*n2x;
     * answer.at(2,6)= n1x*n3x*l8;
     * answer.at(3,2)=-n3x*l8*n1x;
     * answer.at(3,3)= n3x*n3x*ll88;
     * answer.at(3,5)=-n3x*l8*n2x;
     * answer.at(3,6)=-n3x*n3x*ll88;
     *
     * answer.at(5,2)= n2x*n1x;
     * answer.at(5,3)=-n2x*n3x*l8;
     * answer.at(5,5)= n2x*n2x;
     * answer.at(5,6)= n2x*n3x*l8;
     * answer.at(6,2)= n3x*l8*n1x;
     * answer.at(6,3)=-n3x*n3x*ll88;
     * answer.at(6,5)= n3x*l8*n2x;
     * answer.at(6,2)= n3x*n3x*ll88;
     *
     * }
     */

    if ( i == 1 ) {
        //  answer.at(1,1)= n1x*n1x;
        //  answer.at(1,4)= n1x*n2x;
        answer.at(2, 2) = n1x * n1x;
        answer.at(2, 3) = -n1x * n3x * l8;
        answer.at(2, 5) = n1x * n2x;
        answer.at(2, 6) = n1x * n3x * l8;
        answer.at(3, 2) = -n3x * l8 * n1x;
        answer.at(3, 3) = n3x * n3x * ll88;
        answer.at(3, 5) = -n3x * l8 * n2x;
        answer.at(3, 6) = -n3x * n3x * ll88;

        //  answer.at(4,1)= n2x*n1x;
        //  answer.at(4,4)= n2x*n2x;
        answer.at(5, 2) = n2x * n1x;
        answer.at(5, 3) = -n2x * n3x * l8;
        answer.at(5, 5) = n2x * n2x;
        answer.at(5, 6) = n2x * n3x * l8;
        answer.at(6, 2) = n3x * l8 * n1x;
        answer.at(6, 3) = -n3x * n3x * ll88;
        answer.at(6, 5) = n3x * l8 * n2x;
        answer.at(6, 2) = n3x * n3x * ll88;
    }

    /*
     * if (i==2){
     * answer.at(1,2)=n1x*n1x;
     * answer.at(1,5)=n1x*n2x;
     * answer.at(2,1)=n1x*n1x;
     * answer.at(2,4)=n1x*n2x;
     * answer.at(4,2)=n2x*n1x;
     * answer.at(4,5)=n2x*n2x;
     * answer.at(5,1)=n2x*n1x;
     * answer.at(5,4)=n2x*n2x;
     *
     * }
     *
     * if (i==3){
     * answer.at(1,3)=n1x*n1;
     * answer.at(1,6)=n1x*n2;
     * answer.at(3,1)=n1*n1x;
     * answer.at(3,4)=n1*n2x;
     * answer.at(4,3)=n2x*n1;
     * answer.at(4,6)=n2x*n2;
     * answer.at(6,1)=n2*n1x;
     * answer.at(6,4)=n2*n2x;
     * }
     */
}


void
LIBeam2dNL :: computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    int i, j, n;
    double dV;
    GaussPoint *gp;
    IntegrationRule *iRule;
    FloatArray stress;
    FloatMatrix A;
    Material *mat = this->giveMaterial();

    answer.resize(6, 6);
    answer.zero();

    iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    // assemble initial stress matrix
    for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        dV = this->computeVolumeAround(gp);
        stress = ( ( StructuralMaterialStatus * ) mat->giveStatus(gp) )->giveStressVector();
        n = stress.giveSize();
        if ( n ) {
            for ( j = 1; j <= n; j++ ) {
                // loop over each component of strain vector
                this->computeNLBMatrixAt(A, gp, j);
                if ( A.isNotEmpty() ) {
                    A.times(stress.at(j) * dV);
                    answer.add(A);
                }
            }
        }
    }
}


void LIBeam2dNL :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 2);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Line, 1, _2dBeam);
    }
}


void
LIBeam2dNL :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver. This expression is
// valid in both local and global axes.
{
    Material *mat;
    double halfMass;
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    mat        = this->giveMaterial();
    halfMass   = mat->give('d', gp) * this->giveCrossSection()->give(CS_Area) * this->giveLength() / 2.;
    answer.resize(6, 6);
    answer.zero();
    answer.at(1, 1) = halfMass;
    answer.at(2, 2) = halfMass;
    answer.at(4, 4) = halfMass;
    answer.at(5, 5) = halfMass;
}


void
LIBeam2dNL :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
{
    double ksi, n1, n2, n3;
    double l = this->giveLength();

    ksi = aGaussPoint->giveCoordinate(1);
    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;
    n3  = ( 1. - ksi * ksi );

    answer.resize(3, 6);
    answer.zero();

    // us
    answer.at(1, 1) = n1;
    answer.at(1, 4) = n2;
    // w
    answer.at(2, 2) = n1;
    answer.at(2, 3) = -n3 * l / 8.;
    answer.at(2, 5) = n2;
    answer.at(2, 6) = n3 * l / 8.;
    // fi_y
    answer.at(3, 3) = n1;
    answer.at(3, 6) = n2;
}


bool
LIBeam2dNL :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    double sine, cosine;

    sine = sin( this->givePitch() );
    cosine = cos(pitch);

    answer.resize(6, 6);
    answer.zero();

    answer.at(1, 1) =  cosine;
    answer.at(1, 2) =  sine;
    answer.at(2, 1) = -sine;
    answer.at(2, 2) =  cosine;
    answer.at(3, 3) =  1.;
    answer.at(4, 4) =  cosine;
    answer.at(4, 5) =  sine;
    answer.at(5, 4) = -sine;
    answer.at(5, 5) =  cosine;
    answer.at(6, 6) =  1.;

    return true;
}


double LIBeam2dNL :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    double weight  = aGaussPoint->giveWeight();
    return weight * 0.5 * this->giveLength();
}


void
LIBeam2dNL :: computeStrainVectorInLayer(FloatArray &answer, GaussPoint *masterGp,
                                         GaussPoint *slaveGp, TimeStep *tStep)
//
// returns full 3d strain vector of given layer (whose z-coordinate from center-line is
// stored in slaveGp) for given tStep
//
{
    FloatArray masterGpStrain;
    double layerZeta, layerZCoord, top, bottom;

    this->computeStrainVector(masterGpStrain, masterGp, tStep);
    top    = masterGp->giveElement()->giveCrossSection()->give(CS_TopZCoord);
    bottom = masterGp->giveElement()->giveCrossSection()->give(CS_BottomZCoord);
    layerZeta = slaveGp->giveCoordinate(3);
    layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );

    answer.resize(6); // {Exx,Eyy,Ezz,GMyz,GMzx,GMxy}
    answer.zero();

    answer.at(1) = masterGpStrain.at(1) + masterGpStrain.at(2) * layerZCoord;
    answer.at(5) = masterGpStrain.at(3);
}


void
LIBeam2dNL :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(3, D_u, D_w, R_v);
}


int
LIBeam2dNL :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    double ksi, n1, n2;

    ksi = lcoords.at(1);
    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;

    answer.resize(3);
    answer.at(1) = n1 * this->giveNode(1)->giveCoordinate(1) + n2 *this->giveNode(2)->giveCoordinate(1);
    answer.at(3) = n1 * this->giveNode(1)->giveCoordinate(3) + n2 *this->giveNode(2)->giveCoordinate(3);

    return 1;
}


double LIBeam2dNL :: giveLength()
// Returns the length of the receiver.
{
    double dx, dy;
    Node *nodeA, *nodeB;

    if ( length == 0. ) {
        nodeA   = this->giveNode(1);
        nodeB   = this->giveNode(2);
        dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
        dy      = nodeB->giveCoordinate(3) - nodeA->giveCoordinate(3);
        length  = sqrt(dx * dx + dy * dy);
    }

    return length;
}


double LIBeam2dNL :: givePitch()
// Returns the pitch of the receiver.
{
    double xA, xB, yA, yB;
    Node *nodeA, *nodeB;

    if ( pitch == 10. ) {             // 10. : dummy initialization value
        nodeA  = this->giveNode(1);
        nodeB  = this->giveNode(2);
        xA     = nodeA->giveCoordinate(1);
        xB     = nodeB->giveCoordinate(1);
        yA     = nodeA->giveCoordinate(3);
        yB     = nodeB->giveCoordinate(3);
        pitch  = atan2(yB - yA, xB - xA);
    }

    return pitch;
}


IRResultType
LIBeam2dNL :: initializeFrom(InputRecord *ir)
{
    this->NLStructuralElement :: initializeFrom(ir);
    this->computeGaussPoints();
    return IRRT_OK;
}


void
LIBeam2dNL :: computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint)
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
LIBeam2dNL :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
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
}

double
LIBeam2dNL ::   computeEdgeVolumeAround(GaussPoint *aGaussPoint, int iEdge)
{
    if ( iEdge != 1 ) { // edge between nodes 1 2
        _error("computeEdgeVolumeAround: wrong egde number");
    }

    double weight  = aGaussPoint->giveWeight();
    return 0.5 * this->giveLength() * weight;
}

int
LIBeam2dNL :: computeLoadGToLRotationMtrx(FloatMatrix &answer)
{
    /*
     * Returns transformation matrix from global coordinate system to local
     * element coordinate system for element load vector components.
     * If no transformation is necessary, answer is empty matrix (default);
     */

    // f(elemLocal) = T * f (global)
    double sine, cosine;

    answer.resize(3, 3);
    answer.zero();

    sine           = sin( this->givePitch() );
    cosine         = cos(pitch);

    answer.at(1, 1) = cosine;
    answer.at(1, 2) = -sine;
    answer.at(2, 1) = sine;
    answer.at(2, 2) = cosine;
    answer.at(3, 3) = 1.0;

    return 1;
}

int
LIBeam2dNL :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    // returns transformation matrix from
    // edge local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)
    //
    answer.beEmptyMtrx();
    return 0;
}


void
LIBeam2dNL :: computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode)
{
    StructuralElement::computeBodyLoadVectorAt(answer, load, tStep, mode);
    answer.times(this->giveCrossSection()->give(CS_Area));
}


#ifdef __OOFEG
void LIBeam2dNL :: drawRawGeometry(oofegGraphicContext &gc)
{
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    //  if (!go) { // create new one
    WCRec p [ 2 ];   /* poin */
    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = 0.;
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = 0.;
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void LIBeam2dNL :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
{
    GraphicObj *go;
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    double defScale = gc.getDefScale();
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* poin */
    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].y = 0.;
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);

    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].y = 0.;
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}
#endif
} // end namespace oofem
