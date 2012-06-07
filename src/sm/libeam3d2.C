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

#include "libeam3d2.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "structuralcrosssection.h"
#include "structuralmaterial.h"
#include "structuralms.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "flotarry.h"
#include "dof.h"
#include "timestep.h"
#include "contextioerr.h"
#include "util.h"

#ifndef __MAKEDEPEND
 #include <math.h>
#endif

#ifdef __OOFEG
 #include "engngm.h"
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
 #include "conTable.h"
#endif

namespace oofem {
LIBeam3d2 :: LIBeam3d2(int n, Domain *aDomain) : NLStructuralElement(n, aDomain), tc(), tempTc()
{
    numberOfDofMans     = 2;
    referenceNode       = 0;
    length              = 0.;

    tempTcCounter       = 0;
}


void
LIBeam3d2 :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the strain matrix of the receiver.
// eeps = {\eps_x, \gamma_xz, \gamma_xy, \der{phi_x}{x}, \kappa_y, \kappa_z}^T
{
    double l, ksi, n1, n2, n1x, n2x;

    if ( this->nlGeometry ) {
        l = this->giveCurrentLength( domain->giveEngngModel()->giveCurrentStep() );
    } else {
        l = this->giveLength();
    }

    ksi   = aGaussPoint->giveCoordinate(1);

    answer.resize(6, 12);
    answer.zero();

    n1    = 0.5 * ( 1 - ksi );
    n2    = 0.5 * ( 1. + ksi );
    n1x   = -1.0 / l;
    n2x   =  1.0 / l;

    answer.at(1, 1) =  -1. / l;
    answer.at(1, 7) =   1. / l;

    answer.at(2, 3)  = n1x;
    answer.at(2, 5)  = n1;
    answer.at(2, 9)  = n2x;
    answer.at(2, 11) = n2;

    answer.at(3, 2)  = n1x;
    answer.at(3, 6)  = -n1;
    answer.at(3, 8)  = n2x;
    answer.at(3, 12) = -n2;

    answer.at(4, 4)  =  -1. / l;
    answer.at(4, 10) =   1. / l;

    answer.at(5, 5)  = n1x;
    answer.at(5, 11) = n2x;

    answer.at(6, 6)  = n1x;
    answer.at(6, 12) = n2x;
}

void
LIBeam3d2 :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 2);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Line, 1, _3dBeam);
    }
}



void
LIBeam3d2 :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver. This expression is
// valid in both local and global axes.
{
    Material *mat;
    double halfMass;
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    mat        = this->giveMaterial();
    halfMass   = mat->give('d', gp) * this->giveCrossSection()->give(CS_Area) * this->giveLength() / 2.;
    answer.resize(12, 12);
    answer.zero();
    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = halfMass;
    answer.at(7, 7) = answer.at(8, 8) = answer.at(9, 9) = halfMass;
}


void
LIBeam3d2 :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
{
    double ksi, n1, n2;

    ksi = aGaussPoint->giveCoordinate(1);
    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;

    answer.resize(6, 12);
    answer.zero();

    answer.at(1, 1) = n1;
    answer.at(1, 7) = n2;
    answer.at(2, 2) = n1;
    answer.at(2, 8) = n2;
    answer.at(3, 3) = n1;
    answer.at(3, 9) = n2;

    answer.at(4, 4) = n1;
    answer.at(4, 10) = n2;
    answer.at(5, 5) = n1;
    answer.at(5, 11) = n2;
    answer.at(6, 6) = n1;
    answer.at(6, 12) = n2;
}


void
LIBeam3d2 :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
// Returns the stiffness matrix of the receiver, expressed in the global
// axes.
{
    this->StructuralElement :: computeStiffnessMatrix(answer, rMode, tStep);
}


bool
LIBeam3d2 :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    int i, j;

    answer.resize(12, 12);
    answer.zero();

    if ( this->nlGeometry  ==  0 ) {
        FloatMatrix lcs;

        this->giveLocalCoordinateSystem(lcs);
        for ( i = 1; i <= 3; i++ ) {
            for ( j = 1; j <= 3; j++ ) {
                answer.at(i, j) = lcs.at(i, j);
                answer.at(i + 3, j + 3) = lcs.at(i, j);
                answer.at(i + 6, j + 6) = lcs.at(i, j);
                answer.at(i + 9, j + 9) = lcs.at(i, j);
            }
        }
    } else {
        this->updateTempTriad( domain->giveEngngModel()->giveCurrentStep() );

        for ( i = 1; i <= 3; i++ ) {
            for ( j = 1; j <= 3; j++ ) {
                answer.at(i, j) = tempTc.at(j, i);
                answer.at(i + 3, j + 3) = tempTc.at(j, i);
                answer.at(i + 6, j + 6) = tempTc.at(j, i);
                answer.at(i + 9, j + 9) = tempTc.at(j, i);
            }
        }
    }

    return 1;
}


double
LIBeam3d2 :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    double weight  = aGaussPoint->giveWeight();
    return weight * 0.5 * this->giveLength();
}


void
LIBeam3d2 :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(6, D_u, D_v, D_w, R_u, R_v, R_w);
}


int
LIBeam3d2 :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    double ksi, n1, n2;

    ksi = lcoords.at(1);
    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;

    answer.resize(3);
    answer.at(1) = n1 * this->giveNode(1)->giveCoordinate(1) + n2 *this->giveNode(2)->giveCoordinate(1);
    answer.at(2) = n1 * this->giveNode(1)->giveCoordinate(2) + n2 *this->giveNode(2)->giveCoordinate(2);
    answer.at(3) = n1 * this->giveNode(1)->giveCoordinate(3) + n2 *this->giveNode(2)->giveCoordinate(3);

    return 1;
}


int
LIBeam3d2 :: testElementExtension(ElementExtension ext)
{
    return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 );
}


double
LIBeam3d2 :: giveLength()
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


IRResultType
LIBeam3d2 :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // first call parent
    NLStructuralElement :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, referenceNode, IFT_LIBeam3d2_refnode, "refnode"); // Macro
    if ( referenceNode == 0 ) {
        _error("instanciateFrom: wrong reference node specified");
    }

    //  if (this->hasString (initString, "dofstocondense")) {
    //    dofsToCondense = this->ReadIntArray (initString, "dofstocondense");
    //    if (dofsToCondense->giveSize() >= 12)
    //      _error ("instanciateFrom: wrong input data for condensed dofs");
    //  } else {
    //    dofsToCondense = NULL;
    //  }
    // compute initial triad at centre - requires nodal coordinates
    FloatMatrix lcs;
    this->giveLocalCoordinateSystem(lcs);
    this->tc.beTranspositionOf(lcs);

    this->computeGaussPoints();

    return IRRT_OK;
}


void
LIBeam3d2 :: computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint)
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
LIBeam3d2 :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */
    int i;

    if ( iEdge != 1 ) {
        _error("giveEdgeDofMapping: wrong edge number");
    }


    answer.resize(12);
    for ( i = 1; i <= 12; i++ ) {
        answer.at(i) = i;
    }
}


double
LIBeam3d2 :: computeEdgeVolumeAround(GaussPoint *aGaussPoint, int iEdge)
{
    if ( iEdge != 1 ) { // edge between nodes 1 2
        _error("computeEdgeVolumeAround: wrong egde number");
    }

    double weight  = aGaussPoint->giveWeight();
    return 0.5 * this->giveLength() * weight;
}


int
LIBeam3d2 :: computeLoadGToLRotationMtrx(FloatMatrix &answer)
{
    /*
     * Returns transformation matrix from global coordinate system to local
     * element coordinate system for element load vector components.
     * If no transformation is necessary, answer is empty matrix (default);
     */

    // f(elemLocal) = T * f (global)

    FloatMatrix lcs;
    int i, j;

    answer.resize(6, 6);
    answer.zero();

    this->giveLocalCoordinateSystem(lcs);
    for ( i = 1; i <= 3; i++ ) {
        for ( j = 1; j <= 3; j++ ) {
            answer.at(i, j) = lcs.at(i, j);
            answer.at(3 + i, 3 + j) = lcs.at(i, j);
        }
    }
    return 1;
}


void
LIBeam3d2 :: computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode)
{
    NLStructuralElement::computeBodyLoadVectorAt(answer, load, tStep, mode);
    answer.times(this->giveCrossSection()->give(CS_Area));
}


int
LIBeam3d2 :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
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


int
LIBeam3d2 :: giveLocalCoordinateSystem(FloatMatrix &answer)
//
// returns a unit vectors of local coordinate system at element
// stored rowwise (mainly used by some materials with ortho and anisotrophy)
//
{
    FloatArray lx(3), ly(3), lz(3), help(3);
    double length = this->giveLength();
    Node *nodeA, *nodeB, *refNode;
    int i;

    answer.resize(3, 3);
    answer.zero();
    nodeA  = this->giveNode(1);
    nodeB  = this->giveNode(2);
    refNode = this->giveDomain()->giveNode(this->referenceNode);

    for ( i = 1; i <= 3; i++ ) {
        lx.at(i) = ( nodeB->giveCoordinate(i) - nodeA->giveCoordinate(i) ) / length;
        help.at(i) = ( refNode->giveCoordinate(i) - nodeA->giveCoordinate(i) );
    }

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


void
LIBeam3d2 :: updateTempTriad(TimeStep *tStep)
{
    // test if not previously done
    if ( tStep->giveSolutionStateCounter() == tempTcCounter ) {
        return;
    }

    FloatArray u, centreSpin(3);
    FloatMatrix dR(3, 3);

    // ask element's displacement increments
    this->computeVectorOf(EID_MomentumBalance, VM_Incremental, tStep, u);

    // interpolate spin at the centre
    centreSpin.at(1) = 0.5 * ( u.at(4) + u.at(10) );
    centreSpin.at(2) = 0.5 * ( u.at(5) + u.at(11) );
    centreSpin.at(3) = 0.5 * ( u.at(6) + u.at(12) );

    // compute rotation matrix from centreSpin pseudovector
    this->computeRotMtrx(dR, centreSpin);

    // update triad
    tempTc.beProductOf(dR, tc);
    // remember timestamp
    tempTcCounter = tStep->giveSolutionStateCounter();
}


void
LIBeam3d2 :: computeRotMtrx(FloatMatrix &answer, FloatArray &psi)
{
    FloatMatrix S(3, 3), SS(3, 3);
    double psiSize;

    if ( psi.giveSize() != 3 ) {
        _error("computeSMtrx: psi param size mismatch");
    }

    answer.resize(3, 3);
    answer.zero();

    psiSize = psi.computeNorm();
    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = 1.;

    if ( psiSize <= 1.e-40 ) {
        return;
    }

    this->computeSMtrx(S, psi);
    SS.beProductOf(S, S);
    S.times(sin(psiSize) / psiSize);
    SS.times( ( 1. - cos(psiSize) ) / ( psiSize * psiSize ) );

    answer.add(S);
    answer.add(SS);
}


void
LIBeam3d2 :: computeSMtrx(FloatMatrix &answer, FloatArray &vec)
{
    if ( vec.giveSize() != 3 ) {
        _error("computeSMtrx: vec param size mismatch");
    }

    answer.resize(3, 3);

    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = 0.;
    answer.at(1, 2) = -vec.at(3);
    answer.at(1, 3) =  vec.at(2);
    answer.at(2, 1) =  vec.at(3);
    answer.at(2, 3) = -vec.at(1);
    answer.at(3, 1) = -vec.at(2);
    answer.at(3, 2) =  vec.at(1);
}


void
LIBeam3d2 :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray ui, PrevEpsilon;
    FloatMatrix b;

    this->computeVectorOf(EID_MomentumBalance, VM_Incremental, tStep, ui);

    this->computeBmatrixAt(gp, b);
    // increment of strains
    answer.beProductOf(b, ui);

    if ( this->nlGeometry ) {
        // add increment to previous total value
        PrevEpsilon = ( ( StructuralMaterialStatus * ) this->giveMaterial()->giveStatus(gp) )->giveStrainVector();
        if ( PrevEpsilon.giveSize() ) {
            answer.add(PrevEpsilon);
        }
    }
}


double
LIBeam3d2 :: giveCurrentLength(TimeStep *tStep)
// Returns the length of the receiver.
{
    double dx, dy, dz;
    FloatArray u;

    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);

    dx = ( this->giveNode(2)->giveCoordinate(1) + u.at(7) ) -
         ( this->giveNode(1)->giveCoordinate(1) + u.at(1) );
    dy = ( this->giveNode(2)->giveCoordinate(2) + u.at(8) ) -
         ( this->giveNode(1)->giveCoordinate(2) + u.at(2) );
    dz = ( this->giveNode(2)->giveCoordinate(3) + u.at(9) ) -
         ( this->giveNode(1)->giveCoordinate(3) + u.at(3) );

    return sqrt(dx * dx + dy * dy + dz * dz);
}


void LIBeam3d2 :: updateYourself(TimeStep *tStep)
// Updates the receiver at end of step.
{
    NLStructuralElement :: updateYourself(tStep);

    // update triad
    this->updateTempTriad(tStep);
    tc = tempTc;

    // update curvature
    //FloatArray curv;
    //this->computeTempCurv (curv, tStep);
    //kappa = curv;
}


void
LIBeam3d2 :: initForNewStep()
// initializes receiver to new time step or can be used
// if current time step must be restarted
//
// call material->initGpForNewStep() for all GPs.
//
{
    NLStructuralElement :: initForNewStep();
    tempTc = tc;
}



contextIOResultType LIBeam3d2 :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;
    if ( ( iores = NLStructuralElement :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = tc.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}



contextIOResultType LIBeam3d2 :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;
    if ( ( iores = NLStructuralElement :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = tc.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


void
LIBeam3d2 :: FiberedCrossSectionInterface_computeStrainVectorInFiber(FloatArray &answer, GaussPoint *masterGp,
                                                                     GaussPoint *slaveGp, TimeStep *tStep)
{
    FloatArray masterGpStrain;
    double layerYCoord, layerZCoord;

    this->computeStrainVector(masterGpStrain, masterGp, tStep);
    layerZCoord = slaveGp->giveCoordinate(2);
    layerYCoord = slaveGp->giveCoordinate(1);

    answer.resize(6); // {Exx,Eyy,Ezz,GMyz,GMzx,GMxy}
    answer.zero();

    answer.at(1) = masterGpStrain.at(1) + masterGpStrain.at(5) * layerZCoord - masterGpStrain.at(6) * layerYCoord;
    answer.at(5) = masterGpStrain.at(2);
    answer.at(6) = masterGpStrain.at(3);
}

Interface *
LIBeam3d2 :: giveInterface(InterfaceType interface)
{
    if ( interface == FiberedCrossSectionInterfaceType ) {
        return ( FiberedCrossSectionInterface * ) this;
    }

    return NULL;
}



#ifdef __OOFEG
void
LIBeam3d2 :: drawRawGeometry(oofegGraphicContext &gc)
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


void
LIBeam3d2 :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
{
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    double defScale = gc.getDefScale();
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* poin */
    char const *colors[] = {
        "red", "green", "blue"
    };

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

    // draw centre triad
    int i, succ;
    double coeff = this->giveLength() / 3.;
    //this->updateTempTriad (tStep);

    p [ 0 ].x = 0.5 * ( ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale) + ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale) );
    p [ 0 ].y = 0.5 * ( ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale) + ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale) );
    p [ 0 ].z = 0.5 * ( ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale) + ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale) );

    // draw t1
    for ( i = 1; i <= 3; i++ ) {
        p [ 1 ].x = p [ 0 ].x + coeff *tc.at(1, i);
        p [ 1 ].y = p [ 0 ].y + coeff *tc.at(2, i);
        p [ 1 ].z = p [ 0 ].z + coeff *tc.at(3, i);

        EASValsSetColor( ColorGetPixelFromString(oofem_tmpstr(colors [ i - 1 ]), & succ) );

        go = CreateLine3D(p);
        EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
        EMAddGraphicsToModel(ESIModel(), go);
    }
}


void
LIBeam3d2 :: drawScalar(oofegGraphicContext &context)
{
    WCRec p [ 2 ];
    GraphicObj *go;
    FloatArray v;
    IntArray map;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    double defScale;
    int indx;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    if ( context.getInternalVarsDefGeoFlag() ) {
        // use deformed geometry
        defScale = context.getDefScale();
        p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
        p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
        p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);
        p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
        p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
        p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);
    } else {
        p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
        p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
        p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
        p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
        p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
        p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);
    }

    GaussPoint *gp;
    double s;
    gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    if ( giveIPValue(v, gp, context.giveIntVarType(), tStep) == 0 ) {
        return;
    }

    this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
    if ( ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
        return;
    }

    s = v.at(1);
    context.updateFringeTableMinMax(& s, 1);

    go = CreateLine3D(p);
    EASValsSetColor( context.getElementColor() );
    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);

    EASValsSetMType(FILLED_CIRCLE_MARKER);
    go = CreateMarkerWD3D( p, v.at(1) );
    EGWithMaskChangeAttributes(LAYER_MASK | FILL_MASK | MTYPE_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}

#endif
} // end namespace oofem
