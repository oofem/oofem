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

#include "trplanestrain.h"
#include "node.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "mathfem.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
#endif

namespace oofem {
FEI2dTrLin TrPlaneStrain :: interp(1, 2);

TrPlaneStrain :: TrPlaneStrain(int n, Domain *aDomain) :
    StructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(), NodalAveragingRecoveryModelInterface(),
    SPRNodalRecoveryModelInterface(), SpatialLocalizerInterface(),
    DirectErrorIndicatorRCInterface(),
    EIPrimaryUnknownMapperInterface(), ZZErrorEstimatorInterface(), ZZRemeshingCriteriaInterface(),
    MMAShapeFunctProjectionInterface(), HuertaErrorEstimatorInterface(), HuertaRemeshingCriteriaInterface()
    // Constructor.
{
    numberOfDofMans  = 3;
    area = -1;
    numberOfGaussPoints = 1;
}


Interface *
TrPlaneStrain :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return ( ZZNodalRecoveryModelInterface * ) this;
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return ( NodalAveragingRecoveryModelInterface * ) this;
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return ( SPRNodalRecoveryModelInterface * ) this;
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return ( SpatialLocalizerInterface * ) this;
    } else if ( interface == DirectErrorIndicatorRCInterfaceType ) {
        return ( DirectErrorIndicatorRCInterface * ) this;
    } else if ( interface == EIPrimaryUnknownMapperInterfaceType ) {
        return ( EIPrimaryUnknownMapperInterface * ) this;
    } else if ( interface == ZZErrorEstimatorInterfaceType ) {
        return ( ZZErrorEstimatorInterface * ) this;
    } else if ( interface == ZZRemeshingCriteriaInterfaceType ) {
        return ( ZZRemeshingCriteriaInterface * ) this;
    } else if ( interface == MMAShapeFunctProjectionInterfaceType ) {
        return ( MMAShapeFunctProjectionInterface * ) this;
    } else if ( interface == HuertaErrorEstimatorInterfaceType ) {
        return ( HuertaErrorEstimatorInterface * ) this;
    } else if ( interface == HuertaRemeshingCriteriaInterfaceType ) {
        return ( HuertaRemeshingCriteriaInterface * ) this;
    }

    return NULL;
}


void
TrPlaneStrain :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,
                                  int li, int ui)
// Returns the [4x6] strain-displacement matrix {B} of the receiver, eva-
// luated at gp.
{
    FloatMatrix dN;
    this->interp.evaldNdx( dN, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(4, 6);

    answer.at(1, 1) = dN.at(1, 1);
    answer.at(1, 3) = dN.at(2, 1);
    answer.at(1, 5) = dN.at(3, 1);

    answer.at(2, 2) = dN.at(1, 2);
    answer.at(2, 4) = dN.at(2, 2);
    answer.at(2, 6) = dN.at(3, 2);

    answer.at(4, 1) = dN.at(1, 2);
    answer.at(4, 2) = dN.at(1, 1);
    answer.at(4, 3) = dN.at(2, 2);
    answer.at(4, 4) = dN.at(2, 1);
    answer.at(4, 5) = dN.at(3, 2);
    answer.at(4, 6) = dN.at(3, 1);
}


void TrPlaneStrain :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Triangle, numberOfGaussPoints, _PlaneStrain);
    }
}


void
TrPlaneStrain :: computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp.
{
    FloatArray n;
    this->interp.evalN( n, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(2, 6);
    answer.zero();
    answer.at(1, 1) = answer.at(2, 2) = n.at(1);
    answer.at(1, 3) = answer.at(2, 4) = n.at(2);
    answer.at(1, 5) = answer.at(2, 6) = n.at(3);
}


void
TrPlaneStrain :: computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    FloatArray n;
    this->interp.edgeEvalN( n, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    answer.resize(2, 4);
    answer.at(1, 1) = answer.at(2, 2) = n.at(1);
    answer.at(1, 3) = answer.at(2, 4) = n.at(2);
}


void
TrPlaneStrain :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    answer.resize(4);
    if ( iEdge == 1 ) { // edge between nodes 1,2
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 3;
        answer.at(4) = 4;
    } else if ( iEdge == 2 ) { // edge between nodes 2 3
        answer.at(1) = 3;
        answer.at(2) = 4;
        answer.at(3) = 5;
        answer.at(4) = 6;
    } else if ( iEdge == 3 ) { // edge between nodes 3 1
        answer.at(1) = 5;
        answer.at(2) = 6;
        answer.at(3) = 1;
        answer.at(4) = 2;
    } else {
        _error("giveEdgeDofMapping: wrong edge number");
    }
}


double
TrPlaneStrain :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double detJ = this->interp.edgeGiveTransformationJacobian( iEdge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    return detJ *gp->giveWeight();
}


void
TrPlaneStrain :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    this->interp.edgeLocal2global( answer, iEdge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
}


int
TrPlaneStrain :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    // returns transformation matrix from
    // edge local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)
    //

    ///@todo Move this to interpolation class (similar code is repeated in several elements)
    double dx, dy, length;
    Node *nodeA, *nodeB;
    int aNode = 0, bNode = 0;

    answer.resize(2, 2);
    answer.zero();

    if ( iEdge == 1 ) { // edge between nodes 1 2
        aNode = 1;
        bNode = 2;
    } else if ( iEdge == 2 ) { // edge between nodes 2 3
        aNode = 2;
        bNode = 3;
    } else if ( iEdge == 3 ) { // edge between nodes 2 3
        aNode = 3;
        bNode = 1;
    } else {
        _error("computeEdgeVolumeAround: wrong egde number");
    }

    nodeA = this->giveNode(aNode);
    nodeB = this->giveNode(bNode);

    dx = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
    dy = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
    length = sqrt(dx * dx + dy * dy);

    answer.at(1, 1) = dx / length;
    answer.at(1, 2) = -dy / length;
    answer.at(2, 1) = dy / length;
    answer.at(2, 2) = dx / length;

    return 1;
}


double TrPlaneStrain :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double detJ, weight;

    weight = gp->giveWeight();
    detJ = fabs( this->interp.giveTransformationJacobian( * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ) );

    return detJ *weight *this->giveCrossSection()->give(CS_Thickness);
}

IRResultType
TrPlaneStrain :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    this->StructuralElement :: initializeFrom(ir);
    numberOfGaussPoints = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_TrPlaneStrain_nip, "nip"); // Macro
    if ( numberOfGaussPoints != 1 ) {
        numberOfGaussPoints = 1;
    }

    // set - up Gaussian integration points
    this->computeGaussPoints();

    return IRRT_OK;
}


double
TrPlaneStrain :: giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane)
//
// returns receivers characteristic length in gp (for some material models)
// for crack formed in plane with normal normalToCrackPlane.
//
{
    return this->giveLenghtInDir(normalToCrackPlane);
}


void
TrPlaneStrain :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(2, D_u, D_v);
}


int
TrPlaneStrain :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) ) {
        return 4;
    }

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}


void
TrPlaneStrain :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                            InternalStateType type, TimeStep *tStep)
{
    GaussPoint *gp;
    gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    this->giveIPValue(answer, gp, type, tStep);
}


void
TrPlaneStrain :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                           InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}


void
TrPlaneStrain :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(3);
    pap.at(1) = this->giveNode(1)->giveNumber();
    pap.at(2) = this->giveNode(2)->giveNumber();
    pap.at(3) = this->giveNode(3)->giveNumber();
}


void
TrPlaneStrain :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    answer.resize(1);
    if ( ( pap == this->giveNode(1)->giveNumber() ) ||
        ( pap == this->giveNode(2)->giveNumber() ) ||
        ( pap == this->giveNode(3)->giveNumber() ) ) {
        answer.at(1) = pap;
    } else {
        _error("SPRNodalRecoveryMI_giveDofMansDeterminedByPatch: node unknown");
    }
}


int
TrPlaneStrain :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return 1;
}


SPRPatchType
TrPlaneStrain :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dxy;
}


int
TrPlaneStrain :: SpatialLocalizerI_containsPoint(const FloatArray &coords)
{
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}


double
TrPlaneStrain :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
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


double
TrPlaneStrain :: DirectErrorIndicatorRCI_giveCharacteristicSize()
{
    return sqrt(this->computeArea() * 2.0);
}


int
TrPlaneStrain :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode,
                                                                  TimeStep *stepN, const FloatArray &coords,
                                                                  FloatArray &answer)
{
    FloatArray lcoords, u, nv;
    FloatMatrix n;
    int result;

    result = this->computeLocalCoordinates(lcoords, coords);

    this->interp.evalN( nv, lcoords, FEIElementGeometryWrapper(this) );

    n.resize(2, 6);
    n.zero();
    n.at(1, 1) = n.at(2, 2) = nv.at(1);
    n.at(1, 3) = n.at(2, 4) = nv.at(2);
    n.at(1, 5) = n.at(2, 6) = nv.at(3);

    this->computeVectorOf(EID_MomentumBalance, mode, stepN, u);

    answer.beProductOf(n, u);

    return result;
}


void
TrPlaneStrain :: EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer)
{
    giveDofManDofIDMask(1, EID_MomentumBalance, answer);
}


void
TrPlaneStrain :: MMAShapeFunctProjectionInterface_interpolateIntVarAt(FloatArray &answer, FloatArray &coords,
                                                                      coordType ct, nodalValContainerType &list,
                                                                      InternalStateType type, TimeStep *tStep)
{
    FloatArray n, lcoords;
    int vals;

    if ( ct == MMAShapeFunctProjectionInterface :: coordType_local ) {
        lcoords = coords;
    } else {
        computeLocalCoordinates(lcoords, coords);
    }

    this->interp.evalN( n, lcoords, FEIElementGeometryWrapper(this) );

    vals = list.at(1)->giveSize();
    answer.resize(vals);

    for ( int i = 1; i <= vals; i++ ) {
        answer.at(i) = n.at(1) * list.at(1)->at(i) + n.at(2) * list.at(2)->at(i) + n.at(3) * list.at(3)->at(i);
    }
}


void
TrPlaneStrain :: HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId,
                                                                  IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                                  HuertaErrorEstimatorInterface :: SetupMode sMode, TimeStep *tStep,
                                                                  int &localNodeId, int &localElemId, int &localBcId,
                                                                  IntArray &controlNode, IntArray &controlDof,
                                                                  HuertaErrorEstimator :: AnalysisMode aMode)
{
    Element *element = this->HuertaErrorEstimatorI_giveElement();
    int inode, nodes = 3, iside, sides = 3, nd1, nd2;
    FloatArray *corner [ 3 ], midSide [ 3 ], midNode, cor [ 3 ];
    double x = 0.0, y = 0.0;

    static int sideNode [ 3 ] [ 2 ] = { { 1, 2 }, { 2, 3 }, { 3, 1 } };

    if ( sMode == HuertaErrorEstimatorInterface :: NodeMode ||
        ( sMode == HuertaErrorEstimatorInterface :: BCMode && aMode == HuertaErrorEstimator :: HEE_linear ) ) {
        for ( inode = 0; inode < nodes; inode++ ) {
            corner [ inode ] = element->giveNode(inode + 1)->giveCoordinates();
            if ( corner [ inode ]->giveSize() != 3 ) {
                cor [ inode ].resize(3);
                cor [ inode ].at(1) = corner [ inode ]->at(1);
                cor [ inode ].at(2) = corner [ inode ]->at(2);
                cor [ inode ].at(3) = 0.0;

                corner [ inode ] = & ( cor [ inode ] );
            }

            x += corner [ inode ]->at(1);
            y += corner [ inode ]->at(2);
        }

        for ( iside = 0; iside < sides; iside++ ) {
            midSide [ iside ].resize(3);

            nd1 = sideNode [ iside ] [ 0 ] - 1;
            nd2 = sideNode [ iside ] [ 1 ] - 1;

            midSide [ iside ].at(1) = ( corner [ nd1 ]->at(1) + corner [ nd2 ]->at(1) ) / 2.0;
            midSide [ iside ].at(2) = ( corner [ nd1 ]->at(2) + corner [ nd2 ]->at(2) ) / 2.0;
            midSide [ iside ].at(3) = 0.0;
        }

        midNode.resize(3);

        midNode.at(1) = x / nodes;
        midNode.at(2) = y / nodes;
        midNode.at(3) = 0.0;
    }

    this->setupRefinedElementProblem2D(element, refinedElement, level, nodeId, localNodeIdArray, globalNodeIdArray,
                                       sMode, tStep, nodes, corner, midSide, midNode,
                                       localNodeId, localElemId, localBcId,
                                       controlNode, controlDof, aMode, "Quad1PlaneStrain");
}



#ifdef __OOFEG
 #define TR_LENGHT_REDUCT 0.3333

void TrPlaneStrain :: drawRawGeometry(oofegGraphicContext &gc)
{
    WCRec p [ 3 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 2 ].z = 0.;

    go =  CreateTriangle3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void TrPlaneStrain :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
{
    WCRec p [ 3 ];
    GraphicObj *go;
    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    double defScale = gc.getDefScale();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 2 ].z = 0.;

    go =  CreateTriangle3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


void TrPlaneStrain :: drawScalar(oofegGraphicContext &context)
{
    int i, indx, result = 0;
    WCRec p [ 3 ];
    GraphicObj *tr;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    FloatArray v1, v2, v3;
    double s [ 3 ], defScale;
    IntArray map;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    if ( context.giveIntVarMode() == ISM_recovered ) {
        result += this->giveInternalStateAtNode(v1, context.giveIntVarType(), context.giveIntVarMode(), 1, tStep);
        result += this->giveInternalStateAtNode(v2, context.giveIntVarType(), context.giveIntVarMode(), 2, tStep);
        result += this->giveInternalStateAtNode(v3, context.giveIntVarType(), context.giveIntVarMode(), 3, tStep);
    } else if ( context.giveIntVarMode() == ISM_local ) {
        GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
        result += giveIPValue(v1, gp, context.giveIntVarType(), tStep);
        v2 = v1;
        v3 = v1;
        result *= 3;
    }

    if ( result != 3 ) {
        return;
    }

    result = this->giveIntVarCompFullIndx( map, context.giveIntVarType() );

    if ( ( !result ) || ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
        return;
    }

    s [ 0 ] = v1.at(indx);
    s [ 1 ] = v2.at(indx);
    s [ 2 ] = v3.at(indx);

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);

    if ( context.getScalarAlgo() == SA_ISO_SURF ) {
        for ( i = 0; i < 3; i++ ) {
            if ( context.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = context.getDefScale();
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
                p [ i ].z = 0.;
            } else {
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                p [ i ].z = 0.;
            }
        }

        //EASValsSetColor(gc.getYieldPlotColor(ratio));
        context.updateFringeTableMinMax(s, 3);
        tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
        EGWithMaskChangeAttributes(LAYER_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    } else if ( ( context.getScalarAlgo() == SA_ZPROFILE ) || ( context.getScalarAlgo() == SA_COLORZPROFILE ) ) {
        double landScale = context.getLandScale();

        for ( i = 0; i < 3; i++ ) {
            if ( context.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = context.getDefScale();
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
                p [ i ].z = s [ i ] * landScale;
            } else {
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                p [ i ].z = s [ i ] * landScale;
            }
        }

        if ( context.getScalarAlgo() == SA_ZPROFILE ) {
            EASValsSetColor( context.getDeformedElementColor() );
            EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
            EASValsSetFillStyle(FILL_SOLID);
            tr =  CreateTriangle3D(p);
            EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | FILL_MASK | LAYER_MASK, tr);
        } else {
            context.updateFringeTableMinMax(s, 3);
            EASValsSetFillStyle(FILL_SOLID);
            tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
            EGWithMaskChangeAttributes(FILL_MASK | LAYER_MASK, tr);
        }

        EMAddGraphicsToModel(ESIModel(), tr);
    }
}

/*
 * void TrPlaneStrain :: drawInternalState (oofegGraphicContext& gc)
 * //
 * // Draws internal state graphics representation
 * //
 * {
 * WCRec p[3],l[2];
 * GraphicObj *tr;
 * StructuralMaterial *mat = (StructuralMaterial*) this->giveMaterial();
 * GaussPoint* gp;
 * double v1,v2,v3,ratio;
 * int i, nPlastGp;
 * DrawMode mode = gc.getDrawMode();
 * TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
 * double defScale = gc.getDefScale();
 *
 * if (!gc.testElementGraphicActivity(this)) return;
 *
 * // check for yield mode
 * if (mode == yieldState) {
 * // loop over available GPs
 * nPlastGp = 0;
 * for (i=1 ; i<= numberOfGaussPoints ; i++) {
 *  gp = integrationRulesArray[0]-> getIntegrationPoint(i-1) ;
 *  nPlastGp += (mat->giveStatusCharFlag(gp,ms_yield_flag) != 0);
 * }
 * if (nPlastGp == 0) return;
 * // nPlastGp should contain number of yielded gp in element
 * // good way should be select color accordingly
 * ratio = nPlastGp / numberOfGaussPoints;
 * EASValsSetLayer(OOFEG_YIELD_PATTERN_LAYER);
 * for (i=0; i< 3; i++) {
 * if (gc.getInternalVarsDefGeoFlag()) {
 *  // use deformed geometry
 *  p[i].x = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
 *  p[i].y = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
 *  p[i].z = 0.;
 *
 * } else {
 *  p[i].x = (FPNum) this->giveNode(i+1)->giveCoordinate(1);
 *  p[i].y = (FPNum) this->giveNode(i+1)->giveCoordinate(2);
 *  p[i].z = 0.;
 * }
 * }
 *
 * EASValsSetColor(gc.getYieldPlotColor(ratio));
 * tr =  CreateTriangle3D(p);
 * EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | FILL_MASK | LAYER_MASK, tr);
 * EMAddGraphicsToModel(ESIModel(), tr);
 *
 * } else if (mode == crackedState) {
 * // ask if any active crack exist
 * int crackStatus;
 * double ax,ay,bx,by,norm,xc,yc,length;
 * FloatMatrix crackDir;
 *
 * for (i=1 ; i<= numberOfGaussPoints ; i++) {
 *  gp = integrationRulesArray[0]-> getIntegrationPoint(i-1);
 *
 *  if (mat->giveStatusCharFlag (gp,ms_isCracked_flag) == 0) return;
 *  if (mat->GiveStatusCharMtrx(crackDir, gp,ms_crackDirection_flag)) {
 *  for (i = 1; i <= 3; i++) {
 *   crackStatus = (int) mat->giveStatusCharValue(gp, ms_crackStatus_flag,i);
 *   if ((crackStatus != pscm_NONE) && (crackStatus != pscm_CLOSED)) {
 *    // draw a crack
 *    // this element is 2d element in x-y plane
 *    //
 *    // compute perpendicular line to normal in xy plane
 *    ax = crackDir.at(1,i);
 *    ay = crackDir.at(2,i);
 *    if (fabs(ax) > 1.e-6) {
 *     by = 1.;
 *     bx = -ay*by/ax;
 *     norm = sqrt(bx*bx+by*by);
 *     bx = bx/norm;   // normalize to obtain unit vector
 *     by = by/norm;
 *    } else { by = 0.0; bx = 1.0;}
 *    // obtain gp global coordinates - here only one exists
 *    // it is in centre of gravity.
 *    xc = yc = 0.;
 *    for (i=0; i< 3; i++) {
 *     if (gc.getInternalVarsDefGeoFlag()) {
 *      // use deformed geometry
 *      xc += (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
 *      yc += (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
 *     } else {
 *      xc += (FPNum) this->giveNode(i+1)->giveCoordinate(1);
 *      yc += (FPNum) this->giveNode(i+1)->giveCoordinate(2);
 *     }
 *    }
 *    xc = xc / 3.; yc = yc/3.;
 *    length = TR_LENGHT_REDUCT * sqrt(2*this->giveArea()) / 2.0;
 *    l[0].x = (FPNum) xc+bx*length;
 *    l[0].y = (FPNum) yc+by*length;
 *    l[0].z = 0.;
 *    l[1].x = (FPNum) xc-bx*length;
 *    l[1].y = (FPNum) yc-by*length;
 *    l[1].z = 0.;
 *    EASValsSetLayer(OOFEG_CRACK_PATTERN_LAYER);
 *    EASValsSetLineWidth(OOFEG_CRACK_PATTERN_WIDTH);
 *    if ((crackStatus == pscm_SOFTENING)||(crackStatus == pscm_OPEN))
 *     EASValsSetColor(gc.getActiveCrackColor());
 *    else
 *     EASValsSetColor(gc.getCrackPatternColor());
 *    tr = CreateLine3D (l);
 *    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
 *    EMAddGraphicsToModel(ESIModel(), tr);
 *   }
 *  }
 *  //delete crackDir;
 * }
 * }
 * }
 * // check for valid stress-strain mode
 * if (!((mode == sxForce) || (mode == syForce) || (mode == szForce) || (mode == sxyForce))) return ;
 *
 * EASValsSetLayer(OOFEG_STRESS_CONTOUR_LAYER);
 *
 * for (i=0; i< 3; i++) {
 * if (gc.getInternalVarsDefGeoFlag()) {
 *  // use deformed geometry
 *  p[i].x = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
 *  p[i].y = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
 *  p[i].z = 0.;
 *
 * } else {
 *  p[i].x = (FPNum) this->giveNode(i+1)->giveCoordinate(1);
 *  p[i].y = (FPNum) this->giveNode(i+1)->giveCoordinate(2);
 *  p[i].z = 0.;
 * }
 * }
 *
 * int result = 0;
 * result+= this->giveInternalStateAtNode (gc, 1, &v1);
 * result+= this->giveInternalStateAtNode (gc, 2, &v2);
 * result+= this->giveInternalStateAtNode (gc, 3, &v3);
 *
 * if (result == 3) {
 * tr = CreateTriangleWD3D (p,v1,v2,v3);
 * EGWithMaskChangeAttributes(LAYER_MASK, tr);
 * EMAddGraphicsToModel(ESIModel(), tr);
 * }
 * }
 */

#endif
} // end namespace oofem
