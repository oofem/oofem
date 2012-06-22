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

#include "axisymm3d.h"
#include "node.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "engngm.h"
#include "mathfem.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

namespace oofem {
FEI2dTrLin Axisymm3d :: interpolation(1, 2);

Axisymm3d :: Axisymm3d(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(), NodalAveragingRecoveryModelInterface(),
    SPRNodalRecoveryModelInterface()
    // Constructor.
{
    numberOfDofMans = 3;
    area = -1;
    numberOfGaussPoints = 1;
    numberOfFiAndShGaussPoints = 1;
}

Axisymm3d :: ~Axisymm3d()
// destructor
{ }

Interface *
Axisymm3d :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return ( ZZNodalRecoveryModelInterface * ) this;
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return ( NodalAveragingRecoveryModelInterface * ) this;
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return ( SPRNodalRecoveryModelInterface * ) this;
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return ( SpatialLocalizerInterface * ) this;
    }

    return NULL;
}


void
Axisymm3d :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver,
// evaluated at aGaussPoint.
{
    int i;
    FloatArray n(3);

    answer.resize(2, 6);
    answer.zero();
    this->interpolation.evalN(n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    for ( i = 1; i <= 3; i++ ) {
        answer.at(1, 2 * i - 1) = n.at(i);
        answer.at(2, 2 * i - 0) = n.at(i);
    }
}


void
Axisymm3d :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the [6x6] strain-displacement matrix {B} of the receiver, eva-
// luated at aGaussPoint.
{
    int i;
    double x, r;
    int size, ind = 1;
    FloatMatrix dnx;

    this->interpolation.evaldNdx(dnx, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));


    if ( ui == ALL_STRAINS ) {
        size = 6;
        ui = 6;
    } else {
        size = ui - li + 1;
    }

    if ( ( size < 0 ) || ( size > 6 ) ) {
        _error("ComputeBmatrixAt size mismatch");
    }

    answer.resize(size, 6);
    answer.zero();

    if ( ( li <= 1 ) && ( ui >= 1 ) ) {
        answer.at(ind, 1) = dnx.at(1, 1);
        answer.at(ind, 3) = dnx.at(2, 1);
        answer.at(ind, 5) = dnx.at(3, 1);
        ind++;
    }

    if ( ( li <= 2 ) && ( ui >= 2 ) ) {
        answer.at(ind, 2) = dnx.at(1, 2);
        answer.at(ind, 4) = dnx.at(2, 2);
        answer.at(ind, 6) = dnx.at(3, 2);
        ind++;
    }

    if ( ( li <= 3 ) && ( ui >= 3 ) ) {
        FloatArray n(4);
        this->interpolation.evalN(n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

        r = 0.;
        for ( i = 1; i <= numberOfDofMans; i++ ) {
            x  = this->giveNode(i)->giveCoordinate(1);
            r += x * n.at(i);
        }

        answer.at(ind, 1) = n.at(1) / r;
        answer.at(ind, 3) = n.at(2) / r;
        answer.at(ind, 5) = n.at(3) / r;
        ind++;
    }

    if ( ( li <= 4 ) && ( ui >= 4 ) ) {
        ind++;
    }

    if ( ( li <= 5 ) && ( ui >= 5 ) ) {
        ind++;
    }

    if ( ( li <= 6 ) && ( ui >= 6 ) ) {
        answer.at(ind, 1) = dnx.at(1, 2);
        answer.at(ind, 2) = dnx.at(1, 1);
        answer.at(ind, 3) = dnx.at(2, 2);
        answer.at(ind, 4) = dnx.at(2, 1);
        answer.at(ind, 5) = dnx.at(3, 2);
        answer.at(ind, 6) = dnx.at(3, 1);
        ind++;
    }
}


double
Axisymm3d :: giveArea()
// returns the area occupied by the receiver
{
    if ( area > 0 ) {
        return area;         // check if previously computed
    }
    area = fabs( this->interpolation.giveArea(FEIElementGeometryWrapper(this)) );
    return area;
}


double
Axisymm3d :: computeVolumeAround(GaussPoint *aGaussPoint)
{
    int i;
    double determinant, weight, volume, r, x;
    FloatArray n(4);

    this->interpolation.evalN(n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    r = 0.;
    for ( i = 1; i <= numberOfDofMans; i++ ) {
        x  = this->giveNode(i)->giveCoordinate(1);
        r += x * n.at(i);
    }

    determinant = fabs( this->interpolation.giveTransformationJacobian(* aGaussPoint->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this)) );

    weight = aGaussPoint->giveWeight();
    volume = determinant * weight * r;

    return volume;
}


void
Axisymm3d :: computeGaussPoints()
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 2;
        integrationRulesArray = new IntegrationRule * [ 2 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 2);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Triangle, numberOfGaussPoints, _3dMat);
        integrationRulesArray [ 1 ] = new GaussIntegrationRule(2, this, 3, 6);
        integrationRulesArray [ 1 ]->setUpIntegrationPoints(_Triangle, numberOfFiAndShGaussPoints, _3dMat);
    }
}


IRResultType
Axisymm3d :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    this->StructuralElement :: initializeFrom(ir);

    numberOfGaussPoints      = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_Axisymm3d_nip, "nip"); // Macro

    numberOfFiAndShGaussPoints    = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfFiAndShGaussPoints, IFT_Axisymm3d_nipfish, "nipfish"); // Macro

    if ( !( ( numberOfGaussPoints == 1 ) ||
           ( numberOfGaussPoints == 4 ) ||
           ( numberOfGaussPoints == 7 ) ) ) {
        numberOfGaussPoints = 1;
    }

    if ( !( ( numberOfFiAndShGaussPoints == 1 ) ||
           ( numberOfFiAndShGaussPoints == 4 ) ||
           ( numberOfFiAndShGaussPoints == 7 ) ) ) {
        numberOfFiAndShGaussPoints = 1;
    }

    this->computeGaussPoints();

    return IRRT_OK;
}


void
Axisymm3d :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step stepN. The nature of these strains depends
// on the element's type.
{
    int i;
    FloatMatrix b, A;
    FloatArray u, Epsilon, help;
    fMode mode = domain->giveEngngModel()->giveFormulation();

    answer.resize(6);
    answer.zero();

    if ( mode == TL ) {  // Total Lagrange formulation
        this->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, u);
        // linear part of strain tensor (in vector form)

        this->computeBmatrixAt(gp, b, 1, 2);
        Epsilon.beProductOf(b, u);
        answer.at(1) = Epsilon.at(1);
        answer.at(2) = Epsilon.at(2);
        // delete Epsilon; // delete b;

        if ( numberOfFiAndShGaussPoints == 1 ) {
            //
            // if reduced integration in one gp only
            // force the evaluation of eps_fi in this gauss point
            // instead of evaluating in given gp
            //
            GaussPoint *helpGaussPoint;
            helpGaussPoint = integrationRulesArray [ 1 ]->getIntegrationPoint(0);

            this->computeBmatrixAt(helpGaussPoint, b, 3, 6);
        } else {
            _error("ComputeStrainVector: numberOfRotGaussPoints size mismatch");
        }

        Epsilon.beProductOf(b, u);
        answer.at(3) = Epsilon.at(1);
        answer.at(6) = Epsilon.at(4);

        if ( nlGeometry ) {
            for ( i = 1; i <= 6; i++ ) {
                // nonlin part of strain vector
                this->computeNLBMatrixAt(A, gp, i);
                if ( A.isNotEmpty() ) {
                    help.beProductOf(A, u);
                    answer.at(i) += 0.5 * u.dotProduct(help);
                }
            }
        }
    } else if ( mode == AL ) { // actualized Lagrange formulation
        _error("computeStrainVector : unsupported mode");
    }
}



#define NONZERO_COORD_TOL 1.e-2

double
Axisymm3d :: giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane)
//
// returns receivers characteristic length in gp (for some material models)
// for crack formed in plane with normal normalToCrackPlane.
//
{
    // return this -> giveLenghtInDir(normalToCrackPlane) / sqrt (this->numberOfGaussPoints);
    //
    // we must handle special case - crack caused by hoop stres Sigma_z
    // note: z-axis is always principal axis, because there exist only one
    // nonzero shear strain (Sigma_xy), which will cause only rotation of
    // principal axises in x-y plane (around z-axis).
    // so this yields following :
    // normalCrackPlane(3) component  can be only zero
    // (then crack normal lies in x-y plane) or equal to 1.0 (crack normal
    // is perpendicular to to x--y plane) - crack caused by hoop strain.
    if ( fabs( normalToCrackPlane.at(3) ) > NONZERO_COORD_TOL ) {
        double r = 0.;
        int i;
        for ( i = 1; i <= numberOfDofMans; i++ ) {
            r += this->giveNode(i)->giveCoordinate(1);
        }

        r = r / ( ( double ) numberOfDofMans );
        return r;
    } else {
        return this->giveLenghtInDir(normalToCrackPlane);
    }
}


void
Axisymm3d :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    answer.setValues(2, D_u, D_v);
}


int
Axisymm3d :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) ) {
        return 6;
    }

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}


void
Axisymm3d :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                        InternalStateType type, TimeStep *tStep)
{
    GaussPoint *gp;

    if ( numberOfGaussPoints != 1 ) {
        answer.resize(0); // for more gp's need to be refined
        return;
    }

    gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    this->giveIPValue(answer, gp, type, tStep);
}

void
Axisymm3d :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                       InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}


void
Axisymm3d :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(3);
    pap.at(1) = this->giveNode(1)->giveNumber();
    pap.at(2) = this->giveNode(2)->giveNumber();
    pap.at(3) = this->giveNode(3)->giveNumber();
}

void
Axisymm3d :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
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
Axisymm3d :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return this->giveDefaultIntegrationRulePtr()->getNumberOfIntegrationPoints();
}


SPRPatchType
Axisymm3d :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dxy;
}


void
Axisymm3d :: computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint)
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

    FloatArray n(2);
    this->interpolation.edgeEvalN(n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(2, 4);
    answer.zero();

    answer.at(1, 1) = n.at(1);
    answer.at(1, 3) = n.at(2);
    answer.at(2, 2) = n.at(1);
    answer.at(2, 4) = n.at(2);
}


void
Axisymm3d :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
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
Axisymm3d ::   computeEdgeVolumeAround(GaussPoint *aGaussPoint, int iEdge)
{
    FloatArray c(2);
    this->computeEdgeIpGlobalCoords(c, aGaussPoint, iEdge);
    double result = this->interpolation.edgeGiveTransformationJacobian(iEdge, * aGaussPoint->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this));


    return c.at(1) * result * aGaussPoint->giveWeight();
}


void
Axisymm3d ::   computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    this->interpolation.edgeLocal2global(answer, iEdge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));
}


int
Axisymm3d :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    // returns transformation matrix from
    // edge local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)
    //
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

    nodeA   = this->giveNode(aNode);
    nodeB   = this->giveNode(bNode);

    dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
    dy      = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
    length = sqrt(dx * dx + dy * dy);

    answer.at(1, 1) = dx / length;
    answer.at(1, 2) = -dy / length;
    answer.at(2, 1) = dy / length;
    answer.at(2, 2) = dx / length;

    return 1;
}


int
Axisymm3d :: SpatialLocalizerI_containsPoint(const FloatArray &coords)
{
    FloatArray lcoords;
    return this->interpolation.global2local(lcoords, coords, FEIElementGeometryWrapper(this));
}


double
Axisymm3d :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
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


#ifdef __OOFEG
 #define TR_LENGHT_REDUCT 0.3333

void Axisymm3d :: drawRawGeometry(oofegGraphicContext &gc)
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


void Axisymm3d :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
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

void
Axisymm3d :: drawScalar(oofegGraphicContext &context)
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

    this->giveIntVarCompFullIndx( map, context.giveIntVarType() );

    if ( ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
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
 * void Axisymm3d :: drawInternalState (oofegGraphicContext &gc)
 * //
 * // Draws internal state graphics representation
 * //
 * {
 * WCRec p[3],l[2],q[4];
 * GraphicObj *tr;
 * StructuralMaterial *mat = (StructuralMaterial*) this->giveMaterial();
 * GaussPoint* gp;
 * double v1,v2,v3,ratio;
 * int i, nPlastGp;
 * IntegrationRule* iRule = integrationRulesArray[0];
 * TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
 * double defScale = gc.getDefScale();
 *
 * if (!gc.testElementGraphicActivity(this)) return;
 *
 * // check for yield mode
 * if (gc.getDrawMode() == yieldState) {
 * // loop over available GPs
 * nPlastGp = 0;
 * for (i=1 ; i<= iRule->getNumberOfIntegrationPoints() ; i++) {
 *  gp = iRule-> getIntegrationPoint(i) ;
 *  nPlastGp += (mat->giveStatusCharFlag(gp,ms_yield_flag) != 0);
 * }
 * if (nPlastGp == 0) return;
 * // nPlastGp should contain number of yielded gp in element
 * // good way should be select color accordingly
 * ratio = nPlastGp / numberOfGaussPoints;
 * EASValsSetLayer(OOFEG_YIELD_PATTERN_LAYER);
 * if (gc.getInternalVarsDefGeoFlag()) {
 * // use deformed geometry
 * p[0].x = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(1,tStep,DisplacementVector,defScale);
 * p[0].y = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(2,tStep,DisplacementVector,defScale);
 * p[0].z = 0.;
 * p[1].x = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(1,tStep,DisplacementVector,defScale);
 * p[1].y = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(2,tStep,DisplacementVector,defScale);
 * p[1].z = 0.;
 * p[2].x = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(1,tStep,DisplacementVector,defScale);
 * p[2].y = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(2,tStep,DisplacementVector,defScale);
 * p[2].z = 0.;
 *
 * } else {
 * p[0].x = (FPNum) this->giveNode(1)->giveCoordinate(1);
 * p[0].y = (FPNum) this->giveNode(1)->giveCoordinate(2);
 * p[0].z = 0.;
 * p[1].x = (FPNum) this->giveNode(2)->giveCoordinate(1);
 * p[1].y = (FPNum) this->giveNode(2)->giveCoordinate(2);
 * p[1].z = 0.;
 * p[2].x = (FPNum) this->giveNode(3)->giveCoordinate(1);
 * p[2].y = (FPNum) this->giveNode(3)->giveCoordinate(2);
 * p[2].z = 0.;
 * }
 *
 * EASValsSetColor(gc.getYieldPlotColor(ratio));
 * tr =  CreateTriangle3D(p);
 * EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | FILL_MASK | LAYER_MASK, tr);
 * EMAddGraphicsToModel(ESIModel(), tr);
 *
 * } else if (gc.getDrawMode() == crackedState) {
 * // ask if any active crack exist
 * int crackStatus;
 * double ax,ay,bx,by,norm,xc,yc,length;
 * FloatMatrix crackDir;
 *
 * if (numberOfGaussPoints != 1) return;
 * //   for (igp=1 ; igp<= numberOfGaussPoints ; igp++) {
 * {
 *  gp = iRule-> getIntegrationPoint(0);
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
 *    if (gc.getInternalVarsDefGeoFlag()) {
 *     // use deformed geometry
 *     xc = (this->giveNode(1)->giveUpdatedCoordinate(1,tStep,DisplacementVector,defScale)+
 *        this->giveNode(2)->giveUpdatedCoordinate(1,tStep,DisplacementVector,defScale)+
 *        this->giveNode(3)->giveUpdatedCoordinate(1,tStep,DisplacementVector,defScale)) / 3.0 ;
 *     yc = (this->giveNode(1)->giveUpdatedCoordinate(2,tStep,DisplacementVector,defScale) +
 *        this->giveNode(2)->giveUpdatedCoordinate(2,tStep,DisplacementVector,defScale) +
 *        this->giveNode(3)->giveUpdatedCoordinate(2,tStep,DisplacementVector,defScale)) / 3.0 ;
 *
 *    } else {
 *     xc = (this->giveNode(1)->giveCoordinate(1) +
 *        this->giveNode(2)->giveCoordinate(1) +
 *        this->giveNode(3)->giveCoordinate(1)) / 3.0 ;
 *     yc = (this->giveNode(1)->giveCoordinate(2) +
 *        this->giveNode(2)->giveCoordinate(2) +
 *        this->giveNode(3)->giveCoordinate(2)) / 3.0 ;
 *    }
 *    length = TR_LENGHT_REDUCT * sqrt(2*this->giveArea()) / 2.0;
 *    if (i!=3) {
 *     l[0].x = (FPNum) xc+bx*length;
 *     l[0].y = (FPNum) yc+by*length;
 *     l[0].z = 0.;
 *     l[1].x = (FPNum) xc-bx*length;
 *     l[1].y = (FPNum) yc-by*length;
 *     l[1].z = 0.;
 *     EASValsSetLayer(OOFEG_CRACK_PATTERN_LAYER);
 *     EASValsSetLineWidth(OOFEG_CRACK_PATTERN_WIDTH);
 *     if ((crackStatus == pscm_SOFTENING)||(crackStatus == pscm_OPEN))
 *      EASValsSetColor(gc.getActiveCrackColor());
 *     else
 *      EASValsSetColor(gc.getCrackPatternColor());
 *     tr = CreateLine3D (l);
 *     EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
 *     EMAddGraphicsToModel(ESIModel(), tr);
 *    } else {
 *     // special case - crack normal is z axis - draw square
 *     bx=by=1.0;
 *     q[0].x = (FPNum) xc+0.5*bx*length;
 *     q[0].y = (FPNum) yc+0.5*by*length;
 *     q[0].z = 0.;
 *     q[1].x = (FPNum) xc-0.5*bx*length;
 *     q[1].y = (FPNum) yc+0.5*by*length;
 *     q[1].z = 0.;
 *     q[2].x = (FPNum) xc-0.5*bx*length;
 *     q[2].y = (FPNum) yc-0.5*by*length;
 *     q[2].z = 0.;
 *     q[3].x = (FPNum) xc+0.5*bx*length;
 *     q[3].y = (FPNum) yc-0.5*by*length;
 *     q[3].z = 0.;
 *
 *     EASValsSetLayer(OOFEG_CRACK_PATTERN_LAYER);
 *     EASValsSetLineWidth(OOFEG_CRACK_PATTERN_WIDTH);
 *     if ((crackStatus == pscm_SOFTENING)||(crackStatus == pscm_OPEN))
 *      EASValsSetColor(gc.getActiveCrackColor());
 *     else
 *      EASValsSetColor(gc.getCrackPatternColor());
 *     EASValsSetFillStyle (FILL_HOLLOW);
 *     tr = CreateQuad3D (q);
 *     EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK | FILL_MASK, tr);
 *     EMAddGraphicsToModel(ESIModel(), tr);
 *    }
 *   }
 *  }
 *  //delete crackDir;
 * }
 * }
 * }
 * // check for valid stress-strain mode
 * if (!((gc.getDrawMode() == sxForce) || (gc.getDrawMode() == syForce) || (gc.getDrawMode() == szForce) || (gc.getDrawMode() == sxyForce))) return ;
 *
 * EASValsSetLayer(OOFEG_STRESS_CONTOUR_LAYER);
 * if (gc.getInternalVarsDefGeoFlag()) {
 * // use deformed geometry
 * p[0].x = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(1,tStep,DisplacementVector,defScale);
 * p[0].y = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(2,tStep,DisplacementVector,defScale);
 * p[0].z = 0.;
 * p[1].x = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(1,tStep,DisplacementVector,defScale);
 * p[1].y = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(2,tStep,DisplacementVector,defScale);
 * p[1].z = 0.;                    // delete help;
 * p[2].x = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(1,tStep,DisplacementVector,defScale);
 * p[2].y = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(2,tStep,DisplacementVector,defScale);
 * p[2].z = 0.;
 * } else {
 * p[0].x = (FPNum) this->giveNode(1)->giveCoordinate(1);
 * p[0].y = (FPNum) this->giveNode(1)->giveCoordinate(2);
 * p[0].z = 0.;
 * p[1].x = (FPNum) this->giveNode(2)->giveCoordinate(1);
 * p[1].y = (FPNum) this->giveNode(2)->giveCoordinate(2);
 * p[1].z = 0.;
 * p[2].x = (FPNum) this->giveNode(3)->giveCoordinate(1);
 * p[2].y = (FPNum) this->giveNode(3)->giveCoordinate(2);
 * p[2].z = 0.;
 * }
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
