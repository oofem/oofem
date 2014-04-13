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

#include "axisymm3d.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "engngm.h"
#include "mathfem.h"
#include "crosssection.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif


namespace oofem {
REGISTER_Element(Axisymm3d);

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
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >(this);
    }

    return NULL;
}


void
Axisymm3d :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [6x6] strain-displacement matrix {B} of the receiver, eva-
// luated at gp.
{
    double x, r;
    int size, ind = 1;
    FloatMatrix dnx;

    this->interpolation.evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );


    if ( ui == ALL_STRAINS ) {
        size = 6;
        ui = 6;
    } else {
        size = ui - li + 1;
    }

    if ( ( size < 0 ) || ( size > 6 ) ) {
        OOFEM_ERROR("ComputeBmatrixAt size mismatch");
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
        this->interpolation.evalN( n, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

        r = 0.;
        for ( int i = 1; i <= numberOfDofMans; i++ ) {
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


void
Axisymm3d :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the [9x6] displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
// BH matrix  -  9 rows : du/dx, dv/dy, dw/dz = u/r, 0, 0, du/dy,  0, 0, dv/dx
// @todo not checked if correct
{
    FloatArray n;
    FloatMatrix dnx;

    this->interpolation.evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(9, 6);
    answer.zero();

    double r = 0., x;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        x  = this->giveNode(i)->giveCoordinate(1);
        r += x * n.at(i);
    }


    // mode is _3dMat !!!!!! answer.at(4,*), answer.at(5,*), answer.at(7,*), and answer.at(8,*) is zero
    for ( int i = 1; i <= 6; i++ ) {
        answer.at(1, 3 * i - 2) = dnx.at(i, 1);     // du/dx
        answer.at(2, 3 * i - 1) = dnx.at(i, 2);     // dv/dy
        answer.at(6, 3 * i - 2) = dnx.at(i, 2);     // du/dy
        answer.at(9, 3 * i - 1) = dnx.at(i, 1);     // dv/dx
    }

    answer.at(3, 1) = n.at(1) / r;
    answer.at(3, 3) = n.at(2) / r;
    answer.at(3, 5) = n.at(3) / r;
}



double
Axisymm3d :: giveArea()
// returns the area occupied by the receiver
{
    if ( area > 0 ) {
        return area;         // check if previously computed
    }

    area = fabs( this->interpolation.giveArea( FEIElementGeometryWrapper(this) ) );
    return area;
}


double
Axisymm3d :: computeVolumeAround(GaussPoint *gp)
{
    double determinant, weight, volume, r, x;
    FloatArray n;

    this->interpolation.evalN( n, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

    r = 0.;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        x  = this->giveNode(i)->giveCoordinate(1);
        r += x * n.at(i);
    }

    determinant = fabs( this->interpolation.giveTransformationJacobian( * gp->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this) ) );

    weight = gp->giveWeight();
    volume = determinant * weight * r;

    return volume;
}


void
Axisymm3d :: computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(2);
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 2);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
        integrationRulesArray [ 1 ] = new GaussIntegrationRule(2, this, 3, 6);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 1 ], numberOfFiAndShGaussPoints, this);
    }
}


IRResultType
Axisymm3d :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    numberOfGaussPoints = 1;
    result = this->StructuralElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    numberOfFiAndShGaussPoints = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfFiAndShGaussPoints, _IFT_Axisymm3d_nipfish);

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

    return IRRT_OK;
}


void
Axisymm3d :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step tStep. The nature of these strains depends
// on the element's type.
{
    FloatMatrix b, A;
    FloatArray u, Epsilon, help;
    fMode mode = domain->giveEngngModel()->giveFormulation();

    answer.resize(6);
    answer.zero();

    if ( mode == TL ) {  // Total Lagrange formulation
        this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
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
            OOFEM_ERROR("numberOfRotGaussPoints size mismatch");
        }

        Epsilon.beProductOf(b, u);
        answer.at(3) = Epsilon.at(1);
        answer.at(6) = Epsilon.at(4);

        if ( nlGeometry ) {
            OOFEM_ERROR("only supports nlGeometry = 0");
        }
    } else if ( mode == AL ) { // actualized Lagrange formulation
        OOFEM_ERROR("unsupported mode");
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
        for ( int i = 1; i <= numberOfDofMans; i++ ) {
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
    answer = {D_u, D_v};
}


void
Axisymm3d :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                        InternalStateType type, TimeStep *tStep)
{
    GaussPoint *gp;

    if ( numberOfGaussPoints != 1 ) {
        answer.clear(); // for more gp's need to be refined
        return;
    }

    gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    this->giveIPValue(answer, gp, type, tStep);
}

void
Axisymm3d :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                       InternalStateType type, TimeStep *tStep)
{
    answer.clear();
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
        OOFEM_ERROR("node unknown");
    }
}


int
Axisymm3d :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return this->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints();
}


SPRPatchType
Axisymm3d :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dxy;
}


void
Axisymm3d :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
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
    this->interpolation.edgeEvalN( n, iedge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

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
        OOFEM_ERROR("wrong edge number");
    }
}


double
Axisymm3d :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    FloatArray c(2);
    this->computeEdgeIpGlobalCoords(c, gp, iEdge);
    double result = this->interpolation.edgeGiveTransformationJacobian( iEdge, * gp->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this) );


    return c.at(1) * result * gp->giveWeight();
}


void
Axisymm3d :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    this->interpolation.edgeLocal2global( answer, iEdge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
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
        OOFEM_ERROR("wrong egde number");
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
    return this->interpolation.global2local( lcoords, coords, FEIElementGeometryWrapper(this) );
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
        OOFEM_ERROR("coordinates size mismatch");
    }

    if ( size == gsize ) {
        dist = coords.distance(gcoords);
    } else {
        FloatArray helpCoords = coords;

        helpCoords.resizeWithValues(gsize);
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
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, defScale);
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

    indx = context.giveIntVarIndx();

    s [ 0 ] = v1.at(indx);
    s [ 1 ] = v2.at(indx);
    s [ 2 ] = v3.at(indx);

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);

    if ( context.getScalarAlgo() == SA_ISO_SURF ) {
        for ( i = 0; i < 3; i++ ) {
            if ( context.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = context.getDefScale();
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
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
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
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

#endif
} // end namespace oofem
