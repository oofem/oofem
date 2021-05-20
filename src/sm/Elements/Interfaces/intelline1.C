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

#include "sm/Elements/Interfaces/intelline1.h"
#include "sm/CrossSections/structuralinterfacecrosssection.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "lobattoir.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "fei2dlinelin.h"
#include "classfactory.h"
#include "nodalaveragingrecoverymodel.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include <Emarkwd3d.h>
#endif

namespace oofem {
REGISTER_Element(IntElLine1);

FEI2dLineLin IntElLine1 :: interp(1, 1);


IntElLine1 :: IntElLine1(int n, Domain *aDomain) :
    StructuralInterfaceElement(n, aDomain)
{
    numberOfDofMans = 4;
    numberOfGaussPoints = 4;
}


void
IntElLine1 :: computeNmatrixAt(GaussPoint *ip, FloatMatrix &answer)
{
    // Returns the modified N-matrix which multiplied with u give the spatial jump.
    auto N = interp.evalN(ip->giveNaturalCoordinates().at(1));

    answer.resize(2, 8);
    answer.zero();
    answer.at(1, 1) = answer.at(2, 2) = -N.at(1);
    answer.at(1, 3) = answer.at(2, 4) = -N.at(2);

    answer.at(1, 5) = answer.at(2, 6) = N.at(1);
    answer.at(1, 7) = answer.at(2, 8) = N.at(2);
}


void
IntElLine1 :: computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );

//        integrationRulesArray[ 0 ] = std::make_unique<LobattoIntegrationRule>(1,this, 1, 2, false);
//        integrationRulesArray [ 0 ]->SetUpPointsOnLine(2, _2dInterface);

        integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this, 1, 2);
        integrationRulesArray [ 0 ]->SetUpPointsOnLine(this->numberOfGaussPoints, _2dInterface);
    }
}

FloatArrayF<2>
IntElLine1 :: computeCovarBaseVectorAt(IntegrationPoint *ip) const
{
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdxi;
    interp->evaldNdxi( dNdxi, ip->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    FloatArrayF<2> G;
    int numNodes = this->giveNumberOfNodes();
    for ( int i = 1; i <= dNdxi.giveNumberOfRows(); i++ ) {
        double X1_i = 0.5 * ( this->giveNode(i)->giveCoordinate(1) + this->giveNode(i + numNodes / 2)->giveCoordinate(1) ); // (mean) point on the fictious mid surface
        double X2_i = 0.5 * ( this->giveNode(i)->giveCoordinate(2) + this->giveNode(i + numNodes / 2)->giveCoordinate(2) );
        G.at(1) += dNdxi.at(i, 1) * X1_i;
        G.at(2) += dNdxi.at(i, 1) * X2_i;
    }
    return G;
}

double
IntElLine1 :: computeAreaAround(IntegrationPoint *ip)
{
    auto G = this->computeCovarBaseVectorAt(ip);

    double weight = ip->giveWeight();
    double ds = norm(G) * weight;
    if ( this->axisymmode ) {
        int numNodes = this->giveNumberOfNodes();
        auto N = this->interp.evalN(ip->giveNaturalCoordinates().at(1));
        // interpolate radius
        double r = 0.0;
        for ( int i = 1; i <= N.giveSize(); i++ ) {
            double X_i = 0.5 * ( this->giveNode(i)->giveCoordinate(1) + this->giveNode(i + numNodes / 2)->giveCoordinate(1) ); // X-coord of the fictious mid surface
            r += N.at(i) * X_i;
        }
        return ds * r;
    } else { // regular 2d
        double thickness  = this->giveCrossSection()->give(CS_Thickness, ip);
        return ds * thickness;
    }
}


void
IntElLine1 :: initializeFrom(InputRecord &ir)
{
    StructuralInterfaceElement :: initializeFrom(ir);

    this->axisymmode = ir.hasField(_IFT_IntElLine1_axisymmode);

    // Check if node numbering is ok
    int nodeInd1 = this->giveDofManagerNumber(1);
    int arrayInd1 = domain->giveDofManPlaceInArray(nodeInd1);
    DofManager *node1 = domain->giveDofManager(arrayInd1);
    const auto &x1 = node1->giveCoordinates();

//    DofManager *node2 = this->giveDofManager(2);
    int nodeInd2 = this->giveDofManagerNumber(2);
    int arrayInd2 = domain->giveDofManPlaceInArray(nodeInd2);
    DofManager *node2 = domain->giveDofManager(arrayInd2);
    const auto &x2 = node2->giveCoordinates();

//    DofManager *node3 = this->giveDofManager(3);
    int nodeInd3 = this->giveDofManagerNumber(3);
    int arrayInd3 = domain->giveDofManPlaceInArray(nodeInd3);
    DofManager *node3 = domain->giveDofManager(arrayInd3);
    const auto &x3 = node3->giveCoordinates();


    double L2 = distance_square(x1, x2);
    double L3 = distance_square(x1, x3);

    if ( L2 < L3 ) {
        printf("Renumbering element %d\n.\n", this->giveNumber());
        dofManArray = {dofManArray.at(3), dofManArray.at(1), dofManArray.at(4), dofManArray.at(2)};
    }
}


void
IntElLine1 :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v};
}

void
IntElLine1 :: computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    // Transformation matrix to the local coordinate system
    // xy plane
    auto G = this->computeCovarBaseVectorAt(gp);
    G /= norm(G);

    answer.resize(2, 2);
//     answer.at(1, 1) =  G.at(1);//tangent vector
//     answer.at(2, 1) = -G.at(2);
//     answer.at(1, 2) =  G.at(2);
//     answer.at(2, 2) =  G.at(1);
    //normal is -G.at(2), G.at(1), perpendicular to nodes 1 2
    answer.at(1, 1) = -G.at(2);//normal vector
    answer.at(2, 1) = G.at(1);
    answer.at(1, 2) = G.at(1);
    answer.at(2, 2) = G.at(2);
}

FEInterpolation *
IntElLine1 :: giveInterpolation() const
{
    return & interp;
}


Interface *IntElLine1::giveInterface( InterfaceType interface )
{
    if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast<NodalAveragingRecoveryModelInterface *>(this);
    } else {
        return nullptr;
    }
}

void IntElLine1::NodalAveragingRecoveryMI_computeNodalValue( FloatArray &answer, int node, InternalStateType type, TimeStep *tStep )
{
    FloatArray nodeMap, gpValsA, gpValsB, N, locCoord;
    nodeMap.resize(1);

    FEInterpolation *interp = this->giveInterpolation();

    //Since it's linear interpolation, it suffices to take first and last Gauss point
    //and recover nodal values using only them.

    GaussPoint* gpA = integrationRulesArray[0]->getIntegrationPoint(0);
    GaussPoint* gpB = integrationRulesArray[0]->getIntegrationPoint(numberOfGaussPoints-1);

    //Get internal values at Gauss Points A and B
    this->giveIPValue(gpValsA, gpA, type, tStep);
    this->giveIPValue(gpValsB, gpB, type, tStep);

    if ( (node == 1 ) || (node == 3) ) { //nodes 1 and 3 pertain to gpA
        locCoord = gpA->giveNaturalCoordinates();
    } else if ( (node == 2 ) || (node == 4) ) { //nodes 2 and 4 pertain to gpB
        locCoord = gpB->giveNaturalCoordinates();
    }

    nodeMap.at(1) = 1/locCoord.at(1);
    //Evaluate shape functions (associated with Gauss Points) at element nodes
    interp->evalN( N, nodeMap, FEIElementGeometryWrapper(this) );
    answer.resize( gpValsA.giveSize() );
    answer.at(1)= N.at(1)*gpValsA.at(1) + N.at(2)*gpValsB.at(1);
    answer.at(2)= N.at(1)*gpValsA.at(2) + N.at(2)*gpValsB.at(2);
    answer.at(3)= N.at(1)*gpValsA.at(3) + N.at(2)*gpValsB.at(3);
}


#ifdef __OOFEG
void IntElLine1 :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
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
    p [ 0 ].z = 0.0;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = 0.0;
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void IntElLine1 :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    GraphicObj *go;
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* poin */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    double defScale = gc.getDefScale();

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER + 1);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = 0.0;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = 0.0;
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);

    p [ 0 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = 0.0;
    p [ 1 ].x = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = 0.0;
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


void IntElLine1 :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int indx, result = 0;
    IntegrationRule *iRule = this->giveDefaultIntegrationRulePtr();
    FloatArray gcoord(3), v1;
    WCRec p [ 1 ];
    GraphicObj *go;
    double val [ 1 ];

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( gc.giveIntVarMode() == ISM_recovered ) {
        return;
    }

    for ( GaussPoint *gp: *iRule ) {
        result = 0;
        result += giveIPValue(v1, gp, gc.giveIntVarType(), tStep);
        if ( result != 1 ) {
            continue;
        }

        indx = gc.giveIntVarIndx();

        result += this->computeGlobalCoordinates( gcoord, gp->giveNaturalCoordinates() );

        p [ 0 ].x = ( FPNum ) gcoord.at(1);
        p [ 0 ].y = ( FPNum ) gcoord.at(2);
        p [ 0 ].z = 0.;

        val [ 0 ] = v1.at(indx);
        gc.updateFringeTableMinMax(val, 1);
        //if (val[0] > 0.) {

        EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
        EASValsSetMType(FILLED_CIRCLE_MARKER);
        go = CreateMarkerWD3D(p, val [ 0 ]);
        EGWithMaskChangeAttributes(LAYER_MASK | FILL_MASK | MTYPE_MASK, go);
        EMAddGraphicsToModel(ESIModel(), go);
        //}
    }
}

#endif

} // end namespace oofem
