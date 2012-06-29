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

#include "truss1d.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "mathfem.h"
#include "fei1dlin.h"

#ifdef __OOFEG
 #include "engngm.h"
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
FEI1dLin Truss1d :: interp(1); // Initiates the static interpolator


Truss1d :: Truss1d(int n, Domain *aDomain) :
    StructuralElement(n, aDomain),
    ZZNodalRecoveryModelInterface(), NodalAveragingRecoveryModelInterface(),
    SpatialLocalizerInterface(),
    DirectErrorIndicatorRCInterface(),
    EIPrimaryUnknownMapperInterface(), ZZErrorEstimatorInterface(), ZZRemeshingCriteriaInterface(),
    MMAShapeFunctProjectionInterface(), HuertaErrorEstimatorInterface(), HuertaRemeshingCriteriaInterface()
{
    numberOfDofMans = 2;
}


void Truss1d :: computeGaussPoints()
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
Truss1d :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver. This expression is
// valid in both local and global axes.
{
    Material *mat;
    double halfMass;
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    mat = this->giveMaterial();
    halfMass = mat->give('d', gp) * this->giveCrossSection()->give(CS_Area) * this->computeLength() * 0.5;
    answer.resize(2, 2);
    answer.zero();
    answer.at(1, 1) = halfMass;
    answer.at(2, 2) = halfMass;
}


void
Truss1d :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
    FloatMatrix dN;
    this->interp.evaldNdx( dN, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    answer.beTranspositionOf(dN); ///@todo It would be more suitable to follow the column-major version as done in FEI-classes
}


void
Truss1d :: computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver,
// evaluated at gp.
{
    FloatArray n;
    this->interp.evalN( n, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    // Reshape
    answer.resize(1, 2); ///@todo It would be more suitable to follow the column-major order and just do answer.setColumn(...)
    answer.at(1, 1) = n.at(1);
    answer.at(1, 2) = n.at(2);
}


double
Truss1d :: computeVolumeAround(GaussPoint *gp)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    double detJ = fabs( this->interp.giveTransformationJacobian( * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ) );
    return detJ *gp->giveWeight() * this->giveCrossSection()->give(CS_Area);
}


IRResultType
Truss1d :: initializeFrom(InputRecord *ir)
{
    this->StructuralElement :: initializeFrom(ir);
    this->computeGaussPoints();
    return IRRT_OK;
}


void
Truss1d :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(1, D_u);
}


void
Truss1d :: HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId,
                                                            IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                            HuertaErrorEstimatorInterface :: SetupMode sMode, TimeStep *tStep,
                                                            int &localNodeId, int &localElemId, int &localBcId,
                                                            IntArray &controlNode, IntArray &controlDof,
                                                            HuertaErrorEstimator :: AnalysisMode aMode)
{
    Element *element = this->HuertaErrorEstimatorI_giveElement();
    int inode, nodes = 2;
    FloatArray *corner [ 2 ], midNode, cor [ 2 ];
    double x = 0.0;

    if ( sMode == HuertaErrorEstimatorInterface :: NodeMode ||
        ( sMode == HuertaErrorEstimatorInterface :: BCMode && aMode == HuertaErrorEstimator :: HEE_linear ) ) {
        for ( inode = 0; inode < nodes; inode++ ) {
            corner [ inode ] = element->giveNode(inode + 1)->giveCoordinates();
            if ( corner [ inode ]->giveSize() != 3 ) {
                cor [ inode ].resize(3);
                cor [ inode ].at(1) = corner [ inode ]->at(1);
                cor [ inode ].at(2) = 0.0;
                cor [ inode ].at(3) = 0.0;

                corner [ inode ] = & ( cor [ inode ] );
            }

            x += corner [ inode ]->at(1);
        }

        midNode.resize(3);

        midNode.at(1) = x / nodes;
        midNode.at(2) = 0.0;
        midNode.at(3) = 0.0;
    }

    this->setupRefinedElementProblem1D(element, refinedElement, level, nodeId, localNodeIdArray, globalNodeIdArray,
                                       sMode, tStep, nodes, corner, midNode,
                                       localNodeId, localElemId, localBcId,
                                       controlNode, controlDof, aMode, "Truss1d");
}


#ifdef __OOFEG
void Truss1d :: drawRawGeometry(oofegGraphicContext &gc)
{
    GraphicObj *go;
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* point */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = 0.;
    p [ 0 ].z = 0.0;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = 0.;
    p [ 1 ].z = 0.0;
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void Truss1d :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
{
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
    p [ 0 ].y = 0.;
    p [ 0 ].z = 0.;

    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].y = 0.;
    p [ 1 ].z = 0.;
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


void Truss1d :: drawScalar(oofegGraphicContext &context)
{
    int i, indx, result = 0;
    WCRec p [ 2 ];
    GraphicObj *tr;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    FloatArray v1, v2;
    double s [ 2 ], defScale;
    IntArray map;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    if ( context.giveIntVarMode() == ISM_recovered ) {
        result += this->giveInternalStateAtNode(v1, context.giveIntVarType(), context.giveIntVarMode(), 1, tStep);
        result += this->giveInternalStateAtNode(v2, context.giveIntVarType(), context.giveIntVarMode(), 2, tStep);
    } else if ( context.giveIntVarMode() == ISM_local ) {
        GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
        result += giveIPValue(v1, gp, context.giveIntVarType(), tStep);
        v2 = v1;
        result *= 2;
    }

    if ( result != 2 ) {
        return;
    }

    result = this->giveIntVarCompFullIndx( map, context.giveIntVarType() );

    if ( ( !result ) || ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
        return;
    }

    s [ 0 ] = v1.at(indx);
    s [ 1 ] = v2.at(indx);

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);

    if ( ( context.getScalarAlgo() == SA_ISO_SURF ) || ( context.getScalarAlgo() == SA_ISO_LINE ) ) {
        for ( i = 0; i < 2; i++ ) {
            if ( context.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = context.getDefScale();
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
                p [ i ].y = 0.;
                p [ i ].z = 0.;
            } else {
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                p [ i ].y = 0.;
                p [ i ].z = 0.;
            }
        }

        //EASValsSetColor(gc.getYieldPlotColor(ratio));
        tr =  CreateLine3D(p);
        EGWithMaskChangeAttributes(LAYER_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    } else if ( ( context.getScalarAlgo() == SA_ZPROFILE ) || ( context.getScalarAlgo() == SA_COLORZPROFILE ) ) {
        double landScale = context.getLandScale();

        for ( i = 0; i < 2; i++ ) {
            if ( context.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = context.getDefScale();
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
                p [ i ].y = 0.0;
                p [ i ].z = s [ i ] * landScale;
            } else {
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                p [ i ].y = 0.0;
                p [ i ].z = s [ i ] * landScale;
            }
        }

        if ( context.getScalarAlgo() == SA_ZPROFILE ) {
            /*
             * EASValsSetColor(context.getDeformedElementColor());
             * EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
             * tr =  CreateLine3D(p);
             * EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
             */
            WCRec pp [ 4 ];
            pp [ 0 ].x = p [ 0 ].x;
            pp [ 0 ].y = 0.0;
            pp [ 0 ].z = 0.0;
            pp [ 1 ].x = p [ 0 ].x;
            pp [ 1 ].y = 0.0;
            pp [ 1 ].z = p [ 0 ].z;
            pp [ 2 ].x = p [ 1 ].x;
            pp [ 2 ].y = 0.0;
            pp [ 2 ].z = p [ 1 ].z;
            pp [ 3 ].x = p [ 1 ].x;
            pp [ 3 ].y = 0.0;
            pp [ 3 ].z = 0.0;
            tr = CreateQuad3D(pp);
            EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
            EASValsSetColor( context.getDeformedElementColor() );
            //EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
            EASValsSetFillStyle(FILL_HOLLOW);
            EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);
        } else {
            //tr =  CreateTriangleWD3D(p, s[0], s[1], s[2]);
            EASValsSetColor( context.getDeformedElementColor() );
            tr =  CreateLine3D(p);
            EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);
        }
    }
}



#endif


Interface *
Truss1d :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return ( ZZNodalRecoveryModelInterface * ) this;
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return ( NodalAveragingRecoveryModelInterface * ) this;
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


int
Truss1d :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) ) {
        return 1;
    }

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}


void
Truss1d :: ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatArray &answer, GaussPoint *gp, InternalStateType type)
{
    // evaluates N matrix (interpolation estimated stress matrix)
    // according to Zienkiewicz & Zhu paper
    // N(nsigma, nsigma*nnodes)
    // Definition : sigmaVector = N * nodalSigmaVector
    if ( !this->giveIPValueSize(type, gp) ) {
        return;
    }

    this->interp.evalN( answer, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
}


void
Truss1d :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                      InternalStateType type, TimeStep *tStep)
{
    GaussPoint *gp;
    gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    this->giveIPValue(answer, gp, type, tStep);
}


void
Truss1d :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                     InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}


int
Truss1d :: SpatialLocalizerI_containsPoint(const FloatArray &coords)
{
    FloatArray lc;
    return this->computeLocalCoordinates(lc, coords);
}


double
Truss1d :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray lcoords(1), gcoords;
    double dist;
    int size, gsize;

    lcoords.at(1) = 0.0;
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
Truss1d :: DirectErrorIndicatorRCI_giveCharacteristicSize()
{
    return this->computeLength();
}


int
Truss1d :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode,
                                                            TimeStep *stepN, const FloatArray &coords,
                                                            FloatArray &answer)
{
    FloatArray n, u, ksi;
    int result;

    result = this->computeLocalCoordinates(ksi, coords);
    this->interp.evalN( n, ksi, FEIElementGeometryWrapper(this) );
    this->computeVectorOf(EID_MomentumBalance, mode, stepN, u);
    answer.setValues( 1, n.dotProduct(u) );
    return result;
}


void
Truss1d :: EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer)
{
    giveDofManDofIDMask(1, EID_MomentumBalance, answer);
}


void
Truss1d :: MMAShapeFunctProjectionInterface_interpolateIntVarAt(FloatArray &answer, FloatArray &coords,
                                                                coordType ct, nodalValContainerType &list,
                                                                InternalStateType type, TimeStep *tStep)
{
    double vars;
    FloatArray lcoords, n;
    if ( ct == MMAShapeFunctProjectionInterface :: coordType_local ) {
        lcoords = coords;
    } else {
        computeLocalCoordinates(lcoords, coords);
    }

    this->interp.evalN( n, lcoords, FEIElementGeometryWrapper(this) );

    vars = list.at(1)->giveSize();
    answer.resize(vars);
    for ( int i = 1; i <= vars; i++ ) {
        answer.at(i) = n.at(1) * list.at(1)->at(i) + n.at(2) * list.at(2)->at(i);
    }
}
} // end namespace oofem
