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

#include "ltrspace.h"
#include "node.h"
#include "material.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "mathfem.h"

#ifdef __OOFEG
 #include "engngm.h"
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
 #include "rcm2.h"
 #ifndef __MAKEDEPEND
  #include <Etetrawd.h>
 #endif
#endif

namespace oofem {
FEI3dTrLin LTRSpace :: interpolation;

LTRSpace :: LTRSpace(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(), NodalAveragingRecoveryModelInterface(),
    SPRNodalRecoveryModelInterface(), SpatialLocalizerInterface(), DirectErrorIndicatorRCInterface(),
    EIPrimaryUnknownMapperInterface(), ZZErrorEstimatorInterface(), ZZRemeshingCriteriaInterface(),
    MMAShapeFunctProjectionInterface(), HuertaErrorEstimatorInterface(), HuertaRemeshingCriteriaInterface()

{
    numberOfDofMans  = 4;
    numberOfGaussPoints = 1;
}


Interface *
LTRSpace :: giveInterface(InterfaceType interface)
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
LTRSpace :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the [6x12] strain-displacement matrix {B} of the receiver, eva-
// luated at aGaussPoint.
{
    int i;
    FloatMatrix dn(4, 3);
    interpolation.evaldNdx(dn, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(6, 12);
    answer.zero();

    for ( i = 1; i <= 4; i++ ) {
        answer.at(1, 3 * i - 2) = dn.at(i, 1);
        answer.at(2, 3 * i - 1) = dn.at(i, 2);
        answer.at(3, 3 * i - 0) = dn.at(i, 3);

        answer.at(4, 3 * i - 1) = dn.at(i, 3);
        answer.at(4, 3 * i - 0) = dn.at(i, 2);

        answer.at(5, 3 * i - 2) = dn.at(i, 3);
        answer.at(5, 3 * i - 0) = dn.at(i, 1);

        answer.at(6, 3 * i - 2) = dn.at(i, 2);
        answer.at(6, 3 * i - 1) = dn.at(i, 1);
    }
}




void
LTRSpace :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
{
    int i;
    FloatArray n(4);

    this->interpolation.evalN(n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(3, 12);
    answer.zero();

    for ( i = 1; i <= 4; i++ ) {
        answer.at(1, 3 * i - 2)  = n.at(i);
        answer.at(2, 3 * i - 1)  = n.at(i);
        answer.at(3, 3 * i - 0)  = n.at(i);
    }
}

void
LTRSpace :: computeNLBMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint, int i)
//
// Returns the [12x12] nonlinear part of the strain-displacement matrix {B} of the receiver,
// evaluated at aGaussPoint
{
    int j, k, l;
    FloatMatrix dnx;

    interpolation.evaldNdx(dnx, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(12, 12);
    answer.zero();

    if ( i <= 3 ) {
        for ( k = 0; k < 4; k++ ) {
            for ( l = 0; l < 3; l++ ) {
                for ( j = 1; j <= 12; j += 3 ) {
                    answer.at(k * 3 + l + 1, l + j) = dnx.at(k + 1, i) * dnx.at( ( j - 1 ) / 3 + 1, i );
                }
            }
        }
    } else if ( i == 4 )       {
        for ( k = 0; k < 4; k++ ) {
            for ( l = 0; l < 3; l++ ) {
                for ( j = 1; j <= 12; j += 3 ) {
                    answer.at(k * 3 + l + 1, l + j) = dnx.at(k + 1, 2) * dnx.at( ( j - 1 ) / 3 + 1, 3 ) + dnx.at(k + 1, 3) * dnx.at( ( j - 1 ) / 3 + 1, 2 );
                }
            }
        }
    } else if ( i == 5 )        {
        for ( k = 0; k < 4; k++ ) {
            for ( l = 0; l < 3; l++ ) {
                for ( j = 1; j <= 12; j += 3 ) {
                    answer.at(k * 3 + l + 1, l + j) = dnx.at(k + 1, 1) * dnx.at( ( j - 1 ) / 3 + 1, 3 ) + dnx.at(k + 1, 3) * dnx.at( ( j - 1 ) / 3 + 1, 1 );
                }
            }
        }
    } else if ( i == 6 )       {
        for ( k = 0; k < 4; k++ ) {
            for ( l = 0; l < 3; l++ ) {
                for ( j = 1; j <= 12; j += 3 ) {
                    answer.at(k * 3 + l + 1, l + j) = dnx.at(k + 1, 1) * dnx.at( ( j - 1 ) / 3 + 1, 2 ) + dnx.at(k + 1, 2) * dnx.at( ( j - 1 ) / 3 + 1, 1 );
                }
            }
        }
    }
}

double LTRSpace :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
    double determinant, weight, volume;
    determinant = fabs( this->interpolation.giveTransformationJacobian(* aGaussPoint->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this)) );
    weight      = aGaussPoint->giveWeight();
    volume      = determinant * weight;
    return volume;
}

IRResultType
LTRSpace :: initializeFrom(InputRecord *ir)
{
    //const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    //IRResultType result;                            // Required by IR_GIVE_FIELD macro

    this->NLStructuralElement :: initializeFrom(ir);
    numberOfGaussPoints = 1;

    // set up Gaussian integration points
    this->computeGaussPoints();

    return IRRT_OK;
}


void LTRSpace :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 6);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Tetrahedra, numberOfGaussPoints, _3dMat);
    }
}


void
LTRSpace :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver.
{
    GaussPoint *gp;
    double dV, mss1;

    answer.resize(12, 12);
    answer.zero();
    gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    dV = this->computeVolumeAround(gp);
    mss1 = dV * this->giveMaterial()->give('d', gp) / 4.;

    for ( int i = 1; i <= 12; i++ ) {
        answer.at(i, i) = mss1;
    }
}


void
LTRSpace :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(3, D_u, D_v, D_w);
}


double
LTRSpace :: giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane)
{
    return this->giveLenghtInDir(normalToCrackPlane);
}


int
LTRSpace :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) ) {
        return 6;
    }

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}


void
LTRSpace :: ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type)
{
    // evaluates N matrix (interpolation estimated stress matrix)
    // according to Zienkiewicz & Zhu paper
    // N(nsigma, nsigma*nnodes)
    // Definition : sigmaVector = N * nodalSigmaVector

    if ( this->giveIPValueSize(type, aGaussPoint) ) {
        this->interpolation.evalN(answer, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this));
    } else {
        return;
    }
}

void
LTRSpace :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                       InternalStateType type, TimeStep *tStep)
{
    GaussPoint *gp;

    if ( numberOfGaussPoints != 1 ) {
        answer.resize(0); // for more gp's need to be refined
        return;
    }

    gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    giveIPValue(answer, gp, type, tStep);
}

void
LTRSpace :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                      InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}


void
LTRSpace :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(4);
    pap.at(1) = this->giveNode(1)->giveNumber();
    pap.at(2) = this->giveNode(2)->giveNumber();
    pap.at(3) = this->giveNode(3)->giveNumber();
    pap.at(4) = this->giveNode(4)->giveNumber();
}


void
LTRSpace :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int i, found = 0;
    answer.resize(1);

    for ( i = 1; i <= 4; i++ ) {
        if ( this->giveNode(i)->giveNumber() == pap ) {
            found = 1;
        }
    }

    if ( found ) {
        answer.at(1) = pap;
    } else {
        _error("SPRNodalRecoveryMI_giveDofMansDeterminedByPatch: node unknown");
    }
}


int
LTRSpace :: SPRNodalRecoveryMI_giveNumberOfIP()
{ return this->giveDefaultIntegrationRulePtr()->getNumberOfIntegrationPoints(); }


void
LTRSpace :: SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp)
{
    this->computeGlobalCoordinates( coords, * gp->giveCoordinates() );
}


SPRPatchType
LTRSpace :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_3dBiLin;
}


void
LTRSpace :: HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId,
                                                             IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                             HuertaErrorEstimatorInterface :: SetupMode sMode, TimeStep *tStep,
                                                             int &localNodeId, int &localElemId, int &localBcId,
                                                             IntArray &controlNode, IntArray &controlDof,
                                                             HuertaErrorEstimator :: AnalysisMode aMode)
{
    Element *element = this->HuertaErrorEstimatorI_giveElement();
    FloatArray *corner [ 4 ], midSide [ 6 ], midFace [ 4 ], midNode;
    double x = 0.0, y = 0.0, z = 0.0;
    int inode, nodes = 4, iside, sides = 6, iface, faces = 4, nd, nd1, nd2;

    static int sideNode [ 6 ] [ 2 ] = { { 1, 2 }, { 2, 3 }, { 3, 1 }, { 1, 4 }, { 2, 4 }, { 3, 4 } };
    static int faceNode [ 4 ] [ 3 ] = { { 1, 2, 3 }, { 1, 2, 4 }, { 2, 3, 4 }, { 3, 1, 4 } };

    /* ordering of hexa nodes must be compatible with refinedElement connectivity ordering;
     * generally the ordering is: corner side side face side face face center;
     * however the concrete ordering is element dependent (see refineMeshGlobally source if in doubts) */

    int hexaSideNode [ 4 ] [ 3 ] = { { 1, 3, 4 }, { 2, 1, 5 }, { 3, 2, 6 }, { 4, 6, 5 } };
    int hexaFaceNode [ 4 ] [ 3 ] = { { 1, 2, 4 }, { 1, 3, 2 }, { 1, 4, 3 }, { 4, 2, 3 } };

    if ( sMode == HuertaErrorEstimatorInterface :: NodeMode ||
        ( sMode == HuertaErrorEstimatorInterface :: BCMode && aMode == HuertaErrorEstimator :: HEE_linear ) ) {
        for ( inode = 0; inode < nodes; inode++ ) {
            corner [ inode ] = element->giveNode(inode + 1)->giveCoordinates();

            x += corner [ inode ]->at(1);
            y += corner [ inode ]->at(2);
            z += corner [ inode ]->at(3);
        }

        for ( iside = 0; iside < sides; iside++ ) {
            midSide [ iside ].resize(3);

            nd1 = sideNode [ iside ] [ 0 ] - 1;
            nd2 = sideNode [ iside ] [ 1 ] - 1;

            midSide [ iside ].at(1) = ( corner [ nd1 ]->at(1) + corner [ nd2 ]->at(1) ) / 2.0;
            midSide [ iside ].at(2) = ( corner [ nd1 ]->at(2) + corner [ nd2 ]->at(2) ) / 2.0;
            midSide [ iside ].at(3) = ( corner [ nd1 ]->at(3) + corner [ nd2 ]->at(3) ) / 2.0;
        }

        midNode.resize(3);

        midNode.at(1) = x / nodes;
        midNode.at(2) = y / nodes;
        midNode.at(3) = z / nodes;

        for ( iface = 0; iface < faces; iface++ ) {
            x = y = z = 0.0;
            for ( inode = 0; inode < 3; inode++ ) {
                nd = faceNode [ iface ] [ inode ] - 1;
                x += corner [ nd ]->at(1);
                y += corner [ nd ]->at(2);
                z += corner [ nd ]->at(3);
            }

            midFace [ iface ].resize(3);

            midFace [ iface ].at(1) = x / 3;
            midFace [ iface ].at(2) = y / 3;
            midFace [ iface ].at(3) = z / 3;
        }
    }

    this->setupRefinedElementProblem3D(element, refinedElement, level, nodeId, localNodeIdArray, globalNodeIdArray,
                                       sMode, tStep, nodes, corner, midSide, midFace, midNode,
                                       localNodeId, localElemId, localBcId, hexaSideNode, hexaFaceNode,
                                       controlNode, controlDof, aMode, "LSpace");
}

#ifdef __OOFEG
 #define TR_LENGHT_REDUCT 0.3333

void LTRSpace :: drawRawGeometry(oofegGraphicContext &gc)
{
    WCRec p [ 4 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    EASValsSetFillStyle(FILL_SOLID);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 2 ].z = ( FPNum ) this->giveNode(3)->giveCoordinate(3);
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveCoordinate(1);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveCoordinate(2);
    p [ 3 ].z = ( FPNum ) this->giveNode(4)->giveCoordinate(3);

    go =  CreateTetra(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void LTRSpace :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
{
    WCRec p [ 4 ];
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
    EASValsSetFillStyle(FILL_SOLID);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 2 ].z = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 3 ].z = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);

    go =  CreateTetra(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


void LTRSpace :: drawScalar(oofegGraphicContext &context)
{
    int i, indx, result = 0;
    WCRec p [ 4 ];
    GraphicObj *tr;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    FloatArray v [ 4 ];
    double s [ 4 ], defScale = 0.0;
    IntArray map;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    if ( context.giveIntVarMode() == ISM_recovered ) {
        for ( i = 1; i <= 4; i++ ) {
            result += this->giveInternalStateAtNode(v [ i - 1 ], context.giveIntVarType(), context.giveIntVarMode(), i, tStep);
        }

        if ( result != 4 ) {
            return;
        }
    } else if ( context.giveIntVarMode() == ISM_local ) {
        return;
    }

    this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
    if ( ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
        return;
    }

    for ( i = 1; i <= 4; i++ ) {
        s [ i - 1 ] = v [ i - 1 ].at(indx);
    }

    EASValsSetEdgeColor( context.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    if ( context.getScalarAlgo() == SA_ISO_SURF ) {
        for ( i = 0; i < 4; i++ ) {
            if ( context.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
                p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);
            } else {
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(3);
            }
        }

        context.updateFringeTableMinMax(s, 4);
        tr = CreateTetraWD(p, s);
        EGWithMaskChangeAttributes(LAYER_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    }
}

void
LTRSpace :: drawSpecial(oofegGraphicContext &gc)
{
    int i, j, k;
    WCRec q [ 4 ];
    GraphicObj *tr;
    StructuralMaterial *mat = ( StructuralMaterial * ) this->giveMaterial();
    IntegrationRule *iRule = integrationRulesArray [ 0 ];
    GaussPoint *gp;
    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    double defScale = gc.getDefScale();
    FloatArray crackStatuses, cf;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( gc.giveIntVarType() == IST_CrackState ) {
        int crackStatus;
        double xc, yc, zc, length;
        FloatArray crackDir;

        if ( numberOfGaussPoints != 1 ) {
            return;
        }

        //   for (igp=1 ; igp<= numberOfGaussPoints ; igp++) {
        {
            gp = iRule->getIntegrationPoint(0);
            if ( mat->giveIPValue(cf, gp, IST_CrackedFlag, tStep) == 0 ) {
                return;
            }

            if ( ( int ) cf.at(1) == 0 ) {
                return;
            }

            //
            // obtain gp global coordinates - here only one exists
            // it is in centre of gravity.
            xc = yc = zc = 0.;
            for ( i = 0; i < 4; i++ ) {
                if ( gc.getInternalVarsDefGeoFlag() ) {
                    // use deformed geometry
                    xc += ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
                    yc += ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
                    zc += ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);
                } else {
                    xc += ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                    yc += ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                    zc += ( FPNum ) this->giveNode(i + 1)->giveCoordinate(3);
                }
            }

            xc = xc / 4.;
            yc = yc / 4.;
            zc = zc / 4.;
            length = TR_LENGHT_REDUCT * pow(this->computeVolumeAround(gp), 1. / 3.) / 2.0;
            if ( mat->giveIPValue(crackDir, gp, IST_CrackDirs, tStep) ) {
                mat->giveIPValue(crackStatuses, gp, IST_CrackStatuses, tStep);


                for ( i = 1; i <= 3; i++ ) {
                    crackStatus = ( int ) crackStatuses.at(i);
                    if ( ( crackStatus != pscm_NONE ) && ( crackStatus != pscm_CLOSED ) ) {
                        // draw a crack
                        // this element is 3d element

                        if ( i == 1 ) {
                            j = 2;
                            k = 3;
                        } else if ( i == 2 ) {
                            j = 3;
                            k = 1;
                        } else {
                            j = 1;
                            k = 2;
                        }

                        q [ 0 ].x = ( FPNum ) xc + 0.5 * crackDir.at(0 + j) * length + 0.5 * crackDir.at(0 + k) * length;
                        q [ 0 ].y = ( FPNum ) yc + 0.5 * crackDir.at(3 + j) * length + 0.5 * crackDir.at(3 + k) * length;
                        q [ 0 ].z = ( FPNum ) zc + 0.5 * crackDir.at(6 + j) * length + 0.5 * crackDir.at(6 + k) * length;
                        q [ 1 ].x = ( FPNum ) xc + 0.5 * crackDir.at(0 + j) * length - 0.5 * crackDir.at(0 + k) * length;
                        q [ 1 ].y = ( FPNum ) yc + 0.5 * crackDir.at(3 + j) * length - 0.5 * crackDir.at(3 + k) * length;
                        q [ 1 ].z = ( FPNum ) zc + 0.5 * crackDir.at(6 + j) * length - 0.5 * crackDir.at(6 + k) * length;
                        q [ 2 ].x = ( FPNum ) xc - 0.5 * crackDir.at(0 + j) * length - 0.5 * crackDir.at(0 + k) * length;
                        q [ 2 ].y = ( FPNum ) yc - 0.5 * crackDir.at(3 + j) * length - 0.5 * crackDir.at(3 + k) * length;
                        q [ 2 ].z = ( FPNum ) zc - 0.5 * crackDir.at(6 + j) * length - 0.5 * crackDir.at(6 + k) * length;
                        q [ 3 ].x = ( FPNum ) xc - 0.5 * crackDir.at(0 + j) * length + 0.5 * crackDir.at(0 + k) * length;
                        q [ 3 ].y = ( FPNum ) yc - 0.5 * crackDir.at(3 + j) * length + 0.5 * crackDir.at(3 + k) * length;
                        q [ 3 ].z = ( FPNum ) zc - 0.5 * crackDir.at(6 + j) * length + 0.5 * crackDir.at(6 + k) * length;

                        EASValsSetLayer(OOFEG_CRACK_PATTERN_LAYER);
                        EASValsSetLineWidth(OOFEG_CRACK_PATTERN_WIDTH);
                        if ( ( crackStatus == pscm_SOFTENING ) || ( crackStatus == pscm_OPEN ) ) {
                            EASValsSetColor( gc.getActiveCrackColor() );
                        } else {
                            EASValsSetColor( gc.getCrackPatternColor() );
                        }

                        //      EASValsSetFillStyle (FILL_HOLLOW);
                        tr = CreateQuad3D(q);
                        EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
                        EMAddGraphicsToModel(ESIModel(), tr);
                    }
                }
            }
        }
    }
}

#endif

int
LTRSpace :: SpatialLocalizerI_containsPoint(const FloatArray &coords)
{
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}


double
LTRSpace :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray lcoords(4), gcoords;
    double dist;
    int size, gsize;

    lcoords.at(1) = lcoords.at(2) = lcoords.at(3) = lcoords.at(4) = 1. / 4.;
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
LTRSpace :: DirectErrorIndicatorRCI_giveCharacteristicSize() {
    FloatArray lc(4);
    double volume = interpolation.giveTransformationJacobian(lc, FEIElementGeometryWrapper(this));
    return pow(volume * 6.0, 1. / 3.);
}


int
LTRSpace :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode,
                                                             TimeStep *stepN, const FloatArray &coords,
                                                             FloatArray &answer)
{
    FloatArray lcoords, u;
    FloatMatrix n;
    int result;

    result = this->computeLocalCoordinates(lcoords, coords);

    n.resize(3, 12);
    n.zero();

    n.at(1, 1)  = lcoords.at(1);
    n.at(1, 4)  = lcoords.at(2);
    n.at(1, 7)  = lcoords.at(3);
    n.at(1, 10) = lcoords.at(4);

    n.at(2, 2)  = lcoords.at(1);
    n.at(2, 5)  = lcoords.at(2);
    n.at(2, 8)  = lcoords.at(3);
    n.at(2, 11) = lcoords.at(4);

    n.at(3, 3)  = lcoords.at(1);
    n.at(3, 6)  = lcoords.at(2);
    n.at(3, 9)  = lcoords.at(3);
    n.at(3, 12) = lcoords.at(4);

    this->computeVectorOf(EID_MomentumBalance, mode, stepN, u);
    answer.beProductOf(n, u);

    return result;
}

void
LTRSpace :: EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer)
{
    giveDofManDofIDMask(1, EID_MomentumBalance, answer);
}

void
LTRSpace :: MMAShapeFunctProjectionInterface_interpolateIntVarAt(FloatArray &answer, FloatArray &coords,
                                                                 coordType ct, nodalValContainerType &list,
                                                                 InternalStateType type, TimeStep *tStep)
{
    int i, n;
    double l1, l2, l3, l4;
    FloatArray lcoords;
    if ( ct == MMAShapeFunctProjectionInterface :: coordType_local ) {
        lcoords = coords;
    } else {
        computeLocalCoordinates(lcoords, coords);
    }

    l1 = lcoords.at(1);
    l2 = lcoords.at(2);
    l3 = lcoords.at(3);
    l4 = 1.0 - l1 - l2 - l3;
    n = list.at(1)->giveSize();
    answer.resize(n);

    for ( i = 1; i <= n; i++ ) {
        answer.at(i) = l1 * list.at(1)->at(i) + l2 *list.at(2)->at(i) + l3 *list.at(3)->at(i) + l4 *list.at(4)->at(i);
    }
}


void
LTRSpace :: computeEgdeNMatrixAt(FloatMatrix &answer, GaussPoint *aGaussPoint)
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

    answer.resize(3, 6);
    answer.zero();

    answer.at(1, 1) = n.at(1);
    answer.at(1, 4) = n.at(2);
    answer.at(2, 2) = n.at(1);
    answer.at(2, 5) = n.at(2);
    answer.at(3, 3) = n.at(1);
    answer.at(3, 6) = n.at(2);
}

void
LTRSpace :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    answer.resize(6);
    if ( iEdge == 1 ) { // edge between nodes 1 2
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 3;
        answer.at(4) = 4;
        answer.at(5) = 5;
        answer.at(6) = 6;
    } else if ( iEdge == 2 ) { // edge between nodes 2 3
        answer.at(1) = 4;
        answer.at(2) = 5;
        answer.at(3) = 6;
        answer.at(4) = 7;
        answer.at(5) = 8;
        answer.at(6) = 9;
    } else if ( iEdge == 3 ) { // edge between nodes 3 1
        answer.at(1) = 7;
        answer.at(2) = 8;
        answer.at(3) = 9;
        answer.at(4) = 1;
        answer.at(5) = 2;
        answer.at(6) = 3;
    } else if ( iEdge == 4 ) { // edge between nodes 1 4
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 3;
        answer.at(4) = 10;
        answer.at(5) = 11;
        answer.at(6) = 12;
    } else if ( iEdge == 5 ) { // edge between nodes 2 4
        answer.at(1) = 4;
        answer.at(2) = 5;
        answer.at(3) = 6;
        answer.at(4) = 10;
        answer.at(5) = 11;
        answer.at(6) = 12;
    } else if ( iEdge == 6 ) { // edge between nodes 3 4
        answer.at(1) = 7;
        answer.at(2) = 8;
        answer.at(3) = 9;
        answer.at(4) = 10;
        answer.at(5) = 11;
        answer.at(6) = 12;
    } else {
        _error("giveEdgeDofMapping: wrong edge number");
    }
}


double
LTRSpace ::   computeEdgeVolumeAround(GaussPoint *aGaussPoint, int iEdge)
{
    double result = this->interpolation.edgeGiveTransformationJacobian(iEdge, * aGaussPoint->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this));
    return result * aGaussPoint->giveWeight();
}


void
LTRSpace ::   computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    this->interpolation.edgeLocal2global(answer, iEdge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));
}


int
LTRSpace :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    // returns transformation matrix from
    // edge local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)
    //
    _error("computeLoadLEToLRotationMatrix: egde local coordinate system not supported");
    return 1;
}

void
LTRSpace :: computeSurfaceNMatrixAt(FloatMatrix &answer, GaussPoint *sgp)
{
    FloatArray n(3);
    interpolation.surfaceEvalN(n, * sgp->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(3, 9);
    answer.zero();

    answer.at(1, 1) = n.at(1);
    answer.at(1, 4) = n.at(2);
    answer.at(1, 7) = n.at(3);

    answer.at(2, 2) = n.at(1);
    answer.at(2, 5) = n.at(2);
    answer.at(2, 8) = n.at(3);

    answer.at(3, 3) = n.at(1);
    answer.at(3, 6) = n.at(2);
    answer.at(3, 9) = n.at(3);
}

void
LTRSpace :: giveSurfaceDofMapping(IntArray &answer, int iSurf) const
{
    answer.resize(9);
    if ( iSurf == 1 ) {
        answer.at(1) = 1; // node 1
        answer.at(2) = 2;
        answer.at(3) = 3;

        answer.at(4) = 7; // node 3
        answer.at(5) = 8;
        answer.at(6) = 9;

        answer.at(7) = 4; // node 2
        answer.at(8) = 5;
        answer.at(9) = 6;
    } else if ( iSurf == 2 ) {
        answer.at(1) = 1; // node 1
        answer.at(2) = 2;
        answer.at(3) = 3;

        answer.at(4) = 4; // node 2
        answer.at(5) = 5;
        answer.at(6) = 6;

        answer.at(7) = 10; // node 4
        answer.at(8) = 11;
        answer.at(9) = 12;
    } else if ( iSurf == 3 ) {
        answer.at(1) = 4; // node 2
        answer.at(2) = 5;
        answer.at(3) = 6;

        answer.at(4) = 7; // node 3
        answer.at(5) = 8;
        answer.at(6) = 9;

        answer.at(7) = 10; // node 4
        answer.at(8) = 11;
        answer.at(9) = 12;
    } else if ( iSurf == 4 ) {
        answer.at(1) = 1; // node 1
        answer.at(2) = 2;
        answer.at(3) = 3;

        answer.at(4) = 10; // node 4
        answer.at(5) = 11;
        answer.at(6) = 12;

        answer.at(7) = 7; // node 3
        answer.at(8) = 8;
        answer.at(9) = 9;
    } else {
        _error("giveSurfaceDofMapping: wrong surface number");
    }
}


IntegrationRule *
LTRSpace :: GetSurfaceIntegrationRule(int approxOrder)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, this, 1, 1);
    int npoints = iRule->getRequiredNumberOfIntegrationPoints(_Triangle, approxOrder);
    iRule->setUpIntegrationPoints(_Triangle, npoints, _Unknown);
    return iRule;
}


double
LTRSpace :: computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
    double determinant, weight, volume;
    determinant = fabs( interpolation.surfaceGiveTransformationJacobian(iSurf, * gp->giveCoordinates(), FEIElementGeometryWrapper(this)) );

    weight      = gp->giveWeight();
    volume      = 2.0 * determinant * weight;

    return volume;

    // the following works only for 1 GP !!!
    // return interpolation.surfaceGiveTransformationJacobian (iSurf, domain, nodeArray, *gp->giveCoordinates(), 0.0);
}


void
LTRSpace :: computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int isurf)
{
    interpolation.surfaceLocal2global(answer, isurf, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));
}


int
LTRSpace :: computeLoadLSToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *)
{
    _error("computeLoadLSToLRotationMatrix: surface local coordinate system not supported");
    return 1;
}
} // end namespace oofem
