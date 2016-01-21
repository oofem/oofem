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

#include "../sm/Elements/PlaneStress/planstrss.h"
#include "../sm/Materials/structuralms.h"
#include "fei2dquadlin.h"
#include "node.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "strainvector.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
 #include "connectivitytable.h"
 #include "Materials/rcm2.h"
#endif

namespace oofem {
REGISTER_Element(PlaneStress2d);

FEI2dQuadLin PlaneStress2d :: interpolation(1, 2);

PlaneStress2d :: PlaneStress2d(int n, Domain *aDomain) :
    PlaneStressElement(n, aDomain), ZZNodalRecoveryModelInterface(this),
    SPRNodalRecoveryModelInterface(), SpatialLocalizerInterface(this),
    HuertaErrorEstimatorInterface()
    // Constructor.
{
    numberOfDofMans  = 4;
    numberOfGaussPoints = 4;
}

PlaneStress2d :: ~PlaneStress2d()
// Destructor
{ }

FEInterpolation *PlaneStress2d :: giveInterpolation() const { return & interpolation; }

void
PlaneStress2d :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
//
// Returns the [3x8] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
// (epsilon_x,epsilon_y,gamma_xy) = B . r
// r = ( u1,v1,u2,v2,u3,v3,u4,v4)
{
    FloatMatrix dnx;

    this->interpolation.evaldNdx( dnx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );

    answer.resize(3, 8);
    answer.zero();

    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, 2 * i - 1) = dnx.at(i, 1);
        answer.at(2, 2 * i - 0) = dnx.at(i, 2);
    }

#ifdef  PlaneStress2d_reducedShearIntegration
    this->interpolation.evaldNdx( dnx, {0., 0.}, *this->giveCellGeometryWrapper() );
#endif

    for ( int i = 1; i <= 4; i++ ) {
        answer.at(3, 2 * i - 1) = dnx.at(i, 2);
        answer.at(3, 2 * i - 0) = dnx.at(i, 1);
    }
}


void
PlaneStress2d :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
//
// Returns the [4x8] displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
// @todo not checked if correct
{
    FloatMatrix dnx;

    this->interpolation.evaldNdx( dnx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );

    answer.resize(4, 8);

    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, 2 * i - 1) = dnx.at(i, 1);     // du/dx -1
        answer.at(2, 2 * i - 0) = dnx.at(i, 2);     // dv/dy -2
    }

#ifdef  PlaneStress2d_reducedShearIntegration
    this->interpolation.evaldNdx( dnx, {0., 0.}, *this->giveCellGeometryWrapper() );
#endif

    for ( int i = 1; i <= 4; i++ ) {
        answer.at(3, 2 * i - 1) = dnx.at(i, 2);     // du/dy -6
        answer.at(4, 2 * i - 0) = dnx.at(i, 1);     // dv/dx -9
    }
}


IRResultType
PlaneStress2d :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = 4;
    IRResultType result = PlaneStressElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    if ( numberOfGaussPoints != 1 && numberOfGaussPoints != 4 && numberOfGaussPoints != 9 && numberOfGaussPoints != 16 && numberOfGaussPoints != 25 ) {
        numberOfGaussPoints = 4;
        OOFEM_WARNING("Number of Gauss points enforced to 4");
    }

    return IRRT_OK;
}




double
PlaneStress2d :: giveCharacteristicSize(GaussPoint *gp, FloatArray &normalToCrackPlane, ElementCharSizeMethod method)
//
// returns receiver's characteristic size at gp (for some material models)
// for crack formed in plane with normal normalToCrackPlane
// using the selected method
//
{
    if ( method == ECSM_SquareRootOfArea ) {
        // square root of element area
        return this->computeMeanSize();
    }

    if ( method == ECSM_Projection ) {
        // standard projection method
        return this->giveCharacteristicLength(normalToCrackPlane);
    }

    // evaluate average strain and its maximum principal direction
    FloatArray sumstrain, averageNormal;
    for ( GaussPoint *gpi: *this->giveDefaultIntegrationRulePtr() ) {
        StructuralMaterialStatus *matstatus = dynamic_cast< StructuralMaterialStatus * >( gpi->giveMaterialStatus() );
        if ( matstatus ) {
            sumstrain.add( matstatus->giveTempStrainVector() );
        }
    }

    StrainVector sumstrainvec(sumstrain, _PlaneStress);
    sumstrainvec.computeMaxPrincipalDir(averageNormal);

    if ( method == ECSM_ProjectionCentered ) {
        // projection method based on principal direction of average strain
        normalToCrackPlane = averageNormal;
        return this->giveLengthInDir(normalToCrackPlane);
    }

    if ( method == ECSM_Oliver1 || method == ECSM_Oliver1modified ) {
        // method based on derivative of auxiliary function phi at each Gauss point
        // in the maximum principal strain direction determined at
        // ECSM_Oliver1 ... at each Gauss point
        // ECSM_Oliver1modified ... at element center (from average strain)

        // coordinates of the element center
        FloatArray center(2);
        double cx = 0., cy = 0.;
        for ( int i = 1; i <= 4; i++ ) {
            cx += giveNode(i)->giveCoordinate(1);
            cy += giveNode(i)->giveCoordinate(2);
        }

        cx /= 4.;
        cy /= 4.;

        // nodal values of function phi (0 or 1)
        FloatArray phi(4);
        for ( int i = 1; i <= 4; i++ ) {
            if ( ( ( giveNode(i)->giveCoordinate(1) - cx ) * averageNormal.at(1) + ( giveNode(i)->giveCoordinate(2) - cy ) * averageNormal.at(2) ) > 0. ) {
                phi.at(i) = 1.;
            } else {
                phi.at(i) = 0.;
            }
        }

        // gradient of function phi at the current GP
        FloatMatrix dnx;
        this->interpolation.evaldNdx( dnx, gp->giveNaturalCoordinates(), *this->giveCellGeometryWrapper() );
        FloatArray gradPhi(2);
        gradPhi.zero();
        for ( int i = 1; i <= 4; i++ ) {
            gradPhi.at(1) += phi.at(i) * dnx.at(i, 1);
            gradPhi.at(2) += phi.at(i) * dnx.at(i, 2);
        }

        // scalar product of the gradient with crack normal (at GP)
        double dPhidN = 0.;
        if ( method == ECSM_Oliver1modified ) {
            normalToCrackPlane = averageNormal;
        }

        for ( int i = 1; i <= 2; i++ ) {
            dPhidN += gradPhi.at(i) * normalToCrackPlane.at(i);
        }

        if ( dPhidN == 0. ) {
            OOFEM_ERROR("Zero value of dPhidN");
        }

        return 1. / fabs(dPhidN);
    }

    OOFEM_ERROR("invalid method");
    return 0.;
}



Interface *
PlaneStress2d :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >(this);
    } else if ( interface == HuertaErrorEstimatorInterfaceType ) {
        return static_cast< HuertaErrorEstimatorInterface * >(this);
    }

    return NULL;
}


void
PlaneStress2d :: HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId,
                                                                  IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                                  HuertaErrorEstimatorInterface :: SetupMode sMode, TimeStep *tStep,
                                                                  int &localNodeId, int &localElemId, int &localBcId,
                                                                  IntArray &controlNode, IntArray &controlDof,
                                                                  HuertaErrorEstimator :: AnalysisMode aMode)
{
    int inode, nodes = 4, iside, sides = 4, nd1, nd2;
    FloatArray *corner [ 4 ], midSide [ 4 ], midNode, cor [ 4 ];
    double x = 0.0, y = 0.0;

    static int sideNode [ 4 ] [ 2 ] = { { 1, 2 }, { 2, 3 }, { 3, 4 }, { 4, 1 } };

    if ( sMode == HuertaErrorEstimatorInterface :: NodeMode ||
        ( sMode == HuertaErrorEstimatorInterface :: BCMode && aMode == HuertaErrorEstimator :: HEE_linear ) ) {
        for ( inode = 0; inode < nodes; inode++ ) {
            corner [ inode ] = this->giveNode(inode + 1)->giveCoordinates();
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

    this->setupRefinedElementProblem2D(this, refinedElement, level, nodeId, localNodeIdArray, globalNodeIdArray,
                                       sMode, tStep, nodes, corner, midSide, midNode,
                                       localNodeId, localElemId, localBcId,
                                       controlNode, controlDof, aMode, "PlaneStress2d");
}

void PlaneStress2d :: HuertaErrorEstimatorI_computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    computeNmatrixAt(gp->giveSubPatchCoordinates(), answer);
}


#ifdef __OOFEG
 #define TR_LENGHT_REDUCT 0.3333

void PlaneStress2d :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
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
    EASValsSetFillStyle(FILL_HOLLOW);
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

    go =  CreateQuad3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void PlaneStress2d :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    WCRec p [ 4 ];
    GraphicObj *go;
    double defScale = gc.getDefScale();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    EASValsSetFillStyle(FILL_HOLLOW);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, defScale);
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, defScale);
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 2 ].z = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(3, tStep, defScale);
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 3 ].z = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(3, tStep, defScale);

    go =  CreateQuad3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}



void PlaneStress2d :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int i, indx, result = 0;
    WCRec p [ 4 ];
    GraphicObj *tr;
    FloatArray v [ 4 ];
    double s [ 4 ], defScale;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    if ( gc.giveIntVarMode() == ISM_recovered ) {
        for ( i = 1; i <= 4; i++ ) {
            result += this->giveInternalStateAtNode(v [ i - 1 ], gc.giveIntVarType(), gc.giveIntVarMode(), i, tStep);
        }

        if ( result != 4 ) {
            return;
        }

        indx = gc.giveIntVarIndx();

        for ( i = 1; i <= 4; i++ ) {
            s [ i - 1 ] = v [ i - 1 ].at(indx);
        }

        if ( gc.getScalarAlgo() == SA_ISO_SURF ) {
            for ( i = 0; i < 4; i++ ) {
                if ( gc.getInternalVarsDefGeoFlag() ) {
                    // use deformed geometry
                    defScale = gc.getDefScale();
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
            gc.updateFringeTableMinMax(s, 4);
            tr =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
            EGWithMaskChangeAttributes(LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);

            /*  } else if (gc.getScalarAlgo() == SA_ISO_LINE) {
             *
             * EASValsSetColor(context.getActiveCrackColor());
             * EASValsSetLineWidth(OOFEG_ISO_LINE_WIDTH);
             *
             * for (i=0; i< 4; i++) {
             * if (gc.getInternalVarsDefGeoFlag()) {
             * // use deformed geometry
             * defScale = gc.getDefScale();
             * p[i].x = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(1,tStep,defScale);
             * p[i].y = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(2,tStep,defScale);
             * p[i].z = 0.;
             *
             * } else {
             * p[i].x = (FPNum) this->giveNode(i+1)->giveCoordinate(1);
             * p[i].y = (FPNum) this->giveNode(i+1)->giveCoordinate(2);
             * p[i].z = 0.;
             * }
             * }
             *
             * // isoline implementation
             * oofeg_drawIsoLinesOnQuad (p, s);
             */
        } else if ( ( gc.getScalarAlgo() == SA_ZPROFILE ) || ( gc.getScalarAlgo() == SA_COLORZPROFILE ) ) {
            double landScale = gc.getLandScale();

            for ( i = 0; i < 4; i++ ) {
                if ( gc.getInternalVarsDefGeoFlag() ) {
                    // use deformed geometry
                    defScale = gc.getDefScale();
                    p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                    p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                    p [ i ].z = s [ i ] * landScale;
                } else {
                    p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                    p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                    p [ i ].z = s [ i ] * landScale;
                }

                // this fixes a bug in ELIXIR
                if ( fabs(s [ i ]) < 1.0e-6 ) {
                    s [ i ] = 1.0e-6;
                }
            }

            if ( gc.getScalarAlgo() == SA_ZPROFILE ) {
                EASValsSetColor( gc.getDeformedElementColor() );
                EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
                tr =  CreateQuad3D(p);
                EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
            } else {
                gc.updateFringeTableMinMax(s, 4);
                tr =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
                EGWithMaskChangeAttributes(LAYER_MASK, tr);
            }

            EMAddGraphicsToModel(ESIModel(), tr);
        }
    } else if ( gc.giveIntVarMode() == ISM_local ) {
        if ( integrationRulesArray [ 0 ]->giveNumberOfIntegrationPoints() != 4 ) {
            return;
        }

        IntArray ind(4);
        WCRec pp [ 9 ];

        for ( i = 0; i < 4; i++ ) {
            if ( gc.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = gc.getDefScale();
                pp [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                pp [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                pp [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(3, tStep, defScale);
            } else {
                pp [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                pp [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                pp [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(3);
            }
        }

        for ( i = 0; i < 3; i++ ) {
            pp [ i + 4 ].x = 0.5 * ( pp [ i ].x + pp [ i + 1 ].x );
            pp [ i + 4 ].y = 0.5 * ( pp [ i ].y + pp [ i + 1 ].y );
            pp [ i + 4 ].z = 0.5 * ( pp [ i ].z + pp [ i + 1 ].z );
        }

        pp [ 7 ].x = 0.5 * ( pp [ 3 ].x + pp [ 0 ].x );
        pp [ 7 ].y = 0.5 * ( pp [ 3 ].y + pp [ 0 ].y );
        pp [ 7 ].z = 0.5 * ( pp [ 3 ].z + pp [ 0 ].z );

        pp [ 8 ].x = 0.25 * ( pp [ 0 ].x + pp [ 1 ].x + pp [ 2 ].x + pp [ 3 ].x );
        pp [ 8 ].y = 0.25 * ( pp [ 0 ].y + pp [ 1 ].y + pp [ 2 ].y + pp [ 3 ].y );
        pp [ 8 ].z = 0.25 * ( pp [ 0 ].z + pp [ 1 ].z + pp [ 2 ].z + pp [ 3 ].z );

        for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
            const FloatArray& gpCoords = gp->giveNaturalCoordinates();
            if ( ( gpCoords.at(1) > 0. ) && ( gpCoords.at(2) > 0. ) ) {
                ind.at(1) = 0;
                ind.at(2) = 4;
                ind.at(3) = 8;
                ind.at(4) = 7;
            } else if ( ( gpCoords.at(1) < 0. ) && ( gpCoords.at(2) > 0. ) ) {
                ind.at(1) = 4;
                ind.at(2) = 1;
                ind.at(3) = 5;
                ind.at(4) = 8;
            } else if ( ( gpCoords.at(1) < 0. ) && ( gpCoords.at(2) < 0. ) ) {
                ind.at(1) = 5;
                ind.at(2) = 2;
                ind.at(3) = 6;
                ind.at(4) = 8;
            } else {
                ind.at(1) = 6;
                ind.at(2) = 3;
                ind.at(3) = 7;
                ind.at(4) = 8;
            }

            if ( giveIPValue(v [ 0 ], gp, gc.giveIntVarType(), tStep) == 0 ) {
                return;
            }

            indx = gc.giveIntVarIndx();

            for ( i = 1; i <= 4; i++ ) {
                s [ i - 1 ] = v [ 0 ].at(indx);
            }

            for ( i = 0; i < 4; i++ ) {
                p [ i ].x = pp [ ind.at(i + 1) ].x;
                p [ i ].y = pp [ ind.at(i + 1) ].y;
                p [ i ].z = pp [ ind.at(i + 1) ].z;
            }

            gc.updateFringeTableMinMax(s, 4);
            tr =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
            EGWithMaskChangeAttributes(LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);
        }
    }
}




void
PlaneStress2d :: drawSpecial(oofegGraphicContext &gc, TimeStep *tStep)
{
    int i;
    WCRec l [ 2 ];
    GraphicObj *tr;
    double defScale = gc.getDefScale();
    FloatArray crackStatuses, cf;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( gc.giveIntVarType() == IST_CrackState ) {
        // ask if any active crack exist
        int crackStatus;
        double ax, ay, bx, by, norm, xc, yc, length;
        FloatArray crackDir;
        FloatArray gpglobalcoords;

        for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {

            if ( this->giveIPValue(cf, gp, IST_CrackedFlag, tStep) == 0 ) {
                return;
            }

            if ( ( int ) cf.at(1) == 0 ) {
                return;
            }

            if ( this->giveIPValue(crackDir, gp, IST_CrackDirs, tStep) ) {
                this->giveIPValue(crackStatuses, gp, IST_CrackStatuses, tStep);
                for ( i = 1; i <= 3; i++ ) {
                    crackStatus = ( int ) crackStatuses.at(i);
                    if ( ( crackStatus != pscm_NONE ) && ( crackStatus != pscm_CLOSED ) ) {
                        // draw a crack
                        // this element is 2d element in x-y plane
                        //
                        // compute perpendicular line to normal in xy plane
                        ax = crackDir.at(i);
                        ay = crackDir.at(3 + i);
                        if ( fabs(ax) > 1.e-6 ) {
                            by = 1.;
                            bx = -ay * by / ax;
                            norm = sqrt(bx * bx + by * by);
                            bx = bx / norm; // normalize to obtain unit vector
                            by = by / norm;
                        } else {
                            by = 0.0;
                            bx = 1.0;
                        }

                        // obtain gp global coordinates
                        if ( gc.getInternalVarsDefGeoFlag() ) {
                            double ksi, eta, n1, n2, n3, n4;
                            ksi = gp->giveNaturalCoordinate(1);
                            eta = gp->giveNaturalCoordinate(2);

                            n1 = ( 1. + ksi ) * ( 1. + eta ) * 0.25;
                            n2 = ( 1. - ksi ) * ( 1. + eta ) * 0.25;
                            n3 = ( 1. - ksi ) * ( 1. - eta ) * 0.25;
                            n4 = ( 1. + ksi ) * ( 1. - eta ) * 0.25;

                            gpglobalcoords.resize(2);

                            gpglobalcoords.at(1) = ( n1 * this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale) +
                                                    n2 * this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale) +
                                                    n3 * this->giveNode(3)->giveUpdatedCoordinate(1, tStep, defScale) +
                                                    n4 * this->giveNode(4)->giveUpdatedCoordinate(1, tStep, defScale) );
                            gpglobalcoords.at(2) = ( n1 * this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale) +
                                                    n2 * this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale) +
                                                    n3 * this->giveNode(3)->giveUpdatedCoordinate(2, tStep, defScale) +
                                                    n4 * this->giveNode(4)->giveUpdatedCoordinate(2, tStep, defScale) );
                        } else {
                            computeGlobalCoordinates( gpglobalcoords, gp->giveNaturalCoordinates() );
                        }

                        xc = gpglobalcoords.at(1);
                        yc = gpglobalcoords.at(2);
                        length = sqrt( computeVolumeAround(gp) / this->giveCrossSection()->give(CS_Thickness, gp) ) / 2.;
                        l [ 0 ].x = ( FPNum ) xc + bx * length;
                        l [ 0 ].y = ( FPNum ) yc + by * length;
                        l [ 0 ].z = 0.;
                        l [ 1 ].x = ( FPNum ) xc - bx * length;
                        l [ 1 ].y = ( FPNum ) yc - by * length;
                        l [ 1 ].z = 0.;
                        EASValsSetLayer(OOFEG_CRACK_PATTERN_LAYER);
                        EASValsSetLineWidth(OOFEG_CRACK_PATTERN_WIDTH);
                        if ( ( crackStatus == pscm_SOFTENING ) || ( crackStatus == pscm_OPEN ) ) {
                            EASValsSetColor( gc.getActiveCrackColor() );
                        } else {
                            EASValsSetColor( gc.getCrackPatternColor() );
                        }

                        tr = CreateLine3D(l);
                        EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
                        EMAddGraphicsToModel(ESIModel(), tr);
                    }
                }
            }
        }
    }
}

#endif

void
PlaneStress2d :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(4);
    for ( int i = 1; i < 5; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
PlaneStress2d :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i < 5; i++ ) {
        if ( pap == this->giveNode(i)->giveNumber() ) {
            found = 1;
        }
    }

    if ( found ) {
        answer.at(1) = pap;
    } else {
        OOFEM_ERROR("node unknown");
    }
}

int
PlaneStress2d :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return this->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints();
}


SPRPatchType
PlaneStress2d :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dxy;
}

} // end namespace oofem
