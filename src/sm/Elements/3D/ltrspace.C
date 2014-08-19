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

#include "../sm/Elements/3D/ltrspace.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "fei3dtetlin.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
 #include "rcm2.h"

 #include <Etetrawd.h>
#endif

namespace oofem {
REGISTER_Element(LTRSpace);

FEI3dTetLin LTRSpace :: interpolation;

LTRSpace :: LTRSpace(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(this), NodalAveragingRecoveryModelInterface(),
    SPRNodalRecoveryModelInterface(), SpatialLocalizerInterface(this),
    EIPrimaryUnknownMapperInterface(), ZZErrorEstimatorInterface(this),
    MMAShapeFunctProjectionInterface(), HuertaErrorEstimatorInterface()

{
    numberOfDofMans  = 4;
    numberOfGaussPoints = 1;
}


Interface *
LTRSpace :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >(this);
    } else if ( interface == EIPrimaryUnknownMapperInterfaceType ) {
        return static_cast< EIPrimaryUnknownMapperInterface * >(this);
    } else if ( interface == ZZErrorEstimatorInterfaceType ) {
        return static_cast< ZZErrorEstimatorInterface * >(this);
    } else if ( interface == MMAShapeFunctProjectionInterfaceType ) {
        return static_cast< MMAShapeFunctProjectionInterface * >(this);
    } else if ( interface == HuertaErrorEstimatorInterfaceType ) {
        return static_cast< HuertaErrorEstimatorInterface * >(this);
    }

    return NULL;
}


FEInterpolation *
LTRSpace :: giveInterpolation() const
{
    return & interpolation;
}


void
LTRSpace :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [6x12] strain-displacement matrix {B} of the receiver, eva-
// luated at gp.
{
    FloatMatrix dn(4, 3);
    interpolation.evaldNdx( dn, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(6, 12);
    answer.zero();

    for ( int i = 1; i <= 4; i++ ) {
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
LTRSpace :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the [9x24] displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
// BH matrix  -  9 rows : du/dx, dv/dy, dw/dz, dv/dz, du/dz, du/dy, dw/dy, dw/dx, dv/dx
//@todo unchecked
{
    FloatMatrix dnx;

    this->interpolation.evaldNdx( dnx, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(9, 12);
    answer.zero();

    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, 3 * i - 2) = dnx.at(i, 1);     // du/dx
        answer.at(2, 3 * i - 1) = dnx.at(i, 2);     // dv/dy
        answer.at(3, 3 * i - 0) = dnx.at(i, 3);     // dw/dz
        answer.at(4, 3 * i - 1) = dnx.at(i, 3);     // dv/dz
        answer.at(7, 3 * i - 0) = dnx.at(i, 2);     // dw/dy
        answer.at(5, 3 * i - 2) = dnx.at(i, 3);     // du/dz
        answer.at(8, 3 * i - 0) = dnx.at(i, 1);     // dw/dx
        answer.at(6, 3 * i - 2) = dnx.at(i, 2);     // du/dy
        answer.at(9, 3 * i - 1) = dnx.at(i, 1);     // dv/dx
    }

#if 0
    // test if sym(BH) = H*BH == Bsym
    FloatMatrix H, Bsym, Btest;
    H.resize(6, 9);
    H.at(1, 1) = H.at(2, 2) = H.at(3, 3) = H.at(4, 4) = H.at(4, 7) = H.at(5, 5) = H.at(5, 8) = H.at(6, 6) = H.at(6, 9) = 1.0;
    Btest.beProductOf(H, answer);
    computeBmatrixAt(gp, Bsym);
    Btest.printYourself();
    Bsym.printYourself();
#endif
}


double LTRSpace :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double determinant, weight, volume;
    determinant = fabs( this->interpolation.giveTransformationJacobian( * gp->giveNaturalCoordinates(),
                                                                       FEIElementGeometryWrapper(this) ) );
    weight = gp->giveWeight();
    volume = determinant * weight;
    return volume;
}

IRResultType
LTRSpace :: initializeFrom(InputRecord *ir)
{
    //IRResultType result;                            // Required by IR_GIVE_FIELD macro

    numberOfGaussPoints = 1;
    this->NLStructuralElement :: initializeFrom(ir);
    return IRRT_OK;
}


MaterialMode
LTRSpace :: giveMaterialMode()
{
    return _3dMat;
}


void LTRSpace :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 6);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
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
    double density = this->giveStructuralCrossSection()->give('d', gp);
    mss1 = dV * density / 4.;

    for ( int i = 1; i <= 12; i++ ) {
        answer.at(i, i) = mss1;
    }
}


void
LTRSpace :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v, D_w};
}


double
LTRSpace :: giveCharacteristicLength(const FloatArray &normalToCrackPlane)
{
    return this->giveLengthInDir(normalToCrackPlane);
}


void
LTRSpace :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                       InternalStateType type, TimeStep *tStep)
{
    GaussPoint *gp;

    if ( numberOfGaussPoints != 1 ) {
        answer.clear(); // for more gp's need to be refined
        return;
    }

    gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    giveIPValue(answer, gp, type, tStep);
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
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i <= 4; i++ ) {
        if ( this->giveNode(i)->giveNumber() == pap ) {
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
LTRSpace :: SPRNodalRecoveryMI_giveNumberOfIP()
{ return this->giveDefaultIntegrationRulePtr()->giveNumberOfIntegrationPoints(); }


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
            corner [ inode ] = this->giveNode(inode + 1)->giveCoordinates();

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

    this->setupRefinedElementProblem3D(this, refinedElement, level, nodeId, localNodeIdArray, globalNodeIdArray,
                                       sMode, tStep, nodes, corner, midSide, midFace, midNode,
                                       localNodeId, localElemId, localBcId, hexaSideNode, hexaFaceNode,
                                       controlNode, controlDof, aMode, "LSpace");
}

void LTRSpace :: HuertaErrorEstimatorI_computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    computeNmatrixAt(* ( gp->giveSubPatchCoordinates() ), answer);
}

#ifdef __OOFEG
 #define TR_LENGHT_REDUCT 0.3333

void LTRSpace :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
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


void LTRSpace :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
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
    EASValsSetFillStyle(FILL_SOLID);
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

    go =  CreateTetra(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


void LTRSpace :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int i, indx, result = 0;
    WCRec p [ 4 ];
    GraphicObj *tr;
    FloatArray v [ 4 ];
    double s [ 4 ], defScale = 0.0;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( gc.giveIntVarMode() == ISM_recovered ) {
        for ( i = 1; i <= 4; i++ ) {
            result += this->giveInternalStateAtNode(v [ i - 1 ], gc.giveIntVarType(), gc.giveIntVarMode(), i, tStep);
        }

        if ( result != 4 ) {
            return;
        }
    } else if ( gc.giveIntVarMode() == ISM_local ) {
        return;
    }

    indx = gc.giveIntVarIndx();

    for ( i = 1; i <= 4; i++ ) {
        s [ i - 1 ] = v [ i - 1 ].at(indx);
    }

    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    if ( gc.getScalarAlgo() == SA_ISO_SURF ) {
        for ( i = 0; i < 4; i++ ) {
            if ( gc.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(3, tStep, defScale);
            } else {
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(3);
            }
        }

        gc.updateFringeTableMinMax(s, 4);
        tr = CreateTetraWD(p, s);
        EGWithMaskChangeAttributes(LAYER_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    }
}

void
LTRSpace :: drawSpecial(oofegGraphicContext &gc, TimeStep *tStep)
{
    int i, j, k;
    WCRec q [ 4 ];
    GraphicObj *tr;
    IntegrationRule *iRule = integrationRulesArray [ 0 ];
    GaussPoint *gp;
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

        //   for (GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        {
            gp = iRule->getIntegrationPoint(0);
            if ( this->giveIPValue(cf, gp, IST_CrackedFlag, tStep) == 0 ) {
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
                    xc += ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                    yc += ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                    zc += ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(3, tStep, defScale);
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
            if ( this->giveIPValue(crackDir, gp, IST_CrackDirs, tStep) ) {
                this->giveIPValue(crackStatuses, gp, IST_CrackStatuses, tStep);


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


double
LTRSpace :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray lcoords(4), gcoords;
    double dist;
    int size, gsize;

    lcoords.at(1) = lcoords.at(2) = lcoords.at(3) = lcoords.at(4) = 1. / 4.;
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


void
LTRSpace :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType mode,
                                                             TimeStep *tStep, const FloatArray &lcoords,
                                                             FloatArray &answer)
{
    FloatArray u;
    FloatMatrix n;
    n.beNMatrixOf(lcoords, 3);

    this->computeVectorOf(mode, tStep, u);
    answer.beProductOf(n, u);
}


void
LTRSpace :: MMAShapeFunctProjectionInterface_interpolateIntVarAt(FloatArray &answer, FloatArray &coords,
                                                                 coordType ct, nodalValContainerType &list,
                                                                 InternalStateType type, TimeStep *tStep)
{
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
    answer.resize(0);
    answer.add(l1, list[0]);
    answer.add(l2, list[1]);
    answer.add(l3, list[2]);
    answer.add(l4, list[3]);
}


void
LTRSpace :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
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
    this->interpolation.edgeEvalN( n, iedge, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

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
        OOFEM_ERROR("wrong edge number");
    }
}


double
LTRSpace :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double result = this->interpolation.edgeGiveTransformationJacobian( iEdge, * gp->giveNaturalCoordinates(),
                                                                       FEIElementGeometryWrapper(this) );
    return result *gp->giveWeight();
}


void
LTRSpace :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    this->interpolation.edgeLocal2global( answer, iEdge, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
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
    OOFEM_ERROR("egde local coordinate system not supported");
    return 1;
}

void
LTRSpace :: computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *sgp)
{
    FloatArray n(3);
    interpolation.surfaceEvalN( n, iSurf, * sgp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

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
        OOFEM_ERROR("wrong surface number");
    }
}


IntegrationRule *
LTRSpace :: GetSurfaceIntegrationRule(int approxOrder)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, this, 1, 1);
    int npoints = iRule->getRequiredNumberOfIntegrationPoints(_Triangle, approxOrder);
    iRule->SetUpPointsOnTriangle(npoints, _Unknown);
    return iRule;
}


double
LTRSpace :: computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
    double determinant, weight, volume;
    determinant = fabs( interpolation.surfaceGiveTransformationJacobian( iSurf, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );

    weight      = gp->giveWeight();
    volume      = determinant * weight;

    return volume;

    // the following works only for 1 GP !!!
    // return interpolation.surfaceGiveTransformationJacobian (iSurf, domain, nodeArray, *gp->giveNaturalCoordinates(), 0.0);
}


void
LTRSpace :: computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int isurf)
{
    interpolation.surfaceLocal2global( answer, isurf, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
}


int
LTRSpace :: computeLoadLSToLRotationMatrix(FloatMatrix &answer, int, GaussPoint *)
{
    OOFEM_ERROR("surface local coordinate system not supported");
    return 1;
}
} // end namespace oofem
