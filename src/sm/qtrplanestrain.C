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

#include "qtrplanestrain.h"
#include "node.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "crosssection.h"
#include "gaussintegrationrule.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
 #include "oofegutils.h"
 #include "rcm2.h"
#endif


namespace oofem {

REGISTER_Element( QTrPlaneStrain );

FEI2dTrQuad QTrPlaneStrain :: interpolation(1, 2);

QTrPlaneStrain :: QTrPlaneStrain(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain), SpatialLocalizerInterface(),
    DirectErrorIndicatorRCInterface(), EIPrimaryUnknownMapperInterface()
    // Constructor.
{
    numberOfDofMans = 6;
    numberOfGaussPoints = 4;
}


Interface *
QTrPlaneStrain :: giveInterface(InterfaceType interface)
{
    //if (interface == NodalAveragingRecoveryModelInterfaceType) return (NodalAveragingRecoveryModelInterface*) this;
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >( this );
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >( this );
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >( this );
    } else if ( interface == DirectErrorIndicatorRCInterfaceType ) {
        return static_cast< DirectErrorIndicatorRCInterface * >( this );
    } else if ( interface == EIPrimaryUnknownMapperInterfaceType ) {
        return static_cast< EIPrimaryUnknownMapperInterface * >( this );
    }

    return NULL;
}


void
QTrPlaneStrain :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
{
    FloatArray n(6);

    answer.resize(2, 12);
    answer.zero();

    this->interpolation.evalN( n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this) );

    for ( int i = 1; i <= 6; i++ ) {
        answer.at(1, 2 * i - 1) = n.at(i);
        answer.at(2, 2 * i - 0) = n.at(i);
    }
}

IRResultType
QTrPlaneStrain :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    this->StructuralElement :: initializeFrom(ir);
    numberOfGaussPoints = 4;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, _IFT_Element_nip);

    return IRRT_OK;
}



void
QTrPlaneStrain :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the [3x12] strain-displacement matrix {B} of the receiver, eva-
// luated at aGaussPoint.
{
    FloatMatrix dnx;

    this->interpolation.evaldNdx( dnx, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(4, 12);
    answer.zero();

    for ( int i = 1; i <= 6; i++ ) {
        answer.at(1, 2 * i - 1) = dnx.at(i, 1);
        answer.at(2, 2 * i - 0) = dnx.at(i, 2);

        answer.at(4, 2 * i - 1) = dnx.at(i, 2);
        answer.at(4, 2 * i - 0) = dnx.at(i, 1);
    }
}

double
QTrPlaneStrain :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
    double determinant, weight, thickness;
    determinant = fabs( this->interpolation.giveTransformationJacobian( * aGaussPoint->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this) ) );
    weight = aGaussPoint->giveWeight();
    thickness = this->giveCrossSection()->give(CS_Thickness);
    return determinant * weight * thickness;
}

void QTrPlaneStrain :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        this->giveCrossSection()->setupIntegrationPoints( *integrationRulesArray[0], numberOfGaussPoints, this );
    }
}

void
QTrPlaneStrain :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(2, D_u, D_v);
}


int
QTrPlaneStrain :: SpatialLocalizerI_containsPoint(const FloatArray &coords)
{
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}


double
QTrPlaneStrain :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
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

        helpCoords.resizeWithValues(gsize);
        dist = helpCoords.distance(gcoords);
    }

    return dist;
}


#ifdef __OOFEG
 #define TR_LENGHT_REDUCT 0.3333

void QTrPlaneStrain :: drawRawGeometry(oofegGraphicContext &gc)
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


void QTrPlaneStrain :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
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


void QTrPlaneStrain :: drawScalar(oofegGraphicContext &context)
{
    int t, n [ 3 ], i, indx, result = 0;
    WCRec p [ 3 ];
    GraphicObj *tr;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    FloatArray v [ 6 ];
    double s [ 6 ], ss [ 3 ], defScale;
    IntArray map;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    if ( context.giveIntVarMode() == ISM_recovered ) {
        for ( i = 1; i <= 6; i++ ) {
            result += this->giveInternalStateAtNode(v [ i - 1 ], context.giveIntVarType(), context.giveIntVarMode(), i, tStep);
        }
    } else if ( context.giveIntVarMode() == ISM_local ) {
        return;
    }

    if ( result != 6 ) {
        return;
    }

    result = this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
    if ( ( !result ) || ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
        return;
    }

    for ( i = 1; i <= 6; i++ ) {
        s [ i - 1 ] = v [ i - 1 ].at(indx);
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);

    if ( context.getScalarAlgo() == SA_ISO_SURF ) {
        for ( t = 1; t <= 4; t++ ) {
            if ( t == 1 ) {
                n [ 0 ] = 1;
                n [ 1 ] = 4;
                n [ 2 ] = 6;
            } else if ( t == 2 ) {
                n [ 0 ] = 2;
                n [ 1 ] = 5;
                n [ 2 ] = 4;
            } else if ( t == 3 ) {
                n [ 0 ] = 3;
                n [ 1 ] = 6;
                n [ 2 ] = 5;
            } else {
                n [ 0 ] = 4;
                n [ 1 ] = 5;
                n [ 2 ] = 6;
            }


            for ( i = 0; i < 3; i++ ) {
                if ( context.getInternalVarsDefGeoFlag() ) {
                    // use deformed geometry
                    defScale = context.getDefScale();
                    p [ i ].x = ( FPNum ) this->giveNode(n [ i ])->giveUpdatedCoordinate(1, tStep, defScale);
                    p [ i ].y = ( FPNum ) this->giveNode(n [ i ])->giveUpdatedCoordinate(2, tStep, defScale);
                    p [ i ].z = 0.;
                } else {
                    p [ i ].x = ( FPNum ) this->giveNode(n [ i ])->giveCoordinate(1);
                    p [ i ].y = ( FPNum ) this->giveNode(n [ i ])->giveCoordinate(2);
                    p [ i ].z = 0.;
                }
            }

            //EASValsSetColor(gc.getYieldPlotColor(ratio));
            ss [ 0 ] = s [ n [ 0 ] - 1 ];
            ss [ 1 ] = s [ n [ 1 ] - 1 ];
            ss [ 2 ] = s [ n [ 2 ] - 1 ];
            context.updateFringeTableMinMax(ss, 3);
            tr =  CreateTriangleWD3D(p, ss [ 0 ], ss [ 1 ], ss [ 2 ]);
            EGWithMaskChangeAttributes(LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);
        }

        /* } else if (context.getScalarAlgo() == SA_ISO_LINE) {
         *
         * EASValsSetColor(context.getActiveCrackColor());
         * EASValsSetLineWidth(OOFEG_ISO_LINE_WIDTH);
         *
         * for (t=1; t<=4; t++) {
         * if (t==1) {n[0] = 1; n[1]=4; n[2]=6;}
         * else if (t==2) {n[0]=2; n[1]=5; n[2]=4;}
         * else if (t==3) {n[0]=3; n[1]=6; n[2]=5;}
         * else {n[0]=4; n[1]=5; n[2]=6;}
         *
         *
         * for (i=0; i< 3; i++) {
         * if (context.getInternalVarsDefGeoFlag()) {
         * // use deformed geometry
         * defScale = context.getDefScale();
         * p[i].x = (FPNum) this->giveNode(n[i])->giveUpdatedCoordinate(1,tStep,defScale);
         * p[i].y = (FPNum) this->giveNode(n[i])->giveUpdatedCoordinate(2,tStep,defScale);
         * p[i].z = 0.;
         *
         * } else {
         * p[i].x = (FPNum) this->giveNode(n[i])->giveCoordinate(1);
         * p[i].y = (FPNum) this->giveNode(n[i])->giveCoordinate(2);
         * p[i].z = 0.;
         * }
         * }
         * sv[0]=s[n[0]-1];
         * sv[1]=s[n[1]-1];
         * sv[2]=s[n[2]-1];
         *
         * // isoline implementation
         * oofeg_drawIsoLinesOnTriangle (p, sv);
         * } */
    }
}

void
QTrPlaneStrain :: drawSpecial(oofegGraphicContext &gc)
{ }

#endif




int
QTrPlaneStrain :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) || ( type == IST_DamageTensor ) ) {
        return 3;
    }

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}


int
QTrPlaneStrain :: SPRNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) ) {
        return 3;
    }

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}



void
QTrPlaneStrain :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(3);
    pap.at(1) = this->giveNode(1)->giveNumber();
    pap.at(2) = this->giveNode(2)->giveNumber();
    pap.at(3) = this->giveNode(3)->giveNumber();
}

void
QTrPlaneStrain :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    answer.resize(3);
    if ( pap == this->giveNode(1)->giveNumber() ) {
        answer.at(1) = pap;
        answer.at(2) = this->giveNode(4)->giveNumber();
        answer.at(3) = this->giveNode(6)->giveNumber();
    } else if ( pap == this->giveNode(2)->giveNumber() ) {
        answer.at(1) = pap;
        answer.at(2) = this->giveNode(5)->giveNumber();
        answer.at(3) = this->giveNode(4)->giveNumber();
    } else if ( pap == this->giveNode(3)->giveNumber() ) {
        answer.at(1) = pap;
        answer.at(2) = this->giveNode(6)->giveNumber();
        answer.at(3) = this->giveNode(5)->giveNumber();
    } else {
        _error("SPRNodalRecoveryMI_giveDofMansDeterminedByPatch: node unknown");
    }
}

int
QTrPlaneStrain :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return numberOfGaussPoints;
}


SPRPatchType
QTrPlaneStrain :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dquadratic;
}


double
QTrPlaneStrain :: DirectErrorIndicatorRCI_giveCharacteristicSize()
{
    IntegrationRule *iRule = this->giveDefaultIntegrationRulePtr();
    GaussPoint *gp;
    double volume = 0.0;

    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        gp  = iRule->getIntegrationPoint(i);
        volume += this->computeVolumeAround(gp);
    }

    return sqrt( volume * 2.0 / this->giveCrossSection()->give(CS_Thickness) );
}

int
QTrPlaneStrain :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode,
                                                                   TimeStep *stepN, const FloatArray &coords,
                                                                   FloatArray &answer)
{
    FloatArray lcoords, u, nn;
    FloatMatrix n(2, 12);
    int result;

    result = this->computeLocalCoordinates(lcoords, coords);

    this->interpolation.evalN( nn, lcoords, FEIElementGeometryWrapper(this) );

    for ( int i = 1; i <= 6; i++ ) {
        n.at(1, 2 * i - 1) = nn.at(i);
        n.at(2, 2 * i - 0) = nn.at(i);
    }

    this->computeVectorOf(EID_MomentumBalance, mode, stepN, u);
    answer.beProductOf(n, u);

    return result;
}

void
QTrPlaneStrain :: EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer)
{
    giveDofManDofIDMask(1, EID_MomentumBalance, answer);
}


double
QTrPlaneStrain :: giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane)
{
    return this->giveLenghtInDir(normalToCrackPlane) / sqrt((double) gp->giveIntegrationRule()->giveNumberOfIntegrationPoints());
}


} // end namespace oofem
