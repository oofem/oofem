/* $Header: /home/cvs/bp/oofem/sm/src/qtrplstr.C,v 1.5 2003/04/06 14:08:31 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   file QTRPLANSTRSS.CC


#include "qtrplstr.h"
#include "node.h"
#include "material.h"
#include "gausspnt.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "verbose.h"
#include "engngm.h"
#include "crosssection.h"
#include "gaussintegrationrule.h"
#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
 #include "oofegutils.h"
 #include "rcm2.h"
#endif

namespace oofem {
FEI2dTrQuad QTrPlaneStress2d :: interpolation(1, 2);

QTrPlaneStress2d :: QTrPlaneStress2d(int n, Domain *aDomain) :
    StructuralElement(n, aDomain), SpatialLocalizerInterface(),
    DirectErrorIndicatorRCInterface(), EIPrimaryUnknownMapperInterface()
    // Constructor.
{
    numberOfDofMans  = 6;
    numberOfGaussPoints = 4;
}

Interface *
QTrPlaneStress2d :: giveInterface(InterfaceType interface)
{
    //if (interface == NodalAveragingRecoveryModelInterfaceType) return (NodalAveragingRecoveryModelInterface*) this;
    if (interface == ZZNodalRecoveryModelInterfaceType){
        return (ZZNodalRecoveryModelInterface*) this;
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return ( SPRNodalRecoveryModelInterface * ) this;
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return ( SpatialLocalizerInterface * ) this;
    } else if ( interface == DirectErrorIndicatorRCInterfaceType ) {
        return ( DirectErrorIndicatorRCInterface * ) this;
    } else if ( interface == EIPrimaryUnknownMapperInterfaceType ) {
        return ( EIPrimaryUnknownMapperInterface * ) this;
    }

    return NULL;
}




/*
 * void
 * QTrPlaneStress2d :: computeBmatrixAt (GaussPoint *aGaussPoint, FloatMatrix& answer, int li, int ui)
 * // Returns the [3x12] strain-displacement matrix {B} of the receiver, eva-
 * // luated at aGaussPoint.
 * {
 * // FloatMatrix *answer;
 * double x1,x2,x3,y1,y2,y3,b1,b2,b3,c1,c2,c3,area,l1,l2,l3;
 *
 * x1 = this -> giveNode(1) -> giveCoordinate(1);
 * x2 = this -> giveNode(2) -> giveCoordinate(1);
 * x3 = this -> giveNode(3) -> giveCoordinate(1);
 *
 * y1 = this -> giveNode(1) -> giveCoordinate(2);
 * y2 = this -> giveNode(2) -> giveCoordinate(2);
 * y3 = this -> giveNode(3) -> giveCoordinate(2);
 *
 * area = 0.5*(x2*y3+x1*y2+y1*x3-x2*y1-x3*y2-x1*y3);
 *
 * b1 = y2-y3;
 * b2 = y3-y1;
 * b3 = y1-y2;
 *
 * c1 = x3-x2;
 * c2 = x1-x3;
 * c3 = x2-x1;
 *
 * l1 = aGaussPoint -> giveCoordinate(1);
 * l2 = aGaussPoint -> giveCoordinate(2);
 * l3 = 1.-l1-l2;
 *
 * //answer = new FloatMatrix(3,12);
 * answer.resize (3,12);
 * answer.zero();
 *
 * answer.at(1,1)  = 2.*b1*l1+(2.*l1-1.)*b1;
 * answer.at(1,3)  = 2.*b2*l2+(2.*l2-1.)*b2;
 * answer.at(1,5)  = 2.*b3*l3+(2.*l3-1.)*b3;
 * answer.at(1,7)  = 4.*b1*l2+4.*l1*b2;
 * answer.at(1,9)  = 4.*b2*l3+4.*l2*b3;
 * answer.at(1,11) = 4.*b3*l1+4.*l3*b1;
 *
 * answer.at(2,2)  = 2.*c1*l1+(2.*l1-1.)*c1;
 * answer.at(2,4)  = 2.*c2*l2+(2.*l2-1.)*c2;
 * answer.at(2,6)  = 2.*c3*l3+(2.*l3-1.)*c3;
 * answer.at(2,8)  = 4.*c1*l2+4.*l1*c2;
 * answer.at(2,10) = 4.*c2*l3+4.*l2*c3;
 * answer.at(2,12) = 4.*c3*l1+4.*l3*c1;
 *
 * answer.at(3,1)  = 2.*c1*l1+(2.*l1-1.)*c1;
 * answer.at(3,3)  = 2.*c2*l2+(2.*l2-1.)*c2;
 * answer.at(3,5)  = 2.*c3*l3+(2.*l3-1.)*c3;
 * answer.at(3,7)  = 4.*c1*l2+4.*l1*c2;
 * answer.at(3,9)  = 4.*c2*l3+4.*l2*c3;
 * answer.at(3,11) = 4.*c3*l1+4.*l3*c1;
 *
 * answer.at(3,2)  = 2.*b1*l1+(2.*l1-1.)*b1;
 * answer.at(3,4)  = 2.*b2*l2+(2.*l2-1.)*b2;
 * answer.at(3,6)  = 2.*b3*l3+(2.*l3-1.)*b3;
 * answer.at(3,8)  = 4.*b1*l2+4.*l1*b2;
 * answer.at(3,10) = 4.*b2*l3+4.*l2*b3;
 * answer.at(3,12) = 4.*b3*l1+4.*l3*b1;
 *
 * answer.times(1./(2.*area));
 * return ;
 * }
 */

void
QTrPlaneStress2d :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint.
{
    int i;
    FloatArray n(6);

    answer.resize(2, 12);
    answer.zero();

    this->interpolation.evalN(n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this), 0.0);

    for ( i = 1; i <= 6; i++ ) {
        answer.at(1, 2 * i - 1) = n.at(i);
        answer.at(2, 2 * i - 0) = n.at(i);
    }

    return;
}

IRResultType
QTrPlaneStress2d :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    this->StructuralElement :: initializeFrom(ir);
    numberOfGaussPoints = 4;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_QTrPlaneStress2d_nip, "nip"); // Macro

    if ( !( ( numberOfGaussPoints == 1 ) ||
           ( numberOfGaussPoints == 4 ) ||
           ( numberOfGaussPoints == 7 ) ) ) {
        numberOfGaussPoints = 4;
    }

    this->computeGaussPoints();
    return IRRT_OK;
}



void
QTrPlaneStress2d :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the [3x12] strain-displacement matrix {B} of the receiver, eva-
// luated at aGaussPoint.
{
    int i;
    FloatMatrix dnx;

    this->interpolation.evaldNdx(dnx, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this), 0.0);

    answer.resize(3, 12);
    answer.zero();

    for ( i = 1; i <= 6; i++ ) {
        answer.at(1, 2 * i - 1) = dnx.at(i, 1);
        answer.at(2, 2 * i - 0) = dnx.at(i, 2);

        answer.at(3, 2 * i - 1) = dnx.at(i, 2);
        answer.at(3, 2 * i - 0) = dnx.at(i, 1);
    }

    return;
}

double
QTrPlaneStress2d :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
    double determinant, weight, thickness, volume;
    determinant = fabs( this->interpolation.giveTransformationJacobian(* aGaussPoint->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this), 0.0) );
    weight      = aGaussPoint->giveWeight();
    thickness   = this->giveCrossSection()->give(CS_Thickness);
    volume      = determinant * weight * thickness;

    return volume;
}


int
QTrPlaneStress2d :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    this->interpolation.local2global(answer, lcoords, FEIElementGeometryWrapper(this), 0.0);
    return 1;
}


#define POINT_TOL 1.e-3

int
QTrPlaneStress2d :: computeLocalCoordinates(FloatArray &answer, const FloatArray &coords)
{
    return this->interpolation.global2local(answer, coords, FEIElementGeometryWrapper(this), 0.0);
}

void QTrPlaneStress2d :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Triangle, numberOfGaussPoints, _PlaneStress);
    }
}

void
QTrPlaneStress2d ::   giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const {
    // returns DofId mask array for inode element node.
    // DofId mask array determines the dof ordering requsted from node.
    // DofId mask array contains the DofID constants (defined in cltypes.h)
    // describing physical meaning of particular DOFs.

    // IntArray* answer = new IntArray (2);
    answer.resize(2);

    answer.at(1) = D_u;
    answer.at(2) = D_v;

    return;
}



int
QTrPlaneStress2d :: SpatialLocalizerI_containsPoint(const FloatArray &coords) {
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}

double
QTrPlaneStress2d :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
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

void QTrPlaneStress2d :: drawRawGeometry(oofegGraphicContext &gc)
{
    WCRec p [ 3 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(TRUE);
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


void QTrPlaneStress2d :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
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
    EASValsSetEdgeFlag(TRUE);
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


void QTrPlaneStress2d :: drawScalar(oofegGraphicContext &context)
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

    this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
    if ( ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
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
                    p [ i ].x = ( FPNum ) this->giveNode(n [ i ])->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
                    p [ i ].y = ( FPNum ) this->giveNode(n [ i ])->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
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
         * p[i].x = (FPNum) this->giveNode(n[i])->giveUpdatedCoordinate(1,tStep,EID_MomentumBalance,defScale);
         * p[i].y = (FPNum) this->giveNode(n[i])->giveUpdatedCoordinate(2,tStep,EID_MomentumBalance,defScale);
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
QTrPlaneStress2d :: drawSpecial(oofegGraphicContext &gc)
{
    return;
}

#endif




int
QTrPlaneStress2d :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
     if ((type == IST_StressTensor)||(type ==IST_StrainTensor) || ( type == IST_DamageTensor )) return 3;

     GaussPoint *gp = integrationRulesArray[0]-> getIntegrationPoint(0) ;
     return this->giveIPValueSize (type, gp);
}


void
QTrPlaneStress2d :: ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx (FloatMatrix& answer, GaussPoint* aGaussPoint, InternalStateType type)
 {
 // evaluates N matrix (interpolation estimated stress matrix)
 // according to Zienkiewicz & Zhu paper
 // N(nsigma, nsigma*nnodes)
 // Definition : sigmaVector = N * nodalSigmaVector
     int i;
     FloatArray n;

     this->interpolation.evalN(n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this), 0.0);

     if (this->giveIPValueSize(type, aGaussPoint)){
         answer.resize(1,6) ;
     } else {
         return;
     }
     
     for ( i = 1; i <= 6; i++ ) {
        answer.at(1, i)  = n.at(i);
     }

//  double l1,l2,l3;
//  l1 = aGaussPoint -> giveCoordinate(1);
//  l2 = aGaussPoint -> giveCoordinate(2);
//  l3 = 1.0 - l1 - l2;
// 
//  if (this->giveIPValueSize(type, aGaussPoint)){
//      answer.resize(1,6) ;
//  } else {
//      return;
//  }
// 
//  answer.at(1,1)  = (2.*l1-1.)*l1;
//  answer.at(1,2)  = (2.*l2-1.)*l2;
//  answer.at(1,3)  = (2.*l3-1.)*l3;
//  answer.at(1,4)  = 4.*l1*l2;
//  answer.at(1,5)  = 4.*l2*l3;
//  answer.at(1,6)  = 4.*l3*l1;

 return ;
 }
 

int
QTrPlaneStress2d :: SPRNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) ) {
        return 3;
    }

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}



void
QTrPlaneStress2d :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(3);
    pap.at(1) = this->giveNode(1)->giveNumber();
    pap.at(2) = this->giveNode(2)->giveNumber();
    pap.at(3) = this->giveNode(3)->giveNumber();
}

void
QTrPlaneStress2d :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
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
QTrPlaneStress2d :: SPRNodalRecoveryMI_giveNumberOfIP()
{ return numberOfGaussPoints; }


void
QTrPlaneStress2d :: SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp)
{
    if ( this->computeGlobalCoordinates( coords, * gp->giveCoordinates() ) == 0 ) {
        _error("SPRNodalRecoveryMI_computeIPGlobalCoordinates: computeGlobalCoordinates failed");
    }
}

SPRPatchType
QTrPlaneStress2d :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dquadratic;
}



double
QTrPlaneStress2d :: DirectErrorIndicatorRCI_giveCharacteristicSize() {
    IntegrationRule *iRule = this->giveDefaultIntegrationRulePtr();
    GaussPoint *gp;
    double volume = 0.0;

    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp  = iRule->getIntegrationPoint(i);
        volume += this->computeVolumeAround(gp);
    }

    return sqrt( volume * 2.0 / this->giveCrossSection()->give(CS_Thickness) );
}

int
QTrPlaneStress2d :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode,
                                                                     TimeStep *stepN, const FloatArray &coords,
                                                                     FloatArray &answer)
{
    FloatArray lcoords, u, nn;
    FloatMatrix n(2, 12);
    int i, result;

    result = this->computeLocalCoordinates(lcoords, coords);

    this->interpolation.evalN(nn, lcoords, FEIElementGeometryWrapper(this), 0.0);

    for ( i = 1; i <= 6; i++ ) {
        n.at(1, 2 * i - 1) = nn.at(i);
        n.at(2, 2 * i - 0) = nn.at(i);
    }

    this->computeVectorOf(EID_MomentumBalance, mode, stepN, u);
    answer.beProductOf(n, u);

    return result;
}

void
QTrPlaneStress2d :: EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer)
{
    giveDofManDofIDMask(1, EID_MomentumBalance, answer);
}
} // end namespace oofem
