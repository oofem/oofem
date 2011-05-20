/* $Header: /home/cvs/bp/oofem/sm/src/qplanstrss.C,v 1.3.4.1 2004/04/05 15:19:47 bp Exp $ */
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

//   file QPLANSTRSS.CC

#include "qplanestrain.h"
#include "node.h"
#include "material.h"
#include "structuralcrosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "engngm.h"
#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdio.h>
#endif

namespace oofem {
FEI2dQuadQuad QPlaneStrain :: interpolation(1, 2);

QPlaneStrain :: QPlaneStrain(int n, Domain *aDomain) :
    StructuralElement(n, aDomain), ZZNodalRecoveryModelInterface()
    // Constructor.
{
    numberOfDofMans  = 8;
    numberOfGaussPoints = 4;
}

Interface *
QPlaneStrain :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return ( ZZNodalRecoveryModelInterface * ) this;
    }

    return NULL;
}


void
QPlaneStrain :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the [4x16] strain-displacement matrix {B} of the receiver,
// evaluated at aGaussPoint.
{
    int i;
    FloatMatrix dnx;

    this->interpolation.evaldNdx(dnx, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this), 0.0);

    answer.resize(4, 16);
    answer.zero();

    for ( i = 1; i <= 8; i++ ) {
        answer.at(1, 2 * i - 1) = dnx.at(i, 1);
        answer.at(2, 2 * i - 0) = dnx.at(i, 2);

        answer.at(4, 2 * i - 1) = dnx.at(i, 2);
        answer.at(4, 2 * i - 0) = dnx.at(i, 1);
    }

    return;
}

void
QPlaneStrain :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver,
// evaluated at aGaussPoint.
{
    int i;
    FloatArray n(8);

    answer.resize(2, 16);
    answer.zero();

    this->interpolation.evalN(n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this), 0.0);

    for ( i = 1; i <= 8; i++ ) {
        answer.at(1, 2 * i - 1) = n.at(i);
        answer.at(2, 2 * i - 0) = n.at(i);
    }

    return;
}

IRResultType
QPlaneStrain :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    this->Element :: initializeFrom(ir);
    numberOfGaussPoints = 4;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_QPlaneStrain_nip, "nip");

    if ( !( ( numberOfGaussPoints == 1 ) ||
           ( numberOfGaussPoints == 4 ) ||
           ( numberOfGaussPoints == 9 ) ||
           ( numberOfGaussPoints == 16 ) ) ) {
        numberOfGaussPoints = 4;
    }

    this->computeGaussPoints();
    return IRRT_OK;
}

void
QPlaneStrain :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Square, numberOfGaussPoints, _PlaneStrain);
    }
}

double
QPlaneStrain :: computeVolumeAround(GaussPoint *aGaussPoint)
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


void
QPlaneStrain ::   giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const {
    // returns DofId mask array for inode element node.
    // DofId mask array determines the dof ordering requsted from node.
    // DofId mask array contains the DofID constants (defined in cltypes.h)
    // describing physical meaning of particular DOFs.
    //IntArray* answer = new IntArray (2);
    answer.resize(2);

    answer.at(1) = D_u;
    answer.at(2) = D_v;

    return;
}


int
QPlaneStrain :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    this->interpolation.local2global(answer, lcoords, FEIElementGeometryWrapper(this), 0.0);
    return 1;
}


int
QPlaneStrain :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) ) {
        return 3;
    }

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}

void
QPlaneStrain :: ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint, InternalStateType type)
{
    // evaluates N matrix (interpolation estimated stress matrix)
    // according to Zienkiewicz & Zhu paper
    // N(nsigma, nsigma*nnodes)
    // Definition : sigmaVector = N * nodalSigmaVector

    int i;

    if ( this->giveIPValueSize(type, aGaussPoint) ) {
        answer.resize(1, 8);
    } else {
        return;
    }

    FloatArray n;
    this->interpolation.evalN(n, * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this), 0.0);

    for ( i = 1; i <= 8; i++ ) {
        answer.at(1, i) = n.at(i);
    }

    return;
}




#ifdef __OOFEG
void QPlaneStrain :: drawRawGeometry(oofegGraphicContext &gc)
{
    WCRec p [ 4 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(TRUE);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    EASValsSetFillStyle(FILL_HOLLOW);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 2 ].z = 0.;
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveCoordinate(1);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveCoordinate(2);
    p [ 3 ].z = 0.;

    go =  CreateQuad3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void QPlaneStrain :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
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
    EASValsSetEdgeFlag(TRUE);
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    EASValsSetFillStyle(FILL_HOLLOW);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 2 ].z = 0.;
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 3 ].z = 0.;

    go =  CreateQuad3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}

#endif
} // end namespace oofem
