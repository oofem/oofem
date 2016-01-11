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

#include "../sm/Elements/Interfaces/intelline2.h"
#include "../sm/CrossSections/structuralinterfacecrosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "fei2dlinequad.h"
#include "fei2dlinelin.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include <Emarkwd3d.h>
#endif

namespace oofem {
REGISTER_Element(IntElLine2);

FEI2dLineQuad IntElLine2 :: interp(2, 2);
FEI2dLineLin IntElLine2 :: interpLin(1, 1);


IntElLine2 :: IntElLine2(int n, Domain *aDomain) : IntElLine1(n, aDomain)
{
    numberOfDofMans = 6;

    numberOfGaussPoints = 4;
    linear = false;
}


void
IntElLine2 :: computeNmatrixAt(GaussPoint *ip, FloatMatrix &answer)
{

    // Returns the modified N-matrix which multiplied with u give the spatial jump.
    FloatArray N;
    answer.resize(2, 12);
    answer.zero();

    if(linear) {
		interpLin.evalN( N, ip->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

		answer.at(1, 1) = answer.at(2, 2) = -N.at(1);
		answer.at(1, 3) = answer.at(2, 4) = -N.at(2);
//		answer.at(1, 5) = answer.at(2, 6) = -N.at(3);

		answer.at(1, 7) = answer.at(2, 8) = N.at(1);
		answer.at(1, 9) = answer.at(2, 10) = N.at(2);
//		answer.at(1, 11) = answer.at(2, 12) = N.at(3);
    }
    else {
		interp.evalN( N, ip->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

		answer.at(1, 1) = answer.at(2, 2) = -N.at(1);
		answer.at(1, 3) = answer.at(2, 4) = -N.at(2);
		answer.at(1, 5) = answer.at(2, 6) = -N.at(3);

		answer.at(1, 7) = answer.at(2, 8) = N.at(1);
		answer.at(1, 9) = answer.at(2, 10) = N.at(2);
		answer.at(1, 11) = answer.at(2, 12) = N.at(3);
    }
}


void
IntElLine2 :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        //integrationRulesArray[ 0 ].reset( new LobattoIntegrationRule (1,domain, 1, 2) ); ///@todo - should be able to decide
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 2) );
        integrationRulesArray [ 0 ]->SetUpPointsOnLine(numberOfGaussPoints, _2dInterface); ///@todo - should be a parameter with num of ip
    }
}

FEInterpolation *
IntElLine2 :: giveInterpolation() const
{
    return & interp;
}


IRResultType
IntElLine2 :: initializeFrom(InputRecord *ir)
{
	linear = ir->hasField(_IFT_IntElLine2_LinearTraction);
    return IntElLine1 :: initializeFrom(ir);   
}



#ifdef __OOFEG
void IntElLine2 :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
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
    p [ 1 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 1 ].z = 0.0;
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
    p [ 0 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 0 ].z = 0.0;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = 0.0;
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void IntElLine2 :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
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

    p [ 0 ].x = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = 0.0;
    p [ 1 ].x = ( FPNum ) this->giveNode(5)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(5)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = 0.0;
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


void IntElLine2 :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
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
