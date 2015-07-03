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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#include "interfaceelem2dlin.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "fei2dlinelin.h"
#include "../sm/CrossSections/structuralinterfacecrosssection.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include <Emarkwd3d.h>
#endif

namespace oofem {
REGISTER_Element(InterfaceElem2dLin);

FEI2dLineLin InterfaceElem2dLin :: interp(1, 2);


InterfaceElem2dLin :: InterfaceElem2dLin(int n, Domain *aDomain) :
    StructuralElement(n, aDomain)
{
    numberOfDofMans       = 4;
    axisymmode            = false;
}


void
InterfaceElem2dLin :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
    FloatArray n;
    this->interp.evalN (n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(2, 8);
    answer.zero();

    answer.at(1, 2) = answer.at(2, 1) = -n.at(1);
    answer.at(1, 4) = answer.at(2, 3) = -n.at(2);

    answer.at(1, 6) = answer.at(2, 5) = n.at(1);
    answer.at(1, 8) = answer.at(2, 7) = n.at(2);
}


void
InterfaceElem2dLin :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 2) );
        integrationRulesArray [ 0 ]->SetUpPointsOnLine(2, _2dInterface); 
    }
}


double
InterfaceElem2dLin :: computeVolumeAround(GaussPoint *gp)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    double r = 1.0;
    if (this->axisymmode) {
        FloatArray n(2);
        this->interp.evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
        r = n.at(1)*this->giveNode(1)->giveCoordinate(1) + n.at(2)*this->giveNode(2)->giveCoordinate(1);
    }


    double result = this->interp.giveTransformationJacobian(gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this));
    return result * gp->giveWeight() * r;
}


void
InterfaceElem2dLin :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    static_cast< StructuralInterfaceCrossSection* >(this->giveCrossSection())->giveEngTraction_2d(answer, gp, strain, tStep);
}


void
InterfaceElem2dLin :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    static_cast< StructuralInterfaceCrossSection* >(this->giveCrossSection())->give2dStiffnessMatrix_Eng(answer, rMode, gp, tStep);
}


IRResultType
InterfaceElem2dLin :: initializeFrom(InputRecord *ir)
{
    this->axisymmode = ir->hasField(_IFT_InterfaceElem2dLin_axisymmode);
    return StructuralElement :: initializeFrom(ir);
}


void
InterfaceElem2dLin :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v};
}


bool
InterfaceElem2dLin :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    FloatArray grad(2);

    // tangent
    grad.at(1) = this->giveNode(2)->giveCoordinate(1) - this->giveNode(1)->giveCoordinate(1);
    grad.at(2) = this->giveNode(2)->giveCoordinate(2) - this->giveNode(1)->giveCoordinate(2);
    grad.normalize();

    answer.resize(8, 8);
    for ( int i = 0; i < 4; i++ ) {
        answer.at(i * 2 + 1, i * 2 + 1) = grad.at(1);
        answer.at(i * 2 + 1, i * 2 + 2) = grad.at(2);
        answer.at(i * 2 + 2, i * 2 + 1) = -grad.at(2);
        answer.at(i * 2 + 2, i * 2 + 2) = grad.at(1);
    }

    return 1;
}


FEInterpolation *
InterfaceElem2dLin :: giveInterpolation() const
{
    return & interp;
}

#ifdef __OOFEG
void InterfaceElem2dLin :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
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


void InterfaceElem2dLin :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
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


void InterfaceElem2dLin :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int indx, result = 0;
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

    for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
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
