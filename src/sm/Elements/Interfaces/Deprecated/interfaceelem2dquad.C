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

#include "interfaceelem2dquad.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "fei2dlinequad.h"
#include "../sm/CrossSections/structuralinterfacecrosssection.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include <Emarkwd3d.h>
#endif

namespace oofem {
REGISTER_Element(InterfaceElem2dQuad);

FEI2dLineQuad InterfaceElem2dQuad :: interp(1, 2);


InterfaceElem2dQuad :: InterfaceElem2dQuad(int n, Domain *aDomain) :
    StructuralElement(n, aDomain)
{
    numberOfDofMans       = 6;
    axisymmode            = false;
}


void
InterfaceElem2dQuad :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
    ///@todo Use the interpolator everywhere in this file:
    double ksi, n1, n2, n3;

    ksi = gp->giveNaturalCoordinate(1);
    n3  = 1. - ksi * ksi;
    n1  = ( 1. - ksi ) * 0.5 - 0.5 * n3;
    n2  = ( 1. + ksi ) * 0.5 - 0.5 * n3;
    answer.resize(2, 12);
    answer.zero();

    answer.at(1, 2) = answer.at(2, 1) = -n1;
    answer.at(1, 4) = answer.at(2, 3) = -n2;
    answer.at(1, 6) = answer.at(2, 5) = -n3;

    answer.at(1, 8) = answer.at(2, 7) = n1;
    answer.at(1, 10) = answer.at(2, 9) = n2;
    answer.at(1, 12) = answer.at(2, 11) = n3;
}


void
InterfaceElem2dQuad :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        //integrationRulesArray[0].reset( new LobattoIntegrationRule (1,domain, 1, 2) );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 2) );
        integrationRulesArray [ 0 ]->SetUpPointsOnLine(4, _2dInterface);
    }
}


double
InterfaceElem2dQuad :: computeVolumeAround(GaussPoint *gp)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    double weight  = gp->giveWeight();
    double ksi = gp->giveNaturalCoordinate(1);
    double dn1 = ksi - 0.5;
    double dn2 = ksi + 0.5;
    double dn3 = -2.0 * ksi;

    double x1 = this->giveNode(1)->giveCoordinate(1);
    double x2 = this->giveNode(2)->giveCoordinate(1);
    double x3 = this->giveNode(3)->giveCoordinate(1);

    double y1 = this->giveNode(1)->giveCoordinate(2);
    double y2 = this->giveNode(2)->giveCoordinate(2);
    double y3 = this->giveNode(3)->giveCoordinate(2);

    double dx = ( dn1 * x1 ) + ( dn2 * x2 ) + ( dn3 * x3 );
    double dy = ( dn1 * y1 ) + ( dn2 * y2 ) + ( dn3 * y3 );
    double thickness   = this->giveCrossSection()->give(CS_Thickness, gp);

    double r = 1.0;
    if (this->axisymmode) {
      double n3  = 1. - ksi * ksi;
      double n1  = ( 1. - ksi ) * 0.5 - 0.5 * n3;
      double n2  = ( 1. + ksi ) * 0.5 - 0.5 * n3;
      r = n1*this->giveNode(1)->giveCoordinate(1) + n2*this->giveNode(2)->giveCoordinate(1)+ n3*this->giveNode(3)->giveCoordinate(1);
    }

    return sqrt(dx * dx + dy * dy) * weight * thickness * r;
}


void
InterfaceElem2dQuad :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    static_cast< StructuralInterfaceCrossSection* >(this->giveCrossSection())->giveEngTraction_2d(answer, gp, strain, tStep);
}


void
InterfaceElem2dQuad :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    static_cast< StructuralInterfaceCrossSection* >(this->giveCrossSection())->give2dStiffnessMatrix_Eng(answer, rMode, gp, tStep);
}


IRResultType
InterfaceElem2dQuad :: initializeFrom(InputRecord *ir)
{
    this->axisymmode = ir->hasField(_IFT_InterfaceElem2dQuad_axisymmode);
    return StructuralElement :: initializeFrom(ir);
}


void
InterfaceElem2dQuad :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v};
}


bool
InterfaceElem2dQuad :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    FloatArray grad(2);

    //double ksi = gp -> giveCoordinate(1) ;
    double ksi = 0.0; // compute tangent in the middle
    double dn1 = ksi - 0.5;
    double dn2 = ksi + 0.5;
    double dn3 = -2.0 * ksi;

    // tangent
    grad.at(1) = dn1 * this->giveNode(1)->giveCoordinate(1) + dn2 *this->giveNode(2)->giveCoordinate(1) + dn3 *this->giveNode(3)->giveCoordinate(1);
    grad.at(2) = dn1 * this->giveNode(1)->giveCoordinate(2) + dn2 *this->giveNode(2)->giveCoordinate(2) + dn3 *this->giveNode(3)->giveCoordinate(2);
    grad.normalize();

    answer.resize(12, 12);
    for ( int i = 0; i < 6; i++ ) {
        answer.at(i * 2 + 1, i * 2 + 1) = grad.at(1);
        answer.at(i * 2 + 1, i * 2 + 2) = grad.at(2);
        answer.at(i * 2 + 2, i * 2 + 1) = -grad.at(2);
        answer.at(i * 2 + 2, i * 2 + 2) = grad.at(1);
    }

    return 1;
}


FEInterpolation *
InterfaceElem2dQuad :: giveInterpolation() const
{
    return & interp;
}

///@todo Deprecated? Is so, remove it. / Mikael
/*
 * void
 * InterfaceElem2dQuad :: computeStiffnessMatrix (FloatMatrix& answer, MatResponseMode rMode,
 *                                             TimeStep* tStep)
 * // Computes numerically the stiffness matrix of the receiver.
 * {
 * double      dV ;
 * FloatMatrix d, bj, bjl, dbj, t;
 * GaussPoint  *gp ;
 * IntegrationRule* iRule;
 * bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode, this->material);
 *
 * answer.resize (computeNumberOfDofs(),computeNumberOfDofs());
 * answer.zero();
 *
 * iRule = integrationRulesArray[giveDefaultIntegrationRule()];
 *
 * for ( GaussPoint *gp: *iRule ) {
 *  this -> computeBmatrixAt(gp, bjl) ;
 *  this -> computeConstitutiveMatrixAt(d, rMode, gp, tStep);
 *  this -> computeGtoLRotationMatrix(t, gp);
 *  bj.beProductOf(bjl,t);
 *  dbj.beProductOf (d, bj) ;
 *  dV = this->computeVolumeAround (gp);
 *  if (matStiffSymmFlag) answer.plusProductSymmUpper (bj,dbj,dV); else answer.plusProductUnsym (bj,dbj,dV);
 *
 * }
 * }
 */


#ifdef __OOFEG
void InterfaceElem2dQuad :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
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


void InterfaceElem2dQuad :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
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


void InterfaceElem2dQuad :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
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
