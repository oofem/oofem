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

#include "interfaceelem2dquad.h"
#include "node.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "mathfem.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #ifndef __MAKEDEPEND
  #include <Emarkwd3d.h>
 #endif
#endif

namespace oofem {
InterfaceElem2dQuad :: InterfaceElem2dQuad(int n, Domain *aDomain) :
    StructuralElement(n, aDomain)
{
    numberOfDofMans       = 6;
}


void
InterfaceElem2dQuad :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
//
// Returns linear part of geometrical equations of the receiver at gp.
// Returns the linear part of the B matrix
//
{
    double ksi, n1, n2, n3;

    ksi = aGaussPoint->giveCoordinate(1);
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
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        //integrationRulesArray[0] = new LobattoIntegrationRule (1,domain, 1, 2);
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 2);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Line, 4, _2dInterface);
    }
}


int
InterfaceElem2dQuad :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    double ksi, n1, n2, n3;

    ksi = lcoords.at(1);
    n3  = 1. - ksi * ksi;
    n1  = ( 1. - ksi ) * 0.5 - 0.5 * n3;
    n2  = ( 1. + ksi ) * 0.5 - 0.5 * n3;


    answer.resize(2);
    answer.at(1) = n1 * this->giveNode(1)->giveCoordinate(1) + n2 *this->giveNode(2)->giveCoordinate(1) + n3 *this->giveNode(3)->giveCoordinate(1);
    answer.at(2) = n1 * this->giveNode(1)->giveCoordinate(2) + n2 *this->giveNode(2)->giveCoordinate(2) + n3 *this->giveNode(3)->giveCoordinate(2);

    return 1;
}


int
InterfaceElem2dQuad :: computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords)
{
    _error("Not implemented");
    return 0;
}


double
InterfaceElem2dQuad :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    double weight  = aGaussPoint->giveWeight();
    double ksi = aGaussPoint->giveCoordinate(1);
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
    double thickness   = this->giveCrossSection()->give(CS_Thickness);
    return sqrt(dx * dx + dy * dy) * weight * thickness;
}


IRResultType
InterfaceElem2dQuad :: initializeFrom(InputRecord *ir)
{
    this->StructuralElement :: initializeFrom(ir);
    this->computeGaussPoints();
    return IRRT_OK;
}


void
InterfaceElem2dQuad :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(2, D_u, D_v);
}


bool
InterfaceElem2dQuad :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    int i;
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
    for ( i = 0; i < 6; i++ ) {
        answer.at(i * 2 + 1, i * 2 + 1) = grad.at(1);
        answer.at(i * 2 + 1, i * 2 + 2) = grad.at(2);
        answer.at(i * 2 + 2, i * 2 + 1) = -grad.at(2);
        answer.at(i * 2 + 2, i * 2 + 2) = grad.at(1);
    }

    return 1;
}

///@todo Deprecated? Is so, remove it. / Mikael
/*
 * void
 * InterfaceElem2dQuad :: computeStiffnessMatrix (FloatMatrix& answer, MatResponseMode rMode,
 *                                             TimeStep* tStep)
 * // Computes numerically the stiffness matrix of the receiver.
 * {
 * int         j ;
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
 * for (j=0 ; j < iRule->getNumberOfIntegrationPoints() ; j++) {
 *  gp = iRule->getIntegrationPoint(j) ;
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
void InterfaceElem2dQuad :: drawRawGeometry(oofegGraphicContext &gc)
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


void InterfaceElem2dQuad :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
{
    GraphicObj *go;
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* poin */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    double defScale = gc.getDefScale();

    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER + 1);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].z = 0.0;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].z = 0.0;
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);

    p [ 0 ].x = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].z = 0.0;
    p [ 1 ].x = ( FPNum ) this->giveNode(5)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(5)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].z = 0.0;
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


void InterfaceElem2dQuad :: drawScalar(oofegGraphicContext &context)
{
    int i, indx, result = 0;
    GaussPoint *gp;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    FloatArray gcoord(3), v1;
    WCRec p [ 1 ];
    IntArray map;
    GraphicObj *go;
    double val [ 1 ];

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    if ( context.giveIntVarMode() == ISM_recovered ) {
        return;
    }

    for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        result = 0;
        gp  = iRule->getIntegrationPoint(i);
        result += giveIPValue(v1, gp, context.giveIntVarType(), tStep);
        result += this->giveIntVarCompFullIndx( map, context.giveIntVarType() );
        if ( result != 2 ) {
            continue;
        }

        if ( ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
            return;
        }

        result += this->computeGlobalCoordinates( gcoord, * ( gp->giveCoordinates() ) );

        p [ 0 ].x = ( FPNum ) gcoord.at(1);
        p [ 0 ].y = ( FPNum ) gcoord.at(2);
        p [ 0 ].z = 0.;

        val [ 0 ] = v1.at(indx);
        context.updateFringeTableMinMax(val, 1);
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
