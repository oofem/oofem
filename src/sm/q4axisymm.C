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

#include "q4axisymm.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "engngm.h"
#include "mathfem.h"
#include "crosssection.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
//#include "rcm2.h"
//#include "oofegutils.h"
#endif

namespace oofem {
REGISTER_Element(Q4Axisymm);

FEI2dQuadQuad Q4Axisymm :: interp(1, 2);

Q4Axisymm :: Q4Axisymm(int n, Domain *aDomain) :
    StructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(this)
{
    numberOfDofMans = 8;
    numberOfGaussPoints          = 4;
    numberOfFiAndShGaussPoints   = 1;
}


Q4Axisymm :: ~Q4Axisymm()
{ }


FloatArray *
Q4Axisymm :: GiveDerivativeKsi(double ksi, double eta)
{
    FloatArray *n;

    n = new FloatArray(8);

    n->at(1) =  0.25 * ( 1. + eta ) * ( 2.0 * ksi + eta );
    n->at(2) = -0.25 * ( 1. + eta ) * ( -2.0 * ksi + eta );
    n->at(3) = -0.25 * ( 1. - eta ) * ( -2.0 * ksi - eta );
    n->at(4) =  0.25 * ( 1. - eta ) * ( 2.0 * ksi - eta );
    n->at(5) = -ksi * ( 1. + eta );
    n->at(6) = -0.5 * ( 1. - eta * eta );
    n->at(7) = -ksi * ( 1. - eta );
    n->at(8) =  0.5 * ( 1. - eta * eta );

    return n;
}


FloatArray *
Q4Axisymm :: GiveDerivativeEta(double ksi, double eta)
{
    FloatArray *n;

    n = new FloatArray(8);

    n->at(1) =  0.25 * ( 1. + ksi ) * ( 2.0 * eta + ksi );
    n->at(2) =  0.25 * ( 1. - ksi ) * ( 2.0 * eta - ksi );
    n->at(3) = -0.25 * ( 1. - ksi ) * ( -2.0 * eta - ksi );
    n->at(4) = -0.25 * ( 1. + ksi ) * ( -2.0 * eta + ksi );
    n->at(5) =  0.5 * ( 1. - ksi * ksi );
    n->at(6) = -eta * ( 1. - ksi );
    n->at(7) = -0.5 * ( 1. - ksi * ksi );
    n->at(8) = -eta * ( 1. + ksi );

    return n;
}


void
Q4Axisymm :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [6x16] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    ///@todo Not sure how to deal with li and lu. This should be changed, in which case "GiveDerivatice***" and "computeJacobianMatrixAt" will be deprecated.
    //FloatMatrix dN;
    //this->interp.evaldNdx(dN, *gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this));

    int i;
    FloatMatrix jacMtrx, inv;
    FloatArray *nx, *ny;
    double ksi, eta, r, x;
    int size, ind = 1;

    ksi = gp->giveNaturalCoordinate(1);
    eta = gp->giveNaturalCoordinate(2);

    nx = GiveDerivativeKsi(ksi, eta);
    ny = GiveDerivativeEta(ksi, eta);

    this->computeJacobianMatrixAt(jacMtrx, gp);
    inv.beInverseOf(jacMtrx);

    if ( ui == ALL_STRAINS ) {
        size = 6;
        ui = 6;
    } else {
        size = ui - li + 1;
    }

    if ( ( size < 0 ) || ( size > 6 ) ) {
        OOFEM_ERROR("size mismatch");
    }

    answer.resize(size, 16);
    answer.zero();

    if ( ( li <= 1 ) && ( ui >= 1 ) ) {
        for ( i = 1; i <= 8; i++ ) {
            answer.at(ind, 2 * i - 1) = nx->at(i) * inv.at(1, 1) + ny->at(i) * inv.at(1, 2);
        }

        ind++;
    }

    if ( ( li <= 2 ) && ( ui >= 2 ) ) {
        for ( i = 1; i <= 8; i++ ) {
            answer.at(ind, 2 * i - 0) = nx->at(i) * inv.at(2, 1) + ny->at(i) * inv.at(2, 2);
        }

        ind++;
    }

    if ( ( li <= 3 ) && ( ui >= 3 ) ) {
        FloatArray n(8);

        n.at(1) = ( 1. + ksi ) * ( 1. + eta ) * 0.25 * ( ksi + eta - 1. );
        n.at(2) = ( 1. - ksi ) * ( 1. + eta ) * 0.25 * ( -ksi + eta - 1. );
        n.at(3) = ( 1. - ksi ) * ( 1. - eta ) * 0.25 * ( -ksi - eta - 1. );
        n.at(4) = ( 1. + ksi ) * ( 1. - eta ) * 0.25 * ( ksi - eta - 1. );
        n.at(5) = 0.5 * ( 1. - ksi * ksi ) * ( 1. + eta );
        n.at(6) = 0.5 * ( 1. - ksi ) * ( 1. - eta * eta );
        n.at(7) = 0.5 * ( 1. - ksi * ksi ) * ( 1. - eta );
        n.at(8) = 0.5 * ( 1. + ksi ) * ( 1. - eta * eta );

        r = 0.;
        for ( i = 1; i <= 8; i++ ) {
            x  = this->giveNode(i)->giveCoordinate(1);
            r += x * n.at(i);
        }

        for ( i = 1; i <= 8; i++ ) {
            answer.at( ind, ( 2 * i - 1 ) ) = ( n.at(i) ) / r;
        }

        // delete n;
        ind++;
    }

    if ( ( li <= 4 ) && ( ui >= 4 ) ) {
        ind++;
    }

    if ( ( li <= 5 ) && ( ui >= 5 ) ) {
        ind++;
    }


    if ( ( li <= 6 ) && ( ui >= 6 ) ) {
        for ( i = 1; i <= 8; i++ ) {
            answer.at(ind, 2 * i - 1) = nx->at(i) * inv.at(2, 1) + ny->at(i) * inv.at(2, 2);
            answer.at(ind, 2 * i - 0) = nx->at(i) * inv.at(1, 1) + ny->at(i) * inv.at(1, 2);
        }

        ind++;
    }

    delete nx;
    delete ny;
}


void
Q4Axisymm :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the [9x16] displacement gradient matrix {BH} of the receiver,
// evaluated at gp.
// BH matrix  -  9 rows : du/dx, dv/dy, dw/dz = u/r, 0, 0, du/dy,  0, 0, dv/dx
// @todo not checked if correct
{
    FloatArray n;
    FloatMatrix dnx;

    this->interp.evaldNdx( dnx, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(9, 16);
    answer.zero();

    double r = 0., x;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        x  = this->giveNode(i)->giveCoordinate(1);
        r += x * n.at(i);
    }


    // mode is _3dMat !!!!!! answer.at(4,*), answer.at(5,*), answer.at(7,*), and answer.at(8,*) is zero
    for ( int i = 1; i <= 16; i++ ) {
        answer.at(1, 3 * i - 2) = dnx.at(i, 1);     // du/dx
        answer.at(2, 3 * i - 1) = dnx.at(i, 2);     // dv/dy
        answer.at(6, 3 * i - 2) = dnx.at(i, 2);     // du/dy
        answer.at(9, 3 * i - 1) = dnx.at(i, 1);     // dv/dx
    }

    answer.at(3, 1) = n.at(1) / r;
    answer.at(3, 3) = n.at(2) / r;
    answer.at(3, 5) = n.at(3) / r;
    answer.at(3, 7) = n.at(4) / r;
    answer.at(3, 9) = n.at(5) / r;
    answer.at(3, 11) = n.at(6) / r;
    answer.at(3, 13) = n.at(7) / r;
    answer.at(3, 15) = n.at(8) / r;
}


void
Q4Axisymm :: computeJacobianMatrixAt(FloatMatrix &answer, GaussPoint *gp)
// Returns the jacobian matrix  J (x,y)/(ksi,eta)  of the receiver.
{
    double ksi, eta, x, y;
    FloatArray *nx, *ny;

    answer.resize(2, 2);
    answer.zero();

    ksi = gp->giveNaturalCoordinate(1);
    eta = gp->giveNaturalCoordinate(2);

    nx = this->GiveDerivativeKsi(ksi, eta);
    ny = this->GiveDerivativeEta(ksi, eta);

    for ( int i = 1; i <= 8; i++ ) {
        x = this->giveNode(i)->giveCoordinate(1);
        y = this->giveNode(i)->giveCoordinate(2);

        answer.at(1, 1) += nx->at(i) * x;
        answer.at(1, 2) += nx->at(i) * y;
        answer.at(2, 1) += ny->at(i) * x;
        answer.at(2, 2) += ny->at(i) * y;
    }

    delete nx;
    delete ny;
}


IRResultType
Q4Axisymm :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    numberOfGaussPoints          = 4;
    result = this->StructuralElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    numberOfFiAndShGaussPoints   = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfFiAndShGaussPoints, _IFT_Q4Axisymm_nipfish);

    if ( !( ( numberOfGaussPoints == 1 ) ||
           ( numberOfGaussPoints == 4 ) ||
           ( numberOfGaussPoints == 9 ) ||
           ( numberOfGaussPoints == 16 ) ) ) {
        numberOfGaussPoints = 4;
    }

    if ( !( ( numberOfFiAndShGaussPoints == 1 ) ||
           ( numberOfFiAndShGaussPoints == 4 ) ||
           ( numberOfFiAndShGaussPoints == 9 ) ||
           ( numberOfFiAndShGaussPoints == 16 ) ) ) {
        numberOfFiAndShGaussPoints = 1;
    }


    return IRRT_OK;
}


void
Q4Axisymm :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(2);
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 2);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);

        integrationRulesArray [ 1 ] = new GaussIntegrationRule(2, this, 3, 6);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 1 ], numberOfFiAndShGaussPoints, this);
    }
}


double
Q4Axisymm :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    FloatArray n;
    double determinant, r;

    this->interp.evalN( n, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    r = 0;
    for ( int i = 1; i <= 8; i++ ) {
        r += this->giveNode(i)->giveCoordinate(1) * n.at(i);
    }

    determinant = fabs( this->interp.giveTransformationJacobian( * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    return determinant *gp->giveWeight() * r;
}


void
Q4Axisymm :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step tStep. The nature of these strains depends
// on the element's type.
{
    FloatMatrix b, A;
    FloatArray u, Epsilon, help;
    fMode mode = domain->giveEngngModel()->giveFormulation();

    answer.resize(6);
    answer.zero();
    if ( mode == TL ) { // Total Lagrange formulation
        this->computeVectorOf(VM_Total, tStep, u);

        // linear part of strain tensor (in vector form)

        this->computeBmatrixAt(gp, b, 1, 2);
        Epsilon.beProductOf(b, u);
        answer.at(1) = Epsilon.at(1);
        answer.at(2) = Epsilon.at(2);
        // delete Epsilon;  delete b;

        if ( numberOfFiAndShGaussPoints == 1 ) {
            //
            // if reduced integration in one gp only
            // force the evaluation of eps_fi in this gauss point
            // instead of evaluating in given gp
            //
            GaussPoint *helpGaussPoint;
            helpGaussPoint = integrationRulesArray [ 1 ]->getIntegrationPoint(0);

            this->computeBmatrixAt(helpGaussPoint, b, 3, 6);
        } else {
            this->computeBmatrixAt(gp, b, 3, 6);
            //_error("ComputeStrainVector: numberOfFiAndShGaussPoints size mismatch");
        }

        Epsilon.beProductOf(b, u);
        answer.at(3) = Epsilon.at(1);
        answer.at(6) = Epsilon.at(4);
    } else if ( mode == AL ) { // actualized Lagrange formulation
        OOFEM_ERROR("unsupported mode");
    }
}

void
Q4Axisymm :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v};
}



Interface *
Q4Axisymm :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    }

    return NULL;
}



#ifdef __OOFEG

void Q4Axisymm :: drawRawGeometry(oofegGraphicContext &gc)
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
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void Q4Axisymm :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
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
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 2 ].z = 0.;
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(2, tStep, defScale);
    p [ 3 ].z = 0.;

    go =  CreateQuad3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}

#endif
} // end namespace oofem
