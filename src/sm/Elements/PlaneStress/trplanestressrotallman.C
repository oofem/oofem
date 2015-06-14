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

#include "../sm/Elements/PlaneStress/trplanestressrotallman.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "fei2dtrquad.h"
#include "fei2dtrlin.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
#endif

namespace oofem {
REGISTER_Element(TrPlanestressRotAllman);

FEI2dTrQuad TrPlanestressRotAllman :: qinterpolation(1, 2);

TrPlanestressRotAllman :: TrPlanestressRotAllman(int n, Domain *aDomain) :
    TrPlaneStress2d(n, aDomain)
{
    numberOfDofMans  = 3;
    numberOfGaussPoints = 4;
}

Interface *
TrPlanestressRotAllman :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >(this);
    }
    return NULL;
}

void
TrPlanestressRotAllman :: computeLocalNodalCoordinates(std::vector< FloatArray > &lxy)
{
    lxy.resize(6);
    for ( int i = 0; i < 3; i++ ) {
        lxy [ i ] = * this->giveNode(i + 1)->giveCoordinates();
    }
    lxy [ 3 ].resize(2);
    lxy [ 4 ].resize(2);
    lxy [ 5 ].resize(2);
    for ( int i = 1; i <= 2; i++ ) {
        lxy [ 3 ].at(i) = 0.5 * ( lxy [ 0 ].at(i) + lxy [ 1 ].at(i) );
        lxy [ 4 ].at(i) = 0.5 * ( lxy [ 1 ].at(i) + lxy [ 2 ].at(i) );
        lxy [ 5 ].at(i) = 0.5 * ( lxy [ 2 ].at(i) + lxy [ 0 ].at(i) );
    }
}



void
TrPlanestressRotAllman :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp.
{
    FloatArray L, n;
    std::vector< FloatArray > lxy;

    answer.resize(3, 9);
    answer.zero();

    this->computeLocalNodalCoordinates(lxy); // get ready for tranformation into 3d
    this->qinterpolation.evalN( n, iLocCoord, FEIVertexListGeometryWrapper(lxy) );
    this->interp.evalN( L, iLocCoord, FEIVertexListGeometryWrapper(lxy));

    answer.at(1, 1) = answer.at(2, 2) = n.at(1) + n.at(4) / 2. + n.at(6) / 2.;
    answer.at(1, 4) = answer.at(2, 5) = n.at(2) + n.at(4) / 2. + n.at(5) / 2.;
    answer.at(1, 7) = answer.at(2, 8) = n.at(3) + n.at(5) / 2. + n.at(6) / 2.;
    answer.at(1, 3) = n.at(6) * ( lxy [ 0 ].at(2) - lxy [ 2 ].at(2) ) / 8.0 - n.at(4) * ( lxy [ 1 ].at(2) - lxy [ 0 ].at(2) ) / 8.0;
    answer.at(1, 6) = n.at(4) * ( lxy [ 1 ].at(2) - lxy [ 0 ].at(2) ) / 8.0 - n.at(5) * ( lxy [ 2 ].at(2) - lxy [ 1 ].at(2) ) / 8.0;
    answer.at(1, 9) = n.at(5) * ( lxy [ 2 ].at(2) - lxy [ 1 ].at(2) ) / 8.0 - n.at(6) * ( lxy [ 0 ].at(2) - lxy [ 2 ].at(2) ) / 8.0;
    answer.at(2, 3) = -n.at(6) * ( lxy [ 0 ].at(1) - lxy [ 2 ].at(1) ) / 8.0 + n.at(4) * ( lxy [ 1 ].at(1) - lxy [ 0 ].at(1) ) / 8.0;
    answer.at(2, 6) = -n.at(4) * ( lxy [ 1 ].at(1) - lxy [ 0 ].at(1) ) / 8.0 + n.at(5) * ( lxy [ 2 ].at(1) - lxy [ 1 ].at(1) ) / 8.0;
    answer.at(2, 9) = -n.at(5) * ( lxy [ 2 ].at(1) - lxy [ 1 ].at(1) ) / 8.0 + n.at(6) * ( lxy [ 0 ].at(1) - lxy [ 2 ].at(1) ) / 8.0;
    // linear approx for rotations
    answer.at(3, 3) = L.at(1);
    answer.at(3, 6) = L.at(2);
    answer.at(3, 9) = L.at(3);
}

void
TrPlanestressRotAllman :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [3x12] strain-displacement matrix {B} of the receiver, eva-
// luated at gp.
{
    FloatMatrix dnx;
    std::vector< FloatArray > lxy;

    this->computeLocalNodalCoordinates(lxy); // get ready for tranformation into 3d
    this->qinterpolation.evaldNdx( dnx, gp->giveNaturalCoordinates(), FEIVertexListGeometryWrapper(lxy) );

    answer.resize(3, 9);
    answer.zero();

    // epsilon_x
    answer.at(1, 1) = dnx.at(1, 1) + 0.5 * dnx.at(4, 1) + 0.5 * dnx.at(6, 1);
    answer.at(1, 4) = dnx.at(2, 1) + 0.5 * dnx.at(4, 1) + 0.5 * dnx.at(5, 1);
    answer.at(1, 7) = dnx.at(3, 1) + 0.5 * dnx.at(5, 1) + 0.5 * dnx.at(6, 1);
    answer.at(1, 3) =+dnx.at(6, 1) * ( lxy [ 0 ].at(2) - lxy [ 2 ].at(2) ) / 8.0 - dnx.at(4, 1) * ( lxy [ 1 ].at(2) - lxy [ 0 ].at(2) ) / 8.0;
    answer.at(1, 6) =+dnx.at(4, 1) * ( lxy [ 1 ].at(2) - lxy [ 0 ].at(2) ) / 8.0 - dnx.at(5, 1) * ( lxy [ 2 ].at(2) - lxy [ 1 ].at(2) ) / 8.0;
    answer.at(1, 9) =+dnx.at(5, 1) * ( lxy [ 2 ].at(2) - lxy [ 1 ].at(2) ) / 8.0 - dnx.at(6, 1) * ( lxy [ 0 ].at(2) - lxy [ 2 ].at(2) ) / 8.0;

    // epsilon_y
    answer.at(2, 2) = dnx.at(1, 2) + 0.5 * dnx.at(4, 2) + 0.5 * dnx.at(6, 2);
    answer.at(2, 5) = dnx.at(2, 2) + 0.5 * dnx.at(4, 2) + 0.5 * dnx.at(5, 2);
    answer.at(2, 8) = dnx.at(3, 2) + 0.5 * dnx.at(5, 2) + 0.5 * dnx.at(6, 2);
    answer.at(2, 3) =-dnx.at(6, 2) * ( lxy [ 0 ].at(1) - lxy [ 2 ].at(1) ) / 8.0 + dnx.at(4, 2) * ( lxy [ 1 ].at(1) - lxy [ 0 ].at(1) ) / 8.0;
    answer.at(2, 6) =-dnx.at(4, 2) * ( lxy [ 1 ].at(1) - lxy [ 0 ].at(1) ) / 8.0 + dnx.at(5, 2) * ( lxy [ 2 ].at(1) - lxy [ 1 ].at(1) ) / 8.0;
    answer.at(2, 9) =-dnx.at(5, 2) * ( lxy [ 2 ].at(1) - lxy [ 1 ].at(1) ) / 8.0 + dnx.at(6, 2) * ( lxy [ 0 ].at(1) - lxy [ 2 ].at(1) ) / 8.0;

    // gamma_xy (shear)
    answer.at(3, 1) = dnx.at(1, 2) + 0.5 * dnx.at(4, 2) + 0.5 * dnx.at(6, 2);
    answer.at(3, 2) = dnx.at(1, 1) + 0.5 * dnx.at(4, 1) + 0.5 * dnx.at(6, 1);
    answer.at(3, 4) = dnx.at(2, 2) + 0.5 * dnx.at(4, 2) + 0.5 * dnx.at(5, 2);
    answer.at(3, 5) = dnx.at(2, 1) + 0.5 * dnx.at(4, 1) + 0.5 * dnx.at(5, 1);
    answer.at(3, 7) = dnx.at(3, 2) + 0.5 * dnx.at(5, 2) + 0.5 * dnx.at(6, 2);
    answer.at(3, 8) = dnx.at(3, 1) + 0.5 * dnx.at(5, 1) + 0.5 * dnx.at(6, 1);

    answer.at(3, 3) = dnx.at(6, 2) * ( lxy [ 0 ].at(2) - lxy [ 2 ].at(2) ) / 8.0 - dnx.at(4, 2) * ( lxy [ 1 ].at(2) - lxy [ 0 ].at(2) ) / 8.0;
    answer.at(3, 3) += -dnx.at(6, 1) * ( lxy [ 0 ].at(1) - lxy [ 2 ].at(1) ) / 8.0 + dnx.at(4, 1) * ( lxy [ 1 ].at(1) - lxy [ 0 ].at(1) ) / 8.0;
    answer.at(3, 6) = dnx.at(4, 2) * ( lxy [ 1 ].at(2) - lxy [ 0 ].at(2) ) / 8.0 - dnx.at(5, 2) * ( lxy [ 2 ].at(2) - lxy [ 1 ].at(2) ) / 8.0;
    answer.at(3, 6) += -dnx.at(4, 1) * ( lxy [ 1 ].at(1) - lxy [ 0 ].at(1) ) / 8.0 + dnx.at(5, 1) * ( lxy [ 2 ].at(1) - lxy [ 1 ].at(1) ) / 8.0;
    answer.at(3, 9) = dnx.at(5, 2) * ( lxy [ 2 ].at(2) - lxy [ 1 ].at(2) ) / 8.0 - dnx.at(6, 2) * ( lxy [ 0 ].at(2) - lxy [ 2 ].at(2) ) / 8.0;
    answer.at(3, 9) += -dnx.at(5, 1) * ( lxy [ 2 ].at(1) - lxy [ 1 ].at(1) ) / 8.0 + dnx.at(6, 1) * ( lxy [ 0 ].at(1) - lxy [ 2 ].at(1) ) / 8.0;
}


void
TrPlanestressRotAllman :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    // compute standard stiffness matrix
    TrPlaneStress2d :: computeStiffnessMatrix(answer, rMode, tStep);
    // add zero energy mode stabilization
    FloatMatrix ks;
    this->computeStiffnessMatrixZeroEnergyStabilization(ks, rMode, tStep);
    answer.add(ks);
}

void
TrPlanestressRotAllman :: computeStiffnessMatrixZeroEnergyStabilization(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    FloatMatrix b(1, 9);
    FloatMatrix dnx;
    FloatArray lec = {0.333333333333, 0.333333333333, 0.333333333333}; // element center in local coordinates
    std::vector< FloatArray > lxy; 

    this->computeLocalNodalCoordinates(lxy); // get ready for tranformation into 3d
    this->qinterpolation.evaldNdx( dnx, lec, FEIVertexListGeometryWrapper(lxy) );

    // evaluate (dv/dx-du/dy)/2. at element center
    b.at(1, 1) = -1.0 * ( dnx.at(1, 2) + 0.5 * dnx.at(4, 2) + 0.5 * dnx.at(6, 2) );
    b.at(1, 2) = dnx.at(1, 1) + 0.5 * dnx.at(4, 1) + 0.5 * dnx.at(6, 1);
    b.at(1, 4) = -1.0 * ( dnx.at(2, 2) + 0.5 * dnx.at(4, 2) + 0.5 * dnx.at(5, 2) );
    b.at(1, 5) = dnx.at(2, 1) + 0.5 * dnx.at(4, 1) + 0.5 * dnx.at(5, 1);
    b.at(1, 7) = -1.0 * ( dnx.at(3, 2) + 0.5 * dnx.at(5, 2) + 0.5 * dnx.at(6, 2) );
    b.at(1, 8) = dnx.at(3, 1) + 0.5 * dnx.at(5, 1) + 0.5 * dnx.at(6, 1);

    b.at(1, 3) = -dnx.at(6, 2) * ( lxy [ 0 ].at(2) - lxy [ 2 ].at(2) ) / 8.0 + dnx.at(4, 2) * ( lxy [ 1 ].at(2) - lxy [ 0 ].at(2) ) / 8.0;
    b.at(1, 3) += -dnx.at(6, 1) * ( lxy [ 0 ].at(1) - lxy [ 2 ].at(1) ) / 8.0 + dnx.at(4, 1) * ( lxy [ 1 ].at(1) - lxy [ 0 ].at(1) ) / 8.0;
    b.at(1, 6) = -dnx.at(4, 2) * ( lxy [ 1 ].at(2) - lxy [ 0 ].at(2) ) / 8.0 + dnx.at(5, 2) * ( lxy [ 2 ].at(2) - lxy [ 1 ].at(2) ) / 8.0;
    b.at(1, 6) += -dnx.at(4, 1) * ( lxy [ 1 ].at(1) - lxy [ 0 ].at(1) ) / 8.0 + dnx.at(5, 1) * ( lxy [ 2 ].at(1) - lxy [ 1 ].at(1) ) / 8.0;
    b.at(1, 9) = -dnx.at(5, 2) * ( lxy [ 2 ].at(2) - lxy [ 1 ].at(2) ) / 8.0 + dnx.at(6, 2) * ( lxy [ 0 ].at(2) - lxy [ 2 ].at(2) ) / 8.0;
    b.at(1, 9) += -dnx.at(5, 1) * ( lxy [ 2 ].at(1) - lxy [ 1 ].at(1) ) / 8.0 + dnx.at(6, 1) * ( lxy [ 0 ].at(1) - lxy [ 2 ].at(1) ) / 8.0;
    b.times(0.5);
    // add -1.0*sum(r_w)/3.0
    b.at(1, 3) -= 1.0 / 3.0;
    b.at(1, 6) -= 1.0 / 3.0;
    b.at(1, 9) -= 1.0 / 3.0;
    // add alpha*Volume*B^T[G]B to element stiffness matrix
    double G = this->giveStructuralCrossSection()->give( Gxy, this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0) );
    double coeff = G * this->giveArea() * this->giveCrossSection()->give( CS_Thickness, this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0) ) * 1.e-6;
    answer.beTProductOf(b, b);
    answer.times(coeff);
}



void
TrPlanestressRotAllman :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v, R_w};
}



double
TrPlanestressRotAllman :: giveArea()
// returns the area occupied by the receiver
{
    if ( area > 0 ) {
        return area;         // check if previously computed
    }

    return ( area = fabs( this->interp.giveArea( FEIElementGeometryWrapper(this) ) ) );
}

void TrPlanestressRotAllman :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}

void
TrPlanestressRotAllman :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
{
    std::vector< FloatArray > lxy;
    FloatArray l, n;
    IntArray en;
    FEI2dTrQuad qi(1, 2);

    this->computeLocalNodalCoordinates(lxy); // get ready for tranformation into 3d
    qi.edgeEvalN( n, iedge, gp->giveNaturalCoordinates(), FEIVertexListGeometryWrapper(lxy) );
    qi.computeLocalEdgeMapping(en, iedge); // get edge mapping
    this->interp.edgeEvalN( l, iedge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    answer.resize(3, 6);

    answer.at(1, 1) = answer.at(2, 2) = n.at(1) + n.at(3) / 2.0;
    answer.at(1, 4) = answer.at(2, 5) = n.at(2) + n.at(3) / 2.0;
    answer.at(1, 3) = n.at(3) * ( lxy [ en.at(2) - 1 ].at(2) - lxy [ en.at(1) - 1 ].at(2) ) / 8.0;
    answer.at(1, 6) = -answer.at(1, 3);
    answer.at(2, 3) = n.at(3) * ( lxy [ en.at(2) - 1 ].at(1) - lxy [ en.at(1) - 1 ].at(1) ) / 8.0;
    answer.at(2, 6) = -answer.at(2, 3);
    answer.at(3, 3) = l.at(1);
    answer.at(3, 6) = l.at(2);
}

void
TrPlanestressRotAllman :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    answer.resize(6);
    if ( iEdge == 1 ) { // edge between nodes 1,2
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
    } else {
        OOFEM_ERROR("wrong edge number");
    }
}

/*
 * double
 * TrPlanestressRotAllman :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
 * {
 *  // edge with linear geometry -> one can use linear interpolation safely
 *  double detJ = this->interp.edgeGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
 *  return detJ *gp->giveWeight();
 * }
 *
 *
 * void
 * TrPlanestressRotAllman :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
 * {
 *  // edge with linear geometry -> one can use linear interpolation safely
 *  this->interp.edgeLocal2global( answer, iEdge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
 * }
 *
 *
 * int
 * TrPlanestressRotAllman :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
 * {
 *  // returns transformation matrix from
 *  // edge local coordinate system
 *  // to element local c.s
 *  // (same as global c.s in this case)
 *  //
 *  // i.e. f(element local) = T * f(edge local)
 *  //
 *  double dx, dy, length;
 *  Node *nodeA, *nodeB;
 *  int aNode = 0, bNode = 0;
 *
 *  answer.resize(2, 2);
 *  answer.zero();
 *
 *  if ( iEdge == 1 ) { // edge between nodes 1 2
 *      aNode = 1;
 *      bNode = 2;
 *  } else if ( iEdge == 2 ) { // edge between nodes 2 3
 *      aNode = 2;
 *      bNode = 3;
 *  } else if ( iEdge == 3 ) { // edge between nodes 2 3
 *      aNode = 3;
 *      bNode = 1;
 *  } else {
 *      OOFEM_ERROR("wrong egde number");
 *  }
 *
 *  nodeA   = this->giveNode(aNode);
 *  nodeB   = this->giveNode(bNode);
 *
 *  dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
 *  dy      = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
 *  length = sqrt(dx * dx + dy * dy);
 *
 *  answer.at(1, 1) = dx / length;
 *  answer.at(1, 2) = -dy / length;
 *  answer.at(2, 1) = dy / length;
 *  answer.at(2, 2) = dx / length;
 *
 *  return 1;
 * }
 *
 *
 *
 *
 * double TrPlanestressRotAllman :: computeVolumeAround(GaussPoint *gp)
 * // Returns the portion of the receiver which is attached to gp.
 * {
 *  double detJ, weight;
 *
 *  weight = gp->giveWeight();
 *  // safe to use linear interpolation here (geometry is linear)
 *  detJ = fabs( this->interp.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
 *
 *  return detJ *weight *this->giveCrossSection()->give(CS_Thickness), gp;
 * }
 */

IRResultType
TrPlanestressRotAllman :: initializeFrom(InputRecord *ir)
{
    IRResultType result = TrPlaneStress2d :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }
    numberOfGaussPoints = 4;
    return IRRT_OK;
}



int
TrPlanestressRotAllman :: SPRNodalRecoveryMI_giveNumberOfIP()
{
    return this->numberOfGaussPoints;
}
} // end namespace oofem
