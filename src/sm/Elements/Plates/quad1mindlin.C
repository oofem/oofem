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

#include "../sm/Elements/Plates/quad1mindlin.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "load.h"
#include "mathfem.h"
#include "fei2dquadlin.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(Quad1Mindlin);

FEI2dQuadLin Quad1Mindlin :: interp_lin(1, 2);

Quad1Mindlin :: Quad1Mindlin(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(this),
    SPRNodalRecoveryModelInterface()
{
    numberOfGaussPoints = 4;
    numberOfDofMans = 4;
    this->reducedIntegrationFlag = false;
}

FEInterpolation *
Quad1Mindlin :: giveInterpolation() const
{
    return & interp_lin;
}

FEInterpolation *
Quad1Mindlin :: giveInterpolation(DofIDItem id) const
{
    return & interp_lin;
}

void
Quad1Mindlin :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 5) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


void
Quad1Mindlin :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
{
    // Only gravity load
    double dV, load;
    FloatArray force, gravity, n;

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        OOFEM_ERROR("unknown load type");
    }

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(gravity, tStep, mode);

    force.clear();
    if ( gravity.giveSize() ) {
        ///@todo Other/higher integration for lumped mass matrices perhaps?
        for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {

            this->interp_lin.evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
            dV = this->computeVolumeAround(gp) * this->giveCrossSection()->give(CS_Thickness, gp);
            load = this->giveStructuralCrossSection()->give('d', gp) * gravity.at(3) * dV;

            force.add(load, n);
        }

        answer.resize(12);
        answer.zero();

        answer.at(1)  = force.at(1);
        answer.at(4)  = force.at(2);
        answer.at(7)  = force.at(3);
        answer.at(10) = force.at(4);
    } else {
        answer.clear();
    }
}


void
Quad1Mindlin :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [5x9] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    FloatArray n, ns;
    FloatMatrix dn, dns;

    this->interp_lin.evaldNdx( dn, gp->giveNaturalCoordinates(),  FEIElementGeometryWrapper(this) );
    this->interp_lin.evalN( n, gp->giveNaturalCoordinates(),  FEIElementGeometryWrapper(this) );

    answer.resize(5, 12);
    answer.zero();

    // enforce one-point reduced integration if requested
    if ( this->reducedIntegrationFlag ) {
        FloatArray lc(2);
        lc.zero(); // set to element center coordinates

        this->interp_lin.evaldNdx( dns, lc, FEIElementGeometryWrapper(this));
        this->interp_lin.evalN( ns, lc,  FEIElementGeometryWrapper(this));
     } else {
        dns = dn;
        ns = n;
    }

    ///@todo Check sign here
    for ( int i = 0; i < 4; ++i ) {
      answer(0, 2 + i * 3) =  dn(i, 0); // kappa_x = d(fi_y)/dx
      answer(1, 1 + i * 3) = -dn(i, 1); // kappa_y = -d(fi_x)/dy
      answer(2, 2 + i * 3) =  dn(i, 1); // kappa_xy=d(fi_y)/dy-d(fi_x)/dx
      answer(2, 1 + i * 3) = -dn(i, 0);

      answer(3, 0 + i * 3) = dns(i, 0); // gamma_xz = fi_y+dw/dx
      answer(3, 2 + i * 3) = ns(i);
      answer(4, 0 + i * 3) = dns(i, 1);// gamma_yz = -fi_x+dw/dy
      answer(4, 1 + i * 3) = -ns(i);
    }
}


void
Quad1Mindlin :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveGeneralizedStress_Plate(answer, gp, strain, tStep);
}


void
Quad1Mindlin :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->give2dPlateStiffMtrx(answer, rMode, gp, tStep);
}


IRResultType
Quad1Mindlin :: initializeFrom(InputRecord *ir)
{
    this->numberOfGaussPoints = 4;
    this->reducedIntegrationFlag = ir->hasField(_IFT_Quad1Mindlin_ReducedIntegration);
    return NLStructuralElement :: initializeFrom(ir);
}


void
Quad1Mindlin :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_w, R_u, R_v};
}


void
Quad1Mindlin :: computeMidPlaneNormal(FloatArray &answer, const GaussPoint *gp)
{
    FloatArray u, v;
    u.beDifferenceOf( * this->giveNode(2)->giveCoordinates(), * this->giveNode(1)->giveCoordinates() );
    v.beDifferenceOf( * this->giveNode(3)->giveCoordinates(), * this->giveNode(1)->giveCoordinates() );

    answer.beVectorProductOf(u, v);
    answer.normalize();
}


double
Quad1Mindlin :: giveCharacteristicLength(const FloatArray &normalToCrackPlane)
{
    return this->giveCharacteristicLengthForPlaneElements(normalToCrackPlane);
}


double
Quad1Mindlin :: computeVolumeAround(GaussPoint *gp)
{
    double detJ, weight;

    weight = gp->giveWeight();
    detJ = fabs( this->interp_lin.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    return detJ * weight;
}


void
Quad1Mindlin :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver.
{
    double dV, mass = 0.;

    ///@todo Other/higher integration for lumped mass matrices perhaps?
    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {

        dV = this->computeVolumeAround(gp);
        mass += dV * this->giveStructuralCrossSection()->give('d', gp);
    }

    answer.resize(12, 12);
    answer.zero();
    answer.at(1, 1) = mass * 0.25;
    answer.at(4, 4) = mass * 0.25;
    answer.at(7, 7) = mass * 0.25;
    answer.at(10, 10) = mass * 0.25;
}


int
Quad1Mindlin :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    FloatArray help;
    answer.resize(6);
    if ( type == IST_ShellForceTensor || type == IST_ShellStrainTensor ) {
        if ( type == IST_ShellForceTensor ) {
            help = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStressVector();
        } else {
            help = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStrainVector();
        }
        answer.at(1) = 0.0; // nx
        answer.at(2) = 0.0; // ny
        answer.at(3) = 0.0; // nz
        answer.at(4) = help.at(5); // vyz
        answer.at(5) = help.at(4); // vxz
        answer.at(6) = 0.0; // vxy
        return 1;
    } else if ( type == IST_ShellMomentumTensor || type == IST_ShellCurvatureTensor ) {
        if ( type == IST_ShellMomentumTensor ) {
            help = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStressVector();
        } else {
            help = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStrainVector();
        }
        answer.at(1) = help.at(1); // mx
        answer.at(2) = help.at(2); // my
        answer.at(3) = 0.0;        // mz
        answer.at(4) = 0.0;        // myz
        answer.at(5) = 0.0;        // mxz
        answer.at(6) = help.at(3); // mxy
        return 1;
    } else {
        return NLStructuralElement :: giveIPValue(answer, gp, type, tStep);
    }
}

void
Quad1Mindlin :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
{
    IntArray edgeNodes;
    FloatArray n;

    this->interp_lin.edgeEvalN( n, iedge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    this->interp_lin.computeLocalEdgeMapping(edgeNodes, iedge);

    answer.beNMatrixOf(n, 3);
}


void
Quad1Mindlin :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    if ( iEdge == 1 ) { // edge between nodes 1 2
        answer = {1, 2, 3, 4, 5, 6};
    } else if ( iEdge == 2 ) { // edge between nodes 2 3
        answer = {4, 5, 6, 7, 8, 9};
    } else if ( iEdge == 3 ) { // edge between nodes 3 4
        answer = {7, 8, 9, 10, 11, 12};
    } else if ( iEdge == 4 ) { // edge between nodes 4 1
        answer = {10, 11, 12, 1, 2, 3};
    } else {
        OOFEM_ERROR("wrong edge number");
    }
}


double
Quad1Mindlin :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double detJ = this->interp_lin.edgeGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    return detJ *gp->giveWeight();
}


void
Quad1Mindlin :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    this->interp_lin.edgeLocal2global( answer, iEdge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
}


int
Quad1Mindlin :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    double dx, dy, length;
    IntArray edgeNodes;
    Node *nodeA, *nodeB;

    answer.resize(3, 3);
    answer.zero();

    this->interp_lin.computeLocalEdgeMapping(edgeNodes, iEdge);

    nodeA = this->giveNode( edgeNodes.at(1) );
    nodeB = this->giveNode( edgeNodes.at(2) );

    dx = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
    dy = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
    length = sqrt(dx * dx + dy * dy);

    answer.at(1, 1) = 1.0;
    answer.at(2, 2) = dx / length;
    answer.at(2, 3) = -dy / length;
    answer.at(3, 2) = -answer.at(2, 3);
    answer.at(3, 3) = answer.at(2, 2);

    return 1;
}
Interface *
Quad1Mindlin :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    }

    return NULL;
}



void
Quad1Mindlin :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(4);
    for ( int i = 1; i < 5; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
Quad1Mindlin :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    int found = 0;
    answer.resize(1);

    for ( int i = 1; i < 5; i++ ) {
        if ( pap == this->giveNode(i)->giveNumber() ) {
            found = 1;
        }
    }

    if ( found ) {
        answer.at(1) = pap;
    } else {
        OOFEM_ERROR("node unknown");
    }
}
} // end namespace oofem
