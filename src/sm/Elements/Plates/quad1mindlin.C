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

#include "sm/Elements/Plates/quad1mindlin.h"
#include "sm/Materials/structuralms.h"
#include "sm/CrossSections/structuralcrosssection.h"
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

FEI2dQuadLin Quad1Mindlin::interp_lin(1, 2);

Quad1Mindlin::Quad1Mindlin(int n, Domain *aDomain) :
    StructuralElement(n, aDomain), ZZNodalRecoveryModelInterface(this),
    SPRNodalRecoveryModelInterface()
{
    numberOfGaussPoints = 4;
    numberOfDofMans = 4;
}

FEInterpolation *
Quad1Mindlin::giveInterpolation() const
{
    return & interp_lin;
}

FEInterpolation *
Quad1Mindlin::giveInterpolation(DofIDItem id) const
{
    return & interp_lin;
}

void
Quad1Mindlin::computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ] = std::make_unique< GaussIntegrationRule >(1, this, 1, 5);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


void
Quad1Mindlin::computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
{
    // Only gravity load
    FloatArray gravity;

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        OOFEM_ERROR("unknown load type");
    }

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(gravity, tStep, mode);

    if ( gravity.giveSize() ) {
        FloatArrayF< 4 >force;
        ///@todo Other/higher integration for lumped mass matrices perhaps?
        for ( auto &gp: * integrationRulesArray [ 0 ] ) {
            auto n = this->interp_lin.evalN( gp->giveNaturalCoordinates() );
            double dV = this->computeVolumeAround(gp) * this->giveCrossSection()->give(CS_Thickness, gp);
            double load = this->giveStructuralCrossSection()->give('d', gp) * gravity.at(3) * dV;
            force += load * n;
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
Quad1Mindlin::computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [5x9] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    auto dn = this->interp_lin.evaldNdx( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ).second;
    auto n = this->interp_lin.evalN( gp->giveNaturalCoordinates() );

    // enforce one-point reduced integration if requested
    FloatArrayF< 4 >ns;
    FloatMatrixF< 2, 4 >dns;
    if ( this->reducedIntegrationFlag ) {
        FloatArrayF< 2 >lc; // set to element center coordinates
        dns = this->interp_lin.evaldNdx( lc, FEIElementGeometryWrapper(this) ).second;
        ns = this->interp_lin.evalN(lc);
    } else {
        dns = dn;
        ns = n;
    }

    answer.resize(5, 12);
    answer.zero();

    ///@todo Check sign here
    for ( int i = 0; i < 4; ++i ) {
        answer(0, 2 + i * 3) =  dn(0, i);// kappa_x = d(fi_y)/dx
        answer(1, 1 + i * 3) = -dn(1, i); // kappa_y = -d(fi_x)/dy
        answer(2, 2 + i * 3) =  dn(1, i);// kappa_xy=d(fi_y)/dy-d(fi_x)/dx
        answer(2, 1 + i * 3) = -dn(0, i);

        answer(3, 0 + i * 3) = dns(0, i); // gamma_xz = fi_y+dw/dx
        answer(3, 2 + i * 3) = ns [ i ];
        answer(4, 0 + i * 3) = dns(1, i);// gamma_yz = -fi_x+dw/dy
        answer(4, 1 + i * 3) = -ns [ i ];
    }
}


void
Quad1Mindlin::computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    answer = this->giveStructuralCrossSection()->giveGeneralizedStress_Plate(strain, gp, tStep);
}


void
Quad1Mindlin::computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    answer = this->giveStructuralCrossSection()->give2dPlateStiffMtrx(rMode, gp, tStep);
}


void
Quad1Mindlin::initializeFrom(InputRecord &ir)
{
    this->numberOfGaussPoints = 4;
    this->reducedIntegrationFlag = ir.hasField(_IFT_Quad1Mindlin_ReducedIntegration);
    StructuralElement::initializeFrom(ir);
}


void
Quad1Mindlin::giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = { D_w, R_u, R_v };
}


void
Quad1Mindlin::computeMidPlaneNormal(FloatArray &answer, const GaussPoint *gp)
{
    FloatArrayF< 3 >u = this->giveNode(2)->giveCoordinates() - this->giveNode(1)->giveCoordinates();
    FloatArrayF< 3 >v = this->giveNode(3)->giveCoordinates() - this->giveNode(1)->giveCoordinates();

    auto n = cross(u, v);
    answer = n / norm(n);
}


double
Quad1Mindlin::giveCharacteristicLength(const FloatArray &normalToCrackPlane)
{
    return this->giveCharacteristicLengthForPlaneElements(normalToCrackPlane);
}


double
Quad1Mindlin::computeVolumeAround(GaussPoint *gp)
{
    double weight = gp->giveWeight();
    double detJ = fabs(this->interp_lin.giveTransformationJacobian(gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    return detJ * weight;
}


void
Quad1Mindlin::computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver.
{
    ///@todo Other/higher integration for lumped mass matrices perhaps?
    double mass = 0.;
    for ( GaussPoint *gp: * integrationRulesArray [ 0 ] ) {
        double dV = this->computeVolumeAround(gp);
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
Quad1Mindlin::giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    answer.resize(6);
    auto ms = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );
    if ( type == IST_ShellForceTensor || type == IST_ShellStrainTensor ) {
        const auto &s = type == IST_ShellForceTensor ? ms->giveStressVector() : ms->giveStrainVector();
        answer.at(1) = 0.0; // nx
        answer.at(2) = 0.0; // ny
        answer.at(3) = 0.0; // nz
        answer.at(4) = s.at(5); // vyz
        answer.at(5) = s.at(4); // vxz
        answer.at(6) = 0.0; // vxy
        return 1;
    } else if ( type == IST_ShellMomentTensor || type == IST_CurvatureTensor ) {
        const auto &s = type == IST_ShellMomentTensor ? ms->giveStressVector() : ms->giveStrainVector();
        answer.at(1) = s.at(1); // mx
        answer.at(2) = s.at(2); // my
        answer.at(3) = 0.0;     // mz
        answer.at(4) = 0.0;     // myz
        answer.at(5) = 0.0;     // mxz
        answer.at(6) = s.at(3); // mxy
        return 1;
    } else {
        return StructuralElement::giveIPValue(answer, gp, type, tStep);
    }
}

void
Quad1Mindlin::giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    if ( iEdge == 1 ) { // edge between nodes 1 2
        answer = { 1, 2, 3, 4, 5, 6 };
    } else if ( iEdge == 2 ) { // edge between nodes 2 3
        answer = { 4, 5, 6, 7, 8, 9 };
    } else if ( iEdge == 3 ) { // edge between nodes 3 4
        answer = { 7, 8, 9, 10, 11, 12 };
    } else if ( iEdge == 4 ) { // edge between nodes 4 1
        answer = { 10, 11, 12, 1, 2, 3 };
    } else {
        OOFEM_ERROR("wrong edge number");
    }
}


double
Quad1Mindlin::computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double detJ = this->interp_lin.edgeGiveTransformationJacobian(iEdge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    return detJ * gp->giveWeight();
}


int
Quad1Mindlin::computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    const auto &edgeNodes = this->interp_lin.computeLocalEdgeMapping(iEdge);

    auto nodeA = this->giveNode(edgeNodes.at(1) );
    auto nodeB = this->giveNode(edgeNodes.at(2) );

    double dx = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
    double dy = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
    double length = sqrt(dx * dx + dy * dy);

    answer.resize(3, 3);
    answer.zero();
    answer.at(1, 1) = 1.0;
    answer.at(2, 2) = dx / length;
    answer.at(2, 3) = -dy / length;
    answer.at(3, 2) = -answer.at(2, 3);
    answer.at(3, 3) = answer.at(2, 2);

    return 1;
}
Interface *
Quad1Mindlin::giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >( this );
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >( this );
    }

    return nullptr;
}



void
Quad1Mindlin::SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(4);
    for ( int i = 1; i < 5; i++ ) {
        pap.at(i) = this->giveNode(i)->giveNumber();
    }
}

void
Quad1Mindlin::SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
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
