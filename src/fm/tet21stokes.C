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

#include "tet21stokes.h"
#include "node.h"
#include "dof.h"
#include "domain.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "bcgeomtype.h"
#include "generalboundarycondition.h"
#include "load.h"
#include "bodyload.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "fm/Materials/fluiddynamicmaterial.h"
#include "fei3dtetlin.h"
#include "fei3dtetquad.h"
#include "fluidcrosssection.h"
#include "assemblercallback.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(Tet21Stokes);

// Set up interpolation coordinates
FEI3dTetLin Tet21Stokes :: interpolation_lin;
FEI3dTetQuad Tet21Stokes :: interpolation_quad;
// Set up ordering vectors (for assembling)
IntArray Tet21Stokes :: momentum_ordering = {1, 2, 3,  5, 6, 7,  9, 10, 11,  13, 14, 15,  17, 18, 19,
                                             20, 21, 22,  23, 24, 25,  26, 27, 28,  29, 30, 31,  32, 33, 34};
IntArray Tet21Stokes :: conservation_ordering = {4, 8, 12, 16};
IntArray Tet21Stokes :: surf_ordering [ 4 ] = {
    { 1,  2,  3,  9, 10, 11,  5,  6,  7, 23, 24, 25, 20, 21, 22, 17, 18, 19},
    { 1,  2,  3,  5,  6,  7, 13, 14, 15, 17, 18, 19, 29, 30, 31, 26, 27, 28},
    { 5,  6,  7,  9, 10, 11, 13, 14, 15, 20, 21, 22, 32, 33, 34, 29, 30, 31},
    { 1,  2,  3, 13, 14, 15,  9, 10, 11, 26, 27, 28, 32, 33, 34, 23, 24, 25}
};

Tet21Stokes :: Tet21Stokes(int n, Domain *aDomain) : FMElement(n, aDomain), SpatialLocalizerInterface(this)
{
    this->numberOfDofMans = 10;
    this->numberOfGaussPoints = 5;
}

void Tet21Stokes :: computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this, 1, 3);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], this->numberOfGaussPoints, this);
    }
}

int Tet21Stokes :: computeNumberOfDofs()
{
    return 34;
}

void Tet21Stokes :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( inode <= 4 ) {
        answer = {V_u, V_v, V_w, P_f};
    } else {
        answer = {V_u, V_v, V_w};
    }
}

void Tet21Stokes :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
                                             TimeStep *tStep)
{
    // Compute characteristic vector for this element. I.e the load vector(s)
    if ( mtrx == ExternalForcesVector ) {
        this->computeExternalForcesVector(answer, tStep);
    } else if ( mtrx == InternalForcesVector ) {
        this->computeInternalForcesVector(answer, tStep);
    } else {
        OOFEM_ERROR("Unknown Type of characteristic mtrx.");
    }
}

void Tet21Stokes :: giveCharacteristicMatrix(FloatMatrix &answer,
                                             CharType mtrx, TimeStep *tStep)
{
    // Compute characteristic matrix for this element. The only option is the stiffness matrix...
    if ( mtrx == TangentStiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, TangentStiffness, tStep);
    } else {
        OOFEM_ERROR("Unknown Type of characteristic mtrx.");
    }
}

void Tet21Stokes :: computeInternalForcesVector(FloatArray &answer, TimeStep *tStep)
{
    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
    FloatArray a_pressure, a_velocity, devStress, epsp, Nh, dN_V(30);
    FloatMatrix dN, B(6, 30);
    double r_vol, pressure;
    this->computeVectorOfVelocities(VM_Total, tStep, a_velocity);
    this->computeVectorOfPressures(VM_Total, tStep, a_pressure);
    FloatArray momentum, conservation;

    B.zero();
    for ( auto &gp: *integrationRulesArray [ 0 ] ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();

        double detJ = fabs( this->interpolation_quad.evaldNdx( dN, lcoords, FEIElementGeometryWrapper(this) ) );
        this->interpolation_lin.evalN( Nh, lcoords, FEIElementGeometryWrapper(this) );
        double dV = detJ * gp->giveWeight();

        for ( int j = 0, k = 0; j < dN.giveNumberOfRows(); j++, k += 3 ) {
            dN_V(k + 0) = B(0, k + 0) = B(3, k + 1) = B(4, k + 2) = dN(j, 0);
            dN_V(k + 1) = B(1, k + 1) = B(3, k + 0) = B(5, k + 2) = dN(j, 1);
            dN_V(k + 2) = B(2, k + 2) = B(4, k + 0) = B(5, k + 1) = dN(j, 2);
        }

        epsp.beProductOf(B, a_velocity);
        pressure = Nh.dotProduct(a_pressure);
        auto val = mat->computeDeviatoricStress3D(epsp, pressure, gp, tStep);
        devStress = val.first;
        r_vol = val.second;

        momentum.plusProduct(B, devStress, dV);
        momentum.add(-pressure * dV, dN_V);
        conservation.add(r_vol * dV, Nh);
    }

    answer.resize(34);
    answer.zero();
    answer.assemble(momentum, this->momentum_ordering);
    answer.assemble(conservation, this->conservation_ordering);
}

void Tet21Stokes :: computeExternalForcesVector(FloatArray &answer, TimeStep *tStep)
{
    int load_number, load_id;
    Load *load;
    BodyLoad *bLoad;
    bcGeomType ltype;
    FloatArray vec;

    int nLoads = this->boundaryLoadArray.giveSize() / 2;
    answer.clear();

    for ( int i = 1; i <= nLoads; i++ ) {  // For each Neumann boundary condition
        load_number = this->boundaryLoadArray.at(2 * i - 1);
        load_id = this->boundaryLoadArray.at(2 * i);
        load = this->domain->giveLoad(load_number);
        ltype = load->giveBCGeoType();

        if ( ltype == SurfaceLoadBGT ) {
            this->computeBoundarySurfaceLoadVector(vec, static_cast< BoundaryLoad * >(load), load_id, ExternalForcesVector, VM_Total, tStep);
            answer.add(vec);
        }
    }

    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nLoads; i++ ) {
        load  = domain->giveLoad( bodyLoadArray.at(i) );
	if ((bLoad=dynamic_cast<BodyLoad*>(load))) {
	  ltype = load->giveBCGeoType();
	  if ( ltype == BodyLoadBGT && load->giveBCValType() == ForceLoadBVT ) {
            this->computeLoadVector(vec, bLoad, ExternalForcesVector, VM_Total, tStep);
            answer.add(vec);
	  }
        }
    }
}

void Tet21Stokes :: computeLoadVector(FloatArray &answer, BodyLoad *load, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.clear();
        return;
    }

    FloatArray N, gVector, temparray(30);

    load->computeComponentArrayAt(gVector, tStep, VM_Total);
    temparray.zero();
    if ( gVector.giveSize() ) {
        for ( auto &gp: *integrationRulesArray [ 0 ] ) {
            const FloatArray &lcoords = gp->giveNaturalCoordinates();

            double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);
            double detJ = fabs( this->interpolation_quad.giveTransformationJacobian( lcoords, FEIElementGeometryWrapper(this) ) );
            double dA = detJ * gp->giveWeight();

            this->interpolation_quad.evalN( N, lcoords, FEIElementGeometryWrapper(this) );
            for ( int j = 0; j < N.giveSize(); j++ ) {
                temparray(3 * j + 0) += N(j) * rho * gVector(0) * dA;
                temparray(3 * j + 1) += N(j) * rho * gVector(1) * dA;
                temparray(3 * j + 2) += N(j) * rho * gVector(2) * dA;
            }
        }
    }

    answer.resize(34);
    answer.zero();
    answer.assemble(temparray, this->momentum_ordering);
}

  void Tet21Stokes :: computeBoundarySurfaceLoadVector(FloatArray &answer, BoundaryLoad *load, int iSurf, CharType type, ValueModeType mode, TimeStep *tStep, bool global)
{
    if ( type != ExternalForcesVector ) {
        answer.clear();
        return;
    }

    answer.resize(34);
    answer.zero();

    if ( load->giveType() == TransmissionBC ) { // Neumann boundary conditions (traction)

        int numberOfSurfaceIPs = ( int ) ceil( ( load->giveApproxOrder() + 1. ) / 2. ) * 2; ///@todo Check this.

        GaussIntegrationRule iRule(1, this, 1, 1);
        FloatArray N, t, f(18);

        f.zero();
        iRule.SetUpPointsOnTriangle(numberOfSurfaceIPs, _Unknown);

        for ( auto &gp: iRule ) {
            const FloatArray &lcoords = gp->giveNaturalCoordinates();

            this->interpolation_quad.surfaceEvalN( N, iSurf, lcoords, FEIElementGeometryWrapper(this) );
            double dA = gp->giveWeight() * this->interpolation_quad.surfaceGiveTransformationJacobian( iSurf, lcoords, FEIElementGeometryWrapper(this) );

            if ( load->giveFormulationType() == Load :: FT_Entity ) { // load in xi-eta system
                load->computeValueAt(t, tStep, lcoords, VM_Total);
            } else { // Edge load in x-y system
                FloatArray gcoords;
                this->interpolation_quad.surfaceLocal2global( gcoords, iSurf, lcoords, FEIElementGeometryWrapper(this) );
                load->computeValueAt(t, tStep, gcoords, VM_Total);
            }

            // Reshape the vector
            for ( int j = 0; j < N.giveSize(); j++ ) {
                f(3 * j + 0) += N(j) * t(0) * dA;
                f(3 * j + 1) += N(j) * t(1) * dA;
                f(3 * j + 2) += N(j) * t(2) * dA;
            }
        }

        answer.assemble(f, this->surf_ordering [ iSurf - 1 ]);
    } else {
        OOFEM_ERROR("Strange boundary condition type");
    }
}

void Tet21Stokes :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep)
{
    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
    FloatMatrixF<30,30> K;
    FloatMatrixF<30,4> G, Dp;
    FloatMatrixF<4,30> DvT;
    FloatMatrixF<4,4> C;

    for ( auto &gp: *this->integrationRulesArray [ 0 ] ) {
        const auto &lcoords = gp->giveNaturalCoordinates();

        auto Nlin = this->interpolation_lin.evalN( lcoords );
        auto detj_dn = this->interpolation_quad.evaldNdx( lcoords, FEIElementGeometryWrapper(this) );
        auto dN = detj_dn.second;
        auto detJ = detj_dn.first;
        auto dA = std::abs(detJ) * gp->giveWeight();

        auto dN_V = flatten(dN);
        auto B = Bmatrix_3d(dN);

        // Computing the internal forces should have been done first.
        // dsigma_dev/deps_dev  dsigma_dev/dp  deps_vol/deps_dev  deps_vol/dp
        //auto [Ed, Ep, Cd, Cp] = mat->computeTangents2D(mode, gp, tStep);
        auto tangents = mat->computeTangents3D(mode, gp, tStep);

        K.plusProductSymmUpper(B, dot(tangents.dsdd, B), dA);
        G.plusDyadUnsym(dN_V, Nlin, -dA);
        Dp.plusDyadUnsym(Tdot(B, tangents.dsdp), Nlin, dA);
        DvT.plusDyadUnsym(Nlin, Tdot(B, tangents.dedd), dA);
        C.plusDyadSymmUpper(Nlin, tangents.dedp * dA);
    }

    K.symmetrized();
    C.symmetrized();
    auto GTDvT = transpose(G) + DvT;
    auto GDp = G + Dp;

    answer.resize(34, 34);
    answer.zero();
    answer.assemble(K, this->momentum_ordering);
    answer.assemble(GDp, this->momentum_ordering, this->conservation_ordering);
    answer.assemble(GTDvT, this->conservation_ordering, this->momentum_ordering);
    answer.assemble(C, this->conservation_ordering);
}

FEInterpolation *Tet21Stokes :: giveInterpolation() const
{
    return & interpolation_quad;
}

FEInterpolation *Tet21Stokes :: giveInterpolation(DofIDItem id) const
{
    if ( id == P_f ) {
        return & interpolation_lin;
    } else {
        return & interpolation_quad;
    }
}

void Tet21Stokes :: updateYourself(TimeStep *tStep)
{
    Element :: updateYourself(tStep);
}

// Some extension Interfaces to follow:

Interface *Tet21Stokes :: giveInterface(InterfaceType it)
{
    switch ( it ) {
    case NodalAveragingRecoveryModelInterfaceType:
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);

    case SpatialLocalizerInterfaceType:
        return static_cast< SpatialLocalizerInterface * >(this);

    case EIPrimaryUnknownMapperInterfaceType:
        return static_cast< EIPrimaryUnknownMapperInterface * >(this);

    default:
        return FMElement :: giveInterface(it);
    }
}

void Tet21Stokes :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType mode,
                                                                          TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
{
    FloatArray n, n_lin;
    this->interpolation_quad.evalN( n, lcoords, FEIElementGeometryWrapper(this) );
    this->interpolation_lin.evalN( n_lin, lcoords, FEIElementGeometryWrapper(this) );
    answer.resize(4);
    answer.zero();
    for ( int i = 1; i <= n.giveSize(); i++ ) {
        answer(0) += n.at(i) * this->giveNode(i)->giveDofWithID(V_u)->giveUnknown(mode, tStep);
        answer(1) += n.at(i) * this->giveNode(i)->giveDofWithID(V_v)->giveUnknown(mode, tStep);
        answer(2) += n.at(i) * this->giveNode(i)->giveDofWithID(V_w)->giveUnknown(mode, tStep);
    }

    for ( int i = 1; i <= n_lin.giveSize(); i++ ) {
        answer(3) += n_lin.at(i) * this->giveNode(i)->giveDofWithID(P_f)->giveUnknown(mode, tStep);
    }
}

void Tet21Stokes :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_Pressure ) {
        answer.resize(1);
        if ( node <= 4 ) {
            answer.at(1) = this->giveNode(node)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep);
        } else {
            const auto &eNodes = this->interpolation_quad.computeLocalEdgeMapping(node - 4);
            answer.at(1) = 0.5 * (
                this->giveNode( eNodes.at(1) )->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) +
                this->giveNode( eNodes.at(2) )->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) );
        }
    } else {
        answer.clear();
    }
}

} // end namespace oofem
