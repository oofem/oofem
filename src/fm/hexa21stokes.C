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

#include "hexa21stokes.h"
#include "node.h"
#include "dof.h"
#include "domain.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "bcgeomtype.h"
#include "generalboundarycondition.h"
#include "load.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "fluiddynamicmaterial.h"
#include "fei3dhexalin.h"
#include "fei3dhexatriquad.h"
#include "fluidcrosssection.h"
#include "assemblercallback.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(Hexa21Stokes);

// Set up interpolation coordinates
FEI3dHexaLin Hexa21Stokes :: interpolation_lin;
FEI3dHexaTriQuad Hexa21Stokes :: interpolation_quad;
// Set up ordering vectors (for assembling)
IntArray Hexa21Stokes :: momentum_ordering = {
     1,  2,  3,  5,  6,  7,  9, 10, 11, 13, 14, 15, 17, 18, 19, 21, 22, 23, 25, 26, 27, 29, 30, 31, 33, 34, 35,
    37, 38, 39, 41, 42, 43, 45, 46, 47, 49, 50, 51, 53, 54, 55, 57, 58, 59, 61, 62, 63, 65, 66, 67, 69, 70, 71,
    73, 74, 75, 77, 78, 79, 81, 82, 83, 85, 86, 87, 89, 90, 91, 93, 94, 95, 97, 98, 99, 101, 102, 103, 105, 106, 107};
IntArray Hexa21Stokes :: conservation_ordering = {4, 8, 12, 16, 20, 24, 28, 32};
IntArray Hexa21Stokes :: surf_ordering [ 6 ] = {
    { 5,  6,  7,  1,  2,  3, 13, 14, 15,  9, 10, 11, 33, 34, 35, 42, 43, 44, 39, 40, 41, 36, 37, 38, 69, 70, 71},
    {17, 18, 19, 21, 22, 23, 25, 26, 27, 29, 30, 31, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 72, 73, 74},
    { 1,  2,  3, 17, 18, 19, 21, 22, 23,  5,  6,  7, 57, 58, 59, 45, 46, 47, 60, 61, 62, 33, 34, 35, 75, 76, 77},
    { 5,  6,  7,  9, 10, 11, 25, 26, 27, 21, 22, 23, 36, 37, 38, 63, 64, 65, 48, 49, 50, 60, 61, 62, 78, 79, 80},
    { 9, 10, 11, 13, 14, 15, 29, 30, 31, 25, 26, 27, 39, 40, 41, 66, 67, 68, 51, 52, 53, 63, 64, 65, 81, 82, 83},
    {13, 14, 15,  1,  2,  3, 17, 18, 19, 29, 30, 31, 42, 43, 44, 57, 58, 59, 54, 55, 56, 66, 67, 68, 84, 85, 86}
};

Hexa21Stokes :: Hexa21Stokes(int n, Domain *aDomain) : FMElement(n, aDomain), SpatialLocalizerInterface(this)
{
    this->numberOfDofMans = 27;
    this->numberOfGaussPoints = 27;
}

Hexa21Stokes :: ~Hexa21Stokes()
{ }

void Hexa21Stokes :: computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}

int Hexa21Stokes :: computeNumberOfDofs()
{
    return 89;
}

void Hexa21Stokes :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( inode <= 8 ) {
        answer = {V_u, V_v, V_w, P_f};
    } else {
        answer = {V_u, V_v, V_w};
    }
}

void Hexa21Stokes :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
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

void Hexa21Stokes :: giveCharacteristicMatrix(FloatMatrix &answer,
                                              CharType mtrx, TimeStep *tStep)
{
    // Compute characteristic matrix for this element. The only option is the stiffness matrix...
    if ( mtrx == TangentStiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, TangentStiffness, tStep);
    } else {
        OOFEM_ERROR("Unknown Type of characteristic mtrx.");
    }
}

void Hexa21Stokes :: computeInternalForcesVector(FloatArray &answer, TimeStep *tStep)
{
    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
    FloatArray a_pressure, a_velocity, devStress, epsp, Nh, dN_V(81);
    FloatMatrix dN, B(6, 81);
    double r_vol, pressure;
    this->computeVectorOfVelocities(VM_Total, tStep, a_velocity);
    this->computeVectorOfPressures(VM_Total, tStep, a_pressure);
    FloatArray momentum, conservation;

    B.zero();
    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();

        double detJ = fabs( this->interpolation_quad.evaldNdx( dN, lcoords, FEIElementGeometryWrapper(this) ) );
        this->interpolation_lin.evalN( Nh, lcoords, FEIElementGeometryWrapper(this) );
        double dV = detJ * gp->giveWeight();

        for ( int j = 0, k = 0; j < dN.giveNumberOfRows(); j++, k += 3 ) {
            dN_V(k + 0) = B(0, k + 0) = B(5, k + 1) = B(4, k + 2) = dN(j, 0);
            dN_V(k + 1) = B(1, k + 1) = B(5, k + 0) = B(3, k + 2) = dN(j, 1);
            dN_V(k + 2) = B(2, k + 2) = B(4, k + 0) = B(3, k + 1) = dN(j, 2);
        }

        epsp.beProductOf(B, a_velocity);
        pressure = Nh.dotProduct(a_pressure);
        mat->computeDeviatoricStressVector(devStress, r_vol, gp, epsp, pressure, tStep);

        momentum.plusProduct(B, devStress, dV);
        momentum.add(-pressure * dV, dN_V);
        conservation.add(r_vol * dV, Nh);
    }

    answer.resize(89);
    answer.zero();
    answer.assemble(momentum, this->momentum_ordering);
    answer.assemble(conservation, this->conservation_ordering);
}

void Hexa21Stokes :: computeExternalForcesVector(FloatArray &answer, TimeStep *tStep)
{
    FloatArray vec;

    int nLoads = this->boundaryLoadArray.giveSize() / 2;
    answer.clear();

    for ( int i = 1; i <= nLoads; i++ ) {  // For each Neumann boundary condition
        int load_number = this->boundaryLoadArray.at(2 * i - 1);
        int load_id = this->boundaryLoadArray.at(2 * i);
        Load *load = this->domain->giveLoad(load_number);

        if ( load->giveBCGeoType() == SurfaceLoadBGT ) {
            this->computeBoundaryLoadVector(vec, static_cast< BoundaryLoad * >(load), load_id, ExternalForcesVector, VM_Total, tStep);
            answer.add(vec);
        }
    }

    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nLoads; i++ ) {
        Load *load = domain->giveLoad( bodyLoadArray.at(i) );
        if ( load->giveBCGeoType() == BodyLoadBGT && load->giveBCValType() == ForceLoadBVT ) {
            this->computeLoadVector(vec, load, ExternalForcesVector, VM_Total, tStep);
            answer.add(vec);
        }
    }
}

void Hexa21Stokes :: computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.clear();
        return;
    }

    FloatArray N, gVector, temparray(81);

    load->computeComponentArrayAt(gVector, tStep, VM_Total);
    temparray.zero();
    if ( gVector.giveSize() ) {
        for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
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

    answer.resize(89);
    answer.zero();
    answer.assemble(temparray, this->momentum_ordering);
}

void Hexa21Stokes :: computeBoundaryLoadVector(FloatArray &answer, BoundaryLoad *load, int iSurf, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.clear();
        return;
    }

    if ( load->giveType() == TransmissionBC ) { // Neumann boundary conditions (traction)

        int numberOfSurfaceIPs = ( int ) ceil( ( load->giveApproxOrder() + 1. ) / 2. ) * 2; ///@todo Check this.

        GaussIntegrationRule iRule(1, this, 1, 1);
        FloatArray N, t, f(27);

        f.zero();
        iRule.SetUpPointsOnTriangle(numberOfSurfaceIPs, _Unknown);

        for ( GaussPoint *gp: iRule ) {
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

        answer.resize(89);
        answer.zero();
        answer.assemble(f, this->surf_ordering [ iSurf - 1 ]);
    } else {
        OOFEM_ERROR("Strange boundary condition type");
    }
}

void Hexa21Stokes :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep)
{
    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
    FloatMatrix B(6, 81), EdB, K, G, Dp, DvT, C, Ed, dN;
    FloatArray dN_V(81), Nlin, Ep, Cd, tmpA, tmpB;
    double Cp;

    B.zero();

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        // Compute Gauss point and determinant at current element
        const FloatArray &lcoords = gp->giveNaturalCoordinates();

        double detJ = fabs( this->interpolation_quad.evaldNdx( dN, lcoords, FEIElementGeometryWrapper(this) ) );
        double dV = detJ * gp->giveWeight();
        this->interpolation_lin.evalN( Nlin, lcoords, FEIElementGeometryWrapper(this) );

        for ( int j = 0, k = 0; j < dN.giveNumberOfRows(); j++, k += 3 ) {
            dN_V(k + 0) = B(0, k + 0) = B(5, k + 1) = B(4, k + 2) = dN(j, 0);
            dN_V(k + 1) = B(1, k + 1) = B(5, k + 0) = B(3, k + 2) = dN(j, 1);
            dN_V(k + 2) = B(2, k + 2) = B(4, k + 0) = B(3, k + 1) = dN(j, 2);
        }

        // Computing the internal forces should have been done first.
        // dsigma_dev/deps_dev  dsigma_dev/dp  deps_vol/deps_dev  deps_vol/dp
        mat->giveStiffnessMatrices(Ed, Ep, Cd, Cp, mode, gp, tStep);

        EdB.beProductOf(Ed, B);
        K.plusProductSymmUpper(B, EdB, dV);
        G.plusDyadUnsym(dN_V, Nlin, -dV);
        C.plusDyadSymmUpper(Nlin, Cp * dV);

        tmpA.beTProductOf(B, Ep);
        Dp.plusDyadUnsym(tmpA, Nlin, dV);

        tmpB.beTProductOf(B, Cd);
        DvT.plusDyadUnsym(Nlin, tmpB, dV);
    }

    K.symmetrized();
    C.symmetrized();
    FloatMatrix GTDvT, GDp;
    GTDvT.beTranspositionOf(G);
    GTDvT.add(DvT);
    GDp = G;
    GDp.add(Dp);

    answer.resize(89, 89);
    answer.zero();
    answer.assemble(K, this->momentum_ordering);
    answer.assemble(GDp, this->momentum_ordering, this->conservation_ordering);
    answer.assemble(GTDvT, this->conservation_ordering, this->momentum_ordering);
    answer.assemble(C, this->conservation_ordering);
}

FEInterpolation *Hexa21Stokes :: giveInterpolation() const
{
    return & interpolation_quad;
}

FEInterpolation *Hexa21Stokes :: giveInterpolation(DofIDItem id) const
{
    if ( id == P_f ) {
        return & interpolation_lin;
    } else {
        return & interpolation_quad;
    }
}

void Hexa21Stokes :: updateYourself(TimeStep *tStep)
{
    Element :: updateYourself(tStep);
}

// Some extension Interfaces to follow:

Interface *Hexa21Stokes :: giveInterface(InterfaceType it)
{
    switch ( it ) {
    case NodalAveragingRecoveryModelInterfaceType:
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);

    case SpatialLocalizerInterfaceType:
        return static_cast< SpatialLocalizerInterface * >(this);

    default:
        return FMElement :: giveInterface(it);
    }
}


void Hexa21Stokes :: computeField(ValueModeType mode, TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
{
    FloatArray n, n_lin, pressures, velocities;
    this->interpolation_quad.evalN( n, lcoords, FEIElementGeometryWrapper(this) );
    this->interpolation_lin.evalN( n_lin, lcoords, FEIElementGeometryWrapper(this) );
    this->computeVectorOf({P_f}, mode, tStep, pressures);
    this->computeVectorOf({V_u, V_v, V_w}, mode, tStep, velocities);

    answer.resize(4);
    answer.zero();
    for ( int i = 1; i <= n.giveSize(); i++ ) {
        answer(0) += n.at(i) * velocities.at(i*3-2);
        answer(1) += n.at(i) * velocities.at(i*3-1);
        answer(2) += n.at(i) * velocities.at(i*3);
    }

    for ( int i = 1; i <= n_lin.giveSize(); i++ ) {
        answer(3) += n_lin.at(i) * pressures.at(i);
    }
}


void Hexa21Stokes :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    answer.resize(1);
    if ( type == IST_Pressure ) {
        answer.resize(1);
        if ( node <= 8 ) { // Corner nodes
            answer.at(1) = this->giveNode(node)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep);
        } else if ( node <= 20 ) { // Edge nodes
            // Edges are numbered consistently with edge nodes, so node number - 8 = edge number
            IntArray eNodes;
            this->interpolation_quad.computeLocalEdgeMapping(eNodes, node - 8);
            answer.at(1) = 0.5 * (
                this->giveNode( eNodes.at(1) )->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) +
                this->giveNode( eNodes.at(2) )->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) );
        } else if ( node <= 26 ) { // Face nodes
            // Faces are numbered consistently with edge nodes, so node number - 12 = face number
            IntArray fNodes;
            this->interpolation_quad.computeLocalSurfaceMapping(fNodes, node - 20);
            answer.at(1) = 0.25 * (
                this->giveNode( fNodes.at(1) )->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) +
                this->giveNode( fNodes.at(2) )->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) +
                this->giveNode( fNodes.at(3) )->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) +
                this->giveNode( fNodes.at(4) )->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) );
        } else { // Middle node
            answer.at(1) = 0.125 * (
                this->giveNode(1)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) +
                this->giveNode(2)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) +
                this->giveNode(3)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) +
                this->giveNode(4)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) +
                this->giveNode(5)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) +
                this->giveNode(6)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) +
                this->giveNode(7)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) +
                this->giveNode(8)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) );
        }
    } else {
        answer.clear();
    }
}
} // end namespace oofem
