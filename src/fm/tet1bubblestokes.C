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

#include "tet1bubblestokes.h"
#include "node.h"
#include "domain.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "bcgeomtype.h"
#include "load.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "fluiddynamicmaterial.h"
#include "fei3dtetlin.h"
#include "masterdof.h"
#include "fluidcrosssection.h"
#include "assemblercallback.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(Tet1BubbleStokes);

FEI3dTetLin Tet1BubbleStokes :: interp;
// Set up ordering vectors (for assembling)
IntArray Tet1BubbleStokes :: momentum_ordering = {1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15, 17, 18, 19};
IntArray Tet1BubbleStokes :: conservation_ordering = {4, 8, 12, 16};
IntArray Tet1BubbleStokes :: surf_ordering [ 4 ] = {
    {1, 2, 3,  9, 10, 11,  5,  6,  7},
    {1, 2, 3,  5,  6,  7, 13, 14, 15},
    {5, 6, 7,  9, 10, 11, 13, 14, 15},
    {1, 2, 3, 13, 14, 15,  9, 10, 11}
};

Tet1BubbleStokes :: Tet1BubbleStokes(int n, Domain *aDomain) : FMElement(n, aDomain), ZZNodalRecoveryModelInterface(this), SpatialLocalizerInterface(this)
{
    this->numberOfDofMans = 4;
    this->numberOfGaussPoints = 24;

    this->bubble.reset( new ElementDofManager(1, aDomain, this) );
    this->bubble->appendDof( new MasterDof(this->bubble.get(), V_u) );
    this->bubble->appendDof( new MasterDof(this->bubble.get(), V_v) );
    this->bubble->appendDof( new MasterDof(this->bubble.get(), V_w) );
}

Tet1BubbleStokes :: ~Tet1BubbleStokes()
{
}

void Tet1BubbleStokes :: computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], this->numberOfGaussPoints, this);
    }
}

int Tet1BubbleStokes :: computeNumberOfDofs()
{
    return 19;
}

void Tet1BubbleStokes :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {V_u, V_v, V_w, P_f};
}

void Tet1BubbleStokes :: giveInternalDofManDofIDMask(int i, IntArray &answer) const
{
    answer = {V_u, V_v, V_w};
}

double Tet1BubbleStokes :: computeVolumeAround(GaussPoint *gp)
{
    double detJ = fabs( this->interp.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    return detJ *gp->giveWeight();
}

void Tet1BubbleStokes :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
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

void Tet1BubbleStokes :: giveCharacteristicMatrix(FloatMatrix &answer,
                                                  CharType mtrx, TimeStep *tStep)
{
    // Compute characteristic matrix for this element. The only option is the stiffness matrix...
    if ( mtrx == TangentStiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, TangentStiffness, tStep);
    } else {
        OOFEM_ERROR("Unknown Type of characteristic mtrx.");
    }
}

void Tet1BubbleStokes :: computeInternalForcesVector(FloatArray &answer, TimeStep *tStep)
{
    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
    FloatArray a_pressure, a_velocity, devStress, epsp, N, dNv(15);
    double r_vol, pressure;
    FloatMatrix dN, B(6, 15);
    B.zero();

    this->computeVectorOfVelocities(VM_Total, tStep, a_velocity);
    this->computeVectorOfPressures(VM_Total, tStep, a_pressure);

    FloatArray momentum, conservation;

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();

        double detJ = fabs( this->interp.evaldNdx( dN, lcoords, FEIElementGeometryWrapper(this) ) );
        this->interp.evalN( N, lcoords, FEIElementGeometryWrapper(this) );
        double dV = detJ * gp->giveWeight();

        for ( int j = 0, k = 0; j < dN.giveNumberOfRows(); j++, k += 3 ) {
            dNv(k + 0) = B(0, k + 0) = B(5, k + 1) = B(4, k + 2) = dN(j, 0);
            dNv(k + 1) = B(1, k + 1) = B(5, k + 0) = B(3, k + 2) = dN(j, 1);
            dNv(k + 2) = B(2, k + 2) = B(4, k + 0) = B(3, k + 1) = dN(j, 2);
        }

        // Bubble contribution;
        dNv(12) = B(0, 12) = B(5, 13) = B(4, 14) = dN(0, 0) * N(1) * N(2) * N(3) + N(0) * dN(1, 0) * N(2) * N(3) + N(0) * N(1) * dN(2, 0) * N(3) + N(0) * N(1) * N(2) * dN(3, 0);
        dNv(13) = B(1, 13) = B(5, 12) = B(3, 14) = dN(0, 1) * N(1) * N(2) * N(3) + N(0) * dN(1, 1) * N(2) * N(3) + N(0) * N(1) * dN(2, 1) * N(3) + N(0) * N(1) * N(2) * dN(3, 1);
        dNv(14) = B(2, 14) = B(4, 12) = B(3, 13) = dN(0, 2) * N(1) * N(2) * N(3) + N(0) * dN(1, 2) * N(2) * N(3) + N(0) * N(1) * dN(2, 2) * N(3) + N(0) * N(1) * N(2) * dN(3, 2);

        pressure = N.dotProduct(a_pressure);
        epsp.beProductOf(B, a_velocity);

        mat->computeDeviatoricStressVector(devStress, r_vol, gp, epsp, pressure, tStep);

        momentum.plusProduct(B, devStress, dV);
        momentum.add(-pressure * dV, dNv);
        conservation.add(r_vol * dV, N);
    }

    answer.resize(19);
    answer.zero();
    answer.assemble(momentum, this->momentum_ordering);
    answer.assemble(conservation, this->conservation_ordering);
}

void Tet1BubbleStokes :: computeExternalForcesVector(FloatArray &answer, TimeStep *tStep)
{
    FloatArray vec;

    int nLoads = this->boundaryLoadArray.giveSize() / 2;
    answer.clear();
    for ( int i = 1; i <= nLoads; i++ ) {  // For each Neumann boundary condition
        int load_number = this->boundaryLoadArray.at(2 * i - 1);
        int load_id = this->boundaryLoadArray.at(2 * i);
        Load *load = this->domain->giveLoad(load_number);
        bcGeomType ltype = load->giveBCGeoType();

        if ( ltype == SurfaceLoadBGT ) {
            this->computeBoundaryLoadVector(vec, static_cast< BoundaryLoad * >(load), load_id, ExternalForcesVector, VM_Total, tStep);
        } else {
            OOFEM_ERROR("Unsupported boundary condition: %d", load_id);
        }

        answer.add(vec);
    }

    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nLoads; i++ ) {
        Load *load  = domain->giveLoad( bodyLoadArray.at(i) );
        bcGeomType ltype = load->giveBCGeoType();
        if ( ltype == BodyLoadBGT && load->giveBCValType() == ForceLoadBVT ) {
            this->computeLoadVector(vec, load, ExternalForcesVector, VM_Total, tStep);
            answer.add(vec);
        } else {
            OOFEM_ERROR("Unsupported body load: %d", load);
        }
    }
}


void Tet1BubbleStokes :: computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.clear();
        return;
    }

    FloatArray N, gVector, temparray(15);
    double dV, detJ, rho;

    // This is assumed to be the dead weight (thus multiplied by rho)
    load->computeComponentArrayAt(gVector, tStep, VM_Total);
    temparray.zero();
    if ( gVector.giveSize() ) {
        for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
            const FloatArray &lcoords = gp->giveNaturalCoordinates();

            rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);
            detJ = fabs( this->interp.giveTransformationJacobian( lcoords, FEIElementGeometryWrapper(this) ) );
            dV = detJ * gp->giveWeight() * rho;

            this->interp.evalN( N, lcoords, FEIElementGeometryWrapper(this) );
            for ( int j = 0; j < N.giveSize(); j++ ) {
                temparray(3 * j + 0) += N(j) * rho * gVector(0) * dV;
                temparray(3 * j + 1) += N(j) * rho * gVector(1) * dV;
                temparray(3 * j + 2) += N(j) * rho * gVector(2) * dV;
            }

            temparray(12) += N(0) * N(1) * N(2) * N(3) * rho * gVector(0) * dV;
            temparray(13) += N(0) * N(1) * N(2) * N(3) * rho * gVector(1) * dV;
            temparray(14) += N(0) * N(1) * N(2) * N(3) * rho * gVector(1) * dV;
        }
    }

    answer.resize(19);
    answer.zero();
    answer.assemble(temparray, this->momentum_ordering);
}


void Tet1BubbleStokes :: computeBoundaryLoadVector(FloatArray &answer, BoundaryLoad *load, int iSurf, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.clear();
        return;
    }

    if ( load->giveType() == TransmissionBC ) { // Neumann boundary conditions (traction)
        int numberOfIPs = ( int ) ceil( ( load->giveApproxOrder() + 2. ) / 2. );

        GaussIntegrationRule iRule(1, this, 1, 1);
        FloatArray N, t, f(9);

        f.zero();
        iRule.SetUpPointsOnTriangle(numberOfIPs, _Unknown);

        for ( GaussPoint *gp: iRule ) {
            const FloatArray &lcoords = gp->giveNaturalCoordinates();

            this->interp.surfaceEvalN( N, iSurf, lcoords, FEIElementGeometryWrapper(this) );
            double detJ = fabs( this->interp.surfaceGiveTransformationJacobian( iSurf, lcoords, FEIElementGeometryWrapper(this) ) );
            double dA = gp->giveWeight() * detJ;

            if ( load->giveFormulationType() == Load :: FT_Entity ) { // Edge load in xi-eta system
                load->computeValueAt(t, tStep, lcoords, VM_Total);
            } else { // Edge load in x-y system
                FloatArray gcoords;
                this->interp.boundaryLocal2Global( gcoords, iSurf, lcoords, FEIElementGeometryWrapper(this) );
                load->computeValueAt(t, tStep, gcoords, VM_Total);
            }

            // Reshape the vector
            for ( int j = 0; j < N.giveSize(); j++ ) {
                f(3 * j + 0) += N(j) * t(0) * dA;
                f(3 * j + 1) += N(j) * t(1) * dA;
                f(3 * j + 2) += N(j) * t(2) * dA;
            }

            this->interp.surfaceEvalNormal( N, iSurf, lcoords, FEIElementGeometryWrapper(this) );
        }

        answer.resize(19);
        answer.zero();
        answer.assemble(f, this->surf_ordering [ iSurf - 1 ]);
    } else {
        OOFEM_ERROR("Strange boundary condition type");
    }
}

void Tet1BubbleStokes :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep)
{
    // Note: Working with the components; [K, G+Dp; G^T+Dv^T, C] . [v,p]
    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
    FloatMatrix B(6, 15), EdB, K, G, Dp, DvT, C, Ed, dN;
    FloatArray dNv(15), N, Ep, Cd, tmpA, tmpB;
    double Cp;

    B.zero();

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        // Compute Gauss point and determinant at current element
        const FloatArray &lcoords = gp->giveNaturalCoordinates();

        double detJ = fabs( this->interp.evaldNdx( dN, lcoords, FEIElementGeometryWrapper(this) ) );
        double dV = detJ * gp->giveWeight();
        this->interp.evalN( N, lcoords, FEIElementGeometryWrapper(this) );

        for ( int j = 0, k = 0; j < 4; j++, k += 3 ) {
            dNv(k + 0) = B(0, k + 0) = B(5, k + 1) = B(4, k + 2) = dN(j, 0);
            dNv(k + 1) = B(1, k + 1) = B(5, k + 0) = B(3, k + 2) = dN(j, 1);
            dNv(k + 2) = B(2, k + 2) = B(4, k + 0) = B(3, k + 1) = dN(j, 2);
        }

        // Bubble contribution;
        dNv(12) = B(0, 12) = B(5, 13) = B(4, 14) = dN(0, 0) * N(1) * N(2) * N(3) + N(0) * dN(1, 0) * N(2) * N(3) + N(0) * N(1) * dN(2, 0) * N(3) + N(0) * N(1) * N(2) * dN(3, 0);
        dNv(13) = B(1, 13) = B(5, 12) = B(3, 14) = dN(0, 1) * N(1) * N(2) * N(3) + N(0) * dN(1, 1) * N(2) * N(3) + N(0) * N(1) * dN(2, 1) * N(3) + N(0) * N(1) * N(2) * dN(3, 1);
        dNv(14) = B(2, 14) = B(4, 12) = B(3, 13) = dN(0, 2) * N(1) * N(2) * N(3) + N(0) * dN(1, 2) * N(2) * N(3) + N(0) * N(1) * dN(2, 2) * N(3) + N(0) * N(1) * N(2) * dN(3, 2);

        // Computing the internal forces should have been done first.
        // dsigma_dev/deps_dev  dsigma_dev/dp  deps_vol/deps_dev  deps_vol/dp
        mat->giveStiffnessMatrices(Ed, Ep, Cd, Cp, mode, gp, tStep);

        EdB.beProductOf(Ed, B);
        K.plusProductSymmUpper(B, EdB, dV);
        G.plusDyadUnsym(dNv, N, -dV);
        C.plusDyadSymmUpper(N, Cp * dV);

        tmpA.beTProductOf(B, Ep);
        Dp.plusDyadUnsym(tmpA, N, dV);

        tmpB.beTProductOf(B, Cd);
        DvT.plusDyadUnsym(N, tmpB, dV);
    }

    K.symmetrized();
    C.symmetrized();
    FloatMatrix GTDvT, GDp;
    GTDvT.beTranspositionOf(G);
    GTDvT.add(DvT);
    GDp = G;
    GDp.add(Dp);

    answer.resize(19, 19);
    answer.zero();
    answer.assemble(K, this->momentum_ordering);
    answer.assemble(GTDvT, this->conservation_ordering, this->momentum_ordering);
    answer.assemble(GDp, this->momentum_ordering, this->conservation_ordering);
    answer.assemble(C, this->conservation_ordering);
}

FEInterpolation *Tet1BubbleStokes :: giveInterpolation() const
{
    return & interp;
}

FEInterpolation *Tet1BubbleStokes :: giveInterpolation(DofIDItem id) const
{
    return & interp;
}

void Tet1BubbleStokes :: updateYourself(TimeStep *tStep)
{
    Element :: updateYourself(tStep);
}

// Some extension Interfaces to follow:

Interface *Tet1BubbleStokes :: giveInterface(InterfaceType it)
{
    switch ( it ) {
    case ZZNodalRecoveryModelInterfaceType:
        return static_cast< ZZNodalRecoveryModelInterface * >(this);

    case SpatialLocalizerInterfaceType:
        return static_cast< SpatialLocalizerInterface * >(this);

    default:
        return FMElement :: giveInterface(it);
    }
}


void Tet1BubbleStokes :: computeField(ValueModeType mode, TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
{
    FloatArray n, n_lin, pressures, velocities;
    this->interp.evalN( n, lcoords, FEIElementGeometryWrapper(this) );
    this->interp.evalN( n_lin, lcoords, FEIElementGeometryWrapper(this) );
    this->computeVectorOf({P_f}, mode, tStep, pressures);
    this->computeVectorOf({V_u, V_v, V_w}, mode, tStep, velocities);

    answer.resize(4);
    answer.zero();
    for ( int i = 1; i <= n.giveSize(); i++ ) {
        answer(0) += n.at(i) * velocities.at(i*3-2);
        answer(1) += n.at(i) * velocities.at(i*3-1);
        answer(2) += n.at(i) * velocities.at(i*3);
    }
    answer(0) += n.at(1) * n.at(2) * n.at(3) * n.at(4) * velocities.at(13);
    answer(1) += n.at(1) * n.at(2) * n.at(3) * n.at(4) * velocities.at(14);
    answer(2) += n.at(1) * n.at(2) * n.at(3) * n.at(4) * velocities.at(15);

    for ( int i = 1; i <= n_lin.giveSize(); i++ ) {
        answer(3) += n_lin.at(i) * pressures.at(i);
    }
}

} // end namespace oofem
