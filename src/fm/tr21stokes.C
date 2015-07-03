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

#include "tr21stokes.h"
#include "node.h"
#include "dof.h"
#include "domain.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "bcgeomtype.h"
#include "load.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "fluiddynamicmaterial.h"
#include "fei2dtrlin.h"
#include "fei2dtrquad.h"
#include "fluidcrosssection.h"
#include "assemblercallback.h"
#include "classfactory.h"
#include "engngm.h"

namespace oofem {
REGISTER_Element(Tr21Stokes);

// Set up interpolation coordinates
FEI2dTrLin Tr21Stokes :: interpolation_lin(1, 2);
FEI2dTrQuad Tr21Stokes :: interpolation_quad(1, 2);
// Set up ordering vectors (for assembling)
IntArray Tr21Stokes :: momentum_ordering = {1, 2, 4, 5, 7, 8, 10, 11, 12, 13, 14, 15};
IntArray Tr21Stokes :: conservation_ordering = {3, 6, 9};
IntArray Tr21Stokes :: edge_ordering [ 3 ] = {
    {1, 2, 4, 5, 10, 11},
    {4, 5, 7, 8, 12, 13},
    {7, 8, 1, 2, 14, 15}
};

Tr21Stokes :: Tr21Stokes(int n, Domain *aDomain) : FMElement(n, aDomain), ZZNodalRecoveryModelInterface(this), SpatialLocalizerInterface(this)
{
    this->numberOfDofMans = 6;
    if ( aDomain->giveEngngModel()->giveProblemScale() == macroScale ) { // Pick the lowest default value for multiscale computations.
        this->numberOfGaussPoints = 3;
    } else {
        this->numberOfGaussPoints = 4;
    }
}

Tr21Stokes :: ~Tr21Stokes()
{ }

void Tr21Stokes :: computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], this->numberOfGaussPoints, this);
    }
}

int Tr21Stokes :: computeNumberOfDofs()
{
    return 15;
}

void Tr21Stokes :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    if ( inode <= 3 ) {
        answer = {V_u, V_v, P_f};
    } else {
        answer = {V_u, V_v};
    }
}

double Tr21Stokes :: computeVolumeAround(GaussPoint *gp)
{
    double detJ = fabs( this->interpolation_quad.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    return detJ *gp->giveWeight();
}

void Tr21Stokes :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
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

void Tr21Stokes :: giveCharacteristicMatrix(FloatMatrix &answer,
                                            CharType mtrx, TimeStep *tStep)
{
    // Compute characteristic matrix for this element. The only option is the stiffness matrix...
    if ( mtrx == TangentStiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, TangentStiffness, tStep);
    } else {
        OOFEM_ERROR("Unknown Type of characteristic mtrx.");
    }
}

void Tr21Stokes :: computeInternalForcesVector(FloatArray &answer, TimeStep *tStep)
{
    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
    FloatArray a_pressure, a_velocity, devStress, epsp, Nh, dNv(12);
    double r_vol, pressure;
    FloatMatrix dN, B(3, 12);
    B.zero();

    this->computeVectorOfVelocities(VM_Total, tStep, a_velocity);
    this->computeVectorOfPressures(VM_Total, tStep, a_pressure);

    FloatArray momentum, conservation;

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();

        double detJ = fabs( this->interpolation_quad.evaldNdx( dN, lcoords, FEIElementGeometryWrapper(this) ) );
        this->interpolation_lin.evalN( Nh, lcoords, FEIElementGeometryWrapper(this) );
        double dA = detJ * gp->giveWeight();

        for ( int j = 0, k = 0; j < dN.giveNumberOfRows(); j++, k += 2 ) {
            dNv(k)     = B(0, k)     = B(2, k + 1) = dN(j, 0);
            dNv(k + 1) = B(1, k + 1) = B(2, k)     = dN(j, 1);
        }

        pressure = Nh.dotProduct(a_pressure);
        epsp.beProductOf(B, a_velocity);

        mat->computeDeviatoricStressVector(devStress, r_vol, gp, epsp, pressure, tStep);

        momentum.plusProduct(B, devStress, dA);
        momentum.add(-pressure * dA, dNv);
        conservation.add(r_vol * dA, Nh);
    }

    answer.resize(15);
    answer.zero();
    answer.assemble(momentum, this->momentum_ordering);
    answer.assemble(conservation, this->conservation_ordering);
}

void Tr21Stokes :: computeExternalForcesVector(FloatArray &answer, TimeStep *tStep)
{
    FloatArray vec;

    answer.clear();

    int nLoads = this->boundaryLoadArray.giveSize() / 2;
    for ( int i = 1; i <= nLoads; i++ ) {  // For each Neumann boundary condition
        int load_number = this->boundaryLoadArray.at(2 * i - 1);
        int load_id = this->boundaryLoadArray.at(2 * i);
        Load *load = this->domain->giveLoad(load_number);
        bcGeomType ltype = load->giveBCGeoType();

        if ( ltype == EdgeLoadBGT ) {
            this->computeBoundaryLoadVector(vec, static_cast< BoundaryLoad * >(load), load_id, ExternalForcesVector, VM_Total, tStep);
            answer.add(vec);
        }
    }

    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nLoads; i++ ) {
        Load *load = domain->giveLoad( bodyLoadArray.at(i) );
        bcGeomType ltype = load->giveBCGeoType();
        if ( ltype == BodyLoadBGT && load->giveBCValType() == ForceLoadBVT ) {
            this->computeLoadVector(vec, load, ExternalForcesVector, VM_Total, tStep);
            answer.add(vec);
        }
    }
}

void Tr21Stokes :: computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.clear();
        return;
    }

    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
    FloatArray N, gVector, temparray(12);

    load->computeComponentArrayAt(gVector, tStep, VM_Total);
    temparray.zero();
    if ( gVector.giveSize() ) {
        for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
            const FloatArray &lcoords = gp->giveNaturalCoordinates();

            double rho = mat->give('d', gp);
            double detJ = fabs( this->interpolation_quad.giveTransformationJacobian( lcoords, FEIElementGeometryWrapper(this) ) );
            double dA = detJ * gp->giveWeight();

            this->interpolation_quad.evalN( N, lcoords, FEIElementGeometryWrapper(this) );
            for ( int j = 0; j < 6; j++ ) {
                temparray(2 * j)     += N(j) * rho * gVector(0) * dA;
                temparray(2 * j + 1) += N(j) * rho * gVector(1) * dA;
            }
        }
    }

    answer.resize(15);
    answer.zero();
    answer.assemble(temparray, this->momentum_ordering);
}

void Tr21Stokes :: computeBoundaryLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.clear();
        return;
    }

    if ( load->giveType() == TransmissionBC ) { // Neumann boundary conditions (traction)

        int numberOfEdgeIPs = ( int ) ceil( ( load->giveApproxOrder() + 1. ) / 2. ) * 2;

        GaussIntegrationRule iRule(1, this, 1, 1);
        FloatArray N, t, f(6);

        f.zero();
        iRule.SetUpPointsOnLine(numberOfEdgeIPs, _Unknown);

        for ( GaussPoint *gp: iRule ) {
            const FloatArray &lcoords = gp->giveNaturalCoordinates();

            this->interpolation_quad.edgeEvalN( N, boundary, lcoords, FEIElementGeometryWrapper(this) );
            double detJ = fabs( this->interpolation_quad.boundaryGiveTransformationJacobian( boundary, lcoords, FEIElementGeometryWrapper(this) ) );
            double dS = gp->giveWeight() * detJ;

            if ( load->giveFormulationType() == Load :: FT_Entity ) { // Edge load in xi-eta system
                load->computeValueAt(t, tStep, lcoords, VM_Total);
            } else { // Edge load in x-y system
                FloatArray gcoords;
                this->interpolation_quad.boundaryLocal2Global( gcoords, boundary, lcoords, FEIElementGeometryWrapper(this) );
                load->computeValueAt(t, tStep, gcoords, VM_Total);
            }

            // Reshape the vector
            for ( int j = 0; j < 3; j++ ) {
                f(2 * j)     += N(j) * t(0) * dS;
                f(2 * j + 1) += N(j) * t(1) * dS;
            }
        }

        answer.resize(15);
        answer.zero();
        answer.assemble(f, this->edge_ordering [ boundary - 1 ]);
    } else {
        OOFEM_ERROR("Strange boundary condition type");
    }
}

void Tr21Stokes :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep)
{
    // Note: Working with the components; [K, G+Dp; G^T+Dv^T, C] . [v,p]
    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
    FloatMatrix B(3, 12), EdB, K, G, Dp, DvT, C, Ed, dN;
    FloatArray dN_V(12), Nlin, Ep, Cd, tmpA, tmpB;
    double Cp;

    B.zero();

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        // Compute Gauss point and determinant at current element
        const FloatArray &lcoords = gp->giveNaturalCoordinates();

        this->interpolation_lin.evalN( Nlin, lcoords, FEIElementGeometryWrapper(this) );
        double detJ = fabs( this->interpolation_quad.evaldNdx( dN, lcoords, FEIElementGeometryWrapper(this) ) );
        double dA = detJ * gp->giveWeight();

        for ( int j = 0, k = 0; j < 6; j++, k += 2 ) {
            dN_V(k)     = B(0, k)     = B(2, k + 1) = dN(j, 0);
            dN_V(k + 1) = B(1, k + 1) = B(2, k)     = dN(j, 1);
        }

        // Computing the internal forces should have been done first.
        // dsigma_dev/deps_dev  dsigma_dev/dp  deps_vol/deps_dev  deps_vol/dp
        mat->giveStiffnessMatrices(Ed, Ep, Cd, Cp, mode, gp, tStep);

        EdB.beProductOf(Ed, B);
        K.plusProductSymmUpper(B, EdB, dA);
        G.plusDyadUnsym(dN_V, Nlin, -dA);
        C.plusDyadSymmUpper(Nlin, Cp * dA);

        tmpA.beTProductOf(B, Ep);
        Dp.plusDyadUnsym(tmpA, Nlin, dA);

        tmpB.beTProductOf(B, Cd);
        DvT.plusDyadUnsym(Nlin, tmpB, dA);
    }

    K.symmetrized();
    C.symmetrized();
    FloatMatrix GTDvT, GDp;
    GTDvT.beTranspositionOf(G);
    GTDvT.add(DvT);
    GDp = G;
    GDp.add(Dp);

    answer.resize(15, 15);
    answer.zero();
    answer.assemble(K, this->momentum_ordering);
    answer.assemble(GTDvT, this->conservation_ordering, this->momentum_ordering);
    answer.assemble(GDp, this->momentum_ordering, this->conservation_ordering);
    answer.assemble(C, this->conservation_ordering);
}

FEInterpolation *Tr21Stokes :: giveInterpolation() const
{
    return & interpolation_quad;
}

FEInterpolation *Tr21Stokes :: giveInterpolation(DofIDItem id) const
{
    if ( id == P_f ) {
        return & interpolation_lin;
    } else {
        return & interpolation_quad;
    }
}

void Tr21Stokes :: updateYourself(TimeStep *tStep)
{
    Element :: updateYourself(tStep);
}

// Some extension Interfaces to follow:

Interface *Tr21Stokes :: giveInterface(InterfaceType it)
{
    switch ( it ) {
    case NodalAveragingRecoveryModelInterfaceType:
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);

    case ZZNodalRecoveryModelInterfaceType:
        return static_cast< ZZNodalRecoveryModelInterface * >(this);

    case SpatialLocalizerInterfaceType:
        return static_cast< SpatialLocalizerInterface * >(this);

    default:
        return FMElement :: giveInterface(it);
    }
}

void Tr21Stokes :: computeField(ValueModeType mode, TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
{
    FloatArray n, n_lin, pressures, velocities;
    this->interpolation_quad.evalN( n, lcoords, FEIElementGeometryWrapper(this) );
    this->interpolation_lin.evalN( n_lin, lcoords, FEIElementGeometryWrapper(this) );
    this->computeVectorOf({P_f}, mode, tStep, pressures);
    this->computeVectorOf({V_u, V_v, V_w}, mode, tStep, velocities);

    answer.resize(4);
    answer.zero();
    for ( int i = 1; i <= n.giveSize(); i++ ) {
        answer(0) += n.at(i) * velocities.at(i*2-1);
        answer(1) += n.at(i) * velocities.at(i*2);
    }

    for ( int i = 1; i <= n_lin.giveSize(); i++ ) {
        answer(2) += n_lin.at(i) * pressures.at(i);
    }
}


void Tr21Stokes :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_Pressure ) {
        answer.resize(1);
        if ( node == 1 || node == 2 || node == 3 ) {
            answer.at(1) = this->giveNode(node)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep);
        } else {
            double a, b;
            if ( node == 4 ) {
                a = this->giveNode(1)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep);
                b = this->giveNode(2)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep);
            } else if ( node == 5 ) {
                a = this->giveNode(2)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep);
                b = this->giveNode(3)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep);
            } else { /*if ( node == 6 )*/
                a = this->giveNode(3)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep);
                b = this->giveNode(1)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep);
            }

            answer.at(1) = ( a + b ) / 2;
        }
    } else {
        answer.clear();
    }
}

} // end namespace oofem
