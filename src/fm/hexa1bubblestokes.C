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

#include "hexa1bubblestokes.h"
#include "node.h"
#include "domain.h"
#include "equationid.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "bcgeomtype.h"
#include "load.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "fluiddynamicmaterial.h"
#include "fei3dhexalin.h"
#include "masterdof.h"
#include "crosssection.h"
#include "classfactory.h"

namespace oofem {

REGISTER_Element( Hexa1BubbleStokes );

FEI3dHexaLin Hexa1BubbleStokes :: interp;
// Set up ordering vectors (for assembling)
IntArray Hexa1BubbleStokes :: momentum_ordering(27);
IntArray Hexa1BubbleStokes :: conservation_ordering(4);
IntArray Hexa1BubbleStokes :: surf_ordering [ 6 ] = { IntArray(12), IntArray(12), IntArray(12), IntArray(12), IntArray(12), IntArray(12) };
bool Hexa1BubbleStokes :: __initialized = Hexa1BubbleStokes :: initOrdering();

Hexa1BubbleStokes :: Hexa1BubbleStokes(int n, Domain *aDomain) : FMElement(n, aDomain)
{
    this->numberOfDofMans = 8;
    this->numberOfGaussPoints = 27;

    this->bubble = new ElementDofManager(1, aDomain, this);
    this->bubble->appendDof(new MasterDof(1, this->bubble, V_u));
    this->bubble->appendDof(new MasterDof(2, this->bubble, V_v));
    this->bubble->appendDof(new MasterDof(3, this->bubble, V_w));
}

Hexa1BubbleStokes :: ~Hexa1BubbleStokes()
{
    delete this->bubble;
}

IRResultType Hexa1BubbleStokes :: initializeFrom(InputRecord *ir)
{
    this->FMElement :: initializeFrom(ir);
    return IRRT_OK;
}

void Hexa1BubbleStokes :: computeGaussPoints()
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 2 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        this->giveCrossSection()->setupIntegrationPoints( *integrationRulesArray[0], numberOfGaussPoints, this );
    }
}

int Hexa1BubbleStokes :: computeNumberOfDofs(EquationID ut)
{
    if ( ut == EID_MomentumBalance_ConservationEquation ) {
        return 35;
    } else if ( ut == EID_MomentumBalance ) {
        return 27;
    } else if ( ut == EID_ConservationEquation ) {
        return 8;
    }
    return 0;
}

void Hexa1BubbleStokes :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    if ( ut == EID_MomentumBalance ) {
        answer.setValues(3, V_u, V_v, V_w);
    } else if ( ut == EID_ConservationEquation ) {
        answer.setValues(1, P_f);
    } else if ( ut == EID_MomentumBalance_ConservationEquation ) {
        answer.setValues(4, V_u, V_v, V_w, P_f);
    } else {
        answer.resize(0);
    }
}

void Hexa1BubbleStokes :: giveInternalDofManDofIDMask(int i, EquationID eid, IntArray &answer) const
{
    if ( eid == EID_MomentumBalance_ConservationEquation || eid == EID_MomentumBalance ) {
        answer.setValues(3, V_u, V_v, V_w);
    } else {
        answer.resize(0);
    }
}

double Hexa1BubbleStokes :: computeVolumeAround(GaussPoint *gp)
{
    double detJ = fabs(this->interp.giveTransformationJacobian(*gp->giveCoordinates(), FEIElementGeometryWrapper(this)));
    return detJ * gp->giveWeight();
}

void Hexa1BubbleStokes :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
                                            TimeStep *tStep)
{
    // Compute characteristic vector for this element. I.e the load vector(s)
    if ( mtrx == ExternalForcesVector ) {
        this->computeExternalForcesVector(answer, tStep);
    } else if ( mtrx == InternalForcesVector ) {
        this->computeInternalForcesVector(answer, tStep);
    } else {
        OOFEM_ERROR("giveCharacteristicVector: Unknown Type of characteristic mtrx.");
    }
}

void Hexa1BubbleStokes :: giveCharacteristicMatrix(FloatMatrix &answer,
                                            CharType mtrx, TimeStep *tStep)
{
    // Compute characteristic matrix for this element. The only option is the stiffness matrix...
    if ( mtrx == StiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, tStep);
    } else {
        OOFEM_ERROR("giveCharacteristicMatrix: Unknown Type of characteristic mtrx.");
    }
}

void Hexa1BubbleStokes :: computeInternalForcesVector(FloatArray &answer, TimeStep *tStep)
{
    IntegrationRule *iRule = integrationRulesArray [ 0 ];
    FluidDynamicMaterial *mat = static_cast< FluidDynamicMaterial * >( this->giveMaterial() );
    FloatArray a_pressure, a_velocity, devStress, epsp, N, dNv(27);
    double r_vol, pressure;
    FloatMatrix dN, B(6, 27);
    B.zero();
    
    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, a_velocity);
    this->computeVectorOf(EID_ConservationEquation, VM_Total, tStep, a_pressure);

    FloatArray momentum, conservation;

    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(i);
        FloatArray *lcoords = gp->giveCoordinates();

        double detJ = fabs( this->interp.evaldNdx(dN, * lcoords, FEIElementGeometryWrapper(this)) );
        this->interp.evalN(N, * lcoords, FEIElementGeometryWrapper(this));
        double dV = detJ * gp->giveWeight();
                
        for ( int j = 0, k = 0; j < dN.giveNumberOfRows(); j++, k += 3 ) {
            dNv(k + 0) = B(0, k + 0) = B(5, k + 1) = B(4, k + 2) = dN(j, 0);
            dNv(k + 1) = B(1, k + 1) = B(5, k + 0) = B(3, k + 2) = dN(j, 1);
            dNv(k + 2) = B(2, k + 2) = B(4, k + 0) = B(3, k + 1) = dN(j, 2);
        }

        // Bubble contribution;
        double b1 = 0., b2 = 0., b3 = 0.;
        for ( int j = 0; j < 8; ++j ) {
            double x = 16777216.;
            for ( int k = 0; k < 8; ++k ) {
                if ( k != j ) x *= N(k);
            }
            b1 += dN(j, 0)*x;
            b2 += dN(j, 1)*x;
            b3 += dN(j, 2)*x;
        }
        // The bubble we define the bubble function as := prod(N_i) * 8^8 which gives is roughly the same order of magnitude.
        dNv(24) = B(0, 24) = B(5, 25) = B(4, 26) = b1;
        dNv(25) = B(1, 25) = B(5, 24) = B(3, 26) = b2;
        dNv(26) = B(2, 26) = B(4, 24) = B(3, 25) = b3;

        pressure = N.dotProduct(a_pressure);
        epsp.beProductOf(B, a_velocity);

        mat->computeDeviatoricStressVector(devStress, r_vol, gp, epsp, pressure, tStep);

        momentum.plusProduct(B, devStress, dV);
        momentum.add(-pressure*dV, dNv);
        conservation.add(r_vol*dV, N);
    }

    answer.resize(35);
    answer.zero();
    answer.assemble(momentum, this->momentum_ordering);
    answer.assemble(conservation, this->conservation_ordering);
}

void Hexa1BubbleStokes :: computeExternalForcesVector(FloatArray &answer, TimeStep *tStep)
{
    FloatArray vec;

    int nLoads = this->boundaryLoadArray.giveSize() / 2;
    answer.resize(0);
    for ( int i = 1; i <= nLoads; i++ ) {  // For each Neumann boundary condition
        int load_number = this->boundaryLoadArray.at(2 * i - 1);
        int load_id = this->boundaryLoadArray.at(2 * i);
        Load *load = this->domain->giveLoad(load_number);
        bcGeomType ltype = load->giveBCGeoType();

        if ( ltype == SurfaceLoadBGT ) {
            this->computeBoundaryLoadVector(vec, static_cast<BoundaryLoad*>(load), load_id, ExternalForcesVector, VM_Total, tStep);
        } else {
            OOFEM_ERROR2("Hexa1BubbleStokes :: computeLoadVector - Unsupported boundary condition: %d", load_id);
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
            OOFEM_ERROR2("Hexa1BubbleStokes :: computeLoadVector - Unsupported body load: %d", load);
        }
    }
}


void Hexa1BubbleStokes :: computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.resize(0);
        return;
    }

    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    FloatArray N, gVector, temparray(27);

    // This is assumed to be the dead weight (thus multiplied by rho)
    load->computeComponentArrayAt(gVector, tStep, VM_Total);
    temparray.zero();
    if ( gVector.giveSize() ) {
        for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(k);
            FloatArray *lcoords = gp->giveCoordinates();

            double rho = this->giveMaterial()->give('d', gp);
            double detJ = fabs( this->interp.giveTransformationJacobian(* lcoords, FEIElementGeometryWrapper(this)) );
            double dV = detJ * gp->giveWeight() * rho;

            this->interp.evalN(N, * lcoords, FEIElementGeometryWrapper(this));
            for ( int j = 0; j < N.giveSize(); j++ ) {
                temparray(3 * j + 0) += N(j) * rho * gVector(0) * dV;
                temparray(3 * j + 1) += N(j) * rho * gVector(1) * dV;
                temparray(3 * j + 2) += N(j) * rho * gVector(2) * dV;
            }
            temparray(24) += N(0)*N(1)*N(2)*N(3) * rho * gVector(0) * dV;
            temparray(25) += N(0)*N(1)*N(2)*N(3) * rho * gVector(1) * dV;
            temparray(26) += N(0)*N(1)*N(2)*N(3) * rho * gVector(1) * dV;
        }
    }

    answer.resize(35);
    answer.zero();
    answer.assemble( temparray, this->momentum_ordering );
}


void Hexa1BubbleStokes :: computeBoundaryLoadVector(FloatArray &answer, BoundaryLoad *load, int iSurf, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.resize(0);
        return;
    }

    if ( load->giveType() == TransmissionBC ) { // Neumann boundary conditions (traction)
        BoundaryLoad *boundaryLoad = static_cast< BoundaryLoad * >( load );

        int numberOfIPs = ( int ) ceil( ( boundaryLoad->giveApproxOrder() + 2. ) / 2. );

        GaussIntegrationRule iRule(1, this, 1, 1);
        GaussPoint *gp;
        FloatArray N, t, f(9);
        IntArray edge_mapping;

        f.zero();
        iRule.SetUpPointsOnTriangle(numberOfIPs, _Unknown);

        for ( int i = 0; i < iRule.giveNumberOfIntegrationPoints(); i++ ) {
            gp = iRule.getIntegrationPoint(i);
            FloatArray *lcoords = gp->giveCoordinates();

            this->interp.surfaceEvalN(N, iSurf, * lcoords, FEIElementGeometryWrapper(this));
            double detJ = fabs(this->interp.surfaceGiveTransformationJacobian(iSurf, * lcoords, FEIElementGeometryWrapper(this)));
            double dA = gp->giveWeight() * detJ;

            if ( boundaryLoad->giveFormulationType() == BoundaryLoad :: BL_EntityFormulation ) { // Edge load in xi-eta system
                boundaryLoad->computeValueAt(t, tStep, * lcoords, VM_Total);
            } else { // Edge load in x-y system
                FloatArray gcoords;
                this->interp.boundaryLocal2Global(gcoords, iSurf, * lcoords, FEIElementGeometryWrapper(this));
                boundaryLoad->computeValueAt(t, tStep, gcoords, VM_Total);
            }

            // Reshape the vector
            for ( int j = 0; j < N.giveSize(); j++ ) {
                f(3 * j + 0) += N(j) * t(0) * dA;
                f(3 * j + 1) += N(j) * t(1) * dA;
                f(3 * j + 2) += N(j) * t(2) * dA;
            }

            this->interp.surfaceEvalNormal(N, iSurf, * lcoords, FEIElementGeometryWrapper(this));
        }

        answer.resize(35);
        answer.zero();
        answer.assemble(f, this->surf_ordering [ iSurf - 1 ]);
    } else {
        OOFEM_ERROR("Hexa1BubbleStokes :: computeBoundaryLoadVector - Strange boundary condition type");
    }
}

void Hexa1BubbleStokes :: computeStiffnessMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    // Note: Working with the components; [K, G+Dp; G^T+Dv^T, C] . [v,p]
    FluidDynamicMaterial *mat = static_cast< FluidDynamicMaterial * >( this->giveMaterial() );
    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    FloatMatrix B(6, 27), EdB, K, G, Dp, DvT, C, Ed, dN;
    FloatArray dNv(27), N, Ep, Cd, tmpA, tmpB;
    double Cp;
    B.zero();

    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        // Compute Gauss point and determinant at current element
        GaussPoint *gp = iRule->getIntegrationPoint(i);
        FloatArray *lcoords = gp->giveCoordinates();

        double detJ = fabs( this->interp.evaldNdx(dN, * lcoords, FEIElementGeometryWrapper(this)) );
        double dV = detJ * gp->giveWeight();
        this->interp.evalN(N, * lcoords, FEIElementGeometryWrapper(this));
 
        for ( int j = 0, k = 0; j < dN.giveNumberOfRows(); j++, k += 3 ) {
            dNv(k + 0) = B(0, k + 0) = B(5, k + 1) = B(4, k + 2) = dN(j, 0);
            dNv(k + 1) = B(1, k + 1) = B(5, k + 0) = B(3, k + 2) = dN(j, 1);
            dNv(k + 2) = B(2, k + 2) = B(4, k + 0) = B(3, k + 1) = dN(j, 2);
        }

        // Bubble contribution;
        double b1 = 0., b2 = 0., b3 = 0.;
        for ( int j = 0; j < 8; ++j ) {
            double x = 16777216.;
            for ( int k = 0; k < 8; ++k ) {
                if ( k != j ) x *= N(k);
            }
            b1 += dN(j, 0)*x;
            b2 += dN(j, 1)*x;
            b3 += dN(j, 2)*x;
        }
        b1 = b2 = b3 = 0.;
        dNv(24) = B(0, 24) = B(5, 25) = B(4, 26) = b1;
        dNv(25) = B(1, 25) = B(5, 24) = B(3, 26) = b2;
        dNv(26) = B(2, 26) = B(4, 24) = B(3, 25) = b3;

        // Computing the internal forces should have been done first.
        mat->giveDeviatoricStiffnessMatrix(Ed, TangentStiffness, gp, tStep); // dsigma_dev/deps_dev
        mat->giveDeviatoricPressureStiffness(Ep, TangentStiffness, gp, tStep); // dsigma_dev/dp
        mat->giveVolumetricDeviatoricStiffness(Cd, TangentStiffness, gp, tStep); // deps_vol/deps_dev
        mat->giveVolumetricPressureStiffness(Cp, TangentStiffness, gp, tStep); // deps_vol/dp

        EdB.beProductOf(Ed,B);
        K.plusProductSymmUpper(B, EdB, dV);
        G.plusDyadUnsym(dNv, N, -dV);
        C.plusDyadSymmUpper(N, Cp*dV);

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

    answer.resize(35, 35);
    answer.zero();
    answer.assemble(K, this->momentum_ordering);
    answer.assemble(GTDvT, this->conservation_ordering, this->momentum_ordering);
    answer.assemble(GDp, this->momentum_ordering, this->conservation_ordering);
    answer.assemble(C, this->conservation_ordering);
}

FEInterpolation *Hexa1BubbleStokes :: giveInterpolation() const
{
    return &interp;
}

FEInterpolation *Hexa1BubbleStokes :: giveInterpolation(DofIDItem id) const
{
    return &interp;
}

void Hexa1BubbleStokes :: updateYourself(TimeStep *tStep)
{
    Element :: updateYourself(tStep);
}

// Some extension Interfaces to follow:

Interface *Hexa1BubbleStokes :: giveInterface(InterfaceType it)
{
    switch (it) {
        case ZZNodalRecoveryModelInterfaceType:
            return static_cast< ZZNodalRecoveryModelInterface * >(this);
        case SpatialLocalizerInterfaceType:
            return static_cast< SpatialLocalizerInterface * >(this);
        case EIPrimaryUnknownMapperInterfaceType:
            return static_cast< EIPrimaryUnknownMapperInterface * >(this);
        default:
            return FMElement :: giveInterface(it);
    }
}

int Hexa1BubbleStokes :: SpatialLocalizerI_containsPoint(const FloatArray &coords)
{
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}

void Hexa1BubbleStokes :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType mode,
        TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
{
    FloatArray n, n_lin;
    this->interp.evalN(n, lcoords, FEIElementGeometryWrapper(this));
    this->interp.evalN(n_lin, lcoords, FEIElementGeometryWrapper(this));
    answer.resize(4);
    answer.zero();
    for (int i = 1; i <= n.giveSize(); i++) {
        answer(0) += n.at(i)*this->giveNode(i)->giveDofWithID(V_u)->giveUnknown(mode, tStep);
        answer(1) += n.at(i)*this->giveNode(i)->giveDofWithID(V_v)->giveUnknown(mode, tStep);
        answer(2) += n.at(i)*this->giveNode(i)->giveDofWithID(V_w)->giveUnknown(mode, tStep);
    }
    for (int i = 1; i <= n_lin.giveSize(); i++) {
        answer(3) += n_lin.at(i)*this->giveNode(i)->giveDofWithID(P_f)->giveUnknown(mode, tStep);
    }
}

int Hexa1BubbleStokes :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode, TimeStep *tStep, const FloatArray &gcoords, FloatArray &answer)
{
    FloatArray lcoords, n, n_lin;
    bool ok = this->computeLocalCoordinates(lcoords, gcoords);
    if (!ok) {
        answer.resize(0);
        return false;
    }
    this->EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(mode, tStep, lcoords, answer);
    return true;
}

void Hexa1BubbleStokes :: EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer)
{
    answer.setValues(4, V_u, V_v, V_w, P_f);
}

double Hexa1BubbleStokes :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray center;
    FloatArray lcoords;
    lcoords.setValues(3, 0., 0., 0.);
    this->interp.local2global(center, lcoords, FEIElementGeometryWrapper(this));
    return center.distance(coords);
}

} // end namespace oofem
