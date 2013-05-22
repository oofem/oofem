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

#include "tet1bubblestokes.h"
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
#include "fei3dtetlin.h"
#include "masterdof.h"
#include "classfactory.h"

namespace oofem {

REGISTER_Element( Tet1BubbleStokes );

FEI3dTetLin Tet1BubbleStokes :: interp;
// Set up ordering vectors (for assembling)
IntArray Tet1BubbleStokes :: momentum_ordering(15);
IntArray Tet1BubbleStokes :: conservation_ordering(4);
IntArray Tet1BubbleStokes :: edge_ordering [ 6 ] = { IntArray(6), IntArray(6), IntArray(6), IntArray(6), IntArray(6), IntArray(6) };
IntArray Tet1BubbleStokes :: surf_ordering [ 4 ] = { IntArray(9), IntArray(9), IntArray(9), IntArray(9) };
bool Tet1BubbleStokes :: __initialized = Tet1BubbleStokes :: initOrdering();

Tet1BubbleStokes :: Tet1BubbleStokes(int n, Domain *aDomain) : FMElement(n, aDomain)
{
    this->numberOfDofMans = 4;
    this->numberOfGaussPoints = 24;

    this->bubble = new ElementDofManager(1, aDomain, this);
    this->bubble->appendDof(new MasterDof(1, this->bubble, V_u));
    this->bubble->appendDof(new MasterDof(2, this->bubble, V_v));
    this->bubble->appendDof(new MasterDof(3, this->bubble, V_w));
}

Tet1BubbleStokes :: ~Tet1BubbleStokes()
{
    delete this->bubble;
}

IRResultType Tet1BubbleStokes :: initializeFrom(InputRecord *ir)
{
    this->FMElement :: initializeFrom(ir);
    return IRRT_OK;
}

void Tet1BubbleStokes :: computeGaussPoints()
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 2 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Tetrahedra, this->numberOfGaussPoints, _3dFlow);
    }
}

int Tet1BubbleStokes :: computeNumberOfDofs(EquationID ut)
{
    if ( ut == EID_MomentumBalance_ConservationEquation ) {
        return 19;
    } else if ( ut == EID_MomentumBalance ) {
        return 15;
    } else if ( ut == EID_ConservationEquation ) {
        return 4;
    }
    return 0;
}

void Tet1BubbleStokes :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
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

void Tet1BubbleStokes :: giveInternalDofManDofIDMask(int i, EquationID eid, IntArray &answer) const
{
    if ( eid == EID_MomentumBalance_ConservationEquation || eid == EID_MomentumBalance ) {
        answer.setValues(3, V_u, V_v, V_w);
    } else {
        answer.resize(0);
    }
}

double Tet1BubbleStokes :: computeVolumeAround(GaussPoint *gp)
{
    double detJ = fabs(this->interp.giveTransformationJacobian(*gp->giveCoordinates(), FEIElementGeometryWrapper(this)));
    return detJ * gp->giveWeight();
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
        OOFEM_ERROR("giveCharacteristicVector: Unknown Type of characteristic mtrx.");
    }
}

void Tet1BubbleStokes :: giveCharacteristicMatrix(FloatMatrix &answer,
                                            CharType mtrx, TimeStep *tStep)
{
    // Compute characteristic matrix for this element. The only option is the stiffness matrix...
    if ( mtrx == StiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, tStep);
    } else {
        OOFEM_ERROR("giveCharacteristicMatrix: Unknown Type of characteristic mtrx.");
    }
}

void Tet1BubbleStokes :: computeInternalForcesVector(FloatArray &answer, TimeStep *tStep)
{
    IntegrationRule *iRule = integrationRulesArray [ 0 ];
    FluidDynamicMaterial *mat = static_cast< FluidDynamicMaterial * >( this->giveMaterial() );
    FloatArray a_pressure, a_velocity, devStress, epsp, BTs, N, dNv(15);
    double r_vol, pressure;
    FloatMatrix dN, B(6, 15);
    B.zero();
    
    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, a_velocity);
    this->computeVectorOf(EID_ConservationEquation, VM_Total, tStep, a_pressure);

    FloatArray momentum(15), conservation(4);
    momentum.zero();
    conservation.zero();

    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(i);
        FloatArray *lcoords = gp->giveCoordinates();

        double detJ = fabs(this->interp.giveTransformationJacobian(* lcoords, FEIElementGeometryWrapper(this)));
        this->interp.evaldNdx(dN, * lcoords, FEIElementGeometryWrapper(this));
        this->interp.evalN(N, * lcoords, FEIElementGeometryWrapper(this));
        double dV = detJ * gp->giveWeight();

        for ( int j = 0, k = 0; j < 4; j++, k += 3 ) {
            dNv(k + 0) = B(0, k + 0) = B(4, k + 2) = B(5, k + 1) = dN(j, 0);
            dNv(k + 1) = B(1, k + 1) = B(3, k + 2) = B(5, k + 0) = dN(j, 1);
            dNv(k + 2) = B(2, k + 2) = B(3, k + 1) = B(4, k + 0) = dN(j, 2);
        }
        // Bubble contribution;
        dNv(12) = B(0,12) = B(4,14) = B(5,13) = dN(0,0)*N(1)*N(2)*N(3) + N(0)*dN(1,0)*N(2)*N(3) + N(0)*N(1)*dN(2,0)*N(3) + N(0)*N(1)*N(2)*dN(3,0);
        dNv(13) = B(1,13) = B(3,14) = B(5,12) = dN(0,1)*N(1)*N(2)*N(3) + N(0)*dN(1,1)*N(2)*N(3) + N(0)*N(1)*dN(2,1)*N(3) + N(0)*N(1)*N(2)*dN(3,1);
        dNv(14) = B(2,14) = B(3,13) = B(4,12) = dN(0,2)*N(1)*N(2)*N(3) + N(0)*dN(1,2)*N(2)*N(3) + N(0)*N(1)*dN(2,2)*N(3) + N(0)*N(1)*N(2)*dN(3,2);

        pressure = N.dotProduct(a_pressure);
        epsp.beProductOf(B, a_velocity);

        mat->computeDeviatoricStressVector(devStress, r_vol, gp, epsp, pressure, tStep);
        BTs.beTProductOf(B, devStress);

        momentum.add(dV, BTs);
        momentum.add(-pressure*dV, dNv);
        conservation.add(r_vol*dV, N);
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
    answer.resize(0);
    for ( int i = 1; i <= nLoads; i++ ) {  // For each Neumann boundary condition
        int load_number = this->boundaryLoadArray.at(2 * i - 1);
        int load_id = this->boundaryLoadArray.at(2 * i);
        Load *load = this->domain->giveLoad(load_number);
        bcGeomType ltype = load->giveBCGeoType();

        if ( ltype == SurfaceLoadBGT ) {
            this->computeBoundaryLoadVector(vec, load, load_id, ExternalForcesVector, VM_Total, tStep);
        } else if ( ltype == EdgeLoadBGT ) {
            this->computeEdgeLoadVector(vec, load, load_id, ExternalForcesVector, VM_Total, tStep);
        } else {
            OOFEM_ERROR2("Tet1BubbleStokes :: computeLoadVector - Unsupported boundary condition: %d", load_id);
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
            OOFEM_ERROR2("Tet1BubbleStokes :: computeLoadVector - Unsupported body load: %d", load);
        }
    }
}


void Tet1BubbleStokes :: computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.resize(0);
        return;
    }

    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    GaussPoint *gp;
    FloatArray N, gVector, *lcoords, temparray(15);
    double dV, detJ, rho;

    // This is assumed to be the dead weight (thus multiplied by rho)
    load->computeComponentArrayAt(gVector, tStep, VM_Total);
    temparray.zero();
    if ( gVector.giveSize() ) {
        for ( int k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
            gp = iRule->getIntegrationPoint(k);
            lcoords = gp->giveCoordinates();

            rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, gp, tStep);
            detJ = fabs( this->interp.giveTransformationJacobian(* lcoords, FEIElementGeometryWrapper(this)) );
            dV = detJ * gp->giveWeight() * rho;

            this->interp.evalN(N, * lcoords, FEIElementGeometryWrapper(this));
            for ( int j = 0; j < 4; j++ ) {
                temparray(3 * j + 0) += N(j) * rho * gVector(0) * dV;
                temparray(3 * j + 1) += N(j) * rho * gVector(1) * dV;
                temparray(3 * j + 2) += N(j) * rho * gVector(2) * dV;
            }
            temparray(12) += N(0)*N(1)*N(2)*N(3) * rho * gVector(0) * dV;
            temparray(13) += N(0)*N(1)*N(2)*N(3) * rho * gVector(1) * dV;
            temparray(14) += N(0)*N(1)*N(2)*N(3) * rho * gVector(1) * dV;
        }
    }

    answer.resize(19);
    answer.zero();
    answer.assemble( temparray, this->momentum_ordering );
}

void Tet1BubbleStokes :: computeEdgeLoadVector(FloatArray &answer, Load *load, int iEdge, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.resize(0);
        return;
    }

    if ( load->giveType() == TransmissionBC ) { // Neumann boundary conditions (traction)
        BoundaryLoad *boundaryLoad = static_cast< BoundaryLoad * >( load );

        int numberOfEdgeIPs = ( int ) ceil( ( boundaryLoad->giveApproxOrder() + 2. ) / 2. );

        GaussIntegrationRule iRule(1, this, 1, 1);
        GaussPoint *gp;
        FloatArray N, t, f(6);
        IntArray edge_mapping;

        f.zero();
        iRule.setUpIntegrationPoints(_Line, numberOfEdgeIPs, _Unknown);

        for ( int i = 0; i < iRule.getNumberOfIntegrationPoints(); i++ ) {
            gp = iRule.getIntegrationPoint(i);
            FloatArray *lcoords = gp->giveCoordinates();

            this->interp.edgeEvalN(N, iEdge, * lcoords, FEIElementGeometryWrapper(this));
            double detJ = fabs(this->interp.edgeGiveTransformationJacobian(iEdge, * lcoords, FEIElementGeometryWrapper(this)));
            double dS = gp->giveWeight() * detJ;

            if ( boundaryLoad->giveFormulationType() == BoundaryLoad :: BL_EntityFormulation ) { // Edge load in xi-eta system
                boundaryLoad->computeValueAt(t, tStep, * lcoords, VM_Total);
            } else { // Edge load in x-y system
                FloatArray gcoords;
                this->interp.edgeLocal2global(gcoords, iEdge, * lcoords, FEIElementGeometryWrapper(this));
                boundaryLoad->computeValueAt(t, tStep, gcoords, VM_Total);
            }

            // Reshape the vector
            for ( int j = 0; j < N.giveSize(); j++ ) {
                f(3 * j + 0) += N(j) * t(0) * dS;
                f(3 * j + 1) += N(j) * t(1) * dS;
                f(3 * j + 2) += N(j) * t(2) * dS;
            }
        }

        answer.resize(19);
        answer.zero();
        answer.assemble(f, this->edge_ordering [ iEdge - 1 ]);
    } else {
        OOFEM_ERROR("Tet1BubbleStokes :: computeEdgeLoadVector - Strange boundary condition type");
    }
}


void Tet1BubbleStokes :: computeBoundaryLoadVector(FloatArray &answer, Load *load, int iSurf, CharType type, ValueModeType mode, TimeStep *tStep)
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
        iRule.setUpIntegrationPoints(_Triangle, numberOfIPs, _Unknown);

        for ( int i = 0; i < iRule.getNumberOfIntegrationPoints(); i++ ) {
            gp = iRule.getIntegrationPoint(i);
            FloatArray *lcoords = gp->giveCoordinates();

            this->interp.surfaceEvalN(N, iSurf, * lcoords, FEIElementGeometryWrapper(this));
            double detJ = fabs(this->interp.surfaceGiveTransformationJacobian(iSurf, * lcoords, FEIElementGeometryWrapper(this)));
            double dA = gp->giveWeight() * detJ;

            if ( boundaryLoad->giveFormulationType() == BoundaryLoad :: BL_EntityFormulation ) { // Edge load in xi-eta system
                boundaryLoad->computeValueAt(t, tStep, * lcoords, VM_Total);
            } else { // Edge load in x-y system
                FloatArray gcoords;
                this->interp.edgeLocal2global(gcoords, iSurf, * lcoords, FEIElementGeometryWrapper(this));
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

        answer.resize(19);
        answer.zero();
        answer.assemble(f, this->surf_ordering [ iSurf - 1 ]);
    } else {
        OOFEM_ERROR("Tet1BubbleStokes :: computeBoundaryLoadVector - Strange boundary condition type");
    }
}

void Tet1BubbleStokes :: computeStiffnessMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    // Note: Working with the components; [K, G+Dp; G^T+Dv^T, C] . [v,p]
    FluidDynamicMaterial *mat = static_cast< FluidDynamicMaterial * >( this->giveMaterial() );
    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    GaussPoint *gp;
    FloatMatrix B(6, 15), EdB, K(15,15), G, Dp, DvT, C, Ed, dN;
    FloatArray *lcoords, dNv(15), N, Ep, Cd, tmpA, tmpB;
    double Cp;

    K.zero();
    G.zero();
    B.zero();

    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        // Compute Gauss point and determinant at current element
        gp = iRule->getIntegrationPoint(i);
        lcoords = gp->giveCoordinates();

        double detJ = fabs(this->interp.giveTransformationJacobian(* lcoords, FEIElementGeometryWrapper(this)));
        double dV = detJ * gp->giveWeight();

        this->interp.evaldNdx(dN, * lcoords, FEIElementGeometryWrapper(this));
        this->interp.evalN(N, * lcoords, FEIElementGeometryWrapper(this));
 
        for ( int j = 0, k = 0; j < 4; j++, k += 3 ) {
            dNv(k + 0) = B(0, k + 0) = B(4, k + 2) = B(5, k + 1) = dN(j, 0);
            dNv(k + 1) = B(1, k + 1) = B(3, k + 2) = B(5, k + 0) = dN(j, 1);
            dNv(k + 2) = B(2, k + 2) = B(3, k + 1) = B(4, k + 0) = dN(j, 2);
        }
        // Bubble contribution;
        dNv(12) = B(0,12) = B(4,14) = B(5,13) = dN(0,0)*N(1)*N(2)*N(3) + N(0)*dN(1,0)*N(2)*N(3) + N(0)*N(1)*dN(2,0)*N(3) + N(0)*N(1)*N(2)*dN(3,0);
        dNv(13) = B(1,13) = B(3,14) = B(5,12) = dN(0,1)*N(1)*N(2)*N(3) + N(0)*dN(1,1)*N(2)*N(3) + N(0)*N(1)*dN(2,1)*N(3) + N(0)*N(1)*N(2)*dN(3,1);
        dNv(14) = B(2,14) = B(3,13) = B(4,12) = dN(0,2)*N(1)*N(2)*N(3) + N(0)*dN(1,2)*N(2)*N(3) + N(0)*N(1)*dN(2,2)*N(3) + N(0)*N(1)*N(2)*dN(3,2);

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

    answer.resize(19, 19);
    answer.zero();
    answer.assemble(K, this->momentum_ordering);
    answer.assemble(GTDvT, this->conservation_ordering, this->momentum_ordering);
    answer.assemble(GDp, this->momentum_ordering, this->conservation_ordering);
    answer.assemble(C, this->conservation_ordering);
}

FEInterpolation *Tet1BubbleStokes :: giveInterpolation()
{
    return &interp;
}

FEInterpolation *Tet1BubbleStokes :: giveInterpolation(DofIDItem id)
{
    return &interp;
}

void Tet1BubbleStokes :: updateYourself(TimeStep *tStep)
{
    Element :: updateYourself(tStep);
}

// Some extension Interfaces to follow:

Interface *Tet1BubbleStokes :: giveInterface(InterfaceType it)
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

int Tet1BubbleStokes :: SpatialLocalizerI_containsPoint(const FloatArray &coords)
{
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}

void Tet1BubbleStokes :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType mode,
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

int Tet1BubbleStokes :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode, TimeStep *tStep, const FloatArray &gcoords, FloatArray &answer)
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

void Tet1BubbleStokes :: EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer)
{
    answer.setValues(4, V_u, V_v, V_w, P_f);
}

double Tet1BubbleStokes :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray center;
    FloatArray lcoords;
    lcoords.setValues(4, 0.333333, 0.333333, 0.333333, 0.333333);
    this->interp.local2global(center, lcoords, FEIElementGeometryWrapper(this));
    return center.distance(coords);
}

int Tet1BubbleStokes :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}

} // end namespace oofem
