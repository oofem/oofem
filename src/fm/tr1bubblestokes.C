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

#include "tr1bubblestokes.h"
#include "node.h"
#include "domain.h"
#include "equationid.h"
#include "gaussintegrationrule.h"
#include "gausspnt.h"
#include "bcgeomtype.h"
#include "load.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "fluiddynamicmaterial.h"
#include "fei2dtrlin.h"
#include "masterdof.h"

namespace oofem {
// Set up interpolation coordinates
FEI2dTrLin Tr1BubbleStokes :: interp(1, 2);
// Set up ordering vectors (for assembling)
IntArray Tr1BubbleStokes :: ordering(11);
IntArray Tr1BubbleStokes :: edge_ordering [ 3 ] = { IntArray(4), IntArray(4), IntArray(4) };
bool Tr1BubbleStokes :: __initialized = Tr1BubbleStokes :: initOrdering();

Tr1BubbleStokes :: Tr1BubbleStokes(int n, Domain *aDomain) : FMElement(n, aDomain)
{
    this->numberOfDofMans = 3;
    this->numberOfGaussPoints = 7;
    this->computeGaussPoints();
    
    this->bubble = new ElementDofManager(1, aDomain, this);
    this->bubble->appendDof(new MasterDof(1, this->bubble, V_u));
    this->bubble->appendDof(new MasterDof(2, this->bubble, V_v));
}

Tr1BubbleStokes :: ~Tr1BubbleStokes()
{
    delete this->bubble;
}

IRResultType Tr1BubbleStokes :: initializeFrom(InputRecord *ir)
{
    this->FMElement :: initializeFrom(ir);
    return IRRT_OK;
}

void Tr1BubbleStokes :: computeGaussPoints()
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 2 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Triangle, this->numberOfGaussPoints, _2dFlow);
    }
}

int Tr1BubbleStokes :: computeNumberOfDofs(EquationID ut)
{
    if ( ut == EID_MomentumBalance_ConservationEquation ) {
        return 11;
    } else if ( ut == EID_MomentumBalance ) {
        return 8;
    } else if ( ut == EID_ConservationEquation ) {
        return 3;
    }
    return 0;
}

void Tr1BubbleStokes :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    if ( ( ut == EID_MomentumBalance ) || ( ut == EID_AuxMomentumBalance ) ) {
        answer.setValues(2, V_u, V_v);
    } else if ( ut == EID_ConservationEquation ) {
        answer.setValues(1, P_f);
    } else if ( ut == EID_MomentumBalance_ConservationEquation ) {
        answer.setValues(3, V_u, V_v, P_f);
    } else {
        answer.resize(0);
    }
}

void Tr1BubbleStokes :: giveInternalDofManDofIDMask(int i, EquationID eid, IntArray &answer) const
{
    if ( eid == EID_MomentumBalance_ConservationEquation || eid == EID_MomentumBalance ) {
        answer.setValues(2, V_u, V_v);
    } else {
        answer.resize(0);
    }    
}

double Tr1BubbleStokes :: computeVolumeAround(GaussPoint *gp)
{
    double detJ = fabs(this->interp.giveTransformationJacobian(*gp->giveCoordinates(), FEIElementGeometryWrapper(this)));
    return detJ * gp->giveWeight();
}

void Tr1BubbleStokes :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
                                            TimeStep *tStep)
{
    // Compute characteristic vector for this element. I.e the load vector(s)
    if ( mtrx == ExternalForcesVector ) {
        this->computeLoadVector(answer, tStep);
    } else if ( mtrx == InternalForcesVector ) {
        this->computeInternalForcesVector(answer, tStep);
    } else {
        OOFEM_ERROR("giveCharacteristicVector: Unknown Type of characteristic mtrx.");
    }
}

void Tr1BubbleStokes :: giveCharacteristicMatrix(FloatMatrix &answer,
                                            CharType mtrx, TimeStep *tStep)
{
    if ( mtrx == StiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, tStep);
    } else {
        OOFEM_ERROR("giveCharacteristicMatrix: Unknown Type of characteristic mtrx.");
    }
}

void Tr1BubbleStokes :: computeInternalForcesVector(FloatArray &answer, TimeStep *tStep)
{
    IntegrationRule *iRule = integrationRulesArray [ 0 ];
    FluidDynamicMaterial *mat = ( FluidDynamicMaterial * ) this->domain->giveMaterial(this->material);
    FloatArray a_pressure, a_velocity, devStress, epsp, BTs, N, dNv(8);
    double r_vol, pressure;
    FloatMatrix dN, B(3, 8);
    B.zero();
    
    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, a_velocity);
    this->computeVectorOf(EID_ConservationEquation, VM_Total, tStep, a_pressure);

    FloatArray momentum(8), conservation(3);
    momentum.zero();
    conservation.zero();
    GaussPoint *gp;

    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        FloatArray *lcoords = gp->giveCoordinates();

        double detJ = fabs(this->interp.giveTransformationJacobian(* lcoords, FEIElementGeometryWrapper(this)));
        this->interp.evaldNdx(dN, * lcoords, FEIElementGeometryWrapper(this));
        this->interp.evalN(N, * lcoords, FEIElementGeometryWrapper(this));
        double dA = detJ * gp->giveWeight();

        for ( int j = 0, k = 0; j < 3; j++, k += 2 ) {
            dNv(k)     = B(0, k)     = B(2, k + 1) = dN(j, 0);
            dNv(k + 1) = B(1, k + 1) = B(2, k)     = dN(j, 1);
        }
        // Bubble contribution;
        dNv(6) = B(0,6) = B(2,7) = 27. * ( dN(0,0)*N(1)*N(2) + N(0)*dN(1,0)*N(2) + N(0)*N(1)*dN(2,0) ) ;
        dNv(7) = B(1,7) = B(2,6) = 27. * ( dN(0,1)*N(1)*N(2) + N(0)*dN(1,1)*N(2) + N(0)*N(1)*dN(2,1) ) ;

        pressure = N.dotProduct(a_pressure);
        epsp.beProductOf(B, a_velocity);

        mat->computeDeviatoricStressVector(devStress, r_vol, gp, epsp, pressure, tStep);
        BTs.beTProductOf(B, devStress);

        momentum.add(dA, BTs);
        momentum.add(-pressure*dA, dNv);
        conservation.add(r_vol*dA, N);
    }

    FloatArray temp(11);
    temp.zero();
    temp.addSubVector(momentum, 1);
    temp.addSubVector(conservation, 9);
    answer.resize(11);
    answer.zero();
    answer.assemble(temp, this->ordering);
}

void Tr1BubbleStokes :: computeLoadVector(FloatArray &answer, TimeStep *tStep)
{
    int load_number, load_id;
    Load *load;
    bcGeomType ltype;
    FloatArray vec;

    int nLoads = this->boundaryLoadArray.giveSize() / 2;
    answer.resize(11);
    answer.zero();
    for ( int i = 1; i <= nLoads; i++ ) {  // For each Neumann boundary condition
        load_number = this->boundaryLoadArray.at(2 * i - 1);
        load_id = this->boundaryLoadArray.at(2 * i);
        load = this->domain->giveLoad(load_number);
        ltype = load->giveBCGeoType();

        if ( ltype == EdgeLoadBGT ) {
            this->computeEdgeBCSubVectorAt(vec, load, load_id, tStep);
            answer.add(vec);
        }
    }

    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nLoads; i++ ) {
        load  = domain->giveLoad( bodyLoadArray.at(i) );
        ltype = load->giveBCGeoType();
        if ( ltype == BodyLoadBGT && load->giveBCValType() == ForceLoadBVT ) {
            this->computeBodyLoadVectorAt(vec, load, tStep);
            answer.add(vec);
        }
    }
}


void Tr1BubbleStokes :: computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep)
{
    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    GaussPoint *gp;
    FloatArray N, gVector, *lcoords, temparray(8);
    double dA, detJ, rho;

    load->computeComponentArrayAt(gVector, tStep, VM_Total);
    temparray.zero();
    if ( gVector.giveSize() ) {
        for ( int k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
            gp = iRule->getIntegrationPoint(k);
            lcoords = gp->giveCoordinates();

            rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, gp, tStep);
            detJ = fabs( this->interp.giveTransformationJacobian(* lcoords, FEIElementGeometryWrapper(this)) );
            dA = detJ * gp->giveWeight();

            this->interp.evalN(N, * lcoords, FEIElementGeometryWrapper(this));
            for ( int j = 0; j < 3; j++ ) {
                temparray(2 * j)     += N(j) * rho * gVector(0) * dA;
                temparray(2 * j + 1) += N(j) * rho * gVector(1) * dA;
            }
            temparray(6) += N(0)*N(1)*N(2) * rho * gVector(0) * dA;
            temparray(7) += N(0)*N(1)*N(2) * rho * gVector(1) * dA;
        }
    }

    answer.resize(11);
    answer.zero();
    answer.assemble( temparray, this->ordering );
}

void Tr1BubbleStokes :: computeEdgeBCSubVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep)
{
    answer.resize(11);
    answer.zero();

    if ( load->giveType() == TransmissionBC ) { // Neumann boundary conditions (traction)
        BoundaryLoad *boundaryLoad = ( BoundaryLoad * ) load;

        int numberOfEdgeIPs = ( int ) ceil( ( boundaryLoad->giveApproxOrder() + 2. ) / 2. );

        GaussIntegrationRule iRule(1, this, 1, 1);
        GaussPoint *gp;
        FloatArray N, t, f(4);
        IntArray edge_mapping;

        f.zero();
        iRule.setUpIntegrationPoints(_Line, numberOfEdgeIPs, _Unknown);

        for ( int i = 0; i < iRule.getNumberOfIntegrationPoints(); i++ ) {
            gp = iRule.getIntegrationPoint(i);
            FloatArray *lcoords = gp->giveCoordinates();

            this->interp.edgeEvalN(N, * lcoords, FEIElementGeometryWrapper(this));
            double detJ = fabs(this->interp.edgeGiveTransformationJacobian(iEdge, * lcoords, FEIElementGeometryWrapper(this)));
            double dS = gp->giveWeight() * detJ;

            if ( boundaryLoad->giveFormulationType() == BoundaryLoad :: BL_EntityFormulation ) { // Edge load in xi-eta system
                boundaryLoad->computeValueAt(t, tStep, * lcoords, VM_Total);
            } else   { // Edge load in x-y system
                FloatArray gcoords;
                this->interp.edgeLocal2global(gcoords, iEdge, * lcoords, FEIElementGeometryWrapper(this));
                boundaryLoad->computeValueAt(t, tStep, gcoords, VM_Total);
            }

            // Reshape the vector
            for ( int j = 0; j < 2; j++ ) {
                f(2 * j)     += N(j) * t(0) * dS;
                f(2 * j + 1) += N(j) * t(1) * dS;
            }
        }

        answer.assemble(f, this->edge_ordering [ iEdge - 1 ]);
    } else   {
        OOFEM_ERROR("Tr1BubbleStokes :: computeEdgeBCSubVectorAt - Strange boundary condition type");
    }
}

void Tr1BubbleStokes :: computeStiffnessMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    // Note: Working with the components; [K, G+Dp; G^T+Dv^T, C] . [v,p]
    FluidDynamicMaterial *mat = ( FluidDynamicMaterial * ) this->domain->giveMaterial(this->material);
    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    GaussPoint *gp;
    FloatMatrix B(3, 8), EdB, K(8,8), G, Dp, DvT, C, Ed, dN, GT;
    FloatArray *lcoords, dNv(8), N, Ep, Cd, tmpA, tmpB;
    double Cp;

    K.zero();
    G.zero();
    B.zero();

    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        // Compute Gauss point and determinant at current element
        gp = iRule->getIntegrationPoint(i);
        lcoords = gp->giveCoordinates();

        double detJ = fabs(this->interp.giveTransformationJacobian(* lcoords, FEIElementGeometryWrapper(this)));
        double dA = detJ * gp->giveWeight();

        this->interp.evaldNdx(dN, * lcoords, FEIElementGeometryWrapper(this));
        this->interp.evalN(N, * lcoords, FEIElementGeometryWrapper(this));
        for ( int j = 0, k = 0; j < 3; j++, k += 2 ) {
            dNv(k)     = B(0, k)     = B(2, k + 1) = dN(j, 0);
            dNv(k + 1) = B(1, k + 1) = B(2, k)     = dN(j, 1);
        }
        // Bubble contribution;
        dNv(6) = B(0,6) = B(2,7) = 27. * ( dN(0,0)*N(1)*N(2) + N(0)*dN(1,0)*N(2) + N(0)*N(1)*dN(2,0) );
        dNv(7) = B(1,7) = B(2,6) = 27. * ( dN(0,1)*N(1)*N(2) + N(0)*dN(1,1)*N(2) + N(0)*N(1)*dN(2,1) );

        // Computing the internal forces should have been done first.
        mat->giveDeviatoricStiffnessMatrix(Ed, TangentStiffness, gp, tStep); // dsigma_dev/deps_dev
        mat->giveDeviatoricPressureStiffness(Ep, TangentStiffness, gp, tStep); // dsigma_dev/dp
        mat->giveVolumetricDeviatoricStiffness(Cd, TangentStiffness, gp, tStep); // deps_vol/deps_dev
        mat->giveVolumetricPressureStiffness(Cp, TangentStiffness, gp, tStep); // deps_vol/dp

        EdB.beProductOf(Ed,B);
        K.plusProductSymmUpper(B, EdB, dA);
        G.plusDyadUnsym(dNv, N, -dA);
        C.plusDyadSymmUpper(N, N, Cp*dA);

        tmpA.beTProductOf(B, Ep);
        Dp.plusDyadUnsym(tmpA, N, dA);

        tmpB.beTProductOf(B, Cd);
        DvT.plusDyadUnsym(N, tmpB, dA);
    }

    K.symmetrized();
    C.symmetrized();

    GT.beTranspositionOf(G);
    
    FloatMatrix temp(11, 11);
    temp.zero();
    temp.addSubMatrix(K, 1, 1);
    temp.addSubMatrix(GT, 9, 1);
    temp.addSubMatrix(DvT, 9, 1);
    temp.addSubMatrix(G, 1, 9);
    temp.addSubMatrix(Dp, 1, 9);
    temp.addSubMatrix(C, 9, 9);

    answer.resize(11, 11);
    answer.zero();
    answer.assemble(temp, this->ordering);
}

FEInterpolation *Tr1BubbleStokes :: giveInterpolation()
{
    return &interp;
}

FEInterpolation *Tr1BubbleStokes :: giveInterpolation(DofIDItem id)
{
    return &interp;
}

void Tr1BubbleStokes :: updateYourself(TimeStep *tStep)
{
    Element :: updateYourself(tStep);
}

// Some extension Interfaces to follow:

Interface *Tr1BubbleStokes :: giveInterface(InterfaceType it)
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

int Tr1BubbleStokes :: SpatialLocalizerI_containsPoint(const FloatArray &coords)
{
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}

void Tr1BubbleStokes :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType mode,
        TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
{
    FloatArray n, n_lin;
    this->interp.evalN(n, lcoords, FEIElementGeometryWrapper(this));
    this->interp.evalN(n_lin, lcoords, FEIElementGeometryWrapper(this));
    answer.resize(3);
    answer.zero();
    for (int i = 1; i <= n.giveSize(); i++) {
        answer(0) += n.at(i)*this->giveNode(i)->giveDofWithID(V_u)->giveUnknown(EID_MomentumBalance_ConservationEquation, mode, tStep);
        answer(1) += n.at(i)*this->giveNode(i)->giveDofWithID(V_v)->giveUnknown(EID_MomentumBalance_ConservationEquation, mode, tStep);
    }
    for (int i = 1; i <= n_lin.giveSize(); i++) {
        answer(2) += n_lin.at(i)*this->giveNode(i)->giveDofWithID(P_f)->giveUnknown(EID_MomentumBalance_ConservationEquation, mode, tStep);
    }
}

int Tr1BubbleStokes :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode, TimeStep *tStep, const FloatArray &gcoords, FloatArray &answer)
{
    bool ok;
    FloatArray lcoords, n, n_lin;
    ok = this->computeLocalCoordinates(lcoords, gcoords);
    if (!ok) {
        answer.resize(0);
        return false;
    }
    this->EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(mode, tStep, lcoords, answer);
    return true;
}

void Tr1BubbleStokes :: EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer)
{
    answer.setValues(3, V_u, V_v, P_f);
}

double Tr1BubbleStokes :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray center;
    FloatArray lcoords;
    lcoords.setValues(3, 0.333333, 0.333333, 0.333333);
    this->interp.local2global(center, lcoords, FEIElementGeometryWrapper(this));
    return center.distance(coords);
}

int Tr1BubbleStokes :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}

} // end namespace oofem
