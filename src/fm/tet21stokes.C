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

#include "tet21stokes.h"
#include "fmelement.h"
#include "node.h"
#include "domain.h"
#include "equationid.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "bcgeomtype.h"
#include "generalboundarycondition.h"
#include "load.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "fluiddynamicmaterial.h"
#include "fei3dtetlin.h"
#include "fei3dtetquad.h"
#include "crosssection.h"
#include "classfactory.h"

namespace oofem {

REGISTER_Element( Tet21Stokes );

// Set up interpolation coordinates
FEI3dTetLin Tet21Stokes :: interpolation_lin;
FEI3dTetQuad Tet21Stokes :: interpolation_quad;
// Set up ordering vectors (for assembling)
IntArray Tet21Stokes :: momentum_ordering(30);
IntArray Tet21Stokes :: conservation_ordering(4);
IntArray Tet21Stokes :: surf_ordering [ 4 ] = {
    IntArray(18), IntArray(18), IntArray(18), IntArray(18)
};
bool Tet21Stokes :: __initialized = Tet21Stokes :: initOrdering();

Tet21Stokes :: Tet21Stokes(int n, Domain *aDomain) : FMElement(n, aDomain)
{
    this->numberOfDofMans = 10;
    this->numberOfGaussPoints = 4;
}

Tet21Stokes :: ~Tet21Stokes()
{}

IRResultType Tet21Stokes :: initializeFrom(InputRecord *ir)
{
    this->FMElement :: initializeFrom(ir);
    return IRRT_OK;
}

void Tet21Stokes :: computeGaussPoints()
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 2 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        this->giveCrossSection()->setupIntegrationPoints( *integrationRulesArray[0], this->numberOfGaussPoints, this );
    }
}

int Tet21Stokes :: computeNumberOfDofs(EquationID ut)
{
    if ( ut == EID_MomentumBalance_ConservationEquation ) {
        return 34;
    } else if ( ut == EID_MomentumBalance ) {
        return 30;
    } else if ( ut == EID_ConservationEquation ) {
        return 4;
    } else {
        OOFEM_ERROR("Tet21Stokes :: computeNumberOfDofs: Unknown equation id encountered");
    }

    return 0;
}

void Tet21Stokes :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    // Returns the mask for node number inode of this element. The mask tells what quantities
    // are held by each node. Since this element holds velocities (both in x and y direction),
    // in six nodes and pressure in three nodes the answer depends on which node is requested.

    if ( inode <= 4 ) {
        if ( ut == EID_MomentumBalance ) {
            answer.setValues(3, V_u, V_v, V_w);
        } else if ( ut == EID_ConservationEquation ) {
            answer.setValues(1, P_f);
        } else if ( ut == EID_MomentumBalance_ConservationEquation ) {
            answer.setValues(4, V_u, V_v, V_w, P_f);
        } else {
            OOFEM_ERROR("Tet21Stokes :: giveDofManDofIDMask: Unknown equation id encountered");
        }
    } else {
        if ( ut == EID_MomentumBalance || ut == EID_MomentumBalance_ConservationEquation ) {
            answer.setValues(3, V_u, V_v, V_w);
        } else if ( ut == EID_ConservationEquation ) {
            answer.resize(0);
        } else {
            OOFEM_ERROR("Tet21Stokes :: giveDofManDofIDMask: Unknown equation id encountered");
        }
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
        OOFEM_ERROR("Tet21Stokes :: giveCharacteristicVector: Unknown Type of characteristic mtrx.");
    }
}

void Tet21Stokes :: giveCharacteristicMatrix(FloatMatrix &answer,
                                            CharType mtrx, TimeStep *tStep)
{
    // Compute characteristic matrix for this element. The only option is the stiffness matrix...
    if ( mtrx == StiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, tStep);
    } else {
        OOFEM_ERROR("Tet21Stokes :: giveCharacteristicMatrix: Unknown Type of characteristic mtrx.");
    }
}

void Tet21Stokes :: computeInternalForcesVector(FloatArray &answer, TimeStep *tStep)
{
    IntegrationRule *iRule = integrationRulesArray [ 0 ];
    FluidDynamicMaterial *mat = static_cast<FluidDynamicMaterial * >( this->giveMaterial() );
    FloatArray a_pressure, a_velocity, devStress, epsp, Nh, dN_V(30);
    FloatMatrix dN, B(6, 30);
    double r_vol, pressure;
    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, a_velocity);
    this->computeVectorOf(EID_ConservationEquation, VM_Total, tStep, a_pressure);
    FloatArray momentum, conservation;

    B.zero();
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(i);
        FloatArray *lcoords = gp->giveCoordinates();

        double detJ = fabs( this->interpolation_quad.evaldNdx(dN, * lcoords, FEIElementGeometryWrapper(this)) );
        this->interpolation_lin.evalN(Nh, * lcoords, FEIElementGeometryWrapper(this));
        double dV = detJ * gp->giveWeight();

        for ( int j = 0, k = 0; j < dN.giveNumberOfColumns(); j++, k+=3 ) {
            dN_V(k + 0) = B(0, k + 0) = B(3, k + 1) = B(4, k + 2) = dN(0, j);
            dN_V(k + 1) = B(1, k + 1) = B(3, k + 0) = B(5, k + 2) = dN(1, j);
            dN_V(k + 2) = B(2, k + 2) = B(4, k + 0) = B(5, k + 1) = dN(2, j);
        }
        epsp.beProductOf(B, a_velocity);
        pressure = Nh.dotProduct(a_pressure);
        mat->computeDeviatoricStressVector(devStress, r_vol, gp, epsp, pressure, tStep);

        momentum.plusProduct(B, devStress, dV);
        momentum.add(-pressure*dV, dN_V);
        conservation.add(r_vol*dV, Nh);
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
    bcGeomType ltype;
    FloatArray vec;

    int nLoads = this->boundaryLoadArray.giveSize() / 2;
    answer.resize(0);

    for ( int i = 1; i <= nLoads; i++ ) {  // For each Neumann boundary condition
        load_number = this->boundaryLoadArray.at(2 * i - 1);
        load_id = this->boundaryLoadArray.at(2 * i);
        load = this->domain->giveLoad(load_number);
        ltype = load->giveBCGeoType();

        if ( ltype == SurfaceLoadBGT ) {
            this->computeBoundaryLoadVector(vec, load, load_id, ExternalForcesVector, VM_Total, tStep);
            answer.add(vec);
        }
    }

    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nLoads; i++ ) {
        load  = domain->giveLoad( bodyLoadArray.at(i) );
        ltype = load->giveBCGeoType();
        if ( ltype == BodyLoadBGT && load->giveBCValType() == ForceLoadBVT ) {
            this->computeLoadVector(vec, load, ExternalForcesVector, VM_Total, tStep);
            answer.add(vec);
        }
    }
}

void Tet21Stokes :: computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.resize(0);
        return;
    }

    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    FloatArray N, gVector, temparray(30);

    load->computeComponentArrayAt(gVector, tStep, VM_Total);
    temparray.zero();
    if ( gVector.giveSize() ) {
        for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(k);
            FloatArray *lcoords = gp->giveCoordinates();

            double rho = this->giveMaterial()->give('d', gp);
            double detJ = fabs( this->interpolation_quad.giveTransformationJacobian(* lcoords, FEIElementGeometryWrapper(this)) );
            double dA = detJ * gp->giveWeight();

            this->interpolation_quad.evalN(N, * lcoords, FEIElementGeometryWrapper(this));
            for ( int j = 0; j < N.giveSize(); j++ ) {
                temparray(3 * j + 0) += N(j) * rho * gVector(0) * dA;
                temparray(3 * j + 1) += N(j) * rho * gVector(1) * dA;
                temparray(3 * j + 2) += N(j) * rho * gVector(2) * dA;
            }
        }
    }

    answer.resize(34);
    answer.zero();
    answer.assemble( temparray, this->momentum_ordering );

}

void Tet21Stokes :: computeBoundaryLoadVector(FloatArray &answer, Load *load, int iSurf, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.resize(0);
        return;
    }

    answer.resize(34);
    answer.zero();

    if ( load->giveType() == TransmissionBC ) { // Neumann boundary conditions (traction)
        BoundaryLoad *boundaryLoad = static_cast< BoundaryLoad * >( load );

        int numberOfSurfaceIPs = ( int ) ceil( ( boundaryLoad->giveApproxOrder() + 1. ) / 2. ) * 2; ///@todo Check this.

        GaussIntegrationRule iRule(1, this, 1, 1);
        GaussPoint *gp;
        FloatArray N, t, f(18);

        f.zero();
        iRule.SetUpPointsOnTriangle(numberOfSurfaceIPs, _Unknown);

        for ( int i = 0; i < iRule.giveNumberOfIntegrationPoints(); i++ ) {
            gp = iRule.getIntegrationPoint(i);
            FloatArray *lcoords = gp->giveCoordinates();

            this->interpolation_quad.surfaceEvalN(N, iSurf, * lcoords, FEIElementGeometryWrapper(this));
            double dA = gp->giveWeight() * this->interpolation_quad.surfaceGiveTransformationJacobian(iSurf, * lcoords, FEIElementGeometryWrapper(this));

            if ( boundaryLoad->giveFormulationType() == BoundaryLoad :: BL_EntityFormulation ) { // load in xi-eta system
                boundaryLoad->computeValueAt(t, tStep, * lcoords, VM_Total);
            } else { // Edge load in x-y system
                FloatArray gcoords;
                this->interpolation_quad.surfaceLocal2global(gcoords, iSurf, * lcoords, FEIElementGeometryWrapper(this));
                boundaryLoad->computeValueAt(t, tStep, gcoords, VM_Total);
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
        OOFEM_ERROR("Tet21Stokes :: Strange boundary condition type");
    }
}

void Tet21Stokes :: computeStiffnessMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    FluidDynamicMaterial *mat = static_cast< FluidDynamicMaterial * >( this->giveMaterial() );
    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    FloatMatrix B(6, 30), EdB, K, G, Dp, DvT, C, Ed, dN;
    FloatArray dN_V(30), Nlin, Ep, Cd, tmpA, tmpB;
    double Cp;

    B.zero();

    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        // Compute Gauss point and determinant at current element
        GaussPoint *gp = iRule->getIntegrationPoint(i);
        FloatArray *lcoords = gp->giveCoordinates();

        double detJ = fabs( this->interpolation_quad.evaldNdx(dN, * lcoords, FEIElementGeometryWrapper(this)) );
        double dV = detJ * gp->giveWeight();
        this->interpolation_lin.evalN(Nlin, * lcoords, FEIElementGeometryWrapper(this));

//        dN.printYourself();

        for ( int j = 0, k = 0; j < dN.giveNumberOfColumns(); j++, k+=3 ) {
            ///@todo This column-row order of dN is inconsistent with every other interpolator class:
            dN_V(k + 0) = B(0, k + 0) = B(3, k + 1) = B(4, k + 2) = dN(0, j);
            dN_V(k + 1) = B(1, k + 1) = B(3, k + 0) = B(5, k + 2) = dN(1, j);
            dN_V(k + 2) = B(2, k + 2) = B(4, k + 0) = B(5, k + 1) = dN(2, j);
        }

        // Computing the internal forces should have been done first.
        mat->giveDeviatoricStiffnessMatrix(Ed, TangentStiffness, gp, tStep); // dsigma_dev/deps_dev
        mat->giveDeviatoricPressureStiffness(Ep, TangentStiffness, gp, tStep); // dsigma_dev/dp
        mat->giveVolumetricDeviatoricStiffness(Cd, TangentStiffness, gp, tStep); // deps_vol/deps_dev
        mat->giveVolumetricPressureStiffness(Cp, TangentStiffness, gp, tStep); // deps_vol/dp

//        Ed.printYourself();
//        B.printYourself();
        EdB.beProductOf(Ed,B);
        K.plusProductSymmUpper(B, EdB, dV);
        G.plusDyadUnsym(dN_V, Nlin, -dV);
        C.plusDyadSymmUpper(Nlin, Cp*dV);

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

    answer.resize(34, 34);
    answer.zero();
    answer.assemble(K, this->momentum_ordering);
    answer.assemble(GDp, this->momentum_ordering, this->conservation_ordering);
    answer.assemble(GTDvT, this->conservation_ordering, this->momentum_ordering);
    answer.assemble(C, this->conservation_ordering);
//    K.printYourself();
//    GDp.printYourself();
//    GTDvT.printYourself();
//    C.printYourself();
}

FEInterpolation *Tet21Stokes :: giveInterpolation() const
{
    return &interpolation_quad;
}

FEInterpolation *Tet21Stokes :: giveInterpolation(DofIDItem id) const
{
    if (id == P_f) return &interpolation_lin; else return &interpolation_quad;
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

int Tet21Stokes :: SpatialLocalizerI_containsPoint(const FloatArray &coords)
{
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}

void Tet21Stokes :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType mode,
        TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
{
    FloatArray n, n_lin;
    this->interpolation_quad.evalN(n, lcoords, FEIElementGeometryWrapper(this));
    this->interpolation_lin.evalN(n_lin, lcoords, FEIElementGeometryWrapper(this));
    answer.resize(4);
    answer.zero();
    for ( int i = 1; i <= n.giveSize(); i++ ) {
        answer(0) += n.at(i)*this->giveNode(i)->giveDofWithID(V_u)->giveUnknown(mode, tStep);
        answer(1) += n.at(i)*this->giveNode(i)->giveDofWithID(V_v)->giveUnknown(mode, tStep);
        answer(2) += n.at(i)*this->giveNode(i)->giveDofWithID(V_w)->giveUnknown(mode, tStep);
    }
    for ( int i = 1; i <= n_lin.giveSize(); i++ ) {
        answer(3) += n_lin.at(i)*this->giveNode(i)->giveDofWithID(P_f)->giveUnknown(mode, tStep);
    }
}

int Tet21Stokes :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode, TimeStep *tStep, const FloatArray &gcoords, FloatArray &answer)
{
    bool ok;
    FloatArray lcoords, n, n_lin;
    ok = this->computeLocalCoordinates(lcoords, gcoords);
    if ( !ok ) {
        answer.resize(0);
        return false;
    }
    this->EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(mode, tStep, lcoords, answer);
    return true;
}

void Tet21Stokes :: EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer)
{
    answer.setValues(4, V_u, V_v, V_w, P_f);
}

double Tet21Stokes :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray center;
    FloatArray lcoords;
    lcoords.setValues(3, 0.3333333, 0.3333333, 0.3333333);
    this->computeGlobalCoordinates(center, lcoords);
    return center.distance(coords);
}

void Tet21Stokes :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}

int Tet21Stokes :: NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( type == IST_Pressure ) {
        return 1;
    }

    return 0;
}

void Tet21Stokes :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_Pressure ) {
        answer.resize(1);
        if ( node <= 4 ) {
            answer.at(1) = this->giveNode(node)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep);
        } else {
            IntArray eNodes;
            this->interpolation_quad.computeLocalEdgeMapping(eNodes, node - 4);
            answer.at(1) = 0.5 * (
                this->giveNode(eNodes.at(1))->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) +
                this->giveNode(eNodes.at(2))->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) );
        }
    } else {
        answer.resize(0);
    }
}

void Tet21Stokes :: giveIntegratedVelocity(FloatMatrix &answer, TimeStep *tStep )
{
    /*
    * Integrate velocity over element
    */

    IntegrationRule *iRule = integrationRulesArray [ 0 ];
    FloatMatrix v, v_gamma, ThisAnswer, Nmatrix;
    double detJ;
    FloatArray *lcoords, N;
    int i, j, k=0;
    Dof *d;
    GaussPoint *gp;

    v.resize(30,1);
    v.zero();

    FloatArray test;

    for (i=1; i<=this->giveNumberOfDofManagers(); i++) {

        for (j=1; j<=this->giveDofManager(i)->giveNumberOfDofs(); j++) {
            d = this->giveDofManager(i)->giveDof(j);

            if ((d->giveDofID()==V_u) || (d->giveDofID()==V_v) || (d->giveDofID()==V_w) ) {
                k=k+1;
                v.at(k,1)=d->giveUnknown(VM_Total, tStep);
            }
        }
//    	v.printYourself();

    }

    answer.resize(3,1);
    answer.zero();

    Nmatrix.resize(3,30);

    for (i=0; i<iRule->giveNumberOfIntegrationPoints(); i++) {

        gp = iRule->getIntegrationPoint(i);

        lcoords = gp->giveCoordinates();

        this->interpolation_quad.evalN(N, *lcoords, FEIElementGeometryWrapper(this));
        detJ = this->interpolation_quad.giveTransformationJacobian(*lcoords, FEIElementGeometryWrapper(this));

        N.times(detJ*gp->giveWeight());

        for (j=1; j<=6;j++) {
            Nmatrix.at(1,j*3-2)=N.at(j);
            Nmatrix.at(2,j*3-1)=N.at(j);
            Nmatrix.at(3,j*3)=N.at(j);
        }

//        Nmatrix.printYourself();

        ThisAnswer.beProductOf(Nmatrix,v);
        answer.add(ThisAnswer);

    }

}

} // end namespace oofem
