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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "tr21stokes.h"
#include "fmelement.h"
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
#include "fei2dtrquad.h"

namespace oofem {
// Set up interpolation coordinates
FEI2dTrLin Tr21Stokes :: interpolation_lin(1, 2);
FEI2dTrQuad Tr21Stokes :: interpolation_quad(1, 2);
// Set up ordering vectors (for assembling)
IntArray Tr21Stokes :: ordering(15);
IntArray Tr21Stokes :: edge_ordering [ 3 ] = {
    IntArray(6), IntArray(6), IntArray(6)
};
bool Tr21Stokes :: __initialized = Tr21Stokes :: initOrdering();

Tr21Stokes :: Tr21Stokes(int n, Domain *aDomain) : FMElement(n, aDomain)
{
    this->numberOfDofMans = 6;
    if (aDomain->giveEngngModel()->giveProblemScale() == macroScale) // Pick the lowest default value for multiscale computations.
        this->numberOfGaussPoints = 3;
    else 
        this->numberOfGaussPoints = 4;
    this->computeGaussPoints();
}

Tr21Stokes :: ~Tr21Stokes()
{}

IRResultType Tr21Stokes :: initializeFrom(InputRecord *ir)
{
    this->FMElement :: initializeFrom(ir);
    return IRRT_OK;
}

void Tr21Stokes :: computeGaussPoints()
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 2 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Triangle, this->numberOfGaussPoints, _2dFlow);
    }
}

int Tr21Stokes :: computeNumberOfDofs(EquationID ut)
{
    if ( ut == EID_MomentumBalance_ConservationEquation ) {
        return 15;
    } else if ( ut == EID_MomentumBalance ) {
        return 12;
    } else if ( ut == EID_ConservationEquation ) {
        return 3;
    } else {
        OOFEM_ERROR("computeNumberOfDofs: Unknown equation id encountered");
    }

    return 0;
}

void Tr21Stokes :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    // Returns the mask for node number inode of this element. The mask tells what quantities
    // are held by each node. Since this element holds velocities (both in x and y direction),
    // in six nodes and pressure in three nodes the answer depends on which node is requested.

    if ( ( inode == 1 ) || ( inode == 2 ) || ( inode == 3 ) ) {
        if ( ( ut == EID_MomentumBalance ) || ( ut == EID_AuxMomentumBalance ) ) {
            answer.setValues(2, V_u, V_v);
        } else if ( ut == EID_ConservationEquation ) {
            answer.setValues(1, P_f);
        } else if ( ut == EID_MomentumBalance_ConservationEquation ) {
            answer.setValues(3, V_u, V_v, P_f);
        } else {
            _error("giveDofManDofIDMask: Unknown equation id encountered");
        }
    } else if ( ( inode == 4 ) || ( inode == 5 ) || ( inode == 6 ) ) {
        if ( ( ut == EID_MomentumBalance ) || ( ut == EID_AuxMomentumBalance ) ) {
            answer.setValues(2, V_u, V_v);
        } else if ( ut == EID_ConservationEquation ) {
            answer.resize(0);
        } else if ( ut == EID_MomentumBalance_ConservationEquation ) {
            answer.setValues(2, V_u, V_v);
        } else {
            _error("giveDofManDofIDMask: Unknown equation id encountered");
        }
    }
}

double Tr21Stokes :: computeVolumeAround(GaussPoint *gp)
{
    double detJ = fabs(this->interpolation_quad.giveTransformationJacobian(*gp->giveCoordinates(), FEIElementGeometryWrapper(this)));
    return detJ * gp->giveWeight();
}

void Tr21Stokes :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
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

void Tr21Stokes :: giveCharacteristicMatrix(FloatMatrix &answer,
                                            CharType mtrx, TimeStep *tStep)
{
    // Compute characteristic matrix for this element. The only option is the stiffness matrix...
    if ( mtrx == StiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, tStep);
    } else {
        OOFEM_ERROR("giveCharacteristicMatrix: Unknown Type of characteristic mtrx.");
    }
}

void Tr21Stokes :: computeInternalForcesVector(FloatArray &answer, TimeStep *tStep)
{
    IntegrationRule *iRule = integrationRulesArray [ 0 ];
    FluidDynamicMaterial *mat = ( FluidDynamicMaterial * ) this->domain->giveMaterial(this->material);
    FloatArray a_pressure, a_velocity, devStress, epsp, BTs, Nh, dNv(12);
    double r_vol, pressure;
    FloatMatrix dN, B(3, 12);
    B.zero();
    
    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, a_velocity);
    this->computeVectorOf(EID_ConservationEquation, VM_Total, tStep, a_pressure);
    
    FloatArray momentum(12), conservation(3);
    momentum.zero();
    conservation.zero();
    GaussPoint *gp;

    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        FloatArray *lcoords = gp->giveCoordinates();

        double detJ = fabs(this->interpolation_quad.giveTransformationJacobian(* lcoords, FEIElementGeometryWrapper(this)));
        this->interpolation_quad.evaldNdx(dN, * lcoords, FEIElementGeometryWrapper(this));
        this->interpolation_lin.evalN(Nh, * lcoords, FEIElementGeometryWrapper(this));
        double dA = detJ * gp->giveWeight();

        for ( int j = 0, k = 0; j < 6; j++, k += 2 ) {
            dNv(k)     = B(0, k)     = B(2, k + 1) = dN(j, 0);
            dNv(k + 1) = B(1, k + 1) = B(2, k)     = dN(j, 1);
        }

        pressure = Nh.dotProduct(a_pressure);
        epsp.beProductOf(B, a_velocity);

        mat->computeDeviatoricStressVector(devStress, r_vol, gp, epsp, pressure, tStep);
        BTs.beTProductOf(B, devStress);

        momentum.add(dA, BTs);
        momentum.add(-pressure*dA, dNv);
        conservation.add(r_vol*dA, Nh);
    }

    FloatArray temp(15);
    temp.zero();
    temp.addSubVector(momentum, 1);
    temp.addSubVector(conservation, 13);

    answer.resize(15);
    answer.zero();
    answer.assemble(temp, this->ordering);
}

void Tr21Stokes :: computeLoadVector(FloatArray &answer, TimeStep *tStep)
{
    int i, load_number, load_id;
    Load *load;
    bcGeomType ltype;
    FloatArray vec;

    int nLoads = this->boundaryLoadArray.giveSize() / 2;
    answer.resize(15);
    answer.zero();
    for ( i = 1; i <= nLoads; i++ ) {  // For each Neumann boundary condition
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
    for ( i = 1; i <= nLoads; i++ ) {
        load  = domain->giveLoad( bodyLoadArray.at(i) );
        ltype = load->giveBCGeoType();
        if ( ltype == BodyLoadBGT && load->giveBCValType() == ForceLoadBVT ) {
            this->computeBodyLoadVectorAt(vec, load, tStep);
            answer.add(vec);
        }
    }
}


void Tr21Stokes :: computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep)
{
    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    GaussPoint *gp;
    FloatArray N, gVector, *lcoords, temparray(15);
    double dA, detJ, rho;

    load->computeComponentArrayAt(gVector, tStep, VM_Total);
    temparray.zero();
    if ( gVector.giveSize() ) {
        for ( int k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
            gp = iRule->getIntegrationPoint(k);
            lcoords = gp->giveCoordinates();

            rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, gp, tStep);
            detJ = fabs( this->interpolation_quad.giveTransformationJacobian(* lcoords, FEIElementGeometryWrapper(this)) );
            dA = detJ * gp->giveWeight();

            this->interpolation_quad.evalN(N, * lcoords, FEIElementGeometryWrapper(this));
            for ( int j = 0; j < 6; j++ ) {
                temparray(2 * j)     += N(j) * rho * gVector(0) * dA;
                temparray(2 * j + 1) += N(j) * rho * gVector(1) * dA;
            }
        }
    }

    answer.resize(15);
    answer.zero();
    answer.assemble( temparray, this->ordering );
}

void Tr21Stokes :: computeEdgeBCSubVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep)
{
    answer.resize(15);
    answer.zero();

    if ( load->giveType() == TransmissionBC ) { // Neumann boundary conditions (traction)
        BoundaryLoad *boundaryLoad = ( BoundaryLoad * ) load;

        int numberOfEdgeIPs = ( int ) ceil( ( boundaryLoad->giveApproxOrder() + 1. ) / 2. ) * 2;

        GaussIntegrationRule iRule(1, this, 1, 1);
        GaussPoint *gp;
        FloatArray N, t, f(6), f2(6);
        IntArray edge_mapping;

        f.zero();
        iRule.setUpIntegrationPoints(_Line, numberOfEdgeIPs, _Unknown);

        for ( int i = 0; i < iRule.getNumberOfIntegrationPoints(); i++ ) {
            gp = iRule.getIntegrationPoint(i);
            FloatArray *lcoords = gp->giveCoordinates();

            this->interpolation_quad.edgeEvalN(N, * lcoords, FEIElementGeometryWrapper(this));
            double detJ = fabs(this->interpolation_quad.edgeGiveTransformationJacobian(iEdge, * lcoords, FEIElementGeometryWrapper(this)));
            double dS = gp->giveWeight() * detJ;

            if ( boundaryLoad->giveFormulationType() == BoundaryLoad :: BL_EntityFormulation ) { // Edge load in xi-eta system
                boundaryLoad->computeValueAt(t, tStep, * lcoords, VM_Total);
            } else   { // Edge load in x-y system
                FloatArray gcoords;
                this->interpolation_quad.edgeLocal2global(gcoords, iEdge, * lcoords, FEIElementGeometryWrapper(this));
                boundaryLoad->computeValueAt(t, tStep, gcoords, VM_Total);
            }

            // Reshape the vector
            for ( int j = 0; j < 6; j++ ) {
                f(2 * j)     += N(j) * t(0) * dS;
                f(2 * j + 1) += N(j) * t(1) * dS;
            }
        }

        answer.assemble(f, this->edge_ordering [ iEdge - 1 ]);
    } else   {
        OOFEM_ERROR("Tr21Stokes :: computeEdgeBCSubVectorAt - Strange boundary condition type");
    }
}

void Tr21Stokes :: computeStiffnessMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    // Note: Working with the components; [K, G+Dp; G^T+Dv^T, C] . [v,p]
    FluidDynamicMaterial *mat = ( FluidDynamicMaterial * ) this->domain->giveMaterial(this->material);
    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    GaussPoint *gp;
    FloatMatrix B(3, 12), EdB, K(12,12), G, Dp, DvT, C, Ed, dN, GT;
    FloatArray *lcoords, dN_V(12), Nlin, Ep, Cd, tmpA, tmpB;
    double Cp;

    K.zero();
    G.zero();

    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        // Compute Gauss point and determinant at current element
        gp = iRule->getIntegrationPoint(i);
        lcoords = gp->giveCoordinates();

        double detJ = fabs(this->interpolation_quad.giveTransformationJacobian(* lcoords, FEIElementGeometryWrapper(this)));
        double dA = detJ * gp->giveWeight();

        this->interpolation_quad.evaldNdx(dN, * lcoords, FEIElementGeometryWrapper(this));
        this->interpolation_lin.evalN(Nlin, * lcoords, FEIElementGeometryWrapper(this));
        for ( int j = 0, k = 0; j < 6; j++, k += 2 ) {
            dN_V(k)     = B(0, k)     = B(2, k + 1) = dN(j, 0);
            dN_V(k + 1) = B(1, k + 1) = B(2, k)     = dN(j, 1);
        }

        // Computing the internal forces should have been done first.
        mat->giveDeviatoricStiffnessMatrix(Ed, TangentStiffness, gp, tStep); // dsigma_dev/deps_dev
        mat->giveDeviatoricPressureStiffness(Ep, TangentStiffness, gp, tStep); // dsigma_dev/dp
        mat->giveVolumetricDeviatoricStiffness(Cd, TangentStiffness, gp, tStep); // deps_vol/deps_dev
        mat->giveVolumetricPressureStiffness(Cp, TangentStiffness, gp, tStep); // deps_vol/dp

        EdB.beProductOf(Ed,B);
        K.plusProductSymmUpper(B, EdB, dA);
        G.plusDyadUnsym(dN_V, Nlin, -dA);
        C.plusDyadSymmUpper(Nlin, Nlin, Cp*dA);

        tmpA.beTProductOf(B, Ep);
        Dp.plusDyadUnsym(tmpA, Nlin, dA);

        tmpB.beTProductOf(B, Cd);
        DvT.plusDyadUnsym(Nlin, tmpB, dA);
    }

    K.symmetrized();
    C.symmetrized();

    GT.beTranspositionOf(G);

    FloatMatrix temp(15, 15);
    temp.zero();
    temp.addSubMatrix(K, 1, 1);
    temp.addSubMatrix(GT, 13, 1);
    temp.addSubMatrix(DvT, 13, 1);
    temp.addSubMatrix(G, 1, 13);
    temp.addSubMatrix(Dp, 1, 13);
    temp.addSubMatrix(C, 13, 13);

    answer.resize(15, 15);
    answer.zero();
    answer.assemble(temp, this->ordering);
}

FEInterpolation *Tr21Stokes :: giveInterpolation()
{
    return &interpolation_quad;
}

FEInterpolation *Tr21Stokes :: giveInterpolation(DofIDItem id)
{
    if (id == P_f) return &interpolation_lin; else return &interpolation_quad;
}

void Tr21Stokes :: updateYourself(TimeStep *tStep)
{
    Element :: updateYourself(tStep);
}

// Some extension Interfaces to follow:

Interface *Tr21Stokes :: giveInterface(InterfaceType it)
{
    switch (it) {
        case NodalAveragingRecoveryModelInterfaceType:
            return static_cast< NodalAveragingRecoveryModelInterface * >(this);
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

int Tr21Stokes :: SpatialLocalizerI_containsPoint(const FloatArray &coords)
{
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}

void Tr21Stokes :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType mode,
        TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
{
    FloatArray n, n_lin;
    this->interpolation_quad.evalN(n, lcoords, FEIElementGeometryWrapper(this));
    this->interpolation_lin.evalN(n_lin, lcoords, FEIElementGeometryWrapper(this));
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

int Tr21Stokes :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAt(ValueModeType mode, TimeStep *tStep, const FloatArray &gcoords, FloatArray &answer)
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

void Tr21Stokes :: EIPrimaryUnknownMI_givePrimaryUnknownVectorDofID(IntArray &answer)
{
    answer.setValues(3, V_u, V_v, P_f);
}

double Tr21Stokes :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray center;
    FloatArray lcoords;
    lcoords.setValues(3, 0.333333, 0.333333, 0.333333);
    interpolation_quad.local2global(center, lcoords, FEIElementGeometryWrapper(this));
    return center.distance(coords);
}

void Tr21Stokes :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}

int Tr21Stokes :: NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( type == IST_Pressure ) {
        return 1;
    }

    return 0;
}

void Tr21Stokes :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_Pressure ) {
        answer.resize(1);
        if ( node == 1 || node == 2 || node == 3 ) {
            answer.at(1) = this->giveNode(node)->giveDofWithID(P_f)->giveUnknown(EID_ConservationEquation, VM_Total, tStep);
        } else   {
            double a, b;
            if ( node == 4 ) {
                a = this->giveNode(1)->giveDofWithID(P_f)->giveUnknown(EID_ConservationEquation, VM_Total, tStep);
                b = this->giveNode(2)->giveDofWithID(P_f)->giveUnknown(EID_ConservationEquation, VM_Total, tStep);
            } else if ( node == 5 ) {
                a = this->giveNode(2)->giveDofWithID(P_f)->giveUnknown(EID_ConservationEquation, VM_Total, tStep);
                b = this->giveNode(3)->giveDofWithID(P_f)->giveUnknown(EID_ConservationEquation, VM_Total, tStep);
            } else /*if ( node == 6 )*/ {
                a = this->giveNode(3)->giveDofWithID(P_f)->giveUnknown(EID_ConservationEquation, VM_Total, tStep);
                b = this->giveNode(1)->giveDofWithID(P_f)->giveUnknown(EID_ConservationEquation, VM_Total, tStep);
            }

            answer.at(1) = ( a + b ) / 2;
        }
    } else   {
        answer.resize(0);
    }
}

int Tr21Stokes :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}
void Tr21Stokes :: giveGradP(FloatMatrix &answer, TimeStep * tStep )
{
    /*
    * Integrate gradient of P over element
    */

    answer.resize(2,1);
    answer.zero();

/*
    GaussIntegrationRule iRuleEdge(1, this, 1, 1);
    GaussPoint *gpEdge;
    FloatArray Normal, N, *lcoords, p;
    FloatMatrix temp, int_Np_edge;

    iRuleEdge.setUpIntegrationPoints(_Line, this->numberOfGaussPoints, _Unknown);

    this->computeVectorOf(EID_ConservationEquation, VM_Total, tStep, p);

    answer.resize(2,1);
    answer.zero();

    for (int i=1; i<=integrationEdges->giveSize(); i++ ) {

        int iEdge = integrationEdges->at(i);
        givePressureGradientBoundaryIntegral(int_Np_edge, iEdge);

        temp.resize(2,1);
        temp.zero();
        if (iEdge==1) {
            temp.at(1,1)=int_Np_edge.at(1,1)*p.at(1)+int_Np_edge.at(1,2)*p.at(2);
            temp.at(2,1)=int_Np_edge.at(2,1)*p.at(1)+int_Np_edge.at(2,2)*p.at(2);
        } else if (iEdge==2) {
            temp.at(1,1)=int_Np_edge.at(1,1)*p.at(2)+int_Np_edge.at(1,2)*p.at(3);
            temp.at(2,1)=int_Np_edge.at(2,1)*p.at(2)+int_Np_edge.at(2,2)*p.at(3);
        } else if (iEdge==3) {
            temp.at(1,1)=int_Np_edge.at(1,1)*p.at(3)+int_Np_edge.at(1,2)*p.at(1);
            temp.at(2,1)=int_Np_edge.at(2,1)*p.at(3)+int_Np_edge.at(2,2)*p.at(1);
        } else {
            printf ("iEdge != {1,2,3}\n");
        }
        answer.add(temp);
    }
*/
}

void Tr21Stokes :: giveIntegratedVelocity(FloatMatrix &answer, TimeStep *tStep )
{
    /*
    * Integrate velocity over element
    */

    IntegrationRule *iRule = integrationRulesArray [ 0 ];
    FloatMatrix v, v_gamma, ThisAnswer, boundaryV, Nmatrix;
    double detJ;
    FloatArray *lcoords, N;
    int i, j, k=0;
    Dof *d;
    GaussPoint *gp;

    v.resize(12,1);
    v.zero();
    boundaryV.resize(2,1);


    for (i=1; i<=this->giveNumberOfDofManagers(); i++) {
        for (j=1; j<=this->giveDofManager(i)->giveNumberOfDofs(); j++) {
            d = this->giveDofManager(i)->giveDof(j);
            if ((d->giveDofID()==V_u) || (d->giveDofID()==V_v)) {
                k=k+1;
                v.at(k,1)=d->giveUnknown(EID_ConservationEquation, VM_Total, tStep);
            /*} else	if (d->giveDofID()==A_x) {
                boundaryV.at(1,1)=d->giveUnknown(EID_ConservationEquation, VM_Total, tStep);
            } else	if (d->giveDofID()==A_y) {
                boundaryV.at(2,1)=d->giveUnknown(EID_ConservationEquation, VM_Total, tStep);*/
            }
        }
    }

    answer.resize(2,1);
    answer.zero();

    Nmatrix.resize(2,12);

    for (i=0; i<iRule->getNumberOfIntegrationPoints(); i++) {

        gp = iRule->getIntegrationPoint(i);

        lcoords = gp->giveCoordinates();

        this->interpolation_quad.evalN(N, *lcoords, FEIElementGeometryWrapper(this));
        detJ = this->interpolation_quad.giveTransformationJacobian(*lcoords, FEIElementGeometryWrapper(this));

        N.times(detJ*gp->giveWeight());

        for (j=1; j<=6;j++) {
            Nmatrix.at(1,j*2-1)=N.at(j);
            Nmatrix.at(2,j*2)=N.at(j);
        }

        ThisAnswer.beProductOf(Nmatrix,v);
        answer.add(ThisAnswer);

    }

}

void Tr21Stokes :: giveElementFMatrix(FloatMatrix &answer)
{

    IntegrationRule *iRule = integrationRulesArray [ 0 ];
    GaussPoint *gp;
    double detJ;
    FloatArray N, N2, *lcoords;
    IntArray col;
    FloatMatrix temp;

    N2.resize(6);	N2.zero();
    col.resize(2);  col.at(1)=1;  	col.at(2)=2;

    temp.resize(15,2);
    temp.zero();

    for (int i=0; i<iRule->getNumberOfIntegrationPoints(); i++) {
        gp = iRule->getIntegrationPoint(i);
        lcoords = gp->giveCoordinates();

        this->interpolation_quad.evalN(N, *lcoords, FEIElementGeometryWrapper(this));
        detJ = this->interpolation_quad.giveTransformationJacobian(*lcoords, FEIElementGeometryWrapper(this));
        N.times(gp->giveWeight()*detJ);
        //N.printYourself();
        N2.add(N);
    }

    for (int i=1; i<=6; i++) {
        temp.at(i*2-1, 1)=N2.at(i);
        temp.at(i*2, 2)=N2.at(i);
    }

    answer.resize(17,2);
    answer.zero();
    answer.assemble(temp, this->ordering, col);

}
} // end namespace oofem
