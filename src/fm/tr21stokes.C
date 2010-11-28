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
 *               Copyright (C) 1993 - 2010   Borek Patzak
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
#include "domain.h"
#include "equationid.h"
#include "gaussintegrationrule.h"
#include "gausspnt.h"
#include "bcgeomtype.h"
#include "generalbc.h"
#include "load.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "fluiddynamicmaterial.h"

#ifndef __MAKEDEPEND
 #include <math.h>
#endif

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
    numberOfDofMans = 6;
    numberOfGaussPoints = 4;
}

Tr21Stokes :: ~Tr21Stokes()
{}

IRResultType Tr21Stokes :: initializeFrom(InputRecord *ir)
{
    this->FMElement :: initializeFrom(ir);
    this->computeGaussPoints();
    return IRRT_OK;
}

void Tr21Stokes :: computeGaussPoints()
{
    // Set up gausspoints for element

    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 2 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Triangle, numberOfGaussPoints, _2dFlow);
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
            answer.resize(2);
            answer.at(1) = V_u;
            answer.at(2) = V_v;
        } else if ( ut == EID_ConservationEquation ) {
            answer.resize(1);
            answer.at(1) = P_f;
        } else if ( ut == EID_MomentumBalance_ConservationEquation ) {
            answer.resize(3);
            answer.at(1) = V_u;
            answer.at(2) = V_v;
            answer.at(3) = P_f;
        } else {
            _error("giveDofManDofIDMask: Unknown equation id encountered");
        }
    } else if ( ( inode == 4 ) || ( inode == 5 ) || ( inode == 6 ) ) {
        if ( ( ut == EID_MomentumBalance ) || ( ut == EID_AuxMomentumBalance ) ) {
            answer.resize(2);
            answer.at(1) = V_u;
            answer.at(2) = V_v;
        } else if ( ut == EID_ConservationEquation ) {
            answer.resize(0);
        } else if ( ut == EID_MomentumBalance_ConservationEquation ) {
            answer.resize(2);
            answer.at(1) = V_u;
            answer.at(2) = V_v;
        } else {
            _error("giveDofManDofIDMask: Unknown equation id encountered");
        }
    }
}

void Tr21Stokes :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
                                            TimeStep *tStep)
{
    // Compute characteristic vector for this element. I.e the load vector(s)
    if ( mtrx == LoadVector ) {
        this->computeLoadVector(answer, tStep);
    } else if ( mtrx == NodalInternalForcesVector ) {
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
    FloatMatrix dN, B(3, 12);

    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, a_velocity);
    this->computeVectorOf(EID_ConservationEquation, VM_Total, tStep, a_pressure);

    FloatArray momentum(12), conservation(3);
    momentum.zero();
    conservation.zero();
    GaussPoint *gp;
    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        FloatArray *lcoords = gp->giveCoordinates();

        //double detJ = this->interpolation_quad.giveTransformationJacobian(this->domain, dofManArray, *lcoords, 0.0);
        double detJ = this->interpolation_quad.giveTransformationJacobian(* lcoords, FEIElementGeometryWrapper(this), 0.0);
        //this->interpolation_quad.evaldNdx(dN, this->domain, dofManArray, *lcoords, 0.0);
        this->interpolation_quad.evaldNdx(dN, * lcoords, FEIElementGeometryWrapper(this), 0.0);
        //this->interpolation_lin.evalN(Nh, *lcoords, 0.0);
        this->interpolation_lin.evalN(Nh, * lcoords, FEIElementGeometryWrapper(this), 0.0);        // Why FEI?
        double dA = detJ * gp->giveWeight();

        for ( int j = 0, k = 0; j < 6; j++, k += 2 ) {
            dNv(k)   = B(0, k)   = B(2, k + 1) = dN(j, 0);
            dNv(k + 1) = B(1, k + 1) = B(2, k)   = dN(j, 1);
        }

        //FluidDynamicMaterialStatus *ms = static_cast<FluidDynamicMaterialStatus*>(mat->giveStatus(gp));
        //FloatArray devStress = ms->giveDeviatoricStressVector();
        epsp.beProductOf(B, a_velocity);
        mat->computeDeviatoricStressVector(devStress, gp, epsp, tStep);
        double pressure = dotProduct(Nh, a_pressure, 3);
        double w = dotProduct(dNv, a_velocity, 12);
        BTs.beTProductOf(B, devStress);

        momentum += ( BTs - dNv * pressure ) * dA;
        conservation -= Nh * ( w * dA );
    }

    FloatArray temp(15);
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

    for ( i = 1; i <= nLoads; i++ ) {  // For each Neumann boundary condition ....
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
    // FIXME: Untested!

    answer.resize(15, 1);
    answer.zero();

    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    GaussPoint *gp;
    FloatArray N, b, gVector, *lcoords;
    double dA, detJ, rho;

    load->computeComponentArrayAt(gVector, tStep, VM_Total);
    if ( gVector.giveSize() ) {
        for ( int k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
            gp = iRule->getIntegrationPoint(k);
            lcoords = gp->giveCoordinates();

            rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, gp, tStep);
            //detJ = this->interpolation_quad.giveTransformationJacobian(this->domain, this->dofManArray, *lcoords, 0.0);
            detJ = this->interpolation_quad.giveTransformationJacobian(* lcoords, FEIElementGeometryWrapper(this), 0.0);
            dA = detJ * gp->giveWeight();

            //this->interpolation_quad.evalN(N, *lcoords, 0.0);
            this->interpolation_quad.evalN(N, * lcoords, FEIElementGeometryWrapper(this), 0.0);    // Why geometry?
            for ( int j = 0; j < 3; j++ ) {
                answer(2 * j)   += N(j) * rho * b(0) * dA;
                answer(2 * j + 1) += N(j) * rho * b(1) * dA;
            }
        }
    }

    OOFEM_ERROR("computeBodyLoadVectorAt: Not tested yet.");
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

            //this->interpolation_quad.edgeEvalN (N, *(gp->giveCoordinates()), 0.0);
            this->interpolation_quad.edgeEvalN(N, * lcoords, FEIElementGeometryWrapper(this), 0.0);
            double dS = gp->giveWeight() * this->interpolation_quad.edgeGiveTransformationJacobian(iEdge, * lcoords, FEIElementGeometryWrapper(this), 0.0);

            if ( boundaryLoad->giveFormulationType() == BoundaryLoad :: BL_EntityFormulation ) {         // Edge load in xi-eta system
                boundaryLoad->computeValueAt(t, tStep, * lcoords, VM_Total);
            } else   { // Edge load in x-y system
                FloatArray gcoords;
                this->interpolation_quad.edgeLocal2global(gcoords, iEdge, * lcoords, FEIElementGeometryWrapper(this), 0.0);
                boundaryLoad->computeValueAt(t, tStep, gcoords, VM_Total);
            }

            // Reshape the vector
            for ( int j = 0; j < 3; j++ ) {
                f(2 * j)   += N(j) * t(0) * dS;
                f(2 * j + 1) += N(j) * t(1) * dS;
            }
        }

        answer.assemble(f, this->edge_ordering [ iEdge - 1 ]);
    } else   {
        OOFEM_ERROR("Strange boundary condition type");
    }
}

void Tr21Stokes :: computeStiffnessMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    FluidDynamicMaterial *mat = ( FluidDynamicMaterial * ) this->domain->giveMaterial(this->material);
    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    GaussPoint *gp;
    FloatMatrix B(3, 12), BTCB(12, 12), G(12, 3), C, dN, GT, Ia, Ia2, Ib;
    FloatArray *lcoords, dN_V(12), Nlin;

    BTCB.zero();
    G.zero();

    FloatArray a, stress;
    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, a);


    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        // Compute Gauss point and determinant at current element
        gp = iRule->getIntegrationPoint(i);
        lcoords = gp->giveCoordinates();

        double detJ = this->interpolation_quad.giveTransformationJacobian(* lcoords, FEIElementGeometryWrapper(this), 0.0);
        double dA = detJ * gp->giveWeight();

        this->interpolation_quad.evaldNdx(dN, * lcoords, FEIElementGeometryWrapper(this), 0.0);
        this->interpolation_lin.evalN(Nlin, * lcoords, FEIElementGeometryWrapper(this), 0.0);
        for ( int j = 0, k = 0; j < 6; j++, k += 2 ) {
            dN_V(k)   = B(0, k)   = B(2, k + 1) = dN(j, 0);
            dN_V(k + 1) = B(1, k + 1) = B(2, k)   = dN(j, 1);
        }

        // Just so that the stress is stored in the material status.
        FloatArray eps;
        eps.beProductOf(B, a);
        mat->computeDeviatoricStressVector(stress, gp, eps, tStep);
        mat->giveDeviatoricStiffnessMatrix(C, TangentStiffness, gp, tStep);
        // Build B^T*C*B in elemental stiffness matrix
        Ia.beTProductOf(B, C);
        Ia2.beProductOf(Ia, B);
        Ia2.times(dA);
        BTCB.plus(Ia2);

        // Build G in elemental stiffness matrix
        Ib.beDyadicProductOf(dN_V, Nlin);
        Ib.times(-dA);
        G.plus(Ib);
    }

    GT.beTranspositionOf(G);

    FloatMatrix temp(15, 15);
    temp.addSubMatrix(BTCB, 1, 1);
    temp.addSubMatrix(GT, 13, 1);
    temp.addSubMatrix(G, 1, 13);

    answer.resize(15, 15);
    answer.zero();
    answer.assemble(& temp, & this->ordering);
}

void Tr21Stokes :: updateYourself(TimeStep *tStep)
{
    Element :: updateYourself(tStep);
}

// Some extension Interfaces to follow:

Interface *Tr21Stokes :: giveInterface(InterfaceType it)
{
    if ( it == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    } else {
        return FMElement :: giveInterface(it);
    }
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
            } else if ( node == 6 ) {
                a = this->giveNode(3)->giveDofWithID(P_f)->giveUnknown(EID_ConservationEquation, VM_Total, tStep);
                b = this->giveNode(1)->giveDofWithID(P_f)->giveUnknown(EID_ConservationEquation, VM_Total, tStep);
            }

            answer.at(1) = ( a + b ) / 2;
        }
    } else   {
        answer.resize(0);
    }
}
} // end namespace oofem
