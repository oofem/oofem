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


#include "tr1darcy.h"
#include "node.h"
#include "domain.h"
#include "equationid.h"
#include "gaussintegrationrule.h"
#include "gausspnt.h"
#include "bcgeomtype.h"
#include "generalbc.h"
#include "load.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "crosssection.h"
#include "matresponseform.h"
#include "matresponsemode.h"
#include "fei2dtrlin.h"

namespace oofem {
FEI2dTrLin Tr1Darcy :: interpolation_lin(1, 2);

Tr1Darcy :: Tr1Darcy(int n, Domain *aDomain) : TransportElement(n, aDomain)
{
    numberOfDofMans  = 3;
    numberOfGaussPoints = 1;
    this->computeGaussPoints();
}

Tr1Darcy :: ~Tr1Darcy()   
{
}

IRResultType Tr1Darcy :: initializeFrom(InputRecord *ir)
{
    this->TransportElement :: initializeFrom(ir);
    return IRRT_OK;
}

void Tr1Darcy :: computeGaussPoints()
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Triangle, numberOfGaussPoints, _2dHeat);
    }
}

void Tr1Darcy :: computeStiffnessMatrix(FloatMatrix &answer, TimeStep *atTime)
{
    /*
     * Return Ke = integrate(B^T K B)
     */

    FloatMatrix B, BT, K, KB;
    FloatArray *lcoords;
    GaussPoint *gp;

    TransportMaterial *mat = ( TransportMaterial * ) this->domain->giveMaterial(this->material);

    IntegrationRule *iRule = integrationRulesArray [ 0 ];

    answer.resize(3, 3);
    answer.zero();

    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        lcoords = gp->giveCoordinates();

        double detJ = this->interpolation_lin.giveTransformationJacobian( * lcoords, FEIElementGeometryWrapper(this) );
        this->interpolation_lin.evaldNdx( BT, * lcoords, FEIElementGeometryWrapper(this) );
        
        mat->giveCharacteristicMatrix(K, FullForm, TangentStiffness, gp, atTime);

        B.beTranspositionOf(BT);
        KB.beProductOf(K, B);
        answer.plusProductUnsym(B, KB, detJ * gp->giveWeight() ); // Symmetric part is just a single value, not worth it.
    }
}

void Tr1Darcy :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode, TimeStep *tStep)
{
    if ( mtrx == ExternalForcesVector ) {
        this->computeLoadVector(answer, tStep);
    } else if ( mtrx == InternalForcesVector ) {
        this->computeInternalForcesVector(answer, tStep);
    } else {
        _error("giveCharacteristicVector: Unknown Type of characteristic mtrx.");
    }
}

void Tr1Darcy :: computeInternalForcesVector(FloatArray &answer, TimeStep *atTime)
{
    FloatArray *lcoords, w, a, gradP, I;
    FloatMatrix BT;
    GaussPoint *gp;

    TransportMaterial *mat = ( TransportMaterial * ) this->domain->giveMaterial(this->material);
    IntegrationRule *iRule = integrationRulesArray [ 0 ];

    this->computeVectorOf(EID_ConservationEquation, VM_Total, atTime, a);

    answer.resize(3);
    answer.zero();

    for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        lcoords = gp->giveCoordinates();

        double detJ = this->interpolation_lin.giveTransformationJacobian( * lcoords, FEIElementGeometryWrapper(this) );
        this->interpolation_lin.evaldNdx( BT, * lcoords, FEIElementGeometryWrapper(this) );

        gradP.beTProductOf(BT, a);

        mat->giveFluxVector(w, gp, gradP, atTime);

        I.beProductOf(BT, w);
        answer.add(- gp->giveWeight() * detJ, I);
    }
}

void Tr1Darcy :: computeLoadVector(FloatArray &answer, TimeStep *atTime)
{
    // TODO: Implement support for body forces

    FloatArray vec;

    answer.resize(3);
    answer.zero();

    // Compute characteristic vector for Neumann boundary conditions.
    int i, load_number, load_id;
    GeneralBoundaryCondition *load;
    bcGeomType ltype;

    int nLoads = boundaryLoadArray.giveSize() / 2;

    for ( i = 1; i <= nLoads; i++ ) {  // For each Neumann boundary condition ....
        load_number = boundaryLoadArray.at(2 * i - 1);
        load_id = boundaryLoadArray.at(2 * i);
        load = ( GeneralBoundaryCondition * ) domain->giveLoad(load_number);
        ltype = load->giveBCGeoType();

        if ( ltype == EdgeLoadBGT ) {
            this->computeEdgeBCSubVectorAt(vec, ( Load * ) load, load_id, atTime);
        }

        answer.add(vec);
    }

    answer.negated();
}

void Tr1Darcy :: computeEdgeBCSubVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep)
{
    /*
     * Given the load *load, return it's contribution.
     *
     */

    answer.resize(3);
    answer.zero();

    if ( load->giveType() == TransmissionBC ) {                 // Neumann boundary conditions (traction)
        BoundaryLoad *boundaryLoad;
        boundaryLoad = ( BoundaryLoad * ) load;

        int numberOfEdgeIPs;
        numberOfEdgeIPs = ( int ) ceil( ( boundaryLoad->giveApproxOrder() + 1. ) / 2. ) * 2;

        GaussIntegrationRule iRule(1, this, 1, 1);
        GaussPoint *gp;
        FloatArray N, loadValue, reducedAnswer;
        reducedAnswer.resize(3);
        reducedAnswer.zero();
        IntArray mask;

        iRule.setUpIntegrationPoints(_Line, numberOfEdgeIPs, _Unknown);

        for ( int i = 0; i < iRule.getNumberOfIntegrationPoints(); i++ ) {
            gp = iRule.getIntegrationPoint(i);
            FloatArray *lcoords = gp->giveCoordinates();
            this->interpolation_lin.edgeEvalN(N, *lcoords, FEIElementGeometryWrapper(this));
            double dV = this->computeEdgeVolumeAround(gp, iEdge);

            if ( boundaryLoad->giveFormulationType() == BoundaryLoad :: BL_EntityFormulation ) {                // Edge load in xi-eta system
                boundaryLoad->computeValueAt(loadValue, tStep, *lcoords, VM_Total);
            } else {  // Edge load in x-y system
                FloatArray gcoords;
                this->interpolation_lin.edgeLocal2global(gcoords, iEdge, *lcoords, FEIElementGeometryWrapper(this));
                boundaryLoad->computeValueAt(loadValue, tStep, gcoords, VM_Total);
            }

            reducedAnswer.add(loadValue.at(1) * dV, N);
        }

        this->interpolation_lin.computeLocalEdgeMapping(mask, iEdge);
        answer.assemble(reducedAnswer, mask);
    }
}

double Tr1Darcy :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double thickness = 1;
    double detJ = fabs( this->interpolation_lin.edgeGiveTransformationJacobian(iEdge, *gp->giveLocalCoordinates(), FEIElementGeometryWrapper(this)) );
    return detJ * thickness * gp->giveWeight();
}

void Tr1Darcy :: giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep)
{
    /*
     * Compute characteristic matrix for this element. The only option is the stiffness matrix...
     */
    if ( mtrx == StiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, tStep);
    } else {
        _error("giveCharacteristicMatrix: Unknown Type of characteristic mtrx.");
    }
}

void Tr1Darcy :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    /*
     * Returns the mask for node number inode of this element. The mask tells what quantities
     * are held by each node. Since this element holds velocities (both in x and y direction),
     * in six nodes and pressure in three nodes the answer depends on which node is requested.
     */

    if ( ( inode == 1 ) || ( inode == 2 ) || ( inode == 3 ) ) {
        if ( ut == EID_ConservationEquation ) {
            answer.setValues(1, P_f);
        } else {
            _error("giveDofManDofIDMask: Unknown equation id encountered");
        }
    }
}

void
Tr1Darcy :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                       InternalStateType type, TimeStep *tStep)
{
    GaussPoint *gp;
    TransportMaterial *mat = ( TransportMaterial * ) this->domain->giveMaterial(this->material);

    IntegrationRule *iRule = integrationRulesArray [ 0 ];
    gp = iRule->getIntegrationPoint(0);
    mat->giveIPValue(answer, gp, type, tStep);
}

void
Tr1Darcy :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                      InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}

Interface *
Tr1Darcy :: giveInterface(InterfaceType interface)
{
    if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return ( NodalAveragingRecoveryModelInterface * ) this;
    }

    return NULL;
}

int
Tr1Darcy :: computeNumberOfDofs(EquationID ut)
{
    if ( ut == EID_MomentumBalance ) {
        return 6;
    } else if ( ut == EID_ConservationEquation ) {
        return 3;
    } else if ( ut == EID_MomentumBalance_ConservationEquation ) {
        return 9;
    } else {
        _error("computeNumberOfDofs: Unknown equation id encountered");
    }

    return 0;
}
}
