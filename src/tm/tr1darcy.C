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


#include "tr1darcy.h"
#include "fei2dtrlin.h"
#include "gaussintegrationrule.h"
#include "gausspoint.h"
#include "bcgeomtype.h"
#include "generalboundarycondition.h"
#include "transportmaterial.h"
#include "load.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "crosssection.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(Tr1Darcy);

FEI2dTrLin Tr1Darcy :: interpolation_lin(1, 2);

Tr1Darcy :: Tr1Darcy(int n, Domain *aDomain) : TransportElement(n, aDomain)
{
    numberOfDofMans = 3;
}

Tr1Darcy :: ~Tr1Darcy()
{ }

IRResultType Tr1Darcy :: initializeFrom(InputRecord *ir)
{
    this->numberOfGaussPoints = 1;
    return TransportElement :: initializeFrom(ir);
}

FEInterpolation *
Tr1Darcy :: giveInterpolation() const
{
    return & interpolation_lin;
}

void Tr1Darcy :: computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}

void Tr1Darcy :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep)
{
    /*
     * Return Ke = integrate(B^T K B)
     */

    FloatMatrix B, BT, K, KB;

    TransportMaterial *mat = static_cast< TransportMaterial * >( this->giveMaterial() );

    answer.clear();

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        ///@todo Should we make it return the transpose instead?
        double detJ = fabs( this->interpolation_lin.evaldNdx( BT, lcoords, FEIElementGeometryWrapper(this) ) );

        mat->giveCharacteristicMatrix(K, mode, gp, tStep);

        B.beTranspositionOf(BT);
        KB.beProductOf(K, B);
        answer.plusProductUnsym( B, KB, detJ * gp->giveWeight() ); // Symmetric part is just a single value, not worth it.
    }
}

void Tr1Darcy :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode, TimeStep *tStep)
{
    if ( mtrx == ExternalForcesVector ) {
        this->computeExternalForcesVector(answer, tStep, mode);
    } else if ( mtrx == InternalForcesVector ) {
        this->computeInternalForcesVector(answer, tStep);
    } else {
        OOFEM_ERROR("Unknown Type of characteristic mtrx.");
    }
}

void Tr1Darcy :: computeInternalForcesVector(FloatArray &answer, TimeStep *tStep)
{
    FloatArray w, a, gradP, P(1), n;
    FloatMatrix B, BT;

    TransportMaterial *mat = static_cast< TransportMaterial * >( this->giveMaterial() );

    this->computeVectorOf(VM_Total, tStep, a);

    answer.resize(3);
    answer.zero();

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();

        double detJ = fabs( this->interpolation_lin.giveTransformationJacobian( lcoords, FEIElementGeometryWrapper(this) ) );
        this->interpolation_lin.evaldNdx( BT, lcoords, FEIElementGeometryWrapper(this) );
        this->interpolation_lin.evalN( n, lcoords, FEIElementGeometryWrapper(this) );
        B.beTranspositionOf(BT);
        P.at(1) = n.dotProduct(a); // Evaluates the field at this point.

        gradP.beProductOf(B, a);

        mat->giveFluxVector(w, gp, gradP, P, tStep);

        answer.plusProduct(B, w, -gp->giveWeight() * detJ);
    }
}

void Tr1Darcy :: computeExternalForcesVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
    // TODO: Implement support for body forces

    FloatArray vec;

    answer.resize(3);
    answer.zero();

    // Compute characteristic vector for Neumann boundary conditions.
    int load_number, load_id;
    Load *load;
    bcGeomType ltype;

    int nLoads = boundaryLoadArray.giveSize() / 2;

    for ( int i = 1; i <= nLoads; i++ ) {  // For each Neumann boundary condition ....
        load_number = boundaryLoadArray.at(2 * i - 1);
        load_id = boundaryLoadArray.at(2 * i);
        load = domain->giveLoad(load_number);
        ltype = load->giveBCGeoType();

        if ( ltype == EdgeLoadBGT ) {
            this->computeEdgeBCSubVectorAt(vec, load, load_id, tStep, mode, 0);
        }

        answer.add(vec);
    }

    answer.negated();
}

void Tr1Darcy :: computeEdgeBCSubVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep, ValueModeType mode, int indx)
{
    /*
     * Given the load *load, return it's contribution.
     *
     */

    answer.resize(3);
    answer.zero();

    if ( load->giveType() == TransmissionBC ) {                 // Neumann boundary conditions (traction)
        BoundaryLoad *boundaryLoad;
        boundaryLoad = static_cast< BoundaryLoad * >(load);

        int numberOfEdgeIPs;
        numberOfEdgeIPs = ( int ) ceil( ( boundaryLoad->giveApproxOrder() + 1. ) / 2. ) * 2;

        GaussIntegrationRule iRule(1, this, 1, 1);
        FloatArray N, loadValue, reducedAnswer;
        reducedAnswer.resize(3);
        reducedAnswer.zero();
        IntArray mask;

        iRule.SetUpPointsOnLine(numberOfEdgeIPs, _Unknown);

        for ( GaussPoint *gp: iRule ) {
            const FloatArray &lcoords = gp->giveNaturalCoordinates();
            this->interpolation_lin.edgeEvalN( N, iEdge, lcoords, FEIElementGeometryWrapper(this) );
            double dV = this->computeEdgeVolumeAround(gp, iEdge);

            if ( boundaryLoad->giveFormulationType() == Load :: FT_Entity ) {                // Edge load in xi-eta system
                boundaryLoad->computeValueAt(loadValue, tStep, lcoords, mode);
            } else {  // Edge load in x-y system
                FloatArray gcoords;
                this->interpolation_lin.edgeLocal2global( gcoords, iEdge, lcoords, FEIElementGeometryWrapper(this) );
                boundaryLoad->computeValueAt(loadValue, tStep, gcoords, mode);
            }

            reducedAnswer.add(loadValue.at(1) * dV, N);
        }

        this->interpolation_lin.computeLocalEdgeMapping(mask, iEdge);
        answer.assemble(reducedAnswer, mask);
    }
}

double Tr1Darcy :: giveThicknessAt(const FloatArray &gcoords)
{
    return this->giveCrossSection()->give(CS_Thickness, gcoords, this, false);
}


double Tr1Darcy :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double thickness = 1;
    double detJ = fabs( this->interpolation_lin.edgeGiveTransformationJacobian( iEdge, gp->giveSubPatchCoordinates(), FEIElementGeometryWrapper(this) ) );
    return detJ *thickness *gp->giveWeight();
}

void Tr1Darcy :: giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep)
{
    /*
     * Compute characteristic matrix for this element. The only option is the stiffness matrix...
     */
    if ( mtrx == ConductivityMatrix || mtrx == TangentStiffnessMatrix ) {
        this->computeStiffnessMatrix(answer, TangentStiffness, tStep);
    } else {
        OOFEM_ERROR("Unknown Type of characteristic mtrx %d", mtrx);
    }
}

void Tr1Darcy :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {P_f};
}

void
Tr1Darcy :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                       InternalStateType type, TimeStep *tStep)
{
    CrossSection *cs = this->giveCrossSection();
    ///@todo Write support function for getting the closest gp given local c.s. and use that here
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    cs->giveIPValue(answer, gp, type, tStep);
}

Interface *
Tr1Darcy :: giveInterface(InterfaceType interface)
{
    if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    }

    return NULL;
}

int
Tr1Darcy :: computeNumberOfDofs()
{
    return 3;
}
}
