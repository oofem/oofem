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

#include "pfemelement2d.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "engngm.h"
#include "fluiddynamicmaterial.h"
#include "fluidcrosssection.h"
#include "load.h"
#include "timestep.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "contextioerr.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

namespace oofem {
PFEMElement2d :: PFEMElement2d(int n, Domain *aDomain) :
    PFEMElement(n, aDomain)
{
}

PFEMElement2d :: ~PFEMElement2d()
{ }


IRResultType
PFEMElement2d :: initializeFrom(InputRecord *ir)
{
    return this->PFEMElement :: initializeFrom(ir);
}

int
PFEMElement2d :: checkConsistency()
{
    return PFEMElement :: checkConsistency();
}


void
PFEMElement2d :: computeBMatrix(FloatMatrix &answer, GaussPoint *gp)
//
// Returns the [3x6] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
// (epsilon_x,epsilon_y,gamma_xy) = B . r
// r = ( u1,v1,u2,v2,u3,v3)
{
    FloatMatrix dnx;

    this->giveVelocityInterpolation()->evaldNdx( dnx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize( 3, 2 * dnx.giveNumberOfRows() );
    answer.zero();

    for ( int i = 1; i <= dnx.giveNumberOfRows(); i++ ) {
        answer.at(1, 2 * i - 1) = dnx.at(i, 1);
        answer.at(2, 2 * i - 0) = dnx.at(i, 2);
        answer.at(3, 2 * i - 1) = dnx.at(i, 2);
        answer.at(3, 2 * i - 0) = dnx.at(i, 1);
    }
}

void
PFEMElement2d :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *atTime)
{
    FloatMatrix B, D, DB;
    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();

    answer.clear();
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ].get();
    for ( auto &gp : *iRule ) {
        mat->giveDeviatoricStiffnessMatrix(D, mode, gp, atTime);
        this->computeBMatrix(B, gp);
        DB.beProductOf(D, B);
        double dV = this->computeVolumeAround(gp);
        answer.plusProductUnsym(B, DB, dV);
    }
}


void
PFEMElement2d :: computePressureLaplacianMatrix(FloatMatrix &answer, TimeStep *atTime)
{
    answer.zero();

    FloatMatrix dnx;
    FloatMatrix temp;

    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ].get();

    for ( auto &gp : *iRule ) {
        this->givePressureInterpolation()->evaldNdx( dnx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

        temp.resize( dnx.giveNumberOfRows(), dnx.giveNumberOfRows() );
        temp.zero();
        for ( int i = 1; i <= dnx.giveNumberOfRows(); i++ ) {
            for ( int j = 1; j <= dnx.giveNumberOfRows(); j++ ) {
                temp.at(i, j) = dnx.at(i, 1) * dnx.at(j, 1) + dnx.at(i, 2) * dnx.at(j, 2);
            }
        }

        double dV = this->computeVolumeAround(gp);
        answer.add(dV, temp);
    }

    double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity( integrationRulesArray [ 0 ]->getIntegrationPoint(0) );

    answer.times(1.0 / rho);
}


void
PFEMElement2d :: computeGradientMatrix(FloatMatrix &answer, TimeStep *atTime) //G
{
    answer.clear();

    FloatArray N; //(3)
    FloatMatrix dnx; //(3,2)
    FloatMatrix temp;

    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ].get();

    for ( auto &gp : *iRule ) {

        this->givePressureInterpolation()->evalN( N, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
        this->givePressureInterpolation()->evaldNdx( dnx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

        temp.resize( 2 * dnx.giveNumberOfRows(), N.giveSize() );
        temp.zero();
        for ( int i = 1; i <= dnx.giveNumberOfRows(); i++ ) {
            for ( int j = 1; j <= N.giveSize(); j++ ) {
                temp.at(2 * i - 1, j) = N.at(i) * dnx.at(j, 1);
                temp.at(2 * i, j) =  N.at(i) * dnx.at(j, 2);
            }
        }

        double dV = this->computeVolumeAround(gp);
        answer.add(dV, temp);
    }
}


void
PFEMElement2d :: computeDivergenceMatrix(FloatMatrix &answer, TimeStep *atTime) //D = G^T
{
    FloatArray N; //(3)
    FloatMatrix dnx; //(3,2)
    FloatMatrix temp;

    answer.clear();

    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ].get();
    for ( auto &gp : *iRule ) {

        this->giveVelocityInterpolation()->evalN( N, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
        this->giveVelocityInterpolation()->evaldNdx( dnx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

        temp.resize( N.giveSize(), 2 * dnx.giveNumberOfRows() );
        temp.zero();
        for ( int i = 1; i <= N.giveSize(); i++ ) {
            for ( int j = 1; j <= dnx.giveNumberOfRows(); j++ ) {
                temp.at(i, 2 * j - 1) = dnx.at(i, 1) * N.at(j);
                temp.at(i, 2 * j) = dnx.at(i, 2) * N.at(j);
            }
        }

        double dV = this->computeVolumeAround(gp);
        answer.add(dV, temp);
    }
}


void
PFEMElement2d :: computePrescribedRhsVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
    // SHOULD BE MOVED INTO TR1_2D_PFEM

    answer.resize(3);
    answer.zero();
    // computes numericaly the edge load vector of the receiver for given load
    // Each element edge must have unique number assigned to identify it.
    // Integration is done in local edge space (i.e. one dimensional integration is
    // performed on line). This general implementation requires that element must
    // provide following functions:
    // - ComputeEgdeNMatrixAt - returns interpolation matrix of the edge in the
    //   local edge space.
    // - computeEdgeVolumeAround - returns volumeAround in local edge space
    // - GiveEdgeDofMapping - returns integer array specifying local edge dof mapping to
    //   element dofs.
    //
    int numberOfGaussPoints;
    FloatMatrix T;
    FloatArray globalIPcoords;

    numberOfGaussPoints = 3;
    GaussIntegrationRule iRule(1, this, 1, 1);
    iRule.setUpIntegrationPoints(_Line, numberOfGaussPoints, _Unknown);
    FloatArray reducedAnswer, force, ntf, N;
    IntArray mask;
    FloatMatrix Nmtrx;

    FloatArray u_presq, u;
    this->computeVectorOfPrescribed({V_u, V_v}, mode, tStep, u_presq);

    for ( int iEdge = 1; iEdge <= 3; iEdge++ ) {
        reducedAnswer.zero();
        for ( auto &gp : iRule ) {
            this->computeEgdeNVectorAt(N, iEdge, gp);
            this->computeEdgeNMatrixAt(Nmtrx, iEdge, gp);
            double dV = this->computeEdgeVolumeAround(gp, iEdge);
            FloatArray normal;
            this->giveVelocityInterpolation()->boundaryEvalNormal( normal, iEdge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

            FloatArray u_presq_edge(4);
            if ( iEdge == 1 ) {
                u_presq_edge.at(1) = u_presq.at(1);
                u_presq_edge.at(2) = u_presq.at(2);
                u_presq_edge.at(3) = u_presq.at(3);
                u_presq_edge.at(4) = u_presq.at(4);
            } else if ( iEdge == 2 ) {
                u_presq_edge.at(1) = u_presq.at(3);
                u_presq_edge.at(2) = u_presq.at(4);
                u_presq_edge.at(3) = u_presq.at(5);
                u_presq_edge.at(4) = u_presq.at(6);
            } else {
                u_presq_edge.at(1) = u_presq.at(5);
                u_presq_edge.at(2) = u_presq.at(6);
                u_presq_edge.at(3) = u_presq.at(1);
                u_presq_edge.at(4) = u_presq.at(2);
            }

            u.beProductOf(Nmtrx, u_presq_edge);
            double un = normal.dotProduct(u);
            reducedAnswer.add(dV * un, N);
        }

        this->giveEdgeDofMapping(mask, iEdge);
        answer.assemble(reducedAnswer, mask);
    }
}

void
PFEMElement2d :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    // SHOULD BE MOVED INTO TR1_2D_PFEM
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    answer.resize(2);
    if ( iEdge == 1 ) { // edge between nodes 1,2
        answer.at(1) = 1;
        answer.at(2) = 2;
    } else if ( iEdge == 2 ) { // edge between nodes 2 3
        answer.at(1) = 2;
        answer.at(2) = 3;
    } else if ( iEdge == 3 ) { // edge between nodes 3 1
        answer.at(1) = 3;
        answer.at(2) = 1;
    } else {
        OOFEM_ERROR("giveEdgeDofMapping: wrong edge number");
    }
}


void
PFEMElement2d :: computeEdgeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
{
    FloatArray n;
    this->giveVelocityInterpolation()->boundaryEvalN( n, iedge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(n, 2);
}

void
PFEMElement2d :: computeEgdeNVectorAt(FloatArray &answer, int iedge, GaussPoint *gp)
{
    this->giveVelocityInterpolation()->boundaryEvalN(answer, iedge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
}

double
PFEMElement2d :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double detJ = this->giveVelocityInterpolation()->boundaryGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    return detJ *gp->giveWeight();
}


void
PFEMElement2d :: updateYourself(TimeStep *tStep)
{
    PFEMElement :: updateYourself(tStep);
}


void
PFEMElement2d :: printOutputAt(FILE *file, TimeStep *stepN)
// Performs end-of-step operations.
{
    PFEMElement :: printOutputAt(file, stepN);
}


contextIOResultType PFEMElement2d :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( ( iores = PFEMElement :: saveContext(*stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType PFEMElement2d :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;

    if ( ( iores = PFEMElement :: restoreContext(*stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}
} // end namespace oofem
