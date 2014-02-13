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
#include "load.h"
#include "timestep.h"
#include "boundaryload.h"
#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdio.h>
#endif
#include "contextioerr.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

namespace oofem {
PFEMElement2d :: PFEMElement2d(int n, Domain *aDomain) :
    PFEMElement(n, aDomain)

{
    // Constructor.
}

PFEMElement2d :: ~PFEMElement2d()
// Destructor
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

// NOT IN USE
void
PFEMElement2d :: computeNuMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    int i;
    FloatArray n;

    this->giveVelocityInterpolation()->evalN( n, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    answer.resize( 2, 2 * n.giveSize() );
    answer.zero();

    for ( i = 1; i <= n.giveSize(); i++ ) {
        answer.at(1, 2 * i - 1)  = n.at(i);
        answer.at(2, 2 * i - 0)  = n.at(i);
    }

    return;
}

// NOT IN USE
void
PFEMElement2d :: computeNpMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    int i;
    FloatArray n;

    this->givePressureInterpolation()->evalN( n, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    answer.resize( 1, n.giveSize() );
    answer.zero();

    for ( i = 1; i <= n.giveSize(); i++ ) {
        answer.at(1, i)  = n.at(i);
    }

    return;
}

void
PFEMElement2d :: computeBMatrix(FloatMatrix &answer, GaussPoint *gp)
//
// Returns the [3x6] strain-displacement matrix {B} of the receiver,
// evaluated at aGaussPoint.
// (epsilon_x,epsilon_y,gamma_xy) = B . r
// r = ( u1,v1,u2,v2,u3,v3)
{
    int i;
    FloatMatrix dnx;

    this->giveVelocityInterpolation()->evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize( 3, 2 * dnx.giveNumberOfRows() );
    answer.zero();

    for ( i = 1; i <= dnx.giveNumberOfRows(); i++ ) {
        answer.at(1, 2 * i - 1) = dnx.at(i, 1);
        answer.at(2, 2 * i - 0) = dnx.at(i, 2);
        answer.at(3, 2 * i - 1) = dnx.at(i, 2);
        answer.at(3, 2 * i - 0) = dnx.at(i, 1);
    }

    return;
}

void
PFEMElement2d :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *atTime)
{
    //setting sizes????
    FloatMatrix B, D, BTD, BTDB;

    //D.resize(3,3);
    //D.beUnitMatrix();

    int i;
    double dV;
    GaussPoint *gp;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];

    //double mu = static_cast< FluidDynamicMaterial * >( this->giveMaterial() )->giveEffectiveViscosity(iRule->getIntegrationPoint(0), atTime);
    //D.times(mu);

    for ( i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        gp  = iRule->getIntegrationPoint(i);
        ( ( FluidDynamicMaterial * ) this->giveMaterial() )->giveDeviatoricStiffnessMatrix(D, mode, gp, atTime);
        this->computeBMatrix(B, gp);
        BTD.beTProductOf(B, D);
        BTDB.beProductOf(BTD, B);
        dV  = this->computeVolumeAround(gp);
        BTDB.times(dV);
        answer.add(BTDB);
    }

    //????? only for the last gaus point

    //double mu = this->giveMaterial()->giveCharacteristicValue(MRM_Viscosity, gp, atTime);

    //according the paper "PFEM: A powerful tool to solve ..." should not be contained
    //double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, gp, atTime);
    //answer.times(mu/rho);

    //answer.times(mu*rho);

    //  answer.times(mu);
    return;
}

// NOT IN USE
void
PFEMElement2d :: computeVelocityLaplacianMatrix(FloatMatrix &answer, TimeStep *atTime)
{
    answer.zero();

    int i, j;
    FloatMatrix dnx;
    FloatMatrix temp;

    int ip;
    double dV;
    GaussPoint *gp;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];

    for ( ip = 0; ip < iRule->giveNumberOfIntegrationPoints(); ip++ ) {
        gp  = iRule->getIntegrationPoint(ip);
        temp.zero();
        double mu = static_cast< FluidDynamicMaterial * >( this->giveMaterial() )->giveEffectiveViscosity(gp, atTime);

        this->giveVelocityInterpolation()->evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
        temp.resizeWithData( 2 * dnx.giveNumberOfRows(), 2 * dnx.giveNumberOfRows() );

        for ( i = 1; i <= dnx.giveNumberOfRows(); i++ ) {
            for ( j = 1; j <= dnx.giveNumberOfRows(); j++ ) {
#if 0
                temp.at(2 * i - 1, 2 * j - 1) = dnx.at(i, 1) * dnx.at(j, 1);
                temp.at(2 * i, 2 * j) = dnx.at(i, 2) * dnx.at(j, 2);

                /*
                 * //according Ryzshakov's thesis - not working
                 * temp.at(2*i-1, 2*j-1) = (dnx.at(i,1)+dnx.at(i,2))*(dnx.at(j,1)+0.5*dnx.at(j,2));
                 * temp.at(2*i-1, 2*j) = (dnx.at(i,1)+dnx.at(i,2))*0.5*dnx.at(j,1);
                 * temp.at(2*i, 2*j-1) = (dnx.at(i,1)+dnx.at(i,2))*0.5*dnx.at(j,2);
                 * temp.at(2*i, 2*j) = (dnx.at(i,1)+dnx.at(i,2))*(0.5*dnx.at(j,1)+dnx.at(j,2));
                 *
                 */
#else
                //this is acording the kratos implemantation

                temp.at(2 * i - 1, 2 * j - 1) = dnx.at(i, 1) * dnx.at(j, 1) + dnx.at(i, 2) * dnx.at(j, 2);
                temp.at(2 * i, 2 * j) = dnx.at(i, 1) * dnx.at(j, 1) + dnx.at(i, 2) * dnx.at(j, 2);
#endif
            }
        }

        temp.times(mu);
        dV  = this->computeVolumeAround(gp);
        temp.times(dV);
        answer.add(temp);
    }

    return;
}

void
PFEMElement2d :: computePressureLaplacianMatrix(FloatMatrix &answer, TimeStep *atTime)
{
    answer.zero();

    int i, j;
    FloatMatrix dnx;
    FloatMatrix temp;

    int ip;
    double dV;
    GaussPoint *gp;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];

    for ( ip = 0; ip < iRule->giveNumberOfIntegrationPoints(); ip++ ) {
        temp.zero();
        gp  = iRule->getIntegrationPoint(ip);
        this->givePressureInterpolation()->evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
        temp.resizeWithData( dnx.giveNumberOfRows(), dnx.giveNumberOfRows() );
        for ( i = 1; i <= dnx.giveNumberOfRows(); i++ ) {
            for ( j = 1; j <= dnx.giveNumberOfRows(); j++ ) {
                temp.at(i, j) = dnx.at(i, 1) * dnx.at(j, 1) + dnx.at(i, 2) * dnx.at(j, 2);
            }
        }

        dV  = this->computeVolumeAround(gp);
        temp.times(dV);
        answer.add(temp);
    }

    double rho = this->giveMaterial()->give( 'd', integrationRulesArray [ 0 ]->getIntegrationPoint(0) );

    answer.times(1.0 / rho);
    return;
}


// NOT IN USE
void
PFEMElement2d :: computeStabilizedLaplacianMatrix(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(3, 3);
    answer.zero();

    int i, j; // ind;

    FloatArray tau;
    //  FloatArray tau(2);
    FloatMatrix dnx; //(3,2)
    FloatMatrix temp(3, 3);

    int ip;
    double dV;
    GaussPoint *gp;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];

    for ( ip = 0; ip < iRule->giveNumberOfIntegrationPoints(); ip++ ) {
        temp.zero();
        gp  = iRule->getIntegrationPoint(ip);
        this->computeStabilizationParameters(tau, gp, atTime);
        this->givePressureInterpolation()->evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
        //    for(ind = 1; ind<= tau.giveSize(); ind++){
        for ( i = 1; i <= dnx.giveNumberOfRows(); i++ ) {
            for ( j = 1; j <= dnx.giveNumberOfRows(); j++ ) {
                //temp.at(i,j) += dnx.at(i,ind)*dnx.at(j,ind)*tau.at(ind);
                temp.at(i, j) = dnx.at(i, 1) * dnx.at(j, 1) * tau.at(1) + dnx.at(i, 2) * dnx.at(j, 2) * tau.at(2);
            }
        }

        // }
        dV  = this->computeVolumeAround(gp);
        temp.times(dV);
        answer.add(temp);
    }

    double rho = this->giveMaterial()->give( 'd', integrationRulesArray [ 0 ]->getIntegrationPoint(0) );
    answer.times(1.0 / rho);

    return;
}

void
PFEMElement2d :: computeGradientMatrix(FloatMatrix &answer, TimeStep *atTime) //G
{
    answer.zero();

    int i, j;
    FloatArray N; //(3)
    FloatMatrix dnx; //(3,2)
    FloatMatrix temp;

    int ip;
    double dV;
    GaussPoint *gp;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];

    for ( ip = 0; ip < iRule->giveNumberOfIntegrationPoints(); ip++ ) {
        temp.zero();
        gp  = iRule->getIntegrationPoint(ip);

        this->givePressureInterpolation()->evalN( N, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
        this->givePressureInterpolation()->evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

        //?????????????? this is strange ????????????
        // resizing in each loop step
        answer.resizeWithData( 2 * dnx.giveNumberOfRows(), N.giveSize() );
        temp.resizeWithData( 2 * dnx.giveNumberOfRows(), N.giveSize() );

        for ( i = 1; i <= dnx.giveNumberOfRows(); i++ ) {
            for ( j = 1; j <= N.giveSize(); j++ ) {
                temp.at(2 * i - 1, j) = N.at(i) * dnx.at(j, 1);
                temp.at(2 * i, j) =  N.at(i) * dnx.at(j, 2);
            }
        }

        dV  = this->computeVolumeAround(gp);
        temp.times(dV);
        answer.add(temp);
    }

    return;
}

void
PFEMElement2d :: computeDivergenceMatrix(FloatMatrix &answer, TimeStep *atTime) //D = G^T
{
    int i, j;
    FloatArray N; //(3)
    FloatMatrix dnx; //(3,2)
    FloatMatrix temp;

    int ip;
    double dV;
    GaussPoint *gp;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];

    for ( ip = 0; ip < iRule->giveNumberOfIntegrationPoints(); ip++ ) {
        temp.zero();
        gp  = iRule->getIntegrationPoint(ip);

        this->giveVelocityInterpolation()->evalN( N, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
        this->giveVelocityInterpolation()->evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

        //?????????????? this is strange ????????????
        // resizing in each loop step
        answer.resizeWithData( N.giveSize(), 2 * dnx.giveNumberOfRows() );
        temp.resizeWithData( N.giveSize(), 2 * dnx.giveNumberOfRows() );

        for ( i = 1; i <= N.giveSize(); i++ ) {
            for ( j = 1; j <= dnx.giveNumberOfRows(); j++ ) {
                temp.at(i, 2 * j - 1) = dnx.at(i, 1) * N.at(j);
                temp.at(i, 2 * j) = dnx.at(i, 2) * N.at(j);
            }
        }

        dV  = this->computeVolumeAround(gp);
        temp.times(dV);
        answer.add(temp);
    }

    return;
}

// NOT IN USE
void
PFEMElement2d :: computeStabilizationGradientMatrix(FloatMatrix &answer, TimeStep *atTime) //Q
{
    answer.zero();
    FloatArray tau;


    int i, j;
    FloatArray N; //(3)
    FloatMatrix dnx; //(3,2)
    FloatMatrix temp;

    int ip;
    double dV;
    GaussPoint *gp;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];

    for ( ip = 0; ip < iRule->giveNumberOfIntegrationPoints(); ip++ ) {
        temp.zero();
        gp  = iRule->getIntegrationPoint(ip);
        this->computeStabilizationParameters(tau, gp, atTime);

        this->givePressureInterpolation()->evalN( N, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
        this->givePressureInterpolation()->evaldNdx( dnx, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );

        //?????????????? this is strange ????????????
        // resizing in each loop step
        answer.resizeWithData( dnx.giveNumberOfRows(), 2 * N.giveSize() );
        temp.resizeWithData( dnx.giveNumberOfRows(), 2 * N.giveSize() );

        for ( i = 1; i <= dnx.giveNumberOfRows(); i++ ) {
            for ( j = 1; j <= N.giveSize(); j++ ) {
                temp.at(i, 2 * j - 1) = dnx.at(i, 1) * N.at(j) * tau.at(1);
                temp.at(i, 2 * j) = dnx.at(i, 2) * N.at(j) * tau.at(2);
            }
        }

        dV  = this->computeVolumeAround(gp);
        temp.times(dV);
        answer.add(temp);
    }

    return;
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
    int i, numberOfGaussPoints;
    double dV;
    FloatMatrix T;
    FloatArray globalIPcoords;

    numberOfGaussPoints = 3;
    GaussIntegrationRule iRule(1, this, 1, 1);
    iRule.setUpIntegrationPoints(_Line, numberOfGaussPoints, _Unknown);
    GaussPoint *gp;
    FloatArray reducedAnswer, force, ntf, N, temp;
    IntArray mask;
    FloatMatrix Nmtrx;

    FloatArray u_presq, u;
    this->computeVectorOfPrescribed(EID_MomentumBalance, mode, tStep, u_presq);

    for ( int iEdge = 1; iEdge <= 3; iEdge++ ) {
        reducedAnswer.zero();
        for ( i = 0; i < iRule.giveNumberOfIntegrationPoints(); i++ ) {
            gp  = iRule.getIntegrationPoint(i);
            this->computeEgdeNVectorAt(N, iEdge, gp);
            this->computeEdgeNMatrixAt(Nmtrx, iEdge, gp);
            dV  = this->computeEdgeVolumeAround(gp, iEdge);
            FloatArray *lcoords = gp->giveCoordinates();
            FloatArray normal;
            this->giveVelocityInterpolation()->boundaryEvalNormal( normal, iEdge, * lcoords, FEIElementGeometryWrapper(this) );

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
            temp = N;
            temp.times(un);
            reducedAnswer.add(dV, temp);
        }

        this->giveEdgeDofMapping(mask, iEdge);
        answer.assemble(reducedAnswer, mask);
    }

    return;
}

// NOT FINISHED !!! NOT IN USE, BECAUSE NOT NEEDED
void
PFEMElement2d :: computePrescribedPressureRhsVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
    answer.resize(6);
    answer.zero();

    int i;
    int numberOfGaussPoints = 3;
    //double dV;
    FloatMatrix T;
    FloatArray globalIPcoords;

    GaussIntegrationRule iRule(1, this, 1, 1);
    iRule.setUpIntegrationPoints(_Line, numberOfGaussPoints, _Unknown);
    GaussPoint *gp;
    FloatArray reducedAnswer, force, ntf, N, temp;
    IntArray mask;
    FloatMatrix Nmtrx;

    FloatArray p_presq;
    this->computeVectorOf(EID_ConservationEquation, VM_Total, tStep, p_presq);

    double p_norm = p_presq.computeNorm();
    if ( p_norm > 1.e-2 ) {
        //int kl = 1;

        for ( int iEdge = 1; iEdge <= 3; iEdge++ ) {
            reducedAnswer.zero();
            for ( i = 0; i < iRule.giveNumberOfIntegrationPoints(); i++ ) {
                gp  = iRule.getIntegrationPoint(i);
                this->computeEgdeNVectorAt(N, iEdge, gp);
                this->computeEdgeNMatrixAt(Nmtrx, iEdge, gp);
                //dV  = this->computeEdgeVolumeAround(gp, iEdge);
            }
        }
    }
}

void
PFEMElement2d :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
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
        _error("giveEdgeDofMapping: wrong edge number");
    }
}
void
PFEMElement2d :: computeEdgeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
{
    FloatArray n;
    this->giveVelocityInterpolation()->boundaryEvalN( n, iedge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    answer.resize(2, 4);
    answer.at(1, 1) = answer.at(2, 2) = n.at(1);
    answer.at(1, 3) = answer.at(2, 4) = n.at(2);
}

void
PFEMElement2d :: computeEgdeNVectorAt(FloatArray &answer, int iedge, GaussPoint *gp)
{
    FloatArray n;
    this->giveVelocityInterpolation()->boundaryEvalN( n, iedge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    answer.resize(2);
    answer.at(1) = n.at(1);
    answer.at(2) = n.at(2);
}

double
PFEMElement2d :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double detJ = this->giveVelocityInterpolation()->boundaryGiveTransformationJacobian( iEdge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    return detJ *gp->giveWeight();
}


// NOT IN USE
void
PFEMElement2d :: computeStabilizationMassMatrix(FloatMatrix &answer, TimeStep *atTime) //M_hat
{
    /*
     * answer.zero();
     * int i,j;
     * FloatArray N;//(3)
     * FloatArray tau;
     * FloatMatrix temp;
     *
     * int ip;
     * double dV;
     * GaussPoint *gp;
     * IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
     *
     * for ( ip = 0; ip < iRule->getNumberOfIntegrationPoints(); ip++ ) {
     * temp.zero();
     * gp  = iRule->getIntegrationPoint(ip);
     * this->givePressureInterpolation()->evalN(N, * gp->giveCoordinates(), FEIElementGeometryWrapper(this), 0.0);
     * this -> computeStabilizationParameters(tau, gp, atTime);
     *
     * answer.resize(2 * N.giveSize(),2 * N.giveSize());
     * temp.resize(2 * N.giveSize(),2 * N.giveSize());
     *
     * for(i = 1; i <= 3; i++){
     *  for(j = 1; j <= 3; j++){
     *    temp.at(2*i-1, 2*j-1) = N.at(i) * N.at(j) * tau.at(1);
     *    temp.at(2*i, 2*j) = N.at(i) * N.at(j) * tau.at(2);
     *    //temp.at(i,j) = N.at(i) * N.at(j) * tau.at(1);
     *    //		temp.at(i+ N.giveSize(),j+ N.giveSize()) = N.at(i) * N.at(j) * tau.at(2);
     *  }
     * }
     * dV  = this->computeVolumeAround(gp);
     * temp.times(dV);
     * answer.add(temp);
     * }
     */
    return;
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
    //<RESTRICTED_SECTION>

    //</RESTRICTED_SECTION>
}



contextIOResultType PFEMElement2d :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = PFEMElement :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}



contextIOResultType PFEMElement2d :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = PFEMElement :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}
} // end namespace oofem
