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

#include "tr1_2d_pfem.h"
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
 #include "conTable.h"
#endif

namespace oofem {
FEI2dTrLin TR1_2D_PFEM :: velocityInterpolation(1, 2);
FEI2dTrLin TR1_2D_PFEM :: pressureInterpolation(1, 2);

IntArray TR1_2D_PFEM :: ordering(9);
IntArray TR1_2D_PFEM :: edge_ordering [ 3 ] = {
    IntArray(6), IntArray(6), IntArray(6)
};

TR1_2D_PFEM :: TR1_2D_PFEM(int n, Domain *aDomain) :
    PFEMElement2d(n, aDomain)

{
    // Constructor.
    numberOfDofMans  = 3;
}

TR1_2D_PFEM :: TR1_2D_PFEM(int n, Domain *aDomain, int par1, int par2, int par3, int mat, int cs) :
    PFEMElement2d(n, aDomain)

{
    // Constructor.
    numberOfDofMans  = 3;
    IntArray dmans(3), aBodyLoadArry(1);
    dmans.at(1) = par1;
    dmans.at(2) = par2;
    dmans.at(3) = par3;
    aBodyLoadArry.at(1) = 3;
    this->setDofManagers(dmans);
    // CHECK THIS - NOT NICE
    this->setBodyLoads(aBodyLoadArry);
    this->setMaterial(mat);
    this->setCrossSection(cs);
    this->postInitialize();
}
TR1_2D_PFEM :: ~TR1_2D_PFEM()
// Destructor
{ }

int
TR1_2D_PFEM :: computeNumberOfDofs()
{
    return 9;
}

void
TR1_2D_PFEM ::   giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    // returns DofId mask array for inode element node.
    // DofId mask array determines the dof ordering requsted from node.
    // DofId mask array contains the DofID constants (defined in cltypes.h)
    // describing physical meaning of particular DOFs.
    if ( ( ut == EID_MomentumBalance ) || ( ut == EID_AuxMomentumBalance ) ) {
        answer.resize(2);
        answer.at(1) = V_u;
        answer.at(2) = V_v;
    } else if ( ut == EID_ConservationEquation ) {
        answer.resize(1);
        answer.at(1) = P_f;
    } else {
        _error("giveDofManDofIDMask: Unknown equation id encountered");
    }
}
// NOT IN USE
void
TR1_2D_PFEM ::   giveElementDofIDMask(EquationID ut, IntArray &answer) const
{
    this->giveDofManDofIDMask(1, ut, answer);
}


IRResultType
TR1_2D_PFEM :: initializeFrom(InputRecord *ir)
{
    //<RESTRICTED_SECTION>
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    //</RESTRICTED_SECTION>

    this->PFEMElement :: initializeFrom(ir);

    this->computeGaussPoints();
    return IRRT_OK;
}

void
TR1_2D_PFEM :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Triangle, 3, _2dFlow);
    }
}

// NOT IN USE
double
TR1_2D_PFEM :: giveCharacteristicLenght(GaussPoint *gp, const FloatArray &inDirection)
//
// returns receivers characteristic length in gp in given direction (implemented by element.C)
//
{
    return this->giveLenghtInDir(inDirection);
}

// NOT IN USE
void
TR1_2D_PFEM :: computeStabilizationParameters(FloatArray &answer, GaussPoint *gp, TimeStep *atTime)
{
    double mu = static_cast< FluidDynamicMaterial * >( this->giveMaterial() )->giveEffectiveViscosity(gp, atTime);


    answer.resize(2);
    answer.zero();

    FloatArray dirX(3), dirY(3), h(2);
    dirX.zero();
    dirX.at(1) = 1.0;
    dirY.zero();
    dirY.at(2) = 1.0;

    h.at(1) = this->giveCharacteristicLenght(gp, dirX);
    h.at(2) = this->giveCharacteristicLenght(gp, dirY);

    answer.at(1) = 3.0 * h.at(1) * h.at(1) / ( 8.0 * mu );
    answer.at(2) = 3.0 * h.at(2) * h.at(2) / ( 8.0 * mu );

    answer.zero();

    return;
}

void
TR1_2D_PFEM :: computeForceVector(FloatArray &answer, TimeStep *atTime) //F
{
    //copied form tr21stokes.C :: computeLoadVector
    int i, load_number, load_id;
    Load *load;
    bcGeomType ltype;
    FloatArray vec;

    int nLoads = this->boundaryLoadArray.giveSize() / 2;
    answer.resize(6);
    //	answer.resize(15);
    answer.zero();
    for ( i = 1; i <= nLoads; i++ ) {  // For each Neumann boundary condition
        load_number = this->boundaryLoadArray.at(2 * i - 1);
        load_id = this->boundaryLoadArray.at(2 * i);
        load = this->domain->giveLoad(load_number);
        ltype = load->giveBCGeoType();

        if ( ltype == EdgeLoadBGT ) {
            this->computeEdgeBCSubVectorAt(vec, load, load_id, atTime);
            answer.add(vec);
        }
    }

    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( i = 1; i <= nLoads; i++ ) {
        load  = domain->giveLoad( bodyLoadArray.at(i) );
        ltype = load->giveBCGeoType();
        if ( ltype == BodyLoadBGT && load->giveBCValType() == ForceLoadBVT ) {
            this->computeBodyLoadVectorAt(vec, load, atTime);
            answer.add(vec);
        }
    }

    return;
}

//copied from tr21stokes
void
TR1_2D_PFEM :: computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *atTime)
{
    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    GaussPoint *gp;
    FloatArray N, gVector, *lcoords, temparray(6);
    double dA, detJ, rho;

    answer.resize(6);
    answer.zero();

    load->computeComponentArrayAt(gVector, atTime, VM_Total);
    temparray.zero();
    if ( gVector.giveSize() ) {
        for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
            gp = iRule->getIntegrationPoint(k);
            lcoords = gp->giveCoordinates();

            rho = this->giveMaterial()->give( 'd', integrationRulesArray [ 0 ]->getIntegrationPoint(0) );
            //detJ = this->interpolation_quad.giveTransformationJacobian(this->domain, this->dofManArray, *lcoords, 0.0);
            //????? pressure x velocity
            detJ = this->pressureInterpolation.giveTransformationJacobian( * lcoords, FEIElementGeometryWrapper(this) );
            dA = detJ * gp->giveWeight();
            //????? pressure x velocity
            this->pressureInterpolation.evalN( N, * lcoords, FEIElementGeometryWrapper(this) );
            for ( int j = 0; j < 3; j++ ) {
                answer(2 * j)     += N(j) * rho * gVector(0) * dA;
                answer(2 * j + 1) += N(j) * rho * gVector(1) * dA;
                //    temparray(2 * j)     += N(j) * rho * gVector(0) * dA;
                //                temparray(2 * j + 1) += N(j) * rho * gVector(1) * dA;
            }
        }
    }

    //    answer.resize(6);
    //    answer.zero();
    //    answer.assemble( temparray, this->ordering );
}


// NOT IN USE
//copied from tr21stokes
void
TR1_2D_PFEM :: computeEdgeBCSubVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *atTime)
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

        for ( int i = 0; i < iRule.giveNumberOfIntegrationPoints(); i++ ) {
            gp = iRule.getIntegrationPoint(i);
            FloatArray *lcoords = gp->giveCoordinates();

            //this->interpolation_quad.edgeEvalN (N, *(gp->giveCoordinates()), 0.0);
            //????? pressure x velocity
            this->pressureInterpolation.edgeEvalN( N, iEdge, * lcoords, FEIElementGeometryWrapper(this) );
            double dS = gp->giveWeight() * this->pressureInterpolation.edgeGiveTransformationJacobian( iEdge, * lcoords, FEIElementGeometryWrapper(this) );

            if ( boundaryLoad->giveFormulationType() == Load :: FT_Entity ) {         // Edge load in xi-eta system
                boundaryLoad->computeValueAt(t, atTime, * lcoords, VM_Total);
            } else {   // Edge load in x-y system
                FloatArray gcoords;
                //????? pressure x velocity
                this->pressureInterpolation.edgeLocal2global( gcoords, iEdge, * lcoords, FEIElementGeometryWrapper(this) );
                boundaryLoad->computeValueAt(t, atTime, gcoords, VM_Total);
            }

            // Reshape the vector
            for ( int j = 0; j < 3; j++ ) {
                f(2 * j)   += N(j) * t(0) * dS;
                f(2 * j + 1) += N(j) * t(1) * dS;
            }
        }

        answer.assemble(f, this->edge_ordering [ iEdge - 1 ]);
    } else {
        OOFEM_ERROR("Strange boundary condition type");
    }
}


// NOT IN USE
void
TR1_2D_PFEM :: computePFEMSubstitutionMatrix(FloatMatrix &answer, TimeStep *atTime) //S
{
    answer.zero();
    FloatMatrix Q, M_hat, L_tau, temp, temp2;
    this->computeStabilizationGradientMatrix(Q, atTime);
    this->computeStabilizationMassMatrix(M_hat, atTime);
    this->computeStabilizedLaplacianMatrix(L_tau, atTime);
    //this->computeStabilizedLaplacianMatrix(L_tau, gp, atTime);

    temp.beInverseOf(M_hat);
    temp2.beProductOf(Q, temp);
    answer.beProductTOf(temp2, Q);

    answer.times(-1.0);

    answer.add(L_tau);

    return;
}

// NOT IN USE
void
TR1_2D_PFEM :: computeConsistentMassMtrx(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(6, 6);
    answer.zero();
    //double rho = this->giveMaterial()->give('d');
    double rho = this->giveMaterial()->give( 'd', integrationRulesArray [ 0 ]->getIntegrationPoint(0) );

    double ar6 = rho * area / 6.0;
    double ar12 = rho * area / 12.0;

    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = ar6;
    answer.at(4, 4) = answer.at(5, 5) = answer.at(6, 6) = ar6;

    answer.at(1, 3) = answer.at(1, 5) = ar12;
    answer.at(3, 1) = answer.at(3, 5) = ar12;
    answer.at(5, 1) = answer.at(5, 3) = ar12;

    answer.at(2, 4) = answer.at(2, 6) = ar12;
    answer.at(4, 2) = answer.at(4, 6) = ar12;
    answer.at(6, 2) = answer.at(6, 4) = ar12;
}


// NOT IN USE
void
TR1_2D_PFEM :: computeDiagonalMassMtrx(FloatArray &answer, TimeStep *atTime)
{
    int i;
    answer.resize(6);
    answer.zero();

    //double rho = this->giveMaterial()->give('d');
    double rho = this->giveMaterial()->give( 'd', integrationRulesArray [ 0 ]->getIntegrationPoint(0) );
    double mm = rho * this->area / 3.0;
    for ( i = 1; i <= 6; i++ ) {
        answer.at(i) = mm;
    }
}


void
TR1_2D_PFEM :: computeDiagonalMassMtrx(FloatMatrix &answer, TimeStep *atTime)
{
    int i;
    answer.resize(6, 6);
    answer.zero();

    //double rho = this->giveMaterial()->give('d');
    double rho = this->giveMaterial()->give( 'd', integrationRulesArray [ 0 ]->getIntegrationPoint(0) );
    double mm = rho * this->area / 3.0;
    for ( i = 1; i <= 6; i++ ) {
        answer.at(i, i) = mm;
    }
}

// NOT IN USE
void
TR1_2D_PFEM :: computeInvertedStabilizationDiagonalMassMtrx(FloatMatrix &answer, TimeStep *atTime)
{
    int i;
    answer.resize(6, 6);
    answer.zero();
    FloatArray tau;
    this->computeStabilizationParameters(tau, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);

    double mm = 3.0 / this->area;

    for ( i = 1; i <= 3; i++ ) {
        answer.at(2 * i - 1, 2 * i - 1) = mm / tau.at(1);
        answer.at(2 * i, 2 * i) = mm / tau.at(2);
    }
}

// NOT IN USE
int
TR1_2D_PFEM :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    double l1, l2, l3;

    l1 = lcoords.at(1);
    l2 = lcoords.at(2);
    l3 = 1.0 - l1 - l2;

    answer.resize(2);
    answer.at(1) = l1 * this->giveNode(1)->giveCoordinate(1) + l2 *this->giveNode(2)->giveCoordinate(1) +
    l3 *this->giveNode(3)->giveCoordinate(1);
    answer.at(2) = l1 * this->giveNode(1)->giveCoordinate(2) + l2 *this->giveNode(2)->giveCoordinate(2) +
    l3 *this->giveNode(3)->giveCoordinate(2);

    return 1;
}


// NOT IN USE
Interface *
TR1_2D_PFEM :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        //        return ( ZZNodalRecoveryModelInterface * ) this;
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        //        return ( NodalAveragingRecoveryModelInterface * ) this;
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        //        return ( SPRNodalRecoveryModelInterface * ) this;
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        //        return ( SpatialLocalizerInterface * ) this;
    } else if ( interface == EIPrimaryFieldInterfaceType ) {
        //        return ( EIPrimaryFieldInterface * ) this;
    }
    //<RESTRICTED_SECTION>
    else if ( interface == LEPlicElementInterfaceType ) {
        //        return ( LEPlicElementInterface * ) this;
    }

    //</RESTRICTED_SECTION>
    return NULL;
}


#define POINT_TOL 1.e-3

// NOT IN USE
bool
TR1_2D_PFEM :: computeLocalCoordinates(FloatArray &answer, const FloatArray &coords)
{
    Node *node1, *node2, *node3;
    double area, x1, x2, x3, y1, y2, y3;

    node1 = this->giveNode(1);
    node2 = this->giveNode(2);
    node3 = this->giveNode(3);

    x1 = node1->giveCoordinate(1);
    x2 = node2->giveCoordinate(1);
    x3 = node3->giveCoordinate(1);

    y1 = node1->giveCoordinate(2);
    y2 = node2->giveCoordinate(2);
    y3 = node3->giveCoordinate(2);

    area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );

    answer.resize(3);

    answer.at(1) = ( ( x2 * y3 - x3 * y2 ) + ( y2 - y3 ) * coords.at(1) + ( x3 - x2 ) * coords.at(2) ) / 2. / area;
    answer.at(2) = ( ( x3 * y1 - x1 * y3 ) + ( y3 - y1 ) * coords.at(1) + ( x1 - x3 ) * coords.at(2) ) / 2. / area;
    answer.at(3) = ( ( x1 * y2 - x2 * y1 ) + ( y1 - y2 ) * coords.at(1) + ( x2 - x1 ) * coords.at(2) ) / 2. / area;


    for ( int i = 1; i <= 3; i++ ) {
        if ( answer.at(i) < ( 0. - POINT_TOL ) ) {
            return false;
        }

        if ( answer.at(i) > ( 1. + POINT_TOL ) ) {
            return false;
        }
    }

    return true;
}


void
TR1_2D_PFEM :: computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    /* one should call material driver instead */
    FloatArray u(6), eps(3);
    answer.resize(3);


    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);

    eps.at(1) = ( b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5) );
    eps.at(2) = ( c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6) );
    eps.at(3) = ( b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6) + c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5) );
    ( ( FluidDynamicMaterial * ) this->giveMaterial() )->computeDeviatoricStressVector(answer, gp, eps, tStep);
}

int
TR1_2D_PFEM :: checkConsistency()
{
    Node *node1, *node2, *node3;
    double x1, x2, x3, y1, y2, y3;

    node1 = giveNode(1);
    node2 = giveNode(2);
    node3 = giveNode(3);

    // init geometry data
    x1 = node1->giveCoordinate(1);
    x2 = node2->giveCoordinate(1);
    x3 = node3->giveCoordinate(1);

    y1 = node1->giveCoordinate(2);
    y2 = node2->giveCoordinate(2);
    y3 = node3->giveCoordinate(2);

    this->area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );

    b [ 0 ] = ( y2 - y3 ) / ( 2. * area );
    c [ 0 ] = ( x3 - x2 ) / ( 2. * area );
    b [ 1 ] = ( y3 - y1 ) / ( 2. * area );
    c [ 1 ] = ( x1 - x3 ) / ( 2. * area );
    b [ 2 ] = ( y1 - y2 ) / ( 2. * area );
    c [ 2 ] = ( x2 - x1 ) / ( 2. * area );

    return PFEMElement2d :: checkConsistency();
}

double
TR1_2D_PFEM :: computeCriticalTimeStep(TimeStep *tStep)
{
    double deltaT = tStep->giveTimeIncrement();
#if 1
    FloatArray u;
    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
    double vn1 = sqrt( u.at(1) * u.at(1) + u.at(2) * u.at(2) );
    double vn2 = sqrt( u.at(3) * u.at(3) + u.at(4) * u.at(4) );
    double vn3 = sqrt( u.at(5) * u.at(5) + u.at(6) * u.at(6) );

    Node *node1, *node2, *node3;

    node1 = giveNode(1);
    node2 = giveNode(2);
    node3 = giveNode(3);

    double x1, x2, x3, y1, y2, y3;
    double l12, l23, l31;

    x1 = node1->giveCoordinate(1);
    x2 = node2->giveCoordinate(1);
    x3 = node3->giveCoordinate(1);

    y1 = node1->giveCoordinate(2);
    y2 = node2->giveCoordinate(2);
    y3 = node3->giveCoordinate(2);


    l12 = sqrt( ( x2 - x1 ) * ( x2 - x1 ) + ( y2 - y1 ) * ( y2 - y1 ) );
    l23 = sqrt( ( x3 - x2 ) * ( x3 - x2 ) + ( y3 - y2 ) * ( y3 - y2 ) );
    l31 = sqrt( ( x1 - x3 ) * ( x1 - x3 ) + ( y1 - y3 ) * ( y1 - y3 ) );

    double dt1, dt2, dt3;

    if ( vn1 < 1.e-6 ) {
        dt1 = deltaT;
    } else {
        dt1 = min(l12, l31) / vn1;
    }

    if ( vn2 < 1.e-6 ) {
        dt2 = deltaT;
    } else {
        dt2 = min(l12, l23) / vn2;
    }

    if ( vn3 < 1.e-6 ) {
        dt3 = deltaT;
    } else {
        dt3 = min(l23, l31) / vn3;
    }

    double dt_min = min( dt1, min(dt2, dt3) );

    return dt_min;

#else
    FloatArray u;
    double dt1, dt2, dt;
    //double Re = domain->giveEngngModel()->giveUnknownComponent(ReynoldsNumber, VM_Unknown, tStep, domain, NULL);
    double Re = 1.0;

    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);

    double vn1 = sqrt( u.at(1) * u.at(1) + u.at(2) * u.at(2) );
    double vn2 = sqrt( u.at(3) * u.at(3) + u.at(4) * u.at(4) );
    double vn3 = sqrt( u.at(5) * u.at(5) + u.at(6) * u.at(6) );
    double veln = max( vn1, max(vn2, vn3) );

    double l1 = 1.0 / ( sqrt(b [ 0 ] * b [ 0 ] + c [ 0 ] * c [ 0 ]) );
    double l2 = 1.0 / ( sqrt(b [ 1 ] * b [ 1 ] + c [ 1 ] * c [ 1 ]) );
    double l3 = 1.0 / ( sqrt(b [ 2 ] * b [ 2 ] + c [ 2 ] * c [ 2 ]) );

    double ln = min( l1, min(l2, l3) );

    // viscous limit
    dt2 = 0.5 * ln * ln * Re;
    if ( veln != 0.0 ) {
        dt1 = ln / veln;
        dt = dt1 * dt2 / ( dt1 + dt2 );
    } else {
        dt = dt2;
    }

    return dt;

#endif
}




void
TR1_2D_PFEM :: updateYourself(TimeStep *tStep)
{
    PFEMElement :: updateYourself(tStep);
}

// NOT IN USE
int
TR1_2D_PFEM :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    //<RESTRICTED_SECTION>
    if ( type == IST_VOFFraction ) {
        answer.resize(1);
        //        answer.at(1) = this->giveTempVolumeFraction();
        return 1;
    } else
    //</RESTRICTED_SECTION>
    if ( type == IST_Density ) {
        answer.resize(1);
        //answer.at(1) = this->giveMaterial()->giveCharacteristicValue(MRM_Density, aGaussPoint, atTime);
        answer.at(1) = this->giveMaterial()->give( 'd', integrationRulesArray [ 0 ]->getIntegrationPoint(0) );
        return 1;
    } else {
        return PFEMElement :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}

// int
// TR1_2D_PFEM :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type)
// {
//     if ( ( type == IST_VOFFraction ) || ( type == IST_Density ) ) {
//         answer.resize(1);
//         answer.at(1) = 1;
//         return 1;
//     } else {
//         return PFEMElement :: giveIntVarCompFullIndx(answer, type);
//     }
// }


// NOT IN USE
//InternalStateValueType
//TR1_2D_PFEM :: giveIPValueType(InternalStateType type)
//{
//    if ( ( type == IST_VOFFraction ) || ( type == IST_Density ) || ( type == IST_Pressure) ) {
//        return ISVT_SCALAR;
//	} else if (type == IST_Velocity) {
//		return ISVT_VECTOR;
//    } else {
//        return PFEMElement :: giveIPValueType(type);
//    }
//}

void
TR1_2D_PFEM :: printOutputAt(FILE *file, TimeStep *stepN)
// Performs end-of-step operations.
{
    PFEMElement :: printOutputAt(file, stepN);
}



contextIOResultType TR1_2D_PFEM :: saveContext(DataStream *stream, ContextMode mode, void *obj)
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



contextIOResultType TR1_2D_PFEM :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
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

double
TR1_2D_PFEM :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
    double determinant, weight, volume;

    determinant = fabs( this->velocityInterpolation.giveTransformationJacobian( * aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this) ) );


    weight      = aGaussPoint->giveWeight();
    volume      = determinant * weight;

    return volume;
}

#ifdef __OOFEG
int
TR1_2D_PFEM :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                       int node, TimeStep *atTime)
{
    //<RESTRICTED_SECTION>
    if ( type == IST_VOFFraction ) {
        answer.resize(1);
        //        answer.at(1) = this->giveTempVolumeFraction();
        return 1;
    } else
    //</RESTRICTED_SECTION>
    if ( type == IST_Density ) {
        answer.resize(1);
        answer.at(1) = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);
        return 1;
    } else {
        return PFEMElement :: giveInternalStateAtNode(answer, type, mode, node, atTime);
    }
}

void
TR1_2D_PFEM :: drawRawGeometry(oofegGraphicContext &gc)
{
    WCRec p [ 3 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(TRUE);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 2 ].z = 0.;

    go =  CreateTriangle3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}

void TR1_2D_PFEM :: drawScalar(oofegGraphicContext &context)
{
    int i, indx, result = 0;
    WCRec p [ 3 ];
    GraphicObj *tr;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    FloatArray v1, v2, v3;
    double s [ 3 ];
    IntArray map;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    if ( context.giveIntVarMode() == ISM_recovered ) {
        result += this->giveInternalStateAtNode(v1, context.giveIntVarType(), context.giveIntVarMode(), 1, tStep);
        result += this->giveInternalStateAtNode(v2, context.giveIntVarType(), context.giveIntVarMode(), 2, tStep);
        result += this->giveInternalStateAtNode(v3, context.giveIntVarType(), context.giveIntVarMode(), 3, tStep);
    } else if ( context.giveIntVarMode() == ISM_local ) {
        GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
        result += giveIPValue(v1, gp, context.giveIntVarType(), tStep);
        v2 = v1;
        v3 = v1;
        result *= 3;
    }

    if ( result != 3 ) {
        return;
    }

    this->giveIntVarCompFullIndx( map, context.giveIntVarType() );

    if ( ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
        return;
    }

    s [ 0 ] = v1.at(indx);
    s [ 1 ] = v2.at(indx);
    s [ 2 ] = v3.at(indx);

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);

    if ( context.getScalarAlgo() == SA_ISO_SURF ) {
        for ( i = 0; i < 3; i++ ) {
            p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
            p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
            p [ i ].z = 0.;
        }

        //EASValsSetColor(gc.getYieldPlotColor(ratio));
        context.updateFringeTableMinMax(s, 3);
        tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
        EGWithMaskChangeAttributes(LAYER_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    } else if ( ( context.getScalarAlgo() == SA_ZPROFILE ) || ( context.getScalarAlgo() == SA_COLORZPROFILE ) ) {
        double landScale = context.getLandScale();

        for ( i = 0; i < 3; i++ ) {
            p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
            p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
            p [ i ].z = s [ i ] * landScale;
        }

        if ( context.getScalarAlgo() == SA_ZPROFILE ) {
            EASValsSetColor( context.getDeformedElementColor() );
            EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
            EASValsSetFillStyle(FILL_SOLID);
            tr =  CreateTriangle3D(p);
            EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | FILL_MASK | LAYER_MASK, tr);
        } else {
            context.updateFringeTableMinMax(s, 3);
            EASValsSetFillStyle(FILL_SOLID);
            tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
            EGWithMaskChangeAttributes(FILL_MASK | LAYER_MASK, tr);
        }

        EMAddGraphicsToModel(ESIModel(), tr);
    }
}



#endif
} // end namespace oofem
