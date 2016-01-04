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
//last edit: 07/02/2013 by Jan Novak

#include "../sm/Elements/htselement.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "node.h"
#include "mathfem.h"
#include "boundaryload.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(HTSelement);

HTSelement :: HTSelement(int n, Domain *aDomain) : StructuralElement(n, aDomain)
    // Constructor.
{
    this->lambda = 8;
    this->mu = 8;
    this->numberOfStressDofs = 12;
    // numberOfEdges = 3;
    numberOfEdges = 4;
    numberOfDofMans  = 2 * numberOfEdges;
    numberOfDofs = 4 * numberOfEdges;
}


void
HTSelement :: giveDofManDofIDMask(int inode, IntArray &answer) const

{
    if ( inode <= numberOfEdges ) {
        answer.clear();
    } else {
        answer = {D_u_edge_const, D_u_edge_lin, D_v_edge_const, D_v_edge_lin};
    }
}


IRResultType
HTSelement :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // Required by IR_GIVE_FIELD macro
    result = StructuralElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }
    numberOfGaussPoints = 8;
    //IR_GIVE_FIELD(ir, numberOfEdges, _IFT_HTSelement_numberOfEdges, "numberOfEdges");
    //numberOfEdges = 3;

    this->computeCenterOfGravity();
    return IRRT_OK;
}

void
HTSelement :: computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(numberOfEdges);

        for ( int i = 0; i < numberOfEdges; i++ ) {
            integrationRulesArray [ i ].reset( new GaussIntegrationRule(i + 1, this, 1, 100) );
            integrationRulesArray [ i ]->SetUpPointsOnLine(numberOfGaussPoints, _1dMat);
        }
    }
}

void
HTSelement :: computeCenterOfGravity()
{
    double area = 0, sX = 0, sY = 0;
    double aX, aY, bX, bY;
    Node *nodeA, *nodeB;

    for ( int i = 0; i < numberOfEdges; i++ ) {
        nodeA   = this->giveSideNode(i + 1, 1);
        nodeB   = this->giveSideNode(i + 1, 2);

        aX = nodeA->giveCoordinate(1);
        aY = nodeA->giveCoordinate(2);
        bX = nodeB->giveCoordinate(1);
        bY = nodeB->giveCoordinate(2);
        area +=  ( bX - aX ) * ( bY + aY ) / 2;

        sY += ( ( aX * aX ) + ( bX - aX ) * aX + ( bX - aX ) * ( bX - aX ) / 3 ) * ( bY - aY ) / 2.0;
        sX += ( ( aY * aY ) + ( bY - aY ) * aY + ( bY - aY ) * ( bY - aY ) / 3 ) * ( bX - aX ) / 2.0;
    }

    cgX = -sY / area;
    cgY = sX / area;
}

Node *
HTSelement :: giveSideNode(int elementSideNumber, int nodeNumber)
{
    int firstNodeNumber = elementSideNumber;
    int secondNodeNumber = elementSideNumber + 1;
    if ( secondNodeNumber > numberOfEdges ) {
        secondNodeNumber = 1;
    }
    if ( nodeNumber == 1 ) {
        return this->giveNode(firstNodeNumber);
    } else if ( nodeNumber == 2 ) {
        return this->giveNode(secondNodeNumber);
    } else {
        //error("Only two nodes per side");
        return 0;
    }
}


double
HTSelement :: computeVolumeAroundSide(GaussPoint *gp, int elemSideNumber)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    return 0.5 * this->giveSideLength(elemSideNumber) * gp->giveWeight(); // * this->giveCrossSection()->give(CS_Area);
}

double
HTSelement :: giveSideLength(int sideNumber)
{
    Node *nodeA;
    Node *nodeB;
    double dx, dy;
    nodeA   = this->giveSideNode(sideNumber, 1);
    nodeB   = this->giveSideNode(sideNumber, 2);
    dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
    dy      = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
    return sqrt(dx * dx + dy * dy);
}


void
HTSelement :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                     TimeStep *tStep)
// Computes numerically the stiffness matrix of the receiver.
{
    double dV;
    answer.resize(numberOfDofs, numberOfDofs);
    answer.zero();

    FloatMatrix Fedge, Aedge, F, A, N;
    F.resize(numberOfStressDofs, numberOfStressDofs);
    F.zero();
    A.resize(numberOfStressDofs, numberOfDofs);
    A.zero();
    for ( int i = 0; i < numberOfEdges; i++ ) {
        this->computeOutwardNormalMatrix(N, i + 1);
        for ( GaussPoint *gp: *this->giveIntegrationRule(i) ) {
            dV = this->computeVolumeAroundSide(gp, i + 1);
            this->computeFMatrixAt(Fedge, N, gp, i + 1);
            Fedge.times(dV);
            this->computeAMatrixAt(Aedge, N, gp, i + 1);
            Aedge.times(dV);
            F.add(Fedge);
        }
        for ( int k = 0; k < numberOfStressDofs; k++ ) {
            for ( int l = 0; l < 4; l++ ) {
                A(k, l + i * 4) = Aedge(k, l);
            }
        }
    }

    // kondenzace
    FloatMatrix invF;
    invF.beInverseOf(F);

    FloatMatrix FG;
    FG.beProductOf(invF, A);
    answer.beTProductOf(A, FG);
}


void
HTSelement :: computePuVectorAt(FloatArray &answer, FloatMatrix N, FloatArray u, GaussPoint *gp, int sideNumber) //\int{(NSv)^T*u}
{
    FloatMatrix Sv, NSv;
    //answer.resize() ??? otazka zda bude nutne ji preskejlovat, asi alo jo vzhledem k tomu jaky element vracim (ctverec, trojuhelnik)
    answer.zero(); //radsi ano
    this->computeSvMatrixAt(Sv, gp, sideNumber); //the quation is whether (NSv)^t Uv will operate on the same gp
    //if so, this matrix should be calculated just once in F matrix
    NSv.beProductOf(N, Sv);
    //  answer.beTProductOf(NSv, u);
} //end of computePuVectorAt


void
HTSelement :: computePsVectorAt(FloatArray &answer, FloatArray t, GaussPoint *gp) //\int{Ugamma^T*t}
{
    FloatMatrix Ugamma;
    //answer.resize() ??? otazka zda bude nutne ji preskejlovat, asi alo jo vzhledem k tomu jaky element vracim (ctverec, trojuhelnik)
    answer.zero(); //radsi ano
    //vytahnout t - zatizeni na hrane
    this->computeUgammaMatrixAt(Ugamma, gp);
    answer.beTProductOf(Ugamma, t);
} //end of computePsVectorAt

void
HTSelement :: computePrescribedDisplacementLoadVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
    //double dV;

    //  FloatArray PuEdge,Pu(numberOfDofs);
    //FloatMatrix N;
    //IntegrationRule* iRule;

    FloatArray u;
    FloatMatrix K;

    this->computeVectorOf(mode, tStep, u);
    if ( u.containsOnlyZeroes() ) {
        answer.clear();
    } else {
        this->computeStiffnessMatrix(K, TangentStiffness, tStep);
        answer.beProductOf(K, u);
        answer.negated();
    }

#if 0
    Pu.zero();
    answer.resize(numberOfDofs);
    for ( int i = 0; i < numberOfEdges; i++ ) {
        this->computeOutwardNormalMatrix(N, i + 1);
        iRule =  this->giveIntegrationRule(i);
        for ( GaussPoint *gp: *iRule ) {
            dV = this->computeVolumeAroundSide(gp, i + 1);
            this->computePuVectorAt(PuEdge, N, u, gp, i + 1);
            PuEdge.times(dV);
            answer.add(PuEdge);
        }
    }
#endif
}


void
HTSelement :: computeEdgeLoadVectorAt(FloatArray &answer, Load *load,
                                      int iEdge, TimeStep *tStep, ValueModeType mode)
{
    double dV;
    BoundaryLoad *edgeLoad = dynamic_cast< BoundaryLoad * >(load);
    FloatArray force, PsEdge, Ps(numberOfDofs);

    Ps.zero();
    answer.resize(numberOfDofs);
    for ( int i = 0; i < numberOfEdges; i++ ) {
        for ( GaussPoint *gp: *this->giveIntegrationRule(i) ) {
            edgeLoad->computeValueAt(force, tStep, gp->giveNaturalCoordinates(), mode);
            dV = this->computeVolumeAroundSide(gp, i + 1);
            this->computePsVectorAt(PsEdge, force, gp);
            PsEdge.times(dV);
            Ps.add(PsEdge);
        }
        for ( int k = 0; k < numberOfDofs; k++ ) {
            answer(k + i * 4) = Ps(k);
        }
    }
}


void
HTSelement :: computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{
    int nLoads;
    bcGeomType ltype;
    Load *load;
    FloatArray helpLoadVector;

    answer.clear();

    // loop over boundary load array
    nLoads = this->giveBoundaryLoadArray()->giveSize() / 2;
    for ( int i = 1; i <= nLoads; i++ ) {
        int n = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
        int id = boundaryLoadArray.at(i * 2);
        load = domain->giveLoad(n);
        ltype = load->giveBCGeoType();
        if ( ltype == EdgeLoadBGT ) {
            this->computeEdgeLoadVectorAt(helpLoadVector, load, id, tStep, mode);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else {
            OOFEM_ERROR("boundary load %d is of unsupported type (%d)", id, ltype);
        }
    }
    this->computePrescribedDisplacementLoadVectorAt(helpLoadVector, tStep, mode);
    answer.add(helpLoadVector);
}
/*Public functions*/

//element stifness matrix


/*Protected functions*/


void
HTSelement :: computeOutwardNormalMatrix(FloatMatrix &answer, int sideNumber)
{
    double Ax =  this->giveSideNode(sideNumber, 1)->giveCoordinate(1);
    double Bx =  this->giveSideNode(sideNumber, 2)->giveCoordinate(1);
    double Ay =  this->giveSideNode(sideNumber, 1)->giveCoordinate(2);
    double By =  this->giveSideNode(sideNumber, 2)->giveCoordinate(2);

    answer.resize(2, 3);
    answer.zero();
    answer.at(1, 1) =  answer.at(2, 3) = ( By - Ay );
    answer.at(1, 3) =  answer.at(2, 2) = -( Bx - Ax );
    double norm = sqrt( answer.at(1, 1) * answer.at(1, 1) + answer.at(2, 2) * answer.at(2, 2) );
    answer.times(1. / norm);
}



//generalized stress degrees of freedom: v
void
HTSelement :: computeFMatrixAt(FloatMatrix &answer, FloatMatrix N, GaussPoint *gp, int sideNumber) //\int{(NSv)^T Uv}
{
    FloatMatrix Uv, Sv, NSv;
    //answer.resize() ??? otazka zda bude nutne ji preskejlovat, asi alo jo vzhledem k tomu jaky element vracim (ctverec, trojuhelnik)
    answer.zero(); //radsi ano

    this->computeUvMatrixAt(Uv, gp, sideNumber);
    this->computeSvMatrixAt(Sv, gp, sideNumber);

    NSv.beProductOf(N, Sv);
    answer.beTProductOf(NSv, Uv);
} //end of computeFMatrixAt


//addtional displacements degrees of freedom: q (also called G matrix)
void
HTSelement :: computeAMatrixAt(FloatMatrix &answer, FloatMatrix N, GaussPoint *gp, int sideNumber) //\int{(NSv)^T Ug}
{
    //answer.resize() ??? otazka zda bude nutne ji preskejlovat, asi alo jo vzhledem k tomu jaky element vracim (ctverec, trojuhelnik)
    answer.zero(); //radsi ano
    FloatMatrix Sv, Ugamma, NSv;

    this->computeUgammaMatrixAt(Ugamma, gp);
    this->computeSvMatrixAt(Sv, gp, sideNumber);
    NSv.beProductOf(N, Sv);
    answer.beTProductOf(NSv, Ugamma);
} //end of computeAMatrixAt


void
HTSelement :: computeUvMatrixAt(FloatMatrix &answer, GaussPoint *gp, int sideNumber)
{
    FloatArray uv;
    FloatMatrix Uv(numberOfStressDofs, 2);
    uv.resize(2);


    double t = ( gp->giveNaturalCoordinate(1) + 1. ) / 2.;

    double Ax =  ( this->giveSideNode(sideNumber, 1)->giveCoordinate(1) ) - cgX;
    double Bx =  this->giveSideNode(sideNumber, 2)->giveCoordinate(1) - cgX;

    double x = ( Bx - Ax ) * t + Ax;

    double Ay =  this->giveSideNode(sideNumber, 1)->giveCoordinate(2) - cgY;
    double By =  this->giveSideNode(sideNumber, 2)->giveCoordinate(2) - cgY;
    double y = ( By - Ay ) * t + Ay;

    //v1 ... vi generalized stress variables the following matrices are associated with
    this->uv1(uv, x, y);
    Uv(0, 0) = uv(0);
    Uv(0, 1) = uv(1);
    //v2
    this->uv3(uv, x, y);
    Uv(1, 0) = uv(0);
    Uv(1, 1) = uv(1);
    //v3
    this->uv4(uv, x, y);
    Uv(2, 0) = uv(0);
    Uv(2, 1) = uv(1);
    //v4
    this->uv5(uv, x, y);
    Uv(3, 0) = uv(0);
    Uv(3, 1) = uv(1);
    //v5
    this->uv6(uv, x, y);
    Uv(4, 0) = uv(0);
    Uv(4, 1) = uv(1);
    //v6
    this->uv7(uv, x, y);
    Uv(5, 0) = uv(0);
    Uv(5, 1) = uv(1);
    //v7
    this->uv8(uv, x, y);
    Uv(6, 0) = uv(0);
    Uv(6, 1) = uv(1);
    //v8
    this->uv9(uv, x, y);
    Uv(7, 0) = uv(0);
    Uv(7, 1) = uv(1);
    //v9
    this->uv10(uv, x, y);
    Uv(8, 0) = uv(0);
    Uv(8, 1) = uv(1);
    //v10
    this->uv11(uv, x, y);
    Uv(9, 0) = uv(0);
    Uv(9, 1) = uv(1);
    //v11
    this->uv12(uv, x, y);
    Uv(10, 0) = uv(0);
    Uv(10, 1) = uv(1);
    //v12
    this->uv25_4(uv, x, y);
    Uv(11, 0) = uv(0);
    Uv(11, 1) = uv(1);
    //transpose Uv in order to conform with Texeira notation
    answer.beTranspositionOf(Uv);
} //end of computeUvMatrixAt

void
HTSelement :: computeSvMatrixAt(FloatMatrix &answer, GaussPoint *gp, int sideNumber)
{
    FloatArray sv;
    FloatMatrix Sv;
    sv.resize(3);
    Sv.resize(numberOfStressDofs, 3);


    double t = ( gp->giveNaturalCoordinate(1) + 1. ) / 2.;

    double Ax =  this->giveSideNode(sideNumber, 1)->giveCoordinate(1);
    double Bx =  this->giveSideNode(sideNumber, 2)->giveCoordinate(1);

    double x = ( Bx - Ax ) * t + Ax - cgX;

    double Ay =  this->giveSideNode(sideNumber, 1)->giveCoordinate(2);
    double By =  this->giveSideNode(sideNumber, 2)->giveCoordinate(2);
    double y = ( By - Ay ) * t + Ay - cgY;




    //v1
    this->sv1(sv, x, y);
    Sv(0, 0) = sv(0);
    Sv(0, 1) = sv(1);
    Sv(0, 2) = sv(2);
    //v2
    this->sv3(sv, x, y);
    Sv(1, 0) = sv(0);
    Sv(1, 1) = sv(1);
    Sv(1, 2) = sv(2);
    //v3
    this->sv4(sv, x, y);
    Sv(2, 0) = sv(0);
    Sv(2, 1) = sv(1);
    Sv(2, 2) = sv(2);
    //v4
    this->sv5(sv, x, y);
    Sv(3, 0) = sv(0);
    Sv(3, 1) = sv(1);
    Sv(3, 2) = sv(2);
    //v5
    this->sv6(sv, x, y);
    Sv(4, 0) = sv(0);
    Sv(4, 1) = sv(1);
    Sv(4, 2) = sv(2);
    //v6
    this->sv7(sv, x, y);
    Sv(5, 0) = sv(0);
    Sv(5, 1) = sv(1);
    Sv(5, 2) = sv(2);
    //v7
    this->sv8(sv, x, y);
    Sv(6, 0) = sv(0);
    Sv(6, 1) = sv(1);
    Sv(6, 2) = sv(2);
    //v8
    this->sv9(sv, x, y);
    Sv(7, 0) = sv(0);
    Sv(7, 1) = sv(1);
    Sv(7, 2) = sv(2);
    //v9
    this->sv10(sv, x, y);
    Sv(8, 0) = sv(0);
    Sv(8, 1) = sv(1);
    Sv(8, 2) = sv(2);
    //v10
    this->sv11(sv, x, y);
    Sv(9, 0) = sv(0);
    Sv(9, 1) = sv(1);
    Sv(9, 2) = sv(2);
    //v11
    this->sv12(sv, x, y);
    Sv(10, 0) = sv(0);
    Sv(10, 1) = sv(1);
    Sv(10, 2) = sv(2);
    //v12
    this->sv25_4(sv, x, y);
    Sv(11, 0) = sv(0);
    Sv(11, 1) = sv(1);
    Sv(11, 2) = sv(2);

    //transpose Sv in order to conform with Texeira notation
    answer.beTranspositionOf(Sv);
} //end of computeSvMatrixAt


void
HTSelement :: computeUgammaMatrixAt(FloatMatrix &answer, GaussPoint *gp)
{
    FloatMatrix Ugamma(4, 2);
    Ugamma.zero();

    //following constatnts should be sent from elsewhere
    // int d = 2, p = 2;  //@p - polynomial order (number of hierarchical functions)
    //@d - task dimension d={2,3}

    //single edge
    //q_const
    Ugamma(0, 0) = Ugamma(1, 1) = u_gammaConst(gp);
    //q_linear
    Ugamma(2, 0) = Ugamma(3, 1) = u_gammaLin(gp);
    answer.beTranspositionOf(Ugamma);
} //end of computeUgammaMatrixAt



/*Private functions*/

//Uv and Sv basis definition:
//@lambda - lame elastic constant
//@mu     - lame elastic constant
void
HTSelement :: uv1(FloatArray &answer, double x, double y)
{
    answer(0) = x;
    answer(1) = 0;
} //end of uv1

void
HTSelement :: uv2(FloatArray &answer, double x, double y)
{
    answer(0) = y;
    answer(1) = 0;
} //end of uv2

void
HTSelement :: uv3(FloatArray &answer, double x, double y)
{
    answer(0) = 0;
    answer(1) = x;
} //end of uv3

void
HTSelement :: uv4(FloatArray &answer, double x, double y)
{
    answer(0) = 0;
    answer(1) = y;
} //end of uv4

void
HTSelement :: sv1(FloatArray &answer, double x, double y)
{
    answer(0) = 2 * mu + lambda;
    answer(1) = lambda;
    answer(2) = 0;
} //end of sv1

void
HTSelement :: sv2(FloatArray &answer, double x, double y)
{
    answer(0) = 0;
    answer(1) = 0;
    answer(2) = mu;
} //end of sv2

void
HTSelement :: sv3(FloatArray &answer, double x, double y)
{
    answer(0) = 0;
    answer(1) = 0;
    answer(2) = mu;
} //end of sv3

void
HTSelement :: sv4(FloatArray &answer, double x, double y)
{
    answer(0) = lambda;
    answer(1) = 2 * mu + lambda;
    answer(2) = 0;
} //end of sv4

void
HTSelement :: uv5(FloatArray &answer, double x, double y)
{
    answer(0) = x * y;
    answer(1) = -0.5 * y * y * ( lambda + mu ) / ( 2 * mu + lambda );
} //end of uv5

void
HTSelement :: uv6(FloatArray &answer, double x, double y)
{
    answer(0) = -( -2 * y * y * mu - y * y * lambda + x * x * mu ) / ( 2 * mu + lambda );
    answer(1) = 0;
} //end of uv6

void
HTSelement :: uv7(FloatArray &answer, double x, double y)
{
    answer(0) = 0;
    answer(1) = ( 2 * x * x * mu + x * x * lambda - y * y * mu ) / ( 2 * mu + lambda );
} //end of uv7

void
HTSelement :: uv8(FloatArray &answer, double x, double y)
{
    answer(0) = -0.5 * x * x * ( lambda + mu ) / ( 2 * mu + lambda );
    answer(1) = x * y;
} //end of uv8

void
HTSelement :: sv5(FloatArray &answer, double x, double y)
{
    answer(0) = mu * y * ( 4 * mu + 3 * lambda ) / ( 2 * mu + lambda );
    answer(1) = -mu * y;
    answer(2) = mu * x;
} //end of sv5

void
HTSelement :: sv6(FloatArray &answer, double x, double y)
{
    answer(0) = -2 * mu * x;
    answer(1) = -2 * lambda * mu * x / ( 2 * mu + lambda );
    answer(2) = 2 * mu * y;
} //end of sv6

void
HTSelement :: sv7(FloatArray &answer, double x, double y)
{
    answer(0) = -2 * lambda * y * mu / ( 2 * mu + lambda );
    answer(1) = -2 * mu * y;
    answer(2) = 2 * mu * x;
} //end of sv7

void
HTSelement :: sv8(FloatArray &answer, double x, double y)
{
    answer(0) = -mu * x;
    answer(1) = mu * x * ( 4 * mu + 3 * lambda ) / ( 2 * mu + lambda );
    answer(2) = mu * y;
} //end of sv8

void
HTSelement :: uv9(FloatArray &answer, double x, double y)
{
    answer(0) = -( 1 / 3.0 ) * x * x * x * mu / ( 2 * mu + lambda ) + x * y * y;
    answer(1) = -( 1 / 3.0 ) * y * y * y * ( lambda + mu ) / ( 2 * mu + lambda );
} //end of uv9

void
HTSelement :: uv10(FloatArray &answer, double x, double y)
{
    answer(0) = -3 * x * x * y * ( 2 * mu + lambda ) / ( 2 * lambda + 3 * mu ) + y * y * y;
    answer(1) = 3 * x * y * y * ( lambda + mu ) / ( 2 * lambda + 3 * mu );
} //end of uv10

void
HTSelement :: uv11(FloatArray &answer, double x, double y)
{
    answer(0) = -3 * x * x * y * ( -mu - lambda ) / ( 2 * lambda + 3 * mu );
    answer(1) = x * x * x + 3 * x * y * y * ( -lambda - 2 * mu ) / ( 2 * lambda + 3 * mu );
} //end of uv11

void
HTSelement :: uv12(FloatArray &answer, double x, double y)
{
    answer(0) = -( 1 / 3.0 ) * x * x * x * ( lambda + mu ) / ( 2 * mu + lambda );
    answer(1) = x * x * y - ( 1 / 3.0 ) * y * y * y * mu / ( 2 * mu + lambda );
} //end of uv12

void
HTSelement :: sv9(FloatArray &answer, double x, double y)
{
    answer(0) = -mu * ( -4 * mu * y * y - 3 * lambda * y * y + 2 * x * x * mu + lambda * x * x ) / ( 2 * mu + lambda );
    answer(1) = -mu * ( lambda * y * y + 2 * mu * y * y + lambda * x * x ) / ( 2 * mu + lambda );
    answer(2) = 2 * mu * x * y;
} //end of sv9

void
HTSelement :: sv10(FloatArray &answer, double x, double y)
{
    answer(0) = -6 * x * y * mu * ( 4 * mu + 3 * lambda ) / ( 2 * lambda + 3 * mu );
    answer(1) = 6 * mu * x * y * ( 2 * mu + lambda ) / ( 2 * lambda + 3 * mu );
    answer(2) = -3 * mu * ( 2 * x * x * mu + x * x * lambda - 3 * lambda * y * y - 4 * mu * y * y ) / ( 2 * lambda + 3 * mu );
} //end of sv10

void
HTSelement :: sv11(FloatArray &answer, double x, double y)
{
    answer(0) = 6 * mu * x * y * ( 2 * mu + lambda ) / ( 2 * lambda + 3 * mu );
    answer(1) = -6 * x * y * mu * ( 4 * mu + 3 * lambda ) / ( 2 * lambda + 3 * mu );
    answer(2) = 3 * mu * ( -2 * mu * y * y - lambda * y * y + 3 * x * x * lambda + 4 * x * x * mu ) / ( 2 * lambda + 3 * mu );
} //end of sv11

void
HTSelement :: sv12(FloatArray &answer, double x, double y)
{
    answer(0) = -mu * ( x * x * lambda + 2 * x * x * mu + lambda * y * y ) / ( 2 * mu + lambda );
    answer(1) = mu * ( 4 * x * x * mu + 3 * x * x * lambda - 2 * mu * y * y - lambda * y * y ) / ( 2 * mu + lambda );
    answer(2) = 2 * mu * x * y;
} //end of sv12

#if 0
void
HTSelement :: sv13(FloatArray &answer, double x, double y)
{
    answer(0) = ( 6 * x * x * y * y * mu + x * x * x * x * lambda - 2 * y * y * y * y * mu - y * y * y * y * lambda ) / lambda;
    answer(1) = -4 * x * x * x * y * ( lambda + mu ) / lambda;
} //end of sv13

void
HTSelement :: uv13(FloatArray &answer, double x, double y)
{
    answer(0) = x * y * ( y * y * mu + x * x * lambda ) / lambda;
    answer(1) = 0.5 * x * x * ( lambda + mu ) * ( x * x - 3 * y * y ) / lambda;
} //end of uv13

void
HTSelement :: uv15(FloatArray &answer, double x, double y)
{
    answer(0) = -0.5 * y * y * ( lambda + mu ) * ( 3 * x * x - y * y ) / lambda;
    answer(1) = x * y * ( x * x * mu + y * y * lambda ) / lambda;
} //end of uv15

void
HTSelement :: uv16(FloatArray &answer, double x, double y)
{
    answer(0) = -4 * x * y * y * y * ( lambda + mu ) / lambda;
    answer(1) = -( -6 * x * x * y * y * mu - y * y * y * y * lambda + x * x * x * x * lambda + 2 * x * x * x * x * mu ) / lambda;
} //end of uv16

void
HTSelement :: sv13(FloatArray &answer, double x, double y)
{
    answer(0) = 4 * mu * x * ( 6 * y * y * mu + x * x * lambda + 3 * y * y * lambda ) / lambda;
    answer(1) = -4 * mu * x * ( 3 * x * x * lambda + 2 * x * x * mu - 3 * lambda * y * y ) / lambda;
    answer(2) = -4 * mu * y * ( 2 * y * y * mu + y * y * lambda + 3 * x * x * lambda ) / lambda;
} //end of sv13

void
HTSelement :: sv14(FloatArray &answer, double x, double y)
{
    answer(0) = mu * y * ( 2 * y * y * mu + 3 * x * x * lambda + lambda * y * y ) / lambda;
    answer(1) = -mu * y * ( 9 * x * x * lambda + 6 * x * x * mu - lambda * y * y ) / lambda;
    answer(2) = mu * x * ( 3 * x * x * lambda - 3 * y * y * lambda + 2 * x * x * mu ) / lambda;
} //end of sv14

void
HTSelement :: sv15(FloatArray &answer, double x, double y)
{
    answer(0) = mu * x * ( -9 * y * y * lambda - 6 * y * y * mu + lambda * x * x ) / lambda;
    answer(1) = mu * x * ( 3 * y * y * lambda + 2 * x * x * mu + lambda * x * x ) / lambda;
    answer(2) = -mu * y * ( -3 * y * y * lambda - 2 * y * y * mu + 3 * x * x * lambda ) / lambda;
} //end of sv15

void
HTSelement :: sv16(FloatArray &answer, double x, double y)
{
    answer(0) = 4 * y * mu * ( -3 * y * y * lambda - 2 * y * y * mu + 3 * lambda * x * x ) / lambda;
    answer(1) = 4 * mu * y * ( y * y * lambda + 6 * x * x * mu + 3 * lambda * x * x ) / lambda;
    answer(2) = -4 * mu * x * ( x * x * lambda + 2 * x * x * mu + 3 * y * y * lambda ) / lambda;
} //end of sv16

void
HTSelement :: uv17(FloatArray &answer, double x, double y)
{
    answer(0) = x * ( 2 * x * x * x * x * lambda + x * x * x * x * mu - 5 * y * y * y * y * mu - 10 * x * x * y * y * lambda ) / ( 2 * lambda + mu );
    answer(1) = -10 * x * x * y * ( lambda + mu ) * ( x * x - y * y ) / ( 2 * lambda + mu );
} //end of uv17

void
HTSelement :: uv18(FloatArray &answer, double x, double y)
{
    answer(0) = ( 1 / 5.0 ) * y * ( 10 * x * x * y * y * mu - y * y * y * y * lambda - 2 * y * y * y * y * mu + 5 * x * x * x * x * lambda ) / lambda;
    answer(1) = ( 2 / 5.0 ) * x * x * x * ( lambda + mu ) * ( -5 * y * y + x * x ) / lambda;
} //end of uv18

void
HTSelement :: uv19(FloatArray &answer, double x, double y)
{
    answer(0) = -( 2 / 5.0 ) * y * y * y * ( lambda + mu ) * ( 5 * x * x - y * y ) / lambda;
    answer(1) = -( 1 / 5.0 ) * x * ( x * x * x * x * lambda + 2 * x * x * x * x * mu - 5 * y * y * y * y * lambda - 10 * x * x * y * y * mu ) / lambda;
} //end of uv19

void
HTSelement :: uv20(FloatArray &answer, double x, double y)
{
    answer(0) = 10 * x * y * y * ( lambda + mu ) * ( -y * y + x * x ) / ( 2 * lambda + mu );
    answer(1) = -y * ( 5 * x * x * x * x * mu - 2 * y * y * y * y * lambda - y * y * y * y * mu + 10 * x * x * y * y * lambda ) / ( 2 * lambda + mu );
} //end of uv20

void
HTSelement :: sv17(FloatArray &answer, double x, double y)
{
    answer(0) = 5 * mu * ( -6 * x * x * y * y * lambda + 3 * x * x * x * x * lambda + 2 * x * x * x * x * mu - 2 * y * y * y * y * mu - lambda * y * y * y * y ) / ( 2 * lambda + mu );
    answer(1) = -5 * mu * ( 5 * x * x * x * x * lambda - 18 * x * x * y * y * lambda + 4 * x * x * x * x * mu - 12 * x * x * y * y * mu + lambda * y * y * y * y ) / ( 2 * lambda + mu );
    answer(2) = -20 * mu * x * y * ( 3 * x * x * lambda + 2 * x * x * mu - y * y * lambda ) / ( 2 * lambda + mu );
} //end of sv17

void
HTSelement :: sv18(FloatArray &answer, double x, double y)
{
    answer(0) = 4 * mu * x * y * ( 2 * y * y * mu + x * x * lambda + lambda * y * y ) / lambda;
    answer(1) = -4 * mu * x * y * ( 3 * x * x * lambda + 2 * x * x * mu - lambda * y * y ) / lambda;
    answer(2) = mu * ( 3 * x * x * x * x * lambda + 2 * x * x * x * x * mu - y * y * y * y * lambda - 2 * y * y * y * y * mu - 6 * x * x * y * y * lambda ) / lambda;
} //end of sv18

void
HTSelement :: sv19(FloatArray &answer, double x, double y)
{
    answer(0) = 4 * mu * x * y * ( -3 * y * y * lambda - 2 * y * y * mu + lambda * x * x ) / lambda;
    answer(1) = 4 * mu * x * y * ( 2 * x * x * mu + y * y * lambda + lambda * x * x ) / lambda;
    answer(2) = -mu * ( 6 * x * x * y * y * lambda + x * x * x * x * lambda + 2 * x * x * x * x * mu - 3 * y * y * y * y * lambda - 2 * y * y * y * y * mu ) / lambda;
} //end of sv19

void
HTSelement :: sv20(FloatArray &answer, double x, double y)
{
    answer(0) = -5 * mu * ( -18 * x * x * y * y * lambda + 5 * y * y * y * y * lambda - 12 * x * x * y * y * mu + 4 * y * y * y * y * mu + lambda * x * x * x * x ) / ( 2 * lambda + mu );
    answer(1) = -5 * mu * ( 6 * x * x * y * y * lambda - 3 * y * y * y * y * lambda - 2 * y * y * y * y * mu + 2 * x * x * x * x * mu + lambda * x * x * x * x ) / ( 2 * lambda + mu );
    answer(2) = 20 * mu * x * y * ( -3 * y * y * lambda + x * x * lambda - 2 * y * y * mu ) / ( 2 * lambda + mu );
} //end of sv20

void
HTSelement :: uv21(FloatArray &answer, double x, double y)
{
    answer(0) = ( -15 * x * x * y * y * y * y * mu - 15 * x * x * x * x * y * y * lambda + 2 * x * x * x * x * x * x * lambda + x * x * x * x * x * x * mu + y * y * y * y * y * y * lambda + \
                  2 * y * y * y * y * y * y * mu ) / ( 2 * lambda + mu );
    answer(1) = -4 * x * x * x * y * ( lambda + mu ) * ( -5 * y * y + 3 * x * x ) / ( 2 * lambda + mu );
} //end of uv21

void
HTSelement :: uv22(FloatArray &answer, double x, double y)
{
    answer(0) = ( 1 / 3.0 ) * x * y * ( -10 * x * x * y * y * lambda - 3 * y * y * y * y * mu + 6 * x * x * x * x * lambda + 3 * x * x * x * x * mu ) / ( 2 * lambda + mu );
    answer(1) = 0.5 * x * x * ( lambda + mu ) * ( 5 * y * y * y * y + x * x * x * x - 10 * x * x * y * y ) / ( 2 * lambda + mu );
} //end of uv22

void
HTSelement :: uv23(FloatArray &answer, double x, double y)
{
    answer(0) = 0.5 * y * y * ( lambda + mu ) * ( y * y * y * y + 5 * x * x * x * x - 10 * x * x * y * y ) / ( 2 * lambda + mu );
    answer(1) = -( 1 / 3.0 ) * x * y * ( -6 * y * y * y * y * lambda - 3 * y * y * y * y * mu + 3 * x * x * x * x * mu + 10 * x * x * y * y * lambda ) / ( 2 * lambda + mu );
} //end of uv23

void
HTSelement :: uv24(FloatArray &answer, double x, double y)
{
    answer(0) = 4 * x * y * y * y * ( lambda + mu ) * ( 5 * x * x - 3 * y * y ) / ( 2 * lambda + mu );
    answer(1) = ( 2 * y * y * y * y * y * y * lambda + y * y * y * y * y * y * mu + x * x * x * x * x * x * lambda + 2 * x * x * x * x * x * x * mu - 15 * x * x * x * x * y * y * mu - \
                  15 * x * x * y * y * y * y * lambda ) / ( 2 * lambda + mu );
} //end of uv24

void
HTSelement :: sv21(FloatArray &answer, double x, double y)
{
    answer(0) = 6 * mu * x * ( 3 * x * x * x * x * lambda + 2 * x * x * x * x * mu - 10 * y * y * y * y * mu - 10 * x * x * y * y * lambda - 5 * lambda * y * y * y * y ) / ( 2 * lambda + mu );
    answer(1) = -6 * mu * x * ( -30 * x * x * y * y * lambda + 5 * x * x * x * x * lambda - 20 * x * x * y * y * mu + 4 * x * x * x * x * mu + 5 * lambda * y * y * y * y ) / ( 2 * lambda + mu );
    answer(2) = -6 * mu * y * ( -y * y * y * y * lambda - 2 * y * y * y * y * mu + 15 * x * x * x * x * lambda - 10 * x * x * y * y * lambda + 10 * x * x * x * x * mu ) / ( 2 * lambda + mu );
} //end of sv21

void
HTSelement :: sv22(FloatArray &answer, double x, double y)
{
    answer(0) = mu * y * ( -10 * x * x * y * y * lambda - 2 * y * y * y * y * mu + 15 * x * x * x * x * lambda + 10 * x * x * x * x * mu - lambda * y * y * y * y ) / ( 2 * lambda + mu );
    answer(1) = -mu * y * ( -30 * x * x * y * y * lambda + 25 * x * x * x * x * lambda - 20 * x * x * y * y * mu + 20 * x * x * x * x * mu + lambda * y * y * y * y ) / ( 2 * lambda + mu );
    answer(2) = mu * x * ( -30 * x * x * y * y * lambda + 5 * x * x * x * x * lambda + 4 * x * x * x * x * mu + 5 * y * y * y * y * lambda - 20 * x * x * y * y * mu ) / ( 2 * lambda + mu );
} //end of sv22

void
HTSelement :: sv23(FloatArray &answer, double x, double y)
{
    answer(0) = -mu * x * ( 25 * y * y * y * y * lambda - 30 * x * x * y * y * lambda + 20 * y * y * y * y * mu - 20 * x * x * y * y * mu + lambda * x * x * x * x ) / ( 2 * lambda + mu );
    answer(1) = -mu * x * ( 2 * x * x * x * x * mu - 15 * y * y * y * y * lambda - 10 * y * y * y * y * mu + 10 * x * x * y * y * lambda + lambda * x * x * x * x ) / ( 2 * lambda + mu );
    answer(2) = mu * y * ( 5 * x * x * x * x * lambda + 5 * y * y * y * y * lambda + 4 * y * y * y * y * mu - 30 * x * x * y * y * lambda - 20 * x * x * y * y * mu ) / ( 2 * lambda + mu );
} //end of sv23

void
HTSelement :: sv24(FloatArray &answer, double x, double y)
{
    answer(0) = -6 * mu * y * ( -30 * x * x * y * y * lambda + 5 * y * y * y * y * lambda - 20 * x * x * y * y * mu + 4 * y * y * y * y * mu + 5 * lambda * x * x * x * x ) / ( 2 * lambda + mu );
    answer(1) = -6 * mu * y * ( 10 * x * x * x * x * mu - 3 * y * y * y * y * lambda - 2 * y * y * y * y * mu + 10 * x * x * y * y * lambda + 5 * lambda * x * x * x * x ) / ( 2 * lambda + mu );
    answer(2) = 6 * mu * x * ( -15 * y * y * y * y * lambda - 10 * y * y * y * y * mu + x * x * x * x * lambda + 2 * x * x * x * x * mu + 10 * x * x * y * y * lambda ) / ( 2 * lambda + mu );
} //end of sv24
#endif

void
HTSelement :: uv25_4(FloatArray &answer, double x, double y)
{
    answer(0) = -5 * x * x * x * x + 30 * x * x * y * y - 5 * y * y * y * y;
    answer(1) = 20 * x * x * x * y - 20 * x * y * y * y;
} //end of uv25_4

void
HTSelement :: sv25_4(FloatArray &answer, double x, double y)
{
    answer(0) = 2 * mu * ( -20 * x * x * x + 60 * x * y * y );
    answer(1) = 2 * mu * ( 20 * x * x * x - 60 * x * y * y );
    answer(2) = mu * ( 120 * x * x * y - 40 * y * y * y );
} //end of sv25_4


//u_gamma functions
double
HTSelement :: u_gammaConst(GaussPoint *gp)
{
    return 1;
} //end of u_gammaConst

double
HTSelement :: u_gammaLin(GaussPoint *gp)
{
    //  double ksi = 1; must be calculated
    return gp->giveNaturalCoordinate(1);
} //end of u_gammaLin
} // end namespace oofem
