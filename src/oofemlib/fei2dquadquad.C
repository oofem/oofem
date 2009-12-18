/* $Header: /home/cvs/bp/oofem/oofemlib/src/fei2dtrlin.C,v 1.1.4.1 2004/04/05 15:19:43 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

#include "fei2dquadquad.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "node.h"
#include "mathfem.h"

namespace oofem {

void
FEI2dQuadQuad :: evalN(FloatArray &answer, const FloatArray &lcoords, double time)
{

  /* Local Node Numbering

         4----7--- 3
         |         |
         |         |
         |         |
         8         6
         |         |
         |         |
         |         |
         1----5----2


   */

  double ksi, eta;
  
  answer.resize(8);

  ksi = lcoords.at(1);
  eta = lcoords.at(2);
  
  answer.at(1) = ( 1. + ksi ) * ( 1. + eta ) * 0.25 * ( ksi + eta - 1. );
  answer.at(2) = ( 1. - ksi ) * ( 1. + eta ) * 0.25 * ( -ksi + eta - 1. );
  answer.at(3) = ( 1. - ksi ) * ( 1. - eta ) * 0.25 * ( -ksi - eta - 1. );
  answer.at(4) = ( 1. + ksi ) * ( 1. - eta ) * 0.25 * ( ksi - eta - 1. );
  answer.at(5) = 0.5 * ( 1. - ksi * ksi ) * ( 1. + eta );
  answer.at(6) = 0.5 * ( 1. - ksi ) * ( 1. - eta * eta );
  answer.at(7) = 0.5 * ( 1. - ksi * ksi ) * ( 1. - eta );
  answer.at(8) = 0.5 * ( 1. + ksi ) * ( 1. - eta * eta );
  
  return;
}

void
FEI2dQuadQuad :: evaldNdx(FloatMatrix &answer, const FloatArray **coords, const FloatArray &lcoords, double time)
{
    answer.resize(8, 2);
    int i;
    FloatMatrix jacobianMatrix(2, 2), inv(2, 2);
    FloatArray nx(8), ny(8);

    this->giveJacobianMatrixAt(jacobianMatrix, coords, lcoords);
    inv.beInverseOf(jacobianMatrix);

    this->giveDerivativeXi (nx, lcoords);
    this->giveDerivativeEta(ny, lcoords);

    for ( i = 1; i <= 8; i++ ) {
        answer.at(i, 1) = nx.at(i) * inv.at(1, 1) + ny.at(i) * inv.at(1, 2);
        answer.at(i, 2) = nx.at(i) * inv.at(2, 1) + ny.at(i) * inv.at(2, 2);
    }
}

void
FEI2dQuadQuad :: local2global(FloatArray &answer, const FloatArray **coords, const FloatArray &lcoords, double time)
{
    int i;
    FloatArray n(8);
    answer.resize(2);
    answer.zero();

    this->evalN(n, lcoords, time);

    for ( i = 1; i <= 8; i++ ) {
        answer.at(1) += n.at(i) * coords [ i - 1 ]->at(xind);
        answer.at(2) += n.at(i) * coords [ i - 1 ]->at(yind);
    }
}

#define POINT_TOL 1.e-3

int
FEI2dQuadQuad :: global2local(FloatArray &answer, const FloatArray **nc, const FloatArray &coords, double time)
{
    FloatArray lc(2);
    FloatArray r(2), n(8), dksi, deta, delta;
    FloatMatrix p(2, 2);
    double x, y;
    int i, nite = 0;

    // setup initial guess
    answer.resize(2);
    answer.at(1) = answer.at(2) = 0.0;

    // apply Newton-Raphson to solve the problem
    do {
        if ( ( ++nite ) > 10 ) {
            //_error ("computeLocalCoords: no convergence after 10 iterations");
            return 0;
        }

        // compute the residual
	this->evalN(n, answer, time);

	r = answer;
        for ( i = 1; i <= 8; i++ ) {
            r.at(1) -= n.at(i) * nc [ i - 1 ]->at(xind);
            r.at(2) -= n.at(i) * nc [ i - 1 ]->at(yind);
        }

        // check for convergence
        if ( sqrt( dotProduct(r, r, 2) ) < 1.e-10 ) {
            break;
        }

        // compute the corrections
        this->giveDerivativeXi (dksi, answer);
        this->giveDerivativeEta(deta, answer);

        p.zero();
        for ( i = 1; i <= 8; i++ ) {
            x = nc [ i - 1 ]->at(xind);
            y = nc [ i - 1 ]->at(yind);

            p.at(1, 1) += dksi.at(i) * x;
            p.at(1, 2) += deta.at(i) * x;
            p.at(2, 1) += dksi.at(i) * y;
            p.at(2, 2) += deta.at(i) * y;
        }

        // solve for corrections
        p.solveForRhs(r, delta);
        // update guess
        answer.add(delta);
    } while ( 1 );

    for ( i = 1; i <= 2; i++ ) {
      if ( fabs(answer.at(i)) > ( 1. + POINT_TOL ) ) {
	return 0;
      }
      
    }

    return 1;
}


double
FEI2dQuadQuad :: giveTransformationJacobian(const FloatArray **coords, const FloatArray &lcoords, double time)
{
    FloatMatrix jacobianMatrix(2, 2);

    this->giveJacobianMatrixAt(jacobianMatrix, coords, lcoords);
    return jacobianMatrix.giveDeterminant();
}


void
FEI2dQuadQuad :: edgeEvalN(FloatArray &answer, const FloatArray &lcoords, double time)
{
  /*
       1-------3-------2
   */

    double n3, ksi = lcoords.at(1);
    answer.resize(3);

    n3 = 1. - ksi * ksi;
    answer.at(1) = ( 1. - ksi ) * 0.5 - 0.5 * n3;
    answer.at(2) = ( 1. + ksi ) * 0.5 - 0.5 * n3;
    answer.at(3) = n3;
}

void
FEI2dQuadQuad :: edgeEvaldNdx(FloatMatrix &answer, int iedge,
			      const FloatArray **coords, const FloatArray &lcoords, double time)
{
  OOFEM_ERROR("FEI2dQuadQuad :: edgeEvaldNdx: not implemented");
}

void
FEI2dQuadQuad :: edgeLocal2global(FloatArray &answer, int iedge,
				  const FloatArray **coords, const FloatArray &lcoords, double time)
{
  IntArray edgeNodes;
  FloatArray n;
  this->computeLocalEdgeMapping(edgeNodes, iedge);
  this->edgeEvalN(n, lcoords, time);
  
  answer.resize(2);
  answer.at(1) = ( n.at(1) * coords [ edgeNodes.at(1) - 1 ]->at(xind) +
		   n.at(2) * coords [ edgeNodes.at(2) - 1 ]->at(xind) +
		   n.at(3) * coords [ edgeNodes.at(3) - 1 ]->at(xind) );
  answer.at(2) = ( n.at(1) * coords [ edgeNodes.at(1) - 1 ]->at(yind) +
		   n.at(2) * coords [ edgeNodes.at(2) - 1 ]->at(yind) +
		   n.at(3) * coords [ edgeNodes.at(3) - 1 ]->at(yind) );
}


void
FEI2dQuadQuad :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
{
    int aNode = 0, bNode = 0, cNode = 0;
    edgeNodes.resize(3);

    if ( iedge == 1 ) { // edge between nodes 1 2
        aNode = 1;
        bNode = 2;
        cNode = 5;
    } else if ( iedge == 2 ) { // edge between nodes 2 3
        aNode = 2;
        bNode = 3;
        cNode = 6;
    } else if ( iedge == 3 ) { // edge between nodes 2 3
        aNode = 3;
        bNode = 4;
        cNode = 7;
    } else if ( iedge == 4 ) { // edge between nodes 3 4
        aNode = 4;
        bNode = 1;
        cNode = 8;
    } else {
        OOFEM_ERROR2("FEI2dQuadQuad :: computeEdgeMapping: wrong egde number (%d)", iedge);
    }

    edgeNodes.at(1) = aNode;
    edgeNodes.at(2) = bNode;
    edgeNodes.at(3) = cNode;
}

double
FEI2dQuadQuad :: edgeGiveTransformationJacobian(int iedge, const FloatArray **coords, const FloatArray &lcoords, double time)
{
    OOFEM_ERROR("FEI2dQuadQuad :: edgeGiveTransformationJacobian: not implemented");
    return 0.0;
}



void
FEI2dQuadQuad :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray **coords, const FloatArray &lcoords)
// Returns the jacobian matrix  J (x,y)/(ksi,eta)  of the receiver.
// Computes it if it does not exist yet.
{
    int i;
    double x, y;
    FloatArray dxi, deta;

    jacobianMatrix.resize(2, 2);
    jacobianMatrix.zero();

    this->giveDerivativeXi(dxi, lcoords);
    this->giveDerivativeEta(deta, lcoords);

    for ( i = 1; i <= 8; i++ ) {
        x = coords [ i - 1 ]->at(xind);
        y = coords [ i - 1 ]->at(yind);

        jacobianMatrix.at(1, 1) += dxi.at(i) * x;
        jacobianMatrix.at(1, 2) += dxi.at(i) * y;
        jacobianMatrix.at(2, 1) += deta.at(i) * x;
        jacobianMatrix.at(2, 2) += deta.at(i) * y;
    }
}


void
FEI2dQuadQuad :: giveDerivativeXi(FloatArray &n, const FloatArray &lc)
{
    double ksi, eta;
    ksi = lc.at(1); eta = lc.at(2);
    n.resize(8);


    n.at(1) =  0.25 * ( 1. + eta ) * ( 2.0 * ksi + eta );
    n.at(2) = -0.25 * ( 1. + eta ) * ( -2.0 * ksi + eta );
    n.at(3) = -0.25 * ( 1. - eta ) * ( -2.0 * ksi - eta );
    n.at(4) =  0.25 * ( 1. - eta ) * ( 2.0 * ksi - eta );
    n.at(5) = -ksi * ( 1. + eta );
    n.at(6) = -0.5 * ( 1. - eta * eta );
    n.at(7) = -ksi * ( 1. - eta );
    n.at(8) =  0.5 * ( 1. - eta * eta );
}

void
FEI2dQuadQuad :: giveDerivativeEta(FloatArray &n, const FloatArray &lc)
{
    double ksi, eta;
    ksi = lc.at(1); eta = lc.at(2);
    n.resize(8);

    n.at(1) =  0.25 * ( 1. + ksi ) * ( 2.0 * eta + ksi );
    n.at(2) =  0.25 * ( 1. - ksi ) * ( 2.0 * eta - ksi );
    n.at(3) = -0.25 * ( 1. - ksi ) * ( -2.0 * eta - ksi );
    n.at(4) = -0.25 * ( 1. + ksi ) * ( -2.0 * eta + ksi );
    n.at(5) =  0.5 * ( 1. - ksi * ksi );
    n.at(6) = -eta * ( 1. - ksi );
    n.at(7) = -0.5 * ( 1. - ksi * ksi );
    n.at(8) = -eta * ( 1. + ksi );

}

} // end namespace oofem
