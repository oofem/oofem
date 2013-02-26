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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#include "fei3dwedgequad.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "node.h"
#include "mathfem.h"

namespace oofem {

void
FEI3dWedgeQuad :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double x, y, z;
    answer.resize(15);

    x = lcoords.at(1);
    y = lcoords.at(2);
    z = lcoords.at(3);

    answer.at(1) = 0.5 * ((1. - x - y ) * (2. * (1. - x - y) - 1.) * ( 1. - z ) - (1. - x - y) * (1. - z * z));
    answer.at(2) = 0.5 * (x * (2. * x - 1.) * (1. - z ) - x * (1. - z * z));
    answer.at(3) = 0.5 * (y * (2. * y - 1.) * (1. - z ) - y * (1. - z * z));
    answer.at(4) = 0.5 * ((1. - x - y ) * (2. * (1. - x - y) - 1.) * (1. + z) - (1. - x - y) * (1. - z * z));
    answer.at(5) = 0.5 *  (x * (2. * x - 1.) * (1. + z) - x * (1. - z * z));
    answer.at(6) = 0.5 *  (y * (2. * y - 1.) * (1. + z) - y * (1. - z * z));
    answer.at(7) = 2. * (1. - x - y) * x * (1. - z);
    answer.at(8) = 2. * x * y * (1. - z);
    answer.at(9) = 2. * y * (1. -x - y) * (1. - z);
    answer.at(10) = 2.*  x * (1. - x - y) * (1. + z);
    answer.at(11) = 2. * x * y * (1. + z);
    answer.at(12) = 2. * y * (1. - x - y) * (1. + z);
    answer.at(13) = (1. -x - y) * (1. - z * z);
    answer.at(14) = x * (1. - z * z);
    answer.at(15) = y * (1. - z * z);
    return;
}

void
FEI3dWedgeQuad :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    int i;
    FloatMatrix jacobianMatrix(3, 3), inv(3, 3);
    FloatArray dx(15), dy(15), dz(15);
    double u, v, w;

    u = lcoords.at(1);
    v = lcoords.at(2);
    w = lcoords.at(3);

    this->giveJacobianMatrixAt(jacobianMatrix, lcoords, cellgeo);
    inv.beInverseOf(jacobianMatrix);

    this->giveDerivativeKsi(dx, u, v, w);
    this->giveDerivativeEta(dy, u, v, w);
    this->giveDerivativeDzeta(dz, u, v, w);

    answer.resize(15, 3);

    for ( i = 1; i <= 15; i++ ) {
        answer.at(i, 1) = dx.at(i) * inv.at(1, 1) + dy.at(i) * inv.at(1, 2) + dz.at(i) * inv.at(1, 3);
        answer.at(i, 2) = dx.at(i) * inv.at(2, 1) + dy.at(i) * inv.at(2, 2) + dz.at(i) * inv.at(2, 3);
        answer.at(i, 3) = dx.at(i) * inv.at(3, 1) + dy.at(i) * inv.at(3, 2) + dz.at(i) * inv.at(3, 3);
    }
}


void
FEI3dWedgeQuad :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray n(15);
    this -> evalN(n, lcoords, cellgeo);



    answer.resize(3);
    answer.zero();
    for ( int i = 1; i <= 15; i++ ) {
        answer.at(1) += n.at(i) * cellgeo.giveVertexCoordinates(i)->at(1);
        answer.at(2) += n.at(i) * cellgeo.giveVertexCoordinates(i)->at(2);
        answer.at(3) += n.at(i) * cellgeo.giveVertexCoordinates(i)->at(3);
    }
}


int
FEI3dWedgeQuad :: global2local(FloatArray &answer, const FloatArray &coords, const FEICellGeometry &cellgeo)
{
   
     OOFEM_ERROR("FEI3dHexaQuad :: global2local not implemented");
     return 1;
}


double
FEI3dWedgeQuad :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix(3, 3);

    this->giveJacobianMatrixAt(jacobianMatrix, lcoords, cellgeo);
    return jacobianMatrix.giveDeterminant()/2.;
}


void
FEI3dWedgeQuad :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
// Returns the jacobian matrix  J (x,y,z)/(ksi,eta,dzeta)  of the receiver.
// Computes it if it does not exist yet.
{
    int i;
    double x, y, z, u, v, w;
    FloatArray dx(15), dy(15), dz(15);

    jacobianMatrix.resize(3, 3);
    jacobianMatrix.zero();

    u = lcoords.at(1);
    v = lcoords.at(2);
    w = lcoords.at(3);

    this->giveDerivativeKsi(dx, u, v, w);
    this->giveDerivativeEta(dy, u, v,  w);
    this->giveDerivativeDzeta(dz, u, v, w);

    for ( i = 1; i <= 15; i++ ) {
        x = cellgeo.giveVertexCoordinates(i)->at(1);
        y = cellgeo.giveVertexCoordinates(i)->at(2);
        z = cellgeo.giveVertexCoordinates(i)->at(3);

        jacobianMatrix.at(1, 1) += dx.at(i) * x;
        jacobianMatrix.at(1, 2) += dx.at(i) * y;
        jacobianMatrix.at(1, 3) += dx.at(i) * z;
        jacobianMatrix.at(2, 1) += dy.at(i) * x;
        jacobianMatrix.at(2, 2) += dy.at(i) * y;
        jacobianMatrix.at(2, 3) += dy.at(i) * z;
        jacobianMatrix.at(3, 1) += dz.at(i) * x;
        jacobianMatrix.at(3, 2) += dz.at(i) * y;
        jacobianMatrix.at(3, 3) += dz.at(i) * z;
    }
}


void
FEI3dWedgeQuad :: giveDerivativeKsi(FloatArray &dx, double x, double y, double z)
{

  dx.at(1) = 1/2 - ((z - 1)*(x + y - 1))/2 - z*z/2 - ((x + y)*(z - 1))/2; 
  dx.at(2) = z*z/2 - (x*(z - 1))/2 - ((x - 1)*(z - 1))/2 - 1/2;
  dx.at(3) = 0;
  dx.at(4) = ((x + y)*(z + 1))/2 + ((z + 1)*(x + y - 1))/2 - z*z/2 + 1/2;
  dx.at(5) = ((x - 1)*(z + 1))/2 + (x*(z + 1))/2 + z*z/2 - 1/2;
  dx.at(6) = 0;
  dx.at(7) = (z - 1)*(x + y - 1) + x*(z - 1);
  dx.at(8) = -y*(z - 1);
  dx.at(9) = y*(z - 1);
  dx.at(10) = - (z + 1)*(x + y - 1) - x*(z + 1);
  dx.at(11) = y*(z + 1);
  dx.at(12) = -y*(z + 1); 
  dx.at(13) = z*z - 1;
  dx.at(14) = 1 - z*z;
  dx.at(15) = 0;
}

void
FEI3dWedgeQuad :: giveDerivativeEta(FloatArray &dy, double x, double y, double z)
{


  
  dy.at(1) = 1/2 - ((z - 1)*(x + y - 1))/2 - z*z/2 - ((x + y)*(z - 1))/2;
  dy.at(2) = 0;
  dy.at(3) = z*z/2 - (y*(z - 1))/2 - ((y - 1)*(z - 1))/2 - 1/2;
  dy.at(4) = ((x + y)*(z + 1))/2 + ((z + 1)*(x + y - 1))/2 - z*z/2 + 1/2;
  dy.at(5) = 0;
  dy.at(6) = ((y - 1)*(z + 1))/2 + (y*(z + 1))/2 + z*z/2 - 1/2;
  dy.at(7) = x*(z - 1);
  dy.at(8) = -x*(z - 1);
  dy.at(9) = (z - 1)*(x + y - 1) + y*(z - 1);
  dy.at(10) = -x*(z + 1);
  dy.at(11) = x*(z + 1);
  dy.at(12) = - (z + 1)*(x + y - 1) - y*(z + 1);
  dy.at(13) = z*z - 1;
  dy.at(14) = 0;
  dy.at(15) = 1 - z*z;

 
}

void
FEI3dWedgeQuad :: giveDerivativeDzeta(FloatArray &dz,  double x, double y, double z)
{

  dz.at(1) = - ((x + y)*(x + y - 1))/2 - z*(x + y - 1);
  dz.at(2) = x*z - (x*(x - 1))/2;
  dz.at(3) = y*z - (y*(y - 1))/2;
  dz.at(4) = ((x + y)*(x + y - 1))/2 - z*(x + y - 1);
  dz.at(5) = x*z + (x*(x - 1))/2;
  dz.at(6) = y*z + (y*(y - 1))/2;
  dz.at(7) = x*(x + y - 1);
  dz.at(8) = -x*y;
  dz.at(9) = y*(x + y - 1);
  dz.at(10) = -x*(x + y - 1);
  dz.at(11) = x*y;
  dz.at(12) = -y*(x + y - 1);
  dz.at(13) = 2*z*(x + y - 1);
  dz.at(14) = (-2)*x*z; 
  dz.at(15) = (-2)*y*z;

}



void FEI3dWedgeQuad :: edgeEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{ OOFEM_ERROR("FEI3dWedgeQuad :: edgeEvalN not implemented"); }
void FEI3dWedgeQuad :: edgeEvaldNdx(FloatMatrix &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{ OOFEM_ERROR("FEI3dWedgeQuad :: edgeEvaldNdx not implemented"); }
void FEI3dWedgeQuad :: edgeLocal2global(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{ OOFEM_ERROR("FEI3dWedgeQuad :: edgeLocal2global not implemented"); }

void
FEI3dWedgeQuad :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
{OOFEM_ERROR("FEI3dWedgeQuad :: computeLocalEdgeMapping not implemented"); }

void
FEI3dWedgeQuad :: surfaceEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{OOFEM_ERROR("FEI3dWedgeQuad :: computeLocalEdgeMapping not implemented");}

void
FEI3dWedgeQuad :: surfaceLocal2global(FloatArray &answer, int isurf,
                                     const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{OOFEM_ERROR("FEI3dWedgeQuad :: computeLocalEdgeMapping not implemented");}

double
FEI3dWedgeQuad :: surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords,
                                                   const FEICellGeometry &cellgeo)
{OOFEM_ERROR("FEI3dWedgeQuad :: computeLocalEdgeMapping not implemented");return 0;}

void
FEI3dWedgeQuad :: computeLocalSurfaceMapping(IntArray &nodes, int isurf)
{OOFEM_ERROR("FEI3dWedgeQuad :: computeLocalEdgeMapping not implemented");}

void
FEI3dWedgeQuad :: computeGlobalSurfaceMapping(IntArray &surfNodes, IntArray &elemNodes, int iSurf)
{OOFEM_ERROR("FEI3dWedgeQuad :: computeLocalEdgeMapping not implemented");}

double FEI3dWedgeQuad :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
  OOFEM_ERROR("FEI3dWedgeQuad :: edgeGiveTransformationJacobian not implemented");
  return 0.0;
}

} // end namespace oofem
