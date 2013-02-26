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

#include "fei3dwedgelin.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "node.h"
#include "mathfem.h"

namespace oofem {
void
FEI3dWedgeLin :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double x, y, z;
    answer.resize(6);

    x = lcoords.at(1);
    y = lcoords.at(2);
    z = lcoords.at(3);

    answer.at(1)  = 0.5 * ( 1. - x - y ) * ( 1. - z );
    answer.at(2)  = 0.5 * x * (1 - z);
    answer.at(3)  = 0.5 * y * (1 - z);
    answer.at(4)  = 0.5 * (1. - x - y ) * (1 + z);
    answer.at(5)  = 0.5 *  x * ( 1. + z );
    answer.at(6)  = 0.5 *  y * (1. + z );
}

void
FEI3dWedgeLin :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix(3, 3), inv(3, 3);
    FloatArray dx(6), dy(6), dz(6);
    double u, v, w;

    u = lcoords.at(1);
    v = lcoords.at(2);
    w = lcoords.at(3);

    this->giveJacobianMatrixAt(jacobianMatrix, lcoords, cellgeo);
    inv.beInverseOf(jacobianMatrix);

    this->giveDerivativeKsi(dx, v, w);
    this->giveDerivativeEta(dy, u, w);
    this->giveDerivativeDzeta(dz, u, v);

    answer.resize(6, 3);

    for ( int i = 1; i <= 6; i++ ) {
        answer.at(i, 1) = dx.at(i) * inv.at(1, 1) + dy.at(i) * inv.at(1, 2) + dz.at(i) * inv.at(1, 3);
        answer.at(i, 2) = dx.at(i) * inv.at(2, 1) + dy.at(i) * inv.at(2, 2) + dz.at(i) * inv.at(2, 3);
        answer.at(i, 3) = dx.at(i) * inv.at(3, 1) + dy.at(i) * inv.at(3, 2) + dz.at(i) * inv.at(3, 3);
    }
}


void
FEI3dWedgeLin :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double x, y, z;
    FloatArray n(6);

    x = lcoords.at(1);
    y = lcoords.at(2);
    z = lcoords.at(3);


    n.at(1)  = 0.5 * ( 1. - x - y ) * ( 1. - z );
    n.at(2)  = 0.5 * x * (1 - z);
    n.at(3)  = 0.5 * y * (1 - z);
    n.at(4)  = 0.5 * (1. - x - y ) * (1 + z);
    n.at(5)  = 0.5 *  x * ( 1. + z );
    n.at(6)  = 0.5 *  y * (1. + z );


    answer.resize(3);
    answer.zero();
    for ( int i = 1; i <= 6; i++ ) {
        answer.at(1) += n.at(i) * cellgeo.giveVertexCoordinates(i)->at(1);
        answer.at(2) += n.at(i) * cellgeo.giveVertexCoordinates(i)->at(2);
        answer.at(3) += n.at(i) * cellgeo.giveVertexCoordinates(i)->at(3);
    }
}


int
FEI3dWedgeLin :: global2local(FloatArray &answer, const FloatArray &coords, const FEICellGeometry &cellgeo)
{
   
     OOFEM_ERROR("FEI3dHexaQuad :: global2local not implemented");
     return 1;
}


double
FEI3dWedgeLin :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix(3, 3);

    this->giveJacobianMatrixAt(jacobianMatrix, lcoords, cellgeo);
    return jacobianMatrix.giveDeterminant()/2.;
}


void
FEI3dWedgeLin :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
// Returns the jacobian matrix  J (x,y,z)/(ksi,eta,dzeta)  of the receiver.
// Computes it if it does not exist yet.
{
    double x, y, z, u, v, w;
    FloatArray dx(6), dy(6), dz(6);

    jacobianMatrix.resize(3, 3);
    jacobianMatrix.zero();

    u = lcoords.at(1);
    v = lcoords.at(2);
    w = lcoords.at(3);

    this->giveDerivativeKsi(dx, v, w);
    this->giveDerivativeEta(dy, u, w);
    this->giveDerivativeDzeta(dz, u, v);

    for ( int i = 1; i <= 6; i++ ) {
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
FEI3dWedgeLin :: giveDerivativeKsi(FloatArray &dx, double v, double w)
{
    dx.at(1)  = -0.5 * ( 1. - w );
    dx.at(2)  =  0.5 * ( 1. - w );
    dx.at(3)  =  0.;
    dx.at(4)  = -0.5 * ( 1. + w );
    dx.at(5)  =  0.5 * ( 1. + w );
    dx.at(6)  =  0.;
}

void
FEI3dWedgeLin :: giveDerivativeEta(FloatArray &dy, double u, double w)
{
    dy.at(1)  = -0.5 * ( 1. - w );
    dy.at(2)  =  0.;
    dy.at(3)  =  0.5 * ( 1. - w );
    dy.at(4)  =  0.5 * ( 1. + w );
    dy.at(5)  =  0.;
    dy.at(6)  =  0.5 * ( 1. + w );
}

void
FEI3dWedgeLin :: giveDerivativeDzeta(FloatArray &dz, double u, double v)
{
    dz.at(1)  = -0.5 * ( 1. - u - v );
    dz.at(2)  = -0.5 * u;
    dz.at(3)  = -0.5 * v;
    dz.at(4)  =  0.5 * ( 1. - u - v );
    dz.at(5)  =  0.5 * u;
    dz.at(6)  =  0.5 * v;
}



void FEI3dWedgeLin :: edgeEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{ OOFEM_ERROR("FEI3dWedgeLin :: edgeEvalN not implemented"); }
void FEI3dWedgeLin :: edgeEvaldNdx(FloatMatrix &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{ OOFEM_ERROR("FEI3dWedgeLin :: edgeEvaldNdx not implemented"); }
void FEI3dWedgeLin :: edgeLocal2global(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{ OOFEM_ERROR("FEI3dWedgeLin :: edgeLocal2global not implemented"); }

void
FEI3dWedgeLin :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
{OOFEM_ERROR("FEI3dWedgeLin :: computeLocalEdgeMapping not implemented"); }

void
FEI3dWedgeLin :: surfaceEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{OOFEM_ERROR("FEI3dWedgeLin :: computeLocalEdgeMapping not implemented");}

void
FEI3dWedgeLin :: surfaceLocal2global(FloatArray &answer, int isurf,
                                     const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{OOFEM_ERROR("FEI3dWedgeLin :: computeLocalEdgeMapping not implemented");}

double
FEI3dWedgeLin :: surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords,
                                                   const FEICellGeometry &cellgeo)
{OOFEM_ERROR("FEI3dWedgeLin :: computeLocalEdgeMapping not implemented");return 0;}

void
FEI3dWedgeLin :: computeLocalSurfaceMapping(IntArray &nodes, int isurf)
{OOFEM_ERROR("FEI3dWedgeLin :: computeLocalEdgeMapping not implemented");}

void
FEI3dWedgeLin :: computeGlobalSurfaceMapping(IntArray &surfNodes, IntArray &elemNodes, int iSurf)
{OOFEM_ERROR("FEI3dWedgeLin :: computeLocalEdgeMapping not implemented");}

double FEI3dWedgeLin :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dWedgeLin :: edgeGiveTransformationJacobian not implemented");
    return 0.0;
}

} // end namespace oofem
