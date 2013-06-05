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

#include "fei3dwedgequad.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "intarray.h"

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
}

double
FEI3dWedgeQuad :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatMatrix jacobianMatrix, inv, dNduvw, coords;
    this->giveLocalDerivative(dNduvw, lcoords);
    coords.resize(3, 15);
    for ( int i = 1; i <= 15; i++ ) {
        coords.setColumn(*cellgeo.giveVertexCoordinates(i), i);
    }
    jacobianMatrix.beProductOf(coords, dNduvw);
    inv.beInverseOf(jacobianMatrix);

    answer.beProductOf(dNduvw, inv);
    return jacobianMatrix.giveDeterminant();
}


void
FEI3dWedgeQuad :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray n(15);
    this->evalN(n, lcoords, cellgeo);

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
    return jacobianMatrix.giveDeterminant()/2.; ///@todo Should this really be a factor 1/2 here?
}


void
FEI3dWedgeQuad :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
// Returns the jacobian matrix  J (x,y,z)/(ksi,eta,dzeta)  of the receiver.
{
    FloatMatrix dNduvw, coords;
    this->giveLocalDerivative(dNduvw, lcoords);
    coords.resize(3, 15);
    for ( int i = 1; i <= 15; i++ ) {
        coords.setColumn(*cellgeo.giveVertexCoordinates(i), i);
    }
    jacobianMatrix.beProductOf(coords, dNduvw);
}


void
FEI3dWedgeQuad :: giveLocalDerivative(FloatMatrix &dN, const FloatArray &lcoords)
{
    double x, y, z;
    x = lcoords.at(1);
    y = lcoords.at(2);
    z = lcoords.at(3);

    dN.resize(15, 3);

    dN.at(1, 1)  =  1./2. - (z - 1.)*(x + y - 1.) - z*z/2. - ((z - 1.)*(2.*x + 2.*y - 1.))/2.;
    dN.at(2, 1)  =  z*z/2. - x*(z - 1.) - ((2.*x - 1.)*(z - 1.))/2. - 1/2.;
    dN.at(3, 1)  =  0.;
    dN.at(4, 1)  = ((z + 1.)*(2.*x + 2.*y - 1.))/2. + (z + 1.)*(x + y - 1.) - z*z/2. + 1./2.;
    dN.at(5, 1)  = ((2.*x - 1.)*(z + 1.))/2. + x*(z + 1.) + z*z/2. - 1./2.;
    dN.at(6, 1)  =  0.;
    dN.at(7, 1)  = (z - 1.)*(2.*x + 2.*y - 2.) + 2.*x*(z - 1.);
    dN.at(8, 1)  = -2.*y*(z - 1.);
    dN.at(9, 1)  =  2.*y*(z - 1.);
    dN.at(10, 1) = -2.*(z + 1.)*(x + y - 1.) - 2.*x*(z + 1.);
    dN.at(11, 1) =  2.*y*(z + 1.);
    dN.at(12, 1) = -2.*y*(z + 1.);
    dN.at(13, 1) =  z*z - 1.;
    dN.at(14, 1) =  1. - z*z;
    dN.at(15, 1) =  0.;


    dN.at(1, 2) =   1./2. - (z - 1.)*(x + y - 1.) - z*z/2. - ((z - 1.)*(2.*x + 2.*y - 1.))/2.;
    dN.at(2, 2) =   0.;
    dN.at(3, 2) =   z*z/2. - y*(z - 1.) - ((2.*y - 1.)*(z - 1.))/2. - 1./2.;
    dN.at(4, 2) = ((z + 1.)*(2.*x + 2.*y - 1.))/2. + (z + 1.)*(x + y - 1.) - z*z/2. + 1./2.;
    dN.at(5, 2) =   0.;
    dN.at(6, 2) = ((2.*y - 1.)*(z + 1.))/2. + y*(z + 1.) + z*z/2. - 1./2.;
    dN.at(7, 2) =   2.*x*(z - 1.);
    dN.at(8, 2) =  -2.*x*(z - 1.);
    dN.at(9, 2) =   2.*(z - 1.)*(x + y - 1.) + 2.*y*(z - 1.);
    dN.at(10, 2) = -2.*x*(z + 1.);
    dN.at(11, 2) =  2.*x*(z + 1.);
    dN.at(12, 2) = -2.*(z + 1.)*(x + y - 1.) - 2.*y*(z + 1.);
    dN.at(13, 2) =  z*z - 1.;
    dN.at(14, 2) =  0.;
    dN.at(15, 2) =  1. - z*z;


    dN.at(1, 3) =  -z*(x + y - 1.) - ((2.*x + 2.*y - 1.)*(x + y - 1.))/2.;
    dN.at(2, 3) =   x*z - (x*(2.*x - 1.))/2.;
    dN.at(3, 3) =   y*z - (y*(2*y - 1.))/2.;
    dN.at(4, 3) =  ((2.*x + 2.*y - 1.)*(x + y - 1.))/2. - z*(x + y - 1.);
    dN.at(5, 3) =  (x*(2.*x - 1.))/2. + x*z;
    dN.at(6, 3) =  (y*(2.*y - 1.))/2. + y*z;
    dN.at(7, 3) =   x*(2.*x + 2.*y - 2.);
    dN.at(8, 3) =  -2.*x*y;
    dN.at(9, 3) =   2.*y*(x + y - 1.);
    dN.at(10, 3) = -2.*x*(x + y - 1.);
    dN.at(11, 3) =  2.*x*y;
    dN.at(12, 3) = -2.*y*(x + y - 1.);
    dN.at(13, 3) =  2.*z*(x + y - 1.);
    dN.at(14, 3) = -2.*x*z;
    dN.at(15, 3) = -2.*y*z;

}


void FEI3dWedgeQuad :: edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double ksi = lcoords.at(1);
    answer.resize(3);
    answer.at(1) = ksi* ( ksi - 1. ) * 0.5;
    answer.at(2) = ksi* ( 1. + ksi ) * 0.5;
    answer.at(3) = (1. - ksi * ksi);
}


void FEI3dWedgeQuad :: edgeEvaldNdx(FloatMatrix &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dWedgeQuad :: edgeEvaldNdx not implemented");
}


void FEI3dWedgeQuad :: edgeLocal2global(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray nodes;
    FloatArray n;

    this->computeLocalEdgeMapping(nodes, iedge);
    this->edgeEvalN(n, iedge, lcoords, cellgeo);

    answer.resize(0);
    for ( int i = 1; i <= n.giveSize(); ++i ) {
        answer.add( n.at(i), * cellgeo.giveVertexCoordinates(nodes.at(i)));
    }
}


void
FEI3dWedgeQuad :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
{
    if ( iedge == 1 ) {
        edgeNodes.setValues(3, 1, 2, 7);
    } else if ( iedge == 2 ) {
        edgeNodes.setValues(3, 2, 3, 8);
    } else if ( iedge == 3 ) {
        edgeNodes.setValues(3, 3, 1, 9);
    } else if ( iedge == 4 ) {
        edgeNodes.setValues(3, 4, 5, 10);
    } else if ( iedge == 5 ) {
        edgeNodes.setValues(3, 5, 6, 11);
    } else if ( iedge == 6 ) {
        edgeNodes.setValues(3, 6, 4, 12);
    } else if ( iedge == 7 ) {
        edgeNodes.setValues(3, 1, 4, 13);
    } else if ( iedge == 8 ) {
        edgeNodes.setValues(3, 2, 5, 14);
    } else if ( iedge == 9 ) {
        edgeNodes.setValues(3, 3, 6, 15);
    } else {
        OOFEM_ERROR2("FEI3dWedgeQuad :: computeLocalEdgeMapping - Edge %d doesn't exist.\n", iedge);
    }
}


double FEI3dWedgeQuad :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dWedgeQuad :: edgeGiveTransformationJacobian not implemented");
    return 0.0;
}


void
FEI3dWedgeQuad :: surfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double ksi = lcoords.at(1);
    double eta = lcoords.at(2);

    if ( isurf <= 2 ) {
        answer.resize(3);
        answer.at(1) = ksi;
        answer.at(2) = eta;
        answer.at(3) = 1.0 - ksi - eta;
    } else {
        answer.resize(8);
        answer.at(1) = ( 1. + ksi ) * ( 1. + eta ) * 0.25 * ( ksi + eta - 1. );
        answer.at(2) = ( 1. - ksi ) * ( 1. + eta ) * 0.25 * ( -ksi + eta - 1. );
        answer.at(3) = ( 1. - ksi ) * ( 1. - eta ) * 0.25 * ( -ksi - eta - 1. );
        answer.at(4) = ( 1. + ksi ) * ( 1. - eta ) * 0.25 * ( ksi - eta - 1. );
        answer.at(5) = 0.5 * ( 1. - ksi * ksi ) * ( 1. + eta );
        answer.at(6) = 0.5 * ( 1. - ksi ) * ( 1. - eta * eta );
        answer.at(7) = 0.5 * ( 1. - ksi * ksi ) * ( 1. - eta );
        answer.at(8) = 0.5 * ( 1. + ksi ) * ( 1. - eta * eta );
    }
}


void
FEI3dWedgeQuad :: surfaceLocal2global(FloatArray &answer, int isurf,
                                     const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    IntArray nodes;
    FloatArray n;

    this->computeLocalSurfaceMapping(nodes, isurf);
    this->surfaceEvalN(n, isurf, lcoords, cellgeo);

    answer.resize(0);
    for ( int i = 1; i <= n.giveSize(); ++i ) {
        answer.add( n.at(i), * cellgeo.giveVertexCoordinates(nodes.at(i)));
    }
}


void
FEI3dWedgeQuad :: computeLocalSurfaceMapping(IntArray &nodes, int isurf)
{
    if ( isurf == 1 ) {
        nodes.setValues(6, 1, 2, 3, 7, 8, 9);
    } else if ( isurf == 2 ) {
        nodes.setValues(6, 4, 5, 6, 10, 11, 12);
    } else if ( isurf == 3 ) {
        nodes.setValues(8, 1, 2, 5, 4, 7, 14, 10, 13);
    } else if ( isurf == 4 ) {
        nodes.setValues(8, 2, 3, 6, 5, 8, 15, 11, 14);
    } else if ( isurf == 5 ) {
        nodes.setValues(8, 3, 1, 4, 6, 9, 13, 12, 15);
    } else {
        OOFEM_ERROR2("FEI3dWedgeQuad :: computeLocalSurfaceMapping - Surface %d doesn't exist.\n", isurf);
    }
}


double
FEI3dWedgeQuad :: surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dWedgeQuad :: surfaceGiveTransformationJacobian not implemented");
    return 0;
}

} // end namespace oofem
