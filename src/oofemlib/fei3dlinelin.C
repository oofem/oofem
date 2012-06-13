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

#include "fei3dlinelin.h"
#include "mathfem.h"
#include "flotmtrx.h"
#include "flotarry.h"

namespace oofem {

double
FEI3dLineLin :: giveLength(const FEICellGeometry &cellgeo) const
{
    return cellgeo.giveVertexCoordinates(2)->distance(*cellgeo.giveVertexCoordinates(1));
}

void
FEI3dLineLin :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double ksi = lcoords.at(1);
    answer.resize(2);

    answer.at(1) = ( 1. - ksi ) * 0.5;
    answer.at(2) = ( 1. + ksi ) * 0.5;
}

void
FEI3dLineLin :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray vec;
    vec.beDifferenceOf(*cellgeo.giveVertexCoordinates(2), *cellgeo.giveVertexCoordinates(1));

    double l2_inv = 1.0 / vec.computeSquaredNorm();
    answer.resize(2, 3);

    answer.at(1, 1) = -vec.at(1)*l2_inv;
    answer.at(2, 1) =  vec.at(1)*l2_inv;
    answer.at(1, 2) = -vec.at(2)*l2_inv;
    answer.at(2, 2) =  vec.at(2)*l2_inv;
    answer.at(1, 3) = -vec.at(3)*l2_inv;
    answer.at(2, 3) =  vec.at(3)*l2_inv;
}

void
FEI3dLineLin :: evald2Ndx2(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    answer.resize(2,3); ///@todo Check this part
    answer.zero();
}

void
FEI3dLineLin :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double ksi;
    FloatArray n(2);
    ksi = lcoords.at(1);

    answer.resize(0);
    answer.add( ( 1. - ksi ) * 0.5, *cellgeo.giveVertexCoordinates(1) );
    answer.add( ( 1. + ksi ) * 0.5, *cellgeo.giveVertexCoordinates(2) );
}


int
FEI3dLineLin :: global2local(FloatArray &answer, const FloatArray &coords, const FEICellGeometry &cellgeo)
{
    FloatArray vec, x;
    vec.beDifferenceOf(*cellgeo.giveVertexCoordinates(2), *cellgeo.giveVertexCoordinates(1));
    x.beDifferenceOf(coords, *cellgeo.giveVertexCoordinates(1));
    double l2 = vec.computeSquaredNorm();
    double xvec = x.dotProduct(vec);

    answer.setValues(1, 2.0 * xvec / l2 - 1.0);
    return 1;
}


double
FEI3dLineLin :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    return 0.5 * this->giveLength(cellgeo);
}


void
FEI3dLineLin :: edgeEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    this->evalN(answer, lcoords, cellgeo);
}

void
FEI3dLineLin :: edgeEvaldNdx(FloatMatrix &answer, int iedge,
                             const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    double l = this->giveLength(cellgeo);

    answer.resize(2, 1);
    answer.at(1, 1) = -1.0 / l;
    answer.at(2, 1) =  1.0 / l;
}

void
FEI3dLineLin :: edgeLocal2global(FloatArray &answer, int iedge,
                                 const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    return this->local2global(answer, lcoords, cellgeo);
}


double
FEI3dLineLin :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    return this->giveTransformationJacobian(lcoords, cellgeo);
}


void
FEI3dLineLin :: computeLocalEdgeMapping(IntArray &edgeNodes, int iedge)
{
    if ( iedge != 1 ) {
        OOFEM_ERROR2("FEI3dLineLin :: computeEdgeMapping: wrong edge number (%d)", iedge);
    }
    edgeNodes.setValues(2, 1, 2);
}

void
FEI3dLineLin :: surfaceEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dLineLin :: computeEdgeMapping: no surfaces available");
}

double
FEI3dLineLin :: surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dLineLin :: surfaceEvalNormal: no surfaces available");
}

void
FEI3dLineLin :: surfaceLocal2global(FloatArray &answer, int iedge,
                                    const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dLineLin :: computeEdgeMapping: no surfaces available");
}

double
FEI3dLineLin :: surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords,
                                                  const FEICellGeometry &cellgeo)
{
    OOFEM_ERROR("FEI3dLineLin :: computeEdgeMapping: no surfaces available");
    return 0.0;
}

void
FEI3dLineLin :: computeLocalSurfaceMapping(IntArray &surfNodes, int isurf)
{
    OOFEM_ERROR("FEI3dLineLin :: computeEdgeMapping: no surfaces available");
}


void
FEI3dLineLin :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
// Returns the jacobian matrix  J (x,y,z)/(ksi,eta,dzeta)  of the receiver.
// Computes it if it does not exist yet.
{
    ///@todo Not sure about this matrix
    jacobianMatrix.resize(1,1);
    jacobianMatrix.at(1,1) = 1.0;
}

} // end namespace oofem
