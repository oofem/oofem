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

#include "fei3dlinelin.h"
#include "mathfem.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "gaussintegrationrule.h"

namespace oofem {
double
FEI3dLineLin :: giveLength(const FEICellGeometry &cellgeo) const
{
    return distance(cellgeo.giveVertexCoordinates(2), cellgeo.giveVertexCoordinates(1));
}

void
FEI3dLineLin :: evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double ksi = lcoords.at(1);
    answer.resize(2);

    answer.at(1) = ( 1. - ksi ) * 0.5;
    answer.at(2) = ( 1. + ksi ) * 0.5;
}

double
FEI3dLineLin :: evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    ///@todo Not clear what this function should return. Just dNds would make sense if the caller defines a local coordinate system.
    FloatArray vec;
    vec.beDifferenceOf( cellgeo.giveVertexCoordinates(2), cellgeo.giveVertexCoordinates(1) );

    double detJ = vec.computeSquaredNorm() * 0.5;
    double l2_inv = 0.5 / detJ;
    answer.resize(2, 3);

    answer.at(1, 1) = -vec.at(1) * l2_inv;
    answer.at(2, 1) =  vec.at(1) * l2_inv;
    answer.at(1, 2) = -vec.at(2) * l2_inv;
    answer.at(2, 2) =  vec.at(2) * l2_inv;
    answer.at(1, 3) = -vec.at(3) * l2_inv;
    answer.at(2, 3) =  vec.at(3) * l2_inv;

    return detJ;
}

void
FEI3dLineLin :: evald2Ndx2(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    answer.resize(2, 3); ///@todo Check this part
    answer.zero();
}

void
FEI3dLineLin :: local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double ksi = lcoords.at(1);

    answer.beScaled( ( 1. - ksi ) * 0.5, cellgeo.giveVertexCoordinates(1) );
    answer.add( ( 1. + ksi ) * 0.5, cellgeo.giveVertexCoordinates(2) );
}


int
FEI3dLineLin :: global2local(FloatArray &answer, const FloatArray &coords, const FEICellGeometry &cellgeo) const
{
    FloatArray vec, x;
    vec.beDifferenceOf( cellgeo.giveVertexCoordinates(2), cellgeo.giveVertexCoordinates(1) );
    x.beDifferenceOf( coords, cellgeo.giveVertexCoordinates(1) );
    double l2 = vec.computeSquaredNorm();
    double xvec = x.dotProduct(vec);

    answer = FloatArray{2.0 * xvec / l2 - 1.0};
    answer.at(1) = clamp(answer.at(1), -1.0, 1.0);
    return false; // No point to check if point is "inside".
}


double
FEI3dLineLin :: giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    return 0.5 * this->giveLength(cellgeo);
}


void
FEI3dLineLin :: edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    this->evalN(answer, lcoords, cellgeo);
}

void
FEI3dLineLin :: edgeEvaldNdx(FloatMatrix &answer, int iedge,
                             const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    double l_inv = 1.0 / this->giveLength(cellgeo);

    answer.resize(2, 1);
    answer.at(1, 1) = -l_inv;
    answer.at(2, 1) =   l_inv;
}

void
FEI3dLineLin :: edgeLocal2global(FloatArray &answer, int iedge,
                                 const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    return this->local2global(answer, lcoords, cellgeo);
}


double
FEI3dLineLin :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    return this->giveTransformationJacobian(lcoords, cellgeo);
}


IntArray
FEI3dLineLin :: computeLocalEdgeMapping(int iedge) const
{
    if ( iedge != 1 ) {
        OOFEM_ERROR("wrong edge number (%d)", iedge);
    }
    return {1, 2};
}

void
FEI3dLineLin :: surfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    OOFEM_ERROR("no surfaces available");
}

double
FEI3dLineLin :: surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    OOFEM_ERROR("no surfaces available");
    //return 0.0;
}

void
FEI3dLineLin :: surfaceLocal2global(FloatArray &answer, int iedge,
                                    const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    OOFEM_ERROR("no surfaces available");
}

double
FEI3dLineLin :: surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords,
                                                  const FEICellGeometry &cellgeo) const
{
    OOFEM_ERROR("no surfaces available");
    //return 0.0;
}

IntArray
FEI3dLineLin :: computeLocalSurfaceMapping(int isurf) const
{
    OOFEM_ERROR("no surfaces available");
    //return {};
}


void
FEI3dLineLin :: giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
// Returns the jacobian matrix  J (x,y,z)/(ksi,eta,dzeta)  of the receiver.
{
    ///@todo Not sure about this matrix
    jacobianMatrix.resize(1, 1);
    jacobianMatrix.at(1, 1) = 1.0;
}

std::unique_ptr<IntegrationRule>
FEI3dLineLin :: giveIntegrationRule(int order, Element_Geometry_Type egt) const
{
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Line, order + 0);
    iRule->SetUpPointsOnLine(points, _Unknown);
    return std::move(iRule);
}

std::unique_ptr<IntegrationRule>
FEI3dLineLin :: giveBoundaryIntegrationRule(int order, int boundary, Element_Geometry_Type egt) const
{
    ///@todo Not sure about this.
    OOFEM_ERROR("Not supported");
    //return nullptr;
}
} // end namespace oofem
