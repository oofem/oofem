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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "feinterpol3d.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "gaussintegrationrule.h"

namespace oofem {
double FEInterpolation3d::giveVolume(const FEICellGeometry& cellgeo) const
{
    OOFEM_ERROR("Not implemented in subclass.");
}

IntArray FEInterpolation3d :: boundaryEdgeGiveNodes(int boundary, Element_Geometry_Type egt, bool includeHierarchical) const
{
    return this->computeLocalEdgeMapping(boundary);
}

void FEInterpolation3d :: boundaryEdgeEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const 
{
    this->edgeEvalN(answer, boundary, lcoords, cellgeo);
}

double FEInterpolation3d :: boundaryEdgeGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const 
{
    return this->edgeGiveTransformationJacobian(boundary, lcoords, cellgeo);
}

void FEInterpolation3d :: boundaryEdgeLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const 
{
    this->edgeLocal2global(answer, boundary, lcoords, cellgeo);
}

IntArray FEInterpolation3d :: boundaryGiveNodes(int boundary, Element_Geometry_Type egt) const
{
    return this->computeLocalSurfaceMapping(boundary);
}

void FEInterpolation3d :: boundaryEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    this->surfaceEvalN(answer, boundary, lcoords, cellgeo);
}

double FEInterpolation3d :: boundaryEvalNormal(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const 
{
    return this->surfaceEvalNormal(answer, boundary, lcoords, cellgeo);
}

double FEInterpolation3d :: boundaryGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    return this->surfaceGiveTransformationJacobian(boundary, lcoords, cellgeo);
}

void FEInterpolation3d :: boundaryLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const 
{
    return this->surfaceLocal2global(answer, boundary, lcoords, cellgeo);
}

IntArray FEInterpolation3d :: computeEdgeMapping(const IntArray &elemNodes, int iedge) const
{
    const auto &ln = this->computeLocalEdgeMapping(iedge);
    int size = ln.giveSize();
    IntArray edgeNodes(size);
    for ( int i = 1; i <= size; i++ ) {
        edgeNodes.at(i) = elemNodes.at( ln.at(i) );
    }
    return edgeNodes;
}

IntArray FEInterpolation3d :: computeSurfaceMapping(const IntArray &elemNodes, int isurf) const
{
    const auto &ln = this->computeLocalSurfaceMapping(isurf);
    int size = ln.giveSize();
    IntArray surfNodes(size);
    for ( int i = 1; i <= size; i++ ) {
        surfNodes.at(i) = elemNodes.at( ln.at(i) );
    }
    return surfNodes;
}

std::unique_ptr<IntegrationRule> FEInterpolation3d :: giveBoundaryEdgeIntegrationRule(int order, int boundary, Element_Geometry_Type egt) const
{
    auto iRule = std::make_unique<GaussIntegrationRule>(1, nullptr);
    int points = iRule->getRequiredNumberOfIntegrationPoints(_Line, order + this->order);
    iRule->SetUpPointsOnLine(points, _Unknown);
    return iRule;
}

void FEInterpolation3d :: edgeEvaldNdxi(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    OOFEM_ERROR("Not implemented");
}

void FEInterpolation3d :: surfaceEvaldNdx(FloatMatrix &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    OOFEM_ERROR("Not implemented");
}

double FEInterpolation3d :: edgeEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    OOFEM_ERROR("Not implemented");
}

IntArray FEInterpolation3d::boundarySurfaceGiveNodes(int boundary, Element_Geometry_Type egt, bool includeHierarchical) const
{
    return this->computeLocalSurfaceMapping(boundary);
}


double FEInterpolation3d::surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    auto [G1, G2] = this->surfaceEvalBaseVectorsAt(isurf, lcoords, cellgeo);
    auto normal = cross(G1, G2);
    double J = norm(normal);
    answer = normal;
    return J;
}

std::tuple<double, FloatArrayF<3>>
FEInterpolation3d :: surfaceEvalUnitNormal(int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    auto [G1, G2] = this->surfaceEvalBaseVectorsAt(isurf, lcoords, cellgeo);
    auto normal = cross(G1, G2);
    double J = norm(normal);
    return std::make_tuple(J, normal/J);
}


void FEInterpolation3d::surfaceEvaldNdxi(FloatMatrix & answer, const FloatArray & lcoords) const
{
    OOFEM_ERROR("Not implemented");
}

void FEInterpolation3d::surfaceEvald2Ndxi2(FloatMatrix & answer, const FloatArray & lcoords) const
{
    OOFEM_ERROR("Not implemented");
}
  


  std::tuple<FloatArrayF<3>, FloatArrayF<3>> FEInterpolation3d::surfaceEvalBaseVectorsAt(int isurf, const FloatArray & lcoords, const FEICellGeometry & cellgeo) const
{
    //Adapted from FEI3dQuadLin
    // Note: These are not normalized. Returns the two tangent vectors to the surface.
    FloatMatrix dNdxi;
    this->surfaceEvaldNdxi(dNdxi, lcoords);
    //Get nodes which correspond to the surface in question
    auto nodeIndices = this->computeLocalSurfaceMapping(isurf);
    FloatArrayF<3> G1, G2;
    for (int i = 0; i < nodeIndices.giveSize(); ++i) {
      G1 += dNdxi(i, 0) * FloatArrayF<3>(cellgeo.giveVertexCoordinates(nodeIndices(i)));
      G2 += dNdxi(i, 1) * FloatArrayF<3>(cellgeo.giveVertexCoordinates(nodeIndices(i)));
    }
    return std::make_tuple(G1, G2);
}

FloatMatrixF<3,3>
FEInterpolation3d :: surfaceGiveJacobianMatrixAt(int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
    // Jacobian matrix consists of the three curvilinear base vectors. The third is taken as the normal to the surface.
    // Note! The base vectors are not normalized except the third (normal)
    auto [G1, G2] = this->surfaceEvalBaseVectorsAt(isurf, lcoords, cellgeo);
    auto G3 = cross(G1, G2);
    FloatMatrixF<3,3> jacobianMatrix;
    jacobianMatrix.setColumn(G1,0);
    jacobianMatrix.setColumn(G2,1);
    jacobianMatrix.setColumn(G3,2);
    return jacobianMatrix;
}



double
FEInterpolation3d :: surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords,
                                                  const FEICellGeometry &cellgeo) const
{
    auto [J, N] =  surfaceEvalUnitNormal(isurf, lcoords, cellgeo);
    return J;
}
  



  
double FEInterpolation3d:: boundarySurfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const
{
  auto [J, N] = surfaceEvalUnitNormal(isurf, lcoords, cellgeo);
  answer = N;
  return  J;
}
  
  
} // end namespace oofem
