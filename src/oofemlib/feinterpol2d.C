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

#include "feinterpol2d.h"
#include "floatarray.h"
#include "gaussintegrationrule.h"

namespace oofem {
void FEInterpolation2d :: boundaryEdgeGiveNodes(IntArray &answer, int boundary)
{
  this->computeLocalEdgeMapping(answer, boundary);
}

void FEInterpolation2d :: boundaryEdgeEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
  this->edgeEvalN(answer, boundary, lcoords, cellgeo);
}

double FEInterpolation2d :: boundaryEdgeGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    return this->edgeGiveTransformationJacobian(boundary, lcoords, cellgeo);
}

void FEInterpolation2d :: boundaryEdgeLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    this->edgeLocal2global(answer, boundary, lcoords, cellgeo);
}

double FEInterpolation2d :: giveArea(const FEICellGeometry &cellgeo) const
{
    OOFEM_ERROR("Not implemented in subclass.");
    return 0;
}

void FEInterpolation2d :: boundaryGiveNodes(IntArray &answer, int boundary)
{
    this->computeLocalEdgeMapping(answer, boundary);
}

void FEInterpolation2d :: boundaryEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    this->edgeEvalN(answer, boundary, lcoords, cellgeo);
}

double FEInterpolation2d :: boundaryEvalNormal(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    return this->edgeEvalNormal(answer, boundary, lcoords, cellgeo);
}

double FEInterpolation2d :: boundaryGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    return this->edgeGiveTransformationJacobian(boundary, lcoords, cellgeo);
}

void FEInterpolation2d :: boundaryLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    return this->edgeLocal2global(answer, boundary, lcoords, cellgeo);
}

void FEInterpolation2d :: computeEdgeMapping(IntArray &edgeNodes, IntArray &elemNodes, int iedge)
{
    IntArray ln;
    this->computeLocalEdgeMapping(ln, iedge);
    int size = ln.giveSize();
    edgeNodes.resize(size);
    for ( int i = 1; i <= size; i++ ) {
        edgeNodes.at(i) = elemNodes.at( ln.at(i) );
    }
}

double FEInterpolation2d :: edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
    FloatArray normal;
    return this->edgeEvalNormal(normal, iedge, lcoords, cellgeo);
}

void FEInterpolation2d::boundarySurfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
  this->evalN(answer, lcoords, cellgeo);
}

void FEInterpolation2d::boundarySurfaceEvaldNdx(FloatMatrix &answer, int isurf,
					const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
  this->evaldNdx(answer, lcoords, cellgeo);
}

double FEInterpolation2d::boundarySurfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords,
					    const FEICellGeometry &cellgeo)
{
  answer = {0,0,1};
  return this->giveTransformationJacobian(lcoords, cellgeo);
}

void FEInterpolation2d::boundarySurfaceLocal2global(FloatArray &answer, int isurf,
					    const FloatArray &lcoords, const FEICellGeometry &cellgeo)
{
  this->local2global(answer, lcoords, cellgeo);
}

double FEInterpolation2d::boundarySurfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords,
							    const FEICellGeometry &cellgeo)
{
  return this->giveTransformationJacobian(lcoords, cellgeo);
}

void FEInterpolation2d::boundarySurfaceGiveNodes(IntArray &answer, int boundary)
{
  int nnode = this->giveNumberOfNodes();
  answer.resize(nnode);
  for (int i =1; i<=nnode; i++) {
    answer.at(i)=i;
  }
}

  
} // end namespace oofem
