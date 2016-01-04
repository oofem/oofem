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

#ifndef fei3dhexaquad_h
#define fei3dhexaquad_h

#include "feinterpol3d.h"

namespace oofem {
/**
 * Class representing implementation of quadratic hexahedra interpolation class.
 *  ***  numbering like T3d  ***            ***  numbering like T3d  ***
 *  ***  numbering of nodes  ***            ***  numbering of surfs  ***
 *
 *              zeta
 *        1      ^  9         2
 *         +-----|--+--------+                     +-----------------+
 *        /|     |          /|                    /|                /|
 *       / |     |         / |                   / |               / |
 *    12+  |     o      10+  |                  /  |     1        /  |
 *     /   |     |       /   |                 /   |             /   |
 *  4 /  17+  11 |    3 /    +18              /    |       (3)  /    |
 *   +--------+--------+     |               +-----------------+     |
 *   |     |     |     |     |               |     |           |     |
 *   |     |     +-----|--o------> eta       | (6) |           |  4  |
 *   |     |    /   13 |     |               |     |           |     |
 *   |   5 +---/----+--|-----+ 6             |     +-----------|-----+
 * 20+    /   o        +19  /                |    /   5        |    /
 *   |   /   /         |   /                 |   /             |   /
 *   |16+   /          |  +14                |  /       (2)    |  /
 *   | /   /           | /                   | /               | /
 *   |/   L ksi        |/                    |/                |/
 *   +--------+--------+                     +-----------------+
 *  8         15        7
 *
 * @author Ladislav Svoboda
 * @author Mikael Öhman
 */
class OOFEM_EXPORT FEI3dHexaQuad : public FEInterpolation3d
{
public:
    FEI3dHexaQuad() : FEInterpolation3d(2) { }

    virtual integrationDomain giveIntegrationDomain() const { return _Cube; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_hexa_2; }

    virtual double giveCharacteristicLength(const FEICellGeometry &cellgeo) const;

    // Bulk
    virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual int  global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual int giveNumberOfNodes() const { return 20; }
    
    // Edge
    virtual void edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void edgeEvaldNdx(FloatMatrix &answer, int iedge,
                              const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void edgeLocal2global(FloatArray &answer, int iedge,
                                  const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords,
                                                  const FEICellGeometry &cellgeo);
    virtual void computeLocalEdgeMapping(IntArray &edgeNodes, int iedge);

    // Surface
    virtual void surfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    //virtual void surfaceEvaldNdx (FloatMatrix&answer, int isurf, const FloatArray& lcoords, const FEICellGeometry& cellgeo);
    virtual void surfaceLocal2global(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void computeLocalSurfaceMapping(IntArray &nodes, int iSurf);
    virtual double evalNXIntegral(int iEdge, const FEICellGeometry &cellgeo);

    virtual void giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo);

    virtual IntegrationRule *giveIntegrationRule(int order);
    virtual IntegrationRule *giveBoundaryIntegrationRule(int order, int boundary);

protected:
    virtual void giveLocalDerivative(FloatMatrix &dN, const FloatArray &lcoords);
};
} // end namespace oofem
#endif
