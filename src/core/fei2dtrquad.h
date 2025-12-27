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

#ifndef fei2dtrquad_h
#define fei2dtrquad_h

#include "feinterpol2d.h"

namespace oofem {
/**
 * Second order triangular interpolation in 2D (6 nodes).
 */
class OOFEM_EXPORT FEI2dTrQuad : public FEInterpolation2d
{
public:
    FEI2dTrQuad(int ind1, int ind2) : FEInterpolation2d(2, ind1, ind2) { }

    integrationDomain giveIntegrationDomain(const Element_Geometry_Type) const override { return _Triangle; }
    const Element_Geometry_Type giveGeometryType() const override { return EGT_triangle_2; }
    const Element_Geometry_Type giveBoundaryGeometryType(int boundary) const override { return EGT_line_2; } 

    integrationDomain giveBoundaryIntegrationDomain(int ib, const Element_Geometry_Type) const override { return _Line; }
    integrationDomain giveBoundarySurfaceIntegrationDomain(int isurf, const Element_Geometry_Type) const override { return _Triangle; }
    integrationDomain giveBoundaryEdgeIntegrationDomain(int iedge, const Element_Geometry_Type) const override { return _Line; }

    // Bulk
    static FloatArrayF<6> evalN(const FloatArrayF<2> &lcoords);
    static FloatMatrixF<2,6> evaldNdxi(const FloatArrayF<2> &lcoords);
    std::pair<double,FloatMatrixF<2,6>> evaldNdx(const FloatArrayF<2> &lcoords, const FEICellGeometry &cellgeo) const;

    void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void evald2Ndx2(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double giveArea(const FEICellGeometry &cellgeo) const override;
    int giveNumberOfNodes(const Element_Geometry_Type) const override { return 6; }

    bool inside(const FloatArray &lcoords) const override;

    // Edge
    IntArray computeLocalEdgeMapping(int iedge) const override;
    int giveNumberOfEdges(const Element_Geometry_Type) const override { return 3; }
    void edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double edgeEvalNormal(FloatArray &normal, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void edgeEvaldNds(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void edgeLocal2global(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double evalNXIntegral(int iEdge, const FEICellGeometry &cellgeo) const override;

    std::unique_ptr<IntegrationRule> giveIntegrationRule(int order, const Element_Geometry_Type) const override;

    void evaldNdxi(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;

protected:
    double edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo);
};
} // end namespace oofem
#endif // fei2dtrquad_h
