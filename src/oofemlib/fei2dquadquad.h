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

#ifndef fei2dquadquad_h
#define fei2dquadquad_h

#include "feinterpol2d.h"

namespace oofem {
/**
 * Class representing a 2d quadrilateral with quadratic interpolation based on isoparametric coordinates.
 * Local Node Numbering
 * (4)--(7)--(3)
 *  |         |
 *  |         |
 * (8)       (6)
 *  |         |
 *  |         |
 * (1)--(5)--(2)
 */
class OOFEM_EXPORT FEI2dQuadQuad : public FEInterpolation2d
{
public:
    FEI2dQuadQuad(int ind1, int ind2) : FEInterpolation2d(2, ind1, ind2) { }

    integrationDomain giveIntegrationDomain(const Element_Geometry_Type) const override { return _Square; }
    const Element_Geometry_Type giveGeometryType() const override { return EGT_quad_2; }
    integrationDomain giveBoundaryIntegrationDomain(int ib, const Element_Geometry_Type) const override { return _Line; }
    integrationDomain giveBoundarySurfaceIntegrationDomain(int isurf, const Element_Geometry_Type) const override { return _Square; }
    integrationDomain giveBoundaryEdgeIntegrationDomain(int iedge, const Element_Geometry_Type) const override { return _Line; }

    double giveArea(const FEICellGeometry &cellgeo) const override;

    // Bulk
    static FloatArrayF<8> evalN(const FloatArrayF<2> &lcoords) ;
    static FloatMatrixF<2,8> evaldNdxi(const FloatArrayF<2> &lcoords) ;
    std::pair<double, FloatMatrixF<2,8>> evaldNdx(const FloatArrayF<2> &lcoords, const FEICellGeometry &cellgeo) const;

    void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void evaldNdxi(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    int giveNumberOfNodes(const Element_Geometry_Type) const override { return 8; } 

    double giveCharacteristicLength(const FEICellGeometry &cellgeo) const override;

    bool inside(const FloatArray &lcoords) const override;

    // Edge
    IntArray computeLocalEdgeMapping(int iedge) const override;
    void edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void edgeEvaldNds(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void edgeLocal2global(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double edgeEvalNormal(FloatArray &normal, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double evalNXIntegral(int iEdge, const FEICellGeometry &cellgeo) const override;

    std::unique_ptr<IntegrationRule> giveIntegrationRule(int order, const Element_Geometry_Type) const override;
};


/**
 * Class representing a 2d isoparametric quadratic interpolation based on natural coordinates
 * for quadrilateral elements in axisymmetric setting.
 */
class OOFEM_EXPORT FEI2dQuadQuadAxi : public FEI2dQuadQuad
{
public:
    FEI2dQuadQuadAxi(int ind1, int ind2) : FEI2dQuadQuad(ind1, ind2) { }

    double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double boundaryEdgeGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double boundaryGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
};

} // end namespace oofem
#endif // fei2dquadquad_h
