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

#ifndef fei3dtriquad_h
#define fei3dtriquad_h

#include "feinterpol3d.h"

namespace oofem {
/**
 * Second order triangular interpolation in 3D space (6 nodes).
 * @todo This class is entirely unchecked.
 * @author Jim Brouzoulis
 * @author Mikael Öhman
 */
class OOFEM_EXPORT FEI3dTrQuad : public FEInterpolation3d
{
public:
    FEI3dTrQuad() : FEInterpolation3d(2) { }

    integrationDomain giveIntegrationDomain(const Element_Geometry_Type) const override { return _Triangle; }
    const Element_Geometry_Type giveGeometryType() const override { return EGT_triangle_2; }
    const Element_Geometry_Type giveBoundaryGeometryType(int boundary) const override { return EGT_line_2; } 


    integrationDomain giveBoundaryIntegrationDomain(int ib, const Element_Geometry_Type) const override { return _Triangle; }
    integrationDomain giveBoundarySurfaceIntegrationDomain(int isurf, const Element_Geometry_Type) const override { return _Triangle; }
    integrationDomain giveBoundaryEdgeIntegrationDomain(int iedge, const Element_Geometry_Type) const override { return _Line; }

    // Bulk
    void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void evaldNdxi(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    int global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;

    // new methods
    void giveDerivativeXi(FloatArray &n, const FloatArray &lcoords) const;
    void giveDerivativeEta(FloatArray &n, const FloatArray &lcoords) const ;
    void giveLocalNodeCoords(FloatMatrix &answer, const Element_Geometry_Type) const override;
    void surfaceEvaldNdxi(FloatMatrix &answer, const FloatArray &lcoords) const override;

    // Edge
    void edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void edgeEvaldNdxi(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void edgeEvaldNdx(FloatMatrix &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void edgeLocal2global(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    IntArray computeLocalEdgeMapping(int iedge) const override;

    int giveNumberOfEdges(const Element_Geometry_Type) const override { return 3; }

    // Surface
    void surfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void surfaceEvaldNdx(FloatMatrix &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void surfaceLocal2global(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    IntArray computeLocalSurfaceMapping(int iedge) const override;
    void surfaceEvalBaseVectorsAt(FloatArray &G1, FloatArray &G2, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const ;
    void surfaceGiveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const ;

    std::unique_ptr<IntegrationRule> giveIntegrationRule(int order, const Element_Geometry_Type) const override;
    std::unique_ptr<IntegrationRule> giveBoundaryIntegrationRule(int order, int boundary, const Element_Geometry_Type) const override;
    double giveArea(const FEICellGeometry &cellgeo) const;

    int giveNumberOfNodes(const Element_Geometry_Type) const override { return 6; }


    
protected:
    double edgeComputeLength(const IntArray &edgeNodes, const FEICellGeometry &cellgeo) const;
    double computeVolume(const FEICellGeometry &cellgeo) const;

};
} // end namespace oofem
#endif // fei3dtetquad_h
