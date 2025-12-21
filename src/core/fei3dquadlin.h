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

#ifndef fei3dquadlin_h
#define fei3dquadlin_h

#include "feinterpol3d.h"

namespace oofem {
/**
 * First order (linear) quad interpolation in 3D space (4 nodes).
 * @todo This class is entirely unchecked.
 */
class OOFEM_EXPORT FEI3dQuadLin : public FEInterpolation3d
{
public:
    FEI3dQuadLin() : FEInterpolation3d(2) { }

    integrationDomain giveIntegrationDomain(const Element_Geometry_Type) const override { return _Square; }
    const Element_Geometry_Type giveGeometryType() const override { return EGT_quad_1_interface; }
    const Element_Geometry_Type giveBoundaryGeometryType(int boundary) const override { return EGT_line_1; } 

    integrationDomain giveBoundaryIntegrationDomain(int ib, const Element_Geometry_Type) const override { return _Line; }
    integrationDomain giveBoundarySurfaceIntegrationDomain(int isurf, const Element_Geometry_Type) const override { return _Square; }
    integrationDomain giveBoundaryEdgeIntegrationDomain(int iedge, const Element_Geometry_Type) const override { return _Line; }

    // Bulk
    void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void evaldNdxi(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    int global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;

    // Edge
    void edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void edgeEvaldNdxi(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void edgeEvaldNdx(FloatMatrix &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void edgeLocal2global(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    IntArray computeLocalEdgeMapping(int iedge) const override;

    int giveNumberOfEdges(const Element_Geometry_Type) const override { return 4; }

    // Surface
    void surfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void surfaceEvaldNdx(FloatMatrix &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void surfaceEvaldNdxi(FloatMatrix &answer, const FloatArray &lcoords) const override;
    void surfaceEvald2Ndxi2(FloatMatrix &answer, const FloatArray &lcoords) const override;

    double surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void surfaceLocal2global(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    IntArray computeLocalSurfaceMapping(int iedge) const override;
    void surfaceEvalBaseVectorsAt(FloatArray &G1, FloatArray &G2, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const ;
    void surfaceGiveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const ;

    std::unique_ptr<IntegrationRule> giveIntegrationRule(int order, const Element_Geometry_Type) const override;
    std::unique_ptr<IntegrationRule> giveBoundaryIntegrationRule(int order, int boundary, const Element_Geometry_Type) const override;
    double giveArea(const FEICellGeometry &cellgeo) const;

    int giveNumberOfNodes(const Element_Geometry_Type) const override { return 4; }


    
protected:
    double edgeComputeLength(const IntArray &edgeNodes, const FEICellGeometry &cellgeo) const;
    double computeVolume(const FEICellGeometry &cellgeo) const ;

};
} // end namespace oofem
#endif // fei3dquadlin_h
