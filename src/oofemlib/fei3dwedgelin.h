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


#ifndef fei3dwedgelin_h
#define fei3dwedgelin_h

#include "feinterpol3d.h"

namespace oofem {
/**
 * Class representing implementation of linear wedge interpolation class.
 *
 * @author Milan Jirasek
 * @author Mikael Öhman
 */
class OOFEM_EXPORT FEI3dWedgeLin : public FEInterpolation3d
{
public:
    FEI3dWedgeLin() : FEInterpolation3d(1) { }

    integrationDomain giveIntegrationDomain(const Element_Geometry_Type) const override { return _Wedge; }
    const Element_Geometry_Type giveGeometryType() const override { return EGT_wedge_1; }
    integrationDomain giveBoundaryIntegrationDomain(int ib, const Element_Geometry_Type) const override
    {
        if (ib <= 2) return _Triangle;
        else return _Square;
    }
    integrationDomain giveBoundarySurfaceIntegrationDomain(int isurf, const Element_Geometry_Type egt) const override { return this->giveBoundaryIntegrationDomain(isurf, egt); }
    integrationDomain giveBoundaryEdgeIntegrationDomain(int iedge, const Element_Geometry_Type) const override { return _Line; }

    // Bulk
    static FloatArrayF<6> evalN(const FloatArrayF<3> &lcoords);
    static std::pair<double, FloatMatrixF<3,6>> evaldNdx(const FloatArrayF<3> &lcoords, const FEICellGeometry &cellgeo);
    static FloatMatrixF<3,6> evaldNdxi(const FloatArrayF<3> &lcoords);

    void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void evaldNdxi(FloatMatrix & answer, const FloatArray & lcoords, const FEICellGeometry & cellgeo) const override;
    void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    int global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double giveCharacteristicLength(const FEICellGeometry &cellgeo) const;

    // Edge
    void edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void edgeEvaldNdx(FloatMatrix &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void edgeLocal2global(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo)  const override;
    double edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    IntArray computeLocalEdgeMapping(int iedge) const override;

    // Surface
    void surfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    //void surfaceEvaldNdx (FloatMatrix&answer, int iedge, const FloatArray& lcoords, const FEICellGeometry& cellgeo) override;
    void surfaceLocal2global(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    IntArray computeLocalSurfaceMapping(int iSurf) const override;

    std::unique_ptr<IntegrationRule> giveIntegrationRule(int order, const Element_Geometry_Type) const override;
    std::unique_ptr<IntegrationRule> giveBoundaryIntegrationRule(int order, int boundary, const Element_Geometry_Type) const override;
    std::unique_ptr<IntegrationRule> giveSurfaceIntegrationRule(int order, int isurf, const Element_Geometry_Type egt) const
    { return giveBoundaryIntegrationRule(order, isurf, egt); }

    int giveNumberOfNodes(const Element_Geometry_Type) const override { return 6; }

protected:
    double edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo) const;
};
} // end namespace oofem
#endif
