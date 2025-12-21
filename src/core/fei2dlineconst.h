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

#ifndef fei2dlineconst_h
#define fei2dlineconst_h

#include "feinterpol2d.h"

namespace oofem {
/**
 * Class representing a 2d line with linear interpolation.
 * @todo{Some more routines to be implemented here}
 * @author Mikael Öhman
 */
class OOFEM_EXPORT FEI2dLineConst : public FEInterpolation2d
{
public:
    FEI2dLineConst(int ind1, int ind2) : FEInterpolation2d(0, ind1, ind2) { }

    integrationDomain giveIntegrationDomain(const Element_Geometry_Type) const override { return _Line; }
    const Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }
    const Element_Geometry_Type giveBoundaryGeometryType(int boundary) const override { return EGT_point; } 

    integrationDomain giveBoundaryIntegrationDomain(int ib, const Element_Geometry_Type) const override { return _Point; }
    integrationDomain giveBoundarySurfaceIntegrationDomain(int isurf, const Element_Geometry_Type) const override { return _UnknownIntegrationDomain; }
    integrationDomain giveBoundaryEdgeIntegrationDomain(int iedge, const Element_Geometry_Type) const override { return _Line; }

    double giveArea(const FEICellGeometry &cellgeo) const override { return 0.0; }

    void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    int global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo) const override;

    int giveNumberOfEdges(const Element_Geometry_Type) const override { return 1; }
    void giveCellDofMans(IntArray& nodes, IntArray& internalDofMans, Element* elem) const override {nodes={};internalDofMans={1};}


    // "Bulk"
    void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;

    // Edge (same as bulk for this type, so they are all ignored) (perhaps do it the other way around?).
    IntArray boundaryEdgeGiveNodes(int boundary, const Element_Geometry_Type, bool includeHierarchical=false) const override;
    IntArray computeLocalEdgeMapping(int iedge) const override;
    void edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double edgeEvalNormal(FloatArray &normal, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void edgeEvaldNds(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void edgeLocal2global(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override { }
    double evalNXIntegral(int iEdge, const FEICellGeometry &cellgeo) const override;

    int giveNumberOfNodes(const Element_Geometry_Type) const override { return 2; }

    std::unique_ptr<IntegrationRule> giveIntegrationRule(int order, const Element_Geometry_Type) const override;

protected:
    double edgeComputeLength(const IntArray &edgeNodes, const FEICellGeometry &cellgeo) const;
};
} // end namespace oofem
#endif // fei2dlineconst_h
