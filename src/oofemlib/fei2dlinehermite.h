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

#ifndef fei2dlinehermite_h
#define fei2dlinehermite_h

#include "feinterpol2d.h"

namespace oofem {
/**
 * Class representing a 2d line with Hermitian interpolation.
 * The order used is cubic, quadratic, cubic, quadratic.
 * The functions that need geometric information, a linear interpolation is assumed (for geometry). This means functions such as evaldNdx.
 * @author Mikael Öhman
 */
class OOFEM_EXPORT FEI2dLineHermite : public FEInterpolation2d
{
public:
    FEI2dLineHermite(int ind1, int ind2) : FEInterpolation2d(1, ind1, ind2) { }

    integrationDomain giveIntegrationDomain() const override { return _Line; }
    Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }

    integrationDomain giveBoundaryIntegrationDomain(int ib) const override { return _Point; }
    integrationDomain giveBoundarySurfaceIntegrationDomain(int isurf) const override { return _UnknownIntegrationDomain; }
    integrationDomain giveBoundaryEdgeIntegrationDomain(int iedge) const override { return _UnknownIntegrationDomain; }


    double giveArea(const FEICellGeometry &cellgeo) const override { return 0.0; }
    double giveLength(const FEICellGeometry &cellgeo) const;

    void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    int global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo) const override;

    // "Bulk"
    void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;

    // Edge (same as bulk for this type, so they are all ignored) (perhaps do it the other way around?).
    IntArray computeLocalEdgeMapping(int iedge) const override { return {}; }
    void edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override { }
    double edgeEvalNormal(FloatArray &normal, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void edgeEvaldNds(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void edgeEvald2Nds2(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const ;

    void edgeLocal2global(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override { }

    int giveNumberOfNodes() const override { return 2; }

    std::unique_ptr<IntegrationRule> giveIntegrationRule(int order) const override;
};
} // end namespace oofem
#endif // fei2dlinehermite_h
