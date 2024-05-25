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

#ifndef fei3dhexatriquad_h
#define fei3dhexatriquad_h

#include "fei3dhexaquad.h"

namespace oofem {
/**
 * Class representing implementation of tri-quadratic hexahedra interpolation class.
 * @author Mikael Öhman
 */
class OOFEM_EXPORT FEI3dHexaTriQuad : public FEI3dHexaQuad
{
public:
    FEI3dHexaTriQuad() : FEI3dHexaQuad() { }

    integrationDomain giveIntegrationDomain(const Element_Geometry_Type) const override { return _Cube; }
    const Element_Geometry_Type giveGeometryType() const override { return EGT_hexa_27; }
    integrationDomain giveBoundaryIntegrationDomain(int ib, const Element_Geometry_Type) const override { return _Square; }
    integrationDomain giveBoundarySurfaceIntegrationDomain(int isurf, const Element_Geometry_Type) const override { return _Square; }
    integrationDomain giveBoundaryEdgeIntegrationDomain(int iedge, const Element_Geometry_Type) const override { return _Line; }

    // Bulk
    static FloatArrayF<27> evalN(const FloatArrayF<3> &lcoords);
    static std::pair<double, FloatMatrixF<3,27>> evaldNdx(const FloatArrayF<3> &lcoords, const FEICellGeometry &cellgeo);
    static FloatMatrixF<3,27> evaldNdxi(const FloatArrayF<3> &lcoords);

    void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void evaldNdxi(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    int giveNumberOfNodes(const Element_Geometry_Type) const override { return 27; }

    // Surface
    void surfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    //void surfaceEvaldNdx(FloatMatrix&answer, int isurf, const FloatArray& lcoords, const FEICellGeometry& cellgeo) override;
    IntArray computeLocalSurfaceMapping(int iSurf) const override;
    double evalNXIntegral(int iSurf, const FEICellGeometry &cellgeo) const override;

    std::unique_ptr<IntegrationRule> giveIntegrationRule(int order, const Element_Geometry_Type) const override;
    std::unique_ptr<IntegrationRule> giveBoundaryIntegrationRule(int order, int boundary, const Element_Geometry_Type) const override;
};
} // end namespace oofem
#endif
