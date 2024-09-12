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

#ifndef fei1dlin_h
#define fei1dlin_h

#include "feinterpol1d.h"

namespace oofem {
/**
 * Class representing a 1d linear isoparametric interpolation.
 */
class OOFEM_EXPORT FEI1dLin : public FEInterpolation1d
{
protected:
    int cindx = 0;

public:
    FEI1dLin(int coordIndx) : FEInterpolation1d(1), cindx(coordIndx) { }

    integrationDomain giveIntegrationDomain(const Element_Geometry_Type) const override { return _Line; }
    const Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }
    integrationDomain giveBoundaryIntegrationDomain(int ib,const Element_Geometry_Type) const override { return _Point; }
    integrationDomain giveBoundarySurfaceIntegrationDomain(int isurf, const Element_Geometry_Type) const override { return _UnknownIntegrationDomain; }
    integrationDomain giveBoundaryEdgeIntegrationDomain(int iedge, const Element_Geometry_Type) const override { return _UnknownIntegrationDomain; }

    double giveLength(const FEICellGeometry &cellgeo) const override;

    static FloatArrayF<2> evalN(double ksi);
    std::pair<double, FloatMatrixF<1,2>> evaldNdx(const FEICellGeometry &cellgeo) const;

    void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    int  global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;

    int giveNumberOfNodes(Element_Geometry_Type) const override { return 2; }

    IntArray boundaryEdgeGiveNodes(int boundary,Element_Geometry_Type) const override;
    void boundaryEdgeEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double boundaryEdgeGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void boundaryEdgeLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
};
} // end namespace oofem
#endif // fei1dlin_h
