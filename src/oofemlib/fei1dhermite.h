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

#ifndef fei1dquad_h
#define fei1dquad_h

#include "feinterpol1d.h"

namespace oofem {
/**
 * Class representing a 1d Hermitian cubic isoparametric interpolation.
 * @author Mikael Öhman
 */
class OOFEM_EXPORT FEI1dHermite : public FEInterpolation1d
{
protected:
    int cindx;

public:
    FEI1dHermite(int coordIndx) : FEInterpolation1d(2), cindx(coordIndx) { }

    integrationDomain giveIntegrationDomain(const Element_Geometry_Type) const override { return _Line; }
    const Element_Geometry_Type giveGeometryType() const override { return EGT_line_1; }

    double giveLength(const FEICellGeometry &cellgeo) const override;

    std::pair<double, FloatMatrixF<1,4>> evaldNdx(double ksi, const FEICellGeometry &cellgeo) const;
    FloatMatrixF<1,4> evald2Ndx2(double ksi, const FEICellGeometry &cellgeo) const;

    void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void evald2Ndx2(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    int  global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;

    int giveNumberOfNodes(const Element_Geometry_Type) const override { return 2; }
};
} // end namespace oofem
#endif
