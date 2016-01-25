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

    virtual integrationDomain giveIntegrationDomain() const { return _Square; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_quad_2; }

    virtual double giveArea(const FEICellGeometry &cellgeo) const;

    // Bulk
    virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual int  global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo);
    virtual void giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual int giveNumberOfNodes() const { return 8; } 
    /**
     * Returns a characteristic length of the geometry, typically a diagonal or edge length.
     * @param cellgeo underlying cell geometry
     * @return distance from node 1 to 3
     */
    virtual double giveCharacteristicLength(const FEICellGeometry &cellgeo) const;

    // Edge
    virtual void computeLocalEdgeMapping(IntArray &edgeNodes, int iedge);
    virtual void edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void edgeEvaldNds(FloatArray &answer, int iedge,
                              const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void edgeLocal2global(FloatArray &answer, int iedge,
                                  const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double edgeEvalNormal(FloatArray &normal, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double evalNXIntegral(int iEdge, const FEICellGeometry &cellgeo);

    virtual IntegrationRule *giveIntegrationRule(int order);

protected:
    virtual void giveDerivatives(FloatMatrix &dn, const FloatArray &lc);
};
} // end namespace oofem
#endif // fei2dquadquad_h
