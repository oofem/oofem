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

#ifndef fei2dtrquad_h
#define fei2dtrquad_h

#include "feinterpol2d.h"

namespace oofem {
/**
 * Second order triangular interpolation in 2D (6 nodes).
 */
class OOFEM_EXPORT FEI2dTrQuad : public FEInterpolation2d
{
public:
    FEI2dTrQuad(int ind1, int ind2) : FEInterpolation2d(2, ind1, ind2) { }

    virtual integrationDomain giveIntegrationDomain() const { return _Triangle; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_triangle_2; }

    // Bulk
    virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void evald2Ndx2(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual int  global2local(FloatArray &answer, const FloatArray &gcoords, const FEICellGeometry &cellgeo);
    virtual double giveArea(const FEICellGeometry &cellgeo) const;
    virtual void giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual int giveNumberOfNodes() const { return 6; }    
    /**
     * Returns a characteristic length of the geometry, typically a diagonal or edge length.
     * @param cellgeo Underlying cell geometry.
     * @return Square root of area.
     */
    virtual double giveCharacteristicLength(const FEICellGeometry &cellgeo) const;

    // Edge
    virtual void computeLocalEdgeMapping(IntArray &edgeNodes, int iedge);
    virtual int giveNumberOfEdges() const { return 3; }
    virtual void edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double edgeEvalNormal(FloatArray &normal, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void edgeEvaldNds(FloatArray &answer, int iedge,
                              const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void edgeLocal2global(FloatArray &answer, int iedge,
                                  const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double evalNXIntegral(int iEdge, const FEICellGeometry &cellgeo);

    virtual IntegrationRule *giveIntegrationRule(int order);

protected:
    double edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo);
    void giveDerivatives(FloatMatrix &dn, const FloatArray &lc);
};
} // end namespace oofem
#endif // fei2dtrquad_h
