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
protected:

public:
    FEI3dWedgeLin() : FEInterpolation3d(1) { }

    virtual integrationDomain giveIntegrationDomain() const { return _Wedge; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_wedge_1; }

    // Bulk
    virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual int global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double giveCharacteristicLength(const FEICellGeometry &cellgeo) const;

    // Edge
    virtual void edgeEvalN(FloatArray &answer, int iedge, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void edgeEvaldNdx(FloatMatrix &answer, int iedge,
                              const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void edgeLocal2global(FloatArray &answer, int iedge,
                                  const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords,
                                                  const FEICellGeometry &cellgeo);
    virtual void computeLocalEdgeMapping(IntArray &edgeNodes, int iedge);

    // Surface
    virtual void surfaceEvalN(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    //virtual void surfaceEvaldNdx (FloatMatrix&answer, int iedge,
    //               const FloatArray& lcoords, const FEICellGeometry& cellgeo);
    virtual void surfaceLocal2global(FloatArray &answer, int iedge,
                                     const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords,
                                                     const FEICellGeometry &cellgeo);
    virtual void computeLocalSurfaceMapping(IntArray &nodes, int iSurf);

    virtual IntegrationRule *giveIntegrationRule(int order);
    virtual IntegrationRule *giveBoundaryIntegrationRule(int order, int boundary);
    virtual IntegrationRule *giveSurfaceIntegrationRule(int order, int isurf)
        { return giveBoundaryIntegrationRule(order, isurf); }

    virtual int giveNumberOfNodes() const { return 6; }

protected:
    double edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo);
    void giveLocalDerivative(FloatMatrix &dN, const FloatArray &lcoords);
};
} // end namespace oofem
#endif
