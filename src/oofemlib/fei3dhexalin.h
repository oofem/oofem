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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef fei3dhexalin_h
#define fei3dhexalin_h

#include "feinterpol3d.h"

namespace oofem {
/**
 * Class representing implementation of linear hexahedra interpolation class.
 */
class FEI3dHexaLin : public FEInterpolation3d
{
public:
    FEI3dHexaLin() : FEInterpolation3d(1) { }

    // Bulk
    virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual int  global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo);

    // Edge
    virtual void edgeEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void edgeEvaldNdx(FloatMatrix &answer, int iedge,
                              const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void edgeLocal2global(FloatArray &answer, int iedge,
                                  const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double edgeGiveTransformationJacobian(int iedge, const FloatArray &lcoords,
                                                  const FEICellGeometry &cellgeo);
    virtual void computeLocalEdgeMapping(IntArray &edgeNodes, int iedge);

    // Surface
    virtual void surfaceEvalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    //virtual void surfaceEvaldNdx (FloatMatrix&answer, int isurf,
    //               Domain* d, IntArray& nodes, const FloatArray& lcoords);
    virtual void surfaceLocal2global(FloatArray &answer, int isurf,
                                     const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double surfaceEvalNormal(FloatArray &answer, int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double surfaceGiveTransformationJacobian(int isurf, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void computeLocalSurfaceMapping(IntArray &edgeNodes, int iedge);

    virtual void giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo);

protected:
    double edgeComputeLength(IntArray &edgeNodes, const FEICellGeometry &cellgeo);
    void giveLocalDerivative(FloatMatrix &dN, const FloatArray &lcoords);
};
} // end namespace oofem
#endif //  fei3dhexalin_h






