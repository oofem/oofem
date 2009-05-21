/* $Header: /home/cvs/bp/oofem/oofemlib/src/fei1dlin.h,v 1.1 2003/04/06 14:08:24 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

#ifndef fei1dlin_h
#define fei1dlin_h

#include "feinterpol1d.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"

/**
 * Class representing a 1d linear isparametric interpolation.
 */
class FEI1dLin : public FEInterpolation1d
{
protected:
    int cindx;

public:
    FEI1dLin(int coordIndx) : FEInterpolation1d(1) { cindx = coordIndx; }

    /**
     * Evaluates the array of interpolation functions (shape functions) at given point.
     * @param answer contains resulting array of evaluated interpolation functions
     * @param lcoords array containing (local) coordinates
     * @param cellgeo underlying cell geometry
     * @param time time
     */
    virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEIElementGeometry& cellgeo, double time);
    /**
     * Evaluates the matrix of derivatives of interpolation functions (shape functions) at given point.
     * These derivatives are in global coordinate system (where the nodal coordinates are defined)
     * @param matrix contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxj
     * @param lcoords array containing (local) coordinates
     * @param cellgeo underlying cell geometry
     * @param time time
     */
    virtual void evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEIElementGeometry& cellgeo, double time);
    /**
     * Evaluates global coordinates from given local ones
     * These derivatives are in global coordinate system (where the nodal coordinates are defined)
     * @param answer contains resulting global coordinates
     * @param lcoords array containing (local) coordinates
     * @param cellgeo underlying cell geometry
     * @param time time
     */
    virtual void local2global(FloatArray &answer, const FloatArray &lcoords, const FEIElementGeometry& cellgeo, double time);
    /**
     * Evaluates local coordinates from given global ones. Returns nonzero if local coordinates are interpolating,
     * zero if extrapolating (nonzero is returned if point is within the element geometry, zero otherwise).
     * These derivatives are in global coordinate system (where the nodal coordinates are defined)
     * @param answer contains evaluated local coordinates
     * @param lcoords array containing (local) coordinates
     * @param time time
     * @param cellgeo underlying cell geometry
     * @return nonzero is returned if point is within the element geometry, zero otherwise
     */
    virtual int  global2local(FloatArray &answer, const FloatArray &lcoords, const FEIElementGeometry& cellgeo, double time);
    /**
     * Evaluates the jacobian of transformation between local and global coordinates.
     */
    virtual double giveTransformationJacobian(const FloatArray &lcoords, const FEIElementGeometry& cellgeo, double time);

protected:
    double computeLength(const FEIElementGeometry& cellgeo);
};




#endif // fei1dlin_h






