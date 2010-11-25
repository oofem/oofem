/* $Header: /home/cvs/bp/oofem/oofemlib/src/element.h,v 1.27 2003/04/06 14:08:24 bp Exp $ */
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

#ifndef feitspline_h
#define feitspline_h

/*
 * oofem nodes - control points (coordinates ) + dofs
 * oofem elements - NURBS patches as well as integration elements
 *
 *
 * NURBS PATCH:
 * knot vector - store knot coordinates + multiplicity
 * patch integration rule - keep list of elements
 *
 */

#include "feibspline.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "mathfem.h"

namespace oofem {
/* alternatively, it is possible to store for individual control points open local knot vector;
 *       however, this is not enough as I need to know how many knots have been prepended and appended
 *       (see createLocalKnotVector) as these are relevant for finding relevant knot span using BSpline method
 *       (this could be overcome by writing corresponding TSpline method) and for extraction proper
 *       basis function and its derivatives (computed by BSpline methods) */

class TSplineInterpolation : public BSplineInterpolation
{
protected:
    /// localIndexKnotVector[number_of_control_points][nsd][degree+2]
    int ***localIndexKnotVector;
    int totalNumberOfControlPoints;
    /// temporary open local knot vector to enable use of BSpline algorithms (common for all directions)
    /// openLocalKnotVector[3*max_degree+2]
    double *openLocalKnotVector;
public:
    TSplineInterpolation(int nsd) : BSplineInterpolation(nsd) { }
    ~TSplineInterpolation();

    IRResultType initializeFrom(InputRecord *ir);
    void setNumberOfControlPoints(int num) { this->totalNumberOfControlPoints = num; }
    /**
     * Evaluates the array of interpolation functions (shape functions) at given point.
     * @param answer contains resulting array of evaluated interpolation functions
     * @param lcoords array containing (local) coordinates
     * @param cellgeo underlying cell geometry
     * @param time time
     *
     * see also giveNonzeroBasisFunctMask method of TSplineInterpolation
     */
    virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time);
    /**
     * Evaluates the matrix of derivatives of interpolation functions (shape functions) at given point.
     * These derivatives are in global coordinate system (where the nodal coordinates are defined)
     * @param answer contains resulting matrix of derivatives, the member at i,j position contains value of dNi/dxj
     * @param lcoords array containing (local) coordinates
     * @param cellgeo underlying cell geometry
     * @param time time
     */
    virtual void evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time);
    /**
     * Evaluates global coordinates from given local ones
     * @param answer contains resulting global coordinates
     * @param lcoords array containing (local) coordinates
     * @param cellgeo underlying cell geometry
     * @param time time
     */
    virtual void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time);
    virtual int  global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time) {
        OOFEM_ERROR("Not yet inplemented, contact lazy dr for implementation");
        return 0;
    }
    /**
     * Evaluates the jacobian of transformation between local and global coordinates.
     */
    virtual double giveTransformationJacobian(const FloatArray &lcoords, const FEICellGeometry &cellgeo, double time);

    /**
     * Returns indices (zero based) of nonzero basis functions for given knot span
     */
    int giveKnotSpanBasisFuncMask(const IntArray &knotSpan, IntArray &mask);
    /** Returns the number of nonzero basis functions for given knot span */
    int giveNumberOfKnotSpanBasisFunctions(const IntArray &knotSpan);

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "TSplineInterpolation"; }
protected:
    /**
     * Evaluates the middle basis function on local knot vector at u
     * @param u value at which to evaluate
     * @param p degree
     * @param U global knot values
     *             @param I local index knot vector
     * @return N computed middle basis function
     *
     *             @warning u must be in a valid range.
     */
    double basisFunction(double u, int p, const FloatArray &U, const int *I);
    /*
     * Computes the middle basis function and it derivatives on local knot vector at u
     *
     * The result is stored in the ders vector
     *            N(u)  = ders(0)
     * N'(u) = ders(1)
     *            N''(u)= ders(2)
     *
     * @param n the degree of the derivation
     * @param u the parametric value
     * @param p the degree
     *            @param U global knot values
     *            @param I local index knot vector
     * @param ders vector containing the derivatives of the middle basis function.
     *
     * @warning n and u must be in a valid range.
     */
    void dersBasisFunction(int n, double u, int p, const FloatArray &U, const int *I, FloatArray &ders);
    /*
     * Creates local open knot vector.
     *            This is generally done extracting knot values from global knot vector using the local index knot vector
     *            and by prepending p times the first knot and appending p times the last knot.
     *            However, existing knot multiplicity at the start and end must be accounted for.
     *
     * @param p the degree
     * @param U global knot values
     * @param I local index knot vector
     * @param prepent number of prepended entries
     *            @param append number of appended entries
     */
    void createLocalKnotVector(int p, const FloatArray &U, const int *I, int *prepend, int *append);
    /**
     * Returns indices (zero based) of nonzero basis functions for given knot span interval (from start to end)
     */
    int giveKnotSpanBasisFuncMask(const IntArray &startKnotSpan, const IntArray &endKnotSpan, IntArray &mask);
    /** Returns the number of nonzero basis functions at given knot span interval (from start to end) */
    int  giveNumberOfKnotSpanBasisFunctions(const IntArray &startKnotSpan, const IntArray &endKnotSpan);
}; // end of TSplineInterpolation class definition
} // end namespace oofem
#endif //feitspline_h
