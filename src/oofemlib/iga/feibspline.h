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

#ifndef feibspline_h
#define feibspline_h

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

#include "flotarry.h"
#include "flotmtrx.h"
#include "feinterpol.h"


class BSplineInterpolation : public FEInterpolation
{
protected:
    /// number of spatial directions
    int nsd;
    /// degree in each direction
    int *degree;                                   // eg. 2
    /// knotValues[nsd]
    FloatArray *knotValues;                        // eg. 0 1 2 3 4 5
    /// knotMultiplicity[nsd]
    IntArray *knotMultiplicity;                    // eg. 3 1 1 1 2 3
    /// numberOfControlPoints[nsd]
    /// for TSpline this is filled by values corresponding to case when there are no T-junctions
    /// (i.e. TSpline = BSpline)
    int *numberOfControlPoints;
    /// knotVectors[nsd][knot_vector_size]
    double **knotVector;                           // eg. 0 0 0 1 2 3 4 4 5 5 5
    /// nonzero spans in each directions[nsd]
    int *numberOfKnotSpans;                        // eg. 5 (0-1,1-2,2-3,3-4,4-5)
public:
    BSplineInterpolation(int nsd) : FEInterpolation(0) { this->nsd = nsd; }
    ~BSplineInterpolation();
    /**
     * Returns number of spatial dimensions
     */
    int const giveNsd() { return nsd; }
    IRResultType initializeFrom(InputRecord *ir);
    virtual int giveNumberOfKnotSpans(int dim) { return numberOfKnotSpans [ dim - 1 ]; }
    virtual int giveNumberOfControlPoints(int dim) { return numberOfControlPoints [ dim - 1 ]; }
    virtual double **const giveKnotVector() { return this->knotVector; }
    virtual IntArray *const giveKnotMultiplicity(int dim) { return & this->knotMultiplicity [ dim - 1 ]; }
    virtual FloatArray *const giveKnotValues(int dim) { return & this->knotValues [ dim - 1 ]; }
    /**
     * Evaluates the array of interpolation functions (shape functions) at given point.
     * @param answer contains resulting array of evaluated interpolation functions
     * @param lcoords array containing (local) coordinates
     * @param cellgeo underlying cell geometry
     * @param time time
     *
     * see also giveNonzeroBasisFunctMask method
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
    const char *giveClassName() const { return "BSplineInterpolation"; }
    virtual bool hasSubPatchFormulation() { return true; }
protected:
    /**
     * Evaluates the nonvanishing basis functions of 1d BSpline (algorithm A2.2 from NURBS book)
     * @param span knot span index (zero based)
     * @param u value at which to evaluate
     * @param p degree
     * @param U knot vector
     * @param N computed p+1 nonvanishing functions (N_{span-p,p}-N_{span,p})
     *
     *             @warning u and span must be in a valid range.
     */
    void basisFuns(FloatArray &N, int span, double u, int p, const double *U);
    /**
     * Computes nonzero basis functions and their derivatives at u
     *
     * For information on the algorithm, see A2.3 on p72 of
     * the NURBS book. The result is stored in the ders matrix, where
     * ders is of size (n+1,p+1) and the derivative
     *            N(u)  = ders(0,span-p+j) where j=0...p
     * N'(u) = ders(1,span-p+j) where j=0...p
     *            N''(u)= ders(2,span-p+j) where j=0...p
     *
     * @param n the degree of the derivation
     * @param u the parametric value
     *            @param span knot span index (zero based)
     * @param p the degree
     * @param U knot vector
     * @param ders matrix containing the derivatives of the basis functions.
     *
     * @warning n, u and span must be in a valid range.
     */
    void dersBasisFuns(int n, double u, int span, int p, double *const U, FloatMatrix &ders);
    /**
     * Determines the knot span index (Algorithm A2.1 from the NURBS book)
     *
     * Determines the knot span for which there exists non-zero basis
     * functions. The span is the index k for which the parameter
     * u is valid in the (u_k,u_{k+1}] range.
     *
     * @param n number of control points - 1 (number of ctrl pnts = n + 1)
     * @param u the parametric value
     * @param p the degree
     * @param U knot vector
     * @return the span index at u (zero based)
     * @warning u must be in a valid range
     */
    int findSpan(int n, int p, double u, const double *U) const;
    /** Return the range of nonzero basis functions for given knot span and given degree */
    void giveNonzeroBasisFuncInterval(int span, int deg, int &s, int &e) { s = span - deg;
                                                                           e = span; }
}; // enf of BSplineInterpolation class definition


#endif //feibspline_h
