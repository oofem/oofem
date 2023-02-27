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

#ifndef feitspline_h
#define feitspline_h

#include "feibspline.h"
#include "floatarray.h"
#include "floatmatrix.h"

/*
 * alternatively, it is possible to store for individual control points open local knot vector;
 * however, this is not enough as I need to know how many knots have been prepended and appended
 * (see createLocalKnotVector) as these are relevant for finding relevant knot span using BSpline method
 * (this could be overcome by writing corresponding TSpline method) and for extraction proper
 * basis function and its derivatives (computed by BSpline methods)
 */

///@name Input fields for TSplineInterpolation
//@{
#define _IFT_TSplineInterpolation_weights "weights"
#define _IFT_TSplineInterpolation_localIndexKnotVectorU "localindexknotvectoru"
#define _IFT_TSplineInterpolation_localIndexKnotVectorV "localindexknotvectorv"
#define _IFT_TSplineInterpolation_localIndexKnotVectorW "localindexknotvectorw"
//@}

namespace oofem {
/**
 * Interpolation for T-splines.
 */
class OOFEM_EXPORT TSplineInterpolation : public BSplineInterpolation
{
protected:
    /// Local index knot vector of the dimensions [totalNumberOfControlPoints][nsd][degree+2].
    std::vector< std::array<IntArray, 3> > localIndexKnotVector;
    int totalNumberOfControlPoints = 0;
    FloatArray weights;
    /**
     * Temporary open local knot vector to enable use of BSpline algorithms (common for all directions) [3*max_degree+2].
     */
    //FloatArray openLocalKnotVector;
public:
    TSplineInterpolation(int nsd) : BSplineInterpolation(nsd) { }

    void initializeFrom(InputRecord &ir) override;
    void setNumberOfControlPoints(int num) { this->totalNumberOfControlPoints = num; }
    void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;
    int global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override {
        OOFEM_ERROR("Not yet implemented, contact lazy dr for implementation");
        return 0;
    }
    void giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo) const override;

    int giveKnotSpanBasisFuncMask(const IntArray &knotSpan, IntArray &mask) const override;
    int giveNumberOfKnotSpanBasisFunctions(const IntArray &knotSpan) const override;

    const char *giveClassName() const { return "TSplineInterpolation"; }

protected:
    /**
     * Evaluates the middle basis function on local knot vector at u.
     * @param u Value at which to evaluate.
     * @param p Degree.
     * @param U Global knot values.
     * @param I Local index knot vector.
     * @return N Computed middle basis function.
     * @warning Parameter u must be in a valid range.
     */
    double basisFunction(double u, int p, const FloatArray &U, const IntArray &I) const;
    /**
     * Computes the middle basis function and it derivatives on local knot vector at u.
     * The result is stored in the ders vector
     * @f{align*}{
     * N(u)  &= ders(0)
     * N'(u) &= ders(1)
     * N''(u)&= ders(2)
     * @f}
     * @param n Degree of the derivation.
     * @param u Parametric value.
     * @param p Degree.
     * @param U Global knot values.
     * @param I Local index knot vector.
     * @param ders Vector containing the derivatives of the middle basis function.
     *
     * @warning Parameters n and u must be in a valid range.
     */
    void dersBasisFunction(int n, double u, int p, const FloatArray &U, const IntArray &I, FloatArray &ders) const;
    /**
     * Creates local open knot vector.
     * This is generally done extracting knot values from global knot vector using the local index knot vector
     * and by prepending p times the first knot and appending p times the last knot.
     * However, existing knot multiplicity at the start and end must be accounted for.
     * @param p Degree.
     * @param U Global knot values.
     * @param I Local index knot vector
     * @param prepend Number of prepended entries
     * @param append Number of appended entries
     */
    void createLocalKnotVector(FloatArray& answer, int p, const FloatArray &U, const IntArray &I, int &prepend, int &append) const;
    /**
     * Returns indices (zero based) of nonzero basis functions for given knot span interval.
     */
    int giveKnotSpanBasisFuncMask(const IntArray &startKnotSpan, const IntArray &endKnotSpan, IntArray &mask) const;
    /** Returns the number of nonzero basis functions at given knot span interval. */
    int  giveNumberOfKnotSpanBasisFunctions(const IntArray &startKnotSpan, const IntArray &endKnotSpan) const;
}; // end of TSplineInterpolation class definition
} // end namespace oofem
#endif //feitspline_h
