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

#ifndef feibspline_h
#define feibspline_h

#include "feinterpol.h"
#include "floatarray.h"

///@name Input fields for BSplineInterpolation
//@{
#define _IFT_BSplineInterpolation_degree "degree"
#define _IFT_BSplineInterpolation_knotVectorU "knotvectoru"
#define _IFT_BSplineInterpolation_knotVectorV "knotvectorv"
#define _IFT_BSplineInterpolation_knotVectorW "knotvectorw"
#define _IFT_BSplineInterpolation_knotMultiplicityU "knotmultiplicityu"
#define _IFT_BSplineInterpolation_knotMultiplicityV "knotmultiplicityv"
#define _IFT_BSplineInterpolation_knotMultiplicityW "knotmultiplicityw"
//@}

namespace oofem {
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * Interpolation for B-splines.
 */
class OOFEM_EXPORT BSplineInterpolation : public FEInterpolation
{
protected:
    /// Number of spatial directions.
    int nsd;
    /// Degree in each direction.
    int *degree;                                   // eg. 2
    /// Knot values [nsd]
    FloatArray *knotValues;                        // eg. 0 1 2 3 4 5
    /// Knot multiplicity [nsd]
    IntArray *knotMultiplicity;                    // eg. 3 1 1 1 2 3
    /// numberOfControlPoints[nsd]
    /**
     * For TSpline this is filled by values corresponding to case when there are no T-junctions
     * (i.e. TSpline = BSpline)
     */
    int *numberOfControlPoints;
    /// Knot vectors [nsd][knot_vector_size]
    double **knotVector;                           // eg. 0 0 0 1 2 3 4 4 5 5 5
    /// Nonzero spans in each directions [nsd]
    int *numberOfKnotSpans;                        // eg. 5 (0-1,1-2,2-3,3-4,4-5)
public:
    BSplineInterpolation(int nsd) : FEInterpolation(0) {
        this->nsd = nsd;
    }
    virtual ~BSplineInterpolation();

    virtual integrationDomain giveIntegrationDomain() const {
        if ( nsd == 3 ) {
            return _Cube;
        } else if ( nsd == 2 ) {
            return _Square;
        } else if ( nsd == 1 ) {
            return _Line;
        } else {
            return _Unknown_integrationDomain;
        }
    }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_unknown; }

    virtual int giveNsd() { return nsd; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void boundaryEdgeGiveNodes(IntArray &answer, int boundary)
    { OOFEM_ERROR("Functions not supported for this interpolator."); }
    virtual void boundaryEdgeEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
    { OOFEM_ERROR("Functions not supported for this interpolator."); }
    virtual double boundaryEdgeGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
    { OOFEM_ERROR("Functions not supported for this interpolator.");
      return 0.; }
    virtual void boundaryEdgeLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
    { OOFEM_ERROR("Functions not supported for this interpolator."); }

    virtual void boundaryGiveNodes(IntArray &answer, int boundary)
    { OOFEM_ERROR("Not implemented"); }
    virtual void boundaryEvalN(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
    { OOFEM_ERROR("Not implemented"); }
    virtual double boundaryEvalNormal(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
    { OOFEM_ERROR("Not implemented");
      return 0.; }
    virtual double boundaryGiveTransformationJacobian(int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
    { OOFEM_ERROR("boundaryGiveTransformationJacobian - Not implemented");
      return 0.; }
    virtual void boundaryLocal2Global(FloatArray &answer, int boundary, const FloatArray &lcoords, const FEICellGeometry &cellgeo)
    { OOFEM_ERROR("boundaryLocal2Global - Not implemented"); }


    virtual int giveNumberOfKnotSpans(int dim) { return numberOfKnotSpans [ dim - 1 ]; }
    virtual int giveNumberOfControlPoints(int dim) { return numberOfControlPoints [ dim - 1 ]; }
    virtual const double *const *giveKnotVector() {
        return this->knotVector;
    }
    virtual const IntArray *giveKnotMultiplicity(int dim) { return & this->knotMultiplicity [ dim - 1 ]; }
    virtual const FloatArray *giveKnotValues(int dim) { return & this->knotValues [ dim - 1 ]; }
    virtual void evalN(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual double evaldNdx(FloatMatrix &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual void local2global(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual int global2local(FloatArray &answer, const FloatArray &lcoords, const FEICellGeometry &cellgeo) {
        OOFEM_ERROR("Not yet implemented.");
        return 0;
    }
    virtual void giveJacobianMatrixAt(FloatMatrix &jacobianMatrix, const FloatArray &lcoords, const FEICellGeometry &cellgeo);
    virtual int giveKnotSpanBasisFuncMask(const IntArray &knotSpan, IntArray &mask);
    virtual int giveNumberOfKnotSpanBasisFunctions(const IntArray &knotSpan);

    virtual const char *giveClassName() const { return "BSplineInterpolation"; }
    virtual bool hasSubPatchFormulation() { return true; }

    virtual IntegrationRule *giveIntegrationRule(int order)
    { OOFEM_ERROR("Not supported.");
      return NULL; }
    virtual IntegrationRule *giveBoundaryIntegrationRule(int order, int boundary)
    { OOFEM_ERROR("Not supported.");
      return NULL; }
    virtual IntegrationRule *giveBoundaryEdgeIntegrationRule(int order, int boundary)
    { OOFEM_ERROR("Not supported.");
      return NULL; }

protected:
    /**
     * Evaluates the nonvanishing basis functions of 1d BSpline (algorithm A2.2 from NURBS book)
     * @param span Knot span index (zero based).
     * @param u Value at which to evaluate.
     * @param p Degree.
     * @param U Knot vector.
     * @param N Computed p+1 nonvanishing functions (N_{span-p,p}-N_{span,p})
     * @warning Parameter u and span must be in a valid range.
     */
    void basisFuns(FloatArray &N, int span, double u, int p, const double *U);
    /**
     * Computes nonzero basis functions and their derivatives at u.
     *
     * For information on the algorithm, see A2.3 on p72 of
     * the NURBS book. The result is stored in the ders matrix, where
     * ders is of size (n+1,p+1) and the derivative
     * @f{align*}{
     * N(u)   &= \mathit{ders}(0,span-p+j) \quad\text{where } j=0...p \ \
     * N'(u)  &= \mathit{ders}(1,span-p+j) \quad\text{where } j=0...p \ \
     * N''(u) &= \mathit{ders}(2,span-p+j) \quad\text{where } j=0...p
     * @f}
     *
     * @param n Degree of the derivation.
     * @param u Parametric value.
     * @param span Knot span index (zero based).
     * @param p Degree.
     * @param U Knot vector.
     * @param ders Matrix containing the derivatives of the basis functions.
     *
     * @warning Parameters n, u and span must be in a valid range.
     */
    void dersBasisFuns(int n, double u, int span, int p, double *const U, FloatMatrix &ders);
    /**
     * Determines the knot span index (Algorithm A2.1 from the NURBS book)
     *
     * Determines the knot span for which there exists non-zero basis
     * functions. The span is the index k for which the parameter
     * u is valid in the (u_k,u_{k+1}] range.
     *
     * @param n Number of control points - 1.
     * @param u Parametric value.
     * @param p Degree.
     * @param U Knot vector.
     * @return Span index at u (zero based).
     * @warning Parameter u must be in a valid range.
     */
    int findSpan(int n, int p, double u, const double *U) const;
    /**
     * Returns the range of nonzero basis functions for given knot span and given degree.
     */
    void giveNonzeroBasisFuncInterval(int span, int deg, int &s, int &e) {
        s = span - deg;
        e = span;
    }
};
} // end namespace oofem
#endif // feibspline_h
