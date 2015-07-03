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


#ifndef linesearch_h
#define linesearch_h

#include "nummet.h"
#include "nmstatus.h"
#include "floatarray.h"

///@name Input fields for LineSearchNM
//@{
#define _IFT_LineSearchNM_lsearchtol "lsearchtol"
#define _IFT_LineSearchNM_lsearchamp "lsearchamp"
#define _IFT_LineSearchNM_lsearchmaxeta "lsearchmaxeta"
//@}

namespace oofem {
class EngngModel;
class SparseMtrx;
class FloatArray;
class TimeStep;

/**
 * This base class is an abstraction/implementation for numerical method solving
 * line search optimization problem.
 */
class OOFEM_EXPORT LineSearchNM : public NumericalMethod
{
public:
    enum LS_status { ls_ok, ls_failed };

    int max_iter;
    double ls_tolerance;
    double amplifFactor;
    double maxEta, minEta;
    FloatArray eta;
    FloatArray prod;

public:
    /// Constructor
    LineSearchNM(Domain * d, EngngModel * m);

    /**
     * Solves the line search optimization problem in the form of @f$ g(r)=0; r_{new}=r_{old}+\eta\delta r; 0 < \eta < 1 @f$,
     * The aim is to find @f$ \eta @f$ so that the @f$ g(r) @f$ has decreased sufficiently.
     * The total solution vector is updated at exit as well as InternalRhs vector.
     * @param r  Old total solution.
     * @param dr Increment of solution.
     * @param F  Old InternalRhs (real internal forces).
     * @param R  Reference incremental Rhs (incremental load).
     * @param R0 Initial Rhs (initial load).
     * @param eqnmask Equation numbers to mask out from computations.
     * @param lambda Scaling of R.
     * @param etaValue Reached eta value.
     * @param status Linesearch status
     * @param tStep Time step.
     */
    virtual NM_Status solve(FloatArray &r, FloatArray &dr, FloatArray &F, FloatArray &R, FloatArray *R0,
                            IntArray &eqnmask, double lambda, double &etaValue, LS_status &status, TimeStep *tStep);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "LineSearchNM"; }

protected:
    void search(int istep, FloatArray &prod, FloatArray &eta, double amp, double maxeta, double mineta, int &status);
};
} // end namespace oofem
#endif // linesearch_h
