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

#ifndef sparselinsystemnm_h
#define sparselinsystemnm_h

#include "nummet.h"
#include "nmstatus.h"
#include "linsystsolvertype.h"
#include "sparsemtrxtype.h"

namespace oofem {
class EngngModel;
class SparseMtrx;
class FloatArray;

/**
 * This base class is an abstraction for all numerical methods solving sparse
 * linear system of equations. The purpose of this class is to declare
 * the general interface to all numerical methods solving this kind of
 * problem. This interface allows to use any suitable
 * instance of the Numerical method class to the solve problem,
 * and leave the  whole engineering model code,
 * including mapping, unchanged, because all instances of this class
 * provide the common interface.
 */
class OOFEM_EXPORT SparseLinearSystemNM : public NumericalMethod
{
public:
    /// Constructor.
    SparseLinearSystemNM(Domain * d, EngngModel * m);
    /// Destructor.
    virtual ~SparseLinearSystemNM();

    virtual const char *giveClassName() const = 0;
    std :: string errorInfo(const char *func) { return std :: string(this->giveClassName()) + func; }

    /**
     * @return LinSystSolverType value, corresponding to receiver.
     */
    virtual LinSystSolverType giveLinSystSolverType() const = 0;

    /**
     * Solves the given sparse linear system of equations @f$ A\cdot x=b @f$.
     * @param A Coefficient matrix.
     * @param b Right hand side.
     * @param x Solution array.
     * @return Status of the solver.
     */
    virtual NM_Status solve(SparseMtrx &A, FloatArray &b, FloatArray &x) = 0;

    /**
     * Solves the given sparse linear system of equations @f$ A\cdot X=B @f$.
     * Default implementation calls solve multiple times.
     * @param A Coefficient matrix.
     * @param B Right hand side.
     * @param X Solution matrix.
     * @return Status of the solver.
     */
    virtual NM_Status solve(SparseMtrx &A, FloatMatrix &B, FloatMatrix &X);
    /**
     * Returns the recommended sparse matrix type for this solver.
     */
    virtual SparseMtrxType giveRecommendedMatrix(bool symmetric) const = 0;
};
} // end namespace oofem
#endif // sparselinsystemnm_h
