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

#ifndef sparselinsystemnm_h
#define sparselinsystemnm_h

#include "nummet.h"
#include "classtype.h"
#include "nmstatus.h"
#include "linsystsolvertype.h"

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
class SparseLinearSystemNM : public NumericalMethod
{
public:
    /// Constructor.
    SparseLinearSystemNM(int i, Domain *d, EngngModel *m);
    /// Destructor.
    virtual ~SparseLinearSystemNM();

    virtual const char *giveClassName() const { return "SparseLinearSystemNM"; }
    virtual classType giveClassID() const { return SparseLinearSystemNMClass; }

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
    virtual NM_Status solve(SparseMtrx *A, FloatArray *b, FloatArray *x) = 0;

    /**
     * Solves the given sparse linear system of equations @f$ A\cdot X=B @f$.
     * Default implementation calls solve multiple times.
     * @param A Coefficient matrix.
     * @param B Right hand side.
     * @param X Solution matrix.
     * @return Status of the solver.
     */
    virtual NM_Status solve(SparseMtrx *A, FloatMatrix &B, FloatMatrix &X);
};
} // end namespace oofem
#endif // sparselinsystemnm_h
