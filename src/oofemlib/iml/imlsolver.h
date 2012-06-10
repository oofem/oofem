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

#ifndef imlsolver_h
#define imlsolver_h

#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "flotarry.h"

#include "precond.h"

namespace oofem {
class Domain;
class EngngModel;
class FloatMatrix;

/**
 * Implements the solution of linear system of equation in the form @f$ A\cdot x=b @f$ using iterative solvers
 * from IML++ library. Can work with any sparse matrix implementation.
 */
class IMLSolver : public SparseLinearSystemNM
{
private:
    /// Solver type.
    enum IMLSolverType { IML_ST_CG, IML_ST_GMRES };
    /// Preconditioner type.
    enum IMLPrecondType { IML_VoidPrec, IML_DiagPrec, IML_ILU_CompColPrec, IML_ILU_CompRowPrec, IML_ICPrec };

    /// Last mapped Lhs matrix
    SparseMtrx *Lhs;
    /// Last mapped matrix version
    SparseMtrx :: SparseMtrxVersionType lhsVersion;
    /// Preconditioner.
    Preconditioner *M;
    /// IML Solver type.
    IMLSolverType solverType;
    /// IML Preconditioner type.
    IMLPrecondType precondType;
    /// Precond. init flag.
    bool precondInit;
    // Preconditioner attribute string
    // InputRecord precondAttributes;

    /// Tolerance of residual.
    double tol;
    /// Max number of iterations.
    int maxite;


public:
    /// Constructor. Creates new instance of LDLTFactorization, with number i, belonging to domain d and Engngmodel m.
    IMLSolver(int i, Domain *d, EngngModel *m);
    /// Destructor
    virtual ~IMLSolver();

    /**
     * Solves the given linear system iteratively by method described by IMLSolverType.
     * @param A Coefficient matrix.
     * @param b Right hand side.
     * @param x Solution array.
     * @return Status value.
     */
    virtual NM_Status solve(SparseMtrx *A, FloatArray *b, FloatArray *x);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "IMLSolver"; }
    virtual classType giveClassID() const { return IMLSolverClass; }
    virtual LinSystSolverType giveLinSystSolverType() const { return ST_IML; }
};
} // end namespace oofem
#endif // imlsolver_h
