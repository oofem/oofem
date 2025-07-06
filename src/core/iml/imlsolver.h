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

#ifndef imlsolver_h
#define imlsolver_h

#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "convergedreason.h"
#include "floatarray.h"
#include "precond.h"

#include <memory>

///@name Input fields for IMLSolver
//@{
#define _IFT_IMLSolver_Name "iml"
#define _IFT_IMLSolver_stype "stype"
#define _IFT_IMLSolver_lstol "lstol"
#define _IFT_IMLSolver_lsiter "lsiter"
#define _IFT_IMLSolver_lsprecond "lsprecond"
//@}

namespace oofem {
class Domain;
class EngngModel;
class FloatMatrix;

/**
 * Implements the solution of linear system of equation in the form @f$ A\cdot x=b @f$ using iterative solvers
 * from IML++ library. Can work with any sparse matrix implementation.
 */
class OOFEM_EXPORT IMLSolver : public SparseLinearSystemNM
{
private:
    /// Solver type.
    enum IMLSolverType { IML_ST_CG, IML_ST_GMRES };
    /// Preconditioner type.
    enum IMLPrecondType { IML_VoidPrec, IML_DiagPrec, IML_ILU_CompColPrec, IML_ILU_CompRowPrec, IML_ICPrec };

    /// Last mapped Lhs matrix
    SparseMtrx *lhs;
    /// Last mapped matrix version
    SparseMtrx :: SparseMtrxVersionType lhsVersion;
    /// Preconditioner.
    std::unique_ptr<Preconditioner> M;
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
    IMLSolver(Domain * d, EngngModel * m);
    /// Destructor
    virtual ~IMLSolver() {}

    ConvergedReason solve(SparseMtrx &A, FloatArray &b, FloatArray &x) override;

    void initializeFrom(InputRecord &ir) override;
    const char *giveClassName() const override { return "IMLSolver"; }
    LinSystSolverType giveLinSystSolverType() const override { return ST_IML; }
    SparseMtrxType giveRecommendedMatrix(bool symmetric) const override { return symmetric ? SMT_SymCompCol : SMT_CompCol; }
};
} // end namespace oofem
#endif // imlsolver_h
