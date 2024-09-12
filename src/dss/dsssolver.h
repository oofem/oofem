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

#ifndef dsssolver_h
#define dsssolver_h

#include "oofemenv.h"
#include "sparselinsystemnm.h"
#include "convergedreason.h"
#include "sparsemtrx.h"

#define _IFT_DSSSolver_Name "dss"

namespace oofem {
class Domain;
class EngngModel;
class FloatMatrix;
class FloatArray;

/**
 * Implements the solution of linear system of equation in the form Ax=b using direct factorization method.
 * Can work with any sparse matrix implementation. However, the sparse matrix implementation have to support
 * its factorization (canBeFactorized method).
 */
class OOFEM_EXPORT DSSSolver : public SparseLinearSystemNM
{
public:
    /**
     * Constructor.
     * Creates new instance of DSS solver belonging to domain d and Engngmodel m.
     * @param d Domain which solver belongs to.
     * @param m Engineering model which solver belongs to.
     */
    DSSSolver(Domain *d, EngngModel *m);
    /// Destructor.
    virtual ~DSSSolver();

    ConvergedReason solve(SparseMtrx &A, FloatArray &b, FloatArray &x) override;

    const char *giveClassName() const override { return "DSSSolver"; }
    LinSystSolverType giveLinSystSolverType() const override { return ST_DSS; }
    SparseMtrxType giveRecommendedMatrix(bool symmetric) const override { return symmetric ? SMT_DSS_sym_LDL : SMT_DSS_unsym_LU; } ///@todo Check
};
} // end namespace oofem

#endif // dsssolver_h

