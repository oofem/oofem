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

#ifndef dsssolver_h
#define dsssolver_h

#include "sparselinsystemnm.h"
#include "sparsemtrx.h"

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
class DSSSolver : public SparseLinearSystemNM
{
public:
    /**
     * Constructor.
     * Creates new instance of LDLTFactorization, with number i, belonging to domain d and Engngmodel m.
     * @param i Solver number.
     * @param d Domain which solver belongs to.
     * @param m Engineering model which solver belongs to.
     */
    DSSSolver(int i, Domain *d, EngngModel *m);
    /// Destructor.
    virtual ~DSSSolver();

    /**
     * Solves the given linear system by @f$L\cdot D\cdot L^{\mathrm{T}}@f$ factorization.
     * Implementation rely on factorization support provided by mapped sparse matrix.
     * It calls Lhs->factorized()->backSubstitutionWith(*solutionArray). Sets solved flag to 1 if o.k.
     * @param A Coefficient matrix.
     * @param b Right hand side.
     * @param x Solution array.
     * @return NM_Status value.
     */
    virtual NM_Status solve(SparseMtrx *A, FloatArray *b, FloatArray *x);
    virtual IRResultType initializeFrom(InputRecord *ir);

    // identification
    virtual const char *giveClassName() const { return "LDLTFactorization"; }
    virtual classType giveClassID() const { return LDLTFactorizationClass; }
    virtual LinSystSolverType giveLinSystSolverType() const { return ST_DSS; }
};
} // end namespace oofem
#endif // dsssolver_h









