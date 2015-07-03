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


//   ******************************************
//   *** CLASS LDL^T Factorization          ***
//   ******************************************


#ifndef ldltfact_h
#define ldltfact_h

#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "floatarray.h"

namespace oofem {
class Domain;
class EngngModel;
class FloatMatrix;

/**
 * Implements the solution of linear system of equation in the form Ax=b using direct factorization method.
 * Can work with any sparse matrix implementation. However, the sparse matrix implementation have to support
 * its factorization (canBeFactorized method).
 */
class OOFEM_EXPORT LDLTFactorization : public SparseLinearSystemNM
{
private:

public:
    /// Constructor - creates new instance of LDLTFactorization, with number i, belonging to domain d and Engngmodel m.
    LDLTFactorization(Domain * d, EngngModel * m);
    /// Destructor
    virtual ~LDLTFactorization();

    /**
     * Solves the given linear system by LDL^T factorization.
     * Implementation rely on factorization support provided by mapped sparse matrix.
     * It calls Lhs->factorized()->backSubstitutionWith(*solutionArray). Sets solved flag to 1 if o.k.
     * @param A coefficient matrix
     * @param b right hand side
     * @param x solution array
     * @return NM_Status value
     */
    virtual NM_Status solve(SparseMtrx &A, FloatArray &b, FloatArray &x);

    virtual const char *giveClassName() const { return "LDLTFactorization"; }
    virtual LinSystSolverType giveLinSystSolverType() const { return ST_Direct; }
    virtual SparseMtrxType giveRecommendedMatrix(bool symmetric) const { return symmetric ? SMT_Skyline : SMT_SkylineU; }
};
} // end namespace oofem
#endif // ldltfact_h
