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

#ifndef mklpardisosolver_h
#define mklpardisosolver_h

#include "sparselinsystemnm.h"

namespace oofem {

/**
 * Implements the solution of linear system of equation in the form @f$ A\cdot x=b @f$ using solvers
 * from MKL Pardiso. It will not work with Pardiso 5.0.
 * Only works with the CSC (compressed sparse column) sparse matrix implementation.
 * 
 * @author Mikael Öhman
 */
class OOFEM_EXPORT MKLPardisoSolver : public SparseLinearSystemNM
{
public:
    /**
     * Constructor.
     * @param d Domain which solver belongs to.
     * @param m Engineering model which solver belongs to.
     */
    MKLPardisoSolver(Domain * d, EngngModel * m);

    virtual ~MKLPardisoSolver();

    virtual NM_Status solve(SparseMtrx &A, FloatArray &b, FloatArray &x);

    virtual const char *giveClassName() const { return "MKLPardisoSolver"; }
    virtual LinSystSolverType giveLinSystSolverType() const { return ST_MKLPardiso; }
    virtual SparseMtrxType giveRecommendedMatrix(bool symmetric) const { return SMT_CompCol; }
};
} // end namespace oofem
#endif // mklpardisosolver_h
