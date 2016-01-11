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

//   **************************************************
//   *** CLASS GENERALIZED INVERSE ITERATION SOLVER ***
//   **************************************************


#ifndef inverseit_h
#define inverseit_h

#include "sparsegeneigenvalsystemnm.h"
#include "floatarray.h"

namespace oofem {
class Domain;
class EngngModel;

/**
 * This class implements the class NumericalMethod instance Generalized Jacobi
 * Eigen Value Problem Solver
 *
 * DESCRIPTION :
 * Perform solution of eigen value problem in the form
 * K y = (omega)^2 M y
 *
 * TASKS :
 *
 * - solving problem
 *   solveYourselfAt.
 * - returning results (eigen values and associated eigen vectors).
 *
 * Variable description  :
 *
 *       A(N,N)    = STIFFNESS MATRIX (ASSUMED POZITIVE DEFINITE)        *
 *       B(N,N)    = MASS MATRIX (ASSUMED POZITIVE DEFINITE)             *
 *       X(N,N)    = MATRIX STORING EIGENVECTORS ON SOLUTION EXIT        *
 *       EIGV(N)   = VECTOR STORING EIGENVALUES ON SOLUTION EXIT         *
 *       D(N)      = WORKING VECTOR                                      *
 *       N         = ORDER OF WORKING AREA MATRICES A AND B              *
 *       RTOL      = CONVERGENCE TOLERANCE (USUALLY SET TO 10.**-12)     *
 *       NSMAX     = MAXIMUM NUMBER OF SWEEPS ALLOVED                    *
 *                                 (USUALLY SET TO 15)                   *
 *
 * OUTPUT : (after call solveYourselfAt)
 *       A(N,N)    = DIAGONALIZED STIFFNESS MATRIX                       *
 *       B(N,N)    = DIAGONALIZED MASS MATRIX                            *
 *       X(N,N)    = EIGENVECTORS STORED COLUMNWISE                      *
 *       EIGV(N)   = EIGENVALUES                                         *
 *
 *
 */
class OOFEM_EXPORT InverseIteration : public SparseGeneralEigenValueSystemNM
{
private:
    int n, nc, nitem;
    double rtol;
    int solved;

public:
    InverseIteration(Domain * d, EngngModel * m);
    virtual ~InverseIteration();

    virtual NM_Status solve(SparseMtrx &A, SparseMtrx &B, FloatArray &x, FloatMatrix &v, double rtol, int nroot);
    virtual const char *giveClassName() const { return "InverseIteration"; }
};
} // end namespace oofem
#endif // inverseit_h
