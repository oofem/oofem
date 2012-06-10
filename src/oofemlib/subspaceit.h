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

//   ****************************************
//   *** CLASS SUBSPACE ITERATION  SOLVER ***
//   ****************************************


#ifndef subspaceit_h
#define subspaceit_h

#include "sparsegeneigenvalsystemnm.h"
#include "sparsemtrx.h"
#include "flotarry.h"

namespace oofem {
class Domain;
class EngngModel;


class SubspaceIteration : public SparseGeneralEigenValueSystemNM
{
    /*
     * This class implements the class NumericalMethod instance Subspace Iteration
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
     * C     - - - INPUT  DATA - - -
     * C
     * C     A(NWK)  -  STIFFNESS MATRIX
     * C     B(NWM)  -  MASS MARTRIX
     * C     NN  -      SIZE OF PROBLEM
     * C     NNM  -  NN+1
     * C     NROOT  -  REQUIRED NUMBER OF
     * C     RTOL  -  KRITERIUM KONVERGENCE VLASTNICH CISEL
     * C     NC  -  POCET VEKTORU SIMULTANNI ITERACE, DOPORUCUJE SE VOLIT
     * C     NC = MIN (2*NROOT , NROOT+8 )
     * C     NITEM  -  MAXIMALNI POCET ITERACI (OBYC. 16)
     * C
     * C     - - - PRACOVNI POLE - - -
     * C
     * C     TT(NN),W(NN),D(NC),RTOLV(NC),BUP(NC),BLO(NC),BUPC(NC)
     * C     AR(NC,NC)  -  PRACOVNI MATICE - PROJEKCE MATICE  A
     * C     BR(NC,NC)  -  PROJEKCE MATICE  B
     * C
     * C     - - - VYSTUPNI DATA - - -
     * C
     * C     EIGV(NROOT)  -  VLASTNI CISLA
     * C     R(NN,NROOT)  -  VLASTNI VEKTORY
     */

private:
    //SparseMtrx*   a ;
    //SparseMtrx*   b ;
    FloatMatrix *ar;
    FloatMatrix *br;
    //FloatArray*    _eigv;  // only pointer to Engngmethod data, not ownership
    //FloatMatrix*   _r;    // only pointer to Engngmethod data, not ownership
    FloatMatrix *vec;
    int n, nc, nsmax, nitem;
    //double         rtol  ;
    int solved;


public:
    SubspaceIteration(int i, Domain *d, EngngModel *m);
    // constructor
    virtual ~SubspaceIteration();               // destructor


    /**
     * Solves the given sparse generalized eigen value system of equations Ax = o^2 Bx.
     * @param A coefficient matrix
     * @param B coefficient matrix
     * @param x eigen vector(s)
     * @param o eigen value(s)
     * @param rtol tolerance
     * @param nroot number of required eigenvalues
     * @return NM_Status value
     */
    virtual NM_Status solve(SparseMtrx *a, SparseMtrx *b, FloatArray *_eigv, FloatMatrix *r, double rtol, int nroot);


    virtual IRResultType initializeFrom(InputRecord *ir);

    // identification
    virtual const char *giveClassName() const { return "SubspaceIterationSolver"; }
    virtual classType giveClassID() const { return SubspaceIterationSolverClass; }
};
} // end namespace oofem
#endif // subspaceit_h
