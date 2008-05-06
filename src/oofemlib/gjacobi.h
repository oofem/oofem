/* $Header: /home/cvs/bp/oofem/oofemlib/src/gjacobi.h,v 1.4 2003/04/06 14:08:24 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   ***************************************************
//   *** CLASS GENERALIZED JACOBI EIGEN VALUE SOLVER ***
//   ***************************************************


#ifndef gjacobi_h
#define gjacobi_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "nummet.h"
#include "skyline.h"
#include "flotarry.h"

#include "nmstatus.h"

class Domain;
class EngngModel;


class GJacobi : public NumericalMethod
{
    /*
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
private:
    FloatMatrix *a;
    FloatMatrix *b;
    FloatArray *eigv;    // only pointer to caller data, not ownership
    FloatMatrix *x;      // only pointer to caller data, not ownership
    int n, nsmax;
    double rtol;
    int solved;


public:
    GJacobi(int i, Domain *d, EngngModel *m);
    // constructor
    ~GJacobi();                         // destructor

    // solving
    void               solveYourselfAt(TimeStep *);
    void               updateYourself();
    void               updateYourselfExceptLhs();


    /**
     * Solves the given sparse generalized eigen value system of equations Ax = o^2 Bx.
     * @param A coefficient matrix
     * @param B coefficient matrix
     * @param x eigen vector(s)
     * @param o eigen value(s)
     * @return NM_Status value
     */
    virtual NM_Status solve(FloatMatrix *a, FloatMatrix *b, FloatArray *eigv, FloatMatrix *x);

    IRResultType initializeFrom(InputRecord *ir);

    // identification
    const char *giveClassName() const { return "GeneralizedJacobiSolver"; }
    classType giveClassID() const { return GeneralizedJacobiSolverClass; }
protected:
};

#endif // gjacobi_h









