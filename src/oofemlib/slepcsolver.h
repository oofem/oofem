/* $Header: /home/cvs/bp/oofem/oofemlib/src/Attic/petscsolver.h,v 1.1.2.1 2004/04/05 15:19:43 bp Exp $ */
/*
 *
 *****    *****   ******  ******  ***   ***
 **   **  **   **  **      **      ** *** **
 **   **  **   **  ****    ****    **  *  **
 **   **  **   **  **      **      **     **
 **   **  **   **  **      **      **     **
 *****    *****   **      ******  **     **
 *****
 *****
 *****         OOFEM : Object Oriented Finite Element Code
 *****
 *****           Copyright (C) 1993 - 2003   Borek Patzak
 *****
 *****
 *****
 *****   Czech Technical University, Faculty of Civil Engineering,
 *****Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *****
 *****This program is free software; you can redistribute it and/or modify
 *****it under the terms of the GNU General Public License as published by
 *****the Free Software Foundation; either version 2 of the License, or
 *****(at your option) any later version.
 *****
 *****This program is distributed in the hope that it will be useful,
 *****but WITHOUT ANY WARRANTY; without even the implied warranty of
 *****MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *****GNU General Public License for more details.
 *****
 *****You should have received a copy of the GNU General Public License
 *****along with this program; if not, write to the Free Software
 *****Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
#ifndef slepcsolver_h
#define slepcsolver_h

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

#include "sparsegeneigenvalsystemnm.h"
#include "sparsemtrx.h"
#include "flotarry.h"

#ifdef __PETSC_MODULE
 #ifndef __MAKEDEPEND
  #include "petscsparsemtrx.h"
 #endif
#endif

#ifdef __SLEPC_MODULE
 #ifndef __MAKEDEPEND
  #include "slepceps.h"
 #endif
#endif

namespace oofem {
class Domain;
class EngngModel;
class FloatMatrix;

class SLEPcSolver : public SparseGeneralEigenValueSystemNM
{
private:
#ifdef __SLEPC_MODULE
    PetscSparseMtrx *    A;
    PetscSparseMtrx *B;
    /// Eigenvalue solver context
    EPS eps;
    /// flag if context initialized
    bool epsInit;
#endif

public:

    SLEPcSolver(int i, Domain *d, EngngModel *m); //Constructor

    ~SLEPcSolver(); ///Destructor

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


    IRResultType initializeFrom(InputRecord *ir);

    //Identification
    const char *giveClassName() const { return "SLEPcSolver"; }
    classType giveClassID() const { return SlepcSolverClass; }
};
} // end namespace oofem
#endif // slepcsolver_h
