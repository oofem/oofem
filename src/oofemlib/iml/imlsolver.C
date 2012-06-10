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
#ifndef __IML_MODULE
 #include "imlsolver.h"

namespace oofem {
IMLSolver :: IMLSolver(int i, Domain *d, EngngModel *m) : SparseLinearSystemNM(i, d, m)
{
    _error("IMLSolver: can't create, IML support not compiled");
}

IMLSolver :: ~IMLSolver() { }

IRResultType
IMLSolver :: initializeFrom(InputRecord *ir) { return IRRT_OK; }

NM_Status
IMLSolver :: solve(SparseMtrx *A, FloatArray *b, FloatArray *x) { return NM_NoSuccess; }
} // end namespace oofem
#endif

#ifdef __IML_MODULE

 #include <iml/cg.h>
 #include <iml/gmres.h>
 #include "imlsolver.h"
 #include "sparsemtrx.h"
 #include "flotarry.h"
 #include "diagpre.h"
 #include "voidprecond.h"
 #include "compcol.h"
 #include "iluprecond.h"
 #include "icprecond.h"
 #include "verbose.h"
 #include "ilucomprowprecond.h"

 #ifdef TIME_REPORT
  #include "clock.h"
 #endif

namespace oofem {
IMLSolver :: IMLSolver(int i, Domain *d, EngngModel *m) : SparseLinearSystemNM(i, d, m)
{
    Lhs = NULL;
    M = NULL;
    solverType = IML_ST_CG;
    precondType = IML_VoidPrec;
    precondInit = true;
}


IMLSolver :: ~IMLSolver() {
    if ( M ) {
        delete M;
    }
}

IRResultType
IMLSolver :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    int val;

    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_IMLSolver_lstype, "stype"); // Macro
    solverType = ( IMLSolverType ) val;

    tol = 1.e-5;
    IR_GIVE_OPTIONAL_FIELD(ir, tol, IFT_IMLSolver_lstol, "lstol"); // Macro
    maxite = 200;
    IR_GIVE_OPTIONAL_FIELD(ir, maxite, IFT_IMLSolver_lsiter, "lsiter"); // Macro
    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, IFT_IMLSolver_lsprecond, "lsprecond"); // Macro
    precondType = ( IMLPrecondType ) val;

    // create preconditioner
    if ( precondType == IML_DiagPrec ) {
        M = new DiagPreconditioner();
    } else if ( precondType == IML_VoidPrec ) {
        M = new VoidPreconditioner();
    } else if ( precondType == IML_ILU_CompRowPrec ) {
        M = new CompRow_ILUPreconditioner();
    } else if ( precondType == IML_ILU_CompColPrec ) {
        M = new CompCol_ILUPreconditioner();
    } else if ( precondType == IML_ICPrec ) {
        M = new CompCol_ICPreconditioner();
    } else {
        _error("setSparseMtrxAsComponent: unknown preconditioner type");
    }

    // initialize precond attributes
    M->initializeFrom(ir);

    /*
     * IR_GIVE_OPTIONAL_FIELD (ir, precondAttributesRecord, "precondattributes"); // Macro
     */
    return IRRT_OK;
}


NM_Status
IMLSolver :: solve(SparseMtrx *A, FloatArray *b, FloatArray *x)
{
    int result;

    // first check whether Lhs is defined
    if ( !A ) {
        _error("solveYourselfAt: unknown Lhs");
    }

    // and whether Rhs
    if ( !b ) {
        _error("solveYourselfAt: unknown Rhs");
    }

    // and whether previous Solution exist
    if ( !x ) {
        _error("solveYourselfAt: unknown solution array");
    }

    if ( x->giveSize() != b->giveSize() ) {
        _error("solveYourselfAt: size mismatch");
    }


    // check preconditioner
    if ( M ) {
        if ( ( precondInit ) || ( Lhs != A ) || ( this->lhsVersion != A->giveVersion() ) ) {
            M->init(* A);
        }
    } else {
        _error("setSparseMtrxAsComponent: preconditioner creation error");
    }

    Lhs = A;
    this->lhsVersion = A->giveVersion();

 #ifdef TIME_REPORT
    //clock_t tstart = clock();
    oofem_timeval tstart;
    getUtime(tstart);
 #endif


    if ( solverType == IML_ST_CG ) {
        int mi = this->maxite;
        double t = this->tol;
        result = CG(* Lhs, * x, * b, * M, mi, t);
        OOFEM_LOG_INFO("CG(%s): flag=%d, nite %d, achieved tol. %g\n", M->giveClassName(), result, mi, t);
    } else if ( solverType == IML_ST_GMRES ) {
        int mi = this->maxite, restart = 100;
        double t = this->tol;
        FloatMatrix H(restart + 1, restart); // storage for upper Hesenberg
        result = GMRES(* Lhs, * x, * b, * M, H, restart, mi, t);
        OOFEM_LOG_INFO("GMRES(%s): flag=%d, nite %d, achieved tol. %g\n", M->giveClassName(), result, mi, t);
    } else {
        _error("solveYourselfAt: unknown lsover type");
    }

 #ifdef TIME_REPORT
    oofem_timeval ut;
    getRelativeUtime(ut, tstart);
    OOFEM_LOG_INFO( "IMLSolver info: user time consumed by solution: %.2fs\n", ( double ) ( ut.tv_sec + ut.tv_usec / ( double ) OOFEM_USEC_LIM ) );
 #endif


    //solved = 1;
    return NM_Success;
}
} // end namespace oofem
#endif //ifdef __IML_MODULE
