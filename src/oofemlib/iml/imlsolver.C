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

#include <iml/cg.h>
#include <iml/gmres.h>

#include "imlsolver.h"
#include "sparsemtrx.h"
#include "floatarray.h"
#include "diagpre.h"
#include "voidprecond.h"
#include "compcol.h"
#include "iluprecond.h"
#include "icprecond.h"
#include "verbose.h"
#include "ilucomprowprecond.h"
#include "linsystsolvertype.h"
#include "classfactory.h"

#ifdef TIME_REPORT
 #include "timer.h"
#endif

namespace oofem {

REGISTER_SparseLinSolver(IMLSolver, ST_IML)

IMLSolver :: IMLSolver(Domain *d, EngngModel *m) : SparseLinearSystemNM(d, m)
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
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_IMLSolver_lstype);
    solverType = ( IMLSolverType ) val;

    tol = 1.e-5;
    IR_GIVE_OPTIONAL_FIELD(ir, tol, _IFT_IMLSolver_lstol);
    maxite = 200;
    IR_GIVE_OPTIONAL_FIELD(ir, maxite, _IFT_IMLSolver_lsiter);
    val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_IMLSolver_lsprecond);
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
        OOFEM_ERROR("IMLSolver::setSparseMtrxAsComponent: unknown preconditioner type");
    }

    // initialize precond attributes
    M->initializeFrom(ir);

    /*
     * IR_GIVE_OPTIONAL_FIELD (ir, precondAttributesRecord, "precondattributes");
     */
    return IRRT_OK;
}


NM_Status
IMLSolver :: solve(SparseMtrx *A, FloatArray *b, FloatArray *x)
{
    int result;

    // first check whether Lhs is defined
    if ( !A ) {
        OOFEM_ERROR("IMLSolver :: solve: unknown Lhs");
    }

    // and whether Rhs
    if ( !b ) {
        OOFEM_ERROR("IMLSolver :: solve: unknown Rhs");
    }

    // and whether previous Solution exist
    if ( !x ) {
        OOFEM_ERROR("IMLSolver :: solve: unknown solution array");
    }

    if ( x->giveSize() != b->giveSize() ) {
        OOFEM_ERROR("IMLSolver :: solve: size mismatch");
    }


    // check preconditioner
    if ( M ) {
        if ( ( precondInit ) || ( Lhs != A ) || ( this->lhsVersion != A->giveVersion() ) ) {
            M->init(* A);
        }
    } else {
        OOFEM_ERROR("IMLSolver :: solve: preconditioner creation error");
    }

    Lhs = A;
    this->lhsVersion = A->giveVersion();

 #ifdef TIME_REPORT
    Timer timer;
    timer.startTimer();
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
        OOFEM_ERROR("IMLSolver :: solve: unknown lsover type");
    }

 #ifdef TIME_REPORT
    timer.stopTimer();
    OOFEM_LOG_INFO( "IMLSolver info: user time consumed by solution: %.2fs\n", timer.getUtime() );
 #endif


    //solved = 1;
    return NM_Success;
}
} // end namespace oofem
