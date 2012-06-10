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

#include "nrsolver2.h"
#include "verbose.h"
#include "timestep.h"
#include "mathfem.h"
#include "usrdefsub.h"

namespace oofem {
#define nrsolver_SMALL_NUM 1.e-20
#define NRSOLVER_MAX_REL_ERROR_BOUND 1.e10
#define NRSOLVER_MAX_RESTARTS 4
#define NRSOLVER_RESET_STEP_REDUCE 0.25
#define NRSOLVER_DEFAULT_NRM_TICKS 10


NRSolver2 :: NRSolver2(int i, Domain *d, EngngModel *m, EquationID ut) :
    SparseNonLinearSystemNM(i, d, m, ut)
{
    //
    // constructor
    //
    nsmax  = 60;     // default maximum number of sweeps allowed
    rtol   = 10.E-3; // convergence tolerance
    //Psi    = 0.1;       // displacement control on
    solved = 0;
    NR_Mode = NR_OldMode = nrsolverModifiedNRM;
    NR_ModeTick = -1; // do not switch to calm_NR_OldMode
    MANRMSteps = 0;

    linSolver = NULL;
    linesearchSolver = NULL;
    lsFlag = 0; // no line-search
}

NRSolver2 :: ~NRSolver2()
{
    //
    // destructor
    //
    if ( linSolver ) {
        delete linSolver;
    }

    if ( linesearchSolver ) {
        delete linesearchSolver;
    }
}


NM_Status
NRSolver2 :: solve(SparseMtrx *k, FloatArray *R, FloatArray *R0,
                   FloatArray *X, FloatArray *dX, FloatArray *F,
                   double &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                   int &nite, TimeStep *tNow)
//
// this function solve the problem of the unbalanced equilibrium
// using NR scheme
//
//
{
    FloatArray rhs, ddX, RT, XInitial;
    double RRT, forceErr, dispErr = 0.;
    double drr;
    int neq = R->giveSize();
    int irest = 0;
    NM_Status status;

    OOFEM_LOG_INFO("Time       Iteration       ForceError      DisplError\n__________________________________________________________\n");

    l = 1.0;
    XInitial = * X;

    status = NM_None;
    this->giveLinearSolver();

    // compute total load R = R+R0
    RT = * R;
    if ( R0 ) {
        RT.add(*R0);
    }

restart:
    dX->zero();
    engngModel->updateComponent(tNow, NonLinearLhs, domain);

    deltaL = tNow->giveTimeIncrement();

    ddX.resize(neq);
    engngModel->updateComponent(tNow, InternalRhs, domain);
    rhs.beDifferenceOf(RT, *F);

    RRT = RT.computeSquaredNorm();

    nite = 0;

    do {
        nite++;

        if ( nite > 1 ) {
            if ( ( NR_Mode == nrsolverFullNRM ) || ( ( NR_Mode == nrsolverAccelNRM ) && ( nite % MANRMSteps == 0 ) ) ) {
                engngModel->updateComponent(tNow, NonLinearLhs, domain);
            }
        }

        linSolver->solve(k, & rhs, & ddX);
        //
        // update solution
        //
        if ( this->lsFlag && ( nite != 1 ) ) {
            // line search
            LineSearchNM :: LS_status status;
            IntArray prescribedEqs(0);
            double eta;

            this->giveLineSearchSolver()->solve(X, & ddX, F, R, R0, prescribedEqs, 1.0, eta, status, tNow);
            dX->add(ddX);
        } else {
            X->add(ddX);
            dX->add(ddX);
            tNow->incrementStateCounter();     // update solution state counter
            //
            // Convergence check
            //
            engngModel->updateComponent(tNow, InternalRhs, domain);
        }

        rhs.beDifferenceOf(RT, *F);

        //
        // compute forceError
        //
        // err is relative error of unbalanced forces
        forceErr = rhs.computeSquaredNorm();
        // we compute a relative error norm
        if ( ( RRT ) > nrsolver_SMALL_NUM ) {
            forceErr = sqrt( forceErr / ( RRT ) );
        } else {
            forceErr = sqrt(forceErr); // absolute norm
        }

        //
        // compute displacement error
        //
        // err is relative displacement change
        drr = X->computeSquaredNorm();
        if ( drr < nrsolver_SMALL_NUM ) {
            dispErr = 1.;
        } else {
            dispErr = sqrt( ddX.computeSquaredNorm() / drr );
        }

        //
        // Restart if nite >= nsmax of if force or displacement error is bigger
        // than allowed limit (rtol * CALM_MAX_REL_ERROR_BOUND)
        //
        if ( ( nite >= nsmax ) ||
            ( fabs(forceErr) > rtol * NRSOLVER_MAX_REL_ERROR_BOUND ) ||
            ( fabs(dispErr)  > rtol * NRSOLVER_MAX_REL_ERROR_BOUND ) ) {
            irest++;
            if ( irest <= NRSOLVER_MAX_RESTARTS ) {
                // convergence problems
                // there must be step restart followed by decrease of step length
                // status |= NM_ForceRestart;
                // reduce step length

                /*
                 * double time;
                 * time = tNow->giveTime() - tNow->giveTimeIncrement()*(1.0-NRSOLVER_RESET_STEP_REDUCE) ;
                 * deltaL =  deltaL * NRSOLVER_RESET_STEP_REDUCE ;
                 * if (deltaL < minStepLength)  deltaL = minStepLength;
                 *
                 * tNow -> setTime(time);
                 * tNow -> setTimeIncrement(tNow->giveTimeIncrement()*NRSOLVER_RESET_STEP_REDUCE);
                 * tNow->incrementStateCounter();              // update solution state counter
                 */

                // restore previous total displacement vector
                *X = XInitial;
                // reset all changes fro previous equilibrium state
                engngModel->initStepIncrements();
                dX->zero();
#ifdef VERBOSE
                OOFEM_LOG_INFO("NRSolver2 iteration Reset ...\n");
#endif
                NR_OldMode  = NR_Mode;
                NR_Mode     = nrsolverFullNRM;
                NR_ModeTick = NRSOLVER_DEFAULT_NRM_TICKS;
                goto restart;
            } else {
                status = NM_NoSuccess;
                _warning2("NRSolver2 - convergence not reached after %d iterations", nsmax);
                break;
            }
        }

        OOFEM_LOG_INFO("%-10d %-15d %-15e %-15e\n", ( int ) tNow->giveTargetTime(), nite, forceErr, dispErr);
    } while ( ( fabs(forceErr) > rtol ) || ( fabs(dispErr) > rtol ) || ( nite < minIterations ) );

    //
    // end of iteration
    //
    status |= NM_Success;
    solved = 1;
    return status;
}

IRResultType
NRSolver2 :: initializeFrom(InputRecord *ir)
//
//
//
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    nsmax = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, nsmax, IFT_NRSolver_maxiter, "maxiter"); // Macro
    if ( nsmax < 30 ) {
        nsmax = 30;
    }

    minIterations = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, minIterations, IFT_NRSolver_miniterations, "miniter"); // Macro

    minStepLength = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, minStepLength, IFT_NRSolver_minsteplength, "minsteplength"); // Macro

    // read if MANRM method is used
    MANRMSteps = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, MANRMSteps, IFT_NRSolver_manrmsteps, "manrmsteps"); // Macro
    if ( MANRMSteps > 0 ) {
        NR_Mode = NR_OldMode = nrsolverAccelNRM;
    } else {
        NR_Mode = nrsolverModifiedNRM;
    }

    int _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, IFT_NRSolver_lstype, "lstype"); // Macro
    solverType = ( LinSystSolverType ) _val;
    this->giveLinearSolver()->initializeFrom(ir);

    // read relative error tolerances of the solver fo each cc
    IR_GIVE_FIELD(ir, rtol, IFT_NRSolver_rtolv, "rtolv"); // Macro

    this->lsFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, lsFlag, IFT_NRSolver_linesearch, "linesearch"); // Macro

    if ( this->lsFlag ) {
        this->giveLineSearchSolver()->initializeFrom(ir);
    }

    return IRRT_OK;
}


contextIOResultType
NRSolver2 :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    return CIO_OK;
}


contextIOResultType
NRSolver2 :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    return CIO_OK;
}


SparseLinearSystemNM *
NRSolver2 :: giveLinearSolver()
{
    if ( linSolver ) {
        if ( linSolver->giveLinSystSolverType() == solverType ) {
            return linSolver;
        } else {
            delete linSolver;
        }
    }

    linSolver = CreateUsrDefSparseLinSolver(solverType, 1, domain, engngModel);
    if ( linSolver == NULL ) {
        _error("giveLinearSolver: linear solver creation failed");
    }

    return linSolver;
}

LineSearchNM *
NRSolver2 :: giveLineSearchSolver()
{
    if ( linesearchSolver == NULL ) {
        linesearchSolver = new LineSearchNM(1, this->giveDomain(), engngModel);
    }

    return linesearchSolver;
}
} // end namespace oofem
