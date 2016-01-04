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

#include "nrsolver.h"
#include "verbose.h"
#include "timestep.h"
#include "mathfem.h"
// includes for ddc - not very clean (NumMethod knows what is "node" and "dof")
#include "node.h"
#include "element.h"
#include "generalboundarycondition.h"
#include "dof.h"
#include "function.h"
#include "linesearch.h"
#include "classfactory.h"
#include "exportmodulemanager.h"
#include "engngm.h"
#include "parallelcontext.h"
#include "unknownnumberingscheme.h"

#ifdef __PETSC_MODULE
 #include "petscsolver.h"
 #include "petscsparsemtrx.h"
#endif

#include <cstdio>

namespace oofem {
#define nrsolver_ERROR_NORM_SMALL_NUM 1.e-6
#define NRSOLVER_MAX_REL_ERROR_BOUND 1.e20
#define NRSOLVER_MAX_RESTARTS 4
#define NRSOLVER_RESET_STEP_REDUCE 0.25
#define NRSOLVER_DEFAULT_NRM_TICKS 10

REGISTER_SparseNonLinearSystemNM(NRSolver)

NRSolver :: NRSolver(Domain *d, EngngModel *m) :
    SparseNonLinearSystemNM(d, m), prescribedDofs(), prescribedDofsValues()
{
    //
    // constructor
    //
    nsmax = 60;     // default maximum number of sweeps allowed
    deltaL = 1.0;
    solved = 0;
    NR_Mode = NR_OldMode = nrsolverModifiedNRM;
    NR_ModeTick = -1; // do not switch to calm_NR_OldMode
    MANRMSteps = 0;
    numberOfPrescribedDofs = 0;
    prescribedDofsFlag = false;
    prescribedEqsInitFlag = false;
    prescribedDisplacementTF = 0;
    lsFlag = 0; // no line-search
    
    constrainedNRFlag = false; 
    this->forceErrVecOld.resize(0);
    this->forceErrVec.resize(0);
    constrainedNRalpha = 0.5; // default

    smConstraintVersion = 0;
    mCalcStiffBeforeRes = true;
}


NRSolver :: ~NRSolver()
{
}


IRResultType
NRSolver :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // Choosing a big "enough" number. (Alternative: Force input of maxinter)
    nsmax = ( int ) 1e8;
    IR_GIVE_OPTIONAL_FIELD(ir, nsmax, _IFT_NRSolver_maxiter);
    if ( nsmax < 0 ) {
        OOFEM_ERROR("nsmax < 0");
    }

    minIterations = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, minIterations, _IFT_NRSolver_miniterations);

    minStepLength = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, minStepLength, _IFT_NRSolver_minsteplength);

    // read if MANRM method is used
    MANRMSteps = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, MANRMSteps, _IFT_NRSolver_manrmsteps);
    if ( MANRMSteps > 0 ) {
        NR_Mode = NR_OldMode = nrsolverAccelNRM;
    } else {
        NR_Mode = nrsolverModifiedNRM;
    }

    int _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, _IFT_NRSolver_lstype);
    solverType = ( LinSystSolverType ) _val;
    this->giveLinearSolver()->initializeFrom(ir); // make sure that linear solver is initialized from ir as well

    // read relative error tolerances of the solver
    // if rtolv provided set to this tolerance both rtolf and rtold
    rtolf.resize(1);
    rtolf.at(1) = 1.e-3; // Default value.
    rtold.resize(1);
    rtold = FloatArray{0.0}; // Default off (0.0 or negative values mean that residual is ignored)
    IR_GIVE_OPTIONAL_FIELD(ir, rtolf.at(1), _IFT_NRSolver_rtolv);
    IR_GIVE_OPTIONAL_FIELD(ir, rtold.at(1), _IFT_NRSolver_rtolv);

    // read optional force and displacement tolerances
    IR_GIVE_OPTIONAL_FIELD(ir, rtolf.at(1), _IFT_NRSolver_rtolf);
    IR_GIVE_OPTIONAL_FIELD(ir, rtold.at(1), _IFT_NRSolver_rtold);

    prescribedDofs.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, prescribedDofs, _IFT_NRSolver_ddm);
    prescribedDofsValues.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, prescribedDofsValues, _IFT_NRSolver_ddv);
    prescribedDisplacementTF = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, prescribedDisplacementTF, _IFT_NRSolver_ddfunc);

    numberOfPrescribedDofs = prescribedDofs.giveSize() / 2;
    if ( numberOfPrescribedDofs != prescribedDofsValues.giveSize() ) {
        OOFEM_ERROR("direct displacement mask size mismatch");
    }

    if ( numberOfPrescribedDofs ) {
        prescribedDofsFlag = true;
    } else {
        prescribedDofsFlag = false;
    }

    this->lsFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, lsFlag, _IFT_NRSolver_linesearch);

    if ( this->lsFlag ) {
        this->giveLineSearchSolver()->initializeFrom(ir);
    }

    int calcStiffBeforeResFlag = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, calcStiffBeforeResFlag, _IFT_NRSolver_calcstiffbeforeres);
    if ( calcStiffBeforeResFlag == 0 ) {
        mCalcStiffBeforeRes = false;
    }

 
    this->constrainedNRminiter = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->constrainedNRminiter, _IFT_NRSolver_constrainedNRminiter);
    this->constrainedNRFlag = this->constrainedNRminiter != 0;

    return SparseNonLinearSystemNM :: initializeFrom(ir);
}


NM_Status
NRSolver :: solve(SparseMtrx &k, FloatArray &R, FloatArray *R0, FloatArray *iR,
                  FloatArray &X, FloatArray &dX, FloatArray &F,
                  const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                  int &nite, TimeStep *tStep)
//
// this function solve the problem of the unbalanced equilibrium
// using NR scheme
//
//
{
    // residual, iteration increment of solution, total external force
    FloatArray rhs, ddX, RT;
    double RRT;
    int neq = X.giveSize();
    bool converged, errorOutOfRangeFlag;
    ParallelContext *parallel_context = engngModel->giveParallelContext( this->domain->giveNumber() );

    if ( engngModel->giveProblemScale() == macroScale ) {
        OOFEM_LOG_INFO("NRSolver: Iteration");
        if ( rtolf.at(1) > 0.0 ) {
            OOFEM_LOG_INFO(" ForceError");
        }
        if ( rtold.at(1) > 0.0 ) {
            OOFEM_LOG_INFO(" DisplError");
        }
        OOFEM_LOG_INFO("\n----------------------------------------------------------------------------\n");
    }

    l = 1.0;

    NM_Status status = NM_None;
    this->giveLinearSolver();

    // compute total load R = R+R0
    RT = R;
    if ( R0 ) {
        RT.add(* R0);
    }

    RRT = parallel_context->localNorm(RT);
    RRT *= RRT;

    ddX.resize(neq);
    ddX.zero();

    // Fetch the matrix before evaluating internal forces.
    // This is intentional, since its a simple way to drastically increase convergence for nonlinear problems.
    // (This old tangent is just used)
    // This improves convergence for many nonlinear problems, but not all. It may actually
    // cause divergence for some nonlinear problems. Therefore a flag is used to determine if
    // the stiffness should be evaluated before the residual (default yes). /ES

    engngModel->updateComponent(tStep, NonLinearLhs, domain);
    if ( this->prescribedDofsFlag ) {
        if ( !prescribedEqsInitFlag ) {
            this->initPrescribedEqs();
        }
        applyConstraintsToStiffness(k);
    }

    nite = 0;
    for ( nite = 0; ; ++nite ) {
        // Compute the residual
        engngModel->updateComponent(tStep, InternalRhs, domain);
        if (nite || iR == NULL) {
            rhs.beDifferenceOf(RT, F);
        } else {
            rhs = R;
            if (iR) {
                rhs.add(*iR); // add initial guess
            }
        }
        if ( this->prescribedDofsFlag ) {
            this->applyConstraintsToLoadIncrement(nite, k, rhs, rlm, tStep);
        }

        // convergence check
        converged = this->checkConvergence(RT, F, rhs, ddX, X, RRT, internalForcesEBENorm, nite, errorOutOfRangeFlag);

        if ( errorOutOfRangeFlag ) {
            status = NM_NoSuccess;
            OOFEM_WARNING("Divergence reached after %d iterations", nite);
            break;
        } else if ( converged && ( nite >= minIterations ) ) {
            break;
        } else if ( nite >= nsmax ) {
            OOFEM_LOG_DEBUG("Maximum number of iterations reached\n");
            break;
        }

        if ( nite > 0 || !mCalcStiffBeforeRes ) {
            if ( ( NR_Mode == nrsolverFullNRM ) || ( ( NR_Mode == nrsolverAccelNRM ) && ( nite % MANRMSteps == 0 ) ) ) {
                engngModel->updateComponent(tStep, NonLinearLhs, domain);
                applyConstraintsToStiffness(k);
            }
        }

        if ( ( nite == 0 ) && ( deltaL < 1.0 ) ) { // deltaL < 1 means no increment applied, only equilibrate current state
            rhs.zero();
            R.zero();
            ddX = rhs;
        } else {
            linSolver->solve(k, rhs, ddX);
        }

        //
        // update solution
        //
        if ( this->lsFlag && ( nite > 0 ) ) { // Why not nite == 0 ?
            // line search
            LineSearchNM :: LS_status LSstatus;
            double eta;
            this->giveLineSearchSolver()->solve(X, ddX, F, R, R0, prescribedEqs, 1.0, eta, LSstatus, tStep);
        } else if ( this->constrainedNRFlag && ( nite > this->constrainedNRminiter ) ) {
            ///@todo This doesn't check units, it is nonsense and must be corrected / Mikael
            if ( this->forceErrVec.computeSquaredNorm() > this->forceErrVecOld.computeSquaredNorm() ) {
                OOFEM_LOG_INFO("Constraining increment to be %e times full increment...\n", this->constrainedNRalpha);
                ddX.times(this->constrainedNRalpha);
            }
            //this->giveConstrainedNRSolver()->solve(X, & ddX, this->forceErrVec, this->forceErrVecOld, status, tStep);
        }
        X.add(ddX);
        dX.add(ddX);
        tStep->incrementStateCounter(); // update solution state counter
        tStep->incrementSubStepNumber();

        engngModel->giveExportModuleManager()->doOutput(tStep, true);
    }

    status |= NM_Success;
    solved = 1;

    // Modify Load vector to include "quasi reaction"
    if ( R0 ) {
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            R.at( prescribedEqs.at(i) ) = F.at( prescribedEqs.at(i) ) - R0->at( prescribedEqs.at(i) ) - R.at( prescribedEqs.at(i) );
        }
    } else {
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            R.at( prescribedEqs.at(i) ) = F.at( prescribedEqs.at(i) ) - R.at( prescribedEqs.at(i) );
        }
    }

    this->lastReactions.resize(numberOfPrescribedDofs);

#ifdef VERBOSE
    if ( numberOfPrescribedDofs ) {
        // print quasi reactions if direct displacement control used
        OOFEM_LOG_INFO("\n");
        OOFEM_LOG_INFO("NRSolver:     Quasi reaction table                                 \n");
        OOFEM_LOG_INFO("NRSolver:     Node            Dof             Displacement    Force\n");
        double reaction;
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            reaction = R.at( prescribedEqs.at(i) );
            if ( R0 ) {
                reaction += R0->at( prescribedEqs.at(i) );
            }
            lastReactions.at(i) = reaction;
            OOFEM_LOG_INFO("NRSolver:     %-15d %-15d %-+15.5e %-+15.5e\n", prescribedDofs.at(2 * i - 1), prescribedDofs.at(2 * i),
                           X.at( prescribedEqs.at(i) ), reaction);
        }
        OOFEM_LOG_INFO("\n");
    }
#endif

    return status;
}


SparseLinearSystemNM *
NRSolver :: giveLinearSolver()
{
    if ( linSolver ) {
        if ( linSolver->giveLinSystSolverType() == solverType ) {
            return linSolver.get();
        } else {
            linSolver.reset(NULL);
        }
    }

    linSolver.reset( classFactory.createSparseLinSolver(solverType, domain, engngModel) );
    if ( !linSolver ) {
        OOFEM_ERROR("linear solver creation failed for lstype %d", solverType);
    }

    return linSolver.get();
}


LineSearchNM *
NRSolver :: giveLineSearchSolver()
{
    if ( !linesearchSolver ) {
        linesearchSolver.reset( new LineSearchNM(domain, engngModel) );
    }

    return linesearchSolver.get();
}

void
NRSolver :: initPrescribedEqs()
{
    EModelDefaultEquationNumbering dn;
    ParallelContext *parallel_context = engngModel->giveParallelContext( this->domain->giveNumber() );
    int count = 0, ndofman = domain->giveNumberOfDofManagers();
    IntArray localPrescribedEqs(numberOfPrescribedDofs);

    for ( int j = 1; j <= ndofman; j++ ) {
        if ( !parallel_context->isLocal( domain->giveNode(j) ) ) {
            continue;
        }
        int jglobnum = domain->giveNode(j)->giveGlobalNumber();
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            int inode = prescribedDofs.at(2 * i - 1);
            int idofid = prescribedDofs.at(2 * i);
            if ( inode == jglobnum ) {
                localPrescribedEqs.at(++count) = domain->giveNode(j)->giveDofWithID(idofid)->giveEquationNumber(dn);
                continue;
            }
        }
    }

    prescribedEqs.resize(count);
    for ( int i = 1; i <= count; i++ ) {
        prescribedEqs.at(i) = localPrescribedEqs.at(i);
    }

    numberOfPrescribedDofs = count;

    this->prescribedEqsInitFlag = true;
}


void
NRSolver :: applyConstraintsToStiffness(SparseMtrx &k)
{
    if ( this->smConstraintVersion == k.giveVersion() ) {
        return;
    }

#ifdef __PETSC_MODULE
    PetscSparseMtrx *lhs = dynamic_cast< PetscSparseMtrx * >(&k);
    if ( lhs ) {
        Vec diag;
        PetscScalar *ptr;
        int eq;

        lhs->createVecGlobal(& diag);
        MatGetDiagonal(* lhs->giveMtrx(), diag);
        VecGetArray(diag, & ptr);
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            eq = prescribedEqs.at(i) - 1;
            MatSetValue(* ( lhs->giveMtrx() ), eq, eq, ptr [ eq ] * 1.e6, INSERT_VALUES);
        }

        MatAssemblyBegin(* lhs->giveMtrx(), MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(* lhs->giveMtrx(), MAT_FINAL_ASSEMBLY);
        VecRestoreArray(diag, & ptr);
        VecDestroy(& diag);
        if ( numberOfPrescribedDofs ) {
            this->smConstraintVersion = k.giveVersion();
        }

        return;
    }

#endif // __PETSC_MODULE
    for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
        k.at( prescribedEqs.at(i), prescribedEqs.at(i) ) *= 1.e6;
    }

    if ( numberOfPrescribedDofs ) {
        this->smConstraintVersion = k.giveVersion();
    }
}


void
NRSolver :: applyConstraintsToLoadIncrement(int nite, const SparseMtrx &k, FloatArray &R,
                                            referenceLoadInputModeType rlm, TimeStep *tStep)
{
    double factor = engngModel->giveDomain(1)->giveFunction(prescribedDisplacementTF)->evaluateAtTime( tStep->giveTargetTime() );
    if ( ( rlm == rlm_total ) && ( !tStep->isTheFirstStep() ) ) {
        //factor -= engngModel->giveDomain(1)->giveFunction(prescribedDisplacementTF)->
        // at(tStep->givePreviousStep()->giveTime()) ;
        factor -= engngModel->giveDomain(1)->giveFunction(prescribedDisplacementTF)->
        evaluateAtTime( tStep->giveTargetTime() - tStep->giveTimeIncrement() );
    }

    if ( nite == 0 ) {
#if 0
 #ifdef __PETSC_MODULE
        if ( solverType == ST_Petsc ) {
            //Natural2LocalOrdering* n2lpm = engngModel->giveParallelContext(1)->giveN2Lmap();
            //IntArray* map = n2lpm->giveN2Lmap();
            for ( i = 1; i <= prescribedEqs.giveSize(); i++ ) {
                eq = prescribedEqs.at(i);
                R.at(eq) = prescribedDofsValues.at(i) * factor; // local eq
            }

            return;
        }

 #endif
#else
 #ifdef __PETSC_MODULE
        const PetscSparseMtrx *lhs = dynamic_cast< const PetscSparseMtrx * >(&k);
        if ( lhs ) {
            Vec diag;
            PetscScalar *ptr;
            lhs->createVecGlobal(& diag);
            MatGetDiagonal(* ( const_cast< PetscSparseMtrx * >(lhs)->giveMtrx() ), diag);
            VecGetArray(diag, & ptr);

            for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
                int eq = prescribedEqs.at(i) - 1;
                R.at(eq + 1) = ptr [ eq ] * prescribedDofsValues.at(i) * factor;
            }

            VecRestoreArray(diag, & ptr);
            VecDestroy(& diag);
            return;
        }
 #endif
#endif
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            int eq = prescribedEqs.at(i);
            R.at(eq) = k.at(eq, eq) * prescribedDofsValues.at(i) * factor;
        }
    } else {
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            R.at( prescribedEqs.at(i) ) = 0.0;
        }
    }
}


void
NRSolver :: printState(FILE *outputStream)
{
#ifdef VERBOSE
    // print quasi reactions if direct displacement control used
    fprintf(outputStream, "\nQuasi reaction table:\n\n");
    fprintf(outputStream, "  node  dof            force\n");
    fprintf(outputStream, "============================\n");
    if ( lastReactions.giveSize() == 0 ) {
        return;
    }

    double reaction;
    for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
        reaction = lastReactions.at(i);
        fprintf(outputStream, "%6d  %3d   %+11.5e\n", prescribedDofs.at(2 * i - 1), prescribedDofs.at(2 * i), reaction);
    }
    fprintf(outputStream, "============================\n\n");
#endif
}


bool
NRSolver :: checkConvergence(FloatArray &RT, FloatArray &F, FloatArray &rhs,  FloatArray &ddX, FloatArray &X,
                             double RRT, const FloatArray &internalForcesEBENorm,
                             int nite, bool &errorOutOfRange)
{
    double forceErr, dispErr;
    FloatArray dg_forceErr, dg_dispErr, dg_totalLoadLevel, dg_totalDisp;
    bool answer;
    EModelDefaultEquationNumbering dn;
    ParallelContext *parallel_context = engngModel->giveParallelContext( this->domain->giveNumber() );

    /*
     * The force errors are (if possible) evaluated as relative errors.
     * If the norm of applied load vector is zero (one may load by temperature, etc)
     * then the norm of reaction forces is used in relative norm evaluation.
     *
     * Note: This is done only when all dofs are included (nccdg = 0). Not implemented if
     * multiple convergence criteria are used.
     *
     */

    answer = true;
    errorOutOfRange = false;

    // Store the errors associated with the dof groups    
    if ( this->constrainedNRFlag ) {
        this->forceErrVecOld = this->forceErrVec; // copy the old values
        this->forceErrVec.resize( internalForcesEBENorm.giveSize() );
        forceErrVec.zero();
    }

    if ( internalForcesEBENorm.giveSize() > 1 ) { // Special treatment when just one norm is given; No grouping
        int nccdg = this->domain->giveMaxDofID();
        // Keeps tracks of which dof IDs are actually in use;
        IntArray idsInUse(nccdg);
        idsInUse.zero();
        // zero error norms per group
        dg_forceErr.resize(nccdg);
        dg_forceErr.zero();
        dg_dispErr.resize(nccdg);
        dg_dispErr.zero();
        dg_totalLoadLevel.resize(nccdg);
        dg_totalLoadLevel.zero();
        dg_totalDisp.resize(nccdg);
        dg_totalDisp.zero();
        // loop over dof managers
        for ( auto &dofman : domain->giveDofManagers() ) {
            if ( !parallel_context->isLocal(dofman.get()) ) {
                continue;
            }

            // loop over individual dofs
            for ( Dof *dof: *dofman ) {
                if ( !dof->isPrimaryDof() ) {
                    continue;
                }
                int eq = dof->giveEquationNumber(dn);
                int dofid = dof->giveDofID();
                if ( !eq ) {
                    continue;
                }

                dg_forceErr.at(dofid) += rhs.at(eq) * rhs.at(eq);
                dg_dispErr.at(dofid) += ddX.at(eq) * ddX.at(eq);
                dg_totalLoadLevel.at(dofid) += RT.at(eq) * RT.at(eq);
                dg_totalDisp.at(dofid) += X.at(eq) * X.at(eq);
                idsInUse.at(dofid) = 1;
            } // end loop over DOFs
        } // end loop over dof managers

        // loop over elements and their DOFs
        for ( auto &elem : domain->giveElements() ) {
            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

            // loop over element internal Dofs
            for ( int idofman = 1; idofman <= elem->giveNumberOfInternalDofManagers(); idofman++ ) {
                DofManager *dofman = elem->giveInternalDofManager(idofman);
                // loop over individual dofs
                for ( Dof *dof: *dofman ) {
                    if ( !dof->isPrimaryDof() ) {
                        continue;
                    }
                    int eq = dof->giveEquationNumber(dn);
                    int dofid = dof->giveDofID();

                    if ( !eq ) {
                        continue;
                    }

                    dg_forceErr.at(dofid) += rhs.at(eq) * rhs.at(eq);
                    dg_dispErr.at(dofid) += ddX.at(eq) * ddX.at(eq);
                    dg_totalLoadLevel.at(dofid) += RT.at(eq) * RT.at(eq);
                    dg_totalDisp.at(dofid) += X.at(eq) * X.at(eq);
                    idsInUse.at(dofid) = 1;
                } // end loop over DOFs
            } // end loop over element internal dofmans
        } // end loop over elements

        // loop over boundary conditions and their internal DOFs
        for ( auto &bc : domain->giveBcs() ) {
            // loop over element internal Dofs
            for ( int idofman = 1; idofman <= bc->giveNumberOfInternalDofManagers(); idofman++ ) {
                DofManager *dofman = bc->giveInternalDofManager(idofman);
                // loop over individual dofs
                for ( Dof *dof: *dofman ) {
                    if ( !dof->isPrimaryDof() ) {
                        continue;
                    }
                    int eq = dof->giveEquationNumber(dn);
                    int dofid = dof->giveDofID();

                    if ( !eq ) {
                        continue;
                    }

                    dg_forceErr.at(dofid) += rhs.at(eq) * rhs.at(eq);
                    dg_dispErr.at(dofid) += ddX.at(eq) * ddX.at(eq);
                    dg_totalLoadLevel.at(dofid) += RT.at(eq) * RT.at(eq);
                    dg_totalDisp.at(dofid) += X.at(eq) * X.at(eq);
                    idsInUse.at(dofid) = 1;
                } // end loop over DOFs
            } // end loop over element internal dofmans
        } // end loop over elements

        // exchange individual partition contributions (simultaneously for all groups)
        FloatArray collectiveErr(nccdg);
        parallel_context->accumulate(dg_forceErr,       collectiveErr);
        dg_forceErr       = collectiveErr;
        parallel_context->accumulate(dg_dispErr,        collectiveErr);
        dg_dispErr        = collectiveErr;
        parallel_context->accumulate(dg_totalLoadLevel, collectiveErr);
        dg_totalLoadLevel = collectiveErr;
        parallel_context->accumulate(dg_totalDisp,      collectiveErr);
        dg_totalDisp      = collectiveErr;

        OOFEM_LOG_INFO("NRSolver: %-5d", nite);
        //bool zeroNorm = false;
        // loop over dof groups and check convergence individually
        for ( int dg = 1; dg <= nccdg; dg++ ) {
            bool zeroFNorm = false, zeroDNorm = false;
            // Skips the ones which aren't used in this problem (the residual will be zero for these anyway, but it is annoying to print them all)
            if ( !idsInUse.at(dg) ) {
                continue;
            }

            OOFEM_LOG_INFO( "  %s:", __DofIDItemToString( ( DofIDItem ) dg ).c_str() );

            if ( rtolf.at(1) > 0.0 ) {
                //  compute a relative error norm
                if ( ( dg_totalLoadLevel.at(dg) + internalForcesEBENorm.at(dg) ) > nrsolver_ERROR_NORM_SMALL_NUM ) {
                    forceErr = sqrt( dg_forceErr.at(dg) / ( dg_totalLoadLevel.at(dg) + internalForcesEBENorm.at(dg) ) );
                } else {
                    // If both external forces and internal ebe norms are zero, then the residual must be zero.
                    //zeroNorm = true; // Warning about this afterwards.
                    zeroFNorm = true;
                    forceErr = sqrt( dg_forceErr.at(dg) );
                }

                if ( forceErr > rtolf.at(1) * NRSOLVER_MAX_REL_ERROR_BOUND ) {
                    errorOutOfRange = true;
                }
                if ( forceErr > rtolf.at(1) ) {
                    answer = false;
                }
                OOFEM_LOG_INFO(zeroFNorm ? " *%.3e" : "  %.3e", forceErr);

                // Store the errors from the current iteration
                if ( this->constrainedNRFlag ) {
                    forceErrVec.at(dg) = forceErr;
                }       
            }

            if ( rtold.at(1) > 0.0 ) {
                // compute displacement error
                if ( dg_totalDisp.at(dg) >  nrsolver_ERROR_NORM_SMALL_NUM ) {
                    dispErr = sqrt( dg_dispErr.at(dg) / dg_totalDisp.at(dg) );
                } else {
                    ///@todo This is almost always the case for displacement error. nrsolveR_ERROR_NORM_SMALL_NUM is no good.
                    //zeroNorm = true; // Warning about this afterwards.
                    //zeroDNorm = true;
                    dispErr = sqrt( dg_dispErr.at(dg) );
                }
                if ( dispErr  > rtold.at(1) * NRSOLVER_MAX_REL_ERROR_BOUND ) {
                    errorOutOfRange = true;
                }
                if ( dispErr > rtold.at(1) ) {
                    answer = false;
                }
                OOFEM_LOG_INFO(zeroDNorm ? " *%.3e" : "  %.3e", dispErr);
            }
        }
        OOFEM_LOG_INFO("\n");
        //if ( zeroNorm ) OOFEM_WARNING("Had to resort to absolute error measure (marked by *)");
    } else { // No dof grouping
        double dXX, dXdX;

        if ( engngModel->giveProblemScale() == macroScale ) {
            OOFEM_LOG_INFO("NRSolver:     %-15d", nite);
        } else {
            OOFEM_LOG_INFO("  NRSolver:     %-15d", nite);
        }


        forceErr = parallel_context->localNorm(rhs);
        forceErr *= forceErr;
        dXX = parallel_context->localNorm(X);
        dXX *= dXX;                                       // Note: Solutions are always total global values (natural distribution makes little sense for the solution)
        dXdX = parallel_context->localNorm(ddX);
        dXdX *= dXdX;

        if ( rtolf.at(1) > 0.0 ) {
            // we compute a relative error norm
            if ( ( RRT + internalForcesEBENorm.at(1) ) > nrsolver_ERROR_NORM_SMALL_NUM ) {
                forceErr = sqrt( forceErr / ( RRT + internalForcesEBENorm.at(1) ) );
            } else {
                forceErr = sqrt(forceErr);   // absolute norm as last resort
            }
            if ( fabs(forceErr) > rtolf.at(1) * NRSOLVER_MAX_REL_ERROR_BOUND ) {
                errorOutOfRange = true;
            }
            if ( fabs(forceErr) > rtolf.at(1) ) {
                answer = false;
            }
            OOFEM_LOG_INFO(" %-15e", forceErr);

            if ( this->constrainedNRFlag ) {
                // store the errors from the current iteration for use in the next
                forceErrVec.at(1) = forceErr;
            }
        }

        if ( rtold.at(1) > 0.0 ) {
            // compute displacement error
            // err is relative displacement change
            if ( dXX > nrsolver_ERROR_NORM_SMALL_NUM ) {
                dispErr = sqrt(dXdX / dXX);
            } else {
                dispErr = sqrt(dXdX);
            }
            if ( fabs(dispErr)  > rtold.at(1) * NRSOLVER_MAX_REL_ERROR_BOUND ) {
                errorOutOfRange = true;
            }
            if ( fabs(dispErr)  > rtold.at(1) ) {
                answer = false;
            }
            OOFEM_LOG_INFO(" %-15e", dispErr);
        }

        OOFEM_LOG_INFO("\n");
    } // end default case (all dofs contributing)

    return answer;
}
} // end namespace oofem
