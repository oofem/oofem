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


#include "staggeredsolver.h"
#include "EngineeringModels/staticstructuralstaggered.h" ///@todo temporary for development

#include "timestep.h"
#include "classfactory.h"
#include "exportmodulemanager.h"


namespace oofem {
#define nrsolver_ERROR_NORM_SMALL_NUM 1.e-6
#define NRSOLVER_MAX_REL_ERROR_BOUND 1.e20
#define NRSOLVER_MAX_RESTARTS 4
#define NRSOLVER_RESET_STEP_REDUCE 0.25
#define NRSOLVER_DEFAULT_NRM_TICKS 10

REGISTER_SparseNonLinearSystemNM(StaggeredSolver)

StaggeredSolver :: StaggeredSolver(Domain *d, EngngModel *m) : NRSolver(d, m)
{
 
    this->UnknownNumberingSchemeList.resize(0);
    
}


IRResultType
StaggeredSolver :: initializeFrom(InputRecord *ir)
{
    //IRResultType result;                // Required by IR_GIVE_FIELD macro

    NRSolver ::initializeFrom(ir);

    // set the dofIdArrays for each solution field to solve for
    IntArray idList1 = {1, 2};
    IntArray idList2 = {3};
    this->UnknownNumberingSchemeList.resize(2);
    this->UnknownNumberingSchemeList[0].setDofIdArray(idList1);
    this->UnknownNumberingSchemeList[0].setNumber(1);
    this->UnknownNumberingSchemeList[1].setDofIdArray(idList2);
    this->UnknownNumberingSchemeList[1].setNumber(2);
    
    
    return IRRT_OK;
}



NM_Status
StaggeredSolver :: solve(SparseMtrx *k, FloatArray *R, FloatArray *R0,
                  FloatArray *Xtotal, FloatArray *dXtotal, FloatArray *F,
                  const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                  int &nite, TimeStep *tNow)

{
    // residual, iteration increment of solution, total external force
    FloatArray rhs, ddXtotal, RT;
    double RRTtotal;
    int neq = Xtotal->giveSize();
    NM_Status status;
    bool converged, errorOutOfRangeFlag;
#ifdef __PARALLEL_MODE
    ParallelContext *parallel_context = engngModel->giveParallelContext( this->domain->giveNumber() );
#endif
        
    if ( engngModel->giveProblemScale() == macroScale ) {
        OOFEM_LOG_INFO("StaggeredSolver: Iteration");
        if ( rtolf.at(1) > 0.0 ) {
            OOFEM_LOG_INFO(" ForceError");
        }
        if ( rtold.at(1) > 0.0 ) {
            OOFEM_LOG_INFO(" DisplError");
        }
        OOFEM_LOG_INFO("\n----------------------------------------------------------------------------\n");
    }

    l = 1.0;

    status = NM_None;
    this->giveLinearSolver();

    // compute total load R = R+R0
    RT = * R;
    if ( R0 ) {
        RT.add(* R0);
    }

#ifdef __PARALLEL_MODE
    RRTtotal = parallel_context->localNorm(RT);
    RRTtotal *= RRTtotal;
#else
    RRTtotal = RT.computeSquaredNorm();
#endif

    ddXtotal.resize(neq);
    ddXtotal.zero();

    // Fetch the matrix before evaluating internal forces.
    // This is intentional, since its a simple way to drastically increase convergence for nonlinear problems.
    // (This old tangent is just used)
    // This improves convergence for many nonlinear problems, but not all. It may actually
    // cause divergence for some nonlinear problems. Therefore a flag is used to determine if
    // the stiffness should be evaluated before the residual (default yes). /ES

    
    //-------------------------------------------------    
    
    StaticStructuralStaggered *ss = static_cast< StaticStructuralStaggered* > (engngModel);
    int numDofIdGroups = this->UnknownNumberingSchemeList.size();
    FloatArray X, dX, ddX; // unknowns for the current dof group
    FloatArray RRT(numDofIdGroups);
    int di = 1; 
    
    // Create stiffness matrices and internal forces vectors corresponding to each dof group
    if ( !this->stiffnessMatrixList.size() ) {
        this->fIntList.resize(numDofIdGroups);
        this->fExtList.resize(numDofIdGroups);            
        this->stiffnessMatrixList.resize(numDofIdGroups);

        for ( int dG = 0; dG < numDofIdGroups; dG++ ) {
            int neq = ss->giveNumberOfDomainEquations( di, UnknownNumberingSchemeList[dG] );
            this->fIntList[dG].resize(neq);    
            this->fExtList[dG].resize(neq);    
            
            SparseMtrx *mtrx = this->stiffnessMatrixList[dG];
            mtrx = classFactory.createSparseMtrx(ss->giveSparseMtrxType()); 
            
            if ( !mtrx ) {
            OOFEM_ERROR("Couldn't create requested sparse matrix of type %d", ss->giveSparseMtrxType() );
            }
            
            mtrx->buildInternalStructure( engngModel, di, UnknownNumberingSchemeList[dG] );
            //mtrx->buildInternalStructure( engngModel, di, EModelDefaultEquationNumbering() );
            mtrx->printStatistics();
            mtrx->printYourself();
            
        }
    }
    
    
    // Compute external forces 
    for ( int dG = 0; dG < numDofIdGroups; dG++ ) {
        ss->updateExternalForcesForStaggeredSolver(this->fExtList[dG], tNow, domain, this->UnknownNumberingSchemeList[dG]);
        RRT(dG) = this->fExtList[dG].computeSquaredNorm();
        
        this->fExtList[dG].printYourself("f_ext");
    }
    
    //TODO add outer iterations loop
    // Staggered iterations
    for ( int dG = 0; dG < (int)this->UnknownNumberingSchemeList.size(); dG++ ) {
      
        //engngModel->updateComponent(tNow, NonLinearLhs, domain);
        //TODO crashes here when zeroing the matrices
        ss->updateTangentStiffnessForStaggeredSolver(this->stiffnessMatrixList[dG], tNow, domain, this->UnknownNumberingSchemeList[dG]);
        if ( this->prescribedDofsFlag ) {
            if ( !prescribedEqsInitFlag ) {
                this->initPrescribedEqs();
            }
            applyConstraintsToStiffness(k);
        }

        nite = 0;
        do {
            // Compute the residual
            // engngModel->updateComponent(tNow, InternalRhs, domain);
            ss->updateInternalForcesForStaggeredSolver(this->fIntList[dG], tNow, domain, this->UnknownNumberingSchemeList[dG]);
            //rhs.beDifferenceOf(RT, * F); 
            rhs.beDifferenceOf(this->fExtList[dG], this->fIntList[dG]);

            if ( this->prescribedDofsFlag ) {
                this->applyConstraintsToLoadIncrement(nite, k, rhs, rlm, tNow);
            }

            // convergence check
            //converged = this->checkConvergence(RT, * F, rhs, ddXtotal, * Xtotal, RRTtotal, internalForcesEBENorm, nite, errorOutOfRangeFlag, tNow);
            converged = this->checkConvergence(this->fExtList[dG], this->fIntList[dG], rhs, ddX, X, RRT.at(dG), internalForcesEBENorm, nite, errorOutOfRangeFlag, tNow);

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
                    //engngModel->updateComponent(tNow, NonLinearLhs, domain);
                    ss->updateTangentStiffnessForStaggeredSolver(this->stiffnessMatrixList[dG], tNow, domain, this->UnknownNumberingSchemeList[dG]);
                    applyConstraintsToStiffness(this->stiffnessMatrixList[dG]);
                }
            }

            if ( ( nite == 0 ) && ( deltaL < 1.0 ) ) { // deltaL < 1 means no increment applied, only equilibrate current state
                rhs.zero();
                R->zero();
                ddX = rhs;
            } else {
                linSolver->solve(this->stiffnessMatrixList[dG], & rhs, & ddX);
            }

            //
            // update solution
            //
            if ( this->lsFlag && ( nite > 0 ) ) { // Why not nite == 0 ?
                // line search
                LineSearchNM :: LS_status status;
                double eta;
                this->giveLineSearchSolver()->solve(Xtotal, & ddX, F, R, R0, prescribedEqs, 1.0, eta, status, tNow);
            } else if ( this->constrainedNRFlag && ( nite > this->constrainedNRminiter ) ) { 
                if ( this->forceErrVec.computeSquaredNorm() > this->forceErrVecOld.computeSquaredNorm() ) {
                    printf("Constraining increment to be %e times full increment...\n", this->constrainedNRalpha);
                    ddX.times(this->constrainedNRalpha);
                }   
            }
            X.add(ddX);
            dX.add(ddX);
            tNow->incrementStateCounter(); // update solution state counter
            tNow->incrementSubStepNumber();
            nite++; // iteration increment

            engngModel->giveExportModuleManager()->doOutput(tNow, true);
        } while ( true ); // end of iteration
    }
    status |= NM_Success;
    solved = 1;

    // Modify Load vector to include "quasi reaction"
    if ( R0 ) {
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            R->at( prescribedEqs.at(i) ) = F->at( prescribedEqs.at(i) ) - R0->at( prescribedEqs.at(i) ) - R->at( prescribedEqs.at(i) );
        }
    } else {
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            R->at( prescribedEqs.at(i) ) = F->at( prescribedEqs.at(i) ) - R->at( prescribedEqs.at(i) );
        }
    }

    this->lastReactions.resize(numberOfPrescribedDofs);

#ifdef VERBOSE
    if ( numberOfPrescribedDofs ) {
        // print quasi reactions if direct displacement control used
        OOFEM_LOG_INFO("\n");
        OOFEM_LOG_INFO("StaggeredSolver:     Quasi reaction table                                 \n");
        OOFEM_LOG_INFO("StaggeredSolver:     Node            Dof             Displacement    Force\n");
        double reaction;
        for ( int i = 1; i <= numberOfPrescribedDofs; i++ ) {
            reaction = R->at( prescribedEqs.at(i) );
            if ( R0 ) {
                reaction += R0->at( prescribedEqs.at(i) );
            }
            lastReactions.at(i) = reaction;
            OOFEM_LOG_INFO("StaggeredSolver:     %-15d %-15d %-+15.5e %-+15.5e\n", prescribedDofs.at(2 * i - 1), prescribedDofs.at(2 * i),
                           X.at( prescribedEqs.at(i) ), reaction);
        }
        OOFEM_LOG_INFO("\n");
    }
#endif

    return status;
}



} // end namespace oofem

