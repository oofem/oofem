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

#include "timestep.h"
#include "classfactory.h"
#include "exportmodulemanager.h"
#include "engngm.h"
#include "domain.h"
#include "dofmanager.h"
#include "element.h"


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
     IntArray idList1 = {1, 2, 3, 15, 16, 17, 18};
     IntArray idList2 = {18};
     IntArray idList3;
     idList1 = {1, 2, 3, 15, 16, 17};
     idList2 = {17, 18};
     idList3 = {1, 2, 3, 15, 16, 17, 18};
     this->UnknownNumberingSchemeList.resize(3);
     this->UnknownNumberingSchemeList[0].setDofIdArray(idList1);
     this->UnknownNumberingSchemeList[0].setNumber(1);
     this->UnknownNumberingSchemeList[1].setDofIdArray(idList2);
     this->UnknownNumberingSchemeList[1].setNumber(2);
     this->UnknownNumberingSchemeList[2].setDofIdArray(idList3);
     this->UnknownNumberingSchemeList[2].setNumber(3);

//     IntArray idList1 = {1, 2, 3};
//     this->UnknownNumberingSchemeList.resize(1);
//     this->UnknownNumberingSchemeList[0].setDofIdArray(idList1);
//     this->UnknownNumberingSchemeList[0].setNumber(1);    
    
    return IRRT_OK;
}




   
void
StaggeredSolver :: giveTotalLocationArray(IntArray &condensedLocationArray, const UnknownNumberingScheme &s, Domain *d)
{
    IntArray masterDofIDs, nodalArray, ids, locationArray;
    locationArray.clear();
    
    for ( int i = 1; i <= d->giveNumberOfDofManagers(); i++ ) {
        d->giveDofManager(i)->giveCompleteLocationArray(nodalArray, s);
        locationArray.followedBy(nodalArray);
    }
    for ( int el = 1; el <= d->giveNumberOfElements(); el++ ) {
        Element *elem = d->giveElement(el);
        for ( int i = 1; i <= elem->giveNumberOfInternalDofManagers(); i++ ) {
            elem->giveInternalDofManDofIDMask(i, ids);
            elem->giveInternalDofManager(i)->giveLocationArray(ids, nodalArray, s);
            locationArray.followedBy(nodalArray);
        }
    }
    
    
    IntArray nonZeroMask;
    nonZeroMask.findNonzeros(locationArray);

    condensedLocationArray.resize(nonZeroMask.giveSize());
    for ( int i = 1; i <= nonZeroMask.giveSize(); i++ ) {
        condensedLocationArray.at(i) = locationArray.at( nonZeroMask.at(i) );    
    }
}    


NM_Status
StaggeredSolver :: solve(SparseMtrx *k, FloatArray *R, FloatArray *R0,
                  FloatArray *Xtotal, FloatArray *dXtotal, FloatArray *F,
                  const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                  int &nite, TimeStep *tNow)

{
    // residual, iteration increment of solution, total external force
    FloatArray RHS, rhs, ddXtotal, RT;
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
    
    int numDofIdGroups = this->UnknownNumberingSchemeList.size();
    FloatArray RRT(numDofIdGroups);
    
   
    
    // Create stiffness matrices and internal forces vectors corresponding to each dof group
    if ( !this->fIntList.size() ) {
        this->fIntList.resize(numDofIdGroups);
        this->fExtList.resize(numDofIdGroups);
        this->locArrayList.resize(numDofIdGroups);            
        this->stiffnessMatrixList.resize(numDofIdGroups);
        this->X.resize(numDofIdGroups);
        this->dX.resize(numDofIdGroups);
        this->ddX.resize(numDofIdGroups);

        for ( int dG = 0; dG < numDofIdGroups; dG++ ) {
            this->giveTotalLocationArray(this->locArrayList[dG], UnknownNumberingSchemeList[dG], domain);         
            int neq = locArrayList[dG].giveSize();
            this->fIntList[dG].resize(neq);    
            this->fExtList[dG].resize(neq);    
            
            this->X[dG].resize(neq);
            this->X[dG].zero();
            this->dX[dG].resize(neq);
            this->dX[dG].zero();            
            this->ddX[dG].resize(neq);
            this->ddX[dG].zero();
         
        }
    }
    
    
    // Compute external forces 
    for ( int dG = 0; dG < numDofIdGroups; dG++ ) {
        this->fExtList[dG].beSubArrayOf( RT, locArrayList[dG] );        
        RRT(dG) = this->fExtList[dG].computeSquaredNorm();
    }
    
  int nStaggeredIter = 0;
  do {
    // Staggered iterations
    for ( int dG = 0; dG < (int)this->UnknownNumberingSchemeList.size(); dG++ ) {
        printf("\n Solving for dof group %d \n\n", dG+1);
        
        engngModel->updateComponent(tNow, NonLinearLhs, domain);      
        this->stiffnessMatrixList[dG] = k->giveSubMatrix( locArrayList[dG], locArrayList[dG]);

        if ( this->prescribedDofsFlag ) {
            if ( !prescribedEqsInitFlag ) {
                this->initPrescribedEqs();
            }
            applyConstraintsToStiffness(k);
        }

        nite = 0;
        do {
            // Compute the residual
            engngModel->updateComponent(tNow, InternalRhs, domain);
            RHS.beDifferenceOf(RT, * F); 
            
            this->fIntList[dG].beSubArrayOf( *F, locArrayList[dG] );
            rhs.beDifferenceOf(this->fExtList[dG], this->fIntList[dG]);
            
            RHS.zero();
            RHS.assemble(rhs, locArrayList[dG]);
            
            
            if ( this->prescribedDofsFlag ) {
                this->applyConstraintsToLoadIncrement(nite, k, rhs, rlm, tNow);
            }

            // convergence check
            converged = this->checkConvergence(RT, * F, RHS, ddXtotal, * Xtotal, RRTtotal, internalForcesEBENorm, nite, errorOutOfRangeFlag, tNow);
            //converged = this->checkConvergence(fExtList[dG], fIntList[dG], rhs, ddX[dG], X[dG], RRT(dG), internalForcesEBENorm, nite, errorOutOfRangeFlag, tNow);

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
                    engngModel->updateComponent(tNow, NonLinearLhs, domain);
                    this->stiffnessMatrixList[dG] = k->giveSubMatrix( locArrayList[dG], locArrayList[dG]);
                    applyConstraintsToStiffness(this->stiffnessMatrixList[dG]);
                }
            }

            if ( ( nite == 0 ) && ( deltaL < 1.0 ) ) { // deltaL < 1 means no increment applied, only equilibrate current state
                rhs.zero();
                R->zero();
                ddX[dG] = rhs;
            } else {
                status = linSolver->solve(this->stiffnessMatrixList[dG], & rhs, & ddX[dG]);
            }

            //
            // update solution
            //
            if ( this->lsFlag && ( nite > 0 ) ) { // Why not nite == 0 ?
                // line search
                LineSearchNM :: LS_status status;
                double eta;               
                this->giveLineSearchSolver()->solve( &X[dG], &ddX[dG], &fIntList[dG], &fExtList[dG], R0, prescribedEqs, 1.0, eta, status, tNow);
            } else if ( this->constrainedNRFlag && ( nite > this->constrainedNRminiter ) ) { 
                if ( this->forceErrVec.computeSquaredNorm() > this->forceErrVecOld.computeSquaredNorm() ) {
                    printf("Constraining increment to be %e times full increment...\n", this->constrainedNRalpha);
                    ddX[dG].times(this->constrainedNRalpha);
                }   
            }
            X[dG].add(ddX[dG]);
            dX[dG].add(ddX[dG]);

            
            // Update total solution 
            Xtotal->assemble(ddX[dG], locArrayList[dG]);
            dXtotal->assemble(ddX[dG], locArrayList[dG]);
            ddXtotal.zero();
            ddXtotal.assemble(ddX[dG], locArrayList[dG]);
            
            tNow->incrementStateCounter(); // update solution state counter
            tNow->incrementSubStepNumber();
            nite++; // iteration increment

            engngModel->giveExportModuleManager()->doOutput(tNow, true);
        } while ( true ); // end of iteration
    }
    
    
    printf("\n Staggered iteration \n");
    
    // Check convergence of total system
    RHS.beDifferenceOf(RT, * F);
    converged = this->checkConvergence(RT, * F, RHS, ddXtotal, * Xtotal, RRTtotal, internalForcesEBENorm, nStaggeredIter, errorOutOfRangeFlag, tNow);
    if ( converged && ( nStaggeredIter >= minIterations ) ) {
        break;
    }    
    
    nStaggeredIter++;
    
   } while ( true ); // end of iteration 

   
   
   
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

