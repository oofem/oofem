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
#include "generalboundarycondition.h"
#include "dof.h"

namespace oofem {

#define ERROR_NORM_SMALL_NUM 1.e-6
#define MAX_REL_ERROR_BOUND 1.e20

REGISTER_SparseNonLinearSystemNM(StaggeredSolver)

StaggeredSolver :: StaggeredSolver(Domain *d, EngngModel *m) : NRSolver(d, m)
{
    this->UnknownNumberingSchemeList.resize(0);
}


IRResultType
StaggeredSolver :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = NRSolver ::initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    IR_GIVE_FIELD(ir, this->totalIdList, _IFT_StaggeredSolver_DofIdList);
    IR_GIVE_FIELD(ir, this->idPos, _IFT_StaggeredSolver_DofIdListPositions);
    
    this->instanciateYourself();
    
    return IRRT_OK;
}


void 
StaggeredSolver :: instanciateYourself()
{
    int numDofIdGroups = idPos.giveSize()/2;
    this->UnknownNumberingSchemeList.resize(numDofIdGroups);
    IntArray idList; 
    for ( int i = 0; i < numDofIdGroups; i++ ) {
        int sz = idPos.at(i*2 + 2) - idPos.at(i*2 + 1) + 1;
        idList.resize(sz);
        for ( int j = 1; j <= sz; j++) {
            int pos = idPos.at(i*2 + 1) + j - 1;
            idList.at(j) = totalIdList.at(pos);
        }
        this->UnknownNumberingSchemeList[i].setDofIdArray(idList);
        this->UnknownNumberingSchemeList[i].setNumber(i+1);
  
    }    
    
    // Allocate stiffness matrices and internal forces vectors corresponding to each dof group
   
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

   
void
StaggeredSolver :: giveTotalLocationArray(IntArray &condensedLocationArray, const UnknownNumberingScheme &s, Domain *d)
{
    IntArray nodalArray, ids, locationArray;
    locationArray.clear();
    
    for ( auto &dman : d->giveDofManagers() ) {
        dman->giveCompleteLocationArray(nodalArray, s);
        locationArray.followedBy(nodalArray);
    }
    for ( auto &elem : d->giveElements() ) {
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
StaggeredSolver :: solve(SparseMtrx &k, FloatArray &R, FloatArray *R0, FloatArray *iR,
                  FloatArray &Xtotal, FloatArray &dXtotal, FloatArray &F,
                  const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                  int &nite, TimeStep *tStep)

{
    // residual, iteration increment of solution, total external force
    FloatArray RHS, rhs, ddXtotal, RT;
    double RRTtotal;
    int neq = Xtotal.giveSize();
    NM_Status status;
    bool converged, errorOutOfRangeFlag;
    ParallelContext *parallel_context = engngModel->giveParallelContext( this->domain->giveNumber() );
        
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
    RT = R;
    if ( R0 ) {
        RT.add(* R0);
    }

    RRTtotal = parallel_context->localNorm(RT);
    RRTtotal *= RRTtotal;

    ddXtotal.resize(neq);
    ddXtotal.zero();

    // Fetch the matrix before evaluating internal forces.
    // This is intentional, since its a simple way to drastically increase convergence for nonlinear problems.
    // (This old tangent is just used)
    // This improves convergence for many nonlinear problems, but not all. It may actually
    // cause divergence for some nonlinear problems. Therefore a flag is used to determine if
    // the stiffness should be evaluated before the residual (default yes). /ES
    
    //-------------------------------------------------    
  
    // Compute external forces 
    int numDofIdGroups = (int)this->UnknownNumberingSchemeList.size();
    FloatArray RRT(numDofIdGroups);
    for ( int dG = 0; dG < numDofIdGroups; dG++ ) {
        this->fExtList[dG].beSubArrayOf( RT, locArrayList[dG] );        
        RRT(dG) = this->fExtList[dG].computeSquaredNorm();
    }

    for (int nStaggeredIter = 0;; ++nStaggeredIter) {

        // Staggered iterations
        for ( int dG = 0; dG < (int)this->UnknownNumberingSchemeList.size(); dG++ ) {
            printf("\nSolving for dof group %d \n", dG+1);
            
            engngModel->updateComponent(tStep, NonLinearLhs, domain);      
            this->stiffnessMatrixList[dG].reset(k.giveSubMatrix( locArrayList[dG], locArrayList[dG]));

            if ( this->prescribedDofsFlag ) {
                if ( !prescribedEqsInitFlag ) {
                    this->initPrescribedEqs();
                }
                applyConstraintsToStiffness(k);
            }

            for (nite = 0;; ++nite) {
                // Compute the residual
                engngModel->updateComponent(tStep, InternalRhs, domain);
                RHS.beDifferenceOf(RT, F); 
                
                this->fIntList[dG].beSubArrayOf( F, locArrayList[dG] );
                rhs.beDifferenceOf(this->fExtList[dG], this->fIntList[dG]);
                
                RHS.zero();
                RHS.assemble(rhs, locArrayList[dG]);
                
                
                if ( this->prescribedDofsFlag ) {
                    this->applyConstraintsToLoadIncrement(nite, k, rhs, rlm, tStep);
                }

                // convergence check
                IntArray &idArray = UnknownNumberingSchemeList[dG].dofIdArray;
                converged = this->checkConvergenceDofIdArray(RT, F, RHS, ddXtotal, Xtotal, RRTtotal, internalForcesEBENorm, nite, errorOutOfRangeFlag, tStep, idArray);
               

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
                        this->stiffnessMatrixList[dG].reset(k.giveSubMatrix( locArrayList[dG], locArrayList[dG]));
                        applyConstraintsToStiffness(*this->stiffnessMatrixList[dG]);
                    }
                }

                if ( ( nite == 0 ) && ( deltaL < 1.0 ) ) { // deltaL < 1 means no increment applied, only equilibrate current state
                    rhs.zero();
                    R.zero();
                    ddX[dG] = rhs;
                } else {
                    status = linSolver->solve(*this->stiffnessMatrixList[dG], rhs, ddX[dG]);
                }

                //
                // update solution
                //
                if ( this->lsFlag && ( nite > 0 ) ) { // Why not nite == 0 ?
                    // line search
                    LineSearchNM :: LS_status LSstatus;
                    double eta;
                    this->giveLineSearchSolver()->solve( X[dG], ddX[dG], fIntList[dG], fExtList[dG], R0, prescribedEqs, 1.0, eta, LSstatus, tStep);
                } else if ( this->constrainedNRFlag && ( nite > this->constrainedNRminiter ) ) { 
                    if ( this->forceErrVec.computeSquaredNorm() > this->forceErrVecOld.computeSquaredNorm() ) {
                        printf("Constraining increment to be %e times full increment...\n", this->constrainedNRalpha);
                        ddX[dG].times(this->constrainedNRalpha);
                    }
                }
                X[dG].add(ddX[dG]);
                dX[dG].add(ddX[dG]);

                
                // Update total solution (containing all dofs)
                Xtotal.assemble(ddX[dG], locArrayList[dG]);
                dXtotal.assemble(ddX[dG], locArrayList[dG]);
                ddXtotal.zero();
                ddXtotal.assemble(ddX[dG], locArrayList[dG]);
                
                tStep->incrementStateCounter(); // update solution state counter
                tStep->incrementSubStepNumber();

                engngModel->giveExportModuleManager()->doOutput(tStep, true);
            }
        }


        printf("\nStaggered iteration (all dof id's) \n");
        
        // Check convergence of total system
        RHS.beDifferenceOf(RT, F);
        converged = this->checkConvergence(RT, F, RHS, ddXtotal, Xtotal, RRTtotal, internalForcesEBENorm, nStaggeredIter, errorOutOfRangeFlag);
        if ( converged && ( nStaggeredIter >= minIterations ) ) {
            break;
        }
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


//copied from nrsolver

bool
StaggeredSolver :: checkConvergenceDofIdArray(FloatArray &RT, FloatArray &F, FloatArray &rhs, FloatArray &ddX, FloatArray &X,
                          double RRT, const FloatArray &internalForcesEBENorm, int nite, bool &errorOutOfRange, TimeStep *tStep, IntArray &dofIdArray)
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

        OOFEM_LOG_INFO("StaggeredSolver: %-5d", nite);
        //bool zeroNorm = false;
        // loop over dof groups and check convergence individually
        for ( int dg = 1; dg <= nccdg; dg++ ) {
            bool zeroFNorm = false, zeroDNorm = false;
            // Skips the ones which aren't used in this problem (the residual will be zero for these anyway, but it is annoying to print them all)
            if ( !idsInUse.at(dg) ) {
                continue;
            }
            
            if ( dofIdArray.giveSize() && !dofIdArray.contains(dg) ) {
                continue;
            }

            
            OOFEM_LOG_INFO( "  %s:", __DofIDItemToString( ( DofIDItem ) dg ).c_str() );

            if ( rtolf.at(1) > 0.0 ) {
                //  compute a relative error norm
                if ( ( dg_totalLoadLevel.at(dg) + internalForcesEBENorm.at(dg) ) > ERROR_NORM_SMALL_NUM ) {
                    forceErr = sqrt( dg_forceErr.at(dg) / ( dg_totalLoadLevel.at(dg) + internalForcesEBENorm.at(dg) ) );
                } else {
                    // If both external forces and internal ebe norms are zero, then the residual must be zero.
                    //zeroNorm = true; // Warning about this afterwards.
                    zeroFNorm = true;
                    forceErr = sqrt( dg_forceErr.at(dg) );
                }

                if ( forceErr > rtolf.at(1) * MAX_REL_ERROR_BOUND ) {
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
                if ( dg_totalDisp.at(dg) >  ERROR_NORM_SMALL_NUM ) {
                    dispErr = sqrt( dg_dispErr.at(dg) / dg_totalDisp.at(dg) );
                } else {
                    ///@todo This is almost always the case for displacement error. ERROR_NORM_SMALL_NUM is no good.
                    //zeroNorm = true; // Warning about this afterwards.
                    //zeroDNorm = true;
                    dispErr = sqrt( dg_dispErr.at(dg) );
                }
                if ( dispErr  > rtold.at(1) * MAX_REL_ERROR_BOUND ) {
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
            OOFEM_LOG_INFO("StaggeredSolver:     %-15d", nite);
        } else {
            OOFEM_LOG_INFO("  StaggeredSolver:     %-15d", nite);
        }

        forceErr = parallel_context->localNorm(rhs);
        forceErr *= forceErr;
        dXX = parallel_context->localNorm(X);
        dXX *= dXX;                                       // Note: Solutions are always total global values (natural distribution makes little sense for the solution)
        dXdX = parallel_context->localNorm(ddX);
        dXdX *= dXdX;

        if ( rtolf.at(1) > 0.0 ) {
            // we compute a relative error norm
            if ( ( RRT + internalForcesEBENorm.at(1) ) > ERROR_NORM_SMALL_NUM ) {
                forceErr = sqrt( forceErr / ( RRT + internalForcesEBENorm.at(1) ) );
            } else {
                forceErr = sqrt(forceErr);   // absolute norm as last resort
            }
            if ( fabs(forceErr) > rtolf.at(1) * MAX_REL_ERROR_BOUND ) {
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
            if ( dXX > ERROR_NORM_SMALL_NUM ) {
                dispErr = sqrt(dXdX / dXX);
            } else {
                dispErr = sqrt(dXdX);
            }
            if ( fabs(dispErr)  > rtold.at(1) * MAX_REL_ERROR_BOUND ) {
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

