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

#include "calmls.h"
#include "verbose.h"
#include "timestep.h"
#include "flotmtrx.h"
#include "datastream.h"
#include "mathfem.h"
#include "element.h"
// includes for HPC - not very clean (NumMethod knows what is "node" and "dof")
#include "node.h"
#include "dof.h"
#include "usrdefsub.h"
#include "contextioerr.h"

namespace oofem {
#define CALM_RESET_STEP_REDUCE 0.25
#define CALM_MAX_RESTARTS 4
#define CALM_TANGENT_STIFF_TRESHOLD 0.1
#define CALM_DEFAULT_NRM_TICKS 2
#define CALM_MAX_REL_ERROR_BOUND 1.e10

CylindricalALM :: CylindricalALM(int i, Domain *d, EngngModel *m, EquationID ut) :
    SparseNonLinearSystemNM(i, d, m, ut), calm_HPCWeights(), calm_HPCIndirectDofMask(), calm_HPCDmanDofSrcArray(), ccDofGroups()
{
    //
    // constructor
    //

    nsmax  = 60;       // default maximum number of sweeps allowed
    numberOfRequiredIterations = 3;
    //rtol   = 10.E-3  ;   // convergence tolerance
    //Psi    = 0.1;       // displacement control on
    Psi    = 1.0;     // load control on
    solved = 0;
    calm_NR_Mode = calm_NR_OldMode = calm_modifiedNRM;
    calm_NR_ModeTick = -1; // do not swith to calm_NR_OldMode
    calm_MANRMSteps = 0;

    //Bergan_k0 = 0.;    // value used for computing Bergan's parameter
    // of current stiffness.

    deltaL    = -1.0;
    // TangenStiffnessTreshold = 0.10;

    // Variables for Hyper Plane Control
    calm_Control = calm_hpc_off; // HPControl is not default
    linSolver = NULL;
    // linesearch default off
    lsFlag = 0;

    ls_tolerance = 0.40;
    amplifFactor = 2.5;
    maxEta = 8.5;
    minEta = 0.2;

    // number of convergence_criteria dof groups, set to 0 (default behavior)
    nccdg = 0;

    minIterations = 0;

#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
    parallel_context = engngModel->givePetscContext(d->giveNumber(), ut);
#endif
#endif

}


CylindricalALM ::  ~CylindricalALM()
{
    //
    // destructor
    //
    delete linSolver;
}


NM_Status
CylindricalALM :: solve(SparseMtrx *k, FloatArray *R, FloatArray *R0,
                        FloatArray *X, FloatArray *dX, FloatArray *F,
                        double &internalForcesEBENorm, double &ReachedLambda, referenceLoadInputModeType rlm,
                        int &nite, TimeStep *tNow)
{
    FloatArray rhs, deltaXt, deltaX_, dXm1, XInitial;
    FloatArray ddX; // total increment of displacements in iteration
    //double Bergan_k0 = 1.0, bk;
    double XX, RR, RR0, XR, p = 0.0;
    double deltaLambda, Lambda, eta, DeltaLambdam1, DeltaLambda = 0.0;
    double drProduct = 0.0;
    int neq = R->giveSize();
    int irest = 0;
    int HPsize, i, ind;
    double _RR, _XX;
    NM_Status status;
    bool converged, errorOutOfRangeFlag;
    // print iteration header

    OOFEM_LOG_INFO("CALMLS:       Initial step length: %-15e\n", deltaL);
    if ( nccdg == 0 ) {
        OOFEM_LOG_INFO("CALMLS:       Iteration       LoadLevel       ForceError      DisplError    \n");
        OOFEM_LOG_INFO("----------------------------------------------------------------------------\n");
    } else {
        OOFEM_LOG_INFO("Iter  LoadLevel       ");
        for ( i = 1; i <= nccdg; i++ ) {
            OOFEM_LOG_INFO("ForceError(%02d)  DisplError(%02d)  ", i, i);
        }

        OOFEM_LOG_INFO("\n__________________________________________________________\n");
    }

    //
    // Now smarter method is default:
    // after convergence troubles for (calm_NR_ModeTick subsequent
    // steps newly set calm_NR_Mode will be used. After the calm_NR_OldMode will be restored
    //
    if ( calm_NR_ModeTick == 0 ) {
        calm_NR_Mode = calm_NR_OldMode;
    }

    if ( calm_NR_ModeTick > 0 ) {
        calm_NR_ModeTick--;
    }

    XInitial = * X;
    ddX.resize(neq);
    deltaXt.resize(neq);

    status = NM_None;
    this->giveLinearSolver();

    // create HPC Map if needed
    if ( calm_hpc_init ) {
        this->convertHPCMap();
        calm_hpc_init = 0;
    }

    if ( R0 ) {
#ifdef __PARALLEL_MODE
        RR0 = parallel_context->localNorm(*R0);
        RR0 *= RR0;
#else
        RR0 = R0->computeSquaredNorm();
#endif
    } else {
        RR0 = 0.0;
    }

    //
    // A  initial step (predictor)
    //
    //

    //
    // A.2.   We assume positive-definitive (0)Kt (tangent stiffness mtrx).
    //
#ifdef __PARALLEL_MODE
    RR = parallel_context->localNorm(*R);
    RR *= RR;
#else
    RR = R->computeSquaredNorm();
#endif

restart:
    //
    // A.1. calculation of (0)VarRt
    //
    dX->zero();
    //engngModel->updateComponent(tNow, InternalRhs, domain); // By not updating this, one obtains the old equilibrated tangent.
    engngModel->updateComponent(tNow, NonLinearLhs, domain);
    linSolver->solve(k, R, & deltaXt);

    if ( calm_Control == calm_hpc_off ) {
#ifdef __PARALLEL_MODE
        XX = parallel_context->localNorm(deltaXt); XX *= XX;
#else
        XX = deltaXt.computeSquaredNorm();
#endif
        p = sqrt(XX + Psi * Psi * RR);
    } else if ( calm_Control == calm_hpc_on ) {
        HPsize = calm_HPCIndirectDofMask.giveSize();
        _XX = 0;
        _RR = 0;
        for ( i = 1; i <= HPsize; i++ ) {
            if ( ( ind = calm_HPCIndirectDofMask.at(i) ) != 0 ) {
                _XX += deltaXt.at(ind) * deltaXt.at(ind);
                _RR += R->at(ind) * R->at(ind);
            }
        }

#ifdef __PARALLEL_MODE
        FloatArray XXRR(2), collected_XXRR;
        XXRR(0) = _XX;
        XXRR(1) = _RR;

        parallel_context->accumulate(XXRR, collected_XXRR);
        _XX = collected_XXRR(0);
        _RR = collected_XXRR(1);
#endif

        p = sqrt(_XX + Psi * Psi * _RR);
    } else if ( calm_Control == calml_hpc ) {
        HPsize = calm_HPCIndirectDofMask.giveSize();
        p = 0.;
        for ( i = 1; i <= HPsize; i++ ) {
            if ( ( ind = calm_HPCIndirectDofMask.at(i) ) != 0 ) {
                p += deltaXt.at(ind) * calm_HPCWeights.at(i);
            }
        }

#ifdef __PARALLEL_MODE
        p = parallel_context->accumulate(p);
#endif
    }

#ifdef __PARALLEL_MODE
    XR = parallel_context->localDotProduct(deltaXt, *R);
#else
    XR = deltaXt.dotProduct(*R);
#endif
    /* XR is unscaled Bergan's param of current stiffness XR = deltaXt^T k deltaXt
     * this is used to test whether k has negative or positive slope */

    Lambda = ReachedLambda;
    DeltaLambda = deltaLambda = sgn(XR) * deltaL / p;
    Lambda += DeltaLambda;
    //
    // A.3.
    //

    rhs = * R;
    rhs.times(DeltaLambda);
    if ( R0 ) {
        rhs.add(*R0);
    }
    linSolver->solve(k, & rhs, dX);
    X->add(*dX);

    nite = 0;

    // update solution state counter
    tNow->incrementStateCounter();
    engngModel->updateComponent(tNow, InternalRhs, domain);

    do {
        nite++;
        dXm1 = * dX;
        DeltaLambdam1 = DeltaLambda;

        //
        // B  - iteration MNRM is used
        //
        // B.1. is ommited because MNRM is used instead of NRM.
        //
        if ( ( calm_NR_Mode == calm_fullNRM ) || ( ( calm_NR_Mode == calm_accelNRM ) && ( nite % calm_MANRMSteps == 0 ) ) ) {
            //
            // ALM with full NRM
            //
            // we assemble new tangent stiffness and compute new deltaXt
            // internal state of elements is updated by previous calling
            // InternalRhs
            //
            engngModel->updateComponent(tNow, NonLinearLhs, domain);
            //
            // compute deltaXt for i-th iteration
            //
            linSolver->solve(k, R, & deltaXt);
        }

        // B.2.
        //

        rhs =  * R;
        rhs.times(Lambda);
        if ( R0 ) {
            rhs.add(*R0);
        }

        rhs.subtract(*F);
        deltaX_.resize(neq);
        linSolver->solve(k, & rhs, & deltaX_);
        eta = 1.0;
        //
        // B.3.
        //
        if ( this->computeDeltaLambda(deltaLambda, * dX, deltaXt, deltaX_, * R, RR, eta, deltaL, DeltaLambda, neq) ) {
            irest++;
            if ( irest <= CALM_MAX_RESTARTS ) {
                // convergence problems
                // there must be step restart followed by decrease of step length
                // status |= NM_ForceRestart;
                // reduce step length
                deltaL =  deltaL * CALM_RESET_STEP_REDUCE;
                if ( deltaL < minStepLength ) {
                    deltaL = minStepLength;
                }

                // restore previous total displacement vector
                *X = XInitial;
                // reset all changes fro previous equilibrium state
                dX->zero();
                // restore initial stiffness
                engngModel->updateComponent(tNow, NonLinearLhs, domain);

                OOFEM_LOG_INFO("CALMLS:       Iteration Reset ...\n");

                calm_NR_OldMode  = calm_NR_Mode;
                calm_NR_Mode     = calm_fullNRM;
                calm_NR_ModeTick = CALM_DEFAULT_NRM_TICKS;
                goto restart;
            } else {
                status = NM_NoSuccess;
                OOFEM_ERROR("CALMLS :: solve - can't continue further");
            }
        }

        //
        // B.4.+ B.5.
        //


        if ( this->lsFlag && ( nite != 1 ) ) {
            //
            //  LINE SEARCH
            //
            this->do_lineSearch(* X, XInitial, deltaX_, deltaXt,
                                dXm1, * dX, ddX,
                                * R, R0, * F,
                                DeltaLambda, DeltaLambdam1, deltaLambda, Lambda, ReachedLambda,
                                RR, drProduct, tNow);
        } else { // no line search
            //
            // update solution vectors
            //
            ddX.resize(0);
            ddX.add(eta*deltaLambda, deltaXt);
            ddX.add(eta, deltaX_);
            *dX = dXm1;
            dX->add(ddX);
            *X = XInitial;
            X->add(*dX);

#ifdef __PARALLEL_MODE
            drProduct = parallel_context->localNorm(ddX);
            drProduct *= drProduct;
#else
            drProduct = ddX.computeSquaredNorm();
#endif
            tNow->incrementStateCounter();     // update solution state counter
            //
            // B.6.
            //
            DeltaLambda = DeltaLambdam1 + deltaLambda;
            Lambda = ReachedLambda + DeltaLambda;

            tNow->incrementStateCounter();     // update solution state counter
            engngModel->updateComponent(tNow, InternalRhs, domain);
        }

        //
        // B.7.
        //
        // convergence check
        //

        converged = this->checkConvergence(* R, R0, * F, * X, ddX, Lambda, RR0, RR, drProduct,
                                           internalForcesEBENorm, nite, errorOutOfRangeFlag);
        if ( ( nite >= nsmax ) || errorOutOfRangeFlag ) {
            irest++;
            if ( irest <= CALM_MAX_RESTARTS ) {
                // convergence problems
                // there must be step restart followed by decrease of step length
                // status |= NM_ForceRestart;
                // reduce step length
                deltaL =  deltaL * CALM_RESET_STEP_REDUCE;
                if ( deltaL < minStepLength ) {
                    deltaL = minStepLength;
                }

                // restore previous total displacement vector
                *X = XInitial;
                // reset all changes from previous equilibrium state
                engngModel->initStepIncrements();
                dX->zero();
                // restore initial stiffness
                engngModel->updateComponent(tNow, NonLinearLhs, domain);

                OOFEM_LOG_INFO("CALMLS:       Iteration Reset ...\n");

                calm_NR_OldMode  = calm_NR_Mode;
                calm_NR_Mode     = calm_fullNRM;
                calm_NR_ModeTick = CALM_DEFAULT_NRM_TICKS;
                goto restart;
            } else {
                status = NM_NoSuccess;
                OOFEM_WARNING2("CALMLS :: solve - Convergence not reached after %d iterations", nsmax);
                // exit(1);
                break;
            }
        }
    } while ( !converged || ( nite < minIterations ) );

    //
    // update dofs,nodes,Elemms and print result
    //
#ifdef VERBOSE
    // printf ("\nCALM - step iteration finished") ;
#endif

    // Seems to be completely unused / Mikael
#if 0
    /* compute Bergan's parameter of current stiffness */
#ifdef __PARALLEL_MODE
    double RDR = parallel_context->localDotProduct(*R,*dX);
    double DR = parallel_context->localNorm(*dX);
    bk = DeltaLambda * RDR/(DR*DR);
#else
    bk = DeltaLambda * dX->dotProduct(*R) / dX->computeSquaredNorm();
#endif

    if ( tNow->giveNumber() == 1 ) {
        // compute Bergan_k0
        Bergan_k0 = bk;
    } else {
        bk = bk / Bergan_k0;
        //if (fabs(bk) < CALM_TANGENT_STIFF_TRESHOLD) status |= NM_KeepTangent;
    }
#endif

    /* we have computed Bergan's  parameter -> you can adjust deltaL according
     * to this value */
    //
    // there has been restart already - set nite to maxiter
    //
    // if (irest > 0) nite = nsmax;

    if ( nite > numberOfRequiredIterations ) {
        deltaL =  deltaL * numberOfRequiredIterations / nite;
    } else {
        deltaL =  deltaL * sqrt( sqrt( ( double ) numberOfRequiredIterations / ( double ) nite ) );
    }

    if ( deltaL > maxStepLength ) {
        deltaL = maxStepLength;
    }

    if ( deltaL < minStepLength ) {
        deltaL = minStepLength;
    }

    OOFEM_LOG_INFO("CALMLS:       Adjusted step length: %-15e\n", deltaL);

    status = NM_Success;
    solved = 1;
    ReachedLambda = Lambda;

    return status;
}

bool
CylindricalALM :: checkConvergence(const FloatArray &R, const FloatArray *R0, const FloatArray &F,
                                   const FloatArray &X, const FloatArray &ddX,
                                   double Lambda, double RR0, double RR, double drProduct,
                                   double internalForcesEBENorm, int nite, bool &errorOutOfRange)
{
    /*
     * typedef std::set<DofID> __DofIDSet;
     * std::list<__DofIDSet> __ccDofGroups;
     * int nccdg; // number of Convergence Criteria Dof Groups
     */
    int _dg, _idofman, _ielem, _idof, _eq, _ndof, _ng = nccdg, ndofman = domain->giveNumberOfDofManagers();
    int nelem = domain->giveNumberOfElements();
    double forceErr, dispErr, _val;
    DofManager *_idofmanptr;
    Element *_ielemptr;
    Dof *_idofptr;
    FloatArray rhs; // residual of momentum balance eq (unbalanced nodal forces)
    FloatArray dg_forceErr(nccdg), dg_dispErr(nccdg), dg_totalLoadLevel(nccdg), dg_totalDisp(nccdg);
    bool answer;
    EModelDefaultEquationNumbering dn;
#ifdef __PARALLEL_MODE
 #ifdef __PETSC_MODULE
    PetscNatural2LocalOrdering *n2l = parallel_context->giveN2Lmap();
 #endif
#endif

    answer = true;
    errorOutOfRange = false;

    // compute residual vector
    rhs =  R;
    rhs.times(Lambda);
    if ( R0 ) {
        rhs.add(*R0);
    }

    rhs.subtract(F);

    if ( _ng > 0 ) {
        forceErr = dispErr = 0.0;
        // zero error norms per group
        dg_forceErr.zero();
        dg_dispErr.zero();
        dg_totalLoadLevel.zero();
        dg_totalDisp.zero();
        // loop over dof managers
        for ( _idofman = 1; _idofman <= ndofman; _idofman++ ) {
            _idofmanptr = domain->giveDofManager(_idofman);
#ifdef __PARALLEL_MODE
            if ( !_idofmanptr->isLocal() ) {
                continue;
            }

#endif

            _ndof = _idofmanptr->giveNumberOfDofs();
            // loop over individual dofs
            for ( _idof = 1; _idof <= _ndof; _idof++ ) {
                _idofptr = _idofmanptr->giveDof(_idof);
                // loop over dof groups
                for ( _dg = 1; _dg <= _ng; _dg++ ) {
                    // test if dof ID is in active set
                    if ( ccDofGroups.at(_dg - 1).find( _idofptr->giveDofID() ) != ccDofGroups.at(_dg - 1).end() ) {
                        _eq = _idofptr->giveEquationNumber(dn);

                        if ( _eq ) {
#ifdef __PARALLEL_MODE
                            if ( !n2l->giveNewEq(_eq) ) {
                                continue;
                            }
#endif

                            _val = rhs.at(_eq);
                            dg_forceErr.at(_dg) += _val * _val;
                            _val = ddX.at(_eq);
                            dg_dispErr.at(_dg)  += _val * _val;
                            // missing - compute norms of total displacement and load vectors (but only for selected dofs)!
                            if ( R0 ) {
                                _val = R0->at(_eq);
                                dg_totalLoadLevel.at(_dg) += _val * _val;
                            }

                            _val = R.at(_eq);
                            dg_totalLoadLevel.at(_dg) += _val * _val * Lambda * Lambda;
                            _val = X.at(_eq);
                            dg_totalDisp.at(_dg) += _val * _val;
                        }
                    }
                } // end loop over dof groups

            } // end loop over DOFs

        } // end loop over dof managers

        // loop over elements and their DOFs
        for ( _ielem = 1; _ielem <= nelem; _ielem++ ) {
            _ielemptr = domain->giveElement(_ielem);
#ifdef __PARALLEL_MODE
            if ( _ielemptr->giveParallelMode() != Element_local ) {
                continue;
            }

#endif
            // loop over element internal Dofs
            for (_idofman=1; _idofman<=_ielemptr->giveNumberOfInternalDofManagers(); _idofman++) {
                _ndof = _ielemptr->giveInternalDofManager(_idofman)->giveNumberOfDofs();
                // loop over individual dofs
                for ( _idof = 1; _idof <= _ndof; _idof++ ) {
                    _idofptr = _ielemptr->giveInternalDofManager(_idofman)->giveDof(_idof);
                    // loop over dof groups
                    for ( _dg = 1; _dg <= _ng; _dg++ ) {
                        // test if dof ID is in active set
                        if ( ccDofGroups.at(_dg - 1).find( _idofptr->giveDofID() ) != ccDofGroups.at(_dg - 1).end() ) {
                            _eq = _idofptr->giveEquationNumber(dn);

                            if ( _eq ) {
#ifdef __PARALLEL_MODE
                                if ( !n2l->giveNewEq(_eq) ) {
                                    continue;
                                }
#endif

                                _val = rhs.at(_eq);
                                dg_forceErr.at(_dg) += _val * _val;
                                _val = ddX.at(_eq);
                                dg_dispErr.at(_dg)  += _val * _val;
                                // missing - compute norms of total displacement and load vectors (but only for selected dofs)!
                                if ( R0 ) {
                                    _val = R0->at(_eq);
                                    dg_totalLoadLevel.at(_dg) += _val * _val;
                                }

                                _val = R.at(_eq);
                                dg_totalLoadLevel.at(_dg) += _val * _val * Lambda * Lambda;
                                _val = X.at(_eq);
                                dg_totalDisp.at(_dg) += _val * _val;
                            }
                        }
                    } // end loop over dof groups
                } // end loop over DOFs
            } // end loop over internal element dofmans
        } // end loop over elements

#ifdef __PARALLEL_MODE
        // exchange individual partition contributions (simultaneously for all groups)
        FloatArray collectiveErr(_ng);
        parallel_context->accumulate(dg_forceErr, collectiveErr);
        dg_forceErr = collectiveErr;
        parallel_context->accumulate(dg_dispErr, collectiveErr);
        dg_dispErr = collectiveErr;
        parallel_context->accumulate(dg_totalLoadLevel, collectiveErr);
        dg_totalLoadLevel = collectiveErr;
        parallel_context->accumulate(dg_totalDisp, collectiveErr);
        dg_totalDisp = collectiveErr;
#endif

        OOFEM_LOG_INFO("CALMLS:       %-15d %-15e ", nite, Lambda);
        // loop over dof groups
        for ( _dg = 1; _dg <= _ng; _dg++ ) {
            //  compute a relative error norm
            if ( ( dg_totalLoadLevel.at(_dg) ) < calm_SMALL_ERROR_NUM ) {
                dg_forceErr.at(_dg) = sqrt( dg_forceErr.at(_dg) );
            } else {
                dg_forceErr.at(_dg) = sqrt( dg_forceErr.at(_dg) / dg_totalLoadLevel.at(_dg) );
            }

            //
            // compute displacement error
            //
            if ( dg_totalDisp.at(_dg) < calm_SMALL_ERROR_NUM ) {
                dg_dispErr.at(_dg) = sqrt( dg_dispErr.at(_dg) );
            } else {
                dg_dispErr.at(_dg) = sqrt( dg_dispErr.at(_dg) / dg_totalDisp.at(_dg) );
            }

            if ( ( fabs( dg_forceErr.at(_dg) ) > rtolf.at(_dg) * CALM_MAX_REL_ERROR_BOUND ) ||
                ( fabs( dg_dispErr.at(_dg) )  > rtold.at(_dg) * CALM_MAX_REL_ERROR_BOUND ) ) {
                errorOutOfRange = true;
            }

            if ( ( fabs( dg_forceErr.at(_dg) ) > rtolf.at(_dg) ) || ( fabs( dg_dispErr.at(_dg) ) > rtold.at(_dg) ) ) {
                answer = false;
            }


            OOFEM_LOG_INFO("%-15e %-15e ", dg_forceErr.at(_dg), dg_dispErr.at(_dg) );
        }

        OOFEM_LOG_INFO("\n");
    } else {
        //
        // _ng==0 (errors computed for all dofs - this is the default)
        //

        //
        // compute force error(s)
        //
        double dXX;
#ifdef __PARALLEL_MODE
        forceErr = parallel_context->localNorm(rhs); forceErr *= forceErr;
        dXX = parallel_context->localNorm(X); dXX *= dXX;
#else
        // err is relative error of unbalanced forces
        forceErr = rhs.computeSquaredNorm();
        // err is relative displacement change
        dXX = X.computeSquaredNorm();

#endif
        // we compute a relative error norm
        if ( ( RR0 + RR * Lambda * Lambda ) > calm_SMALL_ERROR_NUM ) {
            forceErr = sqrt( forceErr / ( RR0 + RR * Lambda * Lambda ) );
        } else if ( internalForcesEBENorm > calm_SMALL_ERROR_NUM ) {
            forceErr = sqrt(forceErr / internalForcesEBENorm);
        } else {
            forceErr = sqrt(forceErr);
        }

        //
        // compute displacement error
        //
        if ( dXX < calm_SMALL_ERROR_NUM ) {
            dispErr = drProduct;
        } else {
            dispErr = drProduct / dXX;
            dispErr = sqrt(dispErr);
        }

        if ( ( fabs(forceErr) > rtolf.at(1) * CALM_MAX_REL_ERROR_BOUND ) ||
            ( fabs(dispErr)  > rtold.at(1) * CALM_MAX_REL_ERROR_BOUND ) ) {
            errorOutOfRange = true;
        }

        if ( ( fabs(forceErr) > rtolf.at(1) ) || ( fabs(dispErr) > rtold.at(1) ) ) {
            answer = false;
        }

        OOFEM_LOG_INFO("CALMLS:       %-15d %-15e %-15e %-15e\n", nite, Lambda, forceErr, dispErr);
    } // end default case (all dofs conributing)

    return answer;
}



IRResultType
CylindricalALM :: initializeFrom(InputRecord *ir)
//
//
//
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    double oldPsi =  Psi; // default from constructor
    double initialStepLength, forcedInitialStepLength;
    int hpcMode;
    IRResultType val;

    IR_GIVE_OPTIONAL_FIELD(ir, Psi, IFT_CylindricalALM_psi, "psi"); // Macro
    if ( Psi < 0. ) {
        Psi = oldPsi;
    }

    // double oldTangenStiffnessTreshold = TangenStiffnessTreshold;
    // TangenStiffnessTreshold = readDouble (initString,"tstiffnesstreshold");
    // if (TangenStiffnessTreshold < 0.01) TangenStiffnessTreshold  = oldTangenStiffnessTreshold;
    nsmax = 30;
    IR_GIVE_OPTIONAL_FIELD(ir, nsmax, IFT_CylindricalALM_maxiter, "maxiter"); // Macro
    if ( nsmax < 30 ) {
        nsmax = 30;
    }

    minStepLength = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, minStepLength, IFT_CylindricalALM_minsteplength, "minsteplength"); // Macro

    IR_GIVE_FIELD(ir, maxStepLength, IFT_CylindricalALM_steplength, "steplength"); // Macro
    initialStepLength = maxStepLength;
    IR_GIVE_OPTIONAL_FIELD(ir, initialStepLength, IFT_CylindricalALM_initialsteplength, "initialsteplength"); // Macro

    //if (deltaL <= 0.0)  deltaL=maxStepLength;
    // This method (instanciate) is called not only at the beginning but also
    // after restart from engngModel updateAttributes -> and in this case
    // we want to keep restored deltaL)
    if ( ( deltaL <= 0.0 ) || ( deltaL > maxStepLength ) ) {
        deltaL = initialStepLength;
    }

    // using this, one can enforce deltaL from the input file after restart
    forcedInitialStepLength = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, forcedInitialStepLength, IFT_CylindricalALM_forcedinitialsteplength, "forcedinitialsteplength"); // Macro
    if ( forcedInitialStepLength > 0. ) {
        deltaL = forcedInitialStepLength;
    }

    numberOfRequiredIterations = 3;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfRequiredIterations, IFT_CylindricalALM_reqiterations, "reqiterations"); // Macro
    if ( numberOfRequiredIterations < 3 ) {
        numberOfRequiredIterations = 3;
    }

    if ( numberOfRequiredIterations > 1000 ) {
        numberOfRequiredIterations = 1000;
    }

    val = IR_GIVE_OPTIONAL_FIELD(ir, minIterations, IFT_CylindricalALM_miniterations, "miniter"); // Macro
    if ( val == IRRT_OK ) {
        if ( minIterations > 3 && minIterations < 1000 ) {
            numberOfRequiredIterations = minIterations;
        }

        if ( nsmax <= minIterations ) {
            nsmax = minIterations + 1;
        }
    }

    // read if MANRM method is used
    calm_MANRMSteps = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, calm_MANRMSteps, IFT_CylindricalALM_manrmsteps, "manrmsteps"); // Macro
    if ( calm_MANRMSteps > 0 ) {
        calm_NR_Mode = calm_NR_OldMode = calm_accelNRM;
    } else {
        calm_NR_Mode = calm_modifiedNRM;
    }

    // read if HPC is requsted
    hpcMode = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, hpcMode, IFT_CylindricalALM_hpcmode, "hpcmode"); // Macro

    calm_HPCDmanDofSrcArray.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, calm_HPCDmanDofSrcArray, IFT_CylindricalALM_hpc, "hpc"); // Macro

    calm_HPCDmanWeightSrcArray.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, calm_HPCDmanWeightSrcArray, IFT_CylindricalALM_hpcw, "hpcw"); // Macro
    // in calm_HPCIndirectDofMask are stored pairs with following meaning:
    // inode idof
    // example HPC 4 1 2 6 1
    // will yield to HPC for node 4 dof 1 and node 6 dof 1
    // calm_HPCIndirectDofMask must be converted to indirect map
    // -> because dof eqs. are not known now, we derefer this to
    // solveYourselfAt() subroutine. The need for converting is indicated by
    // calm_HPControl = hpc_init
    if ( calm_HPCDmanDofSrcArray.giveSize() != 0 ) {
        int i, nsize;
        if ( hpcMode == 1 ) {
            calm_Control = calm_hpc_on;
        } else if ( hpcMode == 2 ) {
            calm_Control = calml_hpc;
        } else {
            calm_Control = calm_hpc_on; // default is to use hpc_on
        }

        if ( ( calm_HPCDmanDofSrcArray.giveSize() % 2 ) != 0 ) {
            _error("HPC Map size must be even number, it contains pairs <node, nodeDof>");
        }

        nsize = calm_HPCDmanDofSrcArray.giveSize() / 2;
        if ( calm_HPCWeights.giveSize() == 0 ) {
            // no weights -> set to 1.0 by default
            calm_HPCWeights.resize(nsize);
            for ( i = 1; i <= nsize; i++ ) {
                calm_HPCWeights.at(i) = 1.0;
            }
        } else if ( nsize != calm_HPCWeights.giveSize() ) {
            _error("HPC map size and weight array size mismatch");
        }

        calm_hpc_init = 1;
    } else {
        if ( hpcMode ) {
            _error("HPC Map must be specified");
        }
    }

    int _value = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _value, IFT_CylindricalALM_lstype, "lstype"); // Macro
    solverType = ( LinSystSolverType ) _value;

    this->giveLinearSolver()->initializeFrom(ir);

    this->lsFlag = 0; // linesearch default is off
    IR_GIVE_OPTIONAL_FIELD(ir, lsFlag, IFT_CylindricalALM_linesearch, "linesearch"); // Macro

    IR_GIVE_OPTIONAL_FIELD(ir, ls_tolerance, IFT_CylindricalALM_lsearchtol, "lsearchtol"); // Macro
    if ( ls_tolerance < 0.6 ) {
        ls_tolerance = 0.6;
    }

    if ( ls_tolerance > 0.95 ) {
        ls_tolerance = 0.95;
    }

    IR_GIVE_OPTIONAL_FIELD(ir, amplifFactor, IFT_CylindricalALM_lsearchamp, "lsearchamp"); // Macro
    if ( amplifFactor < 1.0 ) {
        amplifFactor = 1.0;
    }

    if ( amplifFactor > 10.0 ) {
        amplifFactor = 10.0;
    }

    IR_GIVE_OPTIONAL_FIELD(ir, maxEta, IFT_CylindricalALM_lsearchmaxeta, "lsearchmaxeta"); // Macro
    if ( maxEta < 1.5 ) {
        maxEta = 1.5;
    }

    if ( maxEta > 15.0 ) {
        maxEta = 15.0;
    }

    /* initialize optional dof groups for convergence criteria evaluation */
    this->nccdg = 0; // default, no dof cc group, all norms evaluated for all dofs
    IR_GIVE_OPTIONAL_FIELD(ir, nccdg, IFT_CylindricalALM_nccdg, "nccdg"); // Macro

    if ( nccdg >= 1 ) {
        int _i, _j;
        IntArray _val;
        char name [ 12 ];
        // create an empty set
        __DofIDSet _set;
        // resize gof group vector
        this->ccDofGroups.resize(nccdg, _set);
        for ( _i = 0; _i < nccdg; _i++ ) {
            sprintf(name, "ccdg%d", _i + 1);
            // read dof group as int array under ccdg# keyword
            IR_GIVE_FIELD(ir, _val, IFT_CylindricalALM_ccdg, name); // Macro
            // convert aray into set
            for ( _j = 1; _j <= _val.giveSize(); _j++ ) {
                ccDofGroups.at(_i).insert( (DofIDItem)_val.at(_j) );
            }
        }

        // read relative error tolerances of the solver fo each cc
        // if common rtolv provided, set to this tolerace both rtolf and rtold
        IR_GIVE_OPTIONAL_FIELD(ir, rtolf, IFT_CylindricalALM_rtolv, "rtolv"); // Macro
        rtold = rtolf;
        // read optional force and displacement tolerances
        IR_GIVE_OPTIONAL_FIELD(ir, rtolf, IFT_CylindricalALM_rtolf, "rtolf"); // Macro
        IR_GIVE_OPTIONAL_FIELD(ir, rtold, IFT_CylindricalALM_rtold, "rtold"); // Macro

        if ( ( rtolf.giveSize() != nccdg ) || ( rtold.giveSize() != nccdg ) ) {
            _error2("Incompatible size of rtolf or rtold params, expected size %d (nccdg)", nccdg);
        }
    } else {
        nccdg = 0;
        double _rtol = 1.e-3; // default tolerance
        rtolf.resize(1);
        rtold.resize(1);
        // read relative error tolerances of the solver
        // if common rtolv provided, set to this tolerace both rtolf and rtold
        IR_GIVE_OPTIONAL_FIELD(ir, _rtol, IFT_CylindricalALM_rtolf, "rtolv"); // Macro
        rtolf.at(1) = rtold.at(1) = _rtol;
        IR_GIVE_OPTIONAL_FIELD(ir, _rtol, IFT_CylindricalALM_rtolf, "rtolf"); // Macro
        rtolf.at(1) = _rtol;
        IR_GIVE_OPTIONAL_FIELD(ir, _rtol, IFT_CylindricalALM_rtold, "rtold"); // Macro
        rtold.at(1) = _rtol;
    }


    //refLoadInputMode = (calm_referenceLoadInputModeType) readInteger  (initString, "refloadmode");
    return IRRT_OK;
}


void CylindricalALM  :: convertHPCMap()
//
// converts HPC map from user input map to HPC indirect Map
//
// indirect map:
// size of indirect map is number of controlled DOFs
// on i-th position this map contain equation number of i-th controlled dof
//
// user input map;
// user input map has size 2*controlled DOFs
// and contain pairs (node, nodeDof);
//
// This is used in order to hide equation numbering from user
//
{
    IntArray indirectMap;
    FloatArray weights;
    int size, i;
    int inode, idof;
    EModelDefaultEquationNumbering dn;

#if defined ( __PARALLEL_MODE ) || defined ( __ENABLE_COMPONENT_LABELS )
    int j, jglobnum, count = 0, ndofman = domain->giveNumberOfDofManagers();
    size = calm_HPCDmanDofSrcArray.giveSize() / 2;
    indirectMap.resize(size);
    weights.resize(size);
    for ( j = 1; j <= ndofman; j++ ) {
        jglobnum = domain->giveNode(j)->giveLabel();
        for ( i = 1; i <= size; i++ ) {
            inode = calm_HPCDmanDofSrcArray.at(2 * i - 1);
            idof  = calm_HPCDmanDofSrcArray.at(2 * i);
            if ( inode == jglobnum ) {
 #ifdef __PARALLEL_MODE
                // HUHU hard wired domain no 1
                if ( parallel_context->isLocal( domain->giveNode(j) ) ) {
                    indirectMap.at(++count) = domain->giveNode(j)->giveDof(idof)->giveEquationNumber(dn);
                    if ( calm_Control == calml_hpc ) {
                        weights.at(count) = calm_HPCDmanWeightSrcArray.at(i);
                    }
                }

 #else
                indirectMap.at(++count) = domain->giveNode(j)->giveDof(idof)->giveEquationNumber(dn);
                if ( calm_Control == calml_hpc ) {
                    weights.at(count) = calm_HPCDmanWeightSrcArray.at(i);
                }

 #endif

                continue;
            }
        }
    }

 #ifndef __PARALLEL_MODE
    if ( count != size ) {
        OOFEM_WARNING("CylindricalALM  :: convertHPCMap: some dofmans/Dofs in HPCarray not recognized");
    }

 #endif

    calm_HPCIndirectDofMask.resize(count);
    if ( calm_Control == calml_hpc ) {
        calm_HPCWeights.resize(count);
    }

    for ( i = 1; i <= count; i++ ) {
        calm_HPCIndirectDofMask.at(i) = indirectMap.at(i);
        calm_HPCWeights.at(i) = weights.at(i);
    }

#else
    size = calm_HPCDmanDofSrcArray.giveSize() / 2;
    calm_HPCIndirectDofMask.resize(size);
    if ( calm_Control == calml_hpc ) {
        calm_HPCWeights.resize(size);
    }

    for ( i = 1; i <= size; i++ ) {
        inode = calm_HPCDmanDofSrcArray.at(2 * i - 1);
        idof  = calm_HPCDmanDofSrcArray.at(2 * i);
        calm_HPCIndirectDofMask.at(i) = domain->giveNode(inode)->giveDof(idof)->giveEquationNumber(dn);

        if ( calm_Control == calml_hpc ) {
            calm_HPCWeights.at(i) = calm_HPCDmanWeightSrcArray.at(i);
        }
    }

#endif
}


contextIOResultType
CylindricalALM :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    // write current deltaL
    if ( !stream->write(& deltaL, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}


contextIOResultType
CylindricalALM :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    // read last deltaL
    if ( !stream->read(& deltaL, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}


SparseLinearSystemNM *
CylindricalALM :: giveLinearSolver()
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


int
CylindricalALM :: computeDeltaLambda(double &deltaLambda, const FloatArray &dX, const FloatArray &deltaXt,
                                     const FloatArray &deltaX_, const FloatArray &R, double RR, double eta,
                                     double deltaL, double DeltaLambda0, int neq)
{
    double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4, a5;
    double _RR, _rr, _a2, _a3, _pr;
    double lam1, lam2, cos1, cos2;
    int i, ind, HPsize = 0;
    FloatArray help;
#ifndef __PARALLEL_MODE
    double XX;
#endif

#ifdef __PARALLEL_MODE
 #ifdef __PETSC_MODULE
    //PetscNatural2LocalOrdering *n2l = parallel_context->giveN2Lmap();
 #endif
#endif
    //
    // B.3.
    //
    if ( ( calm_Control == calm_hpc_off ) || ( calm_Control == calm_hpc_on ) ) {
        if ( calm_Control == calm_hpc_off ) {
            // this two lines are necessary if NRM is used
            // (for MNRM they can be computed at startup A1).
#ifdef __PARALLEL_MODE
            double prod[6];
            prod[0] = parallel_context->localNorm(deltaXt); prod[0] *= prod[0];
            prod[1] = parallel_context->localDotProduct(dX, deltaXt);
            prod[2] = parallel_context->localDotProduct(deltaX_, deltaXt);
            prod[3] = parallel_context->localNorm(dX); prod[3] *= prod[3];
            prod[4] = parallel_context->localDotProduct(dX, deltaX_);
            prod[5] = parallel_context->localNorm(deltaX_); prod[5] *= prod[5];

            a1 = eta * eta * prod[0] + Psi * Psi * RR;
            a2 = RR * Psi * Psi * DeltaLambda0 * 2.0;
            a2 += 2.0 * eta * prod [ 1 ];
            a2 += 2.0 * eta * eta * prod [ 2 ];
            a3 = prod [ 3 ];
            a3 += 2.0 * eta * prod [ 4 ];
            a3 += eta * eta * prod [ 5 ];
            a3 += DeltaLambda0 * DeltaLambda0 * RR * Psi * Psi - deltaL * deltaL;
#else
            XX = deltaXt.computeSquaredNorm();
            a1 = eta * eta * XX + Psi * Psi * RR;
            a2 = RR * Psi * Psi * DeltaLambda0 * 2.0;
            a2 += 2.0 * eta * dX.dotProduct(deltaXt);
            a2 += 2.0 * eta * eta * deltaX_.dotProduct(deltaXt);
            a3 = dX.computeSquaredNorm();
            a3 += 2.0 *eta * dX.dotProduct(deltaX_);
            a3 += eta * eta * deltaX_.computeSquaredNorm();
            a3 += DeltaLambda0 * DeltaLambda0 * RR * Psi * Psi - deltaL * deltaL;
#endif
        } else if ( calm_Control == calm_hpc_on ) {
            HPsize = calm_HPCIndirectDofMask.giveSize();
            _rr = 0.;
            _RR = 0.;
            _a2 = _a3 = 0.;
            for ( i = 1; i <= HPsize; i++ ) {
                if ( ( ind = calm_HPCIndirectDofMask.at(i) ) != 0 ) {
                    _rr += deltaXt.at(ind) * deltaXt.at(ind);
                    _RR += R.at(ind) * R.at(ind);
                    _pr   = ( dX.at(ind) + eta * deltaX_.at(ind) );
                    _a2 += eta * deltaXt.at(ind)  * _pr;
                    _a3 += _pr * _pr;
                }
            }

#ifdef __PARALLEL_MODE
            FloatArray my_(4), col_;
            my_(0) = _rr;
            my_(1) = _RR;
            my_(2) = _a2;
            my_(3) = _a3;
            parallel_context->accumulate(my_, col_);
            a1 = eta * eta * col_(0) + Psi * Psi * col_(1);
            a2 = _RR * Psi * Psi * DeltaLambda0 * 2.0;
            a2 += 2.0 * col_(2);
            a3 = col_(3) - deltaL * deltaL + DeltaLambda0 * DeltaLambda0 * col_(1) * Psi * Psi;
#else
            a1 = eta * eta * _rr + Psi * Psi * _RR;
            a2 = _RR * Psi * Psi * DeltaLambda0 * 2.0;
            a2 += 2.0 * _a2;
            a3 = _a3 - deltaL * deltaL + DeltaLambda0 * DeltaLambda0 * _RR * Psi * Psi;
#endif
            //printf ("\na1=%e, a2=%e, a3=%e", a1,a2,a3);
        }

        // solution of quadratic eqn.
        double discr = a2 * a2 - 4.0 * a1 * a3;
        if ( discr < 0.0 ) {
            _error("computeDeltaLambda: discriminant is negative, solution failed");
        }

        discr = sqrt(discr);
        lam1 = ( -a2 + discr ) / 2. / a1;
        lam2 = ( -a2 - discr ) / 2. / a1;

        // select better lam (according to angle between deltar0 and deltar1(2).
        //
        // up to there rewritten
        //
        if ( calm_Control == calm_hpc_off ) {
#ifdef __PARALLEL_MODE
            double tmp = parallel_context->localNorm(dX); tmp *= tmp;
            a4 = eta*parallel_context->localDotProduct(dX, deltaX_) + tmp;
            a5 = eta*parallel_context->localDotProduct(dX, deltaXt);
#else
            a4 = eta * dX.dotProduct(deltaX_) + dX.computeSquaredNorm();
            a5 = eta * dX.dotProduct(deltaXt);
#endif
        } else {
            a4 = 0.;
            a5 = 0.;
            for ( i = 1; i <= HPsize; i++ ) {
                if ( ( ind = calm_HPCIndirectDofMask.at(i) ) != 0 ) {
                    a4 += eta * dX.at(ind) * deltaX_.at(ind) + dX.at(ind) * dX.at(ind);
                    a5 += eta * dX.at(ind) * deltaXt.at(ind);
                }
            }

#ifdef __PARALLEL_MODE
            FloatArray mya(2), cola;
            mya(0) = a4;
            mya(1) = a5;
            parallel_context->accumulate(mya, cola);
            a4 = cola(0);
            a5 = cola(1);
#endif
        }

        cos1 = ( a4 + a5 * lam1 ) / deltaL / deltaL;
        cos2 = ( a4 + a5 * lam2 ) / deltaL / deltaL;
        if ( cos1 > cos2 ) {
            deltaLambda = lam1;
        } else {
            deltaLambda = lam2;
        }

        //printf ("eta=%e, lam1=%e, lam2=%e", eta, lam1, lam2);
    } else if ( calm_Control == calml_hpc ) {
        // linearized control
        double nom = 0., denom = 0.;

        HPsize = calm_HPCIndirectDofMask.giveSize();
        for ( i = 1; i <= HPsize; i++ ) {
            if ( ( ind = calm_HPCIndirectDofMask.at(i) ) != 0 ) {
                nom += ( dX.at(ind) + eta * deltaX_.at(ind) ) * calm_HPCWeights.at(i);
                denom += eta * deltaXt.at(ind) * calm_HPCWeights.at(i);
            }
        }

#ifdef __PARALLEL_MODE
        FloatArray myv(2), colv(2);
        myv(0) = nom;
        myv(1) = denom;
        parallel_context->accumulate(myv, colv);
        nom = colv(0);
        denom = colv(1);
#endif
        if ( fabs(denom) < calm_SMALL_NUM ) {
            _error("\ncalm: zero denominator in linearized control");
        }

        deltaLambda = ( deltaL - nom ) / denom;
    }

    return 0;
}


void
CylindricalALM :: search(int istep, FloatArray &prod, FloatArray &eta, double amp,
                         double maxetalim, double minetalim, int &ico)
{
    int i, ipos, ineg = 0;
    double etaneg = 1.0;
    double etamax = 0.0;


    // obtain ineg (number of previous line search iteration with negative ratio nearest to origin)
    // as well as max previous step length, etamax

    for ( i = 1; i <= istep; i++ ) {
        etamax = max( etamax, eta.at(i) );
        if ( prod.at(i) >= 0.0 ) {
            continue;
        }

        if ( eta.at(i) >= etaneg ) {
            continue;
        }

        etaneg = eta.at(i);
        ineg = i;
    }

    if ( ineg ) {
        // allow interpolation
        // first find ipos (position of previous s-l with positive ratio that is
        // closest to ineg (but with smaller s-l)
        ipos = 1;
        for ( i = 1; i <= istep; i++ ) {
            if ( prod.at(i) <= 0.0 ) {
                continue;
            }

            if ( eta.at(i) > eta.at(ineg) ) {
                continue;
            }

            if ( eta.at(i) < eta.at(ipos) ) {
                continue;
            }

            ipos = i;
        }

        // interpolate to get step-length
        double etaint = ( prod.at(ineg) * eta.at(ipos) - prod.at(ipos) * eta.at(ineg) ) / ( prod.at(ineg) - prod.at(ipos) );
        // alternatively get eta ensuring reasonable change
        double etaalt = eta.at(ipos) + 0.2 * ( eta.at(ineg) - eta.at(ipos) );
        etaint = max(etaint, etaalt);
        if ( etaint < minetalim ) {
            etaint = minetalim;
            if ( ico == 1 ) {
                ico = 2;
            } else {
                ico = 1;
            }
        }

        eta.at(istep + 1) = etaint;

    } else { // ineq == 0
        // allow extrapolation
        double etamaxstep = amp * etamax;
        // extrapolate between current and previous
        double etaextrap = ( prod.at(istep) * eta.at(istep - 1) - prod.at(istep - 1) * eta.at(istep) ) /
                           ( prod.at(istep) - prod.at(istep - 1) );
        eta.at(istep + 1) = etaextrap;
        // check if in limits
        if ( ( etaextrap <= 0.0 ) || ( etaextrap > etamaxstep ) ) {
            eta.at(istep + 1) = etamaxstep;
        }

        if ( ( eta.at(istep + 1) > maxetalim ) && ( ico == 1 ) ) {
            ico = 2;
            return;
        }

        if ( ( eta.at(istep + 1) > maxetalim ) ) {
            ico = 1;
            eta.at(istep + 1) = maxetalim;
        }
    }
}

void
CylindricalALM :: do_lineSearch(FloatArray &r, const FloatArray &rInitial, const FloatArray &deltaX_, const FloatArray &deltaXt,
                                const FloatArray &dXm1, FloatArray &dX, FloatArray &ddX,
                                const FloatArray &R, const FloatArray *R0, const FloatArray &F,
                                double &DeltaLambda, double &DeltaLambdam1, double &deltaLambda,
                                double &Lambda, double &ReachedLambda, double RR, double &drProduct, TimeStep *tNow)
{
    //
    //  LINE SEARCH
    //

    int neq = r.giveSize();
    int ls_failed, dl_failed = 0;
    int _iter = 0;
    int ico;
    int ls_maxiter = 10;

    DeltaLambda = DeltaLambdam1 + deltaLambda;
    Lambda = ReachedLambda + DeltaLambda;
    double deltaLambdaForEta1 = deltaLambda;

    double d6, d7, d8, d9;

#ifdef __PARALLEL_MODE
 #ifdef __PETSC_MODULE
    //PetscNatural2LocalOrdering *n2l = parallel_context->giveN2Lmap();
 #endif
    d6 = parallel_context->localDotProduct(deltaX_, F);
    d7 = parallel_context->localDotProduct(deltaXt, F);
    d8 = -1.0 * parallel_context->localDotProduct(deltaX_, R);
    d9 = -1.0 * parallel_context->localDotProduct(deltaXt, R);
#else
    d6 = deltaX_.dotProduct(F);
    d7 = deltaXt.dotProduct(F);
    d8 = -1.0 * deltaX_.dotProduct(R);
    d9 = -1.0 * deltaXt.dotProduct(R);
#endif
    double e1, e2, d10 = 0.0, d11 = 0.0;
    double s0, si;
    double prevEta, currEta;

    FloatArray eta(ls_maxiter + 1), prod(ls_maxiter + 1);


    if ( R0 ) {
#ifdef __PARALLEL_MODE
        d10 = -1.0 * parallel_context->localDotProduct(deltaX_, *R0);
        d11 = -1.0 * parallel_context->localDotProduct(deltaXt, *R0);
#else
        d10 = -1.0 * deltaX_.dotProduct(*R0);
        d11 = -1.0 * deltaXt.dotProduct(*R0);
#endif
    }

    // prepare starting product ratios and step lengths
    prod.at(1) = 1.0;
    eta.at(1) = 0.0;
    currEta = eta.at(2) = 1.0;
    // following counter shows how many times the max or min step length has been reached
    ico = 0;

    //
    // begin line search loop
    //
    ls_failed = 1;
    for ( int ils = 2; ils <= ls_maxiter; ils++ ) {
        // update displacements
        drProduct = 0.0; // dotproduct of iterative displacement increment vector

        ddX.resize(0);
        ddX.add(eta.at(ils)*deltaLambda, deltaXt);
        ddX.add(eta.at(ils), deltaX_);
        dX = dXm1;
        dX.add(ddX);

#ifdef __PARALLEL_MODE
        drProduct = parallel_context->localNorm(ddX); drProduct *= drProduct;
#else
        drProduct = ddX.computeSquaredNorm();
#endif

        tNow->incrementStateCounter(); // update solution state counter
        // update internal forces according to new state
        engngModel->updateComponent(tNow, InternalRhs, domain);

#ifdef __PARALLEL_MODE
        e1 = parallel_context->localDotProduct(deltaX_, F);
        e2 = parallel_context->localDotProduct(deltaXt, F);
#else
        e1 = deltaX_.dotProduct(F);
        e2 = deltaXt.dotProduct(F);
#endif

        s0 = d6 + deltaLambda * d7 + Lambda * d8 + deltaLambda * Lambda * d9 + d10 + deltaLambda * d11;
        si = e1 + deltaLambda * e2 + Lambda * d8 + deltaLambda * Lambda * d9 + d10 + deltaLambda * d11;
        prod.at(ils) = si / s0;

        //printf ("\ns0=%e, si=%e, prod=%e", s0, si, prod.at(ils));
        if ( s0 >= 0.0 ) {
            //printf ("solve starting inner product uphill, val=%e",s0);
            ls_failed = 3;
            currEta = 1.0;
            break;
        }

        if ( fabs(si / s0) < ls_tolerance ) {
            ls_failed = 0;
            currEta = eta.at(ils);
            break;
        }

        _iter = 0;

        currEta = eta.at(ils);
        //printf ("\n_ite=%d, deltaLambda=%e, eta=%e", _iter, deltaLambda, currEta);
        do { // solve simultaneously the equations for eta and lambda
            _iter++;
            prevEta = currEta;
            s0 = d6 + deltaLambda * d7 + Lambda * d8 + deltaLambda * Lambda * d9 + d10 + deltaLambda * d11;
            si = e1 + deltaLambda * e2 + Lambda * d8 + deltaLambda * Lambda * d9 + d10 + deltaLambda * d11;
            prod.at(ils) = si / s0;

            // call line-search routine to get new estimate of eta.at(ils+1)
            this->search(ils, prod, eta, amplifFactor, maxEta, minEta, ico);
            if ( ico == 2 ) {
                ls_failed = 2;
                break; // exit the loop
            }

            currEta = eta.at(ils + 1);
            // solve for deltaLambda
            dl_failed = this->computeDeltaLambda(deltaLambda, dXm1, deltaXt, deltaX_, R, RR, currEta, deltaL, DeltaLambdam1, neq);
            if ( dl_failed ) {
                eta.at(ils + 1) = 1.0;
                deltaLambda = deltaLambdaForEta1;
                break;
            }

            DeltaLambda = DeltaLambdam1 + deltaLambda;
            Lambda = ReachedLambda + DeltaLambda;
            //printf ("\n_ite=%d, deltaLambda=%e, eta=%e", _iter, deltaLambda, currEta);
        } while ( ( _iter < 10 ) && ( fabs( ( currEta - prevEta ) / prevEta ) > 0.01 ) );

        if ( ( ls_failed > 1 ) || dl_failed ) {
            break;
        }

        //printf ("\ncalm ls...");
        //printf ("eta = %e, err=%d, ", currEta,ls_failed);
        //printf ("dLambda=%e, lerr=%d ", deltaLambda, dl_failed);
    } // end of line search loop

    if ( ls_failed || dl_failed ) {
        // last resort
        deltaLambda = deltaLambdaForEta1;
        drProduct = 0.0; // dotproduct of iterative displacement increment vector

        ddX.resize(0);
        ddX.add(deltaLambda, deltaXt);
        ddX.add(deltaX_);
        dX = dXm1;
        dX.add(ddX);

#ifdef __PARALLEL_MODE
        drProduct = parallel_context->localNorm(ddX); drProduct *= drProduct;
#endif
        tNow->incrementStateCounter(); // update solution state counter
        engngModel->updateComponent(tNow, InternalRhs, domain);
        DeltaLambda = DeltaLambdam1 + deltaLambda;
        Lambda = ReachedLambda + DeltaLambda;

        OOFEM_LOG_INFO("CALMLS:       Line search - err_id=%d, eta=%e, dlambda=%e\n", ls_failed, 1.0, deltaLambda);
    } else {
        OOFEM_LOG_INFO("CALMLS:       Line search - err_id=%d, eta=%e, dlambda=%e\n", ls_failed, currEta, deltaLambda);
    }
}
} // end namespace oofem
