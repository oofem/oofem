/* $Header: /home/cvs/bp/oofem/oofemlib/src/calmls.C,v 1.2.4.2 2004/04/13 11:28:15 bp Exp $ */
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
//
// file calm.cc
//

#include "calmls.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#include <math.h>
#endif

#include "verbose.h"
#include "ldltfact.h"
#include "imlsolver.h"
#include "timestep.h"
#include "flotmtrx.h"
#include "datastream.h"
#include "mathfem.h"
// includes for HPC - not very clean (NumMethod knows what is "node" and "dof")
#include "node.h"
#include "dof.h"
#include "usrdefsub.h"
#include "contextioerr.h"

#ifdef __PARALLEL_MODE
#ifndef __MAKEDEPEND
#include "mpi.h"
#endif
#endif


#define CALM_RESET_STEP_REDUCE 0.25
#define CALM_MAX_RESTARTS 4
#define CALM_TANGENT_STIFF_TRESHOLD 0.1
#define CALM_DEFAULT_NRM_TICKS 2
#define CALM_MAX_REL_ERROR_BOUND 1.e10

CylindricalALM :: CylindricalALM(int i, Domain *d, EngngModel *m, EquationID ut) :
    SparseNonLinearSystemNM(i, d, m, ut), calm_HPCWeights(), calm_HPCIndirectDofMask(), calm_HPCDmanDofSrcArray()
{
    //
    // constructor
    //
    //k      = NULL ;
    //Rt     = NULL ;
    //R0     = NULL ;
    //DeltaR = NULL ;
    //deltaRt= NULL ;
    //F      = NULL ;

    nsmax  = 60;       // default maximum number of sweeps allowed
    numberOfRequiredIterations = 3;
    //rtol   = 10.E-3  ;   // convergence tolerance
    //Psi    = 0.1;       // displacement control on
    Psi    = 1.0;     // load control on
    solved = 0;
    calm_NR_Mode = calm_NR_OldMode = calm_modifiedNRM;
    calm_NR_ModeTick = -1; // do not swith to calm_NR_OldMode
    calm_MANRMSteps = 0;

    //Bergan_k0 = 0.;    // value used for computing Bergan's parametr
    // of current stiffness.

    deltaL    = -1.0;
    // TangenStiffnessTreshold = 0.10;

    // Variables for Hyper Plane Controll
    calm_Controll = calm_hpc_off; // HPControll is not default
    linSolver = NULL;
    // linesearch default off
    lsFlag = 0;

    ls_tolerance = 0.40;
    amplifFactor = 2.5;
    maxEta = 8.5;
    minEta = 0.2;
}

CylindricalALM ::  ~CylindricalALM() {
    //
    // destructor
    //

    //if (deltaRt) delete deltaRt;

    delete linSolver;
}


NM_Status
CylindricalALM :: solve(SparseMtrx *k, FloatArray *Ri, FloatArray *R0,
                        FloatArray *Rr, FloatArray *r, FloatArray *DeltaR, FloatArray *F,
                        double &ReachedLambda, double rtol, referenceLoadInputModeType rlm,
                        int &nite, TimeStep *tNow)
{
    FloatArray rhs, *R, deltaRt, deltaR_, DeltaRm1;
    FloatArray rInitial;
    //FloatArray F;
    double Bergan_k0 = 1.0;
    double rr, RR, RR0, rR, p = 0.0, bk, forceErr, dispErr, drr;
    double deltaLambda, Lambda, eta, DeltaLambdam1, DeltaLambda = 0.0;
    double __rIterIncr, __rIncr, drProduct = 0.0;
    int neq = Ri->giveSize();
    int irest = 0;
    int HPsize, i, ind;
    double _RR, _rr;
    NM_Status status;
#ifndef __PARALLEL_MODE
    double *p1;
#endif

    OOFEM_LOG_INFO("CALM: Initial step length: %-15e\n", deltaL);
    OOFEM_LOG_INFO("Iteration  LoadLevel       ForceError      DisplError\n__________________________________________________________\n");
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

    rInitial = * r;
    // compute proportional load vector R
    R = Ri;
    /*
     * if (R0 && (refLoadInputMode == rlm_total)) {
     * R = Rt->GiveCopy();
     * R -> substract (R0);
     * } else R = Rt;
     */

    status = NM_None;
    this->giveLinearSolver();
    // NumericalMethod* linSolver =
    //   new LDLTFactorization (this->giveNumber()+1,this->giveDomain(), engngModel);

    // create HPC Map if needed
    if ( calm_hpc_init ) {
        this->convertHPCMap();
        calm_hpc_init = 0;
    }

#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
    // HUHU hard wired domain no 1
    PetscNatural2LocalOrdering *n2l = engngModel->givePetscContext(1, ut)->giveN2Lmap();
    if ( R0 ) {
        double myRR0 = 0.0;
        for ( i = 1; i <= neq; i++ ) {
            if ( n2l->giveNewEq(i) ) {
                myRR0 += R0->at(i) * R0->at(i);
            }
        }

        MPI_Allreduce(& myRR0, & RR0, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    } else {
        RR0 = 0.0;
    }

#endif
#else
    if ( R0 ) {
        RR0 = dotProduct(R0->givePointer(), R0->givePointer(), neq);
    } else {
        RR0 = 0.0;
    }

#endif

    //
    // A  initial step (predictor)
    //
    //
    //
    // A.1. calculation of (0)VarRt
    //
    // A1:
restart:
    //DeltaR -> zero();

    /*
     * linSolver -> setSparseMtrxAsComponent (LinearEquationLhs,k);
     *
     * // if (tNow ->giveNumber() == 1) {
     * deltaRt = new FloatArray (R->giveSize());
     * linSolver -> setFloatArrayAsComponent (LinearEquationRhs,R);
     * linSolver -> setFloatArrayAsComponent (LinearEquationSolution,deltaRt);
     * linSolver -> solveYourselfAt (tNow);
     * linSolver -> updateYourselfExceptLhs ();
     */
    deltaRt.resize( R->giveSize() );
#ifdef __PARALLEL_MODE
#ifdef __VERBOSE_PARALLEL
    //VERBOSEPARALLEL_PRINT("calm:: Solving linear system","", engngModel->giveRank());
#endif
#endif
    linSolver->solve(k, R, & deltaRt);

    // }

    //
    // A.2.   We assume positive-definitive (0)Kt (tangent stiffness mtrx).
    //
    // A2:
#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
    double myRR = 0.0;
    for ( i = 1; i <= neq; i++ ) {
        if ( n2l->giveNewEq(i) ) {
            myRR += R->at(i) * R->at(i);
        }
    }

    MPI_Allreduce(& myRR, & RR, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
#else
    RR = dotProduct(R->givePointer(), R->givePointer(), neq);
#endif

    if ( calm_Controll == calm_hpc_off ) {
#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
        double myrr = 0.0;
        for ( i = 1; i <= neq; i++ ) {
            if ( n2l->giveNewEq(i) ) {
                myrr += deltaRt.at(i) * deltaRt.at(i);
            }
        }

        MPI_Allreduce(& myrr, & rr, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
#else
        rr = dotProduct(deltaRt.givePointer(), deltaRt.givePointer(), neq);
#endif
        p = sqrt(rr + Psi * Psi * RR);
    } else if ( calm_Controll == calm_hpc_on ) {
        HPsize = calm_HPCIndirectDofMask.giveSize();
        _rr = 0;
        _RR = 0;
        for ( i = 1; i <= HPsize; i++ ) {
            if ( ( ind = calm_HPCIndirectDofMask.at(i) ) != 0 ) {
                _rr += deltaRt.at(ind) * deltaRt.at(ind);
                _RR += R->at(ind) * R->at(ind);
            }
        }

#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
        double my_rrRR [ 2 ], colected_rrRR [ 2 ];
        my_rrRR [ 0 ] = _rr;
        my_rrRR [ 1 ] = _RR;
        MPI_Allreduce(my_rrRR, colected_rrRR, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        _rr = colected_rrRR [ 0 ];
        _RR = colected_rrRR [ 1 ];
#endif
#endif

        p = sqrt(_rr + Psi * Psi * _RR);
    } else if ( calm_Controll == calml_hpc ) {
        HPsize = calm_HPCIndirectDofMask.giveSize();
        p = 0.;
        for ( i = 1; i <= HPsize; i++ ) {
            if ( ( ind = calm_HPCIndirectDofMask.at(i) ) != 0 ) {
                p += deltaRt.at(ind) * calm_HPCWeights.at(i);
            }
        }

#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
        double my_p = p;
        MPI_Allreduce(& my_p, & p, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
#endif
    }



#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
    double myrR = 0.0;
    for ( i = 1; i <= neq; i++ ) {
        if ( n2l->giveNewEq(i) ) {
            myrR += deltaRt.at(i) * R->at(i);
        }
    }

    MPI_Allreduce(& myrR, & rR, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
#else
    rR = dotProduct(deltaRt.givePointer(), R->givePointer(), neq);
#endif
    /* rR is unscaled Bergan's param of current stiffness rR = deltaRt^T k deltaRt
     * this is used to test vhether k has negatine or positive slope */

    Lambda = ReachedLambda;
    DeltaLambda = deltaLambda = sgn(rR) * deltaL / p;
    // prevLambda = Lambda ;
    Lambda += DeltaLambda;                 // *
    //
    // A.3.
    //

    // if (stepN->giveNumber() == 1)
    // {
    //  // first load step
    //  ls -> giveRhs()->times(deltaLambda) ;
    // } else
    // {
    //  F = this -> giveInternalForces(deltar,stepN);
    //  ls -> giveRhs()->add (F->GiveCopy()->negated()) ;
    // }
    rhs = * R;
    rhs.times(DeltaLambda);
    if ( R0 ) {
        rhs.add(R0);
    }

    // DeltaR -> zero();
    /*
     * linSolver -> setFloatArrayAsComponent (LinearEquationRhs,&rhs);
     * linSolver -> setFloatArrayAsComponent (LinearEquationSolution,DeltaR);
     * linSolver -> solveYourselfAt (tNow);
     */
    linSolver->solve(k, & rhs, DeltaR);
    r->add(DeltaR);
    //linSolver -> updateYourselfExceptLhs ();

    //delete rhs; rhs = NULL;

    nite = 0;

    // update solution state counter
    tNow->incrementStateCounter();
    engngModel->updateComponent(tNow, InternalRhs, domain);
    //((NonLinearStatic *)engngModel) -> giveInternalForces(F, *DeltaR, tNow);

    do {
        nite++;
        DeltaRm1 = * DeltaR;
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
            // we assemble new tangent stiffness and compute new deltaRt
            // internal state of elements is updated by previous calling
            // of F = engngModel->GiveInternalForces(DeltaR,tNow);
            //
            engngModel->updateComponent(tNow, NonLinearLhs, domain);
            //
            // compute deltaRt for i-th iteration
            //
            /*
             * linSolver -> setSparseMtrxAsComponent (LinearEquationLhs,k);
             *
             * linSolver -> setFloatArrayAsComponent (LinearEquationRhs,R);
             * linSolver -> setFloatArrayAsComponent (LinearEquationSolution,deltaRt);
             * linSolver -> solveYourselfAt (tNow);
             * linSolver -> updateYourselfExceptLhs ();
             */
            linSolver->solve(k, R, & deltaRt);
        }

        // B.2.
        //

        rhs =  * R;
        rhs.times(Lambda);
        if ( R0 ) {
            rhs.add(R0);
        }

        rhs.substract(F);
        deltaR_.resize(neq);
        /*
         * linSolver -> setFloatArrayAsComponent (LinearEquationRhs,&rhs);
         * linSolver -> setFloatArrayAsComponent (LinearEquationSolution,deltaR_);
         * linSolver -> solveYourselfAt (tNow);
         * linSolver -> updateYourselfExceptLhs ();
         */
        linSolver->solve(k, & rhs, & deltaR_);
        eta = 1.0;
        //
        // B.3.
        //
        if ( this->computeDeltaLambda(deltaLambda, * DeltaR, deltaRt, deltaR_, * R, RR, eta, deltaL, DeltaLambda, neq) ) {
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
                r->operator=(rInitial);
                // reset all changes fro previous equilibrium state
                DeltaR->zero();
                // restore initial stiffness
                engngModel->updateComponent(tNow, NonLinearLhs, domain);
                //delete F;
                //delete deltaR_;

                OOFEM_LOG_INFO("calm iteration Reset ...\n");

                calm_NR_OldMode  = calm_NR_Mode;
                calm_NR_Mode     = calm_fullNRM;
                calm_NR_ModeTick = CALM_DEFAULT_NRM_TICKS;
                goto restart;
            } else {
                status = NM_NoSuccess;
                _error("solve: calm - can't continue further");
            }
        }

        //
        // B.4.+ B.5.
        //

        //
        //  LINE SEARCH
        //
        if ( this->lsFlag && ( nite != 1 ) ) {
            int ls_failed, dl_failed = 0;
            int _iter = 0;
            int ico, ii;
            int ls_maxiter = 10;

            DeltaLambda = DeltaLambdam1 + deltaLambda;
            Lambda = ReachedLambda + DeltaLambda;
            double deltaLambdaForEta1 = deltaLambda;

            double d6, d7, d8, d9;

#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
            double myd [ 4 ] = {
                0., 0., 0., 0.
            }, cold [ 4 ];
            for ( i = 1; i <= neq; i++ ) {
                if ( n2l->giveNewEq(i) ) {
                    myd [ 0 ] += deltaR_.at(i) * F->at(i);
                    myd [ 1 ] += deltaRt.at(i) * F->at(i);
                    myd [ 2 ] += deltaR_.at(i) * R->at(i);
                    myd [ 3 ] += deltaRt.at(i) * R->at(i);
                }
            }

            MPI_Allreduce(myd, cold, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            d6 = cold [ 0 ];
            d7 = cold [ 1 ];
            d8 = cold [ 2 ];
            d9 = cold [ 3 ];
#endif
#else
            d6 = dotProduct(deltaR_, * F, neq);
            d7 = dotProduct(deltaRt, * F, neq);
            d8 = -1.0 * dotProduct(deltaR_, * R, neq);
            d9 = -1.0 * dotProduct(deltaRt, * R, neq);
#endif
            double e1, e2, d10 = 0.0, d11 = 0.0;
            double s0, si;
            double prevEta, currEta;

            FloatArray eta(ls_maxiter + 1), prod(ls_maxiter + 1);


            if ( R0 ) {
#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
                double myd [ 2 ] = {
                    0., 0.
                }, cold [ 2 ];
                for ( i = 1; i <= neq; i++ ) {
                    if ( n2l->giveNewEq(i) ) {
                        myd [ 0 ] += deltaR_.at(i) * R0->at(i);
                        myd [ 1 ] += deltaRt.at(i) * R0->at(i);
                    }
                }

                MPI_Allreduce(myd, cold, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                d10 = -1.0 * cold [ 0 ];
                d11 = -1.0 * cold [ 1 ];
#endif
#else
                d10 = -1.0 * dotProduct(deltaR_, * R0, neq);
                d11 = -1.0 * dotProduct(deltaRt, * R0, neq);
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
#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
                double my_drProduct = 0.0;
                for ( ii = 1; ii <= neq; ii++ ) {
                    __rIterIncr = eta.at(ils) * ( deltaLambda * deltaRt.at(ii) + deltaR_.at(ii) );
                    r->at(ii) = rInitial.at(ii) + DeltaRm1.at(ii) + __rIterIncr;
                    DeltaR->at(ii) = DeltaRm1.at(ii) + __rIterIncr;
                    if ( n2l->giveNewEq(ii) ) {
                        my_drProduct += __rIterIncr * __rIterIncr;
                    }
                }

                MPI_Allreduce(& my_drProduct, & drProduct, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
#else
                for ( ii = 1; ii <= neq; ii++ ) {
                    __rIterIncr = eta.at(ils) * ( deltaLambda * deltaRt.at(ii) + deltaR_.at(ii) );
                    r->at(ii) = rInitial.at(ii) + DeltaRm1.at(ii) + __rIterIncr;
                    DeltaR->at(ii) = DeltaRm1.at(ii) + __rIterIncr;
                    drProduct += __rIterIncr * __rIterIncr;
                    //r->at(ii) = rInitial.at(ii) + DeltaRm1.at(ii) + eta.at(ils)*(deltaLambda*deltaRt.at(ii) + deltaR_.at(ii));
                    //DeltaR->at(ii) = DeltaRm1.at(ii) + eta.at(ils)*(deltaLambda*deltaRt.at(ii) + deltaR_.at(ii));
                }

#endif

                tNow->incrementStateCounter();  // update solution state counter
                // update internal forces according to new state
                engngModel->updateComponent(tNow, InternalRhs, domain);

#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
                double mye [ 2 ] = {
                    0., 0.
                }, cole [ 2 ];
                for ( ii = 1; ii <= neq; ii++ ) {
                    if ( n2l->giveNewEq(ii) ) {
                        mye [ 0 ] += deltaR_.at(ii) * F->at(ii);
                        mye [ 1 ] += deltaRt.at(ii) * F->at(ii);
                    }
                }

                MPI_Allreduce(mye, cole, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                e1 = cole [ 0 ];
                e2 = cole [ 1 ];
#endif
#else
                e1 = dotProduct(deltaR_, * F, neq);
                e2 = dotProduct(deltaRt, * F, neq);
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
                    dl_failed = this->computeDeltaLambda(deltaLambda, DeltaRm1, deltaRt, deltaR_, * R, RR, currEta, deltaL, DeltaLambdam1, neq);
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

#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
                double my_drProduct = 0.0;
                for ( ii = 1; ii <= neq; ii++ ) {
                    __rIterIncr = 1.0 * ( deltaLambda * deltaRt.at(ii) + deltaR_.at(ii) );
                    r->at(ii) = rInitial.at(ii) + DeltaRm1.at(ii) + __rIterIncr;
                    DeltaR->at(ii) = DeltaRm1.at(ii) + __rIterIncr;
                    if ( n2l->giveNewEq(ii) ) {
                        my_drProduct += __rIterIncr * __rIterIncr;
                    }
                }

                MPI_Allreduce(& my_drProduct, & drProduct, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
#else
                for ( ii = 1; ii <= neq; ii++ ) {
                    __rIterIncr = 1.0 * ( deltaLambda * deltaRt.at(ii) + deltaR_.at(ii) );
                    r->at(ii) = rInitial.at(ii) + DeltaRm1.at(ii) + __rIterIncr;
                    DeltaR->at(ii) = DeltaRm1.at(ii) + __rIterIncr;
                    drProduct += __rIterIncr * __rIterIncr;
                }

#endif
                tNow->incrementStateCounter();  // update solution state counter
                engngModel->updateComponent(tNow, InternalRhs, domain);
                DeltaLambda = DeltaLambdam1 + deltaLambda;
                Lambda = ReachedLambda + DeltaLambda;

                //printf ("\ncalm fi...eta = %e, err=%d, ", 1.0,ls_failed);
                //printf ("dLambda=%e, lerr=%d", deltaLambda, dl_failed);
                OOFEM_LOG_INFO("LS: err_id=%d, eta=%e, dlambda=%e\n", ls_failed, 1.0, deltaLambda);
            } else {
                //printf ("\ncalm fi...eta = %e, err=%d, ", currEta,ls_failed);
                //printf ("dLambda=%e, lerr=%d", deltaLambda, dl_failed);
                OOFEM_LOG_INFO("LS: err_id=%d, eta=%e, dlambda=%e\n", ls_failed, currEta, deltaLambda);
            }
        } else { // no line search
            //
            // update solution vectors
            //
            drProduct = 0.0; // dotproduct of iterative displacement increment vector
#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
            double my_drProduct = 0.0;
            for ( i = 1; i <= neq; i++ ) {
                __rIterIncr = eta * ( deltaLambda * deltaRt.at(i) + deltaR_.at(i) );
                __rIncr = DeltaRm1.at(i) +  __rIterIncr;

                DeltaR->at(i) = __rIncr;
                r->at(i) = rInitial.at(i) + __rIncr;
                if ( n2l->giveNewEq(i) ) {
                    my_drProduct += __rIterIncr * __rIterIncr;
                }
            }

            MPI_Allreduce(& my_drProduct, & drProduct, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#endif
#else
            for ( i = 1; i <= neq; i++ ) {
                __rIterIncr = eta * ( deltaLambda * deltaRt.at(i) + deltaR_.at(i) );
                __rIncr = DeltaRm1.at(i) +  __rIterIncr;

                DeltaR->at(i) = __rIncr;
                r->at(i) = rInitial.at(i) + __rIncr;

                drProduct += __rIterIncr * __rIterIncr;
            }

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
        // convergency check
        //
        rhs =  * R;
        rhs.times(Lambda);
        if ( R0 ) {
            rhs.add(R0);
        }

        rhs.substract(F);

        //
        // compute forceError
        //
#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
        double myerr [ 2 ] = {
            0., 0.
        }, colerr [ 2 ];
        for ( i = 1; i <= neq; i++ ) {
            if ( n2l->giveNewEq(i) ) {
                myerr [ 0 ] += rhs.at(i) * rhs.at(i);
                myerr [ 1 ] += r->at(i) * r->at(i);
            }
        }

        MPI_Allreduce(myerr, colerr, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        forceErr = colerr [ 0 ];
        drr = colerr [ 1 ];
#endif
#else
        // err is relative error of unbalanced forces
        forceErr = dotProduct(rhs.givePointer(), rhs.givePointer(), neq);
        // err is relative displacement change
        drr = dotProduct(r->givePointer(), r->givePointer(), neq);
#endif
        // we compute a relative error norm
        if ( ( RR0 + RR * Lambda * Lambda ) < calm_SMALL_ERROR_NUM ) {
            sqrt(forceErr);
        } else {
            forceErr = sqrt( forceErr / ( RR0 + RR * Lambda * Lambda ) );
        }

        //
        // compute displacement error
        //
        if ( drr < calm_SMALL_ERROR_NUM ) {
            dispErr = drProduct;
        } else {
            dispErr = drProduct / drr;
            dispErr = sqrt(dispErr);
        }

        //delete rhs; rhs = NULL;
        //delete Delta; Delta = NULL;

        //
        // Restart if nite >= nsmax of if force or displacement error is bigger
        // than allowed limit (rtol * CALM_MAX_REL_ERROR_BOUND)
        //
        if ( ( nite >= nsmax ) ||
            ( fabs(forceErr) > rtol * CALM_MAX_REL_ERROR_BOUND ) ||
            ( fabs(dispErr)  > rtol * CALM_MAX_REL_ERROR_BOUND ) ) {
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
                r->operator=(rInitial);
                // reset all changes fro previous equilibrium state
                engngModel->initStepIncrements();
                DeltaR->zero();
                // restore initial stiffness
                engngModel->updateComponent(tNow, NonLinearLhs, domain);
                //delete F; F = NULL;

                OOFEM_LOG_INFO("calm iteration Reset ...\n");

                calm_NR_OldMode  = calm_NR_Mode;
                calm_NR_Mode     = calm_fullNRM;
                calm_NR_ModeTick = CALM_DEFAULT_NRM_TICKS;
                goto restart;
            } else {
                status = NM_NoSuccess;
                _warning2("CALM - convergence not reached after %d iterations", nsmax);
                // exit(1);
                break;
            }
        }

        OOFEM_LOG_INFO("%-10d %-15e %-15e %-15e\n", nite, Lambda, forceErr, dispErr);
    } while ( ( fabs(forceErr) > rtol ) || ( fabs(dispErr) > rtol ) );

    //delete F;
    //delete linSolver;
    //delete deltaRt; deltaRt = NULL;
    //
    // end of iteration
    //
    // ls ->letSolutionBe(deltar);
    // Lambda += DeltaLambda ;      // *
    //
    // update dofs,nodes,Elemms and print result
    //
#ifdef VERBOSE
    // printf ("\nCALM - step iteration finished") ;
#endif
    /* compute Bergan's parameter of current stiffness */

#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
    double myp [ 2 ] = {
        0., 0.
    }, colp [ 2 ];
    for ( i = 1; i <= neq; i++ ) {
        if ( n2l->giveNewEq(i) ) {
            myp [ 0 ] += R->at(i) * DeltaR->at(i);
            myp [ 1 ] += DeltaR->at(i) * DeltaR->at(i);
        }
    }

    MPI_Allreduce(myp, colp, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    bk = DeltaLambda * colp [ 0 ];
    bk = bk / colp [ 1 ];
#endif
#else
    p1 = DeltaR->givePointer();
    bk = DeltaLambda * dotProduct(R->givePointer(), p1, neq);
    bk = bk / dotProduct(p1, p1, neq);
#endif

    if ( tNow->giveNumber() == 1 ) {
        // compute Bergan_k0
        Bergan_k0 = bk;
    } else {
        bk = bk / Bergan_k0;
        //if (fabs(bk) < CALM_TANGENT_STIFF_TRESHOLD) status |= NM_KeepTangent;
    }

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

    OOFEM_LOG_INFO("CALM: Adjusted step length: %-15e\n", deltaL);

    // delete proportional load vector (if allocated)
    // if (R0 && (refLoadInputMode == rlm_total)) delete R;


    status = NM_Success;
    solved = 1;
    ReachedLambda = Lambda;

    return status;
}


IRResultType
CylindricalALM :: initializeFrom(InputRecord *ir)
//
//
//
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    double oldPsi =  Psi; // defaul from constructor
    double initialStepLength;
    int hpcMode;

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

    numberOfRequiredIterations = 3;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfRequiredIterations, IFT_CylindricalALM_reqiterations, "reqiterations"); // Macro
    if ( numberOfRequiredIterations < 3 ) {
        numberOfRequiredIterations = 3;
    }

    if ( numberOfRequiredIterations > 1000 ) {
        numberOfRequiredIterations = 1000;
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
    // calm_HPControll = hpc_init
    if ( calm_HPCDmanDofSrcArray.giveSize() != 0 ) {
        int i, nsize;
        if ( hpcMode == 1 ) {
            calm_Controll = calm_hpc_on;
        } else if ( hpcMode == 2 )                                                   {
            calm_Controll = calml_hpc;
        } else                                                                                                         {
            calm_Controll = calm_hpc_on; // default is to use hpc_on
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

#if defined(__PARALLEL_MODE) || defined(__ENABLE_COMPONENT_LABELS)
#ifdef __PETSC_MODULE
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
#if defined(__PARALLEL_MODE) && defined (__PETSC_MODULE)
                // HUHU hard wired domain no 1
                if ( engngModel->givePetscContext(1, ut)->giveN2Gmap()->isLocal( domain->giveNode(j) ) ) {
                    indirectMap.at(++count) = domain->giveNode(j)->giveDof(idof)->giveEquationNumber();
                    if ( calm_Controll == calml_hpc ) weights.at(count) = calm_HPCDmanWeightSrcArray.at(i);
                }
#else
                indirectMap.at(++count) = domain->giveNode(j)->giveDof(idof)->giveEquationNumber();
		if ( calm_Controll == calml_hpc ) weights.at(count) = calm_HPCDmanWeightSrcArray.at(i);
#endif

                continue;
            }
        }
    }

#ifndef __PARALLEL_MODE
    if (count != size) OOFEM_WARNING ("CylindricalALM  :: convertHPCMap: some dofmans/Dofs in HPCarray not recognized");
#endif

    calm_HPCIndirectDofMask.resize(count);
    if ( calm_Controll == calml_hpc ) calm_HPCWeights.resize(count);
    for ( i = 1; i <= count; i++ ) {
        calm_HPCIndirectDofMask.at(i) = indirectMap.at(i);
        calm_HPCWeights.at(i) = weights.at(i);
    }

#endif  //__PETSC_MODULE
#else
    size = calm_HPCDmanDofSrcArray.giveSize() / 2;
    calm_HPCIndirectDofMask.resize(size);
    if ( calm_Controll == calml_hpc ) calm_HPCWeights.resize(size);
    for ( i = 1; i <= size; i++ ) {
        inode = calm_HPCDmanDofSrcArray.at(2 * i - 1);
        idof  = calm_HPCDmanDofSrcArray.at(2 * i);
        calm_HPCIndirectDofMask.at(i) = domain->giveNode(inode)->giveDof(idof)->giveEquationNumber();
        if ( calm_Controll == calml_hpc ) calm_HPCWeights.at(i)=calm_HPCDmanWeightSrcArray.at(i);
        
    }
#endif
}


contextIOResultType
CylindricalALM :: saveContext(DataStream *stream, ContextMode mode, void *obj) {
    // write current deltaL
    if ( !stream->write(& deltaL, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

contextIOResultType
CylindricalALM :: restoreContext(DataStream *stream, ContextMode mode, void *obj) {
    // read last deltaL
    if ( !stream->read(& deltaL, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}

SparseLinearSystemNM *
CylindricalALM :: giveLinearSolver() {
    if ( linSolver ) {
        if ( linSolver->giveLinSystSolverType() == solverType ) {
            return linSolver;
        } else {
            delete linSolver;
        }
    }

    linSolver = :: CreateUsrDefSparseLinSolver(solverType, 1, domain, engngModel);
    if ( linSolver == NULL ) {
        _error("giveLinearSolver: linear solver creation failed");
    }

    return linSolver;
}



/*
 * CylindricalALM :: computeDeltaLambda RETURN Values:
 * 0 .. ok
 * 1 .. failed (restart)
 */
int
CylindricalALM :: computeDeltaLambda(double &deltaLambda, FloatArray &DeltaR, FloatArray &deltaRt,
                                     FloatArray &deltaR_, FloatArray &R, double RR, double eta,
                                     double deltaL, double DeltaLambda0, int neq)
{
    double a1 = 0.0, a2 = 0.0, a3 = 0.0, a4, a5;
    double _RR, _rr, _a2, _a3, _pr;
    double lam1, lam2, cos1, cos2;
    int i, ind, HPsize = 0;
    FloatArray help;
#ifndef __PARALLEL_MODE
    double rr, *p1, *p2;
#endif

#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
    // HUHU hard wired domain no 1
    PetscNatural2LocalOrdering *n2l = engngModel->givePetscContext(1, ut)->giveN2Lmap();
#endif
#endif
    //
    // B.3.
    //
    if ( ( calm_Controll == calm_hpc_off ) || ( calm_Controll == calm_hpc_on ) ) {
        if ( calm_Controll == calm_hpc_off ) {
            // this two lines are necesarry if NRM is used
            // (for MNRM they can be computed at startup A1).
#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
            double my_prod [ 6 ] = {
                0., 0., 0., 0., 0., 0.
            }, prod [ 6 ];

            for ( i = 1; i <= neq; i++ ) {
                if ( n2l->giveNewEq(i) ) {
                    my_prod [ 0 ] += deltaRt.at(i) * deltaRt.at(i); //rr
                    my_prod [ 1 ] += DeltaR.at(i) * deltaRt.at(i);
                    my_prod [ 2 ] += deltaR_.at(i) * deltaRt.at(i);
                    my_prod [ 3 ] += DeltaR.at(i) * DeltaR.at(i);
                    my_prod [ 4 ] += DeltaR.at(i) * deltaR_.at(i);
                    my_prod [ 5 ] += deltaR_.at(i) * deltaR_.at(i);
                }
            }

            MPI_Allreduce(my_prod, prod, 6, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // master will receive

            a1 = eta * eta * prod [ 0 ] + Psi * Psi * RR;
            a2 = RR * Psi * Psi * DeltaLambda0 * 2.0;
            a2 += 2.0 * eta * prod [ 1 ];
            a2 += 2.0 * eta * eta * prod [ 2 ];
            a3 = prod [ 3 ];
            a3 += 2.0 * eta * prod [ 4 ];
            a3 += eta * eta * prod [ 5 ];
            a3 += DeltaLambda0 * DeltaLambda0 * RR * Psi * Psi - deltaL * deltaL;
#endif
#else
            p1 = deltaRt.givePointer();
            rr = dotProduct(p1, p1, neq);
            a1 = eta * eta * rr + Psi * Psi * RR;
            a2 = RR * Psi * Psi * DeltaLambda0 * 2.0;
            a2 += 2.0 *eta *dotProduct(DeltaR.givePointer(), p1, neq);
            a2 += 2.0 *eta *eta *dotProduct(deltaR_.givePointer(), p1, neq);
            a3 = dotProduct(DeltaR.givePointer(), DeltaR.givePointer(), neq);
            a3 += 2.0 *eta *dotProduct(DeltaR.givePointer(), deltaR_.givePointer(), neq);
            a3 += eta * eta * dotProduct(deltaR_, deltaR_, neq);
            a3 += DeltaLambda0 * DeltaLambda0 * RR * Psi * Psi - deltaL * deltaL;
#endif
        } else if ( calm_Controll == calm_hpc_on ) {
            HPsize = calm_HPCIndirectDofMask.giveSize();
            _rr = 0.;
            _RR = 0.;
            _a2 = _a3 = 0.;
            for ( i = 1; i <= HPsize; i++ ) {
                if ( ( ind = calm_HPCIndirectDofMask.at(i) ) != 0 ) {
                    _rr += deltaRt.at(ind) * deltaRt.at(ind);
                    _RR += R.at(ind) * R.at(ind);
                    _pr   = ( DeltaR.at(ind) + eta * deltaR_.at(ind) );
                    _a2 += eta * deltaRt.at(ind)  * _pr;
                    _a3 += _pr * _pr;
                }
            }

#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
            double my_ [ 4 ] = {
                _rr, _RR, _a2, _a3
            }, col_ [ 4 ];
            MPI_Allreduce(my_, col_, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // master will receive
            a1 = eta * eta * col_ [ 0 ] + Psi * Psi * col_ [ 1 ];
            a2 = _RR * Psi * Psi * DeltaLambda0 * 2.0;
            a2 += 2.0 * col_ [ 2 ];
            a3 = col_ [ 3 ] - deltaL * deltaL + DeltaLambda0 * DeltaLambda0 * col_ [ 1 ] * Psi * Psi;
#endif
#else

            a1 = eta * eta * _rr + Psi * Psi * _RR;
            a2 = _RR * Psi * Psi * DeltaLambda0 * 2.0;
            a2 += 2.0 * _a2;
            a3 = _a3 - deltaL * deltaL + DeltaLambda0 * DeltaLambda0 * _RR * Psi * Psi;
#endif
            //printf ("\na1=%e, a2=%e, a3=%e", a1,a2,a3);
        }

        // solution of kvdratic eqn.
        double discr = a2 * a2 - 4.0 * a1 * a3;
        if ( discr < 0.0 ) {
            _error("computeDeltaLambda: discriminant is negative, solution failed");
        }

        discr = sqrt(discr);
        lam1 = ( -a2 + discr ) / 2. / a1;
        lam2 = ( -a2 - discr ) / 2. / a1;

        // select better lam (acording to angle between deltar0 and deltar1(2).
        //
        // up to there rewritten
        //
        if ( calm_Controll == calm_hpc_off ) {
#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
            double myp [ 3 ] = {
                0., 0., 0.
            }, colp [ 3 ];
            for ( i = 1; i <= neq; i++ ) {
                if ( n2l->giveNewEq(i) ) {
                    myp [ 0 ] += DeltaR.at(i) * deltaR_.at(i);
                    myp [ 1 ] += DeltaR.at(i) * DeltaR.at(i);
                    myp [ 2 ] += DeltaR.at(i) * deltaRt.at(i);
                }
            }

            MPI_Allreduce(myp, colp, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // master will receive
            a4 = eta * colp [ 0 ] + colp [ 1 ];
            a5 = eta * colp [ 2 ];
#endif
#else
            p1 = DeltaR.givePointer();
            p2 = deltaR_.givePointer();
            a4 = eta * dotProduct(p1, p2, neq) + dotProduct(p1, p1, neq);
            p2 = deltaRt.givePointer();
            a5 = eta * dotProduct(p1, p2, neq);
#endif
        } else {
            a4 = 0.;
            a5 = 0.;
            for ( i = 1; i <= HPsize; i++ ) {
                if ( ( ind = calm_HPCIndirectDofMask.at(i) ) != 0 ) {
                    a4 += eta * DeltaR.at(ind) * deltaR_.at(ind) + DeltaR.at(ind) * DeltaR.at(ind);
                    a5 += eta * DeltaR.at(ind) * deltaRt.at(ind);
                }
            }

#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
            double mya [ 2 ] = {
                a4, a5
            }, cola [ 2 ];
            MPI_Allreduce(mya, cola, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // master will receive
            a4 = cola [ 0 ];
            a5 = cola [ 1 ];
#endif
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
    } else if ( calm_Controll == calml_hpc ) {
        // linearized controll
        double nom = 0., denom = 0.;

        HPsize = calm_HPCIndirectDofMask.giveSize();
        for ( i = 1; i <= HPsize; i++ ) {
            if ( ( ind = calm_HPCIndirectDofMask.at(i) ) != 0 ) {
                nom += ( DeltaR.at(ind) + eta * deltaR_.at(ind) ) * calm_HPCWeights.at(i);
                denom += eta * deltaRt.at(ind) * calm_HPCWeights.at(i);
            }
        }

#ifdef __PARALLEL_MODE
#ifdef __PETSC_MODULE
        double myv [ 2 ] = {
            nom, denom
        }, colv [ 2 ];
        MPI_Allreduce(myv, colv, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); // master will receive
        nom = colv [ 0 ];
        denom = colv [ 1 ];
#endif
#endif
        if ( fabs(denom) < calm_SMALL_NUM ) {
            _error("\ncalm: zero denominator in linearized controll");
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
        // alternativelly get eta ensuring reasonable change
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
        return;
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

    return;
}
