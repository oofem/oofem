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

#include "linesearch.h"
#include "timestep.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "convergedreason.h"
#include "engngm.h"

namespace oofem {
LineSearchNM :: LineSearchNM(Domain *d, EngngModel *m) :
    NumericalMethod(d, m)
{
    max_iter = 10;
    ls_tolerance = 0.80;
    amplifFactor = 2.5;
    maxEta = 4.0;
    minEta = 0.2;
}

ConvergedReason
LineSearchNM :: solve(FloatArray &r, FloatArray &dr, FloatArray &F, FloatArray &R, FloatArray *R0,
                      IntArray &eqnmask, double lambda, double &etaValue, LS_status &status, TimeStep *tStep)
{
    int ico, ils, neq = r.giveSize();
    double s0;

    FloatArray g(neq), rb(neq);
    // Compute inner product at start and stop if positive
    g = R;
    g.times(lambda);
    if ( R0 ) {
        g.add(* R0);
    }

    g.subtract(F);

    for ( auto eq : eqnmask ) {
        g.at( eq ) = 0.0;
    }

    s0 = ( -1.0 ) * g.dotProduct(dr);
    if ( s0 >= 0.0 ) {
        //printf ("\nLineSearchNM::solve starting inner product uphill, val=%e",s0);
        OOFEM_LOG_DEBUG("LS: product uphill, eta=%e\n", 1.0);
        r.add(dr);
        tStep->incrementStateCounter();        // update solution state counter
        engngModel->updateComponent(tStep, InternalRhs, domain);
        etaValue = 1.0;
        status = ls_ok;
        return CR_CONVERGED;
    }

    // keep original total displacement r
    rb = r;

    eta.resize(this->max_iter + 1);
    prod.resize(this->max_iter + 1);
    // prepare starting product ratios and step lengths
    prod.at(1) = 1.0;
    eta.at(1) = 0.0;
    eta.at(2) = 1.0;
    // following counter shows how many times the max or min step length has been reached
    ico = 0;

    // begin line search loop
    for ( ils = 2; ils <= this->max_iter; ils++ ) {
        // update displacements
        r = rb;
        r.add(this->eta.at(ils), dr);

        tStep->incrementStateCounter();        // update solution state counter
        // update internal forces according to new state
        engngModel->updateComponent(tStep, InternalRhs, domain);
        // compute out-of balance forces g in new state
        g = R;
        g.times(lambda);
        if ( R0 ) {
            g.add(* R0);
        }

        g.subtract(F);

        for ( auto eq : eqnmask ) {
            g.at( eq ) = 0.0;
        }

        // compute current inner-product ratio
        double si = ( -1.0 ) * g.dotProduct(dr) / s0;
        prod.at(ils) = si;

        // check if line-search tolerance is satisfied
        if ( fabs(si) < ls_tolerance ) {
            dr.times( this->eta.at(ils) );
            //printf ("\nLineSearchNM::solve tolerance satisfied for eta=%e, ils=%d", eta.at(ils),ils);
            OOFEM_LOG_DEBUG( "LS: ils=%d, eta=%e\n", ils, eta.at(ils) );

            etaValue = eta.at(ils);
            status = ls_ok;
            return CR_CONVERGED;
        }

        // call line-search routine to get new estimate of eta.at(ils)
        this->search(ils, prod, eta, this->amplifFactor, this->maxEta, this->minEta, ico);
        if ( ico == 2 ) {
            break; // exit the loop
        }
    } // end line search loop

    // exceeded no of ls iterations of ls failed
    //if (ico == 2) printf("\nLineSearchNM::solve max or min step length has been reached two times");
    //else printf("\nLineSearchNM::solve reached max number of ls searches");
    OOFEM_LOG_DEBUG( "LS: ils=%d, ico=%d, eta=%e\n", ils, ico, eta.at(ils) );
    /* update F before */
    r = rb;
    r.add(dr);

    tStep->incrementStateCounter();           // update solution state counter
    engngModel->updateComponent(tStep, InternalRhs, domain);
    etaValue = 1.0;
    status = ls_failed;
    return CR_DIVERGED_ITS;
}


void
LineSearchNM :: search(int istep, FloatArray &prod, FloatArray &eta, double amp,
                       double maxetalim, double minetalim, int &ico)
{
    int ineg = 0;
    double etaneg = 1.0;
    double etamax = 0.0;


    // obtain ineg (number of previous line search iteration with negative ratio nearest to origin)
    // as well as max previous step length, etamax

    for ( int i = 1; i <= istep; i++ ) {
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
        int ipos = 1;
        for ( int i = 1; i <= istep; i++ ) {
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
}

void
LineSearchNM :: initializeFrom(InputRecord &ir)
{
    /* default values set in constructor
     * ls_tolerance = 0.80;
     * amplifFactor = 2.5;
     * maxEta = 4.0;
     */
    IR_GIVE_OPTIONAL_FIELD(ir, ls_tolerance, _IFT_LineSearchNM_lsearchtol);
    if ( ls_tolerance < 0.6 ) {
        ls_tolerance = 0.6;
    }

    if ( ls_tolerance > 0.95 ) {
        ls_tolerance = 0.95;
    }

    IR_GIVE_OPTIONAL_FIELD(ir, amplifFactor, _IFT_LineSearchNM_lsearchamp);
    if ( amplifFactor < 1.0 ) {
        amplifFactor = 1.0;
    }

    if ( amplifFactor > 10.0 ) {
        amplifFactor = 10.0;
    }

    IR_GIVE_OPTIONAL_FIELD(ir, maxEta, _IFT_LineSearchNM_lsearchmaxeta);
    if ( maxEta < 1.5 ) {
        maxEta = 1.5;
    }

    if ( maxEta > 15.0 ) {
        maxEta = 15.0;
    }

    //printf ("\nLineSearchNM::initializeFrom: tol=%e, ampl=%e, maxEta=%e\n",
    //    ls_tolerance, amplifFactor,maxEta);
}
} // end namespace oofem
