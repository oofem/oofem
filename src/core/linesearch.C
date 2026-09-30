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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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
#include "convergenceexception.h"
#include "engngm.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace oofem {
LineSearchNM :: LineSearchNM(Domain *d, EngngModel *m) :
    NumericalMethod(d, m)
{
    // Robust Newton defaults for contact.  In particular, a contact-state
    // change should damp the Newton direction; extrapolating beyond the full
    // step is counterproductive when the active projection changes.
    max_iter = 5;
    ls_tolerance = 0.90;
    amplifFactor = 2.5;
    maxEta = 4.0;
    minEta = 0.01;
}

ConvergedReason
LineSearchNM :: solve(FloatArray &r, FloatArray &dr, FloatArray &F, FloatArray &R, FloatArray *R0,
                      IntArray &eqnmask, double lambda, double &etaValue, LS_status &status, TimeStep *tStep)
{
    const int neq = r.giveSize();
    const FloatArray rb(r);
    FloatArray freeCorrection(dr);
    FloatArray prescribedCorrection(neq);
    prescribedCorrection.zero();
    for ( auto eq : eqnmask ) {
        prescribedCorrection.at(eq) = freeCorrection.at(eq);
        freeCorrection.at(eq) = 0.0;
    }

    // OOFEM inserts prescribed increments into the first linear solution.
    // Split that out here rather than forming the initial line-search
    // residual from it directly: prescribed motion is exact, while only the
    // unconstrained Newton correction is scaled by the line search.
    FloatArray trialBase(rb);
    trialBase.add(prescribedCorrection);
    FloatArray g(neq);

    struct TrialMeasure {
        double directional = 0.0;
        double norm = 0.0;
    };

    auto computeResidualMeasure = [&]() {
        g = R;
        g.times(lambda);
        if ( R0 ) {
            g.add(* R0);
        }
        g.subtract(F);
        for ( auto eq : eqnmask ) {
            g.at(eq) = 0.0;
        }
        return TrialMeasure { g.dotProduct(freeCorrection), g.computeNorm() };
    };

    auto evaluateTrial = [&](double step) {
        r = trialBase;
        r.add(step, freeCorrection);
        tStep->incrementStateCounter();
        engngModel->initForNewIteration(domain, tStep, 0, r);
        try {
            engngModel->updateComponent(tStep, InternalRhs, domain);
            return computeResidualMeasure();
        } catch (const ConvergenceException &) {
            const double invalid = std::numeric_limits<double>::quiet_NaN();
            return TrialMeasure { invalid, invalid };
        } catch (const std::domain_error &) {
            const double invalid = std::numeric_limits<double>::quiet_NaN();
            return TrialMeasure { invalid, invalid };
        }
    };

    auto restoreUnperturbedState = [&]() {
        r = rb;
        tStep->incrementStateCounter();
        engngModel->initForNewIteration(domain, tStep, 0, r);
        engngModel->updateComponent(tStep, InternalRhs, domain);
    };

    // The scalar stationarity function is the residual projected onto the
    // fixed Newton direction — the energy criterion of Bonet and Wood.
    const bool hasPrescribedCorrection = prescribedCorrection.computeSquaredNorm() > 0.0;
    const TrialMeasure initial = hasPrescribedCorrection ? evaluateTrial(0.0) : computeResidualMeasure();
    if (!std::isfinite(initial.directional) || !std::isfinite(initial.norm)) {
        // The prescribed part of the increment alone inverted an element.
        // A line search is not allowed to weaken a Dirichlet condition; let
        // the engineering model reduce and retry the complete time step.
        restoreUnperturbedState();
        dr.zero();
        etaValue = 0.0;
        status = ls_failed;
        return CR_DIVERGED_ITS;
    }

    if (freeCorrection.computeSquaredNorm() < 1.e-60) {
        dr = prescribedCorrection;
        r = rb;
        etaValue = 1.0;
        status = ls_ok;
        return CR_CONVERGED;
    }

    const double rInitial = initial.directional;
    const double initialMagnitude = std::abs(rInitial);

    double step = 1.0;
    double bestStep = step;
    double bestDirectionalMagnitude = std::numeric_limits<double>::infinity();
    bool accepted = false;

    for (int iteration = 1; iteration <= max_iter; ++iteration) {
        const TrialMeasure trial = evaluateTrial(step);
        const double rTrial = trial.directional;
        const double magnitude = std::abs(rTrial);
        const bool validTrial = std::isfinite(magnitude) && std::isfinite(trial.norm);
        if (validTrial && magnitude < bestDirectionalMagnitude) {
            bestDirectionalMagnitude = magnitude;
            bestStep = step;
        }

        const double ratio = magnitude / std::max(initialMagnitude, 1.e-30);
        const double normRatio = trial.norm / std::max(initial.norm, 1.e-30);
        OOFEM_LOG_DEBUG("LS: iteration=%d, eta=%e, ratio=%e, norm_ratio=%e\n",
                        iteration, step, ratio, normRatio);
        if (validTrial && (initialMagnitude < 1.e-30 ? trial.norm <= initial.norm
                                                     : ratio <= ls_tolerance)) {
            accepted = true;
            break;
        }

        if (!validTrial) {
            // An inverted trial configuration is rejected and the correction
            // is halved. The minimum line-search bound is deliberately not
            // enforced for this path; several halvings may be needed to
            // regain J > 0.
            step *= 0.5;
            continue;
        }

        // Quadratic energy interpolation.  Unlike OOFEM's former secant
        // search, this never extrapolates beyond the full Newton step.
        double nextStep = 0.5 * step;
        if (initialMagnitude >= 1.e-30 && std::abs(rTrial) > 1.e-30) {
            const double a = rInitial / rTrial;
            const double A = 1.0 + a * (step - 1.0);
            const double B = a * step * step;
            const double discriminant = B * B - 4.0 * A * B;
            if (std::abs(A) > 1.e-30) {
                if (discriminant >= 0.0) {
                    nextStep = (B + std::sqrt(discriminant)) / (2.0 * A);
                    if (nextStep < 0.0) {
                        nextStep = (B - std::sqrt(discriminant)) / (2.0 * A);
                    }
                } else {
                    nextStep = 0.5 * B / A;
                }
            }
        }

        if (!std::isfinite(nextStep) || nextStep <= 0.0 || nextStep > 1.0) {
            nextStep = 0.5 * step;
        }
        if (nextStep < minEta) {
            // Recovery for a vanishing interpolation step.
            nextStep = 0.5;
        }
        step = nextStep;
    }

    // If every sampled state was invalid, keep cutting back instead of
    // silently restoring the original full step (the former OOFEM behavior).
    const bool hasValidTrial = std::isfinite(bestDirectionalMagnitude);
    if (!accepted && !hasValidTrial) {
        restoreUnperturbedState();
        dr.zero();
        etaValue = 0.0;
        status = ls_failed;
        OOFEM_LOG_DEBUG("LS: no admissible trial correction found\n");
        return CR_DIVERGED_ITS;
    }

    double acceptedStep = accepted ? step : bestStep;
    TrialMeasure finalTrial;
    if (!accepted || acceptedStep != step) {
        finalTrial = evaluateTrial(acceptedStep);
    } else {
        finalTrial = computeResidualMeasure();
    }
    if (!std::isfinite(finalTrial.directional) || !std::isfinite(finalTrial.norm)) {
        restoreUnperturbedState();
        dr.zero();
        etaValue = 0.0;
        status = ls_failed;
        return CR_DIVERGED_ITS;
    }

    dr = prescribedCorrection;
    dr.add(acceptedStep, freeCorrection);
    r = rb;
    etaValue = acceptedStep;
    status = ls_ok;
    OOFEM_LOG_DEBUG("LS: accepted eta=%e%s\n", acceptedStep,
                    accepted ? "" : " (best trial)");
    return CR_CONVERGED;
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
LineSearchNM :: initializeFrom(const std::shared_ptr<InputRecord> &ir)
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

    IR_GIVE_OPTIONAL_FIELD(ir, minEta, _IFT_LineSearchNM_lsearchmineta);
    minEta = std::clamp(minEta, 1.e-8, 0.5);

    IR_GIVE_OPTIONAL_FIELD(ir, max_iter, _IFT_LineSearchNM_lsearchmaxiter);
    max_iter = std::clamp(max_iter, 1, 50);

    //printf ("\nLineSearchNM::initializeFrom: tol=%e, ampl=%e, maxEta=%e\n",
    //    ls_tolerance, amplifFactor,maxEta);
}
} // end namespace oofem
