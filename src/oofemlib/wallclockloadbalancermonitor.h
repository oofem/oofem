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

#ifndef wallclockloadbalancer_h
#define wallclockloadbalancer_h

#include "loadbalancer.h"

#define __LB_DEBUG
#ifdef __LB_DEBUG
 #include <list>
 #include "range.h"
 #include "intarray.h"
 #include "floatarray.h"
#endif

///@name Input fields for WallClockLoadBalancerMonitor
//@{
#define _IFT_WallClockLoadBalancerMonitor_Name "wallclock"
#define _IFT_WallClockLoadBalancerMonitor_relwct "relwct"
#define _IFT_WallClockLoadBalancerMonitor_abswct "abswct"
#define _IFT_WallClockLoadBalancerMonitor_minwct "minwct"
#define _IFT_WallClockLoadBalancerMonitor_lbstep "lbstep"
#define _IFT_WallClockLoadBalancerMonitor_perturbedsteps "lbperturbedsteps"
#define _IFT_WallClockLoadBalancerMonitor_perturbfactor "lbperturbfactor"
#define _IFT_WallClockLoadBalancerMonitor_recoveredsteps "lbrecoveredsteps"
#define _IFT_WallClockLoadBalancerMonitor_processingweights "lbprocessingweights"
//@}

namespace oofem {

/**
 * Implementation of simple wall-clock based monitor.
 * It detect imbalance based on wall clock difference required for solution step
 * on particular nodes. When difference in wall clock solution times is greater
 * than a threshold value, the load migration is performed.
 */
class OOFEM_EXPORT WallClockLoadBalancerMonitor : public LoadBalancerMonitor
{
protected:
    /// Declares min abs imbalance to perform relative imbalance check.
    double relWallClockImbalanceTreshold, absWallClockImbalanceTreshold, minAbsWallClockImbalanceTreshold;
    /// The rebalancing done every lbstep.
    int lbstep;
 #ifdef __LB_DEBUG
    /// List of steps with perturbed balancing.
    std :: list< Range >perturbedSteps;
    /// Perturbing factor.
    double perturbFactor;
    /// list of step at which to performed lb recovery.
    IntArray recoveredSteps;
    /// processing weights for lb recovery.
    FloatArray processingWeights;
 #endif
public:
    WallClockLoadBalancerMonitor(EngngModel * em) : LoadBalancerMonitor(em) {
        relWallClockImbalanceTreshold = 0.1;
        absWallClockImbalanceTreshold = 10.0;
        minAbsWallClockImbalanceTreshold = 0.0;
        lbstep = 5;
    }

    LoadBalancerDecisionType decide(TimeStep *);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual const char *giveClassName() const { return "WallClockLoadBalancerMonitor"; }
};

} // end namespace oofem
#endif //wallclockloadbalancer_h