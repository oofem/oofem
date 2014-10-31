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

#include "wallclockloadbalancermonitor.h"
#include "engngm.h"
#include "domain.h"
#include "timestep.h"
#include "element.h"
#include "mathfem.h"
#include "classfactory.h"

#include <mpi.h>

namespace oofem {

REGISTER_LoadBalancerMonitor(WallClockLoadBalancerMonitor);

LoadBalancerMonitor :: LoadBalancerDecisionType
WallClockLoadBalancerMonitor :: decide(TimeStep *tStep)
{
    int nproc = emodel->giveNumberOfProcesses();
    int myrank = emodel->giveRank();
    Domain *d = emodel->giveLoadBalancer()->giveDomain();
    int nelem;
    double *node_solutiontimes = new double [ nproc ];
    double *node_relcomppowers = new double [ nproc ];
    double *node_equivelements = new double [ nproc ];
    double min_st, max_st;
    double relWallClockImbalance;
    double absWallClockImbalance;
    double neqelems, sum_relcomppowers;

    if ( node_solutiontimes == NULL ) {
        OOFEM_ERROR("failed to allocate node_solutiontimes array");
    }

    if ( node_relcomppowers == NULL ) {
        OOFEM_ERROR("failed to allocate node_relcomppowers array");
    }

    if ( node_equivelements == NULL ) {
        OOFEM_ERROR("failed to allocate node_equivelements array");
    }


    // compute wall solution time of my node
    double mySolutionTime = emodel->giveTimer()->getWtime(EngngModelTimer :: EMTT_NetComputationalStepTimer);

#ifdef __LB_DEBUG
    // perturb solution time artificially if requested
    bool perturb = false;
    for ( auto perturbedStep: perturbedSteps ) {
        if ( perturbedStep.test( tStep->giveNumber() ) ) {
            perturb  = true;
            break;
        }
    }

    if ( perturb ) {
        mySolutionTime *= perturbFactor;
        OOFEM_LOG_RELEVANT("[%d] WallClockLoadBalancerMonitor: perturbed solution time by factor=%.2f\n", myrank, perturbFactor);
    }

#endif

    // collect wall clock computational time
    MPI_Allgather(& mySolutionTime, 1, MPI_DOUBLE, node_solutiontimes, 1, MPI_DOUBLE, MPI_COMM_WORLD);

    OOFEM_LOG_RELEVANT("\nLoadBalancer:: individual processor times [sec]: (");
    for ( int i = 0; i < nproc; i++ ) {
        OOFEM_LOG_RELEVANT(" %.3f", node_solutiontimes [ i ]);
    }

    OOFEM_LOG_RELEVANT(")\n");

    // detect imbalance
    min_st = max_st = node_solutiontimes [ 0 ];
    for ( int i = 0; i < nproc; i++ ) {
        min_st = min(min_st, node_solutiontimes [ i ]);
        max_st = max(max_st, node_solutiontimes [ i ]);
    }

    absWallClockImbalance = ( max_st - min_st );
    if ( min_st ) {
        relWallClockImbalance = ( ( max_st - min_st ) / min_st );
    } else {
        relWallClockImbalance = 0.0;
    }

    // update node (processor) weights

    // compute number or equivalent elements (equavalent element has computational weight equal to 1.0)
    nelem = d->giveNumberOfElements();
    neqelems = 0.0;
    for ( int ie = 1; ie <= nelem; ie++ ) {
        if ( d->giveElement(ie)->giveParallelMode() == Element_remote ) {
            continue;
        }

        neqelems += d->giveElement(ie)->predictRelativeComputationalCost();
    }

    // exchange number or equivalent elements
    MPI_Allgather(& neqelems, 1, MPI_DOUBLE, node_equivelements, 1, MPI_DOUBLE, MPI_COMM_WORLD);


    if ( !this->staticNodeWeightFlag ) {
        // compute relative computational powers (solution_time/number_of_equivalent_elements)
        for ( int i = 0; i < nproc; i++ ) {
            node_relcomppowers [ i ] = node_equivelements [ i ] / node_solutiontimes [ i ];
        }

        // normalize computational powers
        sum_relcomppowers = 0.0;
        for ( int i = 0; i < nproc; i++ ) {
            sum_relcomppowers += node_relcomppowers [ i ];
        }

        for ( int i = 0; i < nproc; i++ ) {
            nodeWeights(i) = node_relcomppowers [ i ] / sum_relcomppowers;
        }
    }

    // log equivalent elements on nodes
    OOFEM_LOG_RELEVANT("[%d] LoadBalancer:  node equivalent elements: ", myrank);
    for ( int i = 0; i < nproc; i++ ) {
        OOFEM_LOG_RELEVANT("%6d ", ( int ) node_equivelements [ i ]);
    }

    OOFEM_LOG_RELEVANT("\n");

    // log processor weights
    OOFEM_LOG_RELEVANT("[%d] LoadBalancer: updated proc weights: ", myrank);
    for ( int i = 0; i < nproc; i++ ) {
#ifdef __LB_DEBUG
        OOFEM_LOG_RELEVANT( "%22.15e ", nodeWeights(i) );
#else
        OOFEM_LOG_RELEVANT( "%4.3f ", nodeWeights(i) );
#endif
    }

    OOFEM_LOG_RELEVANT("\n");

    delete[] node_solutiontimes;
    delete[] node_relcomppowers;
    delete[] node_equivelements;

#ifdef __LB_DEBUG
    if ( recoveredSteps.giveSize() ) {
        // recover lb if requested
        int pos;
        if ( ( pos = recoveredSteps.findFirstIndexOf( tStep->giveNumber() ) ) ) {
            double procWeight, sumWeight = 0.0, *procWeights = new double [ nproc ];

            // assign prescribed processing weight
            procWeight = processingWeights.at(pos);
            OOFEM_LOG_RELEVANT("[%d] WallClockLoadBalancerMonitor: processing weight overriden by value=%e\n", myrank, procWeight);

            // exchange processing weights
            MPI_Allgather(& procWeight, 1, MPI_DOUBLE, procWeights, 1, MPI_DOUBLE, MPI_COMM_WORLD);
            for ( int i = 0; i < nproc; i++ ) {
                nodeWeights(i) = procWeights [ i ];
                sumWeight += procWeights [ i ];
            }

            delete[] procWeights;

            if ( fabs(sumWeight - 1.0) > 1.0e-10 ) {
                OOFEM_ERROR("[%d] processing weights do not sum to 1.0 (sum = %e)\n", sumWeight);
            }

            OOFEM_LOG_RELEVANT("[%d] LoadBalancer: wall clock imbalance rel=%.2f\%,abs=%.2fs, recovering load\n", myrank, 100 * relWallClockImbalance, absWallClockImbalance);
            return LBD_RECOVER;
        } else {
            OOFEM_LOG_RELEVANT("[%d] LoadBalancer: wall clock imbalance rel=%.2f\%,abs=%.2fs, continuing\n", myrank, 100 * relWallClockImbalance, absWallClockImbalance);
            return LBD_CONTINUE;
        }
    }

#endif

    // decide
    if ( ( tStep->giveNumber() % this->lbstep == 0 ) &&
        ( ( absWallClockImbalance > this->absWallClockImbalanceTreshold ) ||
         ( ( relWallClockImbalance > this->relWallClockImbalanceTreshold ) && ( absWallClockImbalance > this->minAbsWallClockImbalanceTreshold ) ) ) ) {
        OOFEM_LOG_RELEVANT("[%d] LoadBalancer: wall clock imbalance rel=%.2f\%,abs=%.2fs, recovering load\n", myrank, 100 * relWallClockImbalance, absWallClockImbalance);
        return LBD_RECOVER;
    } else {
        OOFEM_LOG_RELEVANT("[%d] LoadBalancer: wall clock imbalance rel=%.2f\%,abs=%.2fs, continuing\n", myrank, 100 * relWallClockImbalance, absWallClockImbalance);
        return LBD_CONTINUE;
    }
}


IRResultType
WallClockLoadBalancerMonitor :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    result = LoadBalancerMonitor :: initializeFrom(ir);

    IR_GIVE_OPTIONAL_FIELD(ir, relWallClockImbalanceTreshold, _IFT_WallClockLoadBalancerMonitor_relwct);
    IR_GIVE_OPTIONAL_FIELD(ir, absWallClockImbalanceTreshold, _IFT_WallClockLoadBalancerMonitor_abswct);
    IR_GIVE_OPTIONAL_FIELD(ir, minAbsWallClockImbalanceTreshold, _IFT_WallClockLoadBalancerMonitor_minwct);
    IR_GIVE_OPTIONAL_FIELD(ir, lbstep, _IFT_WallClockLoadBalancerMonitor_lbstep);

#ifdef __LB_DEBUG
    perturbedSteps.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, perturbedSteps, _IFT_WallClockLoadBalancerMonitor_perturbedsteps);
    perturbFactor = 1.0;
    IR_GIVE_OPTIONAL_FIELD(ir, perturbFactor, _IFT_WallClockLoadBalancerMonitor_perturbfactor);

    recoveredSteps.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, recoveredSteps, _IFT_WallClockLoadBalancerMonitor_recoveredsteps);
    processingWeights.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, processingWeights, _IFT_WallClockLoadBalancerMonitor_processingweights);
    if ( recoveredSteps.giveSize() != processingWeights.giveSize() ) {
        OOFEM_ERROR("mismatch size of lbrecoveredsteps and lbprocessingweights");
    }

#endif

    return result;
}

}
