/* $Header: /home/cvs/bp/oofem/oofemlib/src/combuff.h,v 1.5 2003/04/06 14:08:23 bp Exp $ */
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
#ifndef loadbalancer_h
#define loadbalancer_h

#ifdef __PARALLEL_MODE
#include "inputrecord.h"
#include "interface.h"
#include "clock.h"
#include "alist.h"

#ifdef __LB_DEBUG
#include "dynalist.h"
#include "range.h"
#endif

namespace oofem {

#define __LB_DEBUG

class Domain;
class EngngModel;
class ProcessCommunicator;
class TimeStep;

#define MIGRATE_LOAD_TAG       9998
/**
 * End-of-data marker, used to identify end of data stream received.
 * The value should not conflict with any classType value
 */
#define LOADBALANCER_END_DATA -1

/**
 * Abstract base class representing general load balancer monitor. The task of the monitor is to
 * detect the imbalance and to make the decision, whether to redistribute the work or to continue
 * with existing partitioning.
 * It provides partition weights, reflecting their relative computational performance. These weights should
 * be continuosly updated to reflect changing work load during solution process.
 */
class LoadBalancerMonitor
{
protected:
    EngngModel *emodel;
    FloatArray nodeWeights;
    bool staticNodeWeightFlag;
public:
    enum LoadBalancerDecisionType { LBD_CONTINUE, LBD_RECOVER };

    LoadBalancerMonitor(EngngModel *em) { emodel = em; }
    virtual ~LoadBalancerMonitor() { }

    ///Initializes receiver acording to object description stored in input record.
    virtual IRResultType initializeFrom(InputRecord *ir);

    /**@name Load evaluation and imbalance detection methods*/
    //@{
    /// returns flag indicating whether rebalancing is necessary; should update node weights as well
    virtual LoadBalancerDecisionType decide(TimeStep *) = 0;
    /// Returns processor weights; the larger weight means more powerfull node, sum of weights should equal to one.
    void giveProcessorWeights(FloatArray &answer) { answer = nodeWeights; }
    //@}

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "LoadBalancerMonitor"; }
};

/**
 * Implementation of simple wall-clock based monitor.
 * It detect imbalance based on waal clock diference required for slotion step
 * on particular nodes. When difference in wall clock solution times is greater
 * than a treshold value, the load migration is performed.
 */
class WallClockLoadBalancerMonitor : public LoadBalancerMonitor
{
protected:
  /// minAbsWallClockImbalanceTreshold declares min abs imbalance to perform relative imbalance check
  double relWallClockImbalanceTreshold, absWallClockImbalanceTreshold, minAbsWallClockImbalanceTreshold;
    /// the rebalancing done every lbstep
    int lbstep;
#ifdef __LB_DEBUG
    // list of steps with perturbed ballancing
    dynaList< Range >perturbedSteps;
    // perturbing factor
    double perturbFactor;
    // list of step at which to performed lb recovery
    IntArray recoveredSteps;
    // processing weights for lb recovery
    FloatArray processingWeights;
#endif
public:
    WallClockLoadBalancerMonitor(EngngModel *em) : LoadBalancerMonitor(em) {
        relWallClockImbalanceTreshold = 0.1;
        absWallClockImbalanceTreshold = 10.0;
        minAbsWallClockImbalanceTreshold = 0.0;
        lbstep = 5;
    }
    LoadBalancerDecisionType decide(TimeStep *);
    ///Initializes receiver acording to object description stored in input record.
    virtual IRResultType initializeFrom(InputRecord *ir);
};


/**
 * Abstract base class representing general load balancer. The task of load balancer is to
 * recover load balance when running in parallel. This is achieved by moving work from busy
 * nodes to other nodes to achieve an equal distribution of work.
 * In general load balancer should repartition the problem domain, taking into account several
 * criteria:
 * - It should take into account different computational requirement of different elements
 * - The new partitioning should minimize the cut (to minimize the communication)
 * - The new partitioning should minimize data movement (the cost of repartitioning) by
 * preserving the locality as much as possible. In other words the new and existing partitioning
 * should be "similar".
 */
class LoadBalancer
{
public:
    /**
     * Describes the state of dofmanager after load balancing
     * on the local partition:
     * DM_NULL  - undefined (undetermined) state, if assigned means internal error
     * DM_Local - local dofman that remains local
     * DM_Remote- local dofman that becames remote (becames local on remote partition)
     * //DM_SharedNew - local shared that became shared
     * DM_Shared- shared dofman that remains shared.
     * //DM_SharedExlude - shared dofman that remains shared,
     *                   possibly with changed partitions, but local partition
     *                   is no more in shared list (should be exluded on remote partitions).
     * //DM_SharedUpdate - Shared dofman that remains shared, the partition list may changed.
     */
    enum DofManMode { DM_NULL, DM_Local, DM_Shared, DM_Remote };
protected:
    Domain *domain;

public:

    LoadBalancer(Domain *d);
    virtual ~LoadBalancer() { }



    /**@name Work transfer calculation methods  */
    //@{
    virtual void calculateLoadTransfer() = 0;
    //@}

    /**@name Work migration methods  */
    //@{
    void migrateLoad(Domain *d);
    //@}

    /// Print receiver statistics
    virtual void  printStatistics() const;

    /**@name Query methods after work transfer calculation */
    //@{
    /// Returns the label of dofmanager after load balancing
    virtual DofManMode giveDofManState(int idofman) = 0;

    /// Returns the partition list of given dofmanager after load balancing
    virtual IntArray *giveDofManPartitions(int idofman) = 0;

    /// Returns the new partition number assigned to local element after LB
    virtual int giveElementPartition(int ielem) = 0;

    //@}
    ///Initializes receiver acording to object description stored in input record.
    virtual IRResultType initializeFrom(InputRecord *ir);

    /// Returns reference to its domain
    Domain *giveDomain() { return domain; }
    /// sets associated Domain
    virtual void         setDomain(Domain *d) { this->domain = d; }
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "LoadBalancer"; }

protected:

    int packMigratingData(Domain *, ProcessCommunicator &pc);
    int unpackMigratingData(Domain *, ProcessCommunicator &pc);
    void deleteRemoteDofManagers(Domain *);
    void deleteRemoteElements(Domain *);
    void initializeWtp(IntArray &wtp);

public:

    class WorkTransferPlugin
    {
protected:
        LoadBalancer *lb;
public:
        WorkTransferPlugin(LoadBalancer *_lb);
        virtual ~WorkTransferPlugin();

        /**
         *  initializes receiver; should be called before any work transfer.
         *  Current implementation assembles for each local element the list
         *  of contributing global element numbers.
         *  This is extracted from IP nonlocal tables;
         */
        virtual void init(Domain *d) = 0;
        /*
         * Migrates necessary local elements to remote processors, where they
         * become remote elements needed to efficiently handle nonlocal dependencies.
         *
         * This involves several steps:
         * - send and receive nonlocElementDependencyArry of migrating regular
         * elements to remote partition
         * - build domain nonlocal element dependency list.
         * - then exclude local elements - what remains are unsatisfied
         * remote dependencies that have to be broadcasted
         * and received from partitions owning relevant elements
         * - transfer of local elements and nodes to remote partitions
         * (remote elements and null dofmans)
         */
        virtual void migrate() = 0;
        /*
         * Called after all wtps migrated their data.
         * Intended to update local data structure.
         * Current implementations rebuilds the nonlocal integration point tables.
         */
        virtual void update() = 0;
    };

protected:
    /// list of work transfer plugins
    AList< WorkTransferPlugin >wtpList;
};


class LoadBalancerElementInterface : public Interface
{
public:
    LoadBalancerElementInterface() { }

    virtual double predictRelativeComputationalCost();
};

} // end namespace oofem
#endif
#endif // loadbalancer_h
