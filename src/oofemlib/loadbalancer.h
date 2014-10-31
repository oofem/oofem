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

#ifndef loadbalancer_h
#define loadbalancer_h

#include "oofemcfg.h"
#include "inputrecord.h"
#include "floatarray.h"

#include <vector>
#include <memory>

///@name Input fields for LoadBalancer
//@{
#define _IFT_LoadBalancer_wtp "wtp"
#define _IFT_LoadBalancerMonitor_nodeWeightMode "nodeweightmode"
#define _IFT_LoadBalancerMonitor_initialnodeweights "nw"
//@}

namespace oofem {
class Domain;
class EngngModel;
class ProcessCommunicator;
class TimeStep;
class IntArray;

 #define MIGRATE_LOAD_TAG       9998

/**
 * Abstract base class representing general load balancer monitor. The task of the monitor is to
 * detect the imbalance and to make the decision, whether to redistribute the work or to continue
 * with existing partitioning.
 * It provides partition weights, reflecting their relative computational performance. These weights should
 * be continuously updated to reflect changing work load during solution process.
 */
class OOFEM_EXPORT LoadBalancerMonitor
{
protected:
    EngngModel *emodel;
    FloatArray nodeWeights;
    bool staticNodeWeightFlag;
public:
    enum LoadBalancerDecisionType { LBD_CONTINUE, LBD_RECOVER };

    LoadBalancerMonitor(EngngModel * em): emodel(em) { }
    virtual ~LoadBalancerMonitor() { }

    /// Initializes receiver according to object description stored in input record.
    virtual IRResultType initializeFrom(InputRecord *ir);

    /**@name Load evaluation and imbalance detection methods*/
    //@{
    /// Returns flag indicating whether rebalancing is necessary; should update node weights as well.
    virtual LoadBalancerDecisionType decide(TimeStep *) = 0;
    /// Returns processor weights; the larger weight means more powerful node, sum of weights should equal to one.
    const FloatArray & giveProcessorWeights() { return nodeWeights; }
    //@}

    /// Returns class name of the receiver.
    virtual const char *giveClassName() const = 0;
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
class OOFEM_EXPORT LoadBalancer
{
public:
    /**
     * Describes the state of dofmanager after load balancing
     * on the local partition.
     */
    enum DofManMode {
        DM_NULL,   ///< Undefined (undetermined) state, if assigned means internal error.
        DM_Local,  ///< Local dofman that remains local.
        DM_Shared, ///< Shared dofman that remains shared.
        DM_Remote, ///< Local dofman that became remote (became local on remote partition).
    };
protected:
    Domain *domain;

public:

    LoadBalancer(Domain * d);
    virtual ~LoadBalancer() { }



    /**@name Work transfer calculation methods  */
    //@{
    virtual void calculateLoadTransfer() = 0;
    //@}

    /**@name Work migration methods */
    //@{
    void migrateLoad(Domain *d);
    //@}

    /// Print receiver statistics
    virtual void  printStatistics() const;

    /**@name Query methods after work transfer calculation */
    //@{
    /// Returns the label of dofmanager after load balancing.
    virtual DofManMode giveDofManState(int idofman) = 0;

    /// Returns the partition list of given dofmanager after load balancing.
    virtual IntArray *giveDofManPartitions(int idofman) = 0;

    /// Returns the new partition number assigned to local element after LB.
    virtual int giveElementPartition(int ielem) = 0;

    //@}
    ///Initializes receiver according to object description stored in input record.
    virtual IRResultType initializeFrom(InputRecord *ir);

    /// Returns reference to its domain.
    Domain *giveDomain() { return domain; }
    /// sets associated Domain
    virtual void setDomain(Domain *d) { this->domain = d; }
    /// Returns class name of the receiver.
    virtual const char *giveClassName() const = 0;

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
        WorkTransferPlugin(LoadBalancer * _lb);
        virtual ~WorkTransferPlugin();

        /**
         *  initializes receiver; should be called before any work transfer.
         *  Current implementation assembles for each local element the list
         *  of contributing global element numbers.
         *  This is extracted from IP nonlocal tables;
         */
        virtual void init(Domain *d) = 0;
        /**
         * Migrates necessary local elements to remote processors, where they
         * become remote elements needed to efficiently handle nonlocal dependencies.
         *
         * This involves several steps:
         * - send and receive nonlocElementDependencyArry of migrating regular
         *   elements to remote partition
         * - build domain nonlocal element dependency list.
         * - then exclude local elements - what remains are unsatisfied
         *   emote dependencies that have to be broadcasted
         *   and received from partitions owning relevant elements
         * - transfer of local elements and nodes to remote partitions
         *   (remote elements and null dofmans)
         */
        virtual void migrate() = 0;
        /**
         * Called after all wtps migrated their data.
         * Intended to update local data structure.
         * Current implementations rebuilds the nonlocal integration point tables.
         */
        virtual void update() = 0;
    };

protected:
    /// List of work transfer plugins.
    std :: vector< std :: unique_ptr< WorkTransferPlugin > >wtpList;
};

/*
class OOFEM_EXPORT LoadBalancerElementInterface : public Interface
{
public:
    LoadBalancerElementInterface() { }

    virtual double predictRelativeComputationalCost();
};
*/
} // end namespace oofem

#endif // loadbalancer_h
