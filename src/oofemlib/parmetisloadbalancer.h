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

#ifndef parmetisloadbalancer_h
#define parmetisloadbalancer_h

#include "loadbalancer.h"
#include "intarray.h"

#include <parmetis.h>
#include <vector>

#define _IFT_ParmetisLoadBalancer_Name "parmetis"

namespace oofem {
/**
 * End-of-data marker, used to identify end of data stream received.
 * The value should not conflict with any globnum id.
 */
 #define PARMETISLB_END_DATA -1
 #define SHARED_DOFMAN_PARTITIONS_TAG 9998

/**
 * ParMetis load balancer.
 */
class OOFEM_EXPORT ParmetisLoadBalancer : public LoadBalancer
{
protected:
    /// Element numbering maps.
    IntArray gToLMap, lToGMap;
    idx_t *elmdist;
    int myGlobNumOffset;
    /// Partition weights (user input).
    real_t *tpwgts;
    /// Array of DofManMode(s).
    IntArray dofManState;
    /// Array of dof man partitions.
    std :: vector< IntArray >dofManPartitions;
    /// Partition vector of the locally-stored elements.
    IntArray elementPart;

public:
    ParmetisLoadBalancer(Domain * d);
    virtual ~ParmetisLoadBalancer();

    virtual void calculateLoadTransfer();

 #if 1
    virtual DofManMode giveDofManState(int idofman);
    virtual IntArray *giveDofManPartitions(int idofman);
    virtual int giveElementPartition(int ielem);
 #endif
protected:
    void handleMasterSlaveDofManLinks();

    void initGlobalParmetisElementNumbering();
    int  giveLocalElementNumber(int globnum) { return gToLMap.at(globnum - myGlobNumOffset); }
    int  giveGlobalElementNumber(int locnum) { return lToGMap.at(locnum); }

    /**
     * Label local partition nodes (the nodes that are local or shared).
     * Labeling consist of assigning corresponding id that characterize the
     * status of local dof manager after balancing the load. Labeling determines
     * which of local nodes remain local, or became local on other partition,
     * or became shared, etc.
     */
    void labelDofManagers();
    int  determineDofManState(int idofman, int myrank, int npart, IntArray *dofManPartitions);

    int packSharedDmanPartitions(ProcessCommunicator &pc);
    int unpackSharedDmanPartitions(ProcessCommunicator &pc);
    void addSharedDofmanPartitions(int _locnum, IntArray _partitions);

    virtual const char *giveClassName() const { return "ParmetisLoadBalancer"; }
};
} // end namespace oofem

#endif // parmetisloadbalancer_h
