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

#ifndef nonlocalmatwtp_h
#define nonlocalmatwtp_h

#include "loadbalancer.h"

#include <map>
#include <list>
#include <vector>
#include <set>

namespace oofem {
/**
 * End-of-data marker, used to identify end of data stream received.
 * The value should not conflict with any classType value and any possible globnum id.
 */
 #define NonlocalMaterialWTP_END_DATA -1

 #define MIGRATE_NONLOCALDEP_TAG 4008
 #define MIGRATE_REMOTE_ELEMENTS_TAG 4009


class GaussPoint;

/**
 * Class implements Work Transfer Plugin that introduces
 * extension to efficiently handle nonlocal dependency of
 * nonlocal materials via remote elements.
 */
class OOFEM_EXPORT NonlocalMaterialWTP : public LoadBalancer :: WorkTransferPlugin
{
protected:
    std :: map< int, IntArray >nonlocElementDependencyMap;
    // local to_Send_List for each partition
    std :: vector< std :: list< int > >toSendList;

public:
    NonlocalMaterialWTP(LoadBalancer * _lb) :
        LoadBalancer :: WorkTransferPlugin(_lb) { }

    /**
     * Initializes receiver; should be called before any work transfer.
     * Current implementation assembles for each local element the list
     * of contributing global element numbers. This is extracted from IP nonlocal tables;
     */
    void init(Domain *d);
    /**
     * Migrates necessary local elements to remote processors, where they
     * become remote elements needed to efficiently handle nonlocal dependencies.
     *
     * This involves several steps:
     * - send and receive nonlocElementDependencyArry of migrating regular
     *   elements to remote partition
     * - build domain nonlocal element dependency list. Then exclude local elements - what remains are unsatisfied remote dependencies that have to be broadcasted and received from partitions owning relevant elements
     * - transfer of local elements and nodes to remote partitions (remote elements and null dofmans)
     */
    void migrate();
    /**
     * Called after all wtps migrated their data. Intended to update local data structure.
     * Current implementations rebuilds the nonlocal integration point tables.
     */
    void update();

protected:
    void giveElementNonlocalDepArry(IntArray &answer, Domain *d, int num);
    int packMigratingElementDependencies(Domain *d, ProcessCommunicator &pc);
    int unpackMigratingElementDependencies(Domain *d, ProcessCommunicator &pc);
    int packRemoteElements(Domain *d, ProcessCommunicator &pc);
    int unpackRemoteElements(Domain *d, ProcessCommunicator &pc);
    /**
     * Plugin for Element::ipEvaluator service to compile nonlocal dependency
     * array for given element (in case of nonlocal material)
     */
    void giveNonlocalDepArryElementPlugin(GaussPoint *gp, std :: set< int > &s);
    /**
     * Rebuilds nonlocal integration tables in element gauss points using
     * information in nonlocElementDependencyMap, this should be significantly
     * faster than rebuild the whole information from scratch.
     */
    void fastRebuildNonlocalTables();
    /**
     * Element plugin for ipEvaluator service to update nonlocal tables of its integration
     * points from given map, containing element dependency.
     */
    void fastElementIPNonlocTableUpdater(GaussPoint *gp, IntArray &map);
};
} // end namespace oofem

#endif // nonlocalmatwtp_h
