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

#include "connectivitytable.h"
#include "domain.h"
#include "element.h"
#include "dofmanager.h"
#include "intarray.h"

namespace oofem {

void
ConnectivityTable :: reset()
{
    nodalConnectivityFlag = 0;
}

void
ConnectivityTable :: instanciateConnectivityTable()
//
// assembles table of nodal connectivities
//
{
    int ndofMan = domain->giveNumberOfDofManagers();
    int nelems = domain->giveNumberOfElements();
    IntArray dofManConnectivity(ndofMan);

    if ( nodalConnectivityFlag ) {
        return;                     // already initialized
    }

    OOFEM_LOG_INFO("ConnectivityTable: initializing\n");

    for ( auto &elem : domain->giveElements() ) {
        int nnodes = elem->giveNumberOfDofManagers();
        for ( int j = 1; j <= nnodes; j++ ) {
            int jnode = elem->giveDofManager(j)->giveNumber();
            dofManConnectivity.at(jnode)++;
        }
    }

    // allocate Nodal connectivity table for domain
    nodalConnectivity.resize(ndofMan);
    for ( int i = 0; i < ndofMan; i++ ) {
        nodalConnectivity[i].resize( dofManConnectivity(i) );
    }

    // build Nodal connectivity table for domain
    dofManConnectivity.zero();

    for ( int i = 1; i <= nelems; i++ ) {
        Element *ielem = domain->giveElement(i);
        int nnodes = ielem->giveNumberOfDofManagers();
        for ( int j = 1; j <= nnodes; j++ ) {
            int jnode = ielem->giveDofManager(j)->giveNumber();
            nodalConnectivity[jnode-1].at( ++dofManConnectivity.at(jnode) ) = i;
        }
    }

    nodalConnectivityFlag = 1;
}

const IntArray *
ConnectivityTable :: giveDofManConnectivityArray(int dofman)
{
    if ( nodalConnectivityFlag == 0 ) {
        this->instanciateConnectivityTable();
    }

    return &this->nodalConnectivity[dofman-1];
}


void
ConnectivityTable :: giveElementNeighbourList(IntArray &answer, IntArray &elemList)
{
    if ( nodalConnectivityFlag == 0 ) {
        this->instanciateConnectivityTable();
    }

    answer.resize(0);

    for ( auto &el_num : elemList ) {
        Element *ielem = domain->giveElement( el_num );
        int nnode = ielem->giveNumberOfDofManagers();
        for ( int j = 1; j <= nnode; j++ ) {
            int jnode = ielem->giveDofManager(j)->giveNumber();
            for ( auto &val : this->nodalConnectivity[jnode-1] ) {
                answer.insertSortedOnce( val );
            }
        }
    }
}


void
ConnectivityTable :: giveNodeNeighbourList(IntArray &answer, IntArray &nodeList)
{
    int nnodes = nodeList.giveSize();
    if ( nodalConnectivityFlag == 0 ) {
        this->instanciateConnectivityTable();
    }

    answer.resize(0);

    for ( int i = 1; i <= nnodes; i++ ) {
        int inode = nodeList.at(i);
        for ( auto &val : this->nodalConnectivity[inode-1] ) {
            answer.insertSortedOnce( val );
        }
    }
}
} // end namespace oofem
