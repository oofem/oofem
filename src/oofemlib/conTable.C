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

#include <cstdio>
#include <cstdlib>
#include <set>

#include "domain.h"
#include "element.h"
#include "dofmanager.h"
#include "conTable.h"

namespace oofem {
ConnectivityTable :: ~ConnectivityTable()
// destructor
{
    this->reset();
}

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
    int i, j, jnode, nnodes;
    Element *ielem;
    IntArray dofManConnectivity(ndofMan);

    if ( nodalConnectivityFlag ) {
        return;                     // already initialized
    }

    OOFEM_LOG_INFO("ConnectivityTable: initializing\n");

    for ( i = 1; i <= nelems; i++ ) {
        ielem = domain->giveElement(i);
        nnodes = ielem->giveNumberOfDofManagers();
        for ( j = 1; j <= nnodes; j++ ) {
            jnode = ielem->giveDofManager(j)->giveNumber();
            dofManConnectivity.at(jnode)++;
        }
    }

    // allocate Nodal connectivity table for domain
    nodalConnectivity.growTo(ndofMan);
    for ( i = 1; i <= ndofMan; i++ ) {
        nodalConnectivity.put( i, new IntArray( dofManConnectivity.at(i) ) );
    }

    // build Nodal connectivity table for domain
    dofManConnectivity.zero();

    for ( i = 1; i <= nelems; i++ ) {
        ielem = domain->giveElement(i);
        nnodes = ielem->giveNumberOfDofManagers();
        for ( j = 1; j <= nnodes; j++ ) {
            jnode = ielem->giveDofManager(j)->giveNumber();
            nodalConnectivity.at(jnode)->at( ++dofManConnectivity.at(jnode) ) = i;
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

    return this->nodalConnectivity.at(dofman);
}


void
ConnectivityTable :: giveElementNeighbourList(IntArray &answer, IntArray &elemList)
{
    int i, j, k, nnode, jnode, nelems = elemList.giveSize();
    Element *ielem;
    if ( nodalConnectivityFlag == 0 ) {
        this->instanciateConnectivityTable();
    }

    std :: set< int >neighbours;

    for ( i = 1; i <= nelems; i++ ) {
        ielem = domain->giveElement( elemList.at(i) );
        nnode = ielem->giveNumberOfDofManagers();
        for ( j = 1; j <= nnode; j++ ) {
            jnode = ielem->giveDofManager(j)->giveNumber();
            for ( k = 1; k <= this->nodalConnectivity.at(jnode)->giveSize(); k++ ) {
                neighbours.insert( this->nodalConnectivity.at(jnode)->at(k) );
            }
        }
    }

    answer.resize( neighbours.size() );
    std :: set< int > :: iterator pos;
    for ( pos = neighbours.begin(), i = 1; pos != neighbours.end(); ++pos, i++ ) {
        answer.at(i) = * pos;
    }
}


void
ConnectivityTable :: giveNodeNeighbourList(IntArray &answer, IntArray &nodeList)
{
    int i, k, inode, nnodes = nodeList.giveSize();
    if ( nodalConnectivityFlag == 0 ) {
        this->instanciateConnectivityTable();
    }

    std :: set< int >neighbours;

    for ( i = 1; i <= nnodes; i++ ) {
        inode = nodeList.at(i);
        for ( k = 1; k <= this->nodalConnectivity.at(inode)->giveSize(); k++ ) {
            neighbours.insert( this->nodalConnectivity.at(inode)->at(k) );
        }
    }

    answer.resize( neighbours.size() );
    std :: set< int > :: iterator pos;
    for ( pos = neighbours.begin(), i = 1; pos != neighbours.end(); ++pos, i++ ) {
        answer.at(i) = * pos;
    }
}
} // end namespace oofem
