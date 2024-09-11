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

ConnectivityTable::ConnectivityTable(Domain * d) : domain(d), nodalConnectivity(), nodalConnectivityFlag(0), elementColoring(), elementColoringFlag(false)  
{
#ifdef _OPENMP
    omp_init_lock(&initLock);
#endif
}

void
ConnectivityTable :: reset()
{
    nodalConnectivityFlag = 0;
    elementColoringFlag = false;
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
#ifdef _OPENMP
    omp_set_lock(&initLock); // if not initialized yet; one thread can proceed with init; others have to wait until init completed
    if ( this->nodalConnectivityFlag ) {
        omp_unset_lock(&initLock); 
        return;
    }
#endif   
//    OOFEM_LOG_INFO("ConnectivityTable: initializing\n");

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
        nodalConnectivity[i].resize( dofManConnectivity[i] );
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
    
 #ifdef _OPENMP
    omp_unset_lock(&initLock);
#endif
 
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
ConnectivityTable :: giveElementNeighbourList(IntArray &answer, const IntArray &elemList)
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

void 
ConnectivityTable::giveElementsWithNodes(IntArray &answer, const IntArray& nodes)
{
    // loop over individual node's connectivity arrays and find intersection (common elements)
    answer.resize(0);
    if (nodes.giveSize()) {
        const IntArray *candidates = this->giveDofManConnectivityArray(nodes.at(1));

        for (int i=1; i<= candidates->giveSize(); i++) {
            IntArray enodes = this->domain->giveElement(candidates->at(i))->giveDofManArray();
            enodes.sort();
            bool found = true;            
            for (auto node : nodes) {
                if (enodes.containsSorted(node) == false) {
                    found = false;
                    break;
                }
            }
            if (found) {
                answer.followedBy(candidates->at(i));
            }
        }
    }
}

void 
ConnectivityTable::buildElementColoring() {
    // greedy coloring algorithm
    // loop over elements
    if (elementColoringFlag) {
        return;
    } else {
        int nelems = domain->giveNumberOfElements();
        elementColoring.resize(nelems);
        elementColoring.zero();
        for (auto &elem : domain->elementList) {
            IntArray neighbors;
            giveElementNeighbourList(neighbors, IntArray({elem->giveNumber()}));
            IntArray neighborsColors;
            for (auto val : neighbors) {
                neighborsColors.insertSortedOnce(elementColoring.at(val));
            }
            // find first unused color code in neighborsColor array
            for (int c=1; true; c++) {
                if (!neighborsColors.containsSorted(c)) {
                    elementColoring.at(elem->giveNumber()) = c;
                    break;
                }
            }
        }
        elementColoringFlag = true;
    }
}

int 
ConnectivityTable::getElementColor(int e) {
    if (!elementColoringFlag) {
        buildElementColoring();
    }
    return elementColoring.at(e);
}

} // end namespace oofem
