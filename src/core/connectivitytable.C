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

#include "connectivitytable.h"
#include "domain.h"
#include "element.h"
#include "dofmanager.h"
#include "feinterpol.h"
#include "intarray.h"
#include <set>


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
        for (int i=1; i<= nodes.giveSize(); i++) {
            const IntArray *candidates = this->giveDofManConnectivityArray(nodes.at(i));
            if (i == 1) {
                answer = *candidates; // first node gives preliminary set of elements
            } else { // other nodes refine the set
                for (int j=1; j<=answer.giveSize(); j++) { // loop over candidates 
                    if (candidates->contains(answer.at(j)) == false) {
                        answer.erase(j);
                        j--;
                    }
                }
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


void 
ConnectivityTable::buildSharedBoundaryEntities(Domain *d) {
    // sets of processed boundary entities per element
    std::vector<std::set<int>> processedBoundaryEdges;
    std::vector<std::set<int>> processedBoundarySurfaces;

    // local vars
    IntArray bnodes, bnodesSorted, neighbors;

    processedBoundaryEdges.resize(domain->giveNumberOfElements());
    processedBoundarySurfaces.resize(domain->giveNumberOfElements());

    // start with a loop over all elements
    for (auto &elem : domain->elementList) {
        // loop over all edges of the element
        // int nsd = elem->giveSpatialDimension();
        // process edges
        int nedges = elem->giveNumberOfEdges();
        for (int i = 1; i <= nedges; i++) {
            if (processedBoundaryEdges[elem->giveNumber()-1].find(i) != processedBoundaryEdges[elem->giveNumber()-1].end()) {
                continue; // edge already processed by neighbor element
            }
            // get element edge
            bnodes = elem->giveBoundaryEdgeNodes(i);
            for ( int j = 1; j <= bnodes.giveSize(); j++ ) {
                bnodes.at(j) = elem->giveDofManager(bnodes.at(j))->giveNumber();
            }
            bnodesSorted = bnodes;
            bnodesSorted.sort();

            // create a boundary entity
            std::unique_ptr<SharedBoundaryEntity> sbe = std::make_unique<SharedBoundaryEntity>();
            sbe->nodes = bnodes;
            SharedBoundaryEntity::elementRec er={elem->giveNumber(), i};
            sbe->elements.push_back(er);
            sbe->geomType = elem->giveEdgeGeometryType(i);
            sbe->spatialDimension = 1;
            std::size_t sbeIndex = sharedBoundaryEntities.size()+1; 
            // store the boundary entity ID on the element
            elem->setSharedEdgeID(i, sbeIndex);

            // now search element neighbors to find one sharing the same boundary nodes
            this->giveElementsWithNodes(neighbors, bnodes);
            // loop over elements sharing the boundary entity
            for (auto neighborelem : neighbors) {
                if (neighborelem == elem->giveNumber()) {
                    continue;
                }
                int neighborboundary = 0;
                IntArray neighborBoundaryNodes;
                for ( int j = 1; j <= domain->giveElement(neighborelem)->giveNumberOfEdges(); j++ ) {
                    bool equal = true;
                    neighborBoundaryNodes = domain->giveElement(neighborelem)->giveBoundaryEdgeNodes(j);// edge or surface?
                    for ( int k = 1; k <= neighborBoundaryNodes.giveSize(); k++ ) {
                        neighborBoundaryNodes.at(k) = domain->giveElement(neighborelem)->giveDofManager(neighborBoundaryNodes.at(k))->giveNumber();
                    }
                    
                    // compare bnodes (sorted) with neighborBoundaryNodes
                    for ( int k = 1; k <= neighborBoundaryNodes.giveSize(); k++ ) {
                        if ( !bnodesSorted.findSorted(neighborBoundaryNodes.at(k))) {
                            equal = false;
                            break;
                        }
                    }
                    if ( equal) {   
                        neighborboundary = j;
                        break;
                    }
                }
                if (neighborboundary) {
                    // store the boundary entity id on the neighbor element
                    domain->giveElement(neighborelem)->setSharedEdgeID(neighborboundary, sbeIndex);
                    SharedBoundaryEntity::elementRec er={neighborelem, neighborboundary};
                    sbe->elements.push_back(er);
                    processedBoundaryEdges[neighborelem-1].insert(neighborboundary);                    
                } else {
                    OOFEM_ERROR("Boundary entity not found in neighbor element");
                }
            } // end loop over neighbors
            sharedBoundaryEntities.emplace_back(std::move(sbe));
        } // end loop over element edges
    } // end loop over elements

    // @TODO now process surfaces

    //@DEBUG print BEs
#if 0
    for (auto& be: this->sharedBoundaryEntities) {
        printf("BE nodes ");
        for (auto n: be->nodes) {
            printf("%4d ", n);
        }
        printf("elements/sides ");
        for (auto& p: be->elements) {
            printf("%4d(%2d) ", p.elementID, p.boundaryID);
        }
        printf("\n");
    }
#endif
}

SharedBoundaryEntity* 
ConnectivityTable::giveBoundaryEntity(int id) {
    SharedBoundaryEntity *be = reinterpret_cast< SharedBoundaryEntity * >( sharedBoundaryEntities[id-1].get() );
    return be;
}


} // end namespace oofem
