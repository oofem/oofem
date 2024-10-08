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

#include "set.h"
#include "error.h"
#include "intarray.h"
#include "inputrecord.h"
#include "domain.h"
#include "element.h"
#include "node.h"
#include "range.h"
#include "feinterpol.h"
#include "feinterpol2d.h"
#include "feinterpol3d.h"
#include "contextioerr.h"
#include "dynamicinputrecord.h"
#include <list>

namespace oofem {
void Set :: initializeFrom(InputRecord &ir)
{
    FEMComponent :: initializeFrom(ir);

    IntArray inputNodes;
    std :: list< Range >inputNodeRanges;
    if ( ir.hasField(_IFT_Set_allNodes) ) { // generate a list with all the node numbers
       this->nodes.enumerate(this->giveDomain()->giveNumberOfDofManagers()); 
    } else {
        IR_GIVE_OPTIONAL_FIELD(ir, inputNodes, _IFT_Set_nodes);
        IR_GIVE_OPTIONAL_FIELD(ir, inputNodeRanges, _IFT_Set_nodeRanges);
        this->computeIntArray(this->nodes, inputNodes, inputNodeRanges);
    }

    if ( ir.hasField(_IFT_Set_allElements) ) { // generate a list with all the element numbers
        this->elements.enumerate(this->giveDomain()->giveNumberOfElements());
        mElementListIsSorted = false;
    } else {
        IntArray inputElements;
        std :: list< Range >inputElementRanges;
        IR_GIVE_OPTIONAL_FIELD(ir, inputElements, _IFT_Set_elements);
        IR_GIVE_OPTIONAL_FIELD(ir, inputElementRanges, _IFT_Set_elementRanges);
        this->computeIntArray(this->elements, inputElements, inputElementRanges);
        mElementListIsSorted = false;
    }


    this->elementBoundaries.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, this->elementBoundaries, _IFT_Set_elementBoundaries);

    this->elementEdges.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, this->elementEdges, _IFT_Set_elementEdges);

    this->elementSurfaces.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, this->elementSurfaces, _IFT_Set_elementSurfaces);

    this->elementInternalNodes.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, this->elementInternalNodes, _IFT_Set_internalElementNodes);
#if 0
    this->inputRec = ir.clone();
#endif
}


void Set :: giveInputRecord(DynamicInputRecord &input)
{
    input.setRecordKeywordField( _IFT_Set_Name, this->giveNumber() );

    if ( this->giveNodeList().giveSize() ) {
        input.setField(this->nodes, _IFT_Set_nodes);
    }
    if ( this->giveElementList().giveSize() ) {
        input.setField(this->elements, _IFT_Set_elements);
    }
    if ( this->giveBoundaryList().giveSize() ) {
        input.setField(this->elementBoundaries, _IFT_Set_elementBoundaries);
    }
    if ( this->giveEdgeList().giveSize() ) {
        input.setField(this->elementEdges, _IFT_Set_elementEdges);
    }
    if ( this->giveSurfaceList().giveSize() ) {
        input.setField(this->elementSurfaces, _IFT_Set_elementSurfaces);
    }
    if ( this->giveInternalElementDofManagerList().giveSize() ) {
        input.setField(this->elementInternalNodes, _IFT_Set_internalElementNodes);
    }
}


void Set :: computeIntArray(IntArray &answer, const IntArray &specified, std :: list< Range >ranges)
{
    // Find the max value;
    int maxIndex = specified.giveSize() == 0 ? 0 : specified.maximum();
    for ( auto &range: ranges ) {
        if ( range.giveEnd() > maxIndex ) {
            maxIndex = range.giveEnd();
        }
    }
    IntArray afflictedNodes(maxIndex);
    afflictedNodes.zero();

    for ( int i = 1; i <= specified.giveSize(); ++i ) {
        afflictedNodes.at( specified.at(i) ) = 1;
    }

    for ( auto &range: ranges ) {
        for ( int i = range.giveStart(); i <= range.giveEnd(); ++i ) {
            afflictedNodes.at(i) = 1;
        }
    }
    answer.findNonzeros(afflictedNodes);
}


const IntArray &Set :: giveElementList() { return elements; }

const IntArray &Set :: giveBoundaryList() { return elementBoundaries; }

const IntArray &Set :: giveEdgeList() { return elementEdges; }

const IntArray &Set :: giveSurfaceList() { return elementSurfaces; }  

const IntArray &Set :: giveInternalElementDofManagerList() { return this->elementInternalNodes;}

const IntArray &Set :: giveNodeList()
{
    // Lazy evaluation, we compute the unique set of nodes if needed (and store it).
    if ( this->totalNodes.giveSize() == 0 ) {
        IntArray afflictedNodes( this->domain->giveNumberOfDofManagers() );
        afflictedNodes.zero();
        for ( int ielem = 1; ielem <= this->elements.giveSize(); ++ielem ) {
            auto e = this->domain->giveElement( this->elements.at(ielem) );
            for ( int inode = 1; inode <= e->giveNumberOfNodes(); ++inode ) {
                afflictedNodes.at( e->giveNode(inode)->giveNumber() ) = 1;
            }
        }
        
        /* boundary entities are obsolete, use edges and/or surfaces instead */
        for ( int ibnd = 1; ibnd <= this->elementBoundaries.giveSize() / 2; ++ibnd ) {
            auto e = this->domain->giveElement( this->elementBoundaries.at(ibnd * 2 - 1) );
            auto boundary = this->elementBoundaries.at(ibnd * 2);
            auto fei = e->giveInterpolation();
            auto bNodes = fei->boundaryGiveNodes(boundary, e->giveGeometryType());
            for ( int inode = 1; inode <= bNodes.giveSize(); ++inode ) {
                afflictedNodes.at( e->giveNode( bNodes.at(inode) )->giveNumber() ) = 1;
            }
        }

        for ( int iedge = 1; iedge <= this->elementEdges.giveSize() / 2; ++iedge ) {
            auto e = this->domain->giveElement( this->elementEdges.at(iedge * 2 - 1) );
            auto edge = this->elementEdges.at(iedge * 2);
            auto fei = e->giveInterpolation();
            auto eNodes = fei->boundaryEdgeGiveNodes(edge, e->giveGeometryType());
            for ( int inode = 1; inode <= eNodes.giveSize(); ++inode ) {
                afflictedNodes.at( e->giveNode( eNodes.at(inode) )->giveNumber() ) = 1;
            }
        }

        for ( int isurf = 1; isurf <= this->elementSurfaces.giveSize() / 2; ++isurf ) {
            auto e = this->domain->giveElement( this->elementSurfaces.at(isurf * 2 - 1) );
            auto surf = this->elementSurfaces.at(isurf * 2);
            auto fei = e->giveInterpolation();
            auto eNodes = fei->boundarySurfaceGiveNodes(surf, e->giveGeometryType());
            for ( int inode = 1; inode <= eNodes.giveSize(); ++inode ) {
                afflictedNodes.at( e->giveNode( eNodes.at(inode) )->giveNumber() ) = 1;
            }
        }

        for ( int inode = 1; inode <= this->nodes.giveSize(); ++inode ) {
            afflictedNodes.at( this->nodes.at(inode) ) = 1;
        }
        totalNodes.findNonzeros(afflictedNodes);
    }
    return this->totalNodes;
}

const IntArray &Set :: giveSpecifiedNodeList() { return this->nodes; }

void Set :: setElementList(IntArray newElements) { this->elements = std :: move(newElements); mElementListIsSorted = false; }

void Set :: setBoundaryList(IntArray newBoundaries) { this->elementBoundaries = std :: move(newBoundaries); }

void Set :: setEdgeList(IntArray newEdges) { this->elementEdges = std :: move(newEdges); }

void Set :: setNodeList(IntArray newNodes) { this->nodes = std :: move(newNodes); }

void Set :: addAllElements()
{
    this->elements.enumerate(this->giveDomain()->giveNumberOfElements());
    mElementListIsSorted = false;
}

bool Set :: hasElement(int number) const
{
    if(!mElementListIsSorted) {
        mElementsSorted = elements;
        mElementsSorted.sort();
        mElementListIsSorted = true;
    }
    return std :: binary_search(mElementsSorted.begin(), mElementsSorted.end(), number);
}

void Set :: clear()
{
    this->elementEdges.clear();
    this->elementBoundaries.clear();
    this->elements.clear();
    this->nodes.clear();
    this->elementInternalNodes.clear();
    this->totalNodes.clear();
}

void Set :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    this->updateLocalNodeNumbering(f);
    this->updateLocalElementNumbering(f);
}

void Set :: updateLocalNodeNumbering(EntityRenumberingFunctor &f)
{
    int localNum;
    IntArray mappedNumbers;

    for ( int i = 1; i <= nodes.giveSize(); i++ ) {
        try {
            localNum = f(nodes.at(i), ERS_DofManager);
        } catch ( std :: out_of_range &e ) {
//#ifndef __MPI_PARALLEL_MODE
            OOFEM_WARNING("Set :: updateLocalNodeNumbering - Node %d not found", nodes.at(i));
//#endif
            continue;       
        }
        mappedNumbers.followedBy(localNum,10);
    }
    nodes = mappedNumbers;
}

void Set :: updateLocalElementNumbering(EntityRenumberingFunctor &f)
{
    int localNum;
    IntArray mappedNumbers;
    for ( int i = 1; i <= elements.giveSize(); i++ ) {
        try {
            localNum = f(elements.at(i), ERS_Element);
        } catch ( std :: out_of_range &e ) {
//#ifndef __MPI_PARALLEL_MODE
            OOFEM_WARNING("Set :: updateLocalElementNumbering - Element %d with indx %d not found", elements.at(i), i);
//#endif
            continue;
        }
        mappedNumbers.followedBy(localNum,10);
    }
    elements = mappedNumbers;

    ///@todo Check order of element number and boundary number (same for edges)
    mappedNumbers.resize(0);
    for ( int i = 1; i <= elementBoundaries.giveSize(); i += 2 ) {
        try {
            localNum = f(elementBoundaries.at(i), ERS_Element);
        } catch ( std :: out_of_range &e ) {
            continue;
        }
        mappedNumbers.followedBy(localNum,10);
        mappedNumbers.followedBy(elementBoundaries.at(i+1));
    }
    elementBoundaries = mappedNumbers;

    mappedNumbers.resize(0);
    for ( int i = 1; i <= elementEdges.giveSize(); i += 2 ) {
        try {
            localNum = f(elementEdges.at(i), ERS_Element);
        } catch ( std :: out_of_range &e ) {
            continue;
        }
        mappedNumbers.followedBy(localNum,10);
        mappedNumbers.followedBy(elementEdges.at(i+1));
    }
    elementEdges = mappedNumbers;
        
    for ( int i = 1; i <= this->elementInternalNodes.giveSize(); i += 2 ) {
        elementInternalNodes.at(i) = f(elementInternalNodes.at(i), ERS_Element);
    }

    mElementListIsSorted = false;
}


void Set :: saveContext(DataStream &stream, ContextMode mode)
{
    FEMComponent :: saveContext(stream, mode);

    if ( ( mode & CM_Definition ) ) {
        contextIOResultType iores;
        if ( ( iores = elements.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = elementBoundaries.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = elementEdges.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = nodes.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = this->elementInternalNodes.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }
}

void Set :: restoreContext(DataStream &stream, ContextMode mode)
{
    FEMComponent :: restoreContext(stream, mode);

    if ( mode & CM_Definition ) {
        contextIOResultType iores;
        if ( ( iores = elements.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = elementBoundaries.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = elementEdges.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = nodes.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = this->elementInternalNodes.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    this->totalNodes.clear();
}
}
