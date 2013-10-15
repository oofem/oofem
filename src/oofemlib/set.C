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

#include <list>

namespace oofem {

IRResultType Set :: initializeFrom(InputRecord* ir)
{
    const char *__proc = "initializeFrom";
    IRResultType result;

    FEMComponent :: initializeFrom(ir);

    IntArray inputNodes;
    std::list< Range > inputNodeRanges;
    IR_GIVE_OPTIONAL_FIELD(ir, inputNodes, _IFT_Set_nodes);
    IR_GIVE_OPTIONAL_FIELD(ir, inputNodeRanges, _IFT_Set_nodeRanges);
    this->computeIntArray(this->nodes, inputNodes, inputNodeRanges);

    IntArray inputElements;
    std::list< Range > inputElementRanges;
    IR_GIVE_OPTIONAL_FIELD(ir, inputElements, _IFT_Set_elements);
    IR_GIVE_OPTIONAL_FIELD(ir, inputElementRanges, _IFT_Set_elementRanges);
    this->computeIntArray(this->elements, inputElements, inputElementRanges);

    this->elementBoundaries.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, this->elementBoundaries, _IFT_Set_elementBoundaries);
    
    this->elementEdges.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, this->elementEdges, _IFT_Set_elementEdges);

    return IRRT_OK;
}

void Set :: computeIntArray(IntArray& answer, const IntArray& specified, std::list< Range > ranges)
{
    // Find the max value;
    int maxIndex = specified.giveSize() == 0 ? 0 : specified.maximum();
    for (std::list< Range >::iterator it = ranges.begin(); it != ranges.end(); ++it) {
        if ( it->giveEnd() > maxIndex ) maxIndex = it->giveEnd();
    }
    IntArray afflictedNodes(maxIndex);
    afflictedNodes.zero();

    for (int i = 1; i <= specified.giveSize(); ++i) {
        afflictedNodes.at(specified.at(i)) = 1;
    }
    
    for (std::list< Range >::iterator it = ranges.begin(); it != ranges.end(); ++it) {
        for (int i = it->giveStart(); i <= it->giveEnd(); ++i) {
            afflictedNodes.at(i) = 1;
        }
    }
    answer.findNonzeros(afflictedNodes);
}


const IntArray& Set :: giveElementList() { return elements; }

const IntArray& Set :: giveBoundaryList() { return elementBoundaries; }

const IntArray& Set :: giveEdgeList() { return elementEdges; }

const IntArray& Set :: giveNodeList()
{
    // Lazy evaluation, we compute the unique set of nodes if needed (and store it).
    if ( this->totalNodes.giveSize() == 0 ) {
        IntArray afflictedNodes(this->domain->giveNumberOfDofManagers());
        afflictedNodes.zero();
        for (int ielem = 1; ielem <= this->elements.giveSize(); ++ielem) {
            Element *e = this->domain->giveElement(this->elements.at(ielem));
            for (int inode = 1; inode <= e->giveNumberOfNodes(); ++inode) {
                afflictedNodes.at(e->giveNode(inode)->giveNumber()) = 1;
            }
        }

        IntArray bNodes;
        for (int ibnd = 1; ibnd <= this->elementBoundaries.giveSize()/2; ++ibnd) {
            Element *e = this->domain->giveElement( this->elementBoundaries.at(ibnd*2-1) );
            int boundary = this->elementBoundaries.at(ibnd*2);
            FEInterpolation *fei = e->giveInterpolation();
            fei->boundaryGiveNodes(bNodes, boundary);
            for (int inode = 1; inode <= bNodes.giveSize(); ++inode) {
                afflictedNodes.at(e->giveNode(bNodes.at(inode))->giveNumber()) = 1;
            }
        }

        IntArray eNodes;
        for (int iedge = 1; iedge <= this->elementEdges.giveSize()/2; ++iedge) {
            Element *e = this->domain->giveElement( this->elementEdges.at(iedge*2-1) );
            int edge = this->elementEdges.at(iedge*2);
            FEInterpolation3d *fei = static_cast< FEInterpolation3d* >( e->giveInterpolation() );
            fei->computeLocalEdgeMapping(eNodes, edge);
            for (int inode = 1; inode <= eNodes.giveSize(); ++inode) {
                afflictedNodes.at(e->giveNode(eNodes.at(inode))->giveNumber()) = 1;
            }
        }

        for (int inode = 1; inode <= this->nodes.giveSize(); ++inode) {
            afflictedNodes.at( this->nodes.at(inode) ) = 1;
        }
        totalNodes.findNonzeros(afflictedNodes);
    }
    return this->totalNodes;
}

void Set :: setElementList(const IntArray& newElements) { this->elements = newElements; }

void Set :: setBoundaryList(const IntArray& newBoundaries) { this->elementBoundaries = newBoundaries; }

void Set :: setEdgeList(const IntArray& newEdges) { this->elementEdges = newEdges; }

void Set :: setNodeList(const IntArray& newNodes) { this->nodes = newNodes; }

void Set :: clear()
{
    this->elementEdges.resize(0);
    this->elementBoundaries.resize(0);
    this->elements.resize(0);
    this->nodes.resize(0);
    this->totalNodes.resize(0);
}

void Set :: updateLocalNumbering(EntityRenumberingFunctor& f)
{
    this->updateLocalNodeNumbering(f);
    this->updateLocalElementNumbering(f);
}

void Set :: updateLocalNodeNumbering(EntityRenumberingFunctor& f)
{
    for (int i = 1; i <= nodes.giveSize(); i++ ) {
        nodes.at(i) = f(nodes.at(i), ERS_DofManager);
    }
}

void Set :: updateLocalElementNumbering(EntityRenumberingFunctor& f)
{
    for (int i = 1; i <= elements.giveSize(); i++ ) {
        elements.at(i) = f(elements.at(i), ERS_Element);
    }
    ///@todo Check order of element number and boundary number (same for edges)
    for (int i = 1; i <= elementBoundaries.giveSize(); i+=2 ) {
        elementBoundaries.at(i) = f(elementBoundaries.at(i), ERS_Element);
    }
    for (int i = 1; i <= elementEdges.giveSize(); i+=2 ) {
        elementEdges.at(i) = f(elementEdges.at(i), ERS_Element);
    }
}


contextIOResultType Set :: saveContext(DataStream* stream, ContextMode mode, void* obj)
{
    contextIOResultType iores;

    if ( ( iores = FEMComponent :: saveContext(stream, mode, obj) ) != CIO_OK ) { THROW_CIOERR(iores); }

    if ( ( mode & CM_Definition ) ) {
        if ( ( iores = elements.storeYourself(stream, mode) ) != CIO_OK ) { THROW_CIOERR(iores); }
        if ( ( iores = elementBoundaries.storeYourself(stream, mode) ) != CIO_OK ) { THROW_CIOERR(iores); }
        if ( ( iores = elementEdges.storeYourself(stream, mode) ) != CIO_OK ) { THROW_CIOERR(iores); }
        if ( ( iores = nodes.storeYourself(stream, mode) ) != CIO_OK ) { THROW_CIOERR(iores); }
    }

    return CIO_OK;
}

contextIOResultType Set :: restoreContext(DataStream* stream, ContextMode mode, void* obj)
{
    contextIOResultType iores;

    if ( ( iores = FEMComponent :: restoreContext(stream, mode, obj) ) != CIO_OK ) { THROW_CIOERR(iores); }

    if ( mode & CM_Definition ) {
        if ( ( iores = elements.restoreYourself(stream, mode) ) != CIO_OK ) { THROW_CIOERR(iores); }
        if ( ( iores = elementBoundaries.restoreYourself(stream, mode) ) != CIO_OK ) { THROW_CIOERR(iores); }
        if ( ( iores = elementEdges.restoreYourself(stream, mode) ) != CIO_OK ) { THROW_CIOERR(iores); }
        if ( ( iores = nodes.restoreYourself(stream, mode) ) != CIO_OK ) { THROW_CIOERR(iores); }
    }

    this->totalNodes.resize(0);

    return CIO_OK;
}

}
