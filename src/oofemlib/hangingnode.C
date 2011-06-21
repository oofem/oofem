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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#include "hangingnode.h"
#include "slavedof.h"
#include "flotarry.h"
#include "intarray.h"
#include "element.h"
#include "feinterpol.h"

namespace oofem {
HangingNode :: HangingNode(int n, Domain *aDomain) : Node(n, aDomain)
{ 
#ifdef __OOFEG
	consistencyChecked = false;
#endif
}

IRResultType HangingNode :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    Node :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, masterElement, IFT_HangingNode_masterElement, "masterelement");
    return IRRT_OK;
}


int HangingNode :: checkConsistency()
{
    Element *e;
    FEInterpolation *fei;
    FloatArray lcoords, masterContribution;
    bool result = true;

#ifdef __OOFEG
    if (consistencyChecked) return result;
    consistencyChecked = true;
#endif

    result = result && Node :: checkConsistency();

    // First check element and interpolation
    if ( !(e = this->giveDomain()->giveElement(this->masterElement)) ) {
        OOFEM_WARNING2("HangingNode :: checkConsistency - Requested element %d doesn't exist.", this->masterElement);
        return false;
    }
    if ( !(fei = e->giveInterpolation()) ) {
        OOFEM_WARNING2("HangingNode :: checkConsistency - Requested element %d doesn't have a interpolator.", this->masterElement);
        return false;
    }
#if 0
#ifdef __PARALLEL_MODE
    // Check if master is in same mode
    if ( parallel_mode != DofManager_local ) {
        for ( int i = 1; i <= countOfMasterNodes; i++ ) {
            if ( e->giveNode(i)->giveParallelMode() != parallel_mode ) {
                OOFEM_WARNING("HangingNode :: checkConsistency - Mismatch in parallel mode of HangingNode and master");
                return false;
            }
        }
    }
#endif
#endif


    fei->global2local(lcoords, coordinates, FEIElementGeometryWrapper(e), 0.0);
    fei->evalN(masterContribution, lcoords, FEIElementGeometryWrapper(e), 0.0);

    // Initialize slave dofs (inside check of consistency of receiver and master dof)
    //AList<Node*> masterNode(e->giveNumberOfNodes()); // TODO: Change to AList, avoiding pointers.
    Node **masterNode = new Node * [ e->giveNumberOfNodes() ];
    for ( int i = 1; i <= e->giveNumberOfNodes(); ++i ) {
        masterNode[i-1] = e->giveNode(i);
        if ( !this->hasSameLCS(e->giveNode(i)) ) {
            OOFEM_WARNING("HangingNode :: checkConsistency - Different lcs for master/slave nodes.");
            result = false;
        }
    }
    for ( int i = 1; i <= numberOfDofs; i++ ) {
        if ( dofArray [ i - 1 ]->giveClassID() == SlaveDofClass ) {
            ( ( SlaveDof * ) dofArray [ i - 1 ] )->initialize(masterContribution.giveSize(), masterNode, NULL, masterContribution);
        }
    }
    delete masterNode;

    return true;
}

} // end namespace oofem
