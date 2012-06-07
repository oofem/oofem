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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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
    // Check local coordinate systems
    for ( int i = 1; i <= e->giveNumberOfNodes(); ++i ) {
        if ( !this->hasSameLCS(e->giveNode(i)) ) {
            OOFEM_WARNING("HangingNode :: checkConsistency - Different lcs for master/slave nodes.");
            result = false;
        }
    }

    fei->global2local(lcoords, coordinates, FEIElementGeometryWrapper(e));

    // Initialize slave dofs (inside check of consistency of receiver and master dof)
    const IntArray &masterNodes = e->giveDofManArray();
    for ( int i = 1; i <= numberOfDofs; i++ ) {
        if ( dofArray [ i - 1 ]->giveClassID() == SlaveDofClass ) {
            SlaveDof *sdof = ( SlaveDof * ) dofArray [ i - 1 ];
            DofIDItem id = sdof->giveDofID();
            fei = e->giveInterpolation(id);
            if (!fei) {
                OOFEM_WARNING3("HangingNode :: checkConsistency - Requested interpolation for dof id %d doesn't exist in element %d.",
                        id, this->masterElement);
                return false;
            }
#if 0// This won't work (yet), as it requires some more general FEI classes, or something similar.
            if (fei->hasMultiField()) {
                FloatMatrix multiContribution;
                IntArray masterDofIDs, masterNodesDup, dofids;
                fei->evalMultiN(multiContribution, dofids, lcoords, FEIElementGeometryWrapper(e), 0.0);
                masterContribution.flatten(multiContribution);
                masterDofIDs.resize(0);
                for (int i = 0; i <= multiContribution.giveNumberOfColumns(); ++i) {
                    masterDofIDs.followedBy(dofids);
                    masterNodesDup.followedBy(masterNodes);
                }
                sdof->initialize(masterContribution.giveSize(), masterNodesDup, &masterDofIDs, masterContribution);
            } else { }
#else
            // Note: There can be more masterNodes than masterContributions, since all the
            // FEI classes are based on that the first nodes correspond to the simpler/linear interpolation.
            // If this assumption is changed in FEIElementGeometryWrapper + friends,
            // masterNode will also need to be modified for each dof accordingly.
            fei->evalN(masterContribution, lcoords, FEIElementGeometryWrapper(e));
            sdof->initialize(masterContribution.giveSize(), masterNodes, NULL, masterContribution);
#endif
        }
    }

    return true;
}

} // end namespace oofem
