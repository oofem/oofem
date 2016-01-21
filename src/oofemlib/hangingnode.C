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

#include "hangingnode.h"
#include "slavedof.h"
#include "floatarray.h"
#include "intarray.h"
#include "element.h"
#include "feinterpol.h"
#include "spatiallocalizer.h"
#include "classfactory.h"

namespace oofem {
REGISTER_DofManager(HangingNode);

HangingNode :: HangingNode(int n, Domain *aDomain) : Node(n, aDomain)
{
#ifdef __OOFEG
    initialized = false;
#endif
}

IRResultType HangingNode :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    this->masterElement = -1;
    IR_GIVE_OPTIONAL_FIELD(ir, this->masterElement, _IFT_HangingNode_masterElement);
    this->masterRegion = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, this->masterRegion, _IFT_HangingNode_masterRegion);
    return Node :: initializeFrom(ir);
}

int HangingNode :: checkConsistency()
{
    Element *e = this->domain->giveElement(this->masterElement);
    int result = Node :: checkConsistency();

#if 0
    // Check if master is in same mode
    if ( parallel_mode != DofManager_local ) {
        for ( int i = 1; i <= countOfMasterNodes; i++ ) {
            if ( e->giveNode(i)->giveParallelMode() != parallel_mode ) {
                OOFEM_WARNING("Mismatch in parallel mode of HangingNode and master");
                return false;
            }
        }
    }
#endif

    // Check local coordinate systems
    for ( int i = 1; i <= e->giveNumberOfNodes(); ++i ) {
        if ( !this->hasSameLCS( e->giveNode(i) ) ) {
            OOFEM_WARNING("Different lcs for master/slave nodes.");
            result = false;
        }
    }
    return result;
}

void HangingNode :: postInitialize()
{
    Node :: postInitialize();

    Element *e;
    FEInterpolation *fei;
    FloatArray lcoords, masterContribution;

#ifdef __OOFEG
    if ( initialized ) {
        return;
    }
    initialized = true;
#endif

    // First check element and interpolation
    if ( masterElement == -1 ) { // Then we find it by taking the closest (probably containing element)
        FloatArray closest;
        SpatialLocalizer *sp = this->domain->giveSpatialLocalizer();
        sp->init();
        // Closest point or containing point? It should be contained, but with numerical errors it might be slightly outside
        // so the closest point is more robust.
        if ( !( e = sp->giveElementClosestToPoint(lcoords, closest, coordinates, this->masterRegion) ) ) {
            OOFEM_ERROR("Couldn't find closest element (automatically).");
        }
        this->masterElement = e->giveNumber();
    } else if ( !( e = this->giveDomain()->giveElement(this->masterElement) ) ) {
        OOFEM_ERROR("Requested element %d doesn't exist.", this->masterElement);
    }
    if ( !( fei = e->giveInterpolation() ) ) {
        OOFEM_ERROR("Requested element %d doesn't have a interpolator.", this->masterElement);
    }

    if ( lcoords.giveSize() == 0 ) { // we don't need to do this again if the spatial localizer was used.
        fei->global2local( lcoords, coordinates, FEIElementGeometryWrapper(e) );
    }

    // Initialize slave dofs (inside check of consistency of receiver and master dof)
    const IntArray &masterNodes = e->giveDofManArray();
    for ( Dof *dof: *this ) {
        SlaveDof *sdof = dynamic_cast< SlaveDof * >(dof);
        if ( sdof ) {
            DofIDItem id = sdof->giveDofID();
            fei = e->giveInterpolation(id);
            if ( !fei ) {
                OOFEM_ERROR("Requested interpolation for dof id %d doesn't exist in element %d.",
                             id, this->masterElement);
            }
#if 0 // This won't work (yet), as it requires some more general FEI classes, or something similar.
            if ( fei->hasMultiField() ) {
                FloatMatrix multiContribution;
                IntArray masterDofIDs, masterNodesDup, dofids;
                fei->evalMultiN(multiContribution, dofids, lcoords, FEIElementGeometryWrapper(e), 0.0);
                masterContribution.flatten(multiContribution);
                masterDofIDs.clear();
                for ( int i = 0; i <= multiContribution.giveNumberOfColumns(); ++i ) {
                    masterDofIDs.followedBy(dofids);
                    masterNodesDup.followedBy(masterNodes);
                }
                sdof->initialize(masterNodesDup, & masterDofIDs, masterContribution);
            } else { }
#else
            // Note: There can be more masterNodes than masterContributions, since all the
            // FEI classes are based on that the first nodes correspond to the simpler/linear interpolation.
            // If this assumption is changed in FEIElementGeometryWrapper + friends,
            // masterNode will also need to be modified for each dof accordingly.
            fei->evalN( masterContribution, lcoords, FEIElementGeometryWrapper(e) );
            sdof->initialize(masterNodes, IntArray(), masterContribution);
#endif
        }
    }
}
} // end namespace oofem
