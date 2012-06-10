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

#include "slavenode.h"
#include "slavedof.h"
#include "flotarry.h"
#include "intarray.h"

namespace oofem {
IRResultType SlaveNode :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    Node :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, masterDofManagers, IFT_SlaveNode_masterDofManagers, "masterdofman");
    IR_GIVE_OPTIONAL_FIELD(ir, masterWeights, IFT_SlaveNode_weights, "weights");

    if ( masterWeights.giveSize() == 0 ) {
        masterWeights.resize(masterDofManagers.giveSize());
        masterWeights.add(1/(double)masterDofManagers.giveSize());
    } else if ( masterDofManagers.giveSize() != masterWeights.giveSize() ) {
        _error("initializeFrom: master dof managers and weights size mismatch.");
    }
    return IRRT_OK;
}


int SlaveNode :: checkConsistency()
{
    bool result = true;

    result = result && Node :: checkConsistency();

    // initialize slave dofs (inside check of consistency of receiver and master dof)
    for ( int i = 1; i <= numberOfDofs; ++i ) {
        if ( dofArray[i-1]->giveClassID() == SlaveDofClass ) {
            ((SlaveDof*)dofArray[i-1])->initialize(masterDofManagers.giveSize(), masterDofManagers, NULL, masterWeights);
        }
    }
    return result;
}


void SlaveNode :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    for ( int i = 1; i <= masterDofManagers.giveSize(); ++i) {
        masterDofManagers.at(i) = f(masterDofManagers.at(i), ERS_DofManager);
    }

    DofManager::updateLocalNumbering(f);
}

} // end namespace oofem
