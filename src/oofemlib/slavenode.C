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

#include "slavenode.h"
#include "slavedof.h"
#include "floatarray.h"
#include "intarray.h"
#include "entityrenumberingscheme.h"
#include "classfactory.h"

namespace oofem {
REGISTER_DofManager(SlaveNode);

IRResultType SlaveNode :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, masterDofManagers, _IFT_SlaveNode_masterDofManagers);
    IR_GIVE_OPTIONAL_FIELD(ir, masterWeights, _IFT_SlaveNode_weights);

    if ( masterWeights.giveSize() == 0 ) {
        masterWeights.resize( masterDofManagers.giveSize() );
        masterWeights.add( 1 / ( double ) masterDofManagers.giveSize() );
    } else if ( masterDofManagers.giveSize() != masterWeights.giveSize() ) {
        OOFEM_WARNING("master dof managers and weights size mismatch.");
        return IRRT_BAD_FORMAT;
    }
    return Node :: initializeFrom(ir);
}


void SlaveNode :: postInitialize()
{
    Node :: postInitialize();

    // initialize slave dofs (inside check of consistency of receiver and master dof)
    for ( Dof *dof: *this ) {
        SlaveDof *sdof = dynamic_cast< SlaveDof * >(dof);
        if ( sdof ) {
            sdof->initialize(masterDofManagers, IntArray(), masterWeights);
        }
    }
}


void SlaveNode :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    for ( int i = 1; i <= masterDofManagers.giveSize(); ++i ) {
        masterDofManagers.at(i) = f(masterDofManagers.at(i), ERS_DofManager);
    }

    DofManager :: updateLocalNumbering(f);
}
} // end namespace oofem
