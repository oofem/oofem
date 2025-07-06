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
 *               Copyright (C) 1993 - 2022   Borek Patzak
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

#include "generalslavenode.h"
#include "slavedof.h"
#include "floatarray.h"
#include "intarray.h"
#include "entityrenumberingscheme.h"
#include "classfactory.h"

namespace oofem {
REGISTER_DofManager(GeneralSlaveNode);

void GeneralSlaveNode::initializeFrom(InputRecord &ir)
{
    IntArray dofTypeMask;
    dofTypeMask.clear();
    IR_GIVE_FIELD(ir, dofTypeMask, _IFT_DofManager_doftypemask);
    //
    Node::initializeFrom(ir);
    //
    int size = 0;
    for ( int i = 1; i <= dofTypeMask.giveSize(); i++ ) {
        if ( dofTypeMask.at(i) != 0 ) {
            size++;
        }
    }

    IntArray masterList;
    FloatArray masterWeights;
    masterSizes = 0;
    if ( size > 0 ) {
        IR_GIVE_FIELD(ir, masterSizes, _IFT_GeneralSlaveNode_masterSizes);
        if ( masterSizes.giveSize() != size ) {
            OOFEM_ERROR("OOFEM ERROR: masterSizes size does not correspond to doftype size");
        }

        IR_GIVE_FIELD(ir, masterList, _IFT_GeneralSlaveNode_masterList);
        IR_GIVE_FIELD(ir, masterWeights, _IFT_GeneralSlaveNode_masterWeights);
    }

    int index = 0;
    for ( int j = 1; j <= masterSizes.giveSize(); j++ ) {
        IntArray dof_masterList( masterSizes.at(j) );
        IntArray dof_dofsList( masterSizes.at(j) );
        FloatArray dof_weightsList( masterSizes.at(j) );

        for ( int i = 1; i <= masterSizes.at(j); i++ ) {
            index++;
            dof_masterList.at(i) =  masterList.at(2 * index - 1);
            dof_dofsList.at(i) = masterList.at(2 * index);
            dof_weightsList.at(i) = masterWeights.at(index);
        }
        dofs_masterList.push_back(dof_masterList);
        dofs_dofsList.push_back(dof_dofsList);
        dofs_weightsList.push_back(dof_weightsList);
    }
}


void
GeneralSlaveNode::postInitialize()
{
    Node::postInitialize();

    // initialize slave dofs (inside check of consistency of receiver and master dof)
    int i = 0;
    for ( Dof *dof: * this ) {
        SlaveDof *sdof = dynamic_cast< SlaveDof * >( dof );
        if ( sdof ) {
            IntArray masterDofManagers = dofs_masterList [ i ];
            IntArray masterDofIDArray =  dofs_dofsList [ i ];
            FloatArray masterWeights = dofs_weightsList [ i ];
            sdof->initialize(masterDofManagers, masterDofIDArray, masterWeights);
            i++;
        }
    }
}



void GeneralSlaveNode::updateLocalNumbering(EntityRenumberingFunctor &f)
{
    for ( int j = 0; j < masterSizes.giveSize(); j++ ) {
        for ( int i = 1; i <= dofs_masterList [ j ].giveSize(); i++ ) {
            dofs_masterList [ j ].at(i) = f(dofs_masterList [ j ].at(i), ERS_DofManager);
        }
    }
    DofManager::updateLocalNumbering(f);
}
} // end namespace oofem
