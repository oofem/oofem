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

#include "generalslavenode.h"
#include "slavedof.h"
#include "floatarray.h"
#include "intarray.h"
#include "entityrenumberingscheme.h"
#include "classfactory.h"
#include "domain.h"
#include "paramkey.h"
#include "parametermanager.h"

namespace oofem {
REGISTER_DofManager(GeneralSlaveNode);

ParamKey GeneralSlaveNode::IPK_GeneralSlaveNode_masterSizes("mastersizes");
ParamKey GeneralSlaveNode::IPK_GeneralSlaveNode_masterList("masterlist");
ParamKey GeneralSlaveNode::IPK_GeneralSlaveNode_masterWeights("masterweights");

void GeneralSlaveNode::initializeFrom(InputRecord &ir, int priority)
{
    ParameterManager &ppm =  domain->dofmanPPM;

    PM_UPDATE_TEMP_PARAMETER(IntArray, ppm, ir, this->number, IPK_DofManager_doftypemask, priority) ;
    PM_UPDATE_PARAMETER(masterSizes, ppm, ir, this->number, IPK_GeneralSlaveNode_masterSizes, priority) ;
    PM_UPDATE_TEMP_PARAMETER(IntArray, ppm, ir, this->number, IPK_GeneralSlaveNode_masterWeights, priority) ;
    PM_UPDATE_TEMP_PARAMETER(IntArray, ppm, ir, this->number, IPK_GeneralSlaveNode_masterList, priority) ;

    Node::initializeFrom(ir, priority);
}


void
GeneralSlaveNode::initializeFinish()
{
    ParameterManager &ppm =  this->giveDomain()->dofmanPPM;
    IntArray masterWeights, masterList;

    Node::initializeFinish();

    PM_ELEMENT_ERROR_IFNOTSET(ppm, this->number, IPK_DofManager_doftypemask) ;
    auto val = ppm.getTempParam(this->number, IPK_DofManager_doftypemask.getIndex());
    IntArray dofTypeMask (std::get<IntArray>(*val)); 


    int size = 0;
    for ( int i = 1; i <= dofTypeMask.giveSize(); i++ ) {
        if ( dofTypeMask.at(i) != 0 ) {
            size++;
        }
    }

    if (size > 0) {
        PM_DOFMAN_ERROR_IFNOTSET(ppm, this->number, IPK_GeneralSlaveNode_masterSizes) ;
        PM_DOFMAN_ERROR_IFNOTSET(ppm, this->number, IPK_GeneralSlaveNode_masterWeights) ;
        PM_DOFMAN_ERROR_IFNOTSET(ppm, this->number, IPK_GeneralSlaveNode_masterList) ;
        if ( masterSizes.giveSize() != size ) {
            OOFEM_ERROR("OOFEM ERROR: masterSizes size does not correspond to doftype size");
        }
        auto val = ppm.getTempParam(this->number, IPK_GeneralSlaveNode_masterWeights.getIndex());
        masterWeights = std::get<IntArray>(*val); 
        auto val2 = ppm.getTempParam(this->number, IPK_GeneralSlaveNode_masterList.getIndex());
        masterList = std::get<IntArray>(*val2); 
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
