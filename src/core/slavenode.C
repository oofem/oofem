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

#include "slavenode.h"
#include "slavedof.h"
#include "floatarray.h"
#include "intarray.h"
#include "entityrenumberingscheme.h"
#include "classfactory.h"
#include "domain.h"
#include "paramkey.h"
#include "parametermanager.h"

namespace oofem {
REGISTER_DofManager(SlaveNode);

ParamKey SlaveNode::IPK_SlaveNode_masterDofManagers("masterdofman");
ParamKey SlaveNode::IPK_SlaveNode_weights("weights");

void SlaveNode :: initializeFrom(InputRecord &ir, int priority)
{
    ParameterManager &ppm =  this->giveDomain()->dofmanPPM;

    Node :: initializeFrom(ir, priority);

    PM_UPDATE_PARAMETER(masterDofManagers, ppm, ir, this->number, IPK_SlaveNode_masterDofManagers, priority) ;
    PM_UPDATE_PARAMETER(masterWeights, ppm, ir, this->number, IPK_SlaveNode_weights, priority) ;
}


void SlaveNode :: postInitialize()
{
    ParameterManager &ppm =  this->giveDomain()->dofmanPPM;

    Node :: postInitialize();
    PM_ELEMENT_ERROR_IFNOTSET(ppm, this->number, IPK_SlaveNode_masterDofManagers) ;

    if ( masterWeights.giveSize() == 0 ) {
        masterWeights.resize( masterDofManagers.giveSize() );
        masterWeights.add( 1 / ( double ) masterDofManagers.giveSize() );
    } else if ( masterDofManagers.giveSize() != masterWeights.giveSize() ) {
        throw ComponentInputException(IPK_SlaveNode_weights.getName(), ComponentInputException::ComponentType::ctDofManager, this->number, "master dof managers and weights size mismatch.");
    }

    // initialize slave dofs (inside check of consistency of receiver and master dof)
    for ( Dof *dof: *this ) {
        SlaveDof *sdof = dynamic_cast< SlaveDof * >(dof);
        if ( sdof ) {
            sdof->initialize(masterDofManagers, IntArray(), masterWeights);
        }
    }
    // clean up
    this->masterWeights.clear();
}


void SlaveNode :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    for ( int i = 1; i <= masterDofManagers.giveSize(); ++i ) {
        masterDofManagers.at(i) = f(masterDofManagers.at(i), ERS_DofManager);
    }

    DofManager :: updateLocalNumbering(f);
}
} // end namespace oofem
