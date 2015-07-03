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

#include "rigidarmnode.h"
#include "slavedof.h"
#include "floatarray.h"
#include "intarray.h"
#include "classfactory.h"
#include "domain.h"
#include "entityrenumberingscheme.h"

namespace oofem {
REGISTER_DofManager(RigidArmNode);

RigidArmNode :: RigidArmNode(int n, Domain *aDomain) : Node(n, aDomain)
{ }


IRResultType
RigidArmNode :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    result = Node :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    IR_GIVE_FIELD(ir, masterDofMngr, _IFT_RigidArmNode_master);

    IR_GIVE_FIELD(ir, masterMask, _IFT_DofManager_mastermask);
    if ( masterMask.giveSize() != this->dofidmask->giveSize() ) {
        OOFEM_WARNING("mastermask size mismatch");
        return IRRT_BAD_FORMAT;
    }

    return IRRT_OK;
}

void
RigidArmNode :: postInitialize()
{
    Node :: postInitialize();

    // auxiliary arrays
    std::map< DofIDItem, IntArray > masterDofID;
    std::map< DofIDItem, FloatArray > masterContribution;

    this->masterNode = dynamic_cast< Node * >( this->domain->giveDofManager(masterDofMngr) );
    if ( !masterNode ) {
        OOFEM_WARNING("master dofManager is not a node");
    }

    int masterNdofs = masterNode->giveNumberOfDofs();

    IntArray masterNodes(masterNdofs);
    for ( int &nodeNum: masterNodes ) {
        nodeNum = this->masterNode->giveNumber();
    }

    this->computeMasterContribution(masterDofID, masterContribution);

    // initialize slave dofs (inside check of consistency of receiver and master dof)
    for ( Dof *dof: *this ) {
        SlaveDof *sdof = dynamic_cast< SlaveDof * >(dof);
        if ( sdof ) {
            DofIDItem id = sdof->giveDofID();
            sdof->initialize(masterNodes, masterDofID [ id ], masterContribution [ id ]);
        }
    }

#if 0
    // check if master in same mode
    if ( parallel_mode != DofManager_local ) {
        if ( ( * masterNode )->giveParallelMode() != parallel_mode ) {
            OOFEM_WARNING("mismatch in parallel mode of RigidArmNode and master", 1);
            result = 0;
        }
    }
#endif
}


int
RigidArmNode :: checkConsistency()
{
    int result;

    result = Node :: checkConsistency();

    // check if receiver has the same coordinate system as master dofManager
    if ( !this->hasSameLCS(this->masterNode) ) {
        OOFEM_WARNING("different lcs for master/slave nodes", 1);
        result = 0;
    }

    // check if created DOFs (dofType) compatible with mastermask
    for ( int i = 1; i <= this->giveNumberOfDofs(); i++ ) {
        if ( this->masterMask.at(i) && this->dofArray [ i - 1 ]->isPrimaryDof() ) {
            OOFEM_ERROR("incompatible mastermask and doftype data");
        }
    }

    return result;
}


void
RigidArmNode :: computeMasterContribution(std::map< DofIDItem, IntArray > &masterDofID, 
                                          std::map< DofIDItem, FloatArray > &masterContribution)
{
    int k;
    IntArray R_uvw(3), uvw(3);
    FloatArray xyz(3);
    int numberOfMasterDofs = masterNode->giveNumberOfDofs();
    //IntArray countOfMasterDofs((int)masterDofID.size());

    // decode of masterMask
    uvw.at(1) = this->dofidmask->findFirstIndexOf(R_u);
    uvw.at(2) = this->dofidmask->findFirstIndexOf(R_v);
    uvw.at(3) = this->dofidmask->findFirstIndexOf(R_w);

    xyz.beDifferenceOf(*this->giveCoordinates(), *masterNode->giveCoordinates());

    if ( hasLocalCS() ) {
        // LCS is stored as global-to-local, so LCS*xyz_glob = xyz_loc
        xyz.rotatedWith(* this->localCoordinateSystem, 'n');
    }

    for ( int i = 1; i <= this->dofidmask->giveSize(); i++ ) {
        Dof *dof = this->giveDofWithID(dofidmask->at(i));
        DofIDItem id = dof->giveDofID();
        masterDofID [ id ].resize(numberOfMasterDofs);
        masterContribution [ id ].resize(numberOfMasterDofs);
        R_uvw.zero();

        switch ( masterMask.at(i) ) {
        case 0: continue;
            break;
        case 1:
            if ( id == D_u ) {
                if ( uvw.at(2) && masterMask.at( uvw.at(2) ) ) {
                    R_uvw.at(3) =  ( ( int ) R_v );
                }

                if ( uvw.at(3) && masterMask.at( uvw.at(3) ) ) {
                    R_uvw.at(2) = -( ( int ) R_w );
                }
            } else if ( id == D_v ) {
                if ( uvw.at(1) && masterMask.at( uvw.at(1) ) ) {
                    R_uvw.at(3) = -( ( int ) R_u );
                }

                if ( uvw.at(3) && masterMask.at( uvw.at(3) ) ) {
                    R_uvw.at(1) =  ( ( int ) R_w );
                }
            } else if ( id == D_w ) {
                if ( uvw.at(1) && masterMask.at( uvw.at(1) ) ) {
                    R_uvw.at(2) =  ( ( int ) R_u );
                }

                if ( uvw.at(2) && masterMask.at( uvw.at(2) ) ) {
                    R_uvw.at(1) = -( ( int ) R_v );
                }
            }

            break;
        default:
            OOFEM_ERROR("unknown value in masterMask");
        }

        //k = ++countOfMasterDofs.at(i);
        k = 1;
        masterDofID [ id ].at(k) = ( int ) id;
        masterContribution [ id ].at(k) = 1.0;

        for ( int j = 1; j <= 3; j++ ) {
            if ( R_uvw.at(j) != 0 ) {
                int sign = R_uvw.at(j) < 0 ? -1 : 1;
                //k = ++countOfMasterDofs.at(i);
                k++;
                masterDofID [ id ].at(k) = sign * R_uvw.at(j);
                masterContribution [ id ].at(k) = sign * xyz.at(j);
            }
        }
        masterDofID [ id ].resizeWithValues(k);
        masterContribution [ id ].resizeWithValues(k);
    }
}


void RigidArmNode :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    Node::updateLocalNumbering (f);
    //update masterNode numbering
    this->masterDofMngr = f( this->masterDofMngr, ERS_DofManager );
}


} // end namespace oofem
