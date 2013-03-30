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

#include "rigidarmnode.h"
#include "slavedof.h"
#include "flotarry.h"
#include "intarray.h"

namespace oofem {
RigidArmNode :: RigidArmNode(int n, Domain *aDomain) : Node(n, aDomain)
{ }


void
RigidArmNode :: allocAuxArrays()
{
    countOfMasterDofs = new IntArray(numberOfDofs);

    masterDofID = new IntArray * [ numberOfDofs ];
    masterContribution = new FloatArray * [ numberOfDofs ];

    for ( int i = 1; i <= numberOfDofs; i++ ) {
        if ( dofArray [ i - 1 ]->isPrimaryDof() ) {
            masterDofID [ i - 1 ] = NULL;
            masterContribution [ i - 1 ] = NULL;
        } else {
            masterDofID [ i - 1 ] = new IntArray;
            masterContribution [ i - 1 ] = new FloatArray;
        }
    }
}

void
RigidArmNode :: deallocAuxArrays()
{
    delete countOfMasterDofs;

    for ( int i = 0; i < numberOfDofs; i++ ) {
        delete masterDofID [ i ];
        delete masterContribution [ i ];
    }

    delete[] masterDofID;
    delete[] masterContribution;
}


IRResultType
RigidArmNode :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    Node :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, masterDofMngr, _IFT_RigidArmNode_master);

    IR_GIVE_FIELD(ir, masterMask, _IFT_DofManager_mastermask); ///@todo Is this really necessary, dofmanager should have read this already.
    if ( masterMask.giveSize() != this->giveNumberOfDofs() ) {
        _error("initializeFrom: mastermask size mismatch");
    }

    return IRRT_OK;
}

void
RigidArmNode :: postInitialize()
{
    Node :: postInitialize();
    
    // allocate auxiliary arrays
    allocAuxArrays();
    
    this->masterNode = dynamic_cast< Node * >( this->domain->giveDofManager(masterDofMngr) );
    if ( !masterNode ) {
        OOFEM_WARNING("RigidArmNode :: postInitialize: master dofManager is not a node");
    }
    
    int masterNdofs = masterNode->giveNumberOfDofs();

    IntArray masterNodes(masterNdofs);
    for ( int i = 1; i <= masterNdofs; i++ ) {
        masterNodes.at(i) = this->masterNode->giveNumber();
    }

    for ( int i = 1; i <= numberOfDofs; i++ ) {
        if ( masterDofID [ i - 1 ] ) {
            countOfMasterDofs->at(i) = 0;
            masterDofID [ i - 1 ]->resize(masterNdofs);
            masterContribution [ i - 1 ]->resize(masterNdofs);
        }
    }

    this->computeMasterContribution();

    // initialize slave dofs (inside check of consistency of receiver and master dof)
    for ( int i = 1; i <= numberOfDofs; i++ ) {
        SlaveDof *sdof = dynamic_cast< SlaveDof* >( dofArray [ i - 1 ] );
        if ( sdof ) {
            sdof->initialize(countOfMasterDofs->at(i), masterNodes, masterDofID [ i - 1 ], *masterContribution [ i - 1 ]);
        }
    }

#if 0
#ifdef __PARALLEL_MODE
    // check if master in same mode
    if ( parallel_mode != DofManager_local ) {
        if ( ( * masterNode )->giveParallelMode() != parallel_mode ) {
        _warning2("checkConsistency: mismatch in parallel mode of RigidArmNode and master", 1);
        result = 0;
        }
    }
#endif
#endif

    // deallocate auxiliary arrays
    deallocAuxArrays();

}


int
RigidArmNode :: checkConsistency()
{
    int result;

    result = Node :: checkConsistency();

    // check if receiver has the same coordinate system as master dofManager
    if ( !this->hasSameLCS(this->masterNode) ) {
        _warning2("checkConsistency: different lcs for master/slave nodes", 1);
        result = 0;
    }

    // check if created DOFs (dofType) compatible with mastermask
    for ( int i = 1; i <= this->giveNumberOfDofs(); i++ ) {
        if ( this->masterMask.at(i) && this->dofArray [ i - 1 ]->isPrimaryDof() ) {
            _error("checkConsistency: incompatible mastermask and doftype data");
        }
    }

    return result;
}


void
RigidArmNode :: computeMasterContribution()
{
    int k, sign;
    IntArray R_uvw(3), uvw(3);
    FloatArray xyz(3);
    DofIDItem id;

    // decode of masterMask
    uvw.at(1) = this->findDofWithDofId(R_u);
    uvw.at(2) = this->findDofWithDofId(R_v);
    uvw.at(3) = this->findDofWithDofId(R_w);

    for ( int i = 1; i <= 3; i++ ) {
        xyz.at(i) = this->giveCoordinate(i) - masterNode->giveCoordinate(i);
    }

    if ( hasLocalCS() ) {
        // LCS is stored as global-to-local, so LCS*xyz_glob = xyz_loc
        xyz.rotatedWith(* this->localCoordinateSystem, 'n');
    }

    for ( int i = 1; i <= numberOfDofs; i++ ) {
        id = this->giveDof(i)->giveDofID();
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
            _error("computeMasterContribution: unknown value in masterMask");
        }

        k = ++ this->countOfMasterDofs->at(i);
        this->masterDofID [ i - 1 ]->at(k) = ( int ) id;
        this->masterContribution [ i - 1 ]->at(k) = 1.0;

        for ( int j = 1; j <= 3; j++ ) {
            if ( R_uvw.at(j) != 0 ) {
                sign = R_uvw.at(j) < 0 ? -1 : 1;

                k = ++ this->countOfMasterDofs->at(i);
                this->masterDofID [ i - 1 ]->at(k) = sign * R_uvw.at(j);
                this->masterContribution [ i - 1 ]->at(k) = sign * xyz.at(j);
            }
        }
    }
}
} // end namespace oofem
