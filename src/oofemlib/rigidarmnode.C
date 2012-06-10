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
    masterMask = new IntArray( this->giveNumberOfDofs() );
    countOfMasterDofs = new IntArray(numberOfDofs);

    masterDofID = new IntArray * [ numberOfDofs ];
    masterContribution = new FloatArray * [ numberOfDofs ];

    for ( int i = 1; i <= numberOfDofs; i++ ) {
        if ( dofArray [ i - 1 ]->giveClassID() == SlaveDofClass ) {
            masterDofID [ i - 1 ] = new IntArray;
            masterContribution [ i - 1 ] = new FloatArray;
        } else {
            masterDofID [ i - 1 ] = NULL;
            masterContribution [ i - 1 ] = NULL;
        }
    }
}

void
RigidArmNode :: deallocAuxArrays()
{
    delete masterMask;
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

    // allocate auxiliary arrays
    allocAuxArrays();

    IR_GIVE_FIELD(ir, masterDofMngr, IFT_RigidArmNode_master, "master"); // Macro

    IR_GIVE_FIELD(ir, * masterMask, IFT_DofManager_mastermask, "mastermask"); // Macro
    if ( masterMask->giveSize() != this->giveNumberOfDofs() ) {
        _error("initializeFrom: mastermask size mismatch");
    }

    return IRRT_OK;
}


int
RigidArmNode :: checkConsistency()
// Checks internal data consistency in node.
// Current implementation checks (when receiver has slave dofs) if receiver has the same
// coordinate system as master dofManager of slave dof.
// If requested, computes natural coordinates on element-like masternodes

{
    int result = 1;
    int i, ndofs;
    Node *master;

    result = result && Node :: checkConsistency();

    // finds master node
    master = dynamic_cast< Node * >( this->domain->giveDofManager(masterDofMngr) );
    if ( !master ) {
        _warning2("checkConsistency: master dofManager is not compatible", 1);
        result = 0;
    }

    // check if receiver has the same coordinate system as master dofManager
    if ( !this->hasSameLCS(master) ) {
        _warning2("checkConsistency: different lcs for master/slave nodes", 1);
        result = 0;
    }

    // check if created DOFs (dofType) compatible with mastermask
    ndofs = master->giveNumberOfDofs();
    for ( i = 1; i <= numberOfDofs; i++ ) {
        if ( masterMask->at(i) && ( dofArray [ i - 1 ]->giveClassID() == MasterDofClass ) ) {
            _error("checkConsistency: incompatible mastermask and doftype data");
        }
    }


    // allocate
    for ( i = 1; i <= numberOfDofs; i++ ) {
        if ( masterDofID [ i - 1 ] ) {
            countOfMasterDofs->at(i) = 0;
            masterDofID [ i - 1 ]->resize(ndofs);
            masterContribution [ i - 1 ]->resize(ndofs);
        }
    }

    IntArray masterNodes(ndofs);
    masterNode = master;
    for ( i = 1; i <= ndofs; i++ ) {
        masterNodes.at(i) = master->giveNumber();
    }

    result = result && computeMasterContribution();

    // initialize slave dofs (inside check of consistency of receiver and master dof)
    for ( i = 1; i <= numberOfDofs; i++ ) {
        if ( dofArray [ i - 1 ]->giveClassID() == SlaveDofClass ) {
            ( ( SlaveDof * ) dofArray [ i - 1 ] )->initialize(countOfMasterDofs->at(i), masterNodes, masterDofID [ i - 1 ], *masterContribution [ i - 1 ]);
        }
    }

    /*
     * #ifdef __PARALLEL_MODE
     *  // check if master in same mode
     *  if ( parallel_mode != DofManager_local ) {
     *      if ( ( * masterNode )->giveParallelMode() != parallel_mode ) {
     *          _warning2("checkConsistency: mismatch in parallel mode of RigidArmNode and master", 1);
     *          result = 0;
     *      }
     *  }
     *
     * #endif
     */

    // deallocate auxiliary arrays
    deallocAuxArrays();

    return result;
}


int
RigidArmNode :: computeMasterContribution()
{
    int i, j, k, sign;
    IntArray R_uvw(3), uvw(3);
    FloatArray xyz(3);
    DofIDItem id;

    // decode of masterMask
    uvw.at(1) = this->findDofWithDofId(R_u);
    uvw.at(2) = this->findDofWithDofId(R_v);
    uvw.at(3) = this->findDofWithDofId(R_w);

    for ( i = 1; i <= 3; i++ ) {
        xyz.at(i) = this->giveCoordinate(i) - masterNode->giveCoordinate(i);
    }

    if ( hasLocalCS() ) {
        // LCS is stored as global-to-local, so LCS*xyz_glob = xyz_loc
        xyz.rotatedWith(* this->localCoordinateSystem, 'n');
    }

    for ( i = 1; i <= numberOfDofs; i++ ) {
        id = this->giveDof(i)->giveDofID();
        R_uvw.zero();

        switch ( masterMask->at(i) ) {
        case 0: continue;
            break;
        case 1:
            if ( id == D_u ) {
                if ( uvw.at(2) && masterMask->at( uvw.at(2) ) ) {
                    R_uvw.at(3) =  ( ( int ) R_v );
                }

                if ( uvw.at(3) && masterMask->at( uvw.at(3) ) ) {
                    R_uvw.at(2) = -( ( int ) R_w );
                }
            } else if ( id == D_v ) {
                if ( uvw.at(1) && masterMask->at( uvw.at(1) ) ) {
                    R_uvw.at(3) = -( ( int ) R_u );
                }

                if ( uvw.at(3) && masterMask->at( uvw.at(3) ) ) {
                    R_uvw.at(1) =  ( ( int ) R_w );
                }
            } else if ( id == D_w ) {
                if ( uvw.at(1) && masterMask->at( uvw.at(1) ) ) {
                    R_uvw.at(2) =  ( ( int ) R_u );
                }

                if ( uvw.at(2) && masterMask->at( uvw.at(2) ) ) {
                    R_uvw.at(1) = -( ( int ) R_v );
                }
            }

            break;
        default:
            _warning("computeMasterContribution: unknown value in masterMask");
            return 0;
        }

        k = ++countOfMasterDofs->at(i);
        masterDofID [ i - 1 ]->at(k) = ( int ) id;
        masterContribution [ i - 1 ]->at(k) = 1.0;

        for ( j = 1; j <= 3; j++ ) {
            if ( R_uvw.at(j) != 0 ) {
                sign = R_uvw.at(j) < 0 ? -1 : 1;

                k = ++countOfMasterDofs->at(i);
                masterDofID [ i - 1 ]->at(k) = sign * R_uvw.at(j);
                masterContribution [ i - 1 ]->at(k) = sign * xyz.at(j);
            }
        }
    }

    return 1;
}
} // end namespace oofem
