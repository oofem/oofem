/* $Header: /home/cvs/bp/oofem/oofemlib/src/dofmanager.C,v 1.18.4.1 2004/04/05 15:19:43 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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
#include "timestep.h"
#include "masterdof.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "intarray.h"

#ifndef __MAKEDEPEND
#include <math.h>
#include <stdlib.h>
#endif

#include "fei1dlin.h"
#include "fei2dtrlin.h"
#include "fei2dquadlin.h"
#include "fei3dhexalin.h"
#include "mathfem.h"



/**
 * Constructor. Creates a hanging node with number n, belonging to aDomain.
 */
HangingNode :: HangingNode(int n, Domain *aDomain) : Node(n, aDomain)
{ }


void
HangingNode :: allocAuxArrays(void)
{
    locoords = new FloatArray(3);
    masterDofMngr = new IntArray;
    masterContribution = new FloatArray;
}

void
HangingNode :: deallocAuxArrays(void)
{
    delete locoords;
    delete masterDofMngr;
    delete masterContribution;

    // only delete
    delete masterNode;
}


/**
 * Gets from the source line from the data file all the data of the receiver.
 */
IRResultType
HangingNode :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    Node :: initializeFrom(ir);

    // allocate auxiliary arrays
    allocAuxArrays();

    IR_GIVE_FIELD(ir, * masterDofMngr, IFT_HangingNode_masters, "masters"); // Macro
    countOfMasterNodes = masterDofMngr->giveSize();

    IR_GIVE_FIELD(ir, typeOfContrib, IFT_HangingNode_type, "type"); // Macro

    if ( typeOfContrib == 0 ) {
        IR_GIVE_FIELD(ir, * masterContribution, IFT_HangingNode_weigths, "weights"); // Macro
    } else   {
        IR_GIVE_FIELD(ir, locoords->at(1), IFT_HangingNode_ksi, "ksi");      // Macro
        IR_GIVE_OPTIONAL_FIELD(ir, locoords->at(2), IFT_HangingNode_eta, "eta"); // Macro
        IR_GIVE_OPTIONAL_FIELD(ir, locoords->at(3), IFT_HangingNode_dzeta, "dzeta"); // Macro
    }

    return IRRT_OK;
}


bool nonzero(double x, double tolerance = 10e-10)
{
    if ( fabs(x) > tolerance ) {
        return 1;
    } else {
        return 0;
    }
}

/**
 * Checks internal data consistency in node.
 * Current implementation checks (when receiver has slave dofs) if receiver has the same
 * coordinate system as master dofManager of slave dof.
 */
int
HangingNode :: checkConsistency()
{
    int result = 1;
    int i, j;
    double suma;

    result = result && Node :: checkConsistency();

    // finds master node
    masterNode = new Node * [ countOfMasterNodes ];
    for ( i = 1; i <= countOfMasterNodes; i++ ) {
        masterNode [ i - 1 ] = dynamic_cast< Node * >( this->domain->giveDofManager( masterDofMngr->at(i) ) );
        if ( !masterNode [ i - 1 ] ) {
            _warning2("checkConsistency: master dofManager is not compatible", 1);
            result = 0;
        }
    }

    // check if receiver has the same coordinate system as master dofManagers
    for ( i = 0; i < countOfMasterNodes; i++ ) {
        if ( !this->hasSameLCS(masterNode [ i ]) ) {
            _warning2("checkConsistency: different lcs for master/slave nodes", 1);
            result = 0;
        }
    }

    result = result && computeMasterContribution();

    // check of master contribution coefficients - SUMA of contributions == 1.0
    if ( masterContribution->giveSize() != countOfMasterNodes ) {
        _warning3("checkConsistency: masterContribution.giveSize()(%f) != countOfMasterDofMngr(%f)", masterContribution->giveSize(), countOfMasterNodes);
        result = 0;
    }

    suma = masterContribution->sum();
    if ( nonzero(suma - 1.0, 10e-9) ) {
        _warning2("checkConsistency: suma of coefficients of master contributions(%.12e) != 1.0", suma);
        result = 0;
    }

    // matching of local and global coordinates
    FloatArray **mnc, coords(3);

    mnc = new FloatArray * [ countOfMasterNodes ];
    for ( i = 0; i < countOfMasterNodes; i++ ) {
        mnc [ i ] = masterNode [ i ]->giveCoordinates();
    }

    if ( typeOfContrib ) {
        coords.zero();
        for ( i = 1; i <= 3; i++ ) {
            for ( j = 1; j <= countOfMasterNodes; j++ ) {
                coords.at(i) += masterContribution->at(j) * mnc [ j - 1 ]->at(i);
            }
        }

        for ( i = 1; i <= 3; i++ ) {
            if ( nonzero( coordinates.at(i) - coords.at(i) ) ) {
                _warning3( "compute_naturalcoord: internal err: coordinates(%12.10f) != coords(%f)", coordinates.at(i), coords.at(i) );
                result = 0;
            }
        }
    }

    // initialize slave dofs (inside check of consistency of receiver and master dof)
    for ( i = 1; i <= numberOfDofs; i++ ) {
        if ( dofArray [ i - 1 ]->giveClassID() == SlaveDofClass ) {
            ( ( SlaveDof * ) dofArray [ i - 1 ] )->initialize(countOfMasterNodes, masterNode, NULL, masterContribution);
        }
    }

    /*
      #ifdef __PARALLEL_MODE
       // check if master in same mode
       if ( parallel_mode != DofManager_local ) {
         for ( i = 1; i <= countOfMasterNodes; i++ ) {
           if ( masterNode [ i - 1 ]->giveParallelMode() != parallel_mode ) {
             _warning2("checkConsistency: mismatch in parallel mode of HangingNode and master", 1);
                result = 0;
           }
         }
       }
	    
      #endif
    */


    // deallocate auxiliary arrays
    deallocAuxArrays();

    return result;
}


int
HangingNode :: computeMasterContribution()
{
    if ( hasLocalCS() && typeOfContrib == 0 ) {
        _error("computeMasterContribution: local cross-section is not supported with this type of contribution");
    }

    switch ( typeOfContrib ) {
    case 211: {
        FEI1dLin(0).evalN(* masterContribution, * locoords, 0.0);
        break;
    }                                                                               // linear truss
    case 312: {
        FEI2dTrLin(0, 0).evalN(* masterContribution, * locoords, 0.0);
        break;
    }                                                                               // linear triangle
    case 412: {
        FEI2dQuadLin(0, 0).evalN(* masterContribution, * locoords, 0.0);
        break;
    }                                                                               // linear rectangle
    case 813: {
        FEI3dHexaLin().evalN(* masterContribution, * locoords, 0.0);
        break;
    }                                                                               // linear hexahedron
    case 321: {
        masterContribution->resize(3);
        masterContribution->at(1) = 0.5 * locoords->at(1) * ( locoords->at(1) - 1.0 );
        masterContribution->at(2) = 0.5 * locoords->at(1) * ( locoords->at(1) + 1.0 );
        masterContribution->at(3) = 1.0 - locoords->at(1) *  locoords->at(1);
        break;
    }                                                                                      // quadratic truss
    case 622: {
        masterContribution->resize(6);
        locoords->at(3) = 1. - locoords->at(1) - locoords->at(2);
        masterContribution->at(1) = ( 2. * locoords->at(1) - 1. ) * locoords->at(1);
        masterContribution->at(2) = ( 2. * locoords->at(2) - 1. ) * locoords->at(2);
        masterContribution->at(3) = ( 2. * locoords->at(3) - 1. ) * locoords->at(3);
        masterContribution->at(4) =  4. * locoords->at(1) * locoords->at(2);
        masterContribution->at(5) =  4. * locoords->at(2) * locoords->at(3);
        masterContribution->at(6) =  4. * locoords->at(3) * locoords->at(1);
        break;
    }                                                                                    // quadratic triangle
    case 822: {
        masterContribution->resize(8);
        masterContribution->at(1) = ( 1. + locoords->at(1) ) * ( 1. + locoords->at(2) ) * 0.25 * ( locoords->at(1) + locoords->at(2) - 1. );
        masterContribution->at(2) = ( 1. - locoords->at(1) ) * ( 1. + locoords->at(2) ) * 0.25 * ( -locoords->at(1) + locoords->at(2) - 1. );
        masterContribution->at(3) = ( 1. - locoords->at(1) ) * ( 1. - locoords->at(2) ) * 0.25 * ( -locoords->at(1) - locoords->at(2) - 1. );
        masterContribution->at(4) = ( 1. + locoords->at(1) ) * ( 1. - locoords->at(2) ) * 0.25 * ( locoords->at(1) - locoords->at(2) - 1. );
        masterContribution->at(5) = 0.5 * ( 1. - locoords->at(1) * locoords->at(1) ) * ( 1. + locoords->at(2) );
        masterContribution->at(6) = 0.5 * ( 1. - locoords->at(1) ) * ( 1. - locoords->at(2) * locoords->at(2) );
        masterContribution->at(7) = 0.5 * ( 1. - locoords->at(1) * locoords->at(1) ) * ( 1. - locoords->at(2) );
        masterContribution->at(8) = 0.5 * ( 1. + locoords->at(1) ) * ( 1. - locoords->at(2) * locoords->at(2) );
        break;
    }                                                                                                           // quadratic rectangle
    case 0: // contribution set in input file
        break;
    default: {
        _warning("computeMasterContribution: unknown type");
        return 0;

        break;
    }
    }

    return 1;
}
