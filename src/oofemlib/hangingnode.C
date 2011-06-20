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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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
 #include <stdlib.h>
#endif

#include "fei1dlin.h"
#include "fei2dtrlin.h"
#include "fei2dquadlin.h"
#include "fei3dhexalin.h"
#include "mathfem.h"

namespace oofem {
HangingNode :: HangingNode(int n, Domain *aDomain) : Node(n, aDomain)
{ 
#ifdef __OOFEG
	consistencyChecked = false;
#endif
}

void
HangingNode :: deallocAuxArrays(void)
{
    delete masterNode;
}


IRResultType
HangingNode :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    Node :: initializeFrom(ir);

    locoords.resize(3);

    // This would be convenient. Renders everything else deprecated / Mikael
    //IR_GIVE_FIELD(ir, masterElement, IFT_HangingNode_masterElement, "masterelement");
    IR_GIVE_FIELD(ir, masterDofMngr, IFT_HangingNode_masters, "masters");
    countOfMasterNodes = masterDofMngr.giveSize();

    IR_GIVE_FIELD(ir, typeOfContrib, IFT_HangingNode_type, "type");

    if ( typeOfContrib == 0 ) {
        IR_GIVE_FIELD(ir, masterContribution, IFT_HangingNode_weigths, "weights");
    } else if ( typeOfContrib < 0 ) {
        //determine ksi, eta, dzeta directly from element-like configuration defined by nodes
    } else {
        IR_GIVE_FIELD(ir, locoords.at(1), IFT_HangingNode_ksi, "ksi");
        IR_GIVE_OPTIONAL_FIELD(ir, locoords.at(2), IFT_HangingNode_eta, "eta");
        IR_GIVE_OPTIONAL_FIELD(ir, locoords.at(3), IFT_HangingNode_dzeta, "dzeta");
    }

    return IRRT_OK;
}


bool nonzero(double x, double tolerance = 10e-10)
{
    return fabs(x) > tolerance;
}


int
HangingNode :: checkConsistency()
{
    int result = 1;
    int i, j;
    double suma;

#ifdef __OOFEG
		if(consistencyChecked) return result;
		consistencyChecked = true;
#endif

    result = result && Node :: checkConsistency();

    // finds master node
    masterNode = new Node * [ countOfMasterNodes ];
    for ( i = 1; i <= countOfMasterNodes; i++ ) {
        masterNode [ i - 1 ] = dynamic_cast< Node * >( this->domain->giveDofManager( masterDofMngr.at(i) ) );
        if ( !masterNode [ i - 1 ] ) {
            _warning2("checkConsistency: master dofManager is not compatible", 1);
            result = 0;
        }
    }

    //compute natural coordinates from the master ones if desired
    if ( typeOfContrib < 0 ) {
        computeNaturalCoordinates();
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
    if ( masterContribution.giveSize() != countOfMasterNodes ) {
        _warning3("checkConsistency: masterContribution.giveSize()(%f) != countOfMasterDofMngr(%f)", masterContribution.giveSize(), countOfMasterNodes);
        result = 0;
    }

    suma = masterContribution.sum();
    // TEMPORARILY TURNED OFF !!!!
    /*
     * if ( nonzero(suma - 1.0, 10e-9) ) {
     *  _warning2("checkConsistency: sum of coefficients of master contributions(%.12e) != 1.0", suma);
     *  result = 0;
     * }
     */

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
                coords.at(i) += masterContribution.at(j) * mnc [ j - 1 ]->at(i);
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
     * #ifdef __PARALLEL_MODE
     * // check if master in same mode
     * if ( parallel_mode != DofManager_local ) {
     *   for ( i = 1; i <= countOfMasterNodes; i++ ) {
     *     if ( masterNode [ i - 1 ]->giveParallelMode() != parallel_mode ) {
     *       _warning2("checkConsistency: mismatch in parallel mode of HangingNode and master", 1);
     *          result = 0;
     *     }
     *   }
     * }
     *
     * #endif
     */


    // deallocate auxiliary arrays
    deallocAuxArrays();

    return result;
}

void
HangingNode :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
  int i;
  for (i=1; i<= masterDofMngr.giveSize(); i++) {
    masterDofMngr.at(i) = f(masterDofMngr.at(i), ERS_DofManager);
  }

  DofManager::updateLocalNumbering(f);
}

int
HangingNode :: computeMasterContribution()
{
    if ( hasLocalCS() && typeOfContrib == 0 ) {
        _error("computeMasterContribution: local cross-section is not supported with this type of contribution");
    }

    // Something like this would be elegant and convenient.
    //Element *e = this->giveDomain()->giveElement(this->masterElement);
    //e->giveDefaultInterpolator()->evalN(masterContribution, locoords, FEIElementGeometryWrapper(e), 0.0);

    switch ( abs(typeOfContrib) ) {
    case 211: {  // linear truss
        FEI1dLin(0).evalN(masterContribution, locoords, FEIVoidCellGeometry(), 0.0);
        break;
    }
    case 312: { // linear triangle
        FEI2dTrLin(0, 0).evalN(masterContribution, locoords, FEIVoidCellGeometry(), 0.0);
        break;
    }
    case 412: { // linear rectangle
        FEI2dQuadLin(0, 0).evalN(masterContribution, locoords, FEIVoidCellGeometry(), 0.0);
        break;
    }
    case 813: { // linear hexahedron
        //FEI3dHexaLin(0, 0).evalN(masterContribution, locoords, FEIVoidCellGeometry(), 0.0);
        double x = locoords.at(1), y = locoords.at(2), z = locoords.at(3);

        masterContribution.resize(8);
        masterContribution.at(1)  = 0.125 * ( 1. - x ) * ( 1. - y ) * ( 1. + z );
        masterContribution.at(2)  = 0.125 * ( 1. - x ) * ( 1. + y ) * ( 1. + z );
        masterContribution.at(3)  = 0.125 * ( 1. + x ) * ( 1. + y ) * ( 1. + z );
        masterContribution.at(4)  = 0.125 * ( 1. + x ) * ( 1. - y ) * ( 1. + z );
        masterContribution.at(5)  = 0.125 * ( 1. - x ) * ( 1. - y ) * ( 1. - z );
        masterContribution.at(6)  = 0.125 * ( 1. - x ) * ( 1. + y ) * ( 1. - z );
        masterContribution.at(7)  = 0.125 * ( 1. + x ) * ( 1. + y ) * ( 1. - z );
        masterContribution.at(8)  = 0.125 * ( 1. + x ) * ( 1. - y ) * ( 1. - z );
        break;
    }
    case 321: { // quadratic truss
        masterContribution.resize(3);
        masterContribution.at(1) = 0.5 * locoords.at(1) * ( locoords.at(1) - 1.0 );
        masterContribution.at(2) = 0.5 * locoords.at(1) * ( locoords.at(1) + 1.0 );
        masterContribution.at(3) = 1.0 - locoords.at(1) *  locoords.at(1);
        break;
    }
    case 622: { // quadratic triangle
        //FEI2dTrQuad(0, 0).evalN(masterContribution, locoords, FEIVoidCellGeometry(), 0.0);
        masterContribution.resize(6);
        locoords.at(3) = 1. - locoords.at(1) - locoords.at(2);
        masterContribution.at(1) = ( 2. * locoords.at(1) - 1. ) * locoords.at(1);
        masterContribution.at(2) = ( 2. * locoords.at(2) - 1. ) * locoords.at(2);
        masterContribution.at(3) = ( 2. * locoords.at(3) - 1. ) * locoords.at(3);
        masterContribution.at(4) =  4. * locoords.at(1) * locoords.at(2);
        masterContribution.at(5) =  4. * locoords.at(2) * locoords.at(3);
        masterContribution.at(6) =  4. * locoords.at(3) * locoords.at(1);
        break;
    }
    case 822: { // quadratic rectangle
        //FEI2dQuadQuad(0, 0).evalN(masterContribution, locoords, FEIVoidCellGeometry(), 0.0);
        masterContribution.resize(8);
        masterContribution.at(1) = ( 1. + locoords.at(1) ) * ( 1. + locoords.at(2) ) * 0.25 * ( locoords.at(1) + locoords.at(2) - 1. );
        masterContribution.at(2) = ( 1. - locoords.at(1) ) * ( 1. + locoords.at(2) ) * 0.25 * ( -locoords.at(1) + locoords.at(2) - 1. );
        masterContribution.at(3) = ( 1. - locoords.at(1) ) * ( 1. - locoords.at(2) ) * 0.25 * ( -locoords.at(1) - locoords.at(2) - 1. );
        masterContribution.at(4) = ( 1. + locoords.at(1) ) * ( 1. - locoords.at(2) ) * 0.25 * ( locoords.at(1) - locoords.at(2) - 1. );
        masterContribution.at(5) = 0.5 * ( 1. - locoords.at(1) * locoords.at(1) ) * ( 1. + locoords.at(2) );
        masterContribution.at(6) = 0.5 * ( 1. - locoords.at(1) ) * ( 1. - locoords.at(2) * locoords.at(2) );
        masterContribution.at(7) = 0.5 * ( 1. - locoords.at(1) * locoords.at(1) ) * ( 1. - locoords.at(2) );
        masterContribution.at(8) = 0.5 * ( 1. + locoords.at(1) ) * ( 1. - locoords.at(2) * locoords.at(2) );
        break;
    }
    case 0: // contribution set in input file
        break;
    default:
        _warning("computeMasterContribution: unknown type");
        return 0;
    }

    return 1;
}

int
HangingNode :: computeNaturalCoordinates()
{
    //Element *e = this->giveDomain()->giveElement(this->masterElement);
    //e->giveDefaultInterpolator()->global2local(locoords, coordinates, FEIElementGeometryWrapper(e), 0.0);
    int i, j;
    //const FloatArray *masterCoords [ countOfMasterNodes ];
    const FloatArray **masterCoords = new const FloatArray * [ countOfMasterNodes ];

    locoords.zero();

    for ( i = 1; i <= countOfMasterNodes; i++ ) { //master nodes of element-like topology
        j = masterDofMngr.at(i);
        masterCoords [ i - 1 ] = new const FloatArray( *this->domain->giveNode( j )->giveCoordinates() );
    }

    //need to extend to other elements
    switch ( typeOfContrib ) {
    case -412: { // linear rectangle
        FEI2dQuadLin(1, 2).global2local(this->locoords, coordinates, FEIVertexListGeometryWrapper(countOfMasterNodes, masterCoords), 0.0);
        break;
    }
    case -813: { // linear hexahedron
        FEI3dHexaLin().global2local(this->locoords, coordinates, FEIVertexListGeometryWrapper(countOfMasterNodes, masterCoords), 0.0);
        break;
    }
    default: {
        _warning("Unknown element-like configuration of master nodes or not implemented element type");

        // clean up
        for ( i = 1; i <= countOfMasterNodes; i++ ) {
            delete masterCoords [ i - 1 ];
        }

        delete [] masterCoords;

        return 0;

        break;
    }
    }

    // clean up
    for ( i = 1; i <= countOfMasterNodes; i++ ) {
        delete masterCoords [ i - 1 ];
    }

    delete [] masterCoords;

    return 1;
}
} // end namespace oofem
