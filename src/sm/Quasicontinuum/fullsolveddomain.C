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
 *  License as published by the Free Software Foundation; eitherc
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

#include "sm/Quasicontinuum/fullsolveddomain.h"

#include "nummet.h"
#include "timestep.h"
#include "element.h"
#include "sparsemtrx.h"
#include "verbose.h"
#include "classfactory.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"

#include "node.h"

namespace oofem {
//REGISTER_Quasicontinuum(QCFullsolveddomain);

QCFullsolveddomain :: QCFullsolveddomain()

// Constructor.
{}

QCFullsolveddomain :: ~QCFullsolveddomain()
// Destructor
{ }

void
QCFullsolveddomain :: initializeFrom(InputRecord &ir)
{
    IR_GIVE_OPTIONAL_FIELD(ir, FullSolvedDomainNodes, _IFT_FullSolvedDomain_nodes);
    IR_GIVE_OPTIONAL_FIELD(ir, FullSolvedDomainElements, _IFT_FullSolvedDomain_elements);
    IR_GIVE_OPTIONAL_FIELD(ir, FullSolvedDomainRadius, _IFT_FullSolvedDomain_radius);
    IR_GIVE_OPTIONAL_FIELD(ir, FullSolvedDomainBox, _IFT_FullSolvedDomain_box);
    // check input format
#ifdef DEBUG
    if ( FullSolvedDomainRadius.giveSize() != 0 && FullSolvedDomainRadius.giveSize() % 4 != 0 ) {
        OOFEM_ERROR("invalid format of FullSolvedDomainRadius");
    }
    if ( FullSolvedDomainBox.giveSize() != 0 && FullSolvedDomainBox.giveSize() % 6 != 0 ) {
        OOFEM_ERROR("invalid format of FullSolvedDomainBox");
    }
#endif
}

void
QCFullsolveddomain :: updateYourself()
{
    // place for adaptivity... //km??
}


bool
QCFullsolveddomain :: isNodeInside(Node *n)
{
    const auto &coordinates = n->giveCoordinates();
    // is tested node in FullSolvedDomainNodes
    if ( FullSolvedDomainNodes.giveSize() != 0 ) {
        for ( int i = 1; i <= FullSolvedDomainNodes.giveSize(); i++ ) {
            if ( n->giveGlobalNumber() == FullSolvedDomainElements.at(i) ) {
                return true;
            }
        }
    }

    // is tested node in FullSolvedDomainElements
    if ( FullSolvedDomainElements.giveSize() != 0 ) {
        for ( int i = 1; i <= FullSolvedDomainElements.giveSize(); i++ ) {
            //if (km??? test zda lezi v i+1 elementu s cislem FullSolvedDomainElements.at(i))
            // return true;
            //}
        }
    }

    // is tested node in FullSolvedDomainRadius
    if ( FullSolvedDomainRadius.giveSize() != 0 ) {
        for ( int i = 0; i <= FullSolvedDomainRadius.giveSize() / 4 - 1; i++ ) {
            FloatArray vector(3);
            vector.at(1) = coordinates.at(1) - FullSolvedDomainRadius.at(4 * i + 1);
            vector.at(2) = coordinates.at(2) - FullSolvedDomainRadius.at(4 * i + 2);
            vector.at(3) = coordinates.at(3) - FullSolvedDomainRadius.at(4 * i + 3);


            if ( vector.computeNorm() <= FullSolvedDomainRadius.at(4 * i + 4) ) {
                return true;
            }
        }
    }

    return false;
}
} // end namespace oofem
