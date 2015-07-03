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

#include "fluidmodel.h"
#include "element.h"
#include "dof.h"
#include "generalboundarycondition.h"
#include "dofmanager.h"

namespace oofem
{
int
FluidModel :: forceEquationNumbering(int id)
{
    // Necessary to number DOFs in special order to guarantee that Skyline matrix factorization to work.
    Domain *domain = this->giveDomain(id);
    TimeStep *currStep = this->giveCurrentStep();

    this->domainNeqs.at(id) = 0;
    this->domainPrescribedNeqs.at(id) = 0;

    // First velocity.
    for ( auto &dman : domain->giveDofManagers() ) {
        for ( Dof *dof: *dman ) {
            DofIDItem type = dof->giveDofID();
            if ( type == V_u || type == V_v || type == V_w ) {
                dof->askNewEquationNumber(currStep);
            }
        }
    }

    for ( auto &elem : domain->giveElements() ) {
        int innodes = elem->giveNumberOfInternalDofManagers();
        for ( int k = 1; k <= innodes; k++ ) {
            DofManager *dman = elem->giveInternalDofManager(k);
            for ( Dof *dof: *dman ) {
                DofIDItem type = dof->giveDofID();
                if ( type == V_u || type == V_v || type == V_w ) {
                    dof->askNewEquationNumber(currStep);
                }
            }
        }
    }

    for ( auto &bc : domain->giveBcs() ) {
        int innodes = bc->giveNumberOfInternalDofManagers();
        for ( int k = 1; k <= innodes; k++ ) {
            DofManager *dman = bc->giveInternalDofManager(k);
            for ( Dof *dof: *dman ) {
                DofIDItem type = dof->giveDofID();
                if ( type == V_u || type == V_v || type == V_w ) {
                    dof->askNewEquationNumber(currStep);
                }
            }
        }
    }

    // Then the rest
    for ( auto &dman : domain->giveDofManagers() ) {
        for ( Dof *dof: *dman ) {
            DofIDItem type = dof->giveDofID();
            if ( !( type == V_u || type == V_v || type == V_w ) ) {
                dof->askNewEquationNumber(currStep);
            }
        }
    }

    for ( auto &elem : domain->giveElements() ) {
        int innodes = elem->giveNumberOfInternalDofManagers();
        for ( int k = 1; k <= innodes; k++ ) {
            DofManager *dman = elem->giveInternalDofManager(k);
            for ( Dof *dof: *dman ) {
                DofIDItem type = dof->giveDofID();
                if ( !( type == V_u || type == V_v || type == V_w ) ) {
                    dof->askNewEquationNumber(currStep);
                }
            }
        }
    }

    for ( auto &bc : domain->giveBcs() ) {
        int innodes = bc->giveNumberOfInternalDofManagers();
        for ( int k = 1; k <= innodes; k++ ) {
            DofManager *dman = bc->giveInternalDofManager(k);
            for ( Dof *dof: *dman ) {
                DofIDItem type = dof->giveDofID();
                if ( !( type == V_u || type == V_v || type == V_w ) ) {
                    dof->askNewEquationNumber(currStep);
                }
            }
        }
    }

    return domainNeqs.at(id);
}
} // end namespace oofem
