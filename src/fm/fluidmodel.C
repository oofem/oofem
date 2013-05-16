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

#include "fluidmodel.h"
#include "element.h"
#include "generalboundarycondition.h"
#include "dofmanager.h"

namespace oofem
{

int
FluidModel :: forceEquationNumbering ( int id )
{
    // Necessary to number DOFs in special order to guarantee that Skyline matrix factorization to work.
    Domain *domain = this->giveDomain ( id );
    TimeStep *currStep = this->giveCurrentStep();
    IntArray loc;

    this->domainNeqs.at ( id ) = 0;
    this->domainPrescribedNeqs.at ( id ) = 0;

    int nnodes = domain->giveNumberOfDofManagers();
    int nelem  = domain->giveNumberOfElements();
    int nbc    = domain->giveNumberOfBoundaryConditions();

    // First velocity.
    for ( int i = 1; i <= nnodes; i++ ) {
        DofManager *dman = domain->giveDofManager ( i );
        int ndofs = dman->giveNumberOfDofs();
        for ( int j = 1; j <= ndofs; j++ ) {
            Dof *jDof = dman->giveDof ( j );
            DofIDItem type = jDof->giveDofID();
            if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                jDof->askNewEquationNumber ( currStep );
            }
        }
    }
    for ( int i = 1; i <= nelem; ++i ) {
        Element *elem = domain->giveElement ( i );
        int innodes = elem->giveNumberOfInternalDofManagers();
        for ( int k = 1; k <= innodes; k++ ) {
            DofManager *dman = elem->giveInternalDofManager ( k );
            int ndofs = dman->giveNumberOfDofs();
            for ( int j = 1; j <= ndofs; j++ ) {
                Dof *jDof = dman->giveDof ( j );
                DofIDItem type = jDof->giveDofID();
                if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                    jDof->askNewEquationNumber ( currStep );
                }
            }
        }
    }
    for ( int i = 1; i <= nbc; ++i ) {
        GeneralBoundaryCondition *bc = domain->giveBc ( i );
        int innodes = bc->giveNumberOfInternalDofManagers();
        for ( int k = 1; k <= innodes; k++ ) {
            DofManager *dman = bc->giveInternalDofManager ( k );
            int ndofs = dman->giveNumberOfDofs();
            for ( int j = 1; j <= ndofs; j++ ) {
                Dof *jDof = dman->giveDof ( j );
                DofIDItem type = jDof->giveDofID();
                if ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) {
                    jDof->askNewEquationNumber ( currStep );
                }
            }
        }
    }
    // Then the rest
    for ( int i = 1; i <= nnodes; i++ ) {
        DofManager *dman = domain->giveDofManager ( i );
        int ndofs = dman->giveNumberOfDofs();
        for ( int j = 1; j <= ndofs; j++ ) {
            Dof *jDof = dman->giveDof ( j );
            DofIDItem type = jDof->giveDofID();
            if ( ! ( ( type == V_u ) || ( type == V_v ) || ( type == V_w ) ) ) {
                jDof->askNewEquationNumber ( currStep );
            }
        }
    }
    for ( int i = 1; i <= nelem; ++i ) {
        Element *elem = domain->giveElement ( i );
        int innodes = elem->giveNumberOfInternalDofManagers();
        for ( int k = 1; k <= innodes; k++ ) {
            DofManager *dman = elem->giveInternalDofManager ( k );
            int ndofs = dman->giveNumberOfDofs();
            for ( int j = 1; j <= ndofs; j++ ) {
                Dof *jDof = dman->giveDof ( j );
                DofIDItem type = jDof->giveDofID();
                if ( ! ( type == V_u || type == V_v || type == V_w ) ) {
                    jDof->askNewEquationNumber ( currStep );
                }
            }
        }
    }
    for ( int i = 1; i <= nbc; ++i ) {
        GeneralBoundaryCondition *bc = domain->giveBc ( i );
        int innodes = bc->giveNumberOfInternalDofManagers();
        for ( int k = 1; k <= innodes; k++ ) {
            DofManager *dman = bc->giveInternalDofManager ( k );
            int ndofs = dman->giveNumberOfDofs();
            for ( int j = 1; j <= ndofs; j++ ) {
                Dof *jDof = dman->giveDof ( j );
                DofIDItem type = jDof->giveDofID();
                if ( ! ( type == V_u || type == V_v || type == V_w ) ) {
                    jDof->askNewEquationNumber ( currStep );
                }
            }
        }
    }

    return domainNeqs.at ( id );
}

} // end namespace oofem
