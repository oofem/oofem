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

#include "simpleslavedof.h"
#include "dofmanager.h"
#include "domain.h"
#include "datastream.h"
#include "entityrenumberingscheme.h"
#include "contextioerr.h"

namespace oofem {
SimpleSlaveDof :: SimpleSlaveDof(DofManager *aNode, int master, DofIDItem id) : Dof(aNode, id)
{
    masterDofMngr = master;
    masterDofIndx = -1;
}


SimpleSlaveDof :: SimpleSlaveDof(DofManager *aNode, DofIDItem id) : Dof(aNode, id)
{
    masterDofMngr = -1;
    masterDofIndx = -1;
}


Dof *SimpleSlaveDof :: giveMasterDof() const
{
    // returns reference to master dof
    // checks dof compatibility and slave to slave references

    return dofManager->giveDomain()->giveDofManager(masterDofMngr)->giveDofWithID(this->dofID);
}

BoundaryCondition *SimpleSlaveDof :: giveBc()
// Returns the boundary condition the receiver is subjected to.
{
    return this->giveMasterDof()->giveBc();
}


int SimpleSlaveDof :: __giveEquationNumber() const
// Returns the number of the equation in the governing system of equations that corres-
// ponds to the receiver. The equation number is 0 if the receiver is
// subjected to a boundary condition, else it is n+1, where n is the
// equation number of the most recently numbered degree of freedom.
{
    return this->giveMasterDof()->__giveEquationNumber();
}

int SimpleSlaveDof :: __givePrescribedEquationNumber()
// Returns the number of the equation in the governing system of equations that corres-
// ponds to the receiver. The equation number is 0 if the receiver is
// subjected to a boundary condition, else it is n+1, where n is the
// equation number of the most recently numbered degree of freedom.
{
    return this->giveMasterDof()->__givePrescribedEquationNumber();
}

InitialCondition *SimpleSlaveDof :: giveIc()
// Returns the initial condition on the receiver. Not used.
{
    return this->giveMasterDof()->giveIc();
}


double SimpleSlaveDof :: giveUnknown(ValueModeType mode, TimeStep *tStep)
{
    return this->giveMasterDof()->giveUnknown(mode, tStep);
}

double SimpleSlaveDof :: giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *tStep)
{
    return this->giveMasterDof()->giveUnknown(field, mode, tStep);
}


bool SimpleSlaveDof :: hasBc(TimeStep *tStep)
// Returns True if the receiver is subjected to a boundary condition, else
// returns False. If necessary, reads the answer in the data file.
{
    return this->giveMasterDof()->hasBc(tStep);
}


bool SimpleSlaveDof :: hasIc()
// Returns True if the receiver is subjected to an initial condition,
// else returns False.
{
    return this->giveMasterDof()->hasIc();
}


bool SimpleSlaveDof :: hasIcOn(ValueModeType u)
// Returns True if the unknown 'u' (e.g., the displacement 'd') of the
// receiver is subjected to an initial condition, else returns False.
{
    return this->giveMasterDof()->hasIcOn(u);
}

int SimpleSlaveDof :: giveBcId()
{
    return this->giveMasterDof()->giveBcId();
}

int SimpleSlaveDof :: giveIcId()
{
    return this->giveMasterDof()->giveIcId();
}


double
SimpleSlaveDof :: giveBcValue(ValueModeType mode, TimeStep *tStep)
{
    return this->giveMasterDof()->giveBcValue(mode, tStep);
}


contextIOResultType SimpleSlaveDof :: saveContext(DataStream &stream, ContextMode mode, void *obj)
//
// saves full node context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;
    if ( ( iores = Dof :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( mode & CM_Definition ) {

        if ( mode & CM_DefinitionGlobal ) {
            int _masterGlobNum = dofManager->giveDomain()->giveDofManager(masterDofMngr)->giveGlobalNumber();
            if ( !stream.write(_masterGlobNum) ) {
                THROW_CIOERR(CIO_IOERR);
            }
        } else {
            if ( !stream.write(masterDofMngr) ) {
                THROW_CIOERR(CIO_IOERR);
            }
        }
    }

    return CIO_OK;
}


contextIOResultType SimpleSlaveDof :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
//
// restores full node context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = Dof :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( mode & CM_Definition ) {
        if ( !stream.read(masterDofMngr) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    this->masterDofIndx = -1;


    return CIO_OK;
}

void
SimpleSlaveDof :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    masterDofMngr = f(masterDofMngr, ERS_DofManager);
}
} // end namespace oofem
