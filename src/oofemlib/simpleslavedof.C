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

#include "simpleslavedof.h"
#include "dofmanager.h"
#include "domain.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {
SimpleSlaveDof :: SimpleSlaveDof(int i, DofManager *aNode, int master, DofIDItem id) : Dof(i, aNode, id)
    // Constructor. Creates a new d.o.f., with number i, belonging
    // to aNode with bc=nbc, ic=nic
{
    masterDofMngr = master;
    masterDofIndx = -1;
}


SimpleSlaveDof :: SimpleSlaveDof(int i, DofManager *aNode, DofIDItem id) : Dof(i, aNode, id)
    // Constructor. Creates a new d.o.f., with number i, belonging
    // to aNode with bc=nbc, ic=nic
{
    masterDofMngr = -1;
    masterDofIndx = -1;
}


Dof *SimpleSlaveDof :: giveMasterDof() const
{
    // returns reference to master dof
    // checks dof compatibility and slave to slave references

    if ( this->masterDofIndx == -1 ) {
        this->masterDofIndx = dofManager->giveDomain()->giveDofManager(masterDofMngr)
                              ->findDofWithDofId(this->dofID);

        if ( this->masterDofIndx ) {
            /*
             * classType masterDofCT = dofManager->giveDomain()->giveDofManager(masterDofMngr)->
             * giveDof(masterDofIndx)->giveClassID();
             * if ((masterDofCT != MasterDofClass)&&(masterDofCT != SharedMasterDofClass)&&(masterDofCT != RemoteMasterDofClass)) {
             *
             * _error ("giveMasterDof: slaveDof to slaveDof reference not allowed");
             * }
             */
        } else {
            _error("giveMasterDof: no dof with dofID in master found");
        }
    }

    return dofManager->giveDomain()->giveDofManager(masterDofMngr)->giveDof(masterDofIndx);
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


double SimpleSlaveDof :: giveUnknown(EquationID type, ValueModeType mode, TimeStep *stepN)
// The key method of class Dof. Returns the value of the unknown 'u'
// (e.g., the displacement) of the receiver, at stepN. This value may,
// or may not, be already available. It may depend on a boundary (if it
// is not a predicted unknown) or initial condition. stepN is not the
// current time step n, it is assumed to be the previous one (n-1).
{
    return this->giveMasterDof()->giveUnknown(type, mode, stepN);
}

double SimpleSlaveDof :: giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *stepN)
{
    return this->giveMasterDof()->giveUnknown(field, mode, stepN);
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


contextIOResultType SimpleSlaveDof :: saveContext(DataStream *stream, ContextMode mode, void *obj)
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
#ifdef __PARALLEL_MODE
        if ( mode & CM_DefinitionGlobal ) {
            int _masterGlobNum = dofManager->giveDomain()->giveDofManager(masterDofMngr)->giveGlobalNumber();
            if ( !stream->write(& _masterGlobNum, 1) ) {
                THROW_CIOERR(CIO_IOERR);
            }
        } else {
            if ( !stream->write(& masterDofMngr, 1) ) {
                THROW_CIOERR(CIO_IOERR);
            }
        }

#else
        if ( !stream->write(& masterDofMngr, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

#endif
    }

    return CIO_OK;
}


contextIOResultType SimpleSlaveDof :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
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
        if ( !stream->read(& masterDofMngr, 1) ) {
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
