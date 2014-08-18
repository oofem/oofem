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

#include "masterdof.h"
#include "dofmanager.h"
#include "domain.h"
#include "timestep.h"
#include "boundarycondition.h"
#include "initialcondition.h"
#include "primaryfield.h"
#include "dictionary.h"
#include "datastream.h"
#include "contextioerr.h"
#include "engngm.h"

namespace oofem {
MasterDof :: MasterDof(DofManager *aNode, int nbc, int nic, DofIDItem id) : Dof(aNode, id)
    // Constructor. Creates a new d.o.f., with number i, belonging
    // to aNode with bc=nbc, ic=nic
{
    equationNumber = 0;                         // means "uninitialized"
    bc             = nbc;
    ic             = nic;
    unknowns       = new Dictionary();
    /*   unknowns       = new Dictionary() ;           // unknown size ?
     * pastUnknowns   = NULL ; */
}

MasterDof :: MasterDof(DofManager *aNode, DofIDItem id) : Dof(aNode, id)
{
    ic = bc = equationNumber = 0;                        // means "uninitialized"
    unknowns = new Dictionary();
}

BoundaryCondition *MasterDof :: giveBc()
// Returns the boundary condition the receiver is subjected to.
{
    if ( bc ) {
        GeneralBoundaryCondition *bcptr = dofManager->giveDomain()->giveBc(bc);
        if ( bcptr->giveType() == DirichletBT ) {
            return static_cast< BoundaryCondition * >(bcptr);
        }
    }

    OOFEM_ERROR("Incompatible BC (%d) applied as Dirichlet/Primary BC", bc);
    return NULL;
}


int MasterDof :: __giveEquationNumber() const
// Returns the number of the equation in the governing system of equations that corres-
// ponds to the receiver. The equationNumber is <0 if the receiver is
// subjected to a boundary condition and then the equationNumber is a Prescribed equation number
// or equationNumber >0.

{
    //if (!equationNumber)
    //  OOFEM_ERROR("Dof has undefined equationNumber");
    return ( equationNumber > 0 ) ? equationNumber : 0;
}

int MasterDof :: __givePrescribedEquationNumber()
// Returns the number of the equation in the governing system of equations that corres-
// ponds to the receiver. The equationNumber is <0 if the receiver is
// subjected to a boundary condition and then the equationNumber is a Prescribed equation number
// or equationNumber >0.
{
    //if (!equationNumber)
    //  OOFEM_ERROR("Dof has undefined equationNumber");
    return ( equationNumber < 0 ) ? -1 * equationNumber : 0;
}

int MasterDof :: askNewEquationNumber(TimeStep *tStep)
// Returns the newly obtained number of the equation in the governing system
// of equations that corres-
// ponds to the receiver. The equation number is 0 if the receiver is
// subjected to a boundary condition, else it is n+1, where n is the
// equation number of the most recently numbered degree of freedom.
{
    EngngModel *model = dofManager->giveDomain()->giveEngngModel();

#ifdef __PARALLEL_MODE
    if ( dofManager->giveParallelMode() == DofManager_null ) {
        equationNumber = 0;
        return 0;
    }

#endif
    if ( this->hasBc(tStep) ) {
        equationNumber = -1 * model->giveNewPrescribedEquationNumber(dofManager->giveDomain()->giveNumber(), this->dofID);
    } else {
        equationNumber = model->giveNewEquationNumber(dofManager->giveDomain()->giveNumber(), this->dofID);
    }

    return equationNumber;
}


InitialCondition *MasterDof :: giveIc()
// Returns the initial condition on the receiver. Not used.
{
    if ( ic ) {
        return  ( dofManager->giveDomain()->giveIc(ic) );
    } else {
        OOFEM_ERROR("does not know yet if has InitCond or not");
        return NULL;
    }
}


double MasterDof :: giveUnknown(ValueModeType mode, TimeStep *tStep)
// The key method of class Dof. Returns the value of the unknown 'u'
// (e.g., the displacement) of the receiver, at tStep. This value may,
// or may not be already available. It may depend on a boundary (if it
// is not a predicted unknown) or initial condition. tStep is not the
// current time step n, it is assumed to be the previous one (n-1).
{
    double value;

#ifdef DEBUG
    //if (type != this->giveUnknownType ())
    // OOFEM_ERROR("Noncompatible Request");
#endif

#ifdef __PARALLEL_MODE
    if ( dofManager->giveParallelMode() == DofManager_null ) {
        return 0.0;
    }

#endif

    // first try if IC apply
    if ( tStep->giveNumber() == dofManager->giveDomain()->giveEngngModel()->giveNumberOfTimeStepWhenIcApply() ) { // step when Ic apply
        if ( this->hasIcOn(mode) ) {
            value = this->giveIc()->give(mode);
        } else if ( this->hasBc(tStep) ) {
            value = this->giveBcValue(mode, tStep);
        } else {
            value = 0.;
        }

        return value;
    }

    //  if ( dofManager->giveDomain()->giveEngngModel()->requiresUnknownsDictionaryUpdate() ) {
    // if this feature is active, engng model must ensure
    // valid data in unknowns dictionary
    // the e-model must ensure that bc and ic values are correctly set in unknowns dictionaries
    // they could not be obtained from bc (which are typically incremental)
    // directly since dictionaries keep the history.

    //    return ( dofManager->giveDomain()->giveEngngModel()->
    //         giveUnknownComponent(mode, tStep, dofManager->giveDomain(), this) );

    //         int hash = dofManager->giveDomain()->giveEngngModel()->giveUnknownDictHashIndx(mode, tStep);
    //         if ( unknowns->includes(hash) ) {
    //             return unknowns->at(hash);
    //         } else {
    //             OOFEM_ERROR(Dof unknowns dictionary does not contain unknown of value mode (%s)", __ValueModeTypeToString(mode));
    //         }
    //  }

    if ( !dofManager->giveDomain()->giveEngngModel()->requiresUnknownsDictionaryUpdate() && this->hasBc(tStep) ) {
        value = this->giveBcValue(mode, tStep);
        return value;
    }

    return ( dofManager->giveDomain()->giveEngngModel()->
            giveUnknownComponent(mode, tStep, dofManager->giveDomain(), this) );
}

double MasterDof :: giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *tStep)
// The key method of class Dof. Returns the value of the unknown field
// (e.g., the displacement) associated to the receiver, at tStep.
{
    double value;

#ifdef DEBUG
#endif

#ifdef __PARALLEL_MODE
    if ( dofManager->giveParallelMode() == DofManager_null ) {
        return 0.0;
    }

#endif

    // first try if IC apply
    if ( tStep->giveNumber() == dofManager->giveDomain()->giveEngngModel()->giveNumberOfTimeStepWhenIcApply() ) { // step when Ic apply
        if ( this->hasIcOn(mode) ) {
            value = this->giveIc()->give(mode);
        } else {
            value = 0.;
        }

        return value;
    }

    // then ask bor BC
    if ( this->hasBc(tStep) ) {    // bound . cond.
        //value = this -> giveBcValue(giveUnknownType(),mode,tStep) ;
        value = this->giveBcValue(mode, tStep);
        return value;
    }

    // try ask field for unknown
    return ( field.giveUnknownValue(this, mode, tStep) );
}


bool MasterDof :: hasBc(TimeStep *tStep)
// Returns True if the receiver is subjected to a boundary condition, else
// returns False. If necessary, reads the answer in the data file.
{
#ifdef __PARALLEL_MODE
    if ( dofManager->giveParallelMode() == DofManager_null ) {
        return true;
    }

#endif

    if ( bc == -1 ) {
        OOFEM_ERROR("does not know yet if has InitCond or not");
    }

    if ( bc ) {
        return this->giveBc()->isImposed(tStep);
    } else {
        return false;
    }
}


bool MasterDof :: hasIc()
// Returns True if the receiver is subjected to an initial condition,
// else returns False.
{
    if ( ic == -1 ) {
        OOFEM_ERROR("does not know yet if has InitCond or not");
    }

    return ic > 0;
}


bool MasterDof :: hasIcOn(ValueModeType u)
// Returns True if the unknown 'u' (e.g., the displacement 'd') of the
// receiver is subjected to an initial condition, else returns False.
{
    if ( this->hasIc() ) {
        return this->giveIc()->hasConditionOn(u);
    } else {
        return false;
    }
}

int MasterDof :: giveBcId()
{
    return this->bc;
}

int MasterDof :: giveIcId()
{
    return this->ic;
}

void MasterDof :: updateYourself(TimeStep *tStep)
// Updates the receiver at end of step.
{
    Dof :: updateYourself(tStep);
    /*   delete pastUnknowns ;
     * pastUnknowns = unknowns ;
     * unknowns     = new Dictionary() ; */
}

void MasterDof :: updateUnknownsDictionary(TimeStep *tStep, ValueModeType mode, double dofValue)
{
    // Updates the receiver's unknown dictionary at end of step.
    // to value dofValue.

    int hash = dofManager->giveDomain()->giveEngngModel()->giveUnknownDictHashIndx(mode, tStep);
    unknowns->at(hash) = dofValue;
}

void MasterDof :: giveUnknownsDictionaryValue(TimeStep *tStep, ValueModeType mode, double &dofValue)
{
    // Updates the receiver's unknown dictionary at end of step.
    // to value dofValue.

    int hash = dofManager->giveDomain()->giveEngngModel()->giveUnknownDictHashIndx(mode, tStep);
    dofValue = unknowns->at(hash);
}

void MasterDof :: printYourself()
{
    printf( "dof %d  of %s %d :\n", dofID, dofManager->giveClassName(), dofManager->giveNumber() );
    printf("equation %d    bc %d \n", equationNumber, bc);
}


contextIOResultType MasterDof :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full node context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;
    if ( stream == NULL ) {
        OOFEM_ERROR("can't write into NULL stream");
    }

    if ( ( iores = Dof :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( mode & CM_Definition ) {
        if ( !stream->write(& bc, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream->write(& ic, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    // store equation number of receiver
    if ( !stream->write(& equationNumber, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( mode & CM_UnknownDictState ) || ( dofManager->giveDomain()->giveEngngModel()->requiresUnknownsDictionaryUpdate() ) ) {
        if ( ( iores = unknowns->saveContext(stream, mode, obj) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


contextIOResultType MasterDof :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full node context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;
    if ( stream == NULL ) {
        OOFEM_ERROR("can't write into NULL stream");
    }

    if ( ( iores = Dof :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( mode & CM_Definition ) {
        if ( !stream->read(& bc, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream->read(& ic, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }


    // read equation number of receiver
    if ( !stream->read(& equationNumber, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( mode & CM_UnknownDictState ) || ( dofManager->giveDomain()->giveEngngModel()->requiresUnknownsDictionaryUpdate() ) ) {
        if ( ( iores = unknowns->restoreContext(stream, mode, obj) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}
} // end namespace oofem
