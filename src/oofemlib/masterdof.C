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

#include "masterdof.h"
#include "dofmanager.h"
#include "domain.h"
#include "timestep.h"
#include "boundary.h"
#include "initial.h"
#include "primaryfield.h"
#include "dictionr.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {
MasterDof :: MasterDof(int i, DofManager *aNode, int nbc, int nic, DofIDItem id) : Dof(i, aNode, id)
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

MasterDof :: MasterDof(int i, DofManager *aNode, DofIDItem id) : Dof(i, aNode, id)
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
            return ( BoundaryCondition * ) bcptr;
        }
    }

    _error2("Incompatible BC (%d) applied as Dirichlet/Primary BC", bc);
    return NULL;
}


int MasterDof :: __giveEquationNumber() const
// Returns the number of the equation in the governing system of equations that corres-
// ponds to the receiver. The equationNumber is <0 if the receiver is
// subjected to a boundary condition and then the equationNumber is a Prescribed equation number
// or equationNumber >0.

{
    //if (!equationNumber)
    //  _error ("giveEquationNumber : Dof has undefined equationNumber");
    return ( equationNumber > 0 ) ? equationNumber : 0;
}

int MasterDof :: __givePrescribedEquationNumber()
// Returns the number of the equation in the governing system of equations that corres-
// ponds to the receiver. The equationNumber is <0 if the receiver is
// subjected to a boundary condition and then the equationNumber is a Prescribed equation number
// or equationNumber >0.
{
    //if (!equationNumber)
    //  _error ("giveEquationNumber: Dof has undefined equationNumber");
    return ( equationNumber < 0 ) ? -1 * equationNumber : 0;
}

int MasterDof :: askNewEquationNumber(TimeStep *tStep)
// Returns the newly obtained number of the equation in the governing system
// of equations that corres-
// ponds to the receiver. The equation number is 0 if the receiver is
// subjected to a boundary condition, else it is n+1, where n is the
// equation number of the most recently numbered degree of freedom.
{
    EngngModel *model;
    model = ( EngngModel * ) ( dofManager->giveDomain()->giveEngngModel() );

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
        OOFEM_ERROR("MasterDof :: giveIc() - does not know yet if has InitCond or not");
    }
}


double MasterDof :: giveUnknown(EquationID type, ValueModeType mode, TimeStep *stepN)
// The key method of class Dof. Returns the value of the unknown 'u'
// (e.g., the displacement) of the receiver, at stepN. This value may,
// or may not be already available. It may depend on a boundary (if it
// is not a predicted unknown) or initial condition. stepN is not the
// current time step n, it is assumed to be the previous one (n-1).
{
    /*   double value ;
     *
     *  if (stepN -> isTheCurrentTimeStep()) {
     *     if (unknowns -> includes(u))                       // already known
     *  value = unknowns -> at(u) ;
     *     else {
     *  if (stepN->giveNumber()==0) {                   // step 0
     *    if (this->hasIcOn(u))                        //   init. cond.
     *       value = this -> giveIc() -> give(u) ;
     *    else                                         //   no init. cond.
     *       value = 0. ;}
     *  else if (this->hasBc() && islower(u))           // bound . cond.
     *    value = this -> giveBc() -> give(u,stepN) ;
     *  else                                            // compute it !
     *    value = this -> computeUnknown(u,stepN) ;
     *  unknowns -> add(u,value) ;}}
     *
     *  else
     *     value = this -> givePastUnknown(u,stepN) ;         // the previous step
     *
     *  return value ;   */
    double value;

#ifdef DEBUG
    //if (type != this->giveUnknownType ())
    // _error ("giveUnknown: Noncompatible Request");
#endif

#ifdef __PARALLEL_MODE
    if ( dofManager->giveParallelMode() == DofManager_null ) {
        return 0.0;
    }

#endif

    // first try if IC apply
    if ( stepN->giveNumber() == dofManager->giveDomain()->giveEngngModel()->giveNumberOfTimeStepWhenIcApply() ) { // step when Ic apply
        if ( this->hasIcOn(mode) ) {
            value = this->giveIc()->give(mode);
        } else if ( this->hasBc(stepN) ) {
            value = this->giveBcValue(mode, stepN);
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
    //         giveUnknownComponent(type, mode, stepN, dofManager->giveDomain(), this) );

    //         int hash = dofManager->giveDomain()->giveEngngModel()->giveUnknownDictHashIndx(type, mode, stepN);
    //         if ( unknowns->includes(hash) ) {
    //             return unknowns->at(hash);
    //         } else {
    //             _error2( "giveUnknown:  Dof unknowns dictionary does not contain unknown of value mode (%s)", __ValueModeTypeToString(mode) );
    //         }
    //  }

    if ( !dofManager->giveDomain()->giveEngngModel()->requiresUnknownsDictionaryUpdate() && this->hasBc(stepN) ) {
        value = this->giveBcValue(mode, stepN);
        return value;
    }

    return ( dofManager->giveDomain()->giveEngngModel()->
            giveUnknownComponent(type, mode, stepN, dofManager->giveDomain(), this) );
}

double MasterDof :: giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *stepN)
// The key method of class Dof. Returns the value of the unknown field
// (e.g., the displacement) associated to the receiver, at stepN.
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
    if ( stepN->giveNumber() == dofManager->giveDomain()->giveEngngModel()->giveNumberOfTimeStepWhenIcApply() ) { // step when Ic apply
        if ( this->hasIcOn(mode) ) {
            value = this->giveIc()->give(mode);
        } else {
            value = 0.;
        }

        return value;
    }

    // then ask bor BC
    if ( this->hasBc(stepN) ) {    // bound . cond.
        //value = this -> giveBcValue(giveUnknownType(),mode,stepN) ;
        value = this->giveBcValue(mode, stepN);
        return value;
    }

    // try ask field for unknown
    return ( field.giveUnknownValue(this, mode, stepN) );
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
        OOFEM_ERROR("MasterDof :: hasBc() - does not know yet if has InitCond or not");
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
        OOFEM_ERROR("MasterDof :: hasIc - does not know yet if has InitCond or not");
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

void MasterDof :: updateUnknownsDictionary(TimeStep *tStep, EquationID type,
                                           ValueModeType mode, double dofValue)
{
    // Updates the receiver's unknown dictionary at end of step.
    // to value dofValue.
#ifdef DEBUG
    // if (type != this->giveUnknownType ())
    //  _error ("giveUnknown: Noncompatible Request");
#endif

    int hash = dofManager->giveDomain()->giveEngngModel()->giveUnknownDictHashIndx(type, mode, tStep);
    unknowns->at(hash) = dofValue;
}

void MasterDof :: giveUnknownsDictionaryValue(TimeStep *tStep, EquationID type,
                                              ValueModeType mode, double &dofValue)
{
    // Updates the receiver's unknown dictionary at end of step.
    // to value dofValue.
#ifdef DEBUG
    // if (type != this->giveUnknownType ())
    //  _error ("giveUnknown: Noncompatible Request");
#endif

    int hash = dofManager->giveDomain()->giveEngngModel()->giveUnknownDictHashIndx(type, mode, tStep);
    dofValue = unknowns->at(hash);
}

void MasterDof :: printYourself()
// Prints the receiver on screen.
{
    printf( "dof %d  of %s %d :\n", number, dofManager->giveClassName(), dofManager->giveNumber() );
    printf("equation %d    bc %d \n", equationNumber, bc);

    // printOutputAt (node->giveDomain()->giveEngngModel()->giveCurrentStep());
}


contextIOResultType MasterDof :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full node context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;
    if ( stream == NULL ) {
        _error("saveContex : can't write into NULL stream");
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
        _error("restoreContex : can't write into NULL stream");
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


#ifdef __PARALLEL_MODE
int
MasterDof :: packUnknowns(CommunicationBuffer &buff, EquationID type, ValueModeType mode, TimeStep *stepN)
{
    return buff.packDouble( this->giveUnknown(type, mode, stepN) );
}

int
MasterDof :: unpackAndUpdateUnknown(CommunicationBuffer &buff, EquationID type,
                                    ValueModeType mode, TimeStep *stepN)
{
    EngngModel :: EngngModel_UpdateMode __umode = EngngModel :: EngngModel_Unknown_Mode;
    double value;
    int result = 0;
    // if dof belonging to shared or remote DofManager, engng model unknowns are updated
    // to accommodate remote contribution or "prescribed" remote values.
    // The unknown dictionary is not updated, it is engng model job to update
    // all unknowns dictionaries.

    if ( dofManager->giveParallelMode() == DofManager_shared ) {
        __umode = EngngModel :: EngngModel_SUMM_Mode;
    } else if ( dofManager->giveParallelMode() == DofManager_remote ) {
        __umode = EngngModel :: EngngModel_SET_Mode;
    } else {
        _error("unpackAndUpdateUnknown: unknown dofManager ParallelMode");
    }

    result = buff.unpackDouble(value);
    dofManager->giveDomain()->giveEngngModel()->
    updateUnknownComponent(type, mode, stepN, this->__giveEquationNumber(),
                           value, __umode);
    return result;
}

#endif
} // end namespace oofem
