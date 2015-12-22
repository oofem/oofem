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

#include "activedof.h"
#include "activebc.h"
#include "domain.h"
#include "dofmanager.h"
#include "contextioerr.h"
#include "activebc.h"
#include "engngm.h"

namespace oofem {
ActiveDof :: ActiveDof(DofManager *aNode, DofIDItem id, int bc) : Dof(aNode, id), equationNumber(0), bc(bc), activeBC(NULL)
{
}

void ActiveDof :: initialize(int cntOfMstrDfMngr, const IntArray &masterNodes, const IntArray *mstrDofID, const FloatArray &mstrContribution)
{
    OOFEM_ERROR("Shouldn't be statically initialized.");
}

ActiveBoundaryCondition *ActiveDof :: giveActiveBoundaryCondition()
{
    if ( !activeBC ) {
        activeBC = dynamic_cast< ActiveBoundaryCondition * >( this->dofManager->giveDomain()->giveBc(bc) );
        if ( !activeBC ) {
            OOFEM_ERROR("No active bc at %d\n", bc);
        }
    }
    return activeBC;
}

int ActiveDof :: giveNumberOfPrimaryMasterDofs()
{
    if ( this->isPrimaryDof() ) {
        return 1;
    }

    int countOfMasterDofs = this->giveActiveBoundaryCondition()->giveNumberOfMasterDofs(this);
    int k = 0;
    for ( int i = 1; i <= countOfMasterDofs; i++ ) {
        k += this->giveMasterDof(i)->giveNumberOfPrimaryMasterDofs();
    }
    return k;
}

bool ActiveDof :: isPrimaryDof()
{
    return this->giveActiveBoundaryCondition()->isPrimaryDof(this);
}

int ActiveDof :: giveNumberOfMasterDofs()
{
    return this->giveActiveBoundaryCondition()->giveNumberOfMasterDofs(this);
}


void ActiveDof :: giveEquationNumbers(IntArray &masterEqNumbers, const UnknownNumberingScheme &s)
{
    if ( this->isPrimaryDof() ) {
        masterEqNumbers.resize(1);
        masterEqNumbers.at(1) = this->giveEquationNumber(s);
        return;
    }

    IntArray mstrEqNmbrs;

    int countOfMasterDofs = this->giveNumberOfMasterDofs();
    masterEqNumbers.preallocate(countOfMasterDofs);
    masterEqNumbers.clear();
    for ( int i = 1; i <= countOfMasterDofs; i++ ) {
        this->giveMasterDof(i)->giveEquationNumbers(mstrEqNmbrs, s);
        masterEqNumbers.followedBy(mstrEqNmbrs);
    }
}

void ActiveDof :: giveDofIDs(IntArray &masterDofIDs)
{
    if ( this->isPrimaryDof() ) {
        masterDofIDs.resize(1);
        masterDofIDs.at(1) = this->giveDofID();
        return;
    }

    IntArray mstrDofIDs;

    int countOfMasterDofs = this->giveNumberOfMasterDofs();
    masterDofIDs.preallocate(countOfMasterDofs);
    masterDofIDs.clear();
    for ( int i = 1; i <= countOfMasterDofs; i++ ) {
        this->giveMasterDof(i)->giveDofIDs(mstrDofIDs);
        masterDofIDs.followedBy(mstrDofIDs);
    }
}

void ActiveDof :: giveMasterDofManArray(IntArray &answer)
{
    if ( this->isPrimaryDof() ) {
        answer.resize(1);
        answer.at(1) = this->giveDofManNumber();
        return;
    }

    IntArray subMasterDofManArray;

    answer.resize( this->giveNumberOfPrimaryMasterDofs() );
    int countOfMasterDofs = this->giveNumberOfMasterDofs();
    for ( int i = 1; i <= countOfMasterDofs; i++ ) {
        this->giveMasterDof(i)->giveMasterDofManArray(subMasterDofManArray);
        answer.followedBy(subMasterDofManArray);
    }
}


void ActiveDof :: giveUnknowns(FloatArray &masterUnknowns, ValueModeType mode, TimeStep *tStep)
{
    FloatArray mstrUnknwns;

    masterUnknowns.resize( this->giveNumberOfPrimaryMasterDofs() );
    int countOfMasterDofs = this->giveNumberOfMasterDofs();
    for ( int k = 1, i = 1; i <= countOfMasterDofs; i++ ) {
        this->giveMasterDof(i)->giveUnknowns(mstrUnknwns, mode, tStep);
        masterUnknowns.copySubVector(mstrUnknwns, k);
        k += mstrUnknwns.giveSize();
    }
}

void ActiveDof :: giveUnknowns(FloatArray &masterUnknowns, PrimaryField &field, ValueModeType mode, TimeStep *tStep)
{
    FloatArray mstrUnknwns;

    masterUnknowns.resize( this->giveNumberOfPrimaryMasterDofs() );
    int countOfMasterDofs = this->giveNumberOfMasterDofs();
    for ( int k = 1, i = 1; i <= countOfMasterDofs; i++ ) {
        this->giveMasterDof(i)->giveUnknowns(mstrUnknwns, field, mode, tStep);
        masterUnknowns.copySubVector(mstrUnknwns, k);
        k += mstrUnknwns.giveSize();
    }
}

void ActiveDof :: computeDofTransformation(FloatArray &primaryMasterContribs)
{
    if ( this->isPrimaryDof() ) {
        primaryMasterContribs.resize(1);
        primaryMasterContribs.at(1) = 1.0;
        return;
    }

    FloatArray masterContribution, subPrimaryMasterContribs;
    this->giveActiveBoundaryCondition()->computeDofTransformation(this, masterContribution);

    primaryMasterContribs.resize( this->giveNumberOfPrimaryMasterDofs() );
    int countOfMasterDofs = this->giveNumberOfMasterDofs();
    for ( int k = 1, i = 1; i <= countOfMasterDofs; i++ ) {
        this->giveMasterDof(i)->computeDofTransformation(subPrimaryMasterContribs);
        subPrimaryMasterContribs.times( masterContribution.at(i) );
        primaryMasterContribs.copySubVector(subPrimaryMasterContribs, k);
        k += subPrimaryMasterContribs.giveSize();
    }
}

double ActiveDof :: giveUnknown(ValueModeType mode, TimeStep *tStep)
{
    return this->giveActiveBoundaryCondition()->giveUnknown(mode, tStep, this);
}

double ActiveDof :: giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *tStep)
{
    return this->giveActiveBoundaryCondition()->giveUnknown(field, mode, tStep, this);
}


int ActiveDof :: __giveEquationNumber() const
{
    return this->equationNumber > 0 ? equationNumber : 0;
}

int ActiveDof :: __givePrescribedEquationNumber()
{
    return this->equationNumber < 0 ? -equationNumber : 0;
}

int ActiveDof :: askNewEquationNumber(TimeStep *tStep)
{
    if ( !this->isPrimaryDof() ) {
        return 0;
    }

    EngngModel *model = dofManager->giveDomain()->giveEngngModel();

    if ( dofManager->giveParallelMode() == DofManager_null ) {
        equationNumber = 0;
        return 0;
    }

    if ( this->hasBc(tStep) ) {
        equationNumber = -model->giveNewPrescribedEquationNumber(dofManager->giveDomain()->giveNumber(), this->dofID);
    } else {
        equationNumber = model->giveNewEquationNumber(dofManager->giveDomain()->giveNumber(), this->dofID);
    }

    return equationNumber;
}

bool ActiveDof :: hasBc(TimeStep *tStep)
{
    return this->giveActiveBoundaryCondition()->hasBc(this, tStep);
}

int ActiveDof :: giveBcId()
{
    return this->bc;
}

void ActiveDof :: setBcId(int bcId)
{
    bc = bcId;
    activeBC = NULL;
}

double ActiveDof :: giveBcValue(ValueModeType mode, TimeStep *tStep)
{
    return this->giveActiveBoundaryCondition()->giveBcValue(this, mode, tStep);
}


// Not sure of initial conditions yet.
bool ActiveDof :: hasIc(TimeStep *tStep) { return false; }
bool ActiveDof :: hasIcOn(ValueModeType type) { return false; }
bool ActiveDof :: hasIc() { return false; }
int ActiveDof :: giveIcId() { return 0; }
InitialCondition *ActiveDof :: giveIc() { return NULL; }


contextIOResultType ActiveDof :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    // Nothing here since the boundary condition deals with all the values.
    return CIO_OK;
}

contextIOResultType ActiveDof :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    return CIO_OK;
}

inline Dof *ActiveDof :: giveMasterDof(int i)
{
    return this->giveActiveBoundaryCondition()->giveMasterDof(this, i);
}

void ActiveDof :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    // No numbering is stored.
}
} // end namespace oofem
