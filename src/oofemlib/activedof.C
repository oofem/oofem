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

#include "activedof.h"
#include "activebc.h"
#include "node.h"
#include "datastream.h"
#include "contextioerr.h"
#include "activebc.h"

namespace oofem {
ActiveDof :: ActiveDof(int n, DofManager *aNode, int bc, DofIDItem id) : Dof(n, aNode, id), bc(bc), activeBC(NULL)
{
}

void ActiveDof :: initialize(int cntOfMstrDfMngr, const IntArray &masterNodes, const IntArray *mstrDofID, const FloatArray &mstrContribution)
{
    OOFEM_ERROR("ActiveDof :: initialize - Shouldn't be statically initialized.");
}

ActiveBoundaryCondition *ActiveDof :: giveActiveBoundaryCondition()
{
    if (!activeBC) {
        activeBC = dynamic_cast<ActiveBoundaryCondition*>(this->dofManager->giveDomain()->giveBc(bc));
        if (!activeBC) {
            OOFEM_ERROR2("ActiveBoundaryCondition :: giveActiveBoundaryCondition - No active bc at %d\n", bc);
        }
    }
    return activeBC;
}

int ActiveDof :: giveNumberOfPrimaryMasterDofs()
{
    int countOfMasterDofs = this->giveActiveBoundaryCondition()->giveNumberOfMasterDofs(this);
    int k = 0;
    for (int i = 1; i <= countOfMasterDofs; i++ ) {
        k += this->giveMasterDof(i)->giveNumberOfPrimaryMasterDofs();
    }
    return k;
}

int ActiveDof :: giveNumberOfMasterDofs()
{
    return this->giveActiveBoundaryCondition()->giveNumberOfMasterDofs(this);
}


void ActiveDof :: giveEquationNumbers(IntArray &masterEqNumbers, const UnknownNumberingScheme &s)
{
    IntArray mstrEqNmbrs;

    masterEqNumbers.resize( this->giveNumberOfPrimaryMasterDofs() );
    int countOfMasterDofs = this->giveNumberOfMasterDofs();
    for (int k = 1, i = 1; i <= countOfMasterDofs; i++ ) {
        this->giveMasterDof(i)->giveEquationNumbers(mstrEqNmbrs, s);
        masterEqNumbers.copySubVector(mstrEqNmbrs, k);
        k += mstrEqNmbrs.giveSize();
    }
}

void ActiveDof :: giveMasterDofManArray(IntArray &answer)
{
    IntArray subMasterDofManArray;

    answer.resize( this->giveNumberOfPrimaryMasterDofs() );
    int countOfMasterDofs = this->giveNumberOfMasterDofs();
    for ( int k = 1, i = 1; i <= countOfMasterDofs; i++ ) {
        this->giveMasterDof(i)->giveMasterDofManArray(subMasterDofManArray);
        answer.copySubVector(subMasterDofManArray, k);
        k += subMasterDofManArray.giveSize();
    }
}


void ActiveDof :: giveUnknowns(FloatArray &masterUnknowns, EquationID type, ValueModeType mode, TimeStep *stepN)
{
    FloatArray mstrUnknwns;

    masterUnknowns.resize( this->giveNumberOfPrimaryMasterDofs() );
    int countOfMasterDofs = this->giveNumberOfMasterDofs();
    for (int k = 1, i = 1; i <= countOfMasterDofs; i++ ) {
        this->giveMasterDof(i)->giveUnknowns(mstrUnknwns, type, mode, stepN);
        masterUnknowns.copySubVector(mstrUnknwns, k);
        k += mstrUnknwns.giveSize();
    }
}

void ActiveDof :: giveUnknowns(FloatArray &masterUnknowns, PrimaryField &field, ValueModeType mode, TimeStep *stepN)
{
    FloatArray mstrUnknwns;

    masterUnknowns.resize( this->giveNumberOfPrimaryMasterDofs() );
    int countOfMasterDofs = this->giveNumberOfMasterDofs();
    for (int k = 1, i = 1; i <= countOfMasterDofs; i++ ) {
        this->giveMasterDof(i)->giveUnknowns(mstrUnknwns, field, mode, stepN);
        masterUnknowns.copySubVector(mstrUnknwns, k);
        k += mstrUnknwns.giveSize();
    }
}

void ActiveDof :: giveBcValues(FloatArray &masterBcValues, ValueModeType mode, TimeStep *stepN)
{
    FloatArray mstrBcVlus;

    masterBcValues.resize( this->giveNumberOfPrimaryMasterDofs() );
    int countOfMasterDofs = this->giveNumberOfMasterDofs();
    for (int k = 1, i = 1; i <= countOfMasterDofs; i++ ) {
        this->giveMasterDof(i)->giveBcValues(mstrBcVlus, mode, stepN);
        masterBcValues.copySubVector(mstrBcVlus, k);
        k += mstrBcVlus.giveSize();
    }
}

void ActiveDof :: computeDofTransformation(FloatArray &primaryMasterContribs)
{
    FloatArray masterContribution, subPrimaryMasterContribs;
    this->giveActiveBoundaryCondition()->computeDofTransformation(this, masterContribution);

    primaryMasterContribs.resize( this->giveNumberOfPrimaryMasterDofs() );
    int countOfMasterDofs = this->giveNumberOfMasterDofs();
    for (int k = 1, i = 1; i <= countOfMasterDofs; i++ ) {
        this->giveMasterDof(i)->computeDofTransformation(subPrimaryMasterContribs);
        subPrimaryMasterContribs.times( masterContribution.at(i) );
        primaryMasterContribs.copySubVector(subPrimaryMasterContribs, k);
        k += primaryMasterContribs.giveSize();
    }
}

double ActiveDof :: giveUnknown(EquationID type, ValueModeType mode, TimeStep *stepN)
{
    return this->giveActiveBoundaryCondition()->giveUnknown(type, mode, stepN, this);
}

double ActiveDof :: giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *stepN)
{
    return this->giveActiveBoundaryCondition()->giveUnknown(field, mode, stepN, this);
}

contextIOResultType ActiveDof :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    // Nothing here since the boundary condition deals with all the values.
    return CIO_OK;
}

contextIOResultType ActiveDof :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    return CIO_OK;
}

inline Dof *ActiveDof :: giveMasterDof(int i)
{
    return this->giveActiveBoundaryCondition()->giveMasterDof(this, i);
    //return dofManager->giveDomain()->giveDofManager( masterDofMans.at(i) )->giveDofWithID( dofIDs.at(i) );
}

void ActiveDof :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    // No numbering is stored.
}

} // end namespace oofem
