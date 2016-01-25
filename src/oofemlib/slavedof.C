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

#include "slavedof.h"
#include "domain.h"
#include "dofmanager.h"
#include "error.h"
#include "datastream.h"
#include "entityrenumberingscheme.h"
#include "contextioerr.h"

namespace oofem {
SlaveDof :: SlaveDof(DofManager *aNode, DofIDItem id) : Dof(aNode, id), masterContribution()
{
    countOfPrimaryMasterDofs = -1;
    countOfMasterDofs = -1;
}


void
SlaveDof :: initialize(const IntArray &masterNodes, const IntArray &mstrDofID, const FloatArray &mstrContribution)
{
    int id;
    bool idSame = false;


    if ( mstrDofID.isEmpty() ) {
        idSame = true;
    } else if ( mstrDofID.giveSize() != mstrContribution.giveSize() ) {
        OOFEM_ERROR("mstrDofID.giveSize %d != mstrContribution.giveSize %d", mstrDofID.giveSize(), mstrContribution.giveSize());
    }


    countOfMasterDofs  = mstrContribution.giveSize();
    masterContribution = mstrContribution;

    masterDofMans.resize(countOfMasterDofs);
    dofIDs.resize(countOfMasterDofs);

    for ( int i = 1; i <= countOfMasterDofs; i++ ) {
        if ( idSame ) {
            id = this->dofID;
        } else {
            id = mstrDofID.at(i);
        }

        masterDofMans.at(i) = masterNodes.at(i);
        dofIDs.at(i) = id;
    }
}

int
SlaveDof :: giveNumberOfPrimaryMasterDofs()
{
    if ( countOfPrimaryMasterDofs > 0 ) {
        return countOfPrimaryMasterDofs;
    } else
    if ( countOfPrimaryMasterDofs == 0 ) {
        OOFEM_ERROR("slaveDof is own master");
    }

    countOfPrimaryMasterDofs = 0;

    int c = 0;
    for ( int i = 1; i <= countOfMasterDofs; i++ ) {
        c += this->giveMasterDof(i)->giveNumberOfPrimaryMasterDofs();
    }

    return countOfPrimaryMasterDofs = c;
}

void
SlaveDof :: giveUnknowns(FloatArray &masterUnknowns, ValueModeType mode, TimeStep *tStep)
{
    FloatArray mstrUnknwns;

    masterUnknowns.resize( this->giveNumberOfPrimaryMasterDofs() );

    for ( int k = 1, i = 1; i <= countOfMasterDofs; i++ ) {
        this->giveMasterDof(i)->giveUnknowns(mstrUnknwns, mode, tStep);
        masterUnknowns.copySubVector(mstrUnknwns, k);
        k += mstrUnknwns.giveSize();
    }
}

void
SlaveDof :: giveUnknowns(FloatArray &masterUnknowns, PrimaryField &field, ValueModeType mode, TimeStep *tStep)
{
    FloatArray mstrUnknwns;

    masterUnknowns.resize( this->giveNumberOfPrimaryMasterDofs() );

    for ( int k = 1, i = 1; i <= countOfMasterDofs; i++ ) {
        this->giveMasterDof(i)->giveUnknowns(mstrUnknwns, field, mode, tStep);
        masterUnknowns.copySubVector(mstrUnknwns, k);
        k += mstrUnknwns.giveSize();
    }
}

void
SlaveDof :: computeDofTransformation(FloatArray &primaryMasterContribs)
{
    FloatArray subPrimaryMasterContribs;

    primaryMasterContribs.resize( this->giveNumberOfPrimaryMasterDofs() );

    for ( int k = 1, i = 1; i <= countOfMasterDofs; i++ ) {
        this->giveMasterDof(i)->computeDofTransformation(subPrimaryMasterContribs);
        subPrimaryMasterContribs.times( masterContribution.at(i) );
        primaryMasterContribs.copySubVector(subPrimaryMasterContribs, k);
        k += subPrimaryMasterContribs.giveSize();
    }
}

void
SlaveDof :: giveEquationNumbers(IntArray &masterEqNumbers, const UnknownNumberingScheme &s)
{
    IntArray mstrEqNmbrs;

    int masterDofs = this->giveNumberOfPrimaryMasterDofs();
    masterEqNumbers.preallocate(masterDofs);
    masterEqNumbers.clear();

    for ( int i = 1; i <= countOfMasterDofs; i++ ) {
        this->giveMasterDof(i)->giveEquationNumbers(mstrEqNmbrs, s);
        masterEqNumbers.followedBy(mstrEqNmbrs);
    }
}


void
SlaveDof :: giveDofIDs(IntArray &masterDofIDs)
{
    IntArray temp;

    int masterDofs = this->giveNumberOfPrimaryMasterDofs();
    masterDofIDs.preallocate(masterDofs);
    masterDofIDs.clear();

    for ( int i = 1; i <= countOfMasterDofs; i++ ) {
        this->giveMasterDof(i)->giveDofIDs(temp);
        masterDofIDs.followedBy(temp);
    }
}


double SlaveDof :: giveUnknown(ValueModeType mode, TimeStep *tStep)
{
    FloatArray masterUnknowns, t;

    this->giveUnknowns(masterUnknowns, mode, tStep);
    this->computeDofTransformation(t);

    return masterUnknowns.dotProduct(t);
}

double SlaveDof :: giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *tStep)
{
    FloatArray masterUnknowns, t;

    giveUnknowns(masterUnknowns, field, mode, tStep);
    computeDofTransformation(t);

    return masterUnknowns.dotProduct(t);
}


int SlaveDof :: __giveEquationNumber() const
{
    OOFEM_ERROR("undefined");
    return 0;
}


int SlaveDof :: __givePrescribedEquationNumber()
{
    OOFEM_ERROR("undefined");
    return 0;
}

contextIOResultType SlaveDof :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( ( iores = Dof :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( mode & CM_Definition ) {
        if ( !stream.write(countOfMasterDofs) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.write(countOfPrimaryMasterDofs) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( ( iores = masterContribution.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        for ( int _idof = 1; _idof <= countOfMasterDofs; _idof++ ) {
            int _idofmanNum;
            if ( mode & CM_DefinitionGlobal ) {
                _idofmanNum = dofManager->giveDomain()->giveDofManager( masterDofMans.at(_idof) )->giveGlobalNumber();
            } else {
                _idofmanNum = masterDofMans.at(_idof);
            }

            if ( !stream.write(_idofmanNum) ) {
                THROW_CIOERR(CIO_IOERR);
            }
        }

        if ( ( iores = dofIDs.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    } // if ( mode & CM_Definition )

    return CIO_OK;
}

contextIOResultType SlaveDof :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( ( iores = Dof :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( mode & CM_Definition ) {
        if ( !stream.read(countOfMasterDofs) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.read(countOfPrimaryMasterDofs) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( ( iores = masterContribution.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        int _idof;

        masterDofMans.resize(countOfMasterDofs);
        for ( _idof = 1; _idof <= countOfMasterDofs; _idof++ ) {
            if ( !stream.read(masterDofMans.at(_idof)) ) {
                THROW_CIOERR(CIO_IOERR);
            }
        }

        if ( ( iores = dofIDs.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    } // if ( mode & CM_Definition )

    return CIO_OK;
}


inline Dof *
SlaveDof :: giveMasterDof(int i)
{
    return dofManager->giveDomain()->giveDofManager( masterDofMans.at(i) )->giveDofWithID( dofIDs.at(i) );
}


void
SlaveDof :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    for ( int i = 1; i <= countOfMasterDofs; i++ ) {
        masterDofMans.at(i) = f(masterDofMans.at(i), ERS_DofManager);
    }
}


void
SlaveDof :: giveMasterDofManArray(IntArray &answer)
{
    IntArray mstrDofManArry;

    int masterDofs = this->giveNumberOfPrimaryMasterDofs();
    answer.preallocate(masterDofs);
    answer.clear();

    for ( int i = 1; i <= countOfMasterDofs; i++ ) {
        this->giveMasterDof(i)->giveMasterDofManArray(mstrDofManArry);
        answer.followedBy(mstrDofManArry);
    }
}
} // end namespace oofem
