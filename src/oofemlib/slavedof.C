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

#include "slavedof.h"
#include "node.h"
#include "datastream.h"
#include "contextioerr.h"


#ifndef __MAKEDEPEND
 #include <stdlib.h>
 #include <ctype.h>
#endif

namespace oofem {
SlaveDof :: SlaveDof(int n, DofManager *aNode, DofIDItem id) : Dof(n, aNode, id), masterContribution()
{
    countOfPrimaryMasterDofs = -1;
}


void
SlaveDof :: initialize(int cntOfMstrDfMngr, const IntArray &masterNodes, const IntArray *mstrDofID, const FloatArray &mstrContribution)
{
    int i, id;
    bool idSame = false;


    if ( mstrDofID == NULL ) {
        idSame = true;
    } else
    if ( mstrDofID->giveSize() < cntOfMstrDfMngr ) {
        _error3("initialize: mstrDofID.giveSize %d != cntOfMstrDfMngr %d", mstrDofID->giveSize(), cntOfMstrDfMngr);
    }

    if ( mstrContribution.giveSize() < cntOfMstrDfMngr ) {
        _error3("initialize: mstrContribution.giveSize %d != cntOfMstrDfMngr %d", mstrContribution.giveSize(), cntOfMstrDfMngr);
    }


    countOfMasterDofs  = cntOfMstrDfMngr;
    masterContribution = mstrContribution;

    masterDofMans.resize(countOfMasterDofs);
    dofIDs.resize(countOfMasterDofs);

    for ( i = 1; i <= countOfMasterDofs; i++ ) {
        if ( idSame ) {
            id = this->dofID;
        } else {
            id = mstrDofID->at(i);
        }

        masterDofMans.at(i) = masterNodes.at(i);
        dofIDs.at(i) = id;
    }
}

int
SlaveDof :: giveNumberOfPrimaryMasterDofs(void)
{
    if ( countOfPrimaryMasterDofs > 0 ) {
        return countOfPrimaryMasterDofs;
    } else
    if ( countOfPrimaryMasterDofs == 0 ) {
        _error2( "giveNumberOfPrimaryDofs: slaveDof number %ld is own master", this->giveNumber() );
    }

    countOfPrimaryMasterDofs = 0;

    long i, c = 0;
    for ( i = 1; i <= countOfMasterDofs; i++ ) {
        if ( !this->giveMasterDof(i)->isPrimaryDof() ) {
            c += this->giveMasterDof(i)->giveNumberOfPrimaryMasterDofs();
        } else {
            c += 1;
        }
    }

    return countOfPrimaryMasterDofs = c;
}

void
SlaveDof :: giveUnknowns(FloatArray &masterUnknowns, EquationID type, ValueModeType mode, TimeStep *stepN)
{
    int i, k;
    FloatArray mstrUnknwns;
    Dof *_md;

    masterUnknowns.resize( this->giveNumberOfPrimaryMasterDofs() );

    for ( k = 1, i = 1; i <= countOfMasterDofs; i++ ) {
        _md = this->giveMasterDof(i);
        if ( !_md->isPrimaryDof() ) {
            _md->giveUnknowns(mstrUnknwns, type, mode, stepN);
            masterUnknowns.copySubVector(mstrUnknwns, k);
            k += mstrUnknwns.giveSize();
        } else {
            masterUnknowns.at(k++) = _md->giveUnknown(type, mode, stepN);
        }
    }
}

void
SlaveDof :: giveUnknowns(FloatArray &masterUnknowns, PrimaryField &field, ValueModeType mode, TimeStep *stepN)
{
    int i, k;
    FloatArray mstrUnknwns;
    Dof *_md;

    masterUnknowns.resize( this->giveNumberOfPrimaryMasterDofs() );

    for ( k = 1, i = 1; i <= countOfMasterDofs; i++ ) {
        _md = this->giveMasterDof(i);
        if ( !_md->isPrimaryDof() ) {
            _md->giveUnknowns(mstrUnknwns, field, mode, stepN);
            masterUnknowns.copySubVector(mstrUnknwns, k);
            k += mstrUnknwns.giveSize();
        } else {
            masterUnknowns.at(k++) = _md->giveUnknown(field, mode, stepN);
        }
    }
}

void
SlaveDof :: giveBcValues(FloatArray &masterBcValues, ValueModeType mode, TimeStep *stepN)
{
    int i, k;
    FloatArray mstrBcVlus;
    Dof *idof;

    masterBcValues.resize( this->giveNumberOfPrimaryMasterDofs() );

    for ( k = 1, i = 1; i <= countOfMasterDofs; i++ ) {
        idof = this->giveMasterDof(i);
        if ( !idof->isPrimaryDof() ) {
            idof->giveBcValues(mstrBcVlus, mode, stepN);
            masterBcValues.copySubVector(mstrBcVlus, k);
            k += mstrBcVlus.giveSize();
        } else {
            if ( idof->hasBc(stepN) ) {
                masterBcValues.at(k++) = idof->giveBcValue(mode, stepN);
            } else {
                masterBcValues.at(k++) = 0.0;
            }
        }
    }
}

void
SlaveDof :: computeDofTransformation(FloatArray &masterContribs)
{
    int i, k;
    FloatArray mstrContrs;
    Dof *idof;

    masterContribs.resize( this->giveNumberOfPrimaryMasterDofs() );

    for ( k = 1, i = 1; i <= countOfMasterDofs; i++ ) {
        idof = this->giveMasterDof(i);
        if ( !idof->isPrimaryDof() ) {
            idof->computeDofTransformation(mstrContrs);
            mstrContrs.times( masterContribution.at(i) );
            masterContribs.copySubVector(mstrContrs, k);
            k += mstrContrs.giveSize();
        } else {
            masterContribs.at(k++) = masterContribution.at(i);
        }
    }
}

void
SlaveDof :: giveEquationNumbers(IntArray &masterEqNumbers, const UnknownNumberingScheme &s)
{
    int i, k;
    IntArray mstrEqNmbrs;
    Dof *idof;

    masterEqNumbers.resize( this->giveNumberOfPrimaryMasterDofs() );

    for ( k = 1, i = 1; i <= countOfMasterDofs; i++ ) {
        idof = this->giveMasterDof(i);
        if ( !idof->isPrimaryDof() ) {
            idof->giveEquationNumbers(mstrEqNmbrs, s);
            masterEqNumbers.copySubVector(mstrEqNmbrs, k);
            k += mstrEqNmbrs.giveSize();
        } else {
            masterEqNumbers.at(k++) = idof->giveEquationNumber(s);
        }
    }
}


double SlaveDof :: giveUnknown(EquationID type, ValueModeType mode, TimeStep *stepN)
{
    FloatArray masterUnknowns, t;

    giveUnknowns(masterUnknowns, type, mode, stepN);
    computeDofTransformation(t);

    return masterUnknowns.dotProduct(t);
}

double SlaveDof :: giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *stepN)
{
    FloatArray masterUnknowns, t;

    giveUnknowns(masterUnknowns, field, mode, stepN);
    computeDofTransformation(t);

    return masterUnknowns.dotProduct(t);
}


contextIOResultType SlaveDof :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( ( iores = Dof :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( mode & CM_Definition ) {
        if ( !stream->write(& countOfMasterDofs, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream->write(& countOfPrimaryMasterDofs, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( ( iores = masterContribution.storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        int _idof, _idofmanNum;
        for ( _idof = 1; _idof <= countOfMasterDofs; _idof++ ) {
#ifdef __PARALLEL_MODE
            if ( mode & CM_DefinitionGlobal ) {
                _idofmanNum = dofManager->giveDomain()->giveDofManager( masterDofMans.at(_idof) )->giveGlobalNumber();
            } else {
                _idofmanNum = masterDofMans.at(_idof);
            }

#else
            _idofmanNum = masterDofMans.at(_idof);
#endif

            if ( !stream->write(& _idofmanNum, 1) ) {
                THROW_CIOERR(CIO_IOERR);
            }
        }

        if ( ( iores = dofIDs.storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    } // if ( mode & CM_Definition )

    return CIO_OK;
}

contextIOResultType SlaveDof :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( ( iores = Dof :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( mode & CM_Definition ) {
        if ( !stream->read(& countOfMasterDofs, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream->read(& countOfPrimaryMasterDofs, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( ( iores = masterContribution.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        int _idof;

        masterDofMans.resize(countOfMasterDofs);
        for ( _idof = 1; _idof <= countOfMasterDofs; _idof++ ) {
            if ( !stream->read(& masterDofMans.at(_idof), 1) ) {
                THROW_CIOERR(CIO_IOERR);
            }
        }

        if ( ( iores = dofIDs.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    } // if ( mode & CM_Definition )

    return CIO_OK;
}


inline Dof *
SlaveDof :: giveMasterDof(int i) {
    return dofManager->giveDomain()->giveDofManager( masterDofMans.at(i) )->giveDofWithID( dofIDs.at(i) );
}


void
SlaveDof :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    int i;
    for ( i = 1; i <= countOfMasterDofs; i++ ) {
        masterDofMans.at(i) = f(masterDofMans.at(i), ERS_DofManager);
    }
}


void
SlaveDof :: giveMasterDofManArray(IntArray &answer)
{
    long i, k;
    IntArray mstrDofManArry;

    answer.resize( this->giveNumberOfPrimaryMasterDofs() );

    for ( k = 1, i = 1; i <= countOfMasterDofs; i++ ) {
        if ( !giveMasterDof(i)->isPrimaryDof() ) {
            this->giveMasterDof(i)->giveMasterDofManArray(mstrDofManArry);
            answer.copySubVector(mstrDofManArry, k);
            k += mstrDofManArry.giveSize();
        } else   {
            answer.at(k++) = this->giveMasterDof(i)->giveDofManNumber();
        }
    }
}
} // end namespace oofem
