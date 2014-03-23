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

#include "primaryfield.h"
#include "spatiallocalizer.h"
#include "dofmanager.h"
#include "dof.h"
#include "element.h"
#include "timestep.h"
#include "datastream.h"
#include "contextioerr.h"

namespace oofem {
PrimaryField :: PrimaryField(EngngModel *a, int idomain,
                             FieldType ft, EquationID ut, int nHist) : Field(ft), solutionVectors(nHist + 1), solStepList(nHist + 1)
{
    this->actualStepNumber = -999;
    this->actualStepIndx = 0;
    this->nHistVectors = nHist;
    this->ut = ut;

    emodel = a;
    domainIndx = idomain;
}

PrimaryField :: ~PrimaryField()
{ }

void
PrimaryField :: initialize(ValueModeType mode, TimeStep *tStep, FloatArray &answer, const UnknownNumberingScheme &s)
{
    int neq =  emodel->giveNumberOfDomainEquations(domainIndx, s);
    answer.resize(neq);
    answer.zero();

    if ( mode == VM_Total ) {
        answer = * ( this->giveSolutionVector(tStep) );
    } else if ( mode == VM_Incremental ) {
        int indxm1 = this->resolveIndx(tStep, -1);
        answer = * ( this->giveSolutionVector(tStep) );
        answer.subtract( * this->giveSolutionVector(indxm1) );
    } else {
        OOFEM_ERROR("unsupported mode %s", __ValueModeTypeToString(mode));
    }
}


double
PrimaryField :: giveUnknownValue(Dof *dof, ValueModeType mode, TimeStep *tStep)
{
    int eq = dof->giveEquationNumber( emodel->giveUnknownNumberingScheme(this->ut) );
    if ( eq == 0 ) {
        OOFEM_ERROR("invalid equation number");
    }

    if ( mode == VM_Total ) {
        return this->giveSolutionVector(tStep)->at(eq);
    } else if ( mode == VM_Incremental ) {
        int indxm1 = this->resolveIndx(tStep, -1);
        return ( this->giveSolutionVector(tStep)->at(eq) - this->giveSolutionVector(indxm1)->at(eq) );
    } else {
        OOFEM_ERROR("unsupported mode");
    }

    return 0.0;
}

int
PrimaryField :: __evaluateAt(FloatArray &answer, DofManager *dman,
                             ValueModeType mode, TimeStep *tStep,
                             IntArray *dofId)
{
    if ( dman->giveDomain() == this->emodel->giveDomain(domainIndx) ) {
        if ( dofId ) {
            dman->giveUnknownVector(answer, * dofId, * this, mode, tStep);
            return 0; // ok
        } else { // all dofs requested
            dman->giveCompleteUnknownVector(answer, mode, tStep);
            return 0; // ok
        }
    } else {
        return this->__evaluateAt(answer, * dman->giveCoordinates(),
                                  mode, tStep, dofId);
    }
}


int
PrimaryField :: __evaluateAt(FloatArray &answer, FloatArray &coords,
                             ValueModeType mode, TimeStep *tStep,
                             IntArray *dofId)
{
    Element *bgelem;
    Domain *domain = emodel->giveDomain(domainIndx);
    SpatialLocalizer *sl = domain->giveSpatialLocalizer();
    // locate background element
    if ( ( bgelem = sl->giveElementContainingPoint(coords) ) == NULL ) {
        //_error("PrimaryField::evaluateAt: point not found in domain\n");
        return 1;
    }

    EIPrimaryFieldInterface *interface = static_cast< EIPrimaryFieldInterface * >( bgelem->giveInterface(EIPrimaryFieldInterfaceType) );
    if ( interface ) {
        if ( dofId ) {
            return interface->EIPrimaryFieldI_evaluateFieldVectorAt(answer, * this, coords, * dofId, mode, tStep);
        } else { // use element default dof id mask
            IntArray elemDofId;
            bgelem->giveElementDofIDMask(this->giveEquationID(), elemDofId);
            return interface->EIPrimaryFieldI_evaluateFieldVectorAt(answer, * this, coords, elemDofId, mode, tStep);
        }
    } else {
        OOFEM_ERROR("background element does not support EIPrimaryFiledInterface\n");
        return 1; // failed
    }
}

int
PrimaryField :: evaluateAt(FloatArray &answer, FloatArray &coords,
                           ValueModeType mode, TimeStep *tStep)
{
    return this->__evaluateAt(answer, coords, mode, tStep, NULL);
}


int
PrimaryField :: evaluateAt(FloatArray &answer, DofManager *dman,
                           ValueModeType mode, TimeStep *tStep)
{
    return this->__evaluateAt(answer, dman, mode, tStep, NULL);
}


FloatArray *
PrimaryField :: giveSolutionVector(TimeStep *tStep)
{
    return this->giveSolutionVector( resolveIndx(tStep, 0) );
}

FloatArray *
PrimaryField :: giveSolutionVector(int i)
{
    FloatArray *answer = NULL;
    if ( ( i >= 1 ) && ( i <= ( nHistVectors + 1 ) ) ) {
        answer = &solutionVectors[i-1];
    } else {
        OOFEM_ERROR("index out of range");
    }

    return answer;
}



int
PrimaryField :: resolveIndx(TimeStep *tStep, int shift)
{
    int tStepo = tStep->giveNumber();
    int relPos = actualStepNumber - tStepo - shift;
    if ( ( relPos >= 0 ) && ( relPos <= nHistVectors ) ) {
        return ( actualStepIndx + relPos ) % ( nHistVectors + 1 ) + 1;
    } else {
        OOFEM_ERROR("History not available for relative step no. %d to step no. %d", shift, tStepo);
    }

    return 0;
}


void
PrimaryField :: advanceSolution(TimeStep *tStep)
{
    TimeStep *newts;
    if ( actualStepNumber == tStep->giveNumber() ) {
        // We're one the correct step already, no need to advance.
        return;
    }
    if ( ( actualStepNumber >= 0 ) && ( actualStepNumber + 1 != tStep->giveNumber() ) ) {
        OOFEM_ERROR("can not advance due to steps skipped");
    }

    actualStepIndx = ( actualStepIndx > 0 ) ? actualStepIndx - 1 : nHistVectors;
    actualStepNumber = tStep->giveNumber();
    if ( ( newts = solStepList.at(actualStepIndx + 1) ) ) {
        * newts = * tStep;
    } else {
        solStepList.put( actualStepIndx + 1, new TimeStep(* tStep) );
    }
}


contextIOResultType
PrimaryField :: saveContext(DataStream *stream, ContextMode mode)
{
    contextIOResultType iores(CIO_IOERR);

    if ( !stream->write(& actualStepNumber, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->write(& actualStepIndx, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    for ( int i = 0; i <= nHistVectors; i++ ) {
        if ( ( iores = solutionVectors[i].storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    TimeStep *iStep;
    int flag;
    for ( int i = 0; i <= nHistVectors; i++ ) {
        if ( ( iStep = solStepList.at(i + 1) ) ) {
            flag = 1;
        } else {
            flag = 0;
        }

        if ( !stream->write(& flag, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( flag ) {
            if ( ( iores = solStepList.at(i + 1)->saveContext(stream, mode) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }
        }
    }

    return CIO_OK;
}

contextIOResultType
PrimaryField :: restoreContext(DataStream *stream, ContextMode mode)
{
    contextIOResultType iores(CIO_IOERR);

    if ( !stream->read(& actualStepNumber, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream->read(& actualStepIndx, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    for ( int i = 0; i <= nHistVectors; i++ ) {
        if ( ( iores = solutionVectors[i].restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    int flag;
    TimeStep *iStep;
    for ( int i = 0; i <= nHistVectors; i++ ) {
        if ( !stream->read(& flag, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( flag ) {
            iStep = new TimeStep(emodel);
            if ( ( iores = iStep->restoreContext(stream, mode) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }
        } else {
            iStep = NULL;
        }

        solStepList.put(i + 1, iStep);
    }

    return CIO_OK;
}
} // end namespace oofem
