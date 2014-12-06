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
#include "engngm.h"
#include "initialcondition.h"
#include "boundarycondition.h"

namespace oofem {
PrimaryField :: PrimaryField(EngngModel *a, int idomain,
                             FieldType ft, int nHist) : Field(ft), solutionVectors(nHist + 1), prescribedVectors(nHist + 1), solStepList(nHist + 1, a)
{
    this->actualStepNumber = -999;
    this->actualStepIndx = 0;
    this->nHistVectors = nHist;

    emodel = a;
    domainIndx = idomain;
}

PrimaryField :: ~PrimaryField()
{ }


void
PrimaryField :: storeDofManager(TimeStep *tStep, DofManager &dman)
{
    for ( Dof *dof: dman ) {
        int eq = dof->giveEqn();
        TimeStep *step = tStep;
        for ( int hist = 0; hist < this->nHistVectors; ++hist ) {
            if ( eq > 0 ) {
                FloatArray *vec = this->giveSolutionVector( resolveIndx(tStep, 0) );
                dof->updateUnknownsDictionary(step, VM_Total, vec->at(eq));
            } else if ( eq < 0 ) {
                FloatArray *vec = this->givePrescribedVector( resolveIndx(tStep, 0) );
                dof->updateUnknownsDictionary(step, VM_Total, vec->at(-eq));
            }
            step = tStep->givePreviousStep();
        }
    }
}

void
PrimaryField :: storeInDofDictionaries(TimeStep *tStep, Domain *d)
{
    for ( auto &dman : d->giveDofManagers() ) {
        this->storeDofManager(tStep, *dman);
    }

    for ( auto &elem : d->giveElements() ) {
        int ndman = elem->giveNumberOfInternalDofManagers();
        for ( int i = 1; i <= ndman; i++ ) {
            this->storeDofManager(tStep, *elem->giveInternalDofManager(i));
        }
    }

    for ( auto &bc : d->giveBcs() ) {
        int ndman = bc->giveNumberOfInternalDofManagers();
        for ( int i = 1; i <= ndman; i++ ) {
            this->storeDofManager(tStep, *bc->giveInternalDofManager(i));
        }
    }
}


void
PrimaryField :: readDofManager(TimeStep *tStep, DofManager &dman)
{
    for ( Dof *dof: dman ) {
        int eq = dof->giveEqn();
        TimeStep *step = tStep;
        for ( int hist = 0; hist < this->nHistVectors; ++hist ) {
            if ( eq > 0 ) {
                FloatArray *vec = this->giveSolutionVector( resolveIndx(tStep, 0));
                vec->at(eq) = dof->giveUnknownsDictionaryValue(step, VM_Total);
            } else if ( eq < 0 ) {
                FloatArray *vec = this->givePrescribedVector( resolveIndx(tStep, 0) );
                vec->at(-eq) = dof->giveUnknownsDictionaryValue(step, VM_Total);
            }
            step = tStep->givePreviousStep();
        }
    }
}


void
PrimaryField :: readFromDofDictionaries(TimeStep *tStep, Domain *d)
{
    for ( auto &dman : d->giveDofManagers() ) {
        this->readDofManager(tStep, *dman);
    }

    for ( auto &elem : d->giveElements() ) {
        int ndman = elem->giveNumberOfInternalDofManagers();
        for ( int i = 1; i <= ndman; i++ ) {
            this->readDofManager(tStep, *elem->giveInternalDofManager(i));
        }
    }

    for ( auto &bc : d->giveBcs() ) {
        int ndman = bc->giveNumberOfInternalDofManagers();
        for ( int i = 1; i <= ndman; i++ ) {
            this->readDofManager(tStep, *bc->giveInternalDofManager(i));
        }
    }
}


void
PrimaryField :: initialize(ValueModeType mode, TimeStep *tStep, FloatArray &answer, const UnknownNumberingScheme &s)
{
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

void
PrimaryField :: applyDefaultInitialCondition(int domain)
{
    int neq = emodel->giveNumberOfDomainEquations( domain, EModelDefaultEquationNumbering() );
    int npeq = emodel->giveNumberOfDomainEquations( domain, EModelDefaultPrescribedEquationNumbering() );
    for ( auto &s : solutionVectors ) {
        s.resize(neq);
        s.zero();
    }
    for ( auto &s : prescribedVectors ) {
        s.resize(npeq);
        s.zero();
    }
}


void
PrimaryField :: applyInitialCondition(InitialCondition &ic)
{
    IntArray loc_s, loc_ps;
    Domain *d = ic.giveDomain();
    Set *set = d->giveSet(ic.giveSetNumber());
    TimeStep *tStep = emodel->giveSolutionStepWhenIcApply();
    int indxm0 = this->resolveIndx(tStep, 0);
    int indxm1 = this->resolveIndx(tStep, -1);
    FloatArray *f0 = this->giveSolutionVector(indxm0);
    FloatArray *f1 = this->giveSolutionVector(indxm1);
    FloatArray *p0 = this->givePrescribedVector(indxm0);
    FloatArray *p1 = this->givePrescribedVector(indxm1);

    // We have to set initial value, and velocity, for this particular primary field.
    for ( int inode : set->giveNodeList() ) {
        DofManager *dman = d->giveDofManager(inode);
        double tot0 = 0., tot1 = 0.;
        if ( ic.hasConditionOn(VM_Total) ) {
            tot0 = ic.give(VM_Total);
        }
        if ( ic.hasConditionOn(VM_Incremental) ) {
            tot1 = tot0 - ic.give(VM_Incremental);
        } else if ( ic.hasConditionOn(VM_Velocity) ) {
            tot1 = tot0 - ic.give(VM_Velocity) * tStep->giveTimeIncrement();
        } else {
            tot1 = tot0;
        }
        for ( auto &dof : *dman ) {
            int eq = dof->giveEqn();
            if ( eq > 0 ) {
                f0->at(eq) = tot0;
                f1->at(eq) = tot1;
            } else if ( eq < 0 ) {
                p0->at(-eq) = tot0;
                p1->at(-eq) = tot1;
            }
        }
    }
}


void
PrimaryField :: applyBoundaryCondition(TimeStep *tStep, int di)
{
    Domain *d = emodel->giveDomain(di);
    FloatArray *f = this->givePrescribedVector(resolveIndx(tStep, 0));
    for ( auto &dman : d->giveDofManagers() ) {
        for ( auto &dof : *dman ) {
            int peq = - dof->giveEqn();
            if ( peq > 0 ) {
                f->at(peq) = dof->giveBcValue(VM_Total, tStep);
            }
        }
    }
}


void
PrimaryField :: applyBoundaryCondition(BoundaryCondition &bc, TimeStep *tStep)
{
    if ( bc.giveSetNumber() == 0 ) {
        return;
    }

    Domain *d = bc.giveDomain();
    Set *set = d->giveSet(bc.giveSetNumber());
    FloatArray *f = this->givePrescribedVector(resolveIndx(tStep, 0));
    for ( int inode : set->giveNodeList() ) {
        DofManager *dman = d->giveDofManager(inode);
        for ( auto &dofid : bc.giveDofIDs() ) {
            Dof *dof = dman->giveDofWithID(dofid);
            int peq = - dof->giveEqn(); // Note, only consider prescribed equations here
            if ( peq > 0 ) {
                f->at(peq) = bc.give(dof, VM_Total, tStep);
            }
        }
    }
}


void
PrimaryField :: update(ValueModeType mode, TimeStep *tStep, FloatArray &vectorToStore)
{
    if ( mode == VM_Total ) {
        * this->giveSolutionVector(tStep) = vectorToStore;
    } else {
        OOFEM_ERROR("unsupported mode %s", __ValueModeTypeToString(mode));
    }
}

double
PrimaryField :: giveUnknownValue(Dof *dof, ValueModeType mode, TimeStep *tStep)
{
    int eq = dof->giveEqn();
    if ( eq == 0 ) {
        OOFEM_ERROR("invalid equation number (slave dof maybe?)");
    }

    if ( mode == VM_Total ) {
        int indxm0 = this->resolveIndx(tStep, 0);
        if ( eq > 0 )
            return this->giveSolutionVector(indxm0)->at(eq);
        else
            return this->givePrescribedVector(indxm0)->at(eq);
    } else if ( mode == VM_Incremental ) {
        int indxm0 = this->resolveIndx(tStep, 0);
        int indxm1 = this->resolveIndx(tStep, -1);
        if ( this->giveSolutionVector(indxm1)->giveSize() == 0 ) ///@todo Clean this up, this is a hack for when we ask for the increment in the first step.
            this->giveSolutionVector(indxm1)->resize(this->giveSolutionVector(indxm0)->giveSize());
        if ( eq > 0 )
            return ( this->giveSolutionVector(indxm0)->at(eq) - this->giveSolutionVector(indxm1)->at(eq) );
        else
            return ( this->givePrescribedVector(indxm0)->at(eq) - this->givePrescribedVector(indxm1)->at(eq) );
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
            bgelem->giveElementDofIDMask(elemDofId);
            return interface->EIPrimaryFieldI_evaluateFieldVectorAt(answer, * this, coords, elemDofId, mode, tStep);
        }
    } else {
        OOFEM_ERROR("background element does not support EIPrimaryFiledInterface");
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
    if ( ( i < 1 ) || ( i > ( nHistVectors + 1 ) ) ) {
        OOFEM_ERROR("index out of range");
    }
    return &solutionVectors[i-1];
}

FloatArray *
PrimaryField :: givePrescribedVector(int i)
{
    if ( ( i < 1 ) || ( i > ( nHistVectors + 1 ) ) ) {
        OOFEM_ERROR("index out of range");
    }
    return &prescribedVectors[i-1];
}


int
PrimaryField :: resolveIndx(TimeStep *tStep, int shift)
{
    int tStepo = tStep->giveNumber();
    int relPos = actualStepNumber - tStepo - shift;
    if ( ( relPos >= 0 ) && ( relPos <= nHistVectors ) ) {
        return ( actualStepIndx + relPos ) % ( nHistVectors + 1 ) + 1;
    } else {
        OOFEM_ERROR("History not available for relative step no. %d to step no. %d (actualStepNumber = %d)", shift, tStepo, actualStepNumber);
    }

    return 0;
}


void
PrimaryField :: advanceSolution(TimeStep *tStep)
{
    if ( actualStepNumber == tStep->giveNumber() ) {
        // We're one the correct step already, no need to advance.
        return;
    }
    if ( ( actualStepNumber >= 0 ) && ( actualStepNumber + 1 != tStep->giveNumber() ) ) {
        OOFEM_ERROR("can not advance due to steps skipped");
    }

    actualStepIndx = ( actualStepIndx > 0 ) ? actualStepIndx - 1 : nHistVectors;
    actualStepNumber = tStep->giveNumber();
    solStepList[actualStepIndx] = * tStep;
    if ( nHistVectors > 1 ) {
        // Copy over the old status to the new nodes
        int curr = this->resolveIndx(tStep, 0);
        int prev = this->resolveIndx(tStep, -1);
        *this->giveSolutionVector(curr) = *this->giveSolutionVector(prev);
        *this->givePrescribedVector(curr) = *this->givePrescribedVector(prev);
    }
}


contextIOResultType
PrimaryField :: saveContext(DataStream &stream, ContextMode mode)
{
    contextIOResultType iores(CIO_IOERR);

    if ( !stream.write(actualStepNumber) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.write(actualStepIndx) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    for ( const auto &vec : solutionVectors ) {
        if ( ( iores = vec.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    for ( const auto &vec : prescribedVectors ) {
        if ( ( iores = vec.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    for ( auto &step : solStepList ) {
        if ( ( iores = step.saveContext(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}

contextIOResultType
PrimaryField :: restoreContext(DataStream &stream, ContextMode mode)
{
    contextIOResultType iores(CIO_IOERR);

    if ( !stream.read(actualStepNumber) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !stream.read(actualStepIndx) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    for ( auto &vec : solutionVectors ) {
        if ( ( iores = vec.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    for ( auto &vec : prescribedVectors ) {
        if ( ( iores = vec.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    for ( int i = 0; i <= nHistVectors; i++ ) {
        solStepList[i] = TimeStep(emodel);
        if ( ( iores = solStepList[i].restoreContext(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}
} // end namespace oofem
