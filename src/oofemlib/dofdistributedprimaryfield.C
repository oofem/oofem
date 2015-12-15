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

#include "dofdistributedprimaryfield.h"
#include "dofmanager.h"
#include "dof.h"
#include "domain.h"
#include "timestep.h"
#include "floatarray.h"
#include "engngm.h"
#include "set.h"
#include "boundarycondition.h"
#include "initialcondition.h"
#include "element.h"
#include "activebc.h"


namespace oofem {
DofDistributedPrimaryField :: DofDistributedPrimaryField(EngngModel *a, int idomain,
                                                         FieldType ft, int nHist) :
    PrimaryField(a, idomain, ft, nHist)
{ }

DofDistributedPrimaryField :: ~DofDistributedPrimaryField()
{ }

double
DofDistributedPrimaryField :: giveUnknownValue(Dof *dof, ValueModeType mode, TimeStep *tStep)
{
#if 0
    double val1 = dof->giveUnknownsDictionaryValue(tStep, VM_Total);
    double val0 = dof->giveUnknownsDictionaryValue(tStep->givePreviousStep(), VM_Total);
    if ( mode == VM_Total ) {
        return this->alpha * val1 + (1.-this->alpha) * val0;
    } else if ( mode == VM_Velocity ) {
        return (val1 - val0) / tStep->giveTimeIncrement();
    } else if ( mode == VM_Incremental ) {
        return val1 - val0;
    } else {
        OOFEM_ERROR("Unknown value mode requested");
        return 0;
    }
#else
    return dof->giveUnknownsDictionaryValue(tStep, mode);
#endif
}

FloatArray *
DofDistributedPrimaryField :: giveSolutionVector(TimeStep *tStep)
{
    //OOFEM_ERROR("This function should not be called. initialize and update should be used instead");
    return PrimaryField :: giveSolutionVector(tStep);
}

void
DofDistributedPrimaryField :: initialize(ValueModeType mode, TimeStep *tStep, FloatArray &answer, const UnknownNumberingScheme &s)
{
    Domain *d = emodel->giveDomain(domainIndx);
    int neq = emodel->giveNumberOfDomainEquations(domainIndx, s);

    answer.resize(neq);
    answer.zero();

    for ( auto &node : d->giveDofManagers() ) {
        for ( Dof *dof: *node ) {
            if ( !dof->isPrimaryDof() ) continue;
            int eqNum = dof->giveEquationNumber(s);
            if ( eqNum ) {
                answer.at(eqNum) = dof->giveUnknownsDictionaryValue(tStep, mode);
            }
        }
    }

    for ( auto &elem : d->giveElements() ) {
        int ndman = elem->giveNumberOfInternalDofManagers();
        for ( int i = 1; i <= ndman; i++ ) {
            for ( auto &dof : *elem->giveInternalDofManager(i) ) {
                if ( !dof->isPrimaryDof() ) continue;
                int eqNum = dof->giveEquationNumber(s);
                if ( eqNum > 0 ) {
                    answer.at(eqNum) = dof->giveUnknownsDictionaryValue(tStep, mode);
                }
            }
        }
    }

    for ( auto &bc : d->giveBcs() ) {
        int ndman = bc->giveNumberOfInternalDofManagers();
        for ( int i = 1; i <= ndman; i++ ) {
            for ( auto &dof : *bc->giveInternalDofManager(i) ) {
                if ( !dof->isPrimaryDof() ) continue;
                int eqNum = dof->giveEquationNumber(s);
                if ( eqNum > 0 ) {
                    answer.at(eqNum) = dof->giveUnknownsDictionaryValue(tStep, mode);
                }
            }
        }
    }
}

// project solutionVector to DoF unknowns dictionary
void
DofDistributedPrimaryField :: update(ValueModeType mode, TimeStep *tStep, const FloatArray &vectorToStore, const UnknownNumberingScheme &s)
{
    Domain *d = emodel->giveDomain(domainIndx);

    for ( auto &node : d->giveDofManagers() ) {
        for ( Dof *dof: *node ) {
            if ( !dof->isPrimaryDof() ) continue;
            int eqNum = dof->giveEquationNumber(s);
            if ( eqNum > 0 ) {
                dof->updateUnknownsDictionary(tStep, mode, vectorToStore.at(eqNum));
            }
            ///@todo This should not be here / Mikael
            if ( mode == VM_Total ) {
                if ( dof->hasBc(tStep) ) {
                    dof->updateUnknownsDictionary(tStep, mode, dof->giveBcValue(VM_Total, tStep));
                }
            }
        }
    }
    
    for ( auto &elem : d->giveElements() ) {
        int ndman = elem->giveNumberOfInternalDofManagers();
        for ( int i = 1; i <= ndman; i++ ) {
            for ( auto &dof : *elem->giveInternalDofManager(i) ) {
                if ( !dof->isPrimaryDof() ) continue;
                int eqNum = dof->giveEquationNumber(s);
                if ( eqNum > 0 ) {
                    dof->updateUnknownsDictionary(tStep, mode, vectorToStore.at(eqNum));
                }
            }
        }
    }

    for ( auto &bc : d->giveBcs() ) {
        int ndman = bc->giveNumberOfInternalDofManagers();
        for ( int i = 1; i <= ndman; i++ ) {
            for ( auto &dof : *bc->giveInternalDofManager(i) ) {
                if ( !dof->isPrimaryDof() ) continue;
                int eqNum = dof->giveEquationNumber(s);
                if ( eqNum > 0 ) {
                    dof->updateUnknownsDictionary(tStep, mode, vectorToStore.at(eqNum));
                }
            }
        }
    }
}


void
DofDistributedPrimaryField :: applyDefaultInitialCondition()
{
    Domain *d = emodel->giveDomain(domainIndx);
    TimeStep *tStep = emodel->giveSolutionStepWhenIcApply();
    // Copy over the old dictionary values to the new step as the initial guess:
    for ( auto &dman : d->giveDofManagers() ) {
        for ( auto &dof : *dman ) {
            dof->updateUnknownsDictionary(tStep, VM_Total, 0.);
            int icid = dof->giveIcId();
            if ( icid > 0 && dof->isPrimaryDof() ) {
                InitialCondition *ic = d->giveIc(icid);
                if ( ic->hasConditionOn(VM_Total) ) {
                    double val = ic->give(VM_Total);
                    dof->updateUnknownsDictionary(tStep, VM_Total, val);
                }
            }
        }
    }

    for ( auto &elem : d->giveElements() ) {
        int ndman = elem->giveNumberOfInternalDofManagers();
        for ( int i = 1; i <= ndman; i++ ) {
            for ( auto &dof : *elem->giveInternalDofManager(i) ) {
                dof->updateUnknownsDictionary(tStep, VM_Total, 0.);
            }
        }
    }

    for ( auto &bc : d->giveBcs() ) {
        int ndman = bc->giveNumberOfInternalDofManagers();
        for ( int i = 1; i <= ndman; i++ ) {
            for ( auto &dof : *bc->giveInternalDofManager(i) ) {
                dof->updateUnknownsDictionary(tStep, VM_Total, 0.);
            }
        }
    }

    // Apply ICs that use sets.
    for ( auto &ic : d->giveIcs() ) {
        this->applyInitialCondition(*ic);
    }
}

void
DofDistributedPrimaryField :: applyInitialCondition(InitialCondition &ic)
{
    if ( ic.giveSetNumber() == 0 ) {
        return;
    }

    Domain *d = ic.giveDomain();
    Set *set = d->giveSet(ic.giveSetNumber());
    TimeStep *tStep = emodel->giveSolutionStepWhenIcApply();
    //TimeStep *prev = tStep->givePreviousStep();

    // We have to set initial value, and velocity, for this particular primary field.
    for ( int inode : set->giveNodeList() ) {
        DofManager *dman = d->giveDofManager(inode);
        double tot0 = 0;
        if ( ic.hasConditionOn(VM_Total) ) {
            tot0 = ic.give(VM_Total);
        }
#if 0
        Not relevant for this time discretization
        double tot1 = 0.
        if ( ic.hasConditionOn(VM_Incremental) ) {
            tot1 = tot0 - ic.give(VM_Incremental);
        } else if ( ic.hasConditionOn(VM_Velocity) ) {
            tot1 = tot0 - ic.give(VM_Velocity) * tStep->giveTimeIncrement();
        } else {
            tot1 = tot0;
        }
#endif
        for ( int dofid : ic.giveDofIDs() ) {
            if ( dman->findDofWithDofId((DofIDItem)dofid) != dman->end() ) {
                Dof *dof = dman->giveDofWithID(dofid);
                //dof->updateUnknownsDictionary(tStep, VM_Total, tot1);
                dof->updateUnknownsDictionary(tStep, VM_Total, tot0);
            }
        }
    }
}


void
DofDistributedPrimaryField :: applyBoundaryCondition(TimeStep *tStep)
{
    Domain *d = emodel->giveDomain(domainIndx);
    for ( auto &dman : d->giveDofManagers() ) {
        for ( auto &dof : *dman ) {
            if ( dof->hasBc(tStep) && dof->isPrimaryDof() ) {
                int bcid = dof->giveBcId();
                double val = static_cast< BoundaryCondition* >(d->giveBc(bcid))->give(dof, VM_Total, tStep->giveTargetTime());
                dof->updateUnknownsDictionary(tStep, VM_Total, val);
            }
        }
    }

    for ( auto &bc : d->giveBcs() ) {
        if ( bc->isImposed(tStep) ) {
            BoundaryCondition *dbc = dynamic_cast< BoundaryCondition* >(bc.get());
            ActiveBoundaryCondition *abc = dynamic_cast< ActiveBoundaryCondition* >(bc.get());
            if ( dbc && dbc->isImposed(tStep) ) {
                this->applyBoundaryCondition(*dbc, tStep);
            } else if ( abc ) {
                for ( auto &dof : *abc->giveInternalDofManager(1) ) {
                    if ( dof->isPrimaryDof() && abc->hasBc(dof, tStep) ) {
                        dof->updateUnknownsDictionary( tStep, VM_Total, abc->giveBcValue(dof, VM_Total, tStep) );
                    }
                }
            }
        }
    }
}


void
DofDistributedPrimaryField :: applyBoundaryCondition(BoundaryCondition &bc, TimeStep *tStep)
{
    if ( bc.giveSetNumber() == 0 ) {
        return;
    }

    Domain *d = bc.giveDomain();
    Set *set = d->giveSet(bc.giveSetNumber());
    for ( int inode : set->giveNodeList() ) {
        DofManager *dman = d->giveDofManager(inode);
        for ( auto &dofid : bc.giveDofIDs() ) {
            if ( !dman->hasDofID((DofIDItem)dofid) ) { ///@todo It's unfortunate that we have to search for the dofid twice.
                continue;
            }
            Dof *dof = dman->giveDofWithID(dofid);
            if ( dof->isPrimaryDof() ) {
                dof->updateUnknownsDictionary( tStep, VM_Total, bc.give(dof, VM_Total, tStep->giveTargetTime()) );
            }
        }
    }
}


void
DofDistributedPrimaryField :: setInitialGuess(DofManager &dman, TimeStep *tStep, TimeStep *prev)
{
    for ( auto &dof : dman ) {
        double val = dof->giveUnknownsDictionaryValue( prev, VM_Total );
        dof->updateUnknownsDictionary( tStep, VM_Total, val );
    }
}

void
DofDistributedPrimaryField :: advanceSolution(TimeStep *tStep)
{
    PrimaryField :: advanceSolution(tStep);
#if 0
    // Copy over the old dictionary values to the new step as the initial guess:
    Domain *d = emodel->giveDomain(1);
    TimeStep *prev = tStep->givePreviousStep();
    for ( auto &dman : d->giveDofManagers() ) {
        this->setInitialGuess(*dman, tStep, prev);
    }

    for ( auto &elem : d->giveElements() ) {
        int ndman = elem->giveNumberOfInternalDofManagers();
        for ( int i = 1; i <= ndman; i++ ) {
            this->setInitialGuess(*elem->giveInternalDofManager(i), tStep, prev);
        }
    }

    for ( auto &bc : d->giveBcs() ) {
        int ndman = bc->giveNumberOfInternalDofManagers();
        for ( int i = 1; i <= ndman; i++ ) {
            this->setInitialGuess(*bc->giveInternalDofManager(i), tStep, prev);
        }
    }

    // Apply dirichlet b.c.s
    //this->applyBoundaryCondition(tStep);
#endif
}


contextIOResultType
DofDistributedPrimaryField :: saveContext(DataStream &stream, ContextMode mode)
{
    // all the job is done by dofs alone
    return CIO_OK;
}

contextIOResultType
DofDistributedPrimaryField :: restoreContext(DataStream &stream, ContextMode mode)
{
    return CIO_OK;
}
} // end namespace oofem
