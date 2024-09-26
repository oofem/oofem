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

#include "dg.h"
#include "timestep.h"
#include "element.h"
#include "dofmanager.h"
#include "dof.h"
#include "dictionary.h"
#include "verbose.h"
#include "classfactory.h"
#include "mathfem.h"
#include "assemblercallback.h"
#include "unknownnumberingscheme.h"
#include "dofdistributedprimaryfield.h"
#include "primaryfield.h"
#include "maskedprimaryfield.h"
#include "nrsolver.h"
#include "activebc.h"
#include "boundarycondition.h"
#include "boundaryload.h"
#include "outputmanager.h"
#include "connectivitytable.h"
#include "mpm.h"

namespace oofem {
REGISTER_EngngModel(DGProblem);

ScalarAdvectionLhsAssembler :: ScalarAdvectionLhsAssembler(double alpha, double deltaT, Variable::VariableQuantity q) : 
    MatrixAssembler(), alpha(alpha), deltaT(deltaT)
{
    this->q = q;
}


void ScalarAdvectionLhsAssembler :: matrixFromElement(FloatMatrix &answer, Element &el, TimeStep *tStep) const
{
    FloatMatrix contrib;
    IntArray loc;
    MPElement *e = dynamic_cast<MPElement*>(&el);
    int ndofs = e->giveNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();

    e->getLocalCodeNumbers (loc, q);

    e->giveCharacteristicMatrix(contrib, MassMatrix, tStep);
    if (contrib.isNotEmpty()) {
        contrib.times(this->deltaT/2.0);
        answer.assemble(contrib, loc, loc);
    }
    e->giveCharacteristicMatrix(contrib, StiffnessMatrix, tStep);
    if (contrib.isNotEmpty()) {
        answer.assemble(contrib, loc, loc);
    }
    // boundary terms
    e->giveCharacteristicMatrix(contrib, InternalFluxVector, tStep);
    if (contrib.isNotEmpty()) {
        answer.assemble(contrib, loc, loc);
    }
}

ScalarAdvectionRhsAssembler :: ScalarAdvectionRhsAssembler(double alpha, double deltaT, Variable::VariableQuantity q) : 
    MatrixAssembler(), alpha(alpha), deltaT(deltaT)
{
    this->q = q;
}


void ScalarAdvectionRhsAssembler :: matrixFromElement(FloatMatrix &answer, Element &el, TimeStep *tStep) const
{
    FloatMatrix contrib;
    IntArray loc;
    MPElement *e = dynamic_cast<MPElement*>(&el);
    int ndofs = e->giveNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();

    e->getLocalCodeNumbers (loc, q);

    e->giveCharacteristicMatrix(contrib, MassMatrix, tStep);
    if (contrib.isNotEmpty()) {
        contrib.times((-1.0)*this->deltaT/2.0);
        answer.assemble(contrib, loc, loc);
    }
    e->giveCharacteristicMatrix(contrib, StiffnessMatrix, tStep);
    if (contrib.isNotEmpty()) {
        answer.assemble(contrib, loc, loc);
    }
    // boundary terms
    e->giveCharacteristicMatrix(contrib, InternalFluxVector, tStep);
    if (contrib.isNotEmpty()) {
        answer.assemble(contrib, loc, loc);
    }

}


DGProblem :: DGProblem(int i, EngngModel *_master = nullptr) : EngngModel(i, _master)
{
    ndomains = 1;
}

NumericalMethod *DGProblem :: giveNumericalMethod(MetaStep *mStep)
{
    if ( !nMethod ) {
        if ( isParallel() ) {
            if ( ( solverType == ST_Petsc ) || ( solverType == ST_Feti ) ) {
                nMethod = classFactory.createSparseLinSolver(solverType, this->giveDomain(1), this);
            }
        } else {
            nMethod = classFactory.createSparseLinSolver(solverType, this->giveDomain(1), this);
        }
        if ( !nMethod ) {
            OOFEM_ERROR("linear solver creation failed for lstype %d", solverType);
        }
    }


    return nMethod.get();
}

  
void
DGProblem :: initializeFrom(InputRecord &ir)
{
    EngngModel :: initializeFrom(ir);

    int val = SMT_Skyline;
    IR_GIVE_OPTIONAL_FIELD(ir, val, _IFT_EngngModel_smtype);
    this->sparseMtrxType = ( SparseMtrxType ) val;

    if ( ir.hasField(_IFT_DGProblem_initt) ) {
        IR_GIVE_FIELD(ir, initT, _IFT_DGProblem_initt);
    }

    if ( ir.hasField(_IFT_DGProblem_deltat) ) {
        IR_GIVE_FIELD(ir, deltaT, _IFT_DGProblem_deltat);
    } else if ( ir.hasField(_IFT_DGProblem_deltatfunction) ) {
        IR_GIVE_FIELD(ir, dtFunction, _IFT_DGProblem_deltatfunction);
    } else if ( ir.hasField(_IFT_DGProblem_prescribedtimes) ) {
        IR_GIVE_FIELD(ir, prescribedTimes, _IFT_DGProblem_prescribedtimes);
    } else {
        throw ValueInputException(ir, "none", "Time step not defined");
    }

    IR_GIVE_FIELD(ir, alpha, _IFT_DGProblem_alpha);
    problemType = "sa"; // compatibility mode @TODO Remove default value
    IR_GIVE_OPTIONAL_FIELD(ir, problemType, _IFT_DGProblem_problemType);
    if (!((problemType == "sa")||(problemType == "tm"))) {
      throw ValueInputException(ir, "none", "Problem type not recognized");
    }
    OOFEM_LOG_RELEVANT("DG: %s formulation\n", problemType.c_str());
    
    this->keepTangent = ir.hasField(_IFT_DGProblem_keepTangent);
    field = std::make_unique<DofDistributedPrimaryField>(this, 1, FT_TransportProblemUnknowns, 2, this->alpha);

    // read field export flag
    exportFields.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, exportFields, _IFT_DGProblem_exportFields);
    if ( exportFields.giveSize() ) {
        FieldManager *fm = this->giveContext()->giveFieldManager();
        for ( int i = 1; i <= exportFields.giveSize(); i++ ) {
            if ( exportFields.at(i) == FT_Displacements ) {
                FieldPtr _displacementField( new MaskedPrimaryField ( ( FieldType ) exportFields.at(i), this->field.get(), {D_u, D_v, D_w} ) );
                fm->registerField( _displacementField, ( FieldType ) exportFields.at(i) );
            } else if ( exportFields.at(i) == FT_Pressure ) {
                FieldPtr _pressureField( new MaskedPrimaryField ( ( FieldType ) exportFields.at(i), this->field.get(), {P_f} ) );
                fm->registerField( _pressureField, ( FieldType ) exportFields.at(i) );
            }
        }
    }
}

double DGProblem :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    return this->field->giveUnknownValue(dof, mode, tStep);
}


double
DGProblem :: giveDeltaT(int n)
{
    if ( this->dtFunction ) {
        return this->giveDomain(1)->giveFunction(this->dtFunction)->evaluateAtTime(n);
    } else if ( this->prescribedTimes.giveSize() > 0 ) {
        return this->giveDiscreteTime(n) - this->giveDiscreteTime(n - 1);
    } else {
        return this->deltaT;
    }
}

double
DGProblem :: giveDiscreteTime(int iStep)
{
    if ( iStep > 0 && iStep <= this->prescribedTimes.giveSize() ) {
        return ( this->prescribedTimes.at(iStep) );
    } else if ( iStep == 0 ) {
        return initT;
    }

    OOFEM_ERROR("invalid iStep");
    return 0.0;
}


TimeStep *DGProblem :: giveNextStep()
{
    if ( !currentStep ) {
        // first step -> generate initial step
        currentStep = std::make_unique<TimeStep>( *giveSolutionStepWhenIcApply() );
    }

    double dt = this->giveDeltaT(currentStep->giveNumber()+1);
    previousStep = std :: move(currentStep);
    currentStep = std::make_unique<TimeStep>(*previousStep, dt);
    currentStep->setIntrinsicTime(previousStep->giveTargetTime() + alpha * dt);
    return currentStep.get();
}


TimeStep* DGProblem :: giveSolutionStepWhenIcApply(bool force)
{
    if ( master && (!force)) {
        return master->giveSolutionStepWhenIcApply();
    } else {
        if ( !stepWhenIcApply ) {
            double dt = this->giveDeltaT(1);
            stepWhenIcApply = std::make_unique<TimeStep>(giveNumberOfTimeStepWhenIcApply(), this, 0, this->initT, dt, 0);
            // The initial step goes from [-dt, 0], so the intrinsic time is at: -deltaT  + alpha*dt
            stepWhenIcApply->setIntrinsicTime(-dt + alpha * dt);
        }

        return stepWhenIcApply.get();
    }
}

void
DGProblem :: constructBoundaryEntities () {
    // This method construct boundary entities of dimension nsd-1 (edges in 2D, surfaces in 3D) from the given mesh entities of dimension nsd.
    // first loop over elements and their boundary entities
    ConnectivityTable *ctable = this->giveDomain(1)->giveConnectivityTable();
    // sets of processed boundary entities per element
    std::vector<std::set<int>> processedBoundaryEntities;
    processedBoundaryEntities.resize(this->giveDomain(1)->giveNumberOfElements());

    for ( auto &elem: this->giveDomain(1)->giveElements() ) {
        Element *e = elem.get();
        for ( int i = 1; i <= e->giveNumberOfBoundarySides(); i++ ) { // should rename as this is effectively number of boundary entities (edges or surfaces) depending on element dimension
            if ( processedBoundaryEntities[e->giveNumber()-1].find(i) != processedBoundaryEntities[e->giveNumber()-1].end() ) {
                continue;  // edge already processed by neighbor element
            }
            IntArray bnodes, neighbors;
            bnodes = e->giveBoundaryNodes(i);
            for ( int j = 1; j <= bnodes.giveSize(); j++ ) {
                bnodes.at(j) = e->giveDofManager(bnodes.at(j))->giveNumber();
            }
            bnodes.sort();
            // now search element neighbors to find one sharing the same boundary nodes
            ctable->giveElementsWithNodes(neighbors, bnodes);
            if ( neighbors.giveSize() > 2 ) {
                OOFEM_ERROR("More than two neighbors found for boundary entity (elem %d, boundary %d)", e->giveNumber(), i);
            } else {
                std::unique_ptr<DGBoundaryEntity> be = std::make_unique<DGBoundaryEntity>();
                be->addElement(e->giveNumber(), i);
                processedBoundaryEntities[e->giveNumber()-1].insert(i);

                if ( neighbors.giveSize() == 2 ) {
                    //boundary entity shared => we are on domain interior
                    int neighborelem = neighbors.at(1) == e->giveNumber() ? neighbors.at(2) : neighbors.at(1);
                    // determine boundary entity number in neighbor element
                    int neighborboundary = 0;
                    for ( int j = 1; j <= this->giveDomain(1)->giveElement(neighborelem)->giveNumberOfBoundarySides(); j++ ) {
                        bool equal = true;
                        IntArray neighborBoundaryNodes = this->giveDomain(1)->giveElement(neighborelem)->giveBoundaryNodes(j);
                        for ( int k = 1; k <= neighborBoundaryNodes.giveSize(); k++ ) {
                            neighborBoundaryNodes.at(k) = this->giveDomain(1)->giveElement(neighborelem)->giveDofManager(neighborBoundaryNodes.at(k))->giveNumber();    
                        }
                        // compare bnodes (sorted) with neighborBoundaryNodes
                        for ( int k = 1; k <= neighborBoundaryNodes.giveSize(); k++ ) {
                            if ( !bnodes.findSorted(neighborBoundaryNodes.at(k))) {
                                equal = false;
                                break;
                            }
                        }
                        if ( equal) {   
                            neighborboundary = j;
                            break;
                        }
                    }
                    if (neighborboundary) {
                        be->addElement(neighborelem, neighborboundary);
                        processedBoundaryEntities[neighborelem-1].insert(neighborboundary);
                    } else {
                        OOFEM_ERROR("Boundary entity not found in neighbor element");
                    }
                }
                this->boundaryEntities.push_back(std::move(be));    
            }
        }
    }
}





void DGProblem :: solveYourselfAt(TimeStep *tStep)
{
    OOFEM_LOG_INFO( "\nSolving [step number %5d, time %e]\n", tStep->giveNumber(), tStep->giveTargetTime() );
    
    Domain *d = this->giveDomain(1);
    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );

    if ( false ) {
        this->constructBoundaryEntities();
        if (true) {
            // print boundary entities
            int id = 1;
            printf("Boundary entities: %ld\n---------------------------\n", this->boundaryEntities.size());
            for ( auto &be: this->boundaryEntities ) {
                if ( be->elements.giveSize() == 1 ) {
                printf("%4d: %4d(%4d) %4s(%4s)\n", id++, be->elements.at(1), be->elementBoundaryIDs.at(1), "-", "-");
                } else {
                    printf("%4d: %4d(%4d) %4d(%4d)\n", id++, be->elements.at(1), be->elementBoundaryIDs.at(1), be->elements.at(2), be->elementBoundaryIDs.at(2));
                }
            }
            printf("---------------------------\n");
            
            printf("Element coloring\n---------------------------\n");
            for (auto &e: this->giveDomain(1)->giveElements()) {
                int c = this->giveDomain(1)->giveConnectivityTable()->getElementColor(e.get()->giveNumber());
                printf("%4d: %2d\n", e->giveNumber(), c);
            }
            printf("---------------------------\n");
        }
        

    }

    // @BP: debug
    
    if ( tStep->isTheFirstStep() ) {
        this->applyIC();

        field->advanceSolution(tStep);
        field->initialize(VM_Total, tStep, solution, EModelDefaultEquationNumbering());

        if ( !lhsMatrix ) {
            lhsMatrix = classFactory.createSparseMtrx(sparseMtrxType);
            lhsMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
        }
        if ( !rhsMatrix ) {
            rhsMatrix = classFactory.createSparseMtrx(sparseMtrxType);
            rhsMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
        }

        OOFEM_LOG_INFO("Assembling system matrices\n");
        this->assemble( *lhsMatrix, tStep, ScalarAdvectionLhsAssembler(this->alpha, tStep->giveTimeIncrement(), unknownQuantity), EModelDefaultEquationNumbering(), d );
        this->assemble( *rhsMatrix, tStep, ScalarAdvectionRhsAssembler(this->alpha, tStep->giveTimeIncrement(), unknownQuantity), EModelDefaultEquationNumbering(), d );
    }
    OOFEM_LOG_INFO("Assembling boundary contributions\n");
    // loop over boundary entities
    // @BP instead using BoundaryEntity, we can set up MPMElement based boundary element representing the boundary entity, 
    // this would allow to use the same assembly code as for the interior
    // the term to evaluate over the boundary is \int (f(s+, s-)\cdot (nv)) ds, where f is numerical flux. So this term has to be parametrized by numerical flux...

    // @BP start of editing
    this->giveNumericalMethod( this->giveCurrentMetaStep() );

    FloatArray *pv = this->field->giveSolutionVector(tStep->givePreviousStep());
    FloatArray *v = this->field->giveSolutionVector(tStep);
    FloatArray rhs;
    rhsMatrix->times(*pv, rhs);
    ConvergedReason status = this->nMethod->solve(*lhsMatrix, rhs, *v);
    tStep->convergedReason = status;
    this->updateSolution(*v, tStep, d); // ?
    // @BP end of editing
}


void
DGProblem :: updateSolution(FloatArray &solutionVector, TimeStep *tStep, Domain *d)
{
    ///@todo NRSolver should report when the solution changes instead of doing it this way.
    this->field->update(VM_Total, tStep, solutionVector, EModelDefaultEquationNumbering());
    ///@todo Need to reset the boundary conditions properly since some "update" is doing strange
    /// things such as applying the (wrong) boundary conditions. This call will be removed when that code can be removed.
    this->field->applyBoundaryCondition(tStep);
}


void
DGProblem :: updateInternalRHS(FloatArray &answer, TimeStep *tStep, Domain *d, FloatArray *eNorm)
{}


void
DGProblem :: updateMatrix(SparseMtrx &mat, TimeStep *tStep, Domain *d)
{}


void
DGProblem :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
{
    OOFEM_ERROR("Unknown component");
}


void
DGProblem :: applyIC()
{
    Domain *domain = this->giveDomain(1);
    OOFEM_LOG_INFO("Applying initial conditions\n");

    this->field->applyDefaultInitialCondition();

    // set initial field IP values (needed by some nonlinear materials)
    TimeStep *s = this->giveSolutionStepWhenIcApply();
    for ( auto &elem : domain->giveElements() ) {
        Element *element = elem.get();
        element->updateInternalState(s);
        element->updateYourself(s);
    }
}


bool
DGProblem :: requiresEquationRenumbering(TimeStep *tStep)
{
    ///@todo This method should be set as the default behavior instead of relying on a user specified flag. Then this function should be removed.
    if ( tStep->isTheFirstStep() ) {
        return true;
    }
    // Check if Dirichlet b.c.s has changed.
    Domain *d = this->giveDomain(1);
    for ( auto &gbc : d->giveBcs() ) {
        ActiveBoundaryCondition *active_bc = dynamic_cast< ActiveBoundaryCondition * >(gbc.get());
        BoundaryCondition *bc = dynamic_cast< BoundaryCondition * >(gbc.get());
        // We only need to consider Dirichlet b.c.s
        if ( bc || ( active_bc && ( active_bc->requiresActiveDofs() || active_bc->giveNumberOfInternalDofManagers() ) ) ) {
            // Check of the dirichlet b.c. has changed in the last step (if so we need to renumber)
            if ( gbc->isImposed(tStep) != gbc->isImposed(tStep->givePreviousStep()) ) {
                return true;
            }
        }
    }
    return false;
}

int
DGProblem :: forceEquationNumbering()
{
    this->lhsMatrix = nullptr;
    this->rhsMatrix = nullptr;
    return EngngModel :: forceEquationNumbering();
}


void
DGProblem :: printOutputAt(FILE *file, TimeStep *tStep)
{
    if ( !this->giveDomain(1)->giveOutputManager()->testTimeStepOutput(tStep) ) {
        return; // Do not print even Solution step header
    }

    EngngModel :: printOutputAt(file, tStep);
}


void
DGProblem :: updateYourself(TimeStep *tStep)
{
    EngngModel :: updateYourself(tStep);
}

void
DGProblem :: saveContext(DataStream &stream, ContextMode mode)
{
    EngngModel :: saveContext(stream, mode);
    field->saveContext(stream);
}


void
DGProblem :: restoreContext(DataStream &stream, ContextMode mode)
{
    EngngModel :: restoreContext(stream, mode);
    field->restoreContext(stream);
}


int
DGProblem :: giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep)
{
    return tStep->giveNumber() % 2;
}


int
DGProblem :: requiresUnknownsDictionaryUpdate()
{
    return true;
}

int
DGProblem :: checkConsistency()
{
  return EngngModel :: checkConsistency();
}


void
DGProblem :: updateDomainLinks()
{
    EngngModel :: updateDomainLinks();
    this->giveNumericalMethod( this->giveCurrentMetaStep() )->setDomain( this->giveDomain(1) );
}


FieldPtr DGProblem::giveField(FieldType key, TimeStep *tStep)
{
    /* Note: the current implementation uses MaskedPrimaryField, that is automatically updated with the model progress, 
        so the returned field always refers to active solution step. 
    */

    if ( tStep != this->giveCurrentStep()) {
        OOFEM_ERROR("Unable to return field representation for non-current time step");
    }
    if ( key == FT_Displacements ) {
      return std::make_shared<MaskedPrimaryField>( key, this->field.get(), IntArray{D_u, D_v, D_w} );
    } else if ( key == FT_Pressure ) {
        return std::make_shared<MaskedPrimaryField>( key, this->field.get(), IntArray{P_f} );
    } else if ( key == FT_Temperature ) {
        return std::make_shared<MaskedPrimaryField>( key, this->field.get(), IntArray{T_f} );
    } else {
        return FieldPtr();
    }
}



} // end namespace oofem
