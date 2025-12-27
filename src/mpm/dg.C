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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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
#include "crosssection.h"
#include "material.h"
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
        answer.assemble(contrib, loc, loc);
    }
    e->giveCharacteristicMatrix(contrib, StiffnessMatrix, tStep);
    if (contrib.isNotEmpty()) {
        contrib.times(this->deltaT/2.0);
        answer.assemble(contrib, loc, loc);
    }
    // boundary terms
    e->giveCharacteristicMatrix(contrib, InternalFluxVector, tStep);
    if (contrib.isNotEmpty()) {
        contrib.times(this->deltaT/2.0);
        answer.assemble(contrib, loc, loc);
    }
}

ScalarAdvectionRhsAssembler :: ScalarAdvectionRhsAssembler(double alpha, double deltaT, Variable::VariableQuantity q) : 
    VectorAssembler(), alpha(alpha), deltaT(deltaT)
{
    this->q = q;
}


void ScalarAdvectionRhsAssembler :: vectorFromElement(FloatArray &answer, Element &el, TimeStep *tStep, ValueModeType mode) const
{
    FloatMatrix contrib, rhsMatrix;
    IntArray loc, dofids;
    MPElement *e = dynamic_cast<MPElement*>(&el);
    int ndofs = e->giveNumberOfDofs();
    answer.resize(ndofs);
    answer.zero();

    rhsMatrix.resize(ndofs, ndofs);
    rhsMatrix.zero();

    e->getLocalCodeNumbers (loc, q);

    e->giveCharacteristicMatrix(contrib, MassMatrix, tStep);
    if (contrib.isNotEmpty()) {
        rhsMatrix.add(contrib);
    }
    e->giveCharacteristicMatrix(contrib, StiffnessMatrix, tStep);
    if (contrib.isNotEmpty()) {
        contrib.times((-1.0)*this->deltaT/2.0);
        rhsMatrix.add(contrib);
    }
    // boundary terms
    e->giveCharacteristicMatrix(contrib, InternalFluxVector, tStep);
    if (contrib.isNotEmpty()) {
        contrib.times((-1.0)*this->deltaT/2.0);
        rhsMatrix.add(contrib);
    }
    FloatArray help, rp;
    e->computeVectorOf(VM_Total, tStep->givePreviousStep(), rp);  

    help.beProductOf(rhsMatrix, rp);
    answer.assemble(help, loc);
}


void 
ClonedDofManager::printOutputAt(FILE *file, TimeStep *tStep) {
   EngngModel *emodel = this->giveDomain()->giveEngngModel();

    fprintf( file, "%-8s%8d (%8d), Master:%d:\n", this->giveClassName(), this->giveLabel(), this->giveNumber(), master );
    for ( Dof *dof: *this ) {
        emodel->printDofOutputAt(file, dof, tStep);
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
    if ( ir.hasField(_IFT_DGProblem_preprocessFEM2DG) ) {
        preprocessFEM2DG = true;
    
        IR_GIVE_OPTIONAL_FIELD(ir, sets2preprocess, _IFT_DGProblem_sets2preprocess);
        IR_GIVE_OPTIONAL_FIELD(ir, targetBoundaryNodeSets, _IFT_DGProblem_targetBoundaryNodeSets);
        if ( sets2preprocess.giveSize() != targetBoundaryNodeSets.giveSize() ) {
            OOFEM_ERROR("Size mismatch in %s and %s attributes", _IFT_DGProblem_sets2preprocess, _IFT_DGProblem_targetBoundaryNodeSets);
        }
        IR_GIVE_FIELD(ir, targetAllNodeSet, _IFT_DGProblem_targetAllNodeSet);
        IR_GIVE_FIELD(ir, targetAllElementSet, _IFT_DGProblem_targetAllElementSet);
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
    Domain *domain = this->giveDomain(1);
    // set maps to assist in generating boundary node sets
    std::vector<std::unordered_multimap<int, int>> setsBoundaryEntities;
    // boundary nodes sets to be generated
    std::vector<std::set<int>> boundaryNodeSets;

    int nnodes = domain->giveNumberOfDofManagers();
    int nelems = domain->giveNumberOfElements();
    int nodeNum = nnodes;
    int elemNum = nelems;

    Timer timer;
    timer.startTimer();

    // initialize setBoundaryEntities
    setsBoundaryEntities.resize(this->sets2preprocess.giveSize());
    boundaryNodeSets.resize(this->sets2preprocess.giveSize());
    for (int iset = 1; iset<=sets2preprocess.giveSize(); iset++) {
        const IntArray& bentities = domain->giveSet(sets2preprocess.at(iset))->giveBoundaryList();
        for (int i = 1; i <= bentities.giveSize()/2; i++) {
            setsBoundaryEntities[iset-1].insert({bentities.at(2*i-1), bentities.at(2*i)});
        }
    }

    std::vector<IntArray> clonedElementNodes(nelems);
    //OOFEM_LOG_INFO ("fem2DG: Decoupling elements\n");
    for ( int ielem = 1; ielem<=nelems; ielem++) {
        // decouple original (FE) elements
        Element *e = domain->giveElement(ielem);
        int numNodes = e->giveNumberOfDofManagers();
        IntArray clonedNodes(numNodes);
        for (int i = 1; i <= numNodes; i++) {
            auto cnode = std::make_unique<ClonedDofManager>(++nodeNum, domain, e->giveDofManagerNumber(i));
            cnode->setCoordinates( domain->giveDofManager(e->giveDofManagerNumber(i))->giveCoordinates() );
            domain->resizeDofManagers(nodeNum);
            domain->setDofManager(nodeNum, std::move(cnode));
            clonedNodes.at(i) = nodeNum;
        }
        // store new element connectivity
        clonedElementNodes.at(ielem-1)=clonedNodes;// e->setDofManagers(clonedNodes);
    }
    //OOFEM_LOG_INFO ("fem2DG: Constructing boundary entities\n");
    for ( int ielem = 1; ielem<=nelems; ielem++) {
        Element *e = domain->giveElement(ielem);
        for ( int i = 1; i <= e->giveNumberOfBoundarySides(); i++ ) { // should rename as this is effectively number of boundary entities (edges or surfaces) depending on element dimension
            if ( processedBoundaryEntities[e->giveNumber()-1].find(i) != processedBoundaryEntities[e->giveNumber()-1].end() ) {
                continue;  // edge already processed by neighbor element
            }
            IntArray bnodes, bnodesSorted, bclonednodes, neighbors;
            bnodes = e->giveBoundaryNodes(i);
            bclonednodes.resize(bnodes.giveSize());
            for ( int j = 1; j <= bnodes.giveSize(); j++ ) {
                bclonednodes.at(j) = clonedElementNodes.at(ielem-1).at(bnodes.at(j));
                bnodes.at(j) = e->giveDofManager(bnodes.at(j))->giveNumber();
            }
            bnodesSorted = bnodes;
            bnodesSorted.sort();
            // now search element neighbors to find one sharing the same boundary nodes
            ctable->giveElementsWithNodes(neighbors, bnodes);
            if ( neighbors.giveSize() > 2 ) {
                OOFEM_ERROR("More than two neighbors found for boundary entity (elem %d, boundary %d)", e->giveNumber(), i);
            } else {
                std::unique_ptr<DGBoundaryEntity> be = std::make_unique<DGBoundaryEntity>();
                be->addElement(e->giveNumber(), i);
                processedBoundaryEntities[e->giveNumber()-1].insert(i);
                IntArray bentityNodes(bclonednodes); // nodes of boundary entity, initially bnodes followed by cloned or neighbor nodes
                int neighborboundary=0; // identified neighbor boundary id

                if ( neighbors.giveSize() == 1 ) {
                    //boundary entity not shared => we are on domain boundary
                    // create cloned nodes on boundary
                    
                    for (int j = 1; j <= bnodes.giveSize(); j++)
                    {
                        auto cnode = std::make_unique<ClonedDofManager>(++nodeNum, domain, bnodes.at(j));
                        cnode->setCoordinates( domain->giveDofManager(bnodes.at(j))->giveCoordinates() );
                        domain->resizeDofManagers(nodeNum);
                        domain->setDofManager(nodeNum, std::move(cnode));
                        bentityNodes.followedBy(nodeNum, bnodes.giveSize());
                    }

                    // update boundary sets
                    for (int iset = 1; iset<=sets2preprocess.giveSize(); iset++) {
                        auto range = setsBoundaryEntities[iset-1].equal_range(e->giveNumber());
                        for (auto it = range.first; it != range.second; ++it) {
                            if (it->second == i) {
                                // add boundary nodes to target set
                                for (int k = 1; k <= bnodes.giveSize(); k++) {
                                    boundaryNodeSets.at(iset-1).insert(bentityNodes.at(bnodes.giveSize()+k));
                                }
                                break;  // only one boundary entity per element
                            }
                        }
                    }

                } else if ( neighbors.giveSize() == 2 ) {
                    //boundary entity shared => we are on domain interior
                    int neighborelem = neighbors.at(1) == e->giveNumber() ? neighbors.at(2) : neighbors.at(1);
                    // neighbor boundary nodes
                    IntArray neighborBoundaryNodes;
                    IntArray neighborClonedBoundaryNodes;
                    for ( int j = 1; j <= this->giveDomain(1)->giveElement(neighborelem)->giveNumberOfBoundarySides(); j++ ) {
                        bool equal = true;
                        neighborBoundaryNodes = this->giveDomain(1)->giveElement(neighborelem)->giveBoundaryNodes(j);
                        int size = neighborBoundaryNodes.giveSize();
                        neighborClonedBoundaryNodes.resize(size);
                        for ( int k = 1; k <= size; k++ ) {
                            neighborClonedBoundaryNodes.at(k) = clonedElementNodes.at(neighborelem-1).at(neighborBoundaryNodes.at(k));
                            neighborBoundaryNodes.at(k) = this->giveDomain(1)->giveElement(neighborelem)->giveDofManager(neighborBoundaryNodes.at(k))->giveNumber();    
                        }
                        // compare bnodes (sorted) with neighborBoundaryNodes
                        for ( int k = 1; k <= neighborBoundaryNodes.giveSize(); k++ ) {
                            if ( !bnodesSorted.findSorted(neighborBoundaryNodes.at(k))) {
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
                        be->addElement(neighborelem, neighborboundary); // delete
                        processedBoundaryEntities[neighborelem-1].insert(neighborboundary); 
                        // add nodes of neighbor element
                        for (int j = 1; j <= bnodes.giveSize(); j++)
                        {
                            bentityNodes.followedBy(neighborClonedBoundaryNodes.at(neighborBoundaryNodes.findFirstIndexOf(bnodes.at(j))), bnodes.giveSize());
                        }
                    } else {
                        OOFEM_ERROR("Boundary entity not found in neighbor element");
                    }
                }
                this->boundaryEntities.push_back(std::move(be));   
                // create boundary element
                    Element_Geometry_Type egt = e->giveInterpolation()->giveBoundaryGeometryType(neighborboundary);
                    std::unique_ptr<Element> belem (this->CreateBoundaryElement(egt, ++elemNum, domain, bentityNodes));
                    domain->resizeElements(elemNum);
                    // @BP TODO fragile, replace by user defined values 
                    belem->setCrossSection(1);
                    belem->setMaterial(1);
                    domain->setElement(elemNum, std::move(belem));
 
            }
        }
    }
    // finally make original elements decoupled
    for ( int ielem = 1; ielem<=nelems; ielem++) {
        Element *e = domain->giveElement(ielem);
        e->setDofManagers(clonedElementNodes.at(ielem-1));
    }
    // update target boundary sets
    for (int iset = 1; iset<=sets2preprocess.giveSize(); iset++) {
        IntArray nodes(boundaryNodeSets.at(iset-1).size());
        int counter = 1;
        for (auto bnode: boundaryNodeSets.at(iset-1)) {
            nodes.at(counter++)=bnode;
        }
        domain->giveSet(targetBoundaryNodeSets.at(iset))->clear();
        domain->giveSet(targetBoundaryNodeSets.at(iset))->setNodeList(nodes);
    }
    // define target all node set
    if ( targetAllNodeSet ) {
        IntArray nodes(nodeNum);
        for (int i = 1; i <= nodeNum; i++) {
            nodes.at(i) = i;
        }
        domain->giveSet(targetAllNodeSet)->clear();
        domain->giveSet(targetAllNodeSet)->setNodeList(nodes);
    }
    // define target all element set
    if ( targetAllElementSet ) {
        IntArray elements(elemNum);
        for (int i = 1; i <= elemNum; i++) {
            elements.at(i) = i;
        }
        domain->giveSet(targetAllElementSet)->clear();
        domain->giveSet(targetAllElementSet)->setElementList(elements);
    } 
    timer.stopTimer();
    //OOFEM_LOG_INFO("fem2DG: generated %d interface elements, %d nodes cloned\n", elemNum-nelems, nodeNum-nnodes); 
    OOFEM_LOG_INFO("fem2DG: FE (%d nodes, %d elements) -> DG (%d nodes, %d elements)\n", nnodes, nelems, nodeNum, elemNum);
    OOFEM_LOG_INFO("fem2DG: done in %.2fs\n", timer.getUtime());
}


std::unique_ptr< Element> 
DGProblem::CreateBoundaryElement( Element_Geometry_Type egt, int elemNum, Domain *domain, IntArray &bentityNodes) const 
{
    std::unique_ptr<Element> belem;
    switch (egt) {
        case EGT_line_1:
            belem = classFactory.createElement("sadgbline1", elemNum, domain);
            break;
        case EGT_quad_1:
            belem = classFactory.createElement("sadgbquad1", elemNum, domain);
            break;
        default:
            OOFEM_ERROR("Unknown boundary element geometry type");
    }
    belem->setDofManagers(bentityNodes);
    return belem;
}

void
DGProblem :: postInitialize()
{
    // set meta step bounds
    int istep = this->giveNumberOfFirstStep(true);
    for ( auto &metaStep: metaStepList ) {
        istep = metaStep.setStepBounds(istep);
    }

    // initialize boundary entities
    if ( preprocessFEM2DG ) {
        this->constructBoundaryEntities();
    }

    for ( auto &domain: domainList ) {
        domain->postInitialize();
    }

#if 0
    if ( preprocessFEM2DG ) {
        // print domain entities
        for ( auto &dman: this->giveDomain(1)->dofManagerList ) {
            ClonedDofManager *cdman = dynamic_cast<ClonedDofManager*>(dman.get());
            Node *ndman = dynamic_cast<Node*>(dman.get());
            std::string name;
            int master=0;
            if (ndman) {
                name = "Node";
            } else if (cdman) {
                name = "ClonedNode";
                master = cdman->giveMasterNumber();
            }
            printf("%10s %4d master %4d coords %6f %6f %6f dofs", name.c_str(), dman->giveNumber(), master, dman->giveCoordinate(1), dman->giveCoordinate(2), dman->giveCoordinate(3));
            for ( Dof *dof: *dman ) {
                printf(" %d", dof->giveDofID());
            }
            printf("\n");
        }
        for (auto &elem: this->giveDomain(1)->elementList) {
            Element *e = elem.get();
            printf("%8s %4d nodes", e->giveClassName(), e->giveNumber());
            for ( int i = 1; i <= e->giveNumberOfDofManagers(); i++ ) {
                printf(" %d", e->giveDofManager(i)->giveNumber());
            }
            printf("\n");
        }

        for (int iset =1; iset<= this->targetBoundaryNodeSets.giveSize(); iset++) {
            Set *set = this->giveDomain(1)->giveSet(this->targetBoundaryNodeSets.at(iset));
            printf("Set %d nodes", set->giveNumber());
            for (int i = 1; i <= (set->giveNodeList()).giveSize(); i++) {
                printf(" %d", (set->giveNodeList()).at(i));
            }
            printf("\n");
        }
    }
#endif
}



void DGProblem :: solveYourselfAt(TimeStep *tStep)
{
    OOFEM_LOG_INFO( "\nSolving [step number %5d, time %e]\n", tStep->giveNumber(), tStep->giveTargetTime() );
    
    Domain *d = this->giveDomain(1);
    int neq = this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering() );
    
    if ( tStep->isTheFirstStep() ) {
        this->applyIC();
        
        if ( !lhsMatrix ) {
            lhsMatrix = classFactory.createSparseMtrx(sparseMtrxType);
            lhsMatrix->buildInternalStructure( this, 1, EModelDefaultEquationNumbering() );
        }
    }
    OOFEM_LOG_INFO("Assembling system matrices\n");
    lhsMatrix->zero();
    this->assemble( *lhsMatrix, tStep, ScalarAdvectionLhsAssembler(this->alpha, tStep->giveTimeIncrement(), unknownQuantity), EModelDefaultEquationNumbering(), d );

    // loop over boundary entities
    // @BP instead using BoundaryEntity, we can set up MPMElement based boundary element representing the boundary entity, 
    // this would allow to use the same assembly code as for the interior
    // the term to evaluate over the boundary is \int (f(s+, s-)\cdot (nv)) ds, where f is numerical flux. So this term has to be parametrized by numerical flux...

    // @BP start of editing
    this->giveNumericalMethod( this->giveCurrentMetaStep() );

    field->advanceSolution(tStep);
    field->initialize(VM_Total, tStep, solution, EModelDefaultEquationNumbering());

    //FloatArray *pv = this->field->giveSolutionVector(tStep->givePreviousStep());
    //FloatArray *v = this->field->giveSolutionVector(tStep);
    FloatArray rhs(neq);  
    this->assembleVector( rhs, tStep, ScalarAdvectionRhsAssembler(this->alpha, tStep->giveTimeIncrement(), unknownQuantity), VM_Total, EModelDefaultEquationNumbering(), d );
    // add dirichlet boundary conditions
    this->assembleDirichletBcRhsVector(rhs, tStep, VM_Total, EModelDefaultEquationNumbering(), d);
    ConvergedReason status = this->nMethod->solve(*lhsMatrix, rhs, solution);
    tStep->convergedReason = status;
    this->updateSolution(solution, tStep, d); // ?
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
    if ( key == FT_Velocity ) {
      return this->giveContext()->giveFieldManager()->giveField(FT_Velocity);
    } else if ( key == FT_Pressure ) {
        return std::make_shared<MaskedPrimaryField>( key, this->field.get(), IntArray{P_f} );
    } else if ( key == FT_Temperature ) {
        return std::make_shared<MaskedPrimaryField>( key, this->field.get(), IntArray{T_f} );
    } else {
        return FieldPtr();
    }
}


/**
 * @BP: This method is used to assemble the Dirichlet boundary conditions into the RHS vector.
 * Note that it assembles the contributions from Lhs and Rhs! 
 * This is needed, as the Rhs is at present evaluated as Matrix * PreviousSolutionVector, and the Dirichlet BCs are not included in previous solution vector!
 */
void
DGProblem::assembleDirichletBcRhsVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode,
                                        const UnknownNumberingScheme &ns, Domain *d) const
{
    IntArray loc, dofids;
    FloatArray rp1, rp, charVec, c1;
    FloatMatrix ke,me,kte,m1;
    FloatMatrix capacity;

    int nelem = d->giveNumberOfElements();

    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
        Element *element = d->giveElement(ielem);

        element->giveElementDofIDMask(dofids);
        element->computeVectorOfPrescribed(dofids, VM_Total, tStep, rp1);
        element->computeVectorOfPrescribed(dofids, VM_Total, tStep->givePreviousStep(), rp);

        if ( !rp.containsOnlyZeroes() || !rp1.containsOnlyZeroes() ) {
            element->giveCharacteristicMatrix(ke, StiffnessMatrix, tStep);
            element->giveCharacteristicMatrix(me, MassMatrix, tStep);
            // internal flux (boundary elements)
            element->giveCharacteristicMatrix(kte, InternalFluxVector, tStep);
            if (ke.isNotEmpty() && kte.isNotEmpty()) {
                ke.add(kte);
            } else if (kte.isNotEmpty()) {
                ke = kte;
            } 
            element->giveLocationArray(loc, ns);

            if ( !rp1.containsOnlyZeroes() ) {
                m1 = ke;
                m1.times(this->deltaT/2.0);
                m1.add(me);
                c1.beProductOf(m1, rp1);
                c1.negated();
                answer.assemble(c1, loc);
            }
/*
            if ( !rp.containsOnlyZeroes() ) {
                m1=ke;
                m1.times((-1.0)*this->deltaT/2.0);
                m1.add(me);
                c1.beProductOf(m1, rp);
                answer.assemble(c1, loc);
            }
*/
        }
    } // end element loop

}


} // end namespace oofem
