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

#include "../sm/EngineeringModels/staticstructuralstaggered.h"
#include "classfactory.h"
#include "primaryfield.h"

namespace oofem {
REGISTER_EngngModel(StaticStructuralStaggered);

StaticStructuralStaggered :: StaticStructuralStaggered(int i, EngngModel *_master) :  StaticStructural(i, _master)
{
    ndomains = 1;

}


StaticStructuralStaggered :: ~StaticStructuralStaggered()
{
    delete field;
    delete stiffnessMatrix;
    delete nMethod;
}


NumericalMethod *StaticStructuralStaggered :: giveNumericalMethod(MetaStep *mStep)
{
    if ( nMethod ) {
        return nMethod;
    }
    nMethod = new NRSolver(this->giveDomain(1), this);
    //nMethod = new StaggeredSolver(this->giveDomain(1), this);
    return nMethod;
}

IRResultType
StaticStructuralStaggered :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    StaticStructural :: initializeFrom(ir);
    
    // set the dofIdArrays for each solution field to solve for
    IntArray idList1 = {1, 2};
    IntArray idList2 = {3};
    this->UnknownNumberingSchemeList.resize(2);
    this->UnknownNumberingSchemeList[0].setDofIdArray(idList1);
    this->UnknownNumberingSchemeList[0].setNumber(1);
    this->UnknownNumberingSchemeList[1].setDofIdArray(idList2);
    this->UnknownNumberingSchemeList[1].setNumber(2);    
    
    return IRRT_OK;
}






void StaticStructuralStaggered :: solveYourselfAt(TimeStep *tStep)
{
    int di = 1;
    // total number of dofs in the problem
    //int neq = this->giveNumberOfDomainEquations( di, EModelDefaultEquationNumbering() );
    int neq = StaticStructural :: giveNumberOfDomainEquations( di, EModelDefaultEquationNumbering() );
    
    
    field->advanceSolution(tStep);

    this->giveNumericalMethod( this->giveCurrentMetaStep() );
    this->initMetaStepAttributes( this->giveCurrentMetaStep() );

    // Fetch vector to fill in from primary field.
    this->solution = field->giveSolutionVector(tStep);


    if(!tStep->isTheFirstStep()) {
        // Old solution as starting guess
        FloatArray *oldSol = field->giveSolutionVector(tStep->givePreviousStep());
        *solution = *oldSol;
    }

    if(solution->giveSize() != neq) {
        printf("Resizing.\n");
        this->solution->resize(neq);
        this->solution->zero();
    }


    // Create "stiffness matrix"
    if ( !this->stiffnessMatrix ) {
        this->stiffnessMatrix = classFactory.createSparseMtrx(sparseMtrxType);
        if ( !this->stiffnessMatrix ) {
            OOFEM_ERROR("Couldn't create requested sparse matrix of type %d", sparseMtrxType);
        }

        this->stiffnessMatrix->buildInternalStructure( this, di, EModelDefaultEquationNumbering() );
    }
    this->internalForces.resize(neq);

    FloatArray incrementOfSolution(neq);

    // Build initial/external load
    FloatArray externalForces(neq);
    externalForces.zero();
    

    // Create stiffness matrices and internal forces vectors corresponding to each dof group
    int numDofIdGroups = this->UnknownNumberingSchemeList.size();
    if ( !this->stiffnessMatrixList.size() ) {
        this->fIntList.resize(numDofIdGroups);
        this->fExtList.resize(numDofIdGroups);            
        this->stiffnessMatrixList.resize(numDofIdGroups);

        for ( int dG = 0; dG < numDofIdGroups; dG++ ) {
	    int neq = giveNumberOfEquations(di, this->UnknownNumberingSchemeList[dG]);
            this->fIntList[dG].resize(neq);    
            this->fExtList[dG].resize(neq);    
            
            SparseMtrx *mtrx = this->stiffnessMatrixList[dG];
            mtrx = classFactory.createSparseMtrx(sparseMtrxType); 
            
            if ( !mtrx ) {
            OOFEM_ERROR("Couldn't create requested sparse matrix of type %d", sparseMtrxType );
            }
            
            mtrx->buildInternalStructure( this, di, UnknownNumberingSchemeList[dG] );
            //mtrx->buildInternalStructure( engngModel, di, EModelDefaultEquationNumbering() );
            //mtrx->printStatistics();
            //mtrx->printYourself();
            
        }
    }    
    
    
    this->assembleVector( externalForces, tStep, ExternalForcesVector, VM_Total,
                         EModelDefaultEquationNumbering(), this->giveDomain(di) );
#ifdef __PARALLEL_MODE
    this->updateSharedDofManagers(externalForces, EModelDefaultEquationNumbering(), LoadExchangeTag);
#endif

    externalForces.printYourself("F_ext_total");
    // Compute external forces 
    for ( int dG = 0; dG < numDofIdGroups; dG++ ) {
        this->updateExternalForcesForStaggeredSolver(this->fExtList[dG], tStep, this->giveDomain(di), this->UnknownNumberingSchemeList[dG]);
	//this->updateInternalForcesForStaggeredSolver(this->fIntList[dG], tStep, this->giveDomain(di), this->UnknownNumberingSchemeList[dG]);
        //RRT(dG) = this->fExtList[dG].computeSquaredNorm();
        
        this->fExtList[dG].printYourself("f_ext");
	this->fIntList[dG].printYourself("f_int");
    }    
    
    
    if ( this->giveProblemScale() == macroScale ) {
        OOFEM_LOG_INFO("\nStaticStructuralStaggered :: solveYourselfAt - Solving step %d, metastep %d, (neq = %d)\n", tStep->giveNumber(), tStep->giveMetaStepNumber(), neq);
    }

    double loadLevel;
    int currentIterations;
    NM_Status status = this->nMethod->solve(this->stiffnessMatrix,
                                            &this->fExtList[0],
                                            NULL,
                                            solution,
                                            & incrementOfSolution,
                                            & ( this->fIntList[0] ),
                                            this->eNorm,
                                            loadLevel, // Only relevant for incrementalBCLoadVector?
                                            SparseNonLinearSystemNM :: rlm_total,
                                            currentIterations,
                                            tStep);    
    
//     NM_Status status = this->nMethod->solve(this->stiffnessMatrix,
//                                             & externalForces,
//                                             NULL,
//                                             solution,
//                                             & incrementOfSolution,
//                                             & ( this->internalForces ),
//                                             this->eNorm,
//                                             loadLevel, // Only relevant for incrementalBCLoadVector?
//                                             SparseNonLinearSystemNM :: rlm_total,
//                                             currentIterations,
//                                             tStep);

    if ( !( status & NM_Success ) ) {
        OOFEM_ERROR("No success in solving problem");
    }
}




int
StaticStructuralStaggered :: forceEquationNumbering()
{
    int numEqn = StructuralEngngModel::forceEquationNumbering();

    delete stiffnessMatrix;
    stiffnessMatrix = NULL;

    return numEqn;
}

int
StaticStructuralStaggered :: giveNewEquationNumber(int domain, DofIDItem id)
{
    int numDofIdGroups = this->UnknownNumberingSchemeList.size();
    for ( int dG = 0; dG < numDofIdGroups; dG++ ) {
	CustomEquationNumbering &uS = this->UnknownNumberingSchemeList[dG];
	if ( uS.dofIdArray.contains( (int)id ) ) {
	    int eq = uS.giveNewEquationNumber();
	    //printf("added eq %d \n", eq);
	    //return eq;
	}
    } 
    //return 0;
    return StaticStructural :: giveNewEquationNumber(domain, id);
}


int
StaticStructuralStaggered :: giveNewPrescribedEquationNumber(int domain, DofIDItem id)
{
    int numDofIdGroups = this->UnknownNumberingSchemeList.size();
    for ( int dG = 0; dG < numDofIdGroups; dG++ ) {
	CustomEquationNumbering &uS = this->UnknownNumberingSchemeList[dG];
	if ( uS.dofIdArray.contains( (int)id ) ) {
	    int eq = uS.giveNewPrescribedEquationNumber();
	    //printf("added pres eq %d \n", eq);
	    //return eq;
	}
    } 
    //return 0;
    return StaticStructural :: giveNewPrescribedEquationNumber(domain, id);
}



int
StaticStructuralStaggered :: giveNumberOfEquations(int id, CustomEquationNumbering &num)
{
    //
    // returns number of equations of current problem
    // this method is implemented here, because some method may add some
    // conditions into the system and this may results into increased number of
    // equations.
    //
    if ( !equationNumberingCompleted ) {
        this->forceEquationNumbering();
    }

    return num.isDefault() ? num.giveNumEquations() : num.giveNumPresEquations();
}








void StaticStructuralStaggered :: terminate(TimeStep *tStep)
{
    StaticStructural :: terminate(tStep);
}

double StaticStructuralStaggered :: giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof)
{
    return this->field->giveUnknownValue(dof, mode, tStep);
}




void StaticStructuralStaggered :: updateInternalForcesForStaggeredSolver(FloatArray &answer, TimeStep* tStep, Domain* d, const CustomEquationNumbering& s)
{
    answer.zero();
    this->assembleVector(answer, tStep, InternalForcesVector, VM_Total, s, d, & this->eNorm);
#ifdef __PARALLEL_MODE
    this->updateSharedDofManagers(answer, s, InternalForcesExchangeTag);
#endif
    internalVarUpdateStamp = tStep->giveSolutionStateCounter(); // Hack for linearstatic

}


void StaticStructuralStaggered :: updateExternalForcesForStaggeredSolver(FloatArray &answer, TimeStep* tStep, Domain* d, const CustomEquationNumbering& s)
{
    answer.zero();
    this->assembleVector( answer, tStep, ExternalForcesVector, VM_Total, s, d );
}


void StaticStructuralStaggered :: updateTangentStiffnessForStaggeredSolver(SparseMtrx *answer, TimeStep* tStep, Domain* d, const CustomEquationNumbering& s)
{
    answer->zero();
    this->assemble(answer, tStep, TangentStiffnessMatrix, s, d);

} 




} // end namespace oofem
