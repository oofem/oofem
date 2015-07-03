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

#include "../sm/Elements/3D/macrolspace.h"
#include "../sm/Materials/micromaterial.h"
#include "../sm/EngineeringModels/structengngmodel.h"
#include "../sm/Elements/3D/lspace.h"
#include "fei3dhexalin.h"
#include "constantfunction.h"
#include "domain.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "node.h"

namespace oofem {
REGISTER_Element(MacroLSpace);

//derived from linear brick element
MacroLSpace :: MacroLSpace(int n, Domain *aDomain) : LSpace(n, aDomain)
{
    this->microMasterNodes.clear();
    this->microBoundaryNodes.clear();
    this->firstCall = true;
    microMaterial = NULL;
    microDomain = NULL;
    microEngngModel = NULL;
    this->iteration = 1;
    this->lastStiffMatrixTimeStep = NULL;
}

MacroLSpace :: ~MacroLSpace() { }


IRResultType MacroLSpace :: initializeFrom(InputRecord *ir)
{
    IRResultType result;              // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->microMasterNodes, _IFT_MacroLspace_microMasterNodes);

    if ( this->microMasterNodes.giveSize() != 8 ) {
        OOFEM_WARNING("Need 8 master nodes from the microproblem defined on macroLspace element");
        return IRRT_BAD_FORMAT;
    }

    IR_GIVE_FIELD(ir, this->microBoundaryNodes, _IFT_MacroLspace_microBoundaryNodes);

    microBoundaryDofManager.resize( 3 * microBoundaryNodes.giveSize() );

#if 0
    if ( ir->hasField(_IFT_MacroLspace_stiffMatrxFileName) ) {
        IR_GIVE_OPTIONAL_FIELD(ir, this->stiffMatrxFileName, _IFT_MacroLspace_stiffMatrxFileName);
        if ( fopen(this->stiffMatrxFileName, "r") != NULL ) { //if the file exist
            stiffMatrxFile = fopen(this->stiffMatrxFileName, "r");
            this->stiffMatrxFileNoneReadingWriting = 1;
        } else {   //or create a new one
            if ( ( stiffMatrxFile = fopen(this->stiffMatrxFileName, "w") ) == NULL ) {
                OOFEM_ERROR("Can not create a new file %s\n", this->stiffMatrxFileName);
            }
            this->stiffMatrxFileNoneReadingWriting = 2;
        }
    }
#endif
    return LSpace :: initializeFrom(ir);;
}



/*Stiffness matrix is taken from the microproblem. No GPs are presented here on macroscale.
 * Stiffness matrix (24,24) goes in order node 1 (u,v,w) to node 2 (u,v,w) ... 8 (u,v,w). Displacements are in global coordinates.
 */
void MacroLSpace :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    //MatResponseMode rMode specifies tangent, secant, or initial matrix
    if ( !this->isActivated(tStep) ) {
        return;
    }


    //called the first time and initiates the microproblem
    if ( this->firstCall ) {
        this->microMaterial = static_cast< MicroMaterial * >( this->giveMaterial() );
        if ( this->microMaterial->microMatIsUsed == true ) {
            OOFEM_ERROR("Micromaterial is already used on another element. Only one micromaterial can be assigned to one macro element");
        }

        this->microDomain = this->microMaterial->problemMicro->giveDomain(1); //from engngm.h
        this->microEngngModel = this->microDomain->giveEngngModel();
        this->microEngngModel->setProblemScale(microScale); //set microScale attribute
        this->microEngngModel->checkProblemConsistency();
        this->microMaterial->init(); //from UnknownNumberingScheme()
        this->microMaterial->setMacroProperties(this->giveDomain(), this, this->microMasterNodes, this->microBoundaryNodes);
        this->firstCall = false;
    }

    //call microproblem
    //activeMStep = microMat->problemMicro->giveMetaStep(1);//->setNumberOfSteps(1);
    //activeMStep->giveMetaStepNumber();
    //microEngngModel->timer.startTimer(EngngModelTimer :: EMTT_AnalysisTimer);
    //microproblem must have the same actual time and zero time increment

    this->microEngngModel->giveCurrentStep()->setTargetTime( tStep->giveTargetTime() ); //adjust total time
    this->microEngngModel->giveCurrentStep()->setIntrinsicTime( tStep->giveIntrinsicTime() ); //adjust intrinsic time
    this->microEngngModel->giveCurrentStep()->setTimeIncrement(0.); //no time increment
    this->microEngngModel->initMetaStepAttributes( microEngngModel->giveCurrentMetaStep() ); //updates numerical method

    OOFEM_LOG_INFO( "\n** Assembling %s stiffness matrix of microproblem %p on macroElement %d, micTimeStep %d, micTime %f\n", __MatResponseModeToString(rMode), this->microMaterial->problemMicro, this->giveNumber(), this->microEngngModel->giveCurrentStep()->giveNumber(), this->microEngngModel->giveCurrentStep()->giveTargetTime() );

    //this->microEngngModel->solveYourselfAt( microEngngModel->giveCurrentStep() );
    //this->microEngngModel->terminate( microEngngModel->giveCurrentStep() );

    if ( tStep != this->lastStiffMatrixTimeStep ) {
        this->microMaterial->giveMacroStiffnessMatrix(answer, this->microEngngModel->giveCurrentStep(), rMode, this->microMasterNodes, this->microBoundaryNodes);
        this->stiffMatrix = answer;
        this->lastStiffMatrixTimeStep = tStep;
        OOFEM_LOG_INFO("** Assembled now\n\n");
    } else {
        answer = this->stiffMatrix;
        OOFEM_LOG_INFO("** Assembled previously in this time step\n\n");
    }
}


//assign values to DOF on the boundary according to definition on macrolspace and actual displacement stage
void MacroLSpace :: changeMicroBoundaryConditions(TimeStep *tStep)
{
    GeneralBoundaryCondition *GeneralBoundaryCond;
    Function *timeFunct;
    DynamicInputRecord ir_func, ir_bc;
    FloatArray n(8), answer, localCoords;
    double displ;
    FloatArray displ_x(8), displ_y(8), displ_z(8);
    int counter;

    //dofManArray has the node order as specified in input file
    for ( int i = 1; i <= this->giveNumberOfNodes(); i++ ) { //8 nodes
        //global displacements
        displ_x.at(i) = this->giveNode(i)->giveDofWithID(D_u)->giveUnknown(VM_Total, tStep);
        displ_y.at(i) = this->giveNode(i)->giveDofWithID(D_v)->giveUnknown(VM_Total, tStep);
        displ_z.at(i) = this->giveNode(i)->giveDofWithID(D_w)->giveUnknown(VM_Total, tStep);
    }

    //overrides the first load-time function to be a constant
    if ( !microDomain->giveNumberOfFunctions() ) {
        microDomain->resizeFunctions(1);
    }

    ir_func.setRecordKeywordField("constantfunction", 1);
    ir_func.setField(1.0, _IFT_ConstantFunction_f);
    if ( ( timeFunct = classFactory.createFunction("constantfunction", 1, microDomain) ) == NULL ) {
        OOFEM_ERROR("Couldn't create constant time function");
    }
    timeFunct->initializeFrom(& ir_func);
    microDomain->setFunction(1, timeFunct);


    /*assign to each boundary node the form "bc 3 # # #", set 0s on free nodes
     * modify bcList = corresponds to "nbc"
     */
    microDomain->clearBoundaryConditions(); //from domain.C
    microDomain->resizeBoundaryConditions( 3 * microBoundaryNodes.giveSize() ); //from domain.C

    counter = 1;
    for ( auto &DofMan : microDomain->giveDofManagers() ) { //go through all nodes on microDomain
        if ( microBoundaryNodes.contains( DofMan->giveGlobalNumber() ) ) { //if the node number is on boundary
            this->evalInterpolation( n, microMaterial->microMasterCoords, * DofMan->giveCoordinates() );
            //n.printYourself();

            for ( Dof *dof: *DofMan ) {
                this->microBoundaryDofManager.at(counter) = DofMan->giveGlobalNumber();
                dof->setBcId(counter);
                DofIDItem id = dof->giveDofID();
                displ = n.dotProduct( id == D_u ? displ_x : ( id == D_v ? displ_y : displ_z ) );
                ir_bc.setRecordKeywordField("boundarycondition", counter);
                ir_bc.setField(1, _IFT_GeneralBoundaryCondition_timeFunct);
                ir_bc.setField(displ, _IFT_BoundaryCondition_PrescribedValue);
                if ( ( GeneralBoundaryCond = classFactory.createBoundaryCondition("boundarycondition", counter, microDomain) ) == NULL ) {
                    OOFEM_ERROR("Couldn't create boundary condition.");
                }
                GeneralBoundaryCond->initializeFrom(& ir_bc);
                microDomain->setBoundaryCondition(counter, GeneralBoundaryCond);
                counter++;
            }
        } else {
            for ( Dof *dof: *DofMan ) {
                dof->setBcId(0);
            }
        }
    }
}

//obtain nodal forces from underlying microScale
//node numbering on element is in the same order as in the input file
//useUpdatedGpRecord=1 is used for printing of reactions
void MacroLSpace :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    //StructuralEngngModel *microStructuralEngngModel;
    FloatArray reactions, localCoords, n(8);
    DofManager *DofMan;
    double reactionForce;

    this->microEngngModel->giveCurrentStep()->setTargetTime( tStep->giveTargetTime() ); //adjust total time
    this->microEngngModel->giveCurrentStep()->setIntrinsicTime( tStep->giveIntrinsicTime() ); //adjust total time
    this->microEngngModel->giveCurrentStep()->setTimeIncrement(0.); //no time increment

    //OOFEM_LOG_INFO("*** useUpdatedGpRecord %d\n", useUpdatedGpRecord);

    if ( useUpdatedGpRecord ) { //printing of data
        answer = this->internalMacroForcesVector;
    } else {
        OOFEM_LOG_INFO( "\n*** Solving reactions %p of macroElement %d, micTimeStep %d, macIteration %d, micTime %f\n", this->microMaterial->problemMicro, this->giveNumber(), this->microEngngModel->giveCurrentStep()->giveNumber(), this->iteration, this->microEngngModel->giveCurrentStep()->giveTargetTime() );

        this->iteration++;
        this->changeMicroBoundaryConditions(tStep);

        this->microEngngModel->solveYourselfAt( this->microEngngModel->giveCurrentStep() );
        this->microEngngModel->updateYourself( this->microEngngModel->giveCurrentStep() );
        //this->microEngngModel->terminate( this->microEngngModel->giveCurrentStep() );
        StructuralEngngModel *microStructuralEngngModel = dynamic_cast< StructuralEngngModel * >(this->microEngngModel);

        //reaction vector contains contributions from unknownNumberingScheme
        microStructuralEngngModel->computeReaction(reactions, this->microEngngModel->giveCurrentStep(), 1);
        //reactions.printYourself();
        answer.resize(24);
        answer.zero();
        /*for ( i = 1; i <= this->giveNumberOfNodes(); i++ ) { //8 nodes
         * DofMan = microDomain->giveDofManager(i);
         */

        for ( int i = 1; i <= this->microBoundaryDofManager.giveSize() / 3; i++ ) { //Number of DoFManagers stored in triplets
            DofMan = microDomain->giveDofManager( this->microBoundaryDofManager.at(3 * i - 2) );
            this->evalInterpolation( n, microMaterial->microMasterCoords, * DofMan->giveCoordinates() );
            for ( int j = 1; j <= DofMan->giveNumberOfDofs(); j++ ) { //3
                reactionForce = reactions.at(3 * i + j - 3);
                for ( int k = 1; k <= 8; k++ ) {
                    answer.at(3 * k + j - 3) += reactionForce * n.at(k);
                }
            }
        }

        this->internalMacroForcesVector = answer;
        OOFEM_LOG_INFO("*** Reactions done\n", this->microMaterial->problemMicro);
    }

    //answer.printYourself();
    //OOFEM_ERROR("STOP");
}

void MacroLSpace :: evalInterpolation(FloatArray &answer, const std::vector< FloatArray > &coords, const FloatArray &gcoords)
{
    FloatArray localCoords;

    //this->interpolation.global2local(localCoords, coords, gcoords, 0.0);//returns even outside the element boundaries
    //this->interpolation.evalN(answer, localCoords, 0.0);
    this->interpolation.global2local( localCoords, gcoords, FEIVertexListGeometryWrapper(coords) ); //returns even outside the element boundaries
    this->interpolation.evalN( answer, localCoords, FEIVertexListGeometryWrapper(coords) );
}


// Updates the receiver at end of time step.
void MacroLSpace :: updateYourself(TimeStep *tStep)
{
    FloatArray answer;

    for ( auto &iRule: integrationRulesArray ) {
        iRule->updateYourself(tStep);
    }

    OOFEM_LOG_INFO("*** Updating macroelement\n");
    //set number of timestep so the microproblem counts from 1 to nsteps in each iteration but not totally
    //this->microEngngModel->giveCurrentStep()->setNumber(1);
    //recalculate microproblem
    //this->giveInternalForcesVector(answer, tStep, 0);

    this->microEngngModel->terminate( this->microEngngModel->giveCurrentStep() ); //perform output, VTK
    this->iteration = 1;
    this->microEngngModel->giveNextStep(); //time step used in ouput file name
}
} // end namespace oofem
