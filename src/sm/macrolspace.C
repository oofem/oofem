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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "macrolspace.h"
#include "micromaterial.h"
#include "lspace.h"
#include "material.h"
#include "domain.h"
#include "usrdefsub.h"
#include "structuralmaterial.h"
#include "oofem_terminate.h"
#include "oofemtxtdatareader.h"
#include "oofemtxtinputrecord.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

#ifdef __OOFEG
 #include "engngm.h"
 #include "metastep.h"
 #include "oofeggraphiccontext.h"
 #include "oofegutils.h"
 #include "conTable.h"
#endif

namespace oofem {
//derived from linear brick element
MacroLSpace :: MacroLSpace(int n, Domain *aDomain) : LSpace(n, aDomain)
{
    this->microMasterNodes.resize(0);
    this->microBoundaryNodes.resize(0);
    this->firstCall = true;
    microMaterial = NULL;
    microDomain = NULL;
    microEngngModel = NULL;
    this->iteration = 1;
    this->lastStiffMatrixTimeStep = NULL;
}

MacroLSpace :: ~MacroLSpace() {}


IRResultType MacroLSpace :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;              // Required by IR_GIVE_FIELD macro
    //IRResultType val;

    this->LSpace :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, this->microMasterNodes, IFT_MacroLspace_microMasterNodes, "micromasternodes");

    if ( this->microMasterNodes.giveSize() != 8 ) {
        OOFEM_ERROR("Need 8 master nodes from the microproblem defined on macroLspace element\n");
    }

    IR_GIVE_FIELD(ir, this->microBoundaryNodes, IFT_MacroLspace_microBoundaryNodes, "microboundarynodes");

    microBoundaryDofManager.resize( 3 * microBoundaryNodes.giveSize() );

    //val = IR_GIVE_OPTIONAL_FIELD2(ir, this->stiffMatrxFileName, IFT_MacroLspace_stiffMatrxFileName, "stiffmatrxfilename"); //Macro

    //     if( ir->hasField(IFT_MacroLspace_stiffMatrxFileName, "stiffmatrxfilename")){
    //       if (fopen(this->stiffMatrxFileName,"r") != NULL){//if the file exist
    //         stiffMatrxFile = fopen(this->stiffMatrxFileName,"r");
    //         this->stiffMatrxFileNoneReadingWriting=1;
    //       }
    //       else {//or create a new one
    //         if((stiffMatrxFile = fopen(this->stiffMatrxFileName,"w")) == NULL)
    //           OOFEM_ERROR2("Can not create a new file %s\n", this->stiffMatrxFileName);
    //         this->stiffMatrxFileNoneReadingWriting=2;
    //       }
    //     }
    return IRRT_OK;
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
        this->microMaterial = ( MicroMaterial * ) this->giveMaterial(); //from element.h
        if ( this->microMaterial->microMatIsUsed == true ) {
            OOFEM_ERROR("Micromaterial is already used on another element. Only one micromaterial can be assigned to one macro element\n");
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
    } else   {
        answer = this->stiffMatrix;
        OOFEM_LOG_INFO("** Assembled previously in this time step\n\n");
    }
}


//assign values to DOF on the boundary according to definition on macrolspace and actual displacement stage
void MacroLSpace :: changeMicroBoundaryConditions(TimeStep *tStep)
{
    //Domain *microDomain = problemMicro->giveDomain(1);
    //EngngModel *microEngngModel = microDomain->giveEngngModel();
    //Domain *domain = this->giveDomain();
    DofManager *DofMan;
    GeneralBoundaryCondition *GeneralBoundaryCond;
    LoadTimeFunction *LoadTimeFunct;
    ///@todo This part should use the DynamicInputRecord instead
    OOFEMTXTInputRecord *ir = new OOFEMTXTInputRecord();
    FloatArray n(8), answer, localCoords;
    char str [ OOFEM_MAX_LINE_LENGTH ];
    double displ;
    FloatArray displ_x(8), displ_y(8), displ_z(8);
    int i, j;
    int counter;

    //dofManArray has the node order as specified in input file
    for ( i = 1; i <= this->giveNumberOfNodes(); i++ ) { //8 nodes
        DofMan = microDomain->giveDofManager(i);
        //global displacements
        displ_x.at(i) = this->giveNode(i)->giveDof(1)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
        displ_y.at(i) = this->giveNode(i)->giveDof(2)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
        displ_z.at(i) = this->giveNode(i)->giveDof(3)->giveUnknown(EID_MomentumBalance, VM_Total, tStep);
    }

    //overrides the first load-time function to be a constant
    if ( !microDomain->giveNumberOfLoadTimeFunctions() ) {
        microDomain->resizeLoadTimeFunctions(1);
    }

    sprintf(str, "ConstantFunction 1 f(t) 1.0");
    ir->setRecordString(str);
    if ( ( LoadTimeFunct = CreateUsrDefLoadTimeFunctionOfType("constantfunction", 1, microDomain) ) == NULL ) {
        OOFEM_ERROR("MacroLSpace :: changeMicroBoundaryConditions - Couldn't create constant time function");
    }
    LoadTimeFunct->initializeFrom(ir);
    microDomain->setLoadTimeFunction(1, LoadTimeFunct);


    /*assign to each boundary node the form "bc 3 # # #", set 0s on free nodes
     * modify bcList = corresponds to "nbc"
     */
    microDomain->clearBoundaryConditions(); //from domain.C
    microDomain->resizeBoundaryConditions( 3 * microBoundaryNodes.giveSize() ); //from domain.C

    counter = 1;
    for ( i = 1; i <= microDomain->giveNumberOfDofManagers(); i++ ) { //go through all nodes on microDomain
        DofMan = microDomain->giveDofManager(i);
        if ( microBoundaryNodes.contains( DofMan->giveGlobalNumber() ) ) { //if the node number is on boundary
            this->evalInterpolation( n, microMaterial->microMasterCoords, * DofMan->giveCoordinates() );
            //n.printYourself();

            for ( j = 1; j <= 3; j++ ) {
                this->microBoundaryDofManager.at(counter) = DofMan->giveGlobalNumber();
                DofMan->giveDof(j)->setBcId(counter);
                displ = n.dotProduct( j == 1 ? displ_x : ( j == 2 ? displ_y : displ_z ) );
                sprintf(str, "boundarycondition %d loadtimefunction 1 prescribedvalue %e", counter, displ);
                //OOFEM_LOG_INFO("%s\n", str);
                ir->setRecordString(str);
                if ( ( GeneralBoundaryCond = CreateUsrDefBoundaryConditionOfType("boundarycondition", counter, microDomain) ) == NULL ) {
                    OOFEM_ERROR("MacroLSpace :: changeMicroBoundaryConditions - Couldn't create boundary condition.");
                }
                GeneralBoundaryCond->initializeFrom(ir);
                microDomain->setBoundaryCondition(counter, GeneralBoundaryCond);
                counter++;
            }
        } else   {
            for ( j = 1; j <= 3; j++ ) {
                DofMan->giveDof(j)->setBcId(0);
            }
        }
    }

    delete ir;
}

//obtain nodal forces from underlying microScale
//node numbering on element is in the same order as in the input file
//useUpdatedGpRecord=1 is used for printing of reactions
void MacroLSpace :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    int i, j, k;
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
    } else    {
        OOFEM_LOG_INFO( "\n*** Solving reactions %p of macroElement %d, micTimeStep %d, macIteration %d, micTime %f\n", this->microMaterial->problemMicro, this->giveNumber(), this->microEngngModel->giveCurrentStep()->giveNumber(), this->iteration, this->microEngngModel->giveCurrentStep()->giveTargetTime() );

        this->iteration++;
        this->changeMicroBoundaryConditions(tStep);

        this->microEngngModel->solveYourselfAt( this->microEngngModel->giveCurrentStep() );
        this->microEngngModel->updateYourself( this->microEngngModel->giveCurrentStep() );
        //this->microEngngModel->terminate( this->microEngngModel->giveCurrentStep() );
        //microStructuralEngngModel = ( StructuralEngngModel * ) &this->microEngngModel;
        StructuralEngngModel *microStructuralEngngModel = dynamic_cast< StructuralEngngModel * >(this->microEngngModel);


        //reaction vector contains contributions from unknownNumberingScheme
        microStructuralEngngModel->computeReaction(reactions, this->microEngngModel->giveCurrentStep(), 1);
        //reactions.printYourself();
        answer.resize(24);
        answer.zero();
        /*for ( i = 1; i <= this->giveNumberOfNodes(); i++ ) { //8 nodes
         * DofMan = microDomain->giveDofManager(i);
         */

        for ( i = 1; i <= this->microBoundaryDofManager.giveSize() / 3; i++ ) { //Number of DoFManagers stored in triplets
            DofMan = microDomain->giveDofManager( this->microBoundaryDofManager.at(3 * i - 2) );
            this->evalInterpolation( n, microMaterial->microMasterCoords, * DofMan->giveCoordinates() );
            for ( j = 1; j <= DofMan->giveNumberOfDofs(); j++ ) { //3
                reactionForce = reactions.at(3 * i + j - 3);
                for ( k = 1; k <= 8; k++ ) {
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

void MacroLSpace :: evalInterpolation(FloatArray &answer, const FloatArray **coords, const FloatArray &gcoords)
{
    FloatArray localCoords;

    //this->interpolation.global2local(localCoords, coords, gcoords, 0.0);//returns even outside the element boundaries
    //this->interpolation.evalN(answer, localCoords, 0.0);
    this->interpolation.global2local(localCoords, gcoords, FEIVertexListGeometryWrapper(8, coords)); //returns even outside the element boundaries
    this->interpolation.evalN(answer, localCoords, FEIVertexListGeometryWrapper(8, coords));
}


// Updates the receiver at end of time step.
void MacroLSpace :: updateYourself(TimeStep *tStep)
{
    int i;
    FloatArray answer;

    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        integrationRulesArray [ i ]->updateYourself(tStep);
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
