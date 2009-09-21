/* $Header: /home/cvs/bp/oofem/sm/src/macrolspace.C,v 1.9 2009/09/20 13:04:00 vs Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   file MACROLSPACE.CC

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

//derived from linear brick element
MacroLSpace :: MacroLSpace(int n, Domain *aDomain) : LSpace(n, aDomain)
{
    this->microDOFs.resize(0);
    this->firstCall = true;
    microMaterial = NULL;
    microDomain = NULL;
    microEngngModel = NULL;
}

MacroLSpace :: ~MacroLSpace() { }


IRResultType MacroLSpace :: initializeFrom(InputRecord *ir)
{
    int i, j;
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;              // Required by IR_GIVE_FIELD macro

    this->LSpace :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, this->microMasterNodes, IFT_MacroLspace_microMasterNodes, "micromasternodes");

    //each node has three DOFs so the array microDOFs must be expanded
    for ( i = 1; i <= this->microMasterNodes.giveSize(); i++ ) {
        for ( j = 1; j <= 3; j++ ) {
            this->microDOFs.followedBy( ( this->microMasterNodes.at(i) - 1 ) * 3 + j );
        }
    }

    if ( microDOFs.giveSize() != 24 || microMasterNodes.giveSize() != 8 ) {
        OOFEM_ERROR("Need 8 nodes and 24 DOFs defined on macroLspace element\n");
    }

    IR_GIVE_FIELD(ir, this->microBoundaryNodes, IFT_MacroLspace_microBoundaryNodes, "microboundarynodes");

    return IRRT_OK;
}



//stiffness matrix [24x24] the rows (columns) goes in order: from node 1 (u,v,w) to node 2 (u,v,w) ... 8 (u,v,w), displacements are in global coordinates
//need to solve microproblem and then to assemble full stiffness matrix for static condensation
void MacroLSpace :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) {
    //int i, j;
    //GaussPoint *gp;
    //double dV;
    FloatMatrix bj, d, dbj;
    //IntegrationRule *iRule;
    //pointer to microproblem

    //Domain *domain = this->giveDomain();
    //MetaStep *activeMStep;
    IntArray condenseWhat;


    //answer.resize( computeNumberOfDofs(EID_MomentumBalance), computeNumberOfDofs(EID_MomentumBalance) );
    //answer.zero();

    if ( !this->isActivated(tStep) ) {
        return;
    }

    //answer.printYourself();

    //called the first time and initiates the microproblem
    if ( this->firstCall ) {
        if ( this->microMaterial != NULL ) {
            OOFEM_ERROR("Micromaterial already called on element. Only one micromaterial can be assigned to one macro element\n");
        }

        this->microMaterial = ( MicroMaterial * ) this->giveMaterial(); //from element.h
        this->microDomain = this->microMaterial->problemMicro->giveDomain(1); //from engngm.h
        this->microEngngModel = this->microDomain->giveEngngModel();
        this->microEngngModel->setProblemScale(microScale); //set microScale attribute
        this->microEngngModel->checkProblemConsistency();
        this->microMaterial->init();//from UnknownNumberingScheme(), obtain all DOFs and set totalNumberOfDomainEquation
        this->microMaterial->setMacroProperties(this->giveDomain(), this);
        this->firstCall = false;
    }



    //call microproblem
    //activeMStep = microMat->problemMicro->giveMetaStep(1);//->setNumberOfSteps(1);
    //activeMStep->giveMetaStepNumber();
    //microEngngModel->timer.startTimer(EngngModelTimer :: EMTT_AnalysisTimer);
    //microproblem must have the same actual time and zero time increment
    this->microEngngModel->giveNextStep();
    this->microEngngModel->giveCurrentStep()->setTime( tStep->giveTime() ); //adjust total time
    this->microEngngModel->giveCurrentStep()->setTimeIncrement(0.); //no time increment

    OOFEM_LOG_INFO( "\n** Solving microproblem %p on macroElement %d, step %d, time %f\n", this->microMaterial->problemMicro, this->giveNumber(), this->microEngngModel->giveCurrentStep()->giveNumber(), this->microEngngModel->giveCurrentStep()->giveTime() );

    this->microEngngModel->initMetaStepAttributes( microEngngModel->giveCurrentStep() );

    //this->changeMicroBoundaryConditions(tStep);
    //this->microEngngModel->solveYourselfAt( microEngngModel->giveCurrentStep() );
    //this->microEngngModel->terminate( microEngngModel->giveCurrentStep() );

    this->microMaterial->giveCondensedStiffnessMatrix(answer, this->microEngngModel->giveCurrentStep(), SecantStiffnessMatrix, this->microMasterNodes, this->microDOFs);

    //   iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    //   for ( j = 0; j <iRule->getNumberOfIntegrationPoints(); j++){
    //   gp = iRule->getIntegrationPoint(j);
    //   this->computeBmatrixAt(gp, bj);//returns with Jacobian in global c.s.
    //   this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);//in global c.s.
    //   dV = this->computeVolumeAround(gp);
    //   dbj.beProductOf(d, bj);
    //   answer.plusProductSymmUpper(bj, dbj, dV);
    //   }

    //answer is already symmetrized and resized from micromaterial model
    //   answer.symmetrized();
    //   answer.resize(24,24);
    OOFEM_LOG_INFO("** Microproblem %p solved with condensed secant stiffness matrix\n\n", this->microMaterial->problemMicro);
}


//assign values to DOF on the boundary according to definition on macrolspace and actual displacement stage
void MacroLSpace :: changeMicroBoundaryConditions(TimeStep *tStep) {
    //Domain *microDomain = problemMicro->giveDomain(1);
    //EngngModel *microEngngModel = microDomain->giveEngngModel();
    Domain *domain = this->giveDomain();
    DofManager *DofMan;
    GeneralBoundaryCondition *GeneralBoundaryCond;
    OOFEMTXTInputRecord *ir = new OOFEMTXTInputRecord();
    FloatArray n(8), answer, localCoords;
    //intArray microDofManArray(8);
    const FloatArray *globalPointCoords [ 8 ];
    char str [ 255 ];
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
        //globalPointCoords[i-1] = new const FloatArray(*this->giveNode(i)->giveCoordinates());
        j = this->microMasterNodes.at(i);
        globalPointCoords [ i - 1 ] = new const FloatArray( * this->microDomain->giveNode(j)->giveCoordinates() );
        //(*globalPointCoords[i-1]).printYourself();
    }

    //displ_x.printYourself();
    //this->interpolation.global2local(answer, globalPointCoords, *this->giveNode(i)->giveCoordinates(), 0.0);
    //answer.printYourself();

    /*assign to each boundary node the form "bc 3 # # #", set 0s on free nodes
     * modify bcList = corresponds to "nbc"
     */
    microDomain->clearBoundaryConditions(); //from domain.C
    microDomain->resizeBoundaryConditions( 3 * microBoundaryNodes.giveSize() ); //from domain.C

    counter = 1;
    for ( i = 1; i <= microDomain->giveNumberOfDofManagers(); i++ ) { //go through all nodes on microDomain
        DofMan = microDomain->giveDofManager(i);
        //         for ( j = 1; j <= DofMan->giveNumberOfDofs(); j++ ) { //j must be 1,2,3 - checked previously
        //             //if(DofMan->giveDof(j)->hasBc(tStep) != 0)
        //             //OOFEM_WARNING3("Deleting all boundary conditions from node %d component %d", i, j);
        //         }

        if ( microBoundaryNodes.contains( DofMan->giveGlobalNumber() ) ) { //if the node number is present in the microBoundaryNodes
            this->interpolation.global2local(localCoords, globalPointCoords, * DofMan->giveCoordinates(), 0.0);
            this->interpolation.evalN(n, localCoords, 0.0);
            //n.printYourself();

            for ( j = 1; j <= 3; j++ ) {
                DofMan->giveDof(j)->setBcId(counter);
                displ = dotProduct(n, j == 1 ? displ_x : ( j == 2 ? displ_y : displ_z ), 8);
                sprintf(str, "boundarycondition %d loadtimefunction 1 prescribedvalue %f", counter, displ);
                //printf("%s\n", str);
                ir->setRecordString(str);
                GeneralBoundaryCond = ( GeneralBoundaryCondition * ) ( GeneralBoundaryCondition(counter, microDomain).ofType("boundarycondition") );
                GeneralBoundaryCond->initializeFrom(ir);
                microDomain->setBoundaryCondition(counter, GeneralBoundaryCond);
                counter++;
            }
        } else {
            for ( j = 1; j <= 3; j++ ) {
                DofMan->giveDof(j)->setBcId(0);
            }
        }
    }





    delete ir;

    for ( i = 1; i <= this->giveNumberOfNodes(); i++ ) { //8 nodes
      delete globalPointCoords[i-1];
    }


    //   for(i=1; i<=microDomain->giveNumberOfDofManagers(); i++){
    //     DofMan=microDomain->giveDofManager(i);
    //     if(DofMan->giveNumberOfDofs() != 3)
    //       OOFEM_ERROR("Need 3 DoFs in each node\n");
    //     for(j=1; j<=3;j++){
    //       //DofMan->giveDof(j)->setPrescribedValue(1.1);
    //     }
    //   }
    //
}

//obtain nodal forces from underlying microScale
//everything is in the default element local coodrinates
void MacroLSpace :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) {
    OOFEM_LOG_INFO( "\n*** giveInternalForcesVector microproblem %p on macroElement %d, step %d, time %f\n", this->microMaterial->problemMicro, this->giveNumber(), this->microEngngModel->giveCurrentStep()->giveNumber(), this->microEngngModel->giveCurrentStep()->giveTime() );

    this->changeMicroBoundaryConditions(tStep);
    this->microEngngModel->solveYourselfAt( this->microEngngModel->giveCurrentStep() );
    this->microEngngModel->terminate( this->microEngngModel->giveCurrentStep() );




    OOFEM_LOG_INFO("\n*** giveInternalForcesVector microproblem %p done\n", this->microMaterial->problemMicro);
}


// Updates the receiver at end of time step.
void MacroLSpace :: updateYourself(TimeStep *tStep)
{
    int i;
    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        integrationRulesArray [ i ]->updateYourself(tStep);
    }

    //set number of timestep so the microproblem counts from 1 to nsteps in each iteration but not totally
    this->microEngngModel->giveCurrentStep()->setNumber(1);
}

