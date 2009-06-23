/* $Header: /home/cvs/bp/oofem/sm/src/lspace.C,v 1.8.4.1 2004/04/05 15:19:47 bp Exp $ */
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


MacroLSpace :: MacroLSpace(int n, Domain *aDomain) : LSpace(n, aDomain)
{
stiffnessMatrixMicro = NULL;//full stiffness matrix of microproblem
}


//stiffness matrix [24x24] the rows (columns) go in order: from node 1 (u,v,w) to node 2 (u,v,w) ... 8 (u,v,w) in global coordinates

void MacroLSpace :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep){
     int i,j;
     GaussPoint *gp;
     double dV;
     FloatMatrix bj,d,dbj;
     IntegrationRule *iRule;
     //pointer to microproblem
     MicroMaterial *microMat = ( MicroMaterial * ) this->giveMaterial();//from element.h
     Domain *microDomain = microMat->problemMicro->giveDomain(1);//from engngm.h
     EngngModel *microEngngModel = microDomain->giveEngngModel ();
     SparseMtrxType sparseMtrxType = (SparseMtrxType) 0;//?
     
     //stiffness matrix is done separately from another scale
     //if(strstr(this->giveClassName (),"Macro")) {
     //this->GiveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
      MetaStep *activeMStep;

  answer.resize( computeNumberOfDofs(EID_MomentumBalance), computeNumberOfDofs(EID_MomentumBalance) );
  answer.zero();
  
  if (!this->isActivated(tStep)) return;
  //answer.printYourself();
  
  //call microproblem
  //pointer to engngm (microMat->problemMicro->)
  microMat->problemMicro->setProblemScale(microScale);//set microScale attribute
  //activeMStep = microMat->problemMicro->giveMetaStep(1);//->setNumberOfSteps(1);
  //activeMStep->giveMetaStepNumber();
  //microEngngModel->timer.startTimer(EngngModelTimer :: EMTT_AnalysisTimer);
  //microproblem must have the same time increment
  microEngngModel->giveNextStep();
  microEngngModel->initMetaStepAttributes(microEngngModel->giveCurrentStep());
  microEngngModel->solveYourselfAt(microEngngModel->giveCurrentStep());
  
//   try {
      //microEngngModel->solveYourself();
//   } catch ( OOFEM_Terminate &c ) {
//         delete microMat->problemMicro;
//   }
  
  //sparseMtrxType = microMat->problemMicro->giveSparseMtrxType();
  
  
   if ( !stiffnessMatrixMicro ) {
             stiffnessMatrixMicro = :: CreateUsrDefSparseMtrx(sparseMtrxType);
   }
   
  stiffnessMatrixMicro->zero();
  stiffnessMatrixMicro->buildInternalStructure(microEngngModel, 1, EID_MomentumBalance);
  
  OOFEM_LOG_INFO("Assembling tangent stiffness matrix of microproblem at address %p\n", microMat->problemMicro);
  //microEngngModel->assemble( stiffnessMatrixMicro, tStep, EID_MomentumBalance, TangentStiffnessMatrix, microDomain );
  //stiffnessMatrixMicro.printYourself();
  //printf("%d",this->giveDomain()->giveEngngModel()->ndomains);
  //microMat->giveClassName();
  
  
//   iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
//   for ( j = 0; j <iRule->getNumberOfIntegrationPoints(); j++){
//   gp = iRule->getIntegrationPoint(j);
//   this->computeBmatrixAt(gp, bj);//returns with Jacobian in global c.s.
//   this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);//in global c.s.
//   dV = this->computeVolumeAround(gp);
//   dbj.beProductOf(d, bj);
//   answer.plusProductSymmUpper(bj, dbj, dV);
//   }
  //stiffnessMatrixMicro.resize(24,24);
  
  answer.symmetrized();
  answer.resize(24,24);
  answer.beUnitMatrix();
}

