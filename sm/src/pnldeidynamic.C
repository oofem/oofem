/* $Header: /home/cvs/bp/oofem/sm/src/pnldeidynamic.C,v 1.6.4.2 2004/05/14 13:45:45 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2000   Borek Patzak                                       



         Czech Technical University, Faculty of Civil Engineering,
     Department of Structural Mechanics, 166 29 Prague, Czech Republic
                                                                               
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.                                                                              
*/

//
// file PNlDEIDynamic.cc
//

#include "pnldeidynamic.h"
#include "nlstructuralelement.h"
#include "nummet.h"
#include "ldltfact.h"
#include "timestep.h"
#include "element.h"
#include "node.h"
#include "elementside.h"
#include "dof.h"
#include "initial.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#include <math.h>
#endif
#include "cltypes.h"
#include "verbose.h"
#include "outputmanager.h"
#include "mathfem.h"

#ifdef __PARALLEL_MODE
#include "problemcomm.h"
//#include "domaincomm.h"
#endif


#define ZERO_REL_MASS  1.E-6

PNlDEIDynamic::  PNlDEIDynamic (int i, EngngModel* _master) : StructuralEngngModel (i,_master), massMatrix(), loadVector(), 
previousIncrementOfDisplacementVector(), displacementVector(), 
velocityVector(), accelerationVector(), internalForces() 
#ifdef __PARALLEL_MODE
//, remoteDofManList()
#endif
{
#ifdef __PARALLEL_MODE
 commMode = ProblemCommunicator::PC__UNKNOWN_MODE;
 nonlocalExt = 0;
 communicator = nonlocCommunicator = NULL;
 commBuff = NULL;
#endif
 ndomains = 1;
 initFlag = 1;
}


PNlDEIDynamic :: ~PNlDEIDynamic ()  {
  //delete massMatrix; 
  //delete loadVector;
  //delete previousIncrementOfDisplacementVector;
  //delete incrementOfDisplacementVector;
  //delete displacementVector;
  //delete velocityVector;
  //delete accelerationVector;
}

NumericalMethod* PNlDEIDynamic :: giveNumericalMethod (TimeStep* atTime)
// only one has reason for DEIDynamic 
//     - SolutionOfLinearEquations

{
/*  if (nMethod) return nMethod ;
  NumericalMethod* nm;
  nm = (NumericalMethod*) new LDLTFactorization (1,domain,this);
  nMethod = nm;
  return nm; */
  return NULL ;   // not necessary here - diagonal matrix is used-simple inversion
}

IRResultType
PNlDEIDynamic :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

 StructuralEngngModel :: initializeFrom (ir);

 IR_GIVE_FIELD (ir, dumpingCoef, IFT_PNlDEIDynamic_dumpcoef, "dumpcoef"); // C = dumpingCoef * M // Macro
 IR_GIVE_FIELD (ir, deltaT, IFT_PNlDEIDynamic_deltat, "deltat"); // Macro

 drFlag = 0;
 IR_GIVE_OPTIONAL_FIELD (ir, drFlag, IFT_PNlDEIDynamic_drflag, "drflag"); // Macro
 if (drFlag) {
   IR_GIVE_FIELD (ir, Tau, IFT_PNlDEIDynamic_tau, "tau"); 
   IR_GIVE_FIELD (ir, pyEstimate, IFT_PNlDEIDynamic_py, "py"); 
 }

#ifdef __PARALLEL_MODE
 if (ir->hasField (IFT_PNlDEIDynamic_nodecutmode, "nodecutmode")) commMode = ProblemCommunicator::PC__NODE_CUT;
 else if (ir->hasField (IFT_PNlDEIDynamic_elementcutmode, "elementcutmode")) commMode = ProblemCommunicator::PC__ELEMENT_CUT;
 else _error ("instanciateFrom: PNlDEIDynamicCommunicatorMode not specified");

 commBuff = new CommunicatorBuff (this->giveNumberOfProcesses());
 communicator = new ProblemCommunicator (this, commBuff, this->giveRank(), 
                                         this->giveNumberOfProcesses(), 
                                         this->commMode);
 
 if (ir->hasField (IFT_PNlDEIDynamic_nonlocalext, "nonlocalext")) {
   nonlocalExt = 1;
   nonlocCommunicator = new ProblemCommunicator (this, commBuff, this->giveRank(), 
                                                 this->giveNumberOfProcesses(), 
                                                 ProblemCommunicator::PC__REMOTE_ELEMENT_MODE);
 }

#endif

  return IRRT_OK;
}





double PNlDEIDynamic ::  giveUnknownComponent (EquationID chc, ValueModeType mode,
                        TimeStep* tStep, Domain* d, Dof* dof)
// returns unknown quantity like displaacement, velocity of equation eq
// This function translates this request to numerical method language
{
 int eq = dof->giveEquationNumber();
 if (eq == 0) _error ("giveUnknownComponent: invalid equation number");

   if (tStep != this->giveCurrentStep ()) {
    _error ("giveUnknownComponent: unknown time step encountered");
  return 0.;
  }
 

 if (chc != EID_MomentumBalance) {
  _error ("giveUnknownComponent: Unknown is of undefined CharType for this problem");
  return 0.;
 }

  switch (mode)
    {
    case VM_Total:
      return displacementVector.at(eq);

  case VM_Incremental:
   return previousIncrementOfDisplacementVector.at(eq);

    case VM_Velocity:
      return velocityVector.at(eq);

    case VM_Acceleration:
      return accelerationVector.at(eq);

    default:
      _error ("giveUnknownComponent: Unknown is of undefined type for this problem");
    }
  return 0.;
}

TimeStep* PNlDEIDynamic :: giveNextStep ()
{
  int istep =0;
  double totalTime = 0;
  StateCounterType counter = 1;

  delete previousStep;
  if (currentStep != NULL) {
    totalTime = currentStep->giveTime() + deltaT;
    istep     = currentStep->giveNumber() + 1   ;
  counter = currentStep->giveSolutionStateCounter() + 1;
 }
  previousStep = currentStep;
  currentStep = new TimeStep (istep,this,1,totalTime,deltaT,counter);
  // time and dt variables are set eq to 0 for staics - has no meaning

  return currentStep;
}



void PNlDEIDynamic :: solveYourself ()
{
  //this -> giveNumericalMethod ();     // can be awoided
 
#ifdef __PARALLEL_MODE
 // force equation numbering before setting up comm maps
  int neq = this -> giveNumberOfEquations (EID_MomentumBalance);
#ifdef __VERBOSE_PARALLEL
  OOFEM_LOG_INFO ("[process rank %d] neq is %d\n", this->giveRank(), neq);
#endif

 // set up communication patterns
 communicator->setUpCommunicationMaps (this);
 if (nonlocalExt) nonlocCommunicator->setUpCommunicationMaps (this);
 // init remote dofman list
 // this->initRemoteDofManList ();
#endif

  StructuralEngngModel::solveYourself();

}



void  PNlDEIDynamic :: solveYourselfAt (TimeStep* tStep) {
//
// creates system of governing eq's and solves them at given time step
//
// first assemble problem at current time step
  Domain* domain = this->giveDomain(1);
  int nDofs,neq ;
  int i,k,j,jj;
  int nman  = domain -> giveNumberOfDofManagers();
  double coeff,maxDt, maxOm= 0.;
  double prevIncrOfDisplacement, incrOfDisplacement;
  
  DofManager* node ;
  Dof* iDof ;
  
  neq = this -> giveNumberOfEquations (EID_MomentumBalance);
  
  if (initFlag) {
#ifdef VERBOSE
    OOFEM_LOG_INFO("Assembling mass matrix\n");
#endif
    
    //
    // assemble mass Matrix
    //
    this->computeMassMtrx (massMatrix, maxOm, tStep);
    //massMatrix.resize(neq);
    //for (i=1; i<=neq; i++) massMatrix.at(i) = 1.0;
    //maxOm = 0.05;
    if (drFlag) {
      // if Dynamic Relaxation assemble amplitude load vector 
      loadRefVector.resize (neq);
      loadRefVector.zero();
      
      this->computeLoadVector (loadRefVector, VM_Total, tStep);
      //this->assembleVectorFromElements (loadRefVector, tStep, EID_MomentumBalance, ElementForceLoadVector, VM_Total, domain) ;
      //this->assembleVectorFromDofManagers(loadRefVector, tStep, EID_MomentumBalance, NodalLoadVector, VM_Total, domain);
      
#ifdef __PARALLEL_MODE
      // compute the processor part of load vector norm pMp 
      this->pMp=0.0;
      double my_pMp=0.0, coeff = 1.0;
      int eqNum, ndofs, ndofman = domain->giveNumberOfDofManagers();
      dofManagerParallelMode dofmanmode;
      DofManager* dman;
      Dof* jdof;
      for (int dm=1; dm <=ndofman; dm++) {
        dman = domain->giveDofManager(dm);
        ndofs = dman->giveNumberOfDofs();
        dofmanmode = dman->giveParallelMode();
        // skip all remote and null dofmanagers
        coeff = 1.0;
        if ((dofmanmode == DofManager_remote) || ((dofmanmode == DofManager_null))) continue;
        else if (dofmanmode == DofManager_shared) coeff = 1./dman->givePartitionsConnectivitySize();
        //coeff *= coeff;

        // for shared nodes we add locally an average= 1/givePartitionsConnectivitySize()*contribution
        for (j=1; j<=ndofs; j++) {
          jdof = dman->giveDof(j);
          if (jdof->isPrimaryDof() && (eqNum = jdof->giveEquationNumber())) {
            my_pMp += coeff*loadRefVector.at(eqNum)*loadRefVector.at(eqNum)/massMatrix.at(eqNum);
          }
        }
      }
      
      // sum up the contributions from processors
      MPI_Allreduce (&my_pMp, &pMp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      for (i=1; i<= neq; i++) {
        pMp += loadRefVector.at(i)*loadRefVector.at(i)/massMatrix.at(i);
      }
#endif
      // solve for rate of loading process (parameter "c") (undamped system assumed)
      if (dumpingCoef < 1.e-3) {
        c = 3.0*this->pyEstimate/pMp/Tau/Tau;
      } else {
        c = this->pyEstimate*Tau*dumpingCoef*dumpingCoef*dumpingCoef/pMp/
          (-3.0/2.0+dumpingCoef*Tau+2.0*exp(-dumpingCoef*Tau)-0.5*exp(-2.0*dumpingCoef*Tau));
      }
    }
    
    
    
    
    initFlag = 0;
  }

  
  if (tStep->giveNumber() == giveNumberOfFirstStep()) {
 
    // 
    // special init step - compute displacements at tstep 0
    //
    displacementVector.                   resize (neq); displacementVector.zero();
    previousIncrementOfDisplacementVector.resize (neq); previousIncrementOfDisplacementVector.zero();
    velocityVector.                       resize (neq); velocityVector.zero();
    accelerationVector.                   resize (neq); accelerationVector.zero();
    
    for (j=1; j<= nman; j++) {
      node = domain->giveDofManager(j);
      nDofs = node->giveNumberOfDofs() ;
      
      for (k=1 ; k<=nDofs ; k++) {
        // ask for initial values obtained from 
        // bc (boundary conditions) and ic (initial conditions)
        // all dofs are expected to be  DisplacementVector type.
        iDof  =  node->giveDof(k);
        if (!iDof->isPrimaryDof()) continue;
        jj = iDof->giveEquationNumber () ;
        if (jj) {
          displacementVector.at(jj) = iDof->giveUnknown(EID_MomentumBalance,VM_Total,tStep) ;
          velocityVector.at(jj)     = iDof->giveUnknown(EID_MomentumBalance,VM_Velocity,tStep) ;
          // accelerationVector = iDof->giveUnknown(AccelerartionVector,tStep) ;
        }
      }
    }
    
    //
    // set-up numerical model 
    //
    
    // try to determine the best deltaT
    // PI = 3.1415926535897932384626383279; // PI =  3.1415926535897931160E0
    maxDt = 2.0/sqrt(maxOm) ;
    if (deltaT > maxDt) {
      // print reduced time step increment and minimum period Tmin
      OOFEM_LOG_RELEVANT("deltaT reduced to %e, Tmin is %e\n",maxDt, maxDt*M_PI);
      deltaT = maxDt;
      tStep -> setTimeIncrement (deltaT);
    } // end of if (tStep->giveNumber() == giveNumberOfFirstStep())
    
    for (j = 1; j <= neq; j++) {
      //   incrementOfDisplacementVector.at(j) = 
      //    velocityVector.at(j)*(deltaT); // becomes previous before used 
      //   displacementVector.at(j) -= incrementOfDisplacementVector.at(j);
      previousIncrementOfDisplacementVector.at(j) =  velocityVector.at(j)*(deltaT);
      displacementVector.at(j) -= previousIncrementOfDisplacementVector.at(j);
    }


#ifdef __PARALLEL_MODE
    //
    // init remote dofs if necessary
    // this->initializeRemoteDofs ();
#endif  
    return;
 }   // end of init step
 
#ifdef VERBOSE
 OOFEM_LOG_INFO("Assembling right hand side\n");
#endif
 
 for (i=1; i<= neq; i++) {
   //   previousIncrementOfDisplacementVector.at(i) = incrementOfDisplacementVector.at(i) ;  
   displacementVector.at(i) += previousIncrementOfDisplacementVector.at(i);
 }

#ifdef __PARALLEL_MODE
 //
 // init remote dofs if necessary
 // this->updateRemoteDofDisplacement ();
#endif  
 
 tStep-> incrementStateCounter();              // update solution state counter

#ifdef __PARALLEL_MODE
 this -> exchangeRemoteElementData ();         // exchange remote element data if necessary
#endif 


 // compute internal forces
 this -> giveInternalForces (internalForces, tStep);
 if (!drFlag) {

   // 
   // assembling the element part of load vector
   //
   this->computeLoadVector (loadVector, VM_Total, tStep);
   //
   // assembling additional parts of right hand side 
   //
   for (k=1; k<=neq; k++) {
     loadVector.at(k) -= internalForces.at(k);
   }
 } else {
   // dynamic relaxation 
   // compute load factor
   pt = 0.0;

#ifdef __PARALLEL_MODE
   double my_pt = 0.0, coeff = 1.0;
   int eqNum, ndofs, ndofman = domain->giveNumberOfDofManagers();
   dofManagerParallelMode dofmanmode;
   DofManager* dman;
   Dof* jdof;
   for (int dm=1; dm <=ndofman; dm++) {
     dman = domain->giveDofManager(dm);
     ndofs = dman->giveNumberOfDofs();
     dofmanmode = dman->giveParallelMode();
     // skip all remote and null dofmanagers
     coeff = 1.0;
     if ((dofmanmode == DofManager_remote) || ((dofmanmode == DofManager_null))) continue;
     else if (dofmanmode == DofManager_shared) coeff = 1./dman->givePartitionsConnectivitySize();
     //coeff *= coeff;

     // for shared nodes we add locally an average= 1/givePartitionsConnectivitySize()*contribution
     for (j=1; j<=ndofs; j++) {
       jdof = dman->giveDof(j);
       if (jdof->isPrimaryDof() && (eqNum = jdof->giveEquationNumber())) {
         my_pt += coeff*internalForces.at(eqNum)*loadRefVector.at(eqNum)/massMatrix.at(eqNum);
       }
     }
   }
   // sum up the contributions from processors
   MPI_Allreduce (&my_pt, &pt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
   for (k=1; k<=neq; k++) {
     pt += internalForces.at(k)*loadRefVector.at(k)/massMatrix.at(k);
   }
#endif
   pt = pt/pMp;
   if (dumpingCoef < 1.e-3) {
     pt += c*(Tau-tStep->giveTime())/Tau;
   } else {
     pt += c*(1.0-exp(dumpingCoef*(tStep->giveTime()-Tau)))/dumpingCoef/Tau;
   }

   loadVector.resize (this->giveNumberOfEquations(EID_MomentumBalance));
   for (k=1; k<=neq; k++) {
     loadVector.at(k) = pt*loadRefVector.at(k)-internalForces.at(k);
   }
   
   
   // compute relative error
   double err = 0.0;
#ifdef __PARALLEL_MODE
   double my_err = 0.0;

   for (int dm=1; dm <=ndofman; dm++) {
     dman = domain->giveDofManager(dm);
     ndofs = dman->giveNumberOfDofs();
     dofmanmode = dman->giveParallelMode();
     // skip all remote and null dofmanagers
     coeff = 1.0;
     if ((dofmanmode == DofManager_remote) || ((dofmanmode == DofManager_null))) continue;
     else if (dofmanmode == DofManager_shared) coeff = 1./dman->givePartitionsConnectivitySize();
     //coeff *= coeff;

     // for shared nodes we add locally an average= 1/givePartitionsConnectivitySize()*contribution
     for (j=1; j<=ndofs; j++) {
       jdof = dman->giveDof(j);
       if (jdof->isPrimaryDof() && (eqNum = jdof->giveEquationNumber())) {
         my_err += coeff*loadVector.at(eqNum)*loadVector.at(eqNum)/massMatrix.at(eqNum);
       }
     }
   }
   // sum up the contributions from processors
   MPI_Allreduce (&my_err, &err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else
   for (k=1; k<=neq; k++) {
     err = loadVector.at(k)*loadVector.at(k)/massMatrix.at(k);
   }
#endif
   err = err / (pMp*pt*pt);
   OOFEM_LOG_RELEVANT("Relative error is %e, loadlevel is %e\n", err, pt);
 }
  
 for (j=1; j<= neq; j++) {
   coeff =  massMatrix.at(j) ;
   loadVector.at(j) +=
     coeff * ((1./(deltaT*deltaT))-dumpingCoef*1./(2.*deltaT)) * 
     previousIncrementOfDisplacementVector.at(j) ;
 }  
 
  //
  // set-up numerical model 
  //
  /* it is not necesary to call numerical method
     approach used here is not good, but effective enough
     inverse of diagonal mass matrix is done here 
  */
  // 
  // call numerical model to solve arised problem - done localy here
  //
#ifdef VERBOSE
 OOFEM_LOG_RELEVANT("Solving [step number %8d, time %15e]\n", tStep->giveNumber(), tStep->giveTime());
#endif

  for (i=1; i<= neq; i++) {
  prevIncrOfDisplacement = previousIncrementOfDisplacementVector.at(i);
  incrOfDisplacement = loadVector.at(i)/
    (massMatrix.at(i)*(1./(deltaT*deltaT) + dumpingCoef/(2.*deltaT))) ;
  accelerationVector.at(i) = (incrOfDisplacement - prevIncrOfDisplacement)/(deltaT*deltaT);
  velocityVector.at(i)     = (incrOfDisplacement + prevIncrOfDisplacement)/(2.*deltaT);
  previousIncrementOfDisplacementVector.at(i) = incrOfDisplacement; // becomes previous
  //  displacementVector.at(i) += incrOfDisplacement; // update total displacements
/*
   incrementOfDisplacementVector.at(i) = loadVector.at(i)/
    (massMatrix.at(i)*(1./(deltaT*deltaT) + dumpingCoef/(2.*deltaT))) ;
   accelerationVector.at(i) = incrementOfDisplacementVector.at(i) - 
    previousIncrementOfDisplacementVector.at(i);
   velocityVector.at(i) = incrementOfDisplacementVector.at(i) + 
    previousIncrementOfDisplacementVector.at(i);
*/
  }
//  accelerationVector.times(1./(deltaT*deltaT));
//  velocityVector.times(1./(2.*deltaT)) ;

#ifdef __PARALLEL_MODE
  //
 // update remote dofs if necessary
 // this->updateRemoteDofs ();
#endif
 
 // update nodes, elements, etc.
 this->updateYourself(this->giveCurrentStep());

} 


void    PNlDEIDynamic :: updateYourself (TimeStep* stepN) 
{
 // updates internal state to reached one
 // all internal variables are directly updated by 
 // numerical method - void function here
 this->updateInternalState(stepN);
 StructuralEngngModel::updateYourself(stepN);
}



/*
void PNlDEIDynamic :: terminate (TimeStep* stepN)
{
  FILE * File;
  int j ;
  File = this -> giveDomain() -> giveOutputStream() ;

  fprintf (File,"\nOutput for time % .3le \n\n",stepN->giveTime());
  // fprintf (File,"\nOutput for time step number %d \n\n",stepN->giveNumber()+1);
  
  int nman   = domain->giveNumberOfDofManagers ();
 
  if (requiresUnknowsDictionaryUpdate()) {
  for( j=1;j<=nman;j++) {
   this->updateDofUnknownsDictionary(domain->giveDofManager(j),stepN) ;
  }
 }

 for( j=1;j<=nman;j++) {
   domain->giveDofManager(j) -> updateYourself(stepN) ;
   //domain->giveDofManager(j)->printOutputAt(File, stepN);
  }
  
  Element* elem;
  
  int nelem = domain->giveNumberOfElements ();
  for (j=1 ; j<=nelem ; j++) {
   elem = domain -> giveElement(j) ;
   elem -> updateYourself(stepN) ;
   //elem -> printOutputAt(File, stepN) ;
  }
   
 domain->giveOutputManager()->doDofManOutput (File, stepN);
 domain->giveOutputManager()->doElementOutput (File, stepN);
 
}
*/

void
PNlDEIDynamic :: computeLoadVector (FloatArray& answer, ValueModeType mode, TimeStep* stepN)
{
 Domain* domain = this->giveDomain(1);
 answer.resize (this->giveNumberOfEquations(EID_MomentumBalance));
 answer.zero();
 // 
 // assembling the nodal part of load vector
 //
 this->assembleVectorFromDofManagers(answer, stepN, EID_MomentumBalance, NodalLoadVector, mode, domain);


 
#ifdef __PARALLEL_MODE
 // nodal load for shared nodes is applied on all partitions; local value should be rescaled
 
 int dm, j, eqNum, ndofs, ndofman = domain->giveNumberOfDofManagers();
 double coeff;
 dofManagerParallelMode dofmanmode;
 DofManager* dman;
 Dof* jdof;
 for (dm=1; dm <=ndofman; dm++) {
   dman = domain->giveDofManager(dm);
   dofmanmode = dman->giveParallelMode();

   if (dofmanmode == DofManager_shared) {
     ndofs = dman->giveNumberOfDofs();
     coeff = 1./dman->givePartitionsConnectivitySize();
     // for shared nodes we add locally an average= 1/givePartitionsConnectivitySize()*contribution
     for (j=1; j<=ndofs; j++) {
       jdof = dman->giveDof(j);
       if (jdof->isPrimaryDof() && (eqNum = jdof->giveEquationNumber())) {
         answer.at(eqNum)*= coeff;
       }
     }
   }
 }
#endif

 // 
 // assembling the element part of load vector
 //
 this->assembleVectorFromElements (answer, stepN, EID_MomentumBalance, ElementForceLoadVector, mode, domain) ;


 // exchange contributions 

#ifdef __PARALLEL_MODE
#ifdef __VERBOSE_PARALLEL
 VERBOSEPARALLEL_PRINT("PNlDEIDynamic :: computeLoadVector","Packing load",this->giveRank());
#endif
 communicator->packAllData ((StructuralEngngModel*)this, &answer, &StructuralEngngModel::packLoad);

#ifdef __VERBOSE_PARALLEL
 VERBOSEPARALLEL_PRINT("PNlDEIDynamic :: computeLoadVector", "Exchange of load started",this->giveRank());
#endif
 communicator->initExchange (LoadExchangeTag);

#ifdef __VERBOSE_PARALLEL
 VERBOSEPARALLEL_PRINT("PNlDEIDynamic :: computeLoadVector","Receiving and unpacking of load started",this->giveRank());
#endif
 
 communicator->unpackAllData ((StructuralEngngModel*)this, &answer, &StructuralEngngModel::unpackLoad);
#endif
}


void
PNlDEIDynamic :: giveInternalForces (FloatArray& answer, TimeStep* stepN)
{
  // computes nodal representation of internal forces (real ones)
  // simply assembles contributions from each element in domain
 Domain* domain = this->giveDomain(1);
  Element* element;
  IntArray loc;
  FloatArray  charVec;
  //FloatArray* answer = new FloatArray(displacementVector->giveSize());
  int nelems;

 answer.resize (displacementVector.giveSize());
 answer.zero();

  nelems = domain-> giveNumberOfElements();
  for (int i = 1; i<= nelems; i++) {
    element = (NLStructuralElement*) domain->giveElement (i);

#ifdef __PARALLEL_MODE
  // skip remote elements (these are used as mirrors of remote eleemnts on other domains
  // when nonlocal constitutive models are used. Their introduction is necessary to
  // allow local averaging on domains without fine grain communication between domains).
  if (element->giveParallelMode () == Element_remote) continue;
#endif
    // if (!element -> hasNLCapability ()) {
    //   error ("giveInternalForces: element with no non-linear capability encountered\n");
    // }
    element -> giveLocationArray (loc, EID_MomentumBalance);
    element -> giveCharacteristicVector (charVec, NodalInternalForcesVector, VM_Total, stepN );
    if (charVec.containsOnlyZeroes ()) continue;
    answer.assemble (charVec, loc) ;
  }    

#ifdef __PARALLEL_MODE
// if (commMode == PNlDEIDynamicCommunicator__NODE_CUT) {
  // exchange Internal forces for  node cut mode
#ifdef __VERBOSE_PARALLEL
  VERBOSEPARALLEL_PRINT("PNlDEIDynamic :: giveInternalForces","Packing internal forces",this->giveRank());
#endif

  communicator->packAllData ((StructuralEngngModel*)this, &answer, &StructuralEngngModel::packInternalForces);

#ifdef __VERBOSE_PARALLEL
  VERBOSEPARALLEL_PRINT("PNlDEIDynamic :: giveInternalForces", "Exchange of internal forces started",this->giveRank());
#endif

  communicator->initExchange (InternalForcesExchangeTag);

#ifdef __VERBOSE_PARALLEL
  VERBOSEPARALLEL_PRINT("PNlDEIDynamic :: giveInternalForces","Receiving and unpacking internal forces started",this->giveRank());
#endif

  communicator->unpackAllData ((StructuralEngngModel*)this, &answer, &StructuralEngngModel::unpackInternalForces);
// }
#endif
 
  internalVarUpdateStamp = stepN ->giveSolutionStateCounter();
  return ;
}

void
PNlDEIDynamic :: computeMassMtrx (FloatArray& massMatrix, double& maxOm, TimeStep* tStep)
{
 Domain* domain = this->giveDomain(1);
  int nelem = domain -> giveNumberOfElements ();
  int neq = this -> giveNumberOfEquations (EID_MomentumBalance);
 int i,j,jj,n;
  double maxOmi, maxOmEl;
 FloatMatrix charMtrx, charMtrx2;
  //FloatArray diagonalStiffMtrx;
  IntArray loc ;
  Element* element;
#ifdef __PARALLEL_MODE
 int result;
#endif

 maxOm = 0.;
 massMatrix.resize (neq); massMatrix.zero();
 //diagonalStiffMtrx.resize (neq); diagonalStiffMtrx.zero();
 for (i = 1; i <= nelem ; i++ ) {
  element = domain -> giveElement(i);
#ifdef __PARALLEL_MODE
  // skip remote elements (these are used as mirrors of remote eleemnts on other domains
  // when nonlocal constitutive models are used. They introduction is necessary to
  // allow local averaging on domains without fine grain communication between domains).
  if (element->giveParallelMode () == Element_remote) continue;
#endif
  element -> giveLocationArray (loc, EID_MomentumBalance);
  element -> giveCharacteristicMatrix (charMtrx, LumpedMassMatrix, tStep );
  //charMtrx.beLumpedOf (fullCharMtrx);

  //delete fullCharMtrx;
  
  // ---> COMMENT REGION IF NO LOCAL VARIANT OF ZERO MASS REPLACENMENT IS NEEDED
  element -> giveCharacteristicMatrix (charMtrx2, StiffnessMatrix, tStep);
  // <--- END REGION
  //
  // assemble it manualy 
  //
#ifdef DEBUG
  if ((n=loc.giveSize()) != charMtrx.giveNumberOfRows()) {
   _error ("solveYourselfAt : dimension mismatch");
  }
#endif
  
  n = loc.giveSize();
  
  // ---> COMMENT REGION IF NO LOCAL VARIANT OF ZERO MASS REPLACENMENT IS NEEDED
  maxOmEl = 0.;
  double maxElmass = -1.0;
  for (j=1 ; j<=n; j++) maxElmass = max(maxElmass,charMtrx.at(j,j));
  if (maxElmass <= 0.0) {
    _warning2 ("solveYourselfAt: Element (%d) with zero (or negative) lumped mass encountered\n", i);
    //charMtrx.printYourself();
  } else {
    for (j=1 ; j<=n; j++) {
      if (charMtrx.at(j,j) > maxElmass*ZERO_REL_MASS) {
        maxOmi =  charMtrx2.at(j,j)/charMtrx.at(j,j) ;
        maxOmEl = (maxOmEl > maxOmi) ? (maxOmEl) : (maxOmi) ; 
        
      }
    }
  
    maxOm = (maxOm > maxOmEl) ? (maxOm) : (maxOmEl) ; 
  
    for (j=1 ; j<=n; j++) {
      jj = loc.at(j) ;
      if ((jj) && (charMtrx.at(j,j) <= maxElmass*ZERO_REL_MASS)) 
        charMtrx.at(j,j) = charMtrx2.at(j,j) / maxOmEl;
    }
    // <--- END REGION
  }
  
  for (j=1 ; j<=n; j++) {
   jj = loc.at(j) ;
   if (jj) { massMatrix.at(jj) += charMtrx.at(j,j) ;
        //diagonalStiffMtrx.at(jj) += charMtrx2.at(j,j);}
       }
    }
 }  
 /*
   // ---> UNCOMMENT REGION IF GLOBAL VARIANT OF ZERO MASS REPLACENMENT IS NEEDED
   
   // if init step - find minimun period of vibration in order to
   // determine maximal admisible time step
   // global variant
   for (i=1; i<=nelem; i++)
   {
   element = domain -> giveElement(i);
   element -> giveLocationArray (loc);
   element -> giveCharacteristicMatrix (charMtrx, StiffnessMatrix, tStep);
   n = loc.giveSize () ;
   for (j=1; j<=n; j++) {
   jj = loc.at(j);
   if (jj) {
   diagonalStiffMtrx.at(jj) += charMtrx.at(j,j);
   }
   }
   //delete charMtrx;
   }
   // find find minimun period of vibration
   // - global variant
   //
   double maxElmass = -1.0;
   for (j=1 ; j<=n; j++) maxElmass = max(maxElmass,charMtrx.at(j,j));
   if (maxElmass <= 0.0) _error ("solveYourselfAt: Element with zero (or negative) lumped mass encountered\n");

   for (j=1; j<= neq; j++) {
   if (massMatrix.at(j) > maxElmass*ZERO_REL_MASS) {
   maxOmi =  diagonalStiffMtrx.at(j)/massMatrix.at(j) ;
   maxOm = (maxOm > maxOmi) ? (maxOm) : (maxOmi) ; 
   }
   }
   // set ZERO MASS members in massMatrix to value which corresponds to 
   // maxOm
   // global variant
   //
   for (i=1; i<= neq; i++) {
   if (massMatrix.at(i) <= maxElmass*ZERO_REL_MASS) {
   massMatrix.at(i) = diagonalStiffMtrx.at(i) / maxOm;
   }
   }
   //delete diagonalStiffMtrx;
   // end global variant
   
   // <--- END REGION
   */
 
#ifdef __PARALLEL_MODE
// if (commMode == PNlDEIDynamicCommunicator__NODE_CUT) {
  // exchange Internal forces for  node and element cut mode
#ifdef __VERBOSE_PARALLEL
  VERBOSEPARALLEL_PRINT("PNlDEIDynamic :: computeMassMtrx","Packing masses",this->giveRank());
#endif

  communicator->packAllData (this, &PNlDEIDynamic::packMasses);

#ifdef __VERBOSE_PARALLEL
  VERBOSEPARALLEL_PRINT("PNlDEIDynamic :: computeMassMtrx","Mass exchangePacking started",this->giveRank());
#endif

  communicator->initExchange (MassExchangeTag);

#ifdef __VERBOSE_PARALLEL
  VERBOSEPARALLEL_PRINT("PNlDEIDynamic :: computeMassMtrx","Receiveng and Unpacking masses",this->giveRank());
#endif

  if (!communicator->unpackAllData (this, &PNlDEIDynamic::unpackMasses))
   _error ("PNlDEIDynamic :: computeMassMtrx: Receiveng and Unpacking masses failed");
// }

 // determine maxOm over all processes
#ifdef __USE_MPI
 double globalMaxOm;

#ifdef __VERBOSE_PARALLEL
 VERBOSEPARALLEL_PRINT("PNlDEIDynamic :: computeMassMtrx","Reduce of maxOm started",this->giveRank());
#endif

 result = MPI_Allreduce (&maxOm, &globalMaxOm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

#ifdef __VERBOSE_PARALLEL
 VERBOSEPARALLEL_PRINT("PNlDEIDynamic :: computeMassMtrx","Reduce of maxOm finished",this->giveRank());
#endif

 if (result != MPI_SUCCESS) _error ("setUpCommunicationMaps: MPI_Allreduce failed");
 maxOm = globalMaxOm;
#else
 WARNING: NOT SUPPORTED MESSAGE PARSING LIBRARY
#endif

#endif
}

/*
#ifdef __PARALLEL_MODE
void
PNlDEIDynamic :: updateRemoteDofs ()
{
#ifdef __PARALLEL_MODE
 if (commMode == PNlDEIDynamicCommunicator__ELEMENT_CUT) {
  // exchange remote dof unknowns for element cut mode
#ifdef __VERBOSE_PARALLEL
  VERBOSEPARALLEL_PRINT("PNlDEIDynamic :: updateRemoteDofs","Packing remote dofs unknowns",domain->giveRank());
#endif

  communicator->packAllData (&PNlDEIDynamic::packRemoteDofsUnknowns);

#ifdef __VERBOSE_PARALLEL
  VERBOSEPARALLEL_PRINT("PNlDEIDynamic :: updateRemoteDofs","Remote dofs unknowns exchange started",domain->giveRank());
#endif

  communicator->initExchange (RemoteDofsUnknwnExchangeTag);

#ifdef __VERBOSE_PARALLEL
  VERBOSEPARALLEL_PRINT("PNlDEIDynamic :: updateRemoteDofs","Receiveng and Unpacking remote dofs unknws",domain->giveRank());
#endif

  if (!communicator->unpackAllData (&PNlDEIDynamic::unpackAndUpdateRemoteDofsUnknowns))
   _error ("PNlDEIDynamic :: updateRemoteDofs: Receiveng and Unpacking remote dofs unknws failed");

 }
#endif
}
#endif
*/

#ifdef __PARALLEL_MODE

int
PNlDEIDynamic :: packMasses (ProcessCommunicator& processComm)
{
 int result = 1;
 int i, size;
 int j, ndofs, eqNum;
 Domain* domain = this->giveDomain(1);
 IntArray const* toSendMap = processComm.giveToSendMap();
 CommunicationBuffer* send_buff = processComm.giveSendBuff();
 DofManager* dman;
 Dof* jdof;

 size = toSendMap->giveSize();
 for (i=1; i<= size; i++) {
  dman = domain->giveDofManager (toSendMap->at(i));
  ndofs = dman->giveNumberOfDofs();
  for (j=1; j<=ndofs; j++) {
   jdof = dman->giveDof(j);
   if (jdof->isPrimaryDof() && (eqNum = jdof->giveEquationNumber())) {
    result &= send_buff->packDouble (massMatrix.at(eqNum));
   }
  }
 }
 return result;
}

int 
PNlDEIDynamic::unpackMasses (ProcessCommunicator& processComm)
{ 
 int result = 1;
 int i, size;
 int j, ndofs, eqNum;
 Domain* domain = this->giveDomain(1);
 dofManagerParallelMode dofmanmode;
 IntArray const* toRecvMap = processComm.giveToRecvMap();
 CommunicationBuffer* recv_buff = processComm.giveRecvBuff();
 DofManager* dman;
 Dof* jdof;
 double value;


 size = toRecvMap->giveSize();
 for (i=1; i<= size; i++) {
  dman = domain->giveDofManager (toRecvMap->at(i));
  ndofs = dman->giveNumberOfDofs();
  dofmanmode = dman->giveParallelMode();
  for (j=1; j<=ndofs; j++) {
   jdof = dman->giveDof(j);
   if (jdof->isPrimaryDof() && (eqNum = jdof->giveEquationNumber())) {
    result &= recv_buff->unpackDouble (value);
    if (dofmanmode == DofManager_shared) 
     massMatrix.at(eqNum) += value;
    else if (dofmanmode == DofManager_remote)
     massMatrix.at(eqNum)  = value;
    else _error ("unpackMasses: unknown dof namager parallel mode");
   }
  }
 }
 return result;
}



/*
int 
PNlDEIDynamic::packRemoteDofsUnknowns (PNlDEIDynamicDomainCommunicator& domainComm)
{ 
 int result = 1;
 int i, size;
 int j, ndofs;
 int eqNum;
 double value;
 IntArray const* toSendMap = domainComm.giveToSendMap();
 CommunicationBuffer* send_buff = domainComm.giveSendBuff();
 TimeStep* tStep = this->giveCurrentStep ();
 DofManager* dofman;
 
 size = toSendMap->giveSize();
 for (i=1; i<= size; i++) {
  dofman = domain->giveDofManager (toSendMap->at(i));
  ndofs = dofman->giveNumberOfDofs ();
  for (j=1; j<=ndofs; j++) {

   
   //  eqNum = dofman->giveDof (j) ->giveEquationNumber ();
   //  if (eqNum) value =  previousIncrementOfDisplacementVector.at(eqNum);
   //  else value  = 0.;
   //  send_buff->packDouble (value);
   

   result &= dofman->giveDof (j)->packUnknowns (*send_buff, DisplacementVector, IncrementalMode, tStep);

  }
 }
 return result;
}


int
PNlDEIDynamic::unpackAndUpdateRemoteDofsUnknowns (PNlDEIDynamicDomainCommunicator& domainComm)
{
 IntArray const* toRecvMap = domainComm.giveToRecvMap();
 CommunicationBuffer* recv_buff = domainComm.giveRecvBuff();
 TimeStep* tStep = this->giveCurrentStep ();
 DofManager* dofman;
 Dof* dof;
 int result = 1;
 int i, size;
 int j, ndofs;
 double prevIncrOfDisplacement, incrOfDisplacement, acceleration, velocity;

 size = toRecvMap->giveSize();
 for (i=1; i<= size; i++) {
  dofman = domain->giveDofManager (toRecvMap->at(i));
  ndofs = dofman->giveNumberOfDofs ();
  for (j=1; j<=ndofs; j++) {
   dof = dofman->giveDof (j);
   prevIncrOfDisplacement = dof->giveUnknown (DisplacementVector, IncrementalMode, tStep);
   // receive displacement increment

   //recv_buff->unpackDouble (incrOfDisplacement);

   result &= dof->unpackAndUpdateUnknown (*recv_buff, DisplacementVector, IncrementalMode, tStep);
   incrOfDisplacement = dof->giveUnknown (DisplacementVector, IncrementalMode, tStep);

   //prevDispl = dof->giveUnknown (DisplacementVector, TotalMode, tStep);
   // compute velocities and accelerations
   acceleration = (incrOfDisplacement - prevIncrOfDisplacement)/(deltaT*deltaT);
   velocity     = (incrOfDisplacement + prevIncrOfDisplacement)/(2.*deltaT);
   dof->updateUnknownsDictionary (tStep, DisplacementVector, IncrementalMode, incrOfDisplacement);
   dof->updateUnknownsDictionary (tStep, DisplacementVector, VelocityMode, velocity);
   dof->updateUnknownsDictionary (tStep, DisplacementVector, AccelerationMode, acceleration);
  }
 }
 return result;
}
*/

int 
PNlDEIDynamic :: estimateMaxPackSize (IntArray& commMap, CommunicationBuffer& buff, int packUnpackType) 
{
 int mapSize = commMap.giveSize();
 int i, j, ndofs, count = 0, pcount = 0;
 IntArray locationArray;
 Domain* domain = this->giveDomain(1);
 DofManager* dman;
 Dof* jdof;

 if (packUnpackType == ProblemCommunicator::PC__ELEMENT_CUT) {
  for (i=1; i<= mapSize; i++) {
   count += domain->giveDofManager (commMap.at(i))->giveNumberOfDofs();
  }
  return (buff.giveDoubleVecPackSize (1) * count);
 } else if (packUnpackType == ProblemCommunicator::PC__NODE_CUT) {
  for (i=1; i<= mapSize; i++) {
   ndofs = (dman = domain->giveDofManager (commMap.at(i)))->giveNumberOfDofs();
   for (j=1; j<=ndofs; j++) {
    jdof = dman->giveDof(j);
    if (jdof->isPrimaryDof() && (jdof->giveEquationNumber())) {
      count++;
    } else {
      pcount++;
    }
   }
  }

  //printf ("\nestimated count is %d\n",count);
  return (buff.giveDoubleVecPackSize (1) * max(count,pcount));
 } else  if (packUnpackType == ProblemCommunicator::PC__REMOTE_ELEMENT_MODE) {

  for (i=1; i<= mapSize; i++) {
   count += domain->giveElement (commMap.at(i))->estimatePackSize (buff);
  }
  return count;

 }

 return 0;
}

/*
int
PNlDEIDynamic :: initRemoteDofManList ()
{
 int i, pos = 0, size = 0, ndofman = domain->giveNumberOfDofManagers();

 if (commMode == PNlDEIDynamicCommunicator__ELEMENT_CUT) {

  for (i=1; i<=ndofman; i++) if (domain->giveDofManager(i)->giveParallelMode () == DofManager_remote) size++;
  this->remoteDofManList.resize (size);
  for (i=1; i<=ndofman; i++) 
   if (domain->giveDofManager(i)->giveParallelMode () == DofManager_remote) {
    this->remoteDofManList.at(++pos) = i;
   }
 }
 return 1;
}


void
PNlDEIDynamic :: initializeRemoteDofs ()
{
 TimeStep* tStep = this->giveCurrentStep ();
 int i,j,ndofs, size = remoteDofManList.giveSize();
 double displ, velocity;
 DofManager* dofman;
 Dof* dof;

 if (commMode == PNlDEIDynamicCommunicator__ELEMENT_CUT) {
  
  for (i=1; i<=size; i++) {
   dofman = domain->giveDofManager (remoteDofManList.at(i));
   ndofs = dofman->giveNumberOfDofs ();
   for (j=1; j<=ndofs; j++) {
    dof = dofman->giveDof (j);
    if (dof->hasIc ()) {
     displ    = dof->giveIc() -> give(TotalMode);
     velocity = dof->giveIc() -> give(VelocityMode);
    } else {
     displ = velocity = 0.;
    }
    dof->updateUnknownsDictionary (tStep, DisplacementVector, IncrementalMode, velocity*(deltaT));
    dof->updateUnknownsDictionary (tStep, DisplacementVector, TotalMode, displ-velocity*(deltaT));
   }
  }
 }
}



void 
PNlDEIDynamic :: updateRemoteDofDisplacement ()
{
 TimeStep* tStep = this->giveCurrentStep ();
 int i,j,ndofs, size = remoteDofManList.giveSize();
 double prevDispl, prevIncrOfDisplacement;
 DofManager* dofman;
 Dof* dof;

 if (commMode == PNlDEIDynamicCommunicator__ELEMENT_CUT) {
  for (i=1; i<=size; i++) {
   dofman = domain->giveDofManager (remoteDofManList.at(i));
   ndofs = dofman->giveNumberOfDofs ();
   for (j=1; j<=ndofs; j++) {
    dof = dofman->giveDof (j);
    prevDispl = dof->giveUnknown (DisplacementVector, TotalMode, tStep);
    prevIncrOfDisplacement = dof->giveUnknown (DisplacementVector, IncrementalMode, tStep);
    dof->updateUnknownsDictionary (tStep, DisplacementVector, TotalMode, prevDispl+prevIncrOfDisplacement);
   }
  }
 }
}
*/

int PNlDEIDynamic::exchangeRemoteElementData ()
{
 int result = 1;

 if (nonlocalExt) {
  
#ifdef __VERBOSE_PARALLEL
  VERBOSEPARALLEL_PRINT("PNlDEIDynamic :: exchangeRemoteElementData","Packing remote element data",this->giveRank());
#endif
  
  result &= nonlocCommunicator->packAllData ((StructuralEngngModel*)this,&StructuralEngngModel::packRemoteElementData);
  
#ifdef __VERBOSE_PARALLEL
  VERBOSEPARALLEL_PRINT("PNlDEIDynamic :: exchangeRemoteElementData","Remote element data exchange started",this->giveRank());
#endif
  
  result &= nonlocCommunicator->initExchange (RemoteElementsExchangeTag);
  
#ifdef __VERBOSE_PARALLEL
  VERBOSEPARALLEL_PRINT("PNlDEIDynamic :: exchangeRemoteElementData","Receiveng and Unpacking remote element data",this->giveRank());
#endif
  
  if (!(result&=nonlocCommunicator->unpackAllData ((StructuralEngngModel*)this,&StructuralEngngModel::unpackRemoteElementData)))
   _error ("PNlDEIDynamic :: exchangeRemoteElementData: Receiveng and Unpacking remote element data");
  // }
  
  return result;
 } // if (nonlocalext)

 return 1;
}
 


#endif



contextIOResultType PNlDEIDynamic :: saveContext (FILE* stream, void *obj)
// 
// saves state variable - displacement vector
//
{
 contextIOResultType iores;
 int closeFlag = 0;

  if (stream==NULL) {
  if (!this->giveContextFile(&stream, this->giveCurrentStep()->giveNumber(), 
                this->giveCurrentStep()->giveVersion(), contextMode_write)) 
   THROW_CIOERR(CIO_IOERR); // override 
  closeFlag = 1;
 }
 if ((iores = StructuralEngngModel :: saveContext (stream)) != CIO_OK) THROW_CIOERR(iores);

 if ((iores = previousIncrementOfDisplacementVector.storeYourself(stream)) != CIO_OK) THROW_CIOERR(iores);
 if ((iores = displacementVector.storeYourself(stream)) != CIO_OK) THROW_CIOERR(iores);
 if ((iores = velocityVector.storeYourself(stream)) != CIO_OK) THROW_CIOERR(iores);
 if ((iores = accelerationVector.storeYourself(stream)) != CIO_OK) THROW_CIOERR(iores);

 if (fwrite(&deltaT,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
 
  if (closeFlag) fclose (stream); // ensure consistent records
  return CIO_OK;
}



contextIOResultType PNlDEIDynamic :: restoreContext (FILE* stream, void *obj)
// 
// restore state variable - displacement vector
//
{
  contextIOResultType iores;
 int closeFlag = 0;
 int istep, iversion;

 this->resolveCorrespondingStepNumber (istep, iversion, obj);

 if (stream == NULL) {
  if (!this->giveContextFile(&stream, istep, iversion, contextMode_read)) 
   THROW_CIOERR(CIO_IOERR); // override 
  closeFlag = 1;
 }
 
 // save element context
 if ((iores = StructuralEngngModel :: restoreContext (stream, obj)) != CIO_OK) THROW_CIOERR(iores);  
 
 if ((iores = previousIncrementOfDisplacementVector.restoreYourself(stream)) != CIO_OK) THROW_CIOERR(iores);  
 if ((iores = displacementVector.restoreYourself(stream)) != CIO_OK) THROW_CIOERR(iores);  
 if ((iores = velocityVector.restoreYourself(stream)) != CIO_OK) THROW_CIOERR(iores);
 if ((iores = accelerationVector.restoreYourself(stream)) != CIO_OK) THROW_CIOERR(iores);

 if (fread (&deltaT,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);

  if (closeFlag) fclose (stream); // ensure consistent records
 return CIO_OK;
}



void 
PNlDEIDynamic::printDofOutputAt (FILE* stream, Dof* iDof, TimeStep* atTime) 
{
 static char dofchar[]="dva";
 static ValueModeType dofmodes[]={VM_Total, VM_Velocity, VM_Acceleration};
 
 iDof->printMultipleOutputAt(stream, atTime, dofchar, EID_MomentumBalance, dofmodes, 3);
}

void
PNlDEIDynamic:: terminate(TimeStep* tStep)
{
 StructuralEngngModel :: terminate (tStep);
 this->printReactionForces (tStep, 1);
}


void
PNlDEIDynamic :: printOutputAt (FILE* File,TimeStep* stepN) 
{
  //FILE* File = this -> giveDomain() -> giveOutputStream() ;

 if (!this->giveDomain(1)->giveOutputManager()->testTimeStepOutput (stepN)) return;  // do not print even Solution step header
 
  fprintf (File,"\n\nOutput for time % .3e, solution step number %d\n",stepN->giveTime(),stepN->giveNumber());  
  if (drFlag) {
    fprintf (File,"Reached load level : %e\n\n", this->pt );
  }

 this->giveDomain(1)->giveOutputManager()->doDofManOutput  (File, stepN);
 this->giveDomain(1)->giveOutputManager()->doElementOutput (File, stepN);
 
}
