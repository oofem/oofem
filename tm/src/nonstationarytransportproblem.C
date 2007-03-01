/* $Header: /home/cvs/bp/oofem/tm/src/nonstationarytransportproblem.C,v 1.2.4.1 2004/04/05 15:19:53 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2002   Borek Patzak                                       



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


#include "nonstationarytransportproblem.h"
#include "nummet.h"
#include "ldltfact.h"
#include "imlsolver.h"
#include "timestep.h"
#include "metastep.h"
#include "element.h"
#include "dofmanager.h"
#include "elementside.h"
#include "dof.h"
#include "cltypes.h"
#include "verbose.h"
#include "conTable.h"
#include "transportelement.h"
#include "usrdefsub.h"
#include "datastream.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif


NumericalMethod* NonStationaryTransportProblem :: giveNumericalMethod (TimeStep* atTime)
// only one has reason for LinearStatic 
//     - SolutionOfLinearEquations

{
  if (nMethod) return nMethod ;
 
 SparseLinearSystemNM* nm;
 if (solverType == ST_Direct) {
  nm = (SparseLinearSystemNM*) new LDLTFactorization (1,this->giveDomain(1),this);
  nMethod = nm;
  return nm;
 } else {
  nm = (SparseLinearSystemNM*) new IMLSolver (1,this->giveDomain(1),this);
  nMethod = nm;
  return nm;
 }
}

IRResultType
NonStationaryTransportProblem :: initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

 EngngModel::initializeFrom (ir);
 int val = 0;
 IR_GIVE_OPTIONAL_FIELD (ir, val, IFT_NonStationaryTransportProblem_lstype, "lstype"); // Macro
 solverType = (LinSystSolverType) val;

 val = 0;
 IR_GIVE_OPTIONAL_FIELD (ir, val, IFT_NonStationaryTransportProblem_smtype, "smtype"); // Macro
 sparseMtrxType = (SparseMtrxType) val;

 IR_GIVE_FIELD (ir, deltaT, IFT_NonStationaryTransportProblem_deltat, "deltat"); // Macro
 IR_GIVE_OPTIONAL_FIELD (ir, dtTimeFunction, IFT_NonStationaryTransportProblem_dtf, "dtf"); // Macro

 IR_GIVE_FIELD (ir, alpha, IFT_NonStationaryTransportProblem_alpha, "alpha"); // Macro
 /* The following done in updateAttributes
   if (this->giveNumericalMethod (giveCurrentStep())) nMethod -> instanciateFrom (ir);
 */
 // read lumped capacity stabilization flag
 if (ir->hasField(IFT_NonStationaryTransportProblem_lumpedcapa, "lumpedcapa")) lumpedCapacityStab = 1;

 // read field export flag
 exportFieldFlag = 0;
 if (ir->hasField(IFT_NonStationaryTransportProblem_exportfields, "exportfields")) {
  IntArray atomicFieldID;
  IR_GIVE_FIELD (ir, atomicFieldID, IFT_NonStationaryTransportProblem_atomicfields, "atomicfields"); // Macro
  // export flux fields
  FieldManager* fm = this->giveContext()->giveFieldManager();
  for (int i=1; i<=atomicFieldID.giveSize(); i++) {
   fm->registerField (&FluxField, (FieldType) atomicFieldID.at(i));
  }
 }

 return IRRT_OK;
}
 


double NonStationaryTransportProblem ::  giveUnknownComponent (EquationID chc, ValueModeType mode, 
                                TimeStep* tStep, Domain* d, Dof* dof)
// returns unknown quantity like displaacement, velocity of equation eq
// This function translates this request to numerical method language
{
 int eq = dof->giveEquationNumber();
 if (eq == 0) _error ("giveUnknownComponent: invalid equation number");

/*
  if (tStep != this->giveCurrentStep ()) {
    _error ("giveUnknownComponent: unknown time step encountered");
  return 0.;
  }
*/
  //if (chc != HeMaCVector) {// heat and mass concetration vector
 if (chc != EID_ConservationEquation) {// heat and mass concetration vector
  _error ("giveUnknownComponent: Unknown is of undefined CharType for this problem");
  return 0.;
  }
/*
 if ((tStep == this->giveCurrentStep()) && (solutionVector.isNotEmpty()))  {
  if (mode == UnknownMode_Total)  return solutionVector.at(eq);
  else if (mode == UnknownMode_Incremental) return solutionVector.at(eq)-previousSolutionVector.at(eq);
  else  _error ("giveUnknownComponent: Unknown is of undefined type for this problem");
 } else if ((tStep == this->givePreviousStep()) && (previousSolutionVector.isNotEmpty()))  {
  if (mode == UnknownMode_Incremental)  return previousSolutionVector.at(eq);
  else  _error ("giveUnknownComponent: Unknown is of undefined type for this problem");
 } else  _error ("giveUnknownComponent: Unknown is of undefined type for this problem");
*/
 return FluxField.giveUnknownValue (dof, mode, tStep) ;
}  


TimeStep* 
NonStationaryTransportProblem :: giveSolutionStepWhenIcApply()
{
 if (stepWhenIcApply == NULL) {
  stepWhenIcApply = new TimeStep (giveNumberOfTimeStepWhenIcApply(),this,0,
                  -deltaT,deltaT,0);
 }
 return stepWhenIcApply;
}


LoadTimeFunction*
NonStationaryTransportProblem :: giveDtTimeFunction ()
   // Returns the load-time function of the receiver.
{
   if (!dtTimeFunction || !ndomains) return NULL;
   return giveDomain(1) -> giveLoadTimeFunction(dtTimeFunction) ;
}

double
NonStationaryTransportProblem :: giveDeltaT (int n)
{
 if (giveDtTimeFunction()) return deltaT * giveDtTimeFunction()->__at(n);
 return deltaT;
}



TimeStep* 
NonStationaryTransportProblem :: giveNextStep ()
{
  int istep = this->giveNumberOfFirstStep();
  double totalTime = 0;
  StateCounterType counter = 1;
 delete previousStep;

  if (currentStep != NULL) {
   istep =  currentStep->giveNumber() + 1   ;
   totalTime = currentStep->giveTime() + giveDeltaT(istep);
   //istep =  currentStep->giveNumber() + 1   ;
   //totalTime = currentStep->giveTime() + deltaT;
  counter = currentStep->giveSolutionStateCounter() + 1;
 } else {
  // first step -> generate initial step
  currentStep = new TimeStep (*giveSolutionStepWhenIcApply());
 }
  previousStep = currentStep;
  currentStep = new TimeStep (istep,this, 1, totalTime, this->giveDeltaT(istep), counter);
  // time and dt variables are set eq to 0 for staics - has no meaning
  return currentStep;
}


void  NonStationaryTransportProblem :: solveYourselfAt (TimeStep* tStep) {
//
// creates system of governing eq's and solves them at given time step
//
// first assemble problem at current time step
  int neq =  this -> giveNumberOfEquations (EID_ConservationEquation);

 if (initFlag) {
  lhs = ::CreateUsrDefSparseMtrx(sparseMtrxType); 
  if (lhs==NULL) _error ("solveYourselfAt: sparse matrix creation failed");
  lhs->buildInternalStructure (this, 1, EID_ConservationEquation);

  bcRhs.resize(neq); bcRhs.zero();
  initFlag = 0;
 }

  if (tStep->giveNumber() == giveNumberOfFirstStep()) {
  TimeStep *stepWhenIcApply = tStep->givePreviousStep();

  this->applyIC (stepWhenIcApply);

  this->assembleVectorFromElements(bcRhs, stepWhenIcApply, EID_ConservationEquation, ElementBCTransportVector, VM_Total, this->giveDomain(1));
  this->assembleDirichletBcRhsVector (bcRhs, stepWhenIcApply, EID_ConservationEquation, VM_Total, NSTP_MidpointLhs, this->giveDomain(1));
  this->assembleVectorFromElements(bcRhs, stepWhenIcApply, EID_ConservationEquation, ElementInternalSourceVector, VM_Total, this->giveDomain(1));
  this->assembleVectorFromDofManagers(bcRhs, stepWhenIcApply, EID_ConservationEquation, NodalLoadVector, VM_Total, this->giveDomain(1)) ;

#ifdef VERBOSE
  OOFEM_LOG_INFO("Assembling conductivity and capacity matrices\n");
#endif
 
  this -> assemble (lhs, stepWhenIcApply, EID_ConservationEquation, LHSBCMatrix, this->giveDomain(1));
  lhs->times(alpha);
  this -> assemble (lhs, stepWhenIcApply, EID_ConservationEquation, NSTP_MidpointLhs, this->giveDomain(1));
 }

 FluxField.advanceSolution(tStep);
 FloatArray* solutionVector = FluxField.giveSolutionVector(tStep);
 solutionVector->resize(neq); solutionVector->zero();

#ifdef VERBOSE
  OOFEM_LOG_INFO("Assembling rhs\n");
#endif

  // 
  // assembling the element part of load vector
  //
 rhs = bcRhs; rhs.times(1.-alpha);
 bcRhs.zero();

  this->assembleVectorFromElements(bcRhs, tStep, EID_ConservationEquation, ElementBCTransportVector, VM_Total, this->giveDomain(1));
  this->assembleDirichletBcRhsVector (bcRhs, tStep, EID_ConservationEquation, VM_Total, NSTP_MidpointLhs, this->giveDomain(1));
  this->assembleVectorFromElements(bcRhs, tStep, EID_ConservationEquation, ElementInternalSourceVector, VM_Total, this->giveDomain(1));
  // 
  // assembling the nodal part of load vector
  //
  this->assembleVectorFromDofManagers(bcRhs, tStep, EID_ConservationEquation, NodalLoadVector, VM_Total, this->giveDomain(1)) ;
  for (int i=1; i<=neq; i++) rhs.at(i)+=bcRhs.at(i)*alpha;
  //
  // add the rhs part depending on previous solution
  //
  assembleAlgorithmicPartOfRhs (rhs, EID_ConservationEquation, tStep->givePreviousStep());
  //
  // set-up numerical model
  //
  this->giveNumericalMethod(tStep);

  // 
  // call numerical model to solve arised problem
  //
#ifdef VERBOSE
  OOFEM_LOG_INFO("Solving ...\n");
#endif


  //nMethod -> solveYourselfAt(tStep);
  nMethod -> solve (lhs, &rhs, FluxField.giveSolutionVector(tStep));
  // update solution state counter
  tStep->incrementStateCounter();   

  // update nodes, elements, etc.
  this->updateYourself(this->giveCurrentStep());

} 

void    
NonStationaryTransportProblem :: updateYourself (TimeStep* stepN) 
{
 //this->updateInternalState(stepN);
 EngngModel::updateYourself(stepN);
 //previousSolutionVector = solutionVector;
}


contextIOResultType 
NonStationaryTransportProblem :: saveContext (DataStream* stream, ContextMode mode, void *obj)
// 
// saves state variable - displacement vector
//
{
 contextIOResultType iores;
 int closeFlag = 0;
 FILE* file;

 if (stream==NULL) {
  if (!this->giveContextFile(&file, this->giveCurrentStep()->giveNumber(), 
                this->giveCurrentStep()->giveVersion(), contextMode_write)) 
   THROW_CIOERR(CIO_IOERR); // override 
  stream = new FileDataStream (file);
  closeFlag = 1;
 }

  if ((iores = EngngModel :: saveContext (stream,mode)) != CIO_OK) THROW_CIOERR(iores);
  if ((iores = FluxField.saveContext(stream,mode)) != CIO_OK) THROW_CIOERR(iores);

 if (closeFlag) {fclose (file); delete stream; stream=NULL;}// ensure consistent records
  return CIO_OK;
}



contextIOResultType 
NonStationaryTransportProblem :: restoreContext (DataStream* stream, ContextMode mode, void *obj)
// 
// restore state variable - displacement vector
//
{
 contextIOResultType iores;
 int closeFlag = 0;
 int istep, iversion;
 FILE* file;

 this->resolveCorrespondingStepNumber (istep, iversion, obj);

 if (stream == NULL) {
  if (!this->giveContextFile(&file, istep, iversion, contextMode_read))
   THROW_CIOERR(CIO_IOERR); // override 
  stream = new FileDataStream (file);
  closeFlag = 1;
 }

 if ((iores = EngngModel :: restoreContext (stream, mode, obj)) != CIO_OK) THROW_CIOERR(iores);
 if ((iores = FluxField.restoreContext(stream, mode)) != CIO_OK) THROW_CIOERR(iores);
 
 if (closeFlag) {fclose (file); delete stream; stream=NULL;}// ensure consistent records
 return CIO_OK;
}


int
NonStationaryTransportProblem::checkConsistency ()
{
// check internal consistency
// if success returns nonzero
 int i, nelem;
 Element* ePtr;
 TransportElement* sePtr;
 Domain* domain = this->giveDomain(1);

 nelem = domain->giveNumberOfElements();
 // check for proper element type

 for (i=1; i<= nelem; i++) {
  ePtr = domain->giveElement(i);
  sePtr = dynamic_cast<TransportElement*>(ePtr);
  if (sePtr == NULL) {
    _warning2 ("Element %d has no TransportElement base",i);
   return 0;
  }
 }

 EngngModel :: checkConsistency ();

 return 1;
}


void
NonStationaryTransportProblem::updateDomainLinks ()
{
 EngngModel::updateDomainLinks();
 this->giveNumericalMethod(giveCurrentStep())->setDomain (this->giveDomain(1));
}

void 
NonStationaryTransportProblem::giveElementCharacteristicMatrix (FloatMatrix& answer, int num, 
                                CharType type, TimeStep* tStep, Domain *domain) 
{
  // we don't directlt call element ->GiveCharacteristicMatrix() function, because some
   // engngm classes may require special modification of base types supported on
   // element class level

 if ((type == NSTP_MidpointLhs) || (type == NSTP_MidpointRhs)) {
  
    Element* element;
    // IntArray loc ;
    FloatMatrix charMtrx1, charMtrx2;

  element = domain -> giveElement(num);
  // element -> giveLocationArray (loc);
  element -> giveCharacteristicMatrix (answer, ConductivityMatrix, tStep );
  element -> giveCharacteristicMatrix (charMtrx2, CapacityMatrix, tStep);

  if (lumpedCapacityStab) {
   int i,j,size = charMtrx2.giveNumberOfRows();
   double s;
   for (i=1; i<=size; i++) {
    s = 0.0;
    for (j=1; j<=size; j++) {
     s+=charMtrx2.at(i,j);
     charMtrx2.at(i,j) = 0.0;
    }
    charMtrx2.at(i,i) = s;
   }
  }

  if (type == NSTP_MidpointLhs) {
   answer.times(this->alpha);
   charMtrx2.times(1./tStep->giveTimeIncrement());
  } else {
   answer.times(this->alpha-1.0);
   charMtrx2.times(1./tStep->giveTimeIncrement());
  }
  answer.plus (charMtrx2);
  return ;

 } else 
  EngngModel::giveElementCharacteristicMatrix(answer, num, type, tStep, domain);
}


void
NonStationaryTransportProblem::assembleAlgorithmicPartOfRhs (FloatArray& answer, EquationID ut, TimeStep* tStep)
{
  int i ;
  IntArray loc ;
  FloatMatrix charMtrx, bcMtrx ;
  FloatArray  unknownVec, contrib;
  Element *element ;
  
  Domain* domain = this->giveDomain(1);
  int nelem = domain -> giveNumberOfElements ();

  for (i = 1; i <= nelem ; i++ ) {
    element = domain -> giveElement(i);
#ifdef __PARALLEL_MODE
    // skip remote elements (these are used as mirrors of remote eleemnts on other domains
    // when nonlocal constitutive models are used. They introduction is necessary to
    // allow local averaging on domains without fine grain communication between domains).
    if (element->giveParallelMode () == Element_remote) continue;
#endif
    element -> giveLocationArray (loc, ut);
    this -> giveElementCharacteristicMatrix (charMtrx, i, NSTP_MidpointRhs, tStep, domain);
    element -> giveCharacteristicMatrix (bcMtrx, LHSBCMatrix, tStep);
    bcMtrx.times(this->alpha-1.0);
    if (bcMtrx.isNotEmpty()) charMtrx.plus(bcMtrx);
    if (charMtrx.isNotEmpty()) {
      element -> computeVectorOf (EID_ConservationEquation, VM_Total, tStep, unknownVec);
      contrib.beProductOf(charMtrx, unknownVec);
      answer.assemble (contrib, loc) ;
    }  
  }
}


void
NonStationaryTransportProblem::printDofOutputAt (FILE* stream, Dof* iDof, TimeStep* atTime) 
{
 iDof->printSingleOutputAt(stream, atTime, 'f', EID_ConservationEquation, VM_Total);
}



void 
NonStationaryTransportProblem :: applyIC (TimeStep* stepWhenIcApply)
{
 Domain* domain = this->giveDomain(1);
  int neq =  this -> giveNumberOfEquations (EID_ConservationEquation);
 FloatArray* solutionVector;
  
#ifdef VERBOSE
 OOFEM_LOG_INFO("Applying initial conditions\n");
#endif
 int nDofs,j,k,jj;
 int nman  = domain -> giveNumberOfDofManagers();
 DofManager *node;
 Dof  *iDof;
  
 FluxField.advanceSolution(stepWhenIcApply);
 solutionVector = FluxField.giveSolutionVector(stepWhenIcApply);
 solutionVector->resize(neq); solutionVector->zero();
    
 for (j=1; j<= nman; j++) {
  node = domain->giveDofManager(j);
  nDofs = node->giveNumberOfDofs() ;
  
  for (k=1 ; k<=nDofs ; k++) {
   // ask for initial values obtained from 
   // bc (boundary conditions) and ic (initial conditions)
   iDof  =  node->giveDof(k);
   if (!iDof->isPrimaryDof()) continue;
   jj = iDof->giveEquationNumber () ;
   if (jj) {
    solutionVector->at(jj) = iDof->giveUnknown(EID_ConservationEquation,VM_Total,stepWhenIcApply);
   }
  }
 }
 
 /* Not relevant in linear case
 // update element state according to given ic
 int nelem = domain -> giveNumberOfElements ();
 TransportElement* element;

 for (j = 1; j <= nelem ; j++ ) {
  element = (TransportElement*) domain -> giveElement(j);
  element -> updateInternalState (stepWhenIcApply);
  element -> updateYourself(stepWhenIcApply);
 }
 */
}


void
NonStationaryTransportProblem :: assembleDirichletBcRhsVector (FloatArray& answer, TimeStep* tStep, EquationID ut, 
                                                               ValueModeType mode, CharType lhsType, Domain* d)
{
  int ielem ;
  IntArray loc ;
  Element *element ;
  FloatArray  rp, charVec;
  FloatMatrix s,bcMtrx;
  
  int nelem = d -> giveNumberOfElements ();

  for (ielem = 1; ielem <= nelem ; ielem++ ) {
    element = d -> giveElement(ielem);
    
    element -> computeVectorOfPrescribed(EID_ConservationEquation, mode,tStep, rp) ;
    if (rp.containsOnlyZeroes())
      continue;
    else {
      this->giveElementCharacteristicMatrix (s, ielem, lhsType, tStep, d);
      element -> giveCharacteristicMatrix (bcMtrx, LHSBCMatrix, tStep);
      s.plus(bcMtrx);
      charVec.beProductOf (s, rp); charVec.negated() ;   
      
      element -> giveLocationArray (loc, ut);
      answer.assemble (charVec, loc) ;
    }
  } // end element loop
}

#ifdef __PETSC_MODULE
void 
NonStationaryTransportProblem::initPetscContexts ()
{
  int i;
  PetscContext *petscContext;

  petscContextList -> growTo(ndomains) ;
  for (i=0; i < this->ndomains ; i++) {
    petscContext =  new PetscContext (this, EID_ConservationEquation);
    petscContextList->put(i+1,petscContext) ;
    
  }
}
#endif
