/* $Header: /home/cvs/bp/oofem/oofemlib/src/engngm.C,v 1.45.4.3 2004/05/28 10:36:48 bp Exp $ */
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
// file: engngm.cc
//

#include "nummet.h"
#include "engngm.h"
#include "timestep.h"
#include "metastep.h"
#include "element.h"
#include "cltypes.h"
#include "mathfem.h"
#include "clock.h"
#include "datastream.h"
//#include "linearstatic.h"
//#include "nlinearstatic.h"
//#include "eigenvaluedynamic.h"
//#include "deidynamic.h"
//#include "nldeidynamic.h"
//#include "diidynamic.h"
//#include "creeplinearstatic.h"
//#include "stationaryflow.h"

#include "femcmpnn.h"
#include "dofmanager.h"
#include "node.h"
#include "elementside.h"
#include "dof.h"
#include "timestep.h"
#include "skyline.h"
#include "verbose.h"
#include "outputmanager.h"
#include "oofemtxtdatareader.h"
#include "util.h"
#include "oofemdef.h"
#include "sloangraph.h"
#include "logger.h"
#include "errorestimator.h"
#ifndef __MAKEDEPEND
#include <limits.h>
#include <string.h>
#include <stdio.h>
#include <stdarg.h>
// include unistd.h; needed for access
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#elif _MSC_VER
#include <io.h>
#endif
#endif

#ifdef __OOFEG  
#include "oofeggraphiccontext.h"
#endif

#ifdef TIME_REPORT
#ifndef __MAKEDEPEND
#include <time.h>
#endif
#endif

#ifdef __PETSC_MODULE
#ifndef __MAKEDEPEND
#include "petscvec.h"
#endif 
#endif

EngngModel :: EngngModel (int i, EngngModel* _master) : domainNeqs(), domainPrescribedNeqs()
// constructor
{
 number = i;
  //nMethod = NULL; 
  currentStep = NULL ; 
  previousStep = NULL; 
 stepWhenIcApply = NULL;
  numberOfSteps = 0;
  numberOfEquations = 0;
  numberOfPrescribedEquations = 0;
 renumberFlag      = 0;
 equationNumberingCompleted = 0;
 ndomains = 0;
 nMetaSteps = 0;

 //dataInputFileName     = NULL ;
 dataOutputFileName    = new char [MAX_FILENAME_LENGTH];
 //inputStream           = NULL ;
 outputStream          = NULL ;

 domainList            = new AList<Domain>(0);
 metaStepList          = new AList<MetaStep>(0);

 contextOutputMode     =  NOCONTEXT ;
 contextOutputStep     = 0 ;
 pMode                 = _processor;     // for giveContextFile()

 exportModuleManager   = new ExportModuleManager(this);
 master                = _master; // master mode by default
 // create context if in master mode; otherwise request context from master
 if (master) context = master->giveContext (); else context = new EngngModelContext();

#ifdef __PARALLEL_MODE
 initParallel ();
 parallelFlag = 1;
 loadBallancingFlag = false;
 force_load_rebalance_in_first_step = false; 
 lb = NULL;
 lbm = NULL;
#else
 parallelFlag=0;
#endif

#ifdef __PETSC_MODULE
 petscContextList            = new AList<PetscContext>(0);
#endif

}

EngngModel :: EngngModel (int i, char* s, EngngModel* _master) : domainNeqs(), domainPrescribedNeqs()
// constructor
{
 number = i;
  //nMethod = NULL; 
  currentStep = NULL ; 
  previousStep = NULL; 
 stepWhenIcApply = NULL;
  numberOfSteps = 0;
  numberOfEquations = 0;
  numberOfPrescribedEquations = 0;
 renumberFlag      = 0;
 equationNumberingCompleted = 0;
 ndomains = 0;
 nMetaSteps = 0;

   //dataInputFileName = new char[strlen(s)+1] ;
   //strcpy (dataInputFileName,s) ;
   dataOutputFileName = new char [MAX_FILENAME_LENGTH];

 //inputStream           = NULL ;
 outputStream          = NULL ;

 domainList            = new AList<Domain>(0);
 metaStepList          = new AList<MetaStep>(0);

 exportModuleManager   = new ExportModuleManager(this);
 master                = _master; // master mode by default
 // create context if in master mode; otherwise request context from master
 if (master) context = master->giveContext (); else context = new EngngModelContext();

#ifndef __PARALLEL_MODE
 //dataInputFileName = new char[strlen(s)+1] ;
  //strcpy (dataInputFileName,s) ;
 dataOutputFileName = new char [MAX_FILENAME_LENGTH];
 parallelFlag=0;
#else
 parallelFlag=1;
 initParallel ();
 //dataInputFileName = new char[strlen(s)+10] ;
 //sprintf (dataInputFileName, "%s.%d", s, rank);
 dataOutputFileName = new char [MAX_FILENAME_LENGTH];
 loadBallancingFlag = false;
 force_load_rebalance_in_first_step = false;
 lb = NULL;
 lbm = NULL;
#endif

#ifdef __PETSC_MODULE
 petscContextList            = new AList<PetscContext>(0);
#endif

}


EngngModel ::  ~EngngModel () 
// destructor
{
  delete currentStep; 
  delete previousStep; 
 delete stepWhenIcApply;

  //delete nMethod;

 //delete dataInputFileName ;
 delete[] dataOutputFileName;
 delete domainList;
 delete metaStepList;

#ifdef __PETSC_MODULE
 delete petscContextList;
#endif

 if (exportModuleManager) delete  exportModuleManager;
 // master deletes the context
 if (master==NULL) delete context;

 //fclose (inputStream) ;
 if (outputStream)
  fclose(outputStream) ;

#ifdef __PARALLEL_MODE
 if (loadBallancingFlag) {
   if (lb) delete lb;
   if (lbm) delete lbm;
 }
#endif
}


int EngngModel :: instanciateYourself (DataReader* dr, InputRecord* ir, char* dataOutputFileName, char* desc) 
// simple input - only number of steps variable is read
{
 //char line [OOFEM_MAX_LINE_LENGTH+1];
 Domain* domain;
 int i;
 bool inputReaderFinish=true;

 strcpy (this->dataOutputFileName, dataOutputFileName);

 if (this->giveProblemMode() ==   _postProcessor) {
  // modify output file name to prevent output to be lost
  strcat (dataOutputFileName,".oofeg");
 }
 if ((outputStream = fopen (dataOutputFileName,"w"))== NULL) {
  _error2 ("instanciateYourself: Can't open output file %s",dataOutputFileName);
 };
 
 fprintf (outputStream,"%s",PRG_HEADER);
 this->startTime = ::getTime ();
 //this->startClock= this-> getClock();
 fprintf (outputStream,"\nStarting analysis on: %s\n",ctime(&this->startTime));
 
 fprintf (outputStream,"%s\n",desc);
 
#  ifdef VERBOSE
 OOFEM_LOG_DEBUG ("Reading all data from input file \n") ;
#  endif
#ifdef __PARALLEL_MODE
 fprintf (outputStream,"Problem rank is %d/%d on %s\n\n",this->rank, this->numProcs, this->processor_name);
#endif

 // create domains
 domainNeqs.resize(this->ndomains);
 domainPrescribedNeqs.resize(this->ndomains);
 domainList -> growTo(ndomains) ;
 for (i=0; i < this->ndomains ; i++) {
  domain =  new Domain (i+1,0,this);
  domainList->put(i+1,domain) ;
 }

#ifdef __PETSC_MODULE
 this->initPetscContexts();
 /*
 petscContextList -> growTo(ndomains) ;
 for (i=0; i < this->ndomains ; i++) {
  petscContext =  new PetscContext (this);
  petscContextList->put(i+1,petscContext) ;
 }
 */
#endif

 // instanciate receiver 
 this->initializeFrom (ir);
 exportModuleManager->initializeFrom (ir);

 if (this->nMetaSteps == 0) {
   inputReaderFinish = false;
   this->instanciateDefaultMetaStep (ir);
 } else this->instanciateMetaSteps (dr);

 // instanciate export module manager
 exportModuleManager->instanciateYourself(dr, ir);
 this->instanciateDomains (dr);
 
 exportModuleManager->initialize();
 
 // check emodel input record if no default metastep, since all has been read
 if (inputReaderFinish) ir->finish();

 return 1;
}

IRResultType 
EngngModel::initializeFrom (InputRecord* ir)
{
 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

 IR_GIVE_FIELD (ir, numberOfSteps, IFT_EngngModel_nsteps, "nsteps"); // Macro
  if (numberOfSteps <= 0) _error ("instanciateFrom: nsteps not specified, bad format");

 contextOutputStep =  0;
 IR_GIVE_OPTIONAL_FIELD (ir, contextOutputStep, IFT_EngngModel_contextoutputstep, "contextoutputstep"); // Macro
 if (contextOutputStep) this->setUDContextOutputMode(contextOutputStep);

 renumberFlag = 0; IR_GIVE_OPTIONAL_FIELD (ir, renumberFlag, IFT_EngngModel_renumberflag, "renumber"); // Macro
 nMetaSteps   = 0; IR_GIVE_OPTIONAL_FIELD (ir, nMetaSteps, IFT_EngngModel_nmsteps, "nmsteps"); // Macro

#ifdef __PARALLEL_MODE
 IR_GIVE_OPTIONAL_FIELD (ir, parallelFlag, IFT_EngngModel_parallelflag, "parallelflag"); // Macro
// fprintf (stderr, "Parallel mode is %d\n", parallelFlag);
 
 /* Load ballancing support */
 int _val = 0;
 IR_GIVE_OPTIONAL_FIELD (ir, _val, IFT_NonLinearStatic_loadBallancingFlag, "lbflag"); // Macro
 loadBallancingFlag = _val;

 _val = 0;
 IR_GIVE_OPTIONAL_FIELD (ir, _val, IFT_NonLinearStatic_forceloadBallancingFlag, "forcelb1"); // Macro
 force_load_rebalance_in_first_step = _val;

#endif
  return IRRT_OK;
}

int
EngngModel::instanciateDomains (DataReader* dr) 
{
 int i, result = 1;
 // read problem domains
 for (i=0; i < this->ndomains ; i++) {
  result &= domainList->at(i+1)->instanciateYourself(dr);
 }
 return result;
}


int
EngngModel::instanciateMetaSteps (DataReader* dr) 
{
 int i, result = 1;
 MetaStep* mstep;

 // creat meta steps
 metaStepList -> growTo(nMetaSteps) ;
 for (i=0; i < this->nMetaSteps ; i++) {
  mstep =  new MetaStep (i+1,this);
  metaStepList->put(i+1,mstep) ;
 }

 // read problem domains
 for (i=0; i < this->nMetaSteps ; i++) {
  InputRecord* ir = dr->giveInputRecord (DataReader::IR_mstepRec, i+1);
  result &= metaStepList->at(i+1)->initializeFrom(ir);
 }

 // set meta step bounds
 int istep = this->giveNumberOfFirstStep();
 for (i=0; i < this->nMetaSteps ; i++) {
  istep = metaStepList->at(i+1)->setStepBounds (istep);
 }
 
 this->numberOfSteps = istep-1;
 OOFEM_LOG_RELEVANT ("Total number of solution steps %d\n", numberOfSteps);
 return result;
}


int 
EngngModel::instanciateDefaultMetaStep (InputRecord* ir)
{
 MetaStep* mstep;

 if (numberOfSteps == 0) _error ("instanciateDefaultMetaStep: nsteps cannot be zero");
 // create default meta steps
 this->nMetaSteps = 1;
 metaStepList -> growTo(nMetaSteps) ;
 mstep =  new MetaStep (1, this, numberOfSteps, *ir);
 metaStepList->put(1,mstep) ;

 // set meta step bounds
 int istep = this->giveNumberOfFirstStep()-1;
 metaStepList->at(1)->setStepBounds (istep+1);
 
 OOFEM_LOG_RELEVANT ("Total number of solution steps %d\n",  numberOfSteps);
 return 1;
}



/*
IntArray*  EngngModel :: GiveBanWidthVector () {
//
// Returns maximal column height for assembled characteristics matrix
// this method is implemented here, because some method may add some 
// conditions in to system and this may results into increased number of
// equations. 
//
  int neq;
  
  neq = this -> giveNumberOfEquations ();
  IntArray* mht = new IntArray (neq);
  IntArray* loc;
  int j,js,ieq,maxle ;

  for ( j =1 ; j<=neq; j++) mht->at(j)=INT_MAX ;          // initialize column height
  int nelem = domain -> giveNumberOfElements() ;
  for (int i = 1 ; i <= nelem ; i++ ) {
    loc = domain -> giveElement(i) -> giveLocationArray () ;
    js = loc -> giveSize() ;
    maxle = INT_MAX;
    for ( j = 1 ; j <= js ; j++ ) {
      ieq = loc->at(j);
      if(ieq != 0) maxle = min (maxle, ieq) ;
    }
    for ( j = 1 ; j <= js ; j++ ) {
      ieq = loc->at(j);
      if (ieq != 0) {
    mht->at(ieq) = min(maxle,mht->at(ieq));
      }
    }
  }
  return mht ;
}
*/

int  EngngModel :: giveNumberOfEquations (EquationID) {
//
// returns number of equations of current problem
// this method is implemented here, because some method may add some 
// conditions in to system and this may results into increased number of
// equations. 
//
  if (equationNumberingCompleted) return  numberOfEquations;

  return  this->forceEquationNumbering() ;
}

int  EngngModel :: giveNumberOfPrescribedEquations (EquationID) {
//
// returns number of equations of current problem
// this method is implemented here, because some method may add some 
// conditions in to system and this may results into increased number of
// equations. 
//
  if (equationNumberingCompleted) return  numberOfPrescribedEquations;

  return  this->forceEquationNumbering() ;
}

int  EngngModel :: giveNumberOfDomainEquations (int id, EquationID) {
//
// returns number of equations of current problem
// this method is implemented here, because some method may add some 
// conditions in to system and this may results into increased number of
// equations. 
//
  if (equationNumberingCompleted) return  domainNeqs.at(id);
 this->forceEquationNumbering();
  return  domainNeqs.at(id) ;
}

int  EngngModel :: giveNumberOfPrescribedDomainEquations (int id, EquationID) {
//
// returns number of equations of current problem
// this method is implemented here, because some method may add some 
// conditions in to system and this may results into increased number of
// equations. 
//
  if (equationNumberingCompleted) return  domainPrescribedNeqs.at(id);
 this->forceEquationNumbering();
  return  domainPrescribedNeqs.at(id) ;
}

int
EngngModel :: forceEquationNumbering (int id)
{
// forces equation renumbering for current time step
// intended mainly for problems with changes of static system
// during solution
// OUTPUT:
// sets this->numberOfEquations and this->numberOfPrescribedEquations and returns this value

 int i,j,ndofs,nnodes,nelem;
 DofManager *inode;
 Domain* domain = this->giveDomain(id);
 TimeStep* currStep = this->giveCurrentStep();
  IntArray loc;

 this->domainNeqs.at(id) = 0;
 this->domainPrescribedNeqs.at(id) = 0;
 

 nnodes = domain -> giveNumberOfDofManagers() ;
 
 if (!this->renumberFlag) {

  for ( i=1; i <= nnodes ; i++) {
   inode = domain -> giveDofManager (i) ;
   ndofs = inode->giveNumberOfDofs ();
   for (j=1;j<= ndofs; j++) inode->giveDof(j)->askNewEquationNumber(currStep);
  }

 } else {
  
  // invoke profile reduction 
  int initialProfile, optimalProfile;
  //clock_t time_0 = this->getClock(), time_1;
  oofem_timeval tstart;
  ::getUtime(tstart);
  OOFEM_LOG_INFO ("Renumbering ... ");
  
  SloanGraph graph (domain);
  graph.initialize();
  graph.tryParameters(0,0);
  initialProfile = graph.giveOptimalProfileSize();
  graph.tryParameters(2,1);
  graph.tryParameters(1,0);
  graph.tryParameters(5,1);
  graph.tryParameters(10,1);
  optimalProfile = graph.giveOptimalProfileSize();
  //FILE* renTableFile = fopen ("rentab.dat","w");
  //graph.writeOptimalRenumberingTable (renTableFile);
  IntArray* reverseRenTable = graph.giveOptimalRenumberingTable();
  //time_1 = this->getClock();
  //long nsec = (time_1 - time_0) / CLOCKS_PER_SEC;

  oofem_timeval ut;
  ::getRelativeUtime (ut, tstart);
  
  OOFEM_LOG_INFO ("done in %.2fs\n", (double)(ut.tv_sec+ut.tv_usec/(double)OOFEM_USEC_LIM));
  OOFEM_LOG_INFO ("Nominal profile %d (old) %d (new)\n", initialProfile, optimalProfile);
//  printf ("\nUser time consumed by renumbering: %lds", nsec);
  
  // undefine all dofs equation numbers
  for ( i=1; i <= nnodes ; i++) {
   inode = domain->giveDofManager (reverseRenTable->at(i));
   ndofs = inode->giveNumberOfDofs ();
   for (j=1;j<= ndofs; j++) inode->giveDof(j)->askNewEquationNumber(currStep);
  }
 }
 
 // invalidate element local copies of location arrays
 nelem = domain -> giveNumberOfElements() ;
 for ( i=1 ; i <= nelem ; i++ ){
  domain -> giveElement (i) -> invalidateLocationArray () ;
  }

/*
 // for equation numbering according to nodes numbers uncomment this section
 for ( i=1; i <= nnodes ; i++) {
 //  domain -> giveNode (i) -> giveCompleteLocationArray (loc) ;
 for equation numbering according to elements numbers uncomment this section
 nelem = domain -> giveNumberOfElements() ;
 for ( i=1 ; i <= nelem ; i++ ){
  loc = domain -> giveElement (i) -> giveLocationArray () ;
  js = loc.giveSize ();
  for (j=1 ; j <= js ; j++ ) 
  numberOfEquations = max (  numberOfEquations , loc.at(j));
  }
*/
 
  return  domainNeqs.at(id) ;
}


int
EngngModel :: forceEquationNumbering ()
{
 int i;
 // set numberOfEquations counter to zero
 this->numberOfEquations = 0;
 this->numberOfPrescribedEquations = 0;
 
 for (i = 1; i<= this->ndomains; i++) {
  domainNeqs.at(i) = 0;
  this->numberOfEquations+= this->forceEquationNumbering (i);
 }

 equationNumberingCompleted = 1;

 for (i = 1; i<= this->ndomains; i++) {
   //this->numberOfPrescribedEquations+=giveNumberOfPrescribedDomainEquations(i);
   this->numberOfPrescribedEquations+= domainPrescribedNeqs.at(i) ;
 }

#ifdef __PETSC_MODULE
 for (i = 1; i<= petscContextList->giveSize(); i++) 
	 this->petscContextList->at(i)->init(i);
#endif


 return this->numberOfEquations;
}



void 
EngngModel :: solveYourself ()
{
 int imstep, jstep;
 int smstep=1, sjstep=1;
 MetaStep* activeMStep;
 FILE* out = this->giveOutputStream ();
//#ifdef TIME_REPORT
 this->timer.startTimer(EngngModelTimer::EMTT_AnalysisTimer);
//#endif

 if (this->currentStep) {
  smstep = this->currentStep->giveMetaStepNumber();
  sjstep = this->giveMetaStep(smstep)->giveStepRelativeNumber(this->currentStep->giveNumber()) + 1;
 }


 for (imstep = smstep; imstep<= nMetaSteps; imstep++, sjstep = 1) {
   activeMStep = this->giveMetaStep(imstep);
   for (jstep = sjstep; jstep <= activeMStep->giveNumberOfSteps(); jstep++) {
//#ifdef TIME_REPORT
     this->timer.startTimer(EngngModelTimer::EMTT_SolutionStepTimer);
//#endif
     this->giveNextStep();
     // update state ccording to new meta step
     if (jstep == sjstep) this->initMetaStepAttributes (this->giveCurrentStep());
     // renumber equations if necessary
     if (this->requiresEquationRenumbering(this->giveCurrentStep())) this->forceEquationNumbering();
     this->solveYourselfAt(this->giveCurrentStep());
     this->terminate (this->giveCurrentStep());
     
//#ifdef TIME_REPORT
     this->timer.stopTimer(EngngModelTimer::EMTT_SolutionStepTimer);
     long int _steptime = this->timer.getUtime (EngngModelTimer::EMTT_SolutionStepTimer);
     OOFEM_LOG_INFO ("\nEngngModel info: user time consumed by solution step %d: %lds\n", 
                     this->giveCurrentStep()->giveNumber(), _steptime);
     
     fprintf (out,"\nUser time consumed by solution step %d: %ld [s]\n\n",
              this->giveCurrentStep()->giveNumber(), _steptime);
//#endif

#ifdef __PARALLEL_MODE
     if (loadBallancingFlag) this->ballanceLoad(this->giveCurrentStep());
#endif

   }
 }
}

void
EngngModel :: initMetaStepAttributes (TimeStep* tStep)
{
  MetaStep* mstep = this->giveMetaStep(tStep->giveMetaStepNumber());

  // update attributes
  this->updateAttributes (tStep); // virtual function
  // finish data acquiring
  mstep->giveAttributesRecord()->finish();
}

void
EngngModel :: updateAttributes (TimeStep* atTime)
{
 MetaStep* mstep = this->giveMetaStep(atTime->giveMetaStepNumber());
 InputRecord* ir = mstep->giveAttributesRecord();
 
 if(this->giveNumericalMethod (atTime)) 
  this->giveNumericalMethod (atTime) -> initializeFrom (ir);
}


void
EngngModel :: updateYourself (TimeStep* stepN)
{
 int idomain, ndomains = this->giveNumberOfDomains();
  int j, nnodes;
 Domain* domain;
 
 for (idomain = 1; idomain <= ndomains; idomain++) {

  domain= this->giveDomain(idomain);
  
#  ifdef VERBOSE
  VERBOSE_PRINT0("Updating domain ",domain->giveNumber())
#  endif
   
  nnodes = domain->giveNumberOfDofManagers ();
  for( j=1;j<=nnodes;j++) {
   domain->giveDofManager(j)->updateYourself(stepN) ;
   //domain->giveDofManager(j)->printOutputAt(File, stepN);
  }
  
#  ifdef VERBOSE
  VERBOSE_PRINT0("Updated nodes & sides ",nnodes)
#  endif
   
   
   Element* elem;
  
  int nelem = domain->giveNumberOfElements ();
  for (j=1 ; j<=nelem ; j++) {
   elem = domain -> giveElement(j) ;
#ifdef __PARALLEL_MODE
   // skip remote elements (these are used as mirrors of remote eleemnts on other domains
   // when nonlocal constitutive models are used. They introduction is necessary to
   // allow local averaging on domains without fine grain communication between domains).
   if (elem->giveParallelMode () == Element_remote) continue;
#endif
   elem -> updateYourself(stepN);
   //elem -> printOutputAt(File, stepN) ;
  }
  
#  ifdef VERBOSE
  VERBOSE_PRINT0("Updated Elements ",nelem)
#  endif
   
 }
}

void 
EngngModel :: terminate (TimeStep* stepN)
{
 this->doStepOutput(stepN);
 this->saveStepContext(stepN);
}


void
EngngModel :: doStepOutput (TimeStep* stepN)
{
 FILE* File = this->giveOutputStream();
 
 // print output
 this->printOutputAt (File, stepN);
 // export using export manager
 exportModuleManager->doOutput(stepN);
}

void
EngngModel :: saveStepContext (TimeStep* stepN)
{

 // save context if required
 // default - save only if ALWAYS is set ( see cltypes.h )
 
 if ((this->giveContextOutputMode() == ALWAYS) ||
     (this->giveContextOutputMode() == REQUIRED)) {
   this->saveContext(NULL, CM_State);
 }
 else if (this->giveContextOutputMode() == USERDEFINED) {
  if (stepN->giveNumber()%this->giveContextOutputStep() == 0) 
    this->saveContext(NULL, CM_State);
 }
}


void
EngngModel :: printOutputAt (FILE * File,TimeStep* stepN) 
{
  //FILE* File = this -> giveDomain() -> giveOutputStream() ;
 int domCount = 0, idomain;
 Domain* domain;

  // fprintf (File,"\nOutput for time step number %d \n\n",stepN->giveNumber());
 for (idomain = 1; idomain <= this->ndomains; idomain++) {
  domain= this->giveDomain(idomain);
  domCount += domain->giveOutputManager()->testTimeStepOutput (stepN);
 }
 if (domCount == 0) return;  // do not print even Solution step header
 
 fprintf (File,"\n==============================================================");
 fprintf (File,"\nOutput for time % .8e ",stepN->giveTime()*this->giveVariableScale(VST_Time));
 fprintf (File,"\n==============================================================\n");
 for (idomain = 1; idomain <= this->ndomains; idomain++) {
  
  domain= this->giveDomain(idomain);
  fprintf (File,"Output for domain %3d\n",domain->giveNumber());
  
  domain->giveOutputManager()->doDofManOutput  (File, stepN);
  domain->giveOutputManager()->doElementOutput (File, stepN);
 }
}

void EngngModel :: printYourself ()
{
  printf ("\nEngineeringModel: instance %s\n",this->giveClassName());
  printf ("number of steps: %d\n",this->giveNumberOfSteps());
  printf ("number of eq's : %d\n",numberOfEquations);
}

/*
void EngngModel :: assemble (SparseMtrx* answer, TimeStep* tStep, CharType type) 
//
// assembles matrix answer by  calling  
// element(i) -> giveCharacteristicMatrix ( type, tStep );
// for each element in domain 
// and assembling every contribution to answer
//
//
{
  Element* element;
  IntArray* loc ;

  if (answer == NULL)  _error("assemble: NULL pointer encountered.");

  int nelem = domain -> giveNumberOfElements ();
  FloatMatrix* charMtrx;
  for (int i = 1; i <= nelem ; i++ ) {
    element = domain -> giveElement(i);
    loc = element -> giveLocationArray ();
    charMtrx = element -> GiveCharacteristicMatrix ( type, tStep );
    if (charMtrx) answer ->  assemble (charMtrx, loc) ;
  delete charMtrx;
  }
}  
*/


void EngngModel :: assemble (SparseMtrx* answer, TimeStep* tStep, EquationID ut, CharType type, Domain* domain) 
//
// assembles matrix answer by  calling  
// element(i) -> giveCharacteristicMatrix ( type, tStep );
// for each element in domain 
// and assembling every contribution to answer
//
//
{
 int ielem;
  IntArray loc ;
  FloatMatrix mat;
  Element *element;

  if (answer == NULL)  _error("assemble: NULL pointer encountered.");

  int nelem = domain -> giveNumberOfElements ();
  for ( ielem = 1; ielem <= nelem ; ielem++ ) {
    element = domain -> giveElement(ielem);
#ifdef __PARALLEL_MODE
    // skip remote elements (these are used as mirrors of remote eleemnts on other domains
    // when nonlocal constitutive models are used. They introduction is necessary to
    // allow local averaging on domains without fine grain communication between domains).
    if (element->giveParallelMode () == Element_remote) continue;
#endif
    element -> giveLocationArray (loc, ut);
    this->giveElementCharacteristicMatrix(mat, ielem, type, tStep, domain );
 
    if (mat.isNotEmpty()) {
      if (answer -> assemble (loc, mat) == 0) 
	_error("assemble: sparse matrix assemble error");
    }  
  }
  answer->assembleBegin();
  answer->assembleEnd();
}

void EngngModel :: assemble (SparseMtrx* answer, TimeStep* tStep, EquationID r_id, EquationID c_id, CharType type, Domain* domain) 
//
// assembles matrix answer by  calling  
// element(i) -> giveCharacteristicMatrix ( type, tStep );
// for each element in domain 
// and assembling every contribution to answer
//
//
{
 int ielem;
  IntArray r_loc, c_loc ;
  FloatMatrix mat;
  Element *element;

  if (answer == NULL)  _error("assemble: NULL pointer encountered.");

  int nelem = domain -> giveNumberOfElements ();
  for ( ielem = 1; ielem <= nelem ; ielem++ ) {
    element = domain -> giveElement(ielem);
#ifdef __PARALLEL_MODE
    // skip remote elements (these are used as mirrors of remote eleemnts on other domains
    // when nonlocal constitutive models are used. They introduction is necessary to
    // allow local averaging on domains without fine grain communication between domains).
    if (element->giveParallelMode () == Element_remote) continue;
#endif
    element -> giveLocationArray (r_loc, r_id);
    element -> giveLocationArray (c_loc, c_id);
    
    this->giveElementCharacteristicMatrix(mat, ielem, type, tStep, domain );
 
    if (mat.isNotEmpty()) {
      if (answer -> assemble (r_loc, c_loc, mat) == 0) 
	_error("assemble: sparse matrix assemble error");
    }
    answer->assembleBegin();
    answer->assembleEnd();
  }
  //answer->assembleBegin();
  //answer->assembleEnd();
}


/*
void EngngModel :: assemble (FloatArray& answer, TimeStep* tStep, CharType type, Domain* domain) 
//
// assembles matrix answer by  calling  
// element(i) -> giveCharacteristicMatrix ( type, tStep );
// or node(i) -> computeLoadVectorAt (tStep); 
// for each element in domain 
// and assembling every contribution to answer
//
//
{
  int i ;
  IntArray loc ;
  FloatArray charVec ;
 ValueModeType mode;
 Element *element ;
 DofManager *node ;
 
  int nnode = domain -> giveNumberOfDofManagers();
  int nelem = domain -> giveNumberOfElements ();

  switch (type) {
  case NodalLoadVector_Total:
  case NodalLoadVector_Incremental:
  
  if (type == NodalLoadVector_Total) mode = TotalMode; else mode = IncrementalMode;
    for (i = 1; i <= nnode ; i++ ) {
      node = domain -> giveDofManager(i);
      node -> giveCompleteLocationArray (loc);
      node -> computeLoadVectorAt (charVec, tStep, mode);
      if(charVec.giveSize()) answer.assemble (charVec, loc) ;
    }
  
    break ;

  case ElementForceLoadVector_Total:
  case ElementForceLoadVector_Incremental:
  case ElementNonForceLoadVector_Total:
  case ElementNonForceLoadVector_Incremental:
  case NodalInternalForcesVector_Total:
  case NodalInternalForcesVector_Incremental:
 case ElementHEMOLoadVector_Total:

    for (i = 1; i <= nelem ; i++ ) {
      element = domain -> giveElement(i);
#ifdef __PARALLEL_MODE
   // skip remote elements (these are used as mirrors of remote eleemnts on other domains
   // when nonlocal constitutive models are used. They introduction is necessary to
   // allow local averaging on domains without fine grain communication between domains).
   if (element->giveParallelMode () == Element_remote) continue;
#endif
      element -> giveLocationArray (loc);
      this -> giveElementCharacteristicVector (charVec, i, type, tStep, domain);
      if(charVec.giveSize()) answer.assemble (charVec, loc) ;
    }  
    break ;

 case ElementForceLoadVectorOfPrescribed_Total:
  case ElementForceLoadVectorOfPrescribed_Incremental:

  CharType mtype;
  if (type == ElementForceLoadVectorOfPrescribed_Total) mtype = ElementForceLoadVector_Total;
  else mtype = ElementForceLoadVector_Incremental;

  for (i=1; i<= nelem; i++) {
      element = domain -> giveElement(i);
#ifdef __PARALLEL_MODE
   // skip remote elements (these are used as mirrors of remote eleemnts on other domains
   // when nonlocal constitutive models are used. They introduction is necessary to
   // allow local averaging on domains without fine grain communication between domains).
   if (element->giveParallelMode () == Element_remote) continue;
#endif
      element -> givePrescribedLocationArray (loc);
      this -> giveElementCharacteristicVector (charVec, i, mtype, tStep, domain);
      if(charVec.giveSize()) answer.assemble (charVec, loc) ;
  }
  break ;

 case NodalLoadVectorOfPrescribed_Total:
  case NodalLoadVectorOfPrescribed_Incremental:
  
  if (type == NodalLoadVector_Total) mode = TotalMode; else mode = IncrementalMode;
    for (i = 1; i <= nnode ; i++ ) {
      node = domain -> giveDofManager(i);
      node -> giveCompletePrescribedLocationArray (loc);
      node -> computeLoadVectorAt (charVec, tStep, mode);
      if(charVec.giveSize()) answer.assemble (charVec, loc) ;
    }
  
    break ;

  default:
     _error("assemble: Unknown Type of characteristic mtrx.");
  }
}
*/

void EngngModel :: assembleVectorFromDofManagers (FloatArray& answer, TimeStep* tStep, EquationID ut, 
                                                  CharType type, ValueModeType mode, Domain* domain) 
//
// assembles matrix answer by  calling  
// node(i) -> computeLoadVectorAt (tStep); 
// for each element in domain 
// and assembling every contribution to answer
//
//
{
  int i ;
  IntArray loc ;
  FloatArray charVec ;
  DofManager *node ;
  
  int nnode = domain -> giveNumberOfDofManagers();
  
  for (i = 1; i <= nnode ; i++ ) {
    node = domain -> giveDofManager(i);
    node -> computeLoadVectorAt (charVec, tStep, mode);
    if(charVec.giveSize()) {
      node -> giveCompleteLocationArray (loc);
      answer.assemble (charVec, loc) ;
    }
  }
}


void EngngModel :: assemblePrescribedVectorFromDofManagers (FloatArray& answer, TimeStep* tStep, EquationID ut,
                                                            CharType type, ValueModeType mode, Domain* domain) 
//
// assembles matrix answer by  calling  
// node(i) -> computeLoadVectorAt (tStep); 
// for each element in domain 
// and assembling every contribution to answer
//
//
{
  int i ;
  IntArray loc ;
  FloatArray charVec ;
  DofManager *node ;
  
  int nnode = domain -> giveNumberOfDofManagers();
  
  for (i = 1; i <= nnode ; i++ ) {
    node = domain -> giveDofManager(i);
    node -> computeLoadVectorAt (charVec, tStep, mode);
    if(charVec.giveSize()) {
      node -> giveCompletePrescribedLocationArray (loc);
      answer.assemble (charVec, loc) ;
    }
  }
}

void EngngModel :: assembleVectorFromElements (FloatArray& answer, TimeStep* tStep, EquationID ut,
                                               CharType type, ValueModeType mode, Domain* domain) 
//
// assembles matrix answer by  calling  
// element(i) -> giveCharacteristicMatrix ( type, tStep );
// for each element in domain 
// and assembling every contribution to answer
//
//
{
  int i ;
  IntArray loc ;
  FloatArray charVec ;
  Element *element ;
 
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
  this -> giveElementCharacteristicVector (charVec, i, type, mode, tStep, domain);
  if(charVec.giveSize()) answer.assemble (charVec, loc) ;
 }  
}

void EngngModel :: assemblePrescribedVectorFromElements (FloatArray& answer, TimeStep* tStep, EquationID ut, 
                                                         CharType type, ValueModeType mode, Domain* domain) 
//
// assembles matrix answer by  calling  
// element(i) -> giveCharacteristicMatrix ( type, tStep );
// for each element in domain 
// and assembling every contribution to answer
//
//
{
  int i ;
  IntArray loc ;
  FloatArray charVec ;
 Element *element ;
 
  int nelem = domain -> giveNumberOfElements ();

 for (i = 1; i <= nelem ; i++ ) {
  element = domain -> giveElement(i);
#ifdef __PARALLEL_MODE
  // skip remote elements (these are used as mirrors of remote eleemnts on other domains
  // when nonlocal constitutive models are used. They introduction is necessary to
  // allow local averaging on domains without fine grain communication between domains).
  if (element->giveParallelMode () == Element_remote) continue;
#endif
  element -> givePrescribedLocationArray (loc, ut);
  if (loc.containsOnlyZeroes()) continue;
  this -> giveElementCharacteristicVector (charVec, i, type, mode, tStep, domain);
  if(charVec.giveSize()) answer.assemble (charVec, loc) ;
 }  
}




#ifdef __PETSC_MODULE
void 
EngngModel :: petsc_assembleVectorFromDofManagers (Vec answer, TimeStep* tStep, EquationID ut,
                                                   CharType type, ValueModeType mode, Domain* domain) 
//
// assembles matrix answer by  calling  
// node(i) -> computeLoadVectorAt (tStep); 
// for each element in domain 
// and assembling every contribution to answer
//
//
{
  int i, ni ;
  IntArray loc;
#ifdef __PARALLEL_MODE
  IntArray gloc ;
#endif
  FloatArray charVec ;
  DofManager *node ;
 
  int nnode = domain -> giveNumberOfDofManagers();
  
  for (i = 1; i <= nnode ; i++ ) {
    node = domain -> giveDofManager(i);
    node -> computeLoadVectorAt (charVec, tStep, mode);
    if((ni = charVec.giveSize())) {
      node -> giveCompleteLocationArray (loc);
#ifdef __PARALLEL_MODE
      this->givePetscContext(domain->giveNumber(),ut)->giveN2Gmap()->map2New(gloc, loc, 0);
      VecSetValues(answer, ni, gloc.givePointer(), charVec.givePointer(), ADD_VALUES);
#else
      loc.add(-1);
      VecSetValues(answer, ni, loc.givePointer(), charVec.givePointer(), ADD_VALUES);
#endif
    }
  }
}

void 
EngngModel :: petsc_assemblePrescribedVectorFromDofManagers (Vec answer, TimeStep* tStep, EquationID ut, 
                                                             CharType type, ValueModeType mode, Domain* domain) 
//
// assembles matrix answer by  calling  
// node(i) -> computeLoadVectorAt (tStep); 
// for each element in domain 
// and assembling every contribution to answer
//
//
{
  int i, ni ;
  IntArray loc ;
#ifdef __PARALLEL_MODE
  IntArray gloc ;
#endif
  FloatArray charVec ;
  DofManager *node ;
  
  int nnode = domain -> giveNumberOfDofManagers();
  
  for (i = 1; i <= nnode ; i++ ) {
    node = domain -> giveDofManager(i);
    node -> computeLoadVectorAt (charVec, tStep, mode);
    if((ni=charVec.giveSize())) {
      node -> giveCompletePrescribedLocationArray (loc);
#ifdef __PARALLEL_MODE
      this->givePetscContext(domain->giveNumber(),ut)->giveN2Gmap()->map2New(gloc, loc, 0); // ????
      VecSetValues(answer, ni, gloc.givePointer(), charVec.givePointer(), ADD_VALUES);
#else
      loc.add(-1);
      VecSetValues(answer, ni, loc.givePointer(), charVec.givePointer(), ADD_VALUES);
#endif
    }
  }
}

void 
EngngModel :: petsc_assembleVectorFromElements (Vec answer, TimeStep* tStep, EquationID ut,
                                                CharType type, ValueModeType mode, Domain* domain) 
//
// assembles matrix answer by  calling  
// element(i) -> giveCharacteristicMatrix ( type, tStep );
// for each element in domain 
// and assembling every contribution to answer
//
//
{
  int i,ni ;
  IntArray loc ;
#ifdef __PARALLEL_MODE
  IntArray gloc ;
#endif
  FloatArray charVec ;
  Element *element ;
 
  int nelem = domain -> giveNumberOfElements ();
  
  for (i = 1; i <= nelem ; i++ ) {
    element = domain -> giveElement(i);
#ifdef __PARALLEL_MODE
    // skip remote elements (these are used as mirrors of remote eleemnts on other domains
    // when nonlocal constitutive models are used. They introduction is necessary to
    // allow local averaging on domains without fine grain communication between domains).
    if (element->giveParallelMode () == Element_remote) continue;
#endif
    this -> giveElementCharacteristicVector (charVec, i, type, mode, tStep, domain);
    if((ni=charVec.giveSize())) {
      element -> giveLocationArray (loc, ut);
#ifdef __PARALLEL_MODE
      this -> givePetscContext(domain->giveNumber(),ut)->giveN2Gmap()->map2New(gloc, loc, 0);
      VecSetValues(answer, ni, gloc.givePointer(), charVec.givePointer(), ADD_VALUES);
#else
      loc.add(-1);
      VecSetValues(answer, ni, loc.givePointer(), charVec.givePointer(), ADD_VALUES);
#endif

    }  
  }
}

void 
EngngModel :: petsc_assemblePrescribedVectorFromElements (Vec answer, TimeStep* tStep, EquationID ut,
							  CharType type, ValueModeType mode, Domain* domain) 
//
// assembles matrix answer by  calling  
// element(i) -> giveCharacteristicMatrix ( type, tStep );
// for each element in domain 
// and assembling every contribution to answer
//
//
{
  int i,ni ;
  IntArray loc ;
#ifdef __PARALLEL_MODE
  IntArray gloc ;
#endif
  FloatArray charVec ;
  Element *element ;
  
  int nelem = domain -> giveNumberOfElements ();
  
  for (i = 1; i <= nelem ; i++ ) {
    element = domain -> giveElement(i);
#ifdef __PARALLEL_MODE
    // skip remote elements (these are used as mirrors of remote eleemnts on other domains
    // when nonlocal constitutive models are used. They introduction is necessary to
    // allow local averaging on domains without fine grain communication between domains).
    if (element->giveParallelMode () == Element_remote) continue;
#endif
    this -> giveElementCharacteristicVector (charVec, i, type, mode, tStep, domain);
    if((ni=charVec.giveSize())) {
      element -> givePrescribedLocationArray (loc, ut);
#ifdef __PARALLEL_MODE
      this -> givePetscContext(domain->giveNumber(),ut)->giveN2Gmap()->map2New(gloc, loc, 0);  // ??
      VecSetValues(answer, ni, gloc.givePointer(), charVec.givePointer(), ADD_VALUES);
#else
      loc.add(-1);
      VecSetValues(answer, ni, loc.givePointer(), charVec.givePointer(), ADD_VALUES);
#endif
    }
  }  
}

#endif

void 
EngngModel ::  updateComponent (TimeStep* tStep, NumericalCmpn cmpn, Domain* d)
//
// updates some componet, which is used by numerical method
// to newly reached state
//
{
 _error("updateComponent: Unknown Type of component.");
}


void
EngngModel :: initStepIncrements ()
//
// resets all temp data in elements and their gps so only
// non temp variables remain.
// this function can be used if called before this->terminate()
// for current step restart, because it causes complete reset
// of temp variables.
//
{
 int j;
 Element *elem;
 Domain* domain;

 for (int idomain = 1; idomain <= this->ndomains; idomain++) {
  domain= this->giveDomain(idomain);

  int nelem = domain->giveNumberOfElements ();
  for (j=1 ; j<=nelem ; j++) {
   elem = domain -> giveElement(j) ;
#ifdef __PARALLEL_MODE
   // skip remote elements (these are used as mirrors of remote eleemnts on other domains
   // when nonlocal constitutive models are used. They introduction is necessary to
   // allow local averaging on domains without fine grain communication between domains).
   if (elem->giveParallelMode () == Element_remote) continue;
#endif
   elem -> initForNewStep();
  }
 }
 return;
}





contextIOResultType EngngModel :: saveContext(DataStream* stream, ContextMode mode, void *obj )
//
// this procedure is used mainly for two reasons:
//
// - to save context of current model in order to be able to backtrace 
//   computational history (useful in nonlinear analysis, when you want
//   explore another, previously detected load path
//
// - to save context of current model in order to enable
//   possibility of post-processing.
//
// saving context means: (if needed may be enhanced)
//
// - save EngngModel state varialbles as displacement, velocity, .. vectors
// - save Elements stress, strain and material history.
//
// This version saves only Element and Material properties.
//
{
 contextIOResultType iores;
 int i, serNum, closeFlag = 0;
 Element *element;
 Domain* domain;
 ErrorEstimator* ee;
 FILE* file;
 

  OOFEM_LOG_INFO ("Storing context\n");
  if (stream==NULL) {
    if (!this->giveContextFile(&file, this->giveCurrentStep()->giveNumber(), 
                               this->giveCurrentStep()->giveVersion(), contextMode_write)) 
      THROW_CIOERR(CIO_IOERR); // override 
    stream = new FileDataStream (file);
    closeFlag = 1;
 }
 
 // store solution step
 if ((iores = giveCurrentStep ()->saveContext (stream, mode)) != CIO_OK) THROW_CIOERR(iores);
 
 // store numberOfEquations and domainNeqs array
 if (!stream->write (&numberOfEquations,1)) THROW_CIOERR(CIO_IOERR);
 if ((iores = domainNeqs.storeYourself (stream, mode)) != CIO_OK) THROW_CIOERR(iores);
 // store numberOfPrescribedEquations and domainNeqs array
 if (!stream->write (&numberOfPrescribedEquations,1)) THROW_CIOERR(CIO_IOERR);
 if ((iores = domainPrescribedNeqs.storeYourself (stream,mode)) != CIO_OK) THROW_CIOERR(iores);
 // store renumber flag
 if (!stream->write (&renumberFlag,1)) THROW_CIOERR(CIO_IOERR);

 
 for (int idomain = 1; idomain <= this->ndomains; idomain++) {
   domain= this->giveDomain(idomain);
   
   // save domain serial number
   serNum = domain->giveSerialNumber();
   if (!stream->write (&serNum, 1)) THROW_CIOERR(CIO_IOERR);
   
   // nodes & sides and corresponding dofs
   int nnodes = domain->giveNumberOfDofManagers ();
   for (i = 1; i <= nnodes ; i++ ) {
     if ((iores = domain -> giveDofManager(i)->saveContext(stream, mode)) != CIO_OK) THROW_CIOERR(iores);
   }
   
   // elements and corresponding integration points
   int nelem = domain -> giveNumberOfElements ();
   for ( i = 1; i <= nelem ; i++ ) {
     element = domain -> giveElement(i);
#ifdef __PARALLEL_MODE
     // skip remote elements (these are used as mirrors of remote eleemnts on other domains
     // when nonlocal constitutive models are used. They introduction is necessary to
     // allow local averaging on domains without fine grain communication between domains).
     if (element->giveParallelMode () == Element_remote) continue;
#endif
     if ((iores = element->saveContext(stream, mode)) != CIO_OK) THROW_CIOERR(iores);
   }
   
   // store error estimator data
   ee = this->giveDomainErrorEstimator(idomain);
   if (ee) 
     if ((iores = ee->saveContext(stream, mode)) != CIO_OK) THROW_CIOERR(iores);
   
 }


 // store nMethod
 NumericalMethod* nmethod = this->giveNumericalMethod(giveCurrentStep());
 if (nmethod) 
  if ((iores = nmethod->saveContext(stream, mode)) != CIO_OK) THROW_CIOERR(iores);


 if (closeFlag) {fclose(file); delete (stream); stream=NULL;} // ensure consistent records
 return CIO_OK;
}  
  

contextIOResultType EngngModel :: restoreContext(DataStream* stream, ContextMode mode, void *obj)
//
// this procedure is used mainly for two reasons:
//
// - to restore context of current model in order to be able to backtrace 
//   computational history (useful in nonlinear analysis, when you want
//   explore another, previously detected load path
//
// - to restore context of current model in order to enable
//   possibility of post-processing.
//
// restoring context means: (if needed may be enhanced)
//
// - rest. EngngModel state varialbles as displacement, velocity, .. vectors
// - rest. Elements stress, strain and material history.
//
// This version loads only Element and Material properties.
//
// This function is inverse to the saveContext() member function
// 
// WARNING obj is cast into int pointer  to time step for which to seek Context.
//
{
  contextIOResultType iores;
  int i, serNum, closeFlag = 0, istep, iversion;
  bool domainUpdated = false;
  Element* element;
  Domain* domain;
  ErrorEstimator* ee;
  FILE* file;

  this -> resolveCorrespondingStepNumber (istep, iversion, obj);
  OOFEM_LOG_RELEVANT ("Restoring context for tStep %d.%d\n",istep,iversion);

 if (stream == NULL) {
   if (!this->giveContextFile(&file, istep, iversion, contextMode_read)) THROW_CIOERR(CIO_IOERR); // override 
   stream = new FileDataStream (file);
   closeFlag = 1;
 }

 // restore solution step
 if (currentStep == NULL) currentStep = new TimeStep (istep,this,0, 0.,0.,0) ;
 if ((iores = currentStep -> restoreContext(stream, mode)) != CIO_OK) THROW_CIOERR(iores);
 // this->updateAttributes (currentStep);

 int pmstep = currentStep->giveMetaStepNumber();
 if (nMetaSteps) {
  if (!this->giveMetaStep(pmstep)->isStepValid(istep-1)) pmstep--;
 }

 if (previousStep) delete previousStep;
 previousStep = new TimeStep (istep-1, this, pmstep, currentStep->giveTime()-currentStep->giveTimeIncrement(),
                currentStep->giveTimeIncrement(), currentStep->giveSolutionStateCounter () - 1);

 // restore numberOfEquations and domainNeqs array
 if (!stream->read (&numberOfEquations,1)) THROW_CIOERR(CIO_IOERR);
 if ((iores = domainNeqs.restoreYourself (stream, mode)) != CIO_OK) THROW_CIOERR(iores);
 // restore numberOfPrescribedEquations and domainNeqs array
 if (!stream->read (&numberOfPrescribedEquations,1)) THROW_CIOERR(CIO_IOERR);
 if ((iores = domainPrescribedNeqs.restoreYourself (stream, mode)) != CIO_OK) THROW_CIOERR(iores);
 // restore renumber flag
 if (!stream->read (&renumberFlag,1)) THROW_CIOERR(CIO_IOERR);

  for (int idomain = 1; idomain <= this->ndomains; idomain++) {
    domain= this->giveDomain(idomain);
    
    domainUpdated = false;
    // restore domain serial number
    if (!stream->read (&serNum, 1)) THROW_CIOERR(CIO_IOERR);
    if (serNum != domain->giveSerialNumber()) {
      
      OOFEM_LOG_INFO ("deleting old domain\n");
      domainList->remove(idomain);
      
      // read corresponding domain
      OOFEM_LOG_INFO ("restoring domain %d.%d\n", idomain, serNum);
      Domain* dNew = new Domain (idomain, serNum, this);
      DataReader* domainDr = this->GiveDomainDataReader (1, serNum, contextMode_read);
      if (!dNew -> instanciateYourself(domainDr)) _error ("initializeAdaptive: domain Instanciation failed");
      delete domainDr;
      
      domain = dNew;
      domainList->put(idomain, dNew);
      domainUpdated = true;
    }
    
    //this->forceEquationNumbering (); - Equation numbers are stored too in dof contexts
    
    // nodes & sides and corresponding dofs
    int nnodes = domain->giveNumberOfDofManagers ();
    for (i = 1; i <= nnodes ; i++ ) {
      if ((iores = domain -> giveDofManager(i)->restoreContext(stream,mode)) != CIO_OK) THROW_CIOERR(iores);
    }
    
    int nelem = domain -> giveNumberOfElements ();
    for ( i = 1; i <= nelem ; i++ ) {
      element = domain -> giveElement(i);
#ifdef __PARALLEL_MODE
      // skip remote elements (these are used as mirrors of remote eleemnts on other domains
      // when nonlocal constitutive models are used. They introduction is necessary to
      // allow local averaging on domains without fine grain communication between domains).
      if (element->giveParallelMode () == Element_remote) continue;
#endif
      if ((iores = element->restoreContext(stream,mode)) != CIO_OK) THROW_CIOERR(iores) ;
    }
    
    // restore error estimator data
    ee = this->giveDomainErrorEstimator(idomain);
    if (domainUpdated) ee->setDomain (domain);
    if (ee) 
      if ((iores = ee->restoreContext(stream,mode)) != CIO_OK) THROW_CIOERR(iores);
    
  }
 
 // restore nMethod
 NumericalMethod* nmethod = this->giveNumericalMethod(giveCurrentStep());
 if (nmethod) 
 if ((iores = nmethod->restoreContext(stream,mode)) != CIO_OK) THROW_CIOERR(iores);

 this->updateDomainLinks();
 this->updateAttributes (currentStep);
 this->initStepIncrements();

  if (closeFlag) {fclose (file); delete stream; stream = NULL;} // ensure consistent records
  return CIO_OK;
}  


void
EngngModel :: resolveCorrespondingStepNumber (int& istep, int& iversion, void* obj)
{
// 
// returns corresponding step number
//
 if (obj == NULL) {
  istep = 1;
  iversion = 0;
  return;
 }
 istep = *(int*) obj;
 iversion = *(((int*) obj)+1);

 if (istep > this->giveNumberOfSteps()) {
  istep = this->giveNumberOfSteps();
 }
 if (istep <= 0) istep = 1;
 return;
}

 
int
EngngModel :: giveContextFile(FILE** contextFile, int stepNumber, int stepVersion, ContextFileMode cmode, int errLevel)
//
//
// assigns context file of given step number to stream
// returns nonzero on success
//
{

  char fname[MAX_FILENAME_LENGTH];
  sprintf (fname, "%s.%d.%d.osf", this->dataOutputFileName, stepNumber, stepVersion);
 
 if (cmode ==  contextMode_read) { 
  *contextFile = fopen(fname,"rb"); // open for reading
 } else {
  *contextFile = fopen(fname,"wb"); // open for writting, 
 }  
//  rewind (*contextFile); // seek the beginning
// // overwrite if exist
// else *contextFile = fopen(fname,"r+"); // open for reading and writting
 
 if((*contextFile == NULL) && errLevel > 0) {
   _error2 ("giveContextFile : can't open %s",fname);
 }

  return 1 ;
}

bool
EngngModel :: testContextFile (int stepNumber, int stepVersion)
{
  char fname[MAX_FILENAME_LENGTH];
  sprintf (fname, "%s.%d.%d.osf", this->dataOutputFileName, stepNumber, stepVersion);
#ifdef HAVE_ACCESS
 if (access (fname, R_OK) == 0) return true;
 else return false;
#elif _MSC_VER
	if (_access(fname, 4) == 0) return true;
  else return false;
#else
  return true;
#endif

}
 
DataReader*
EngngModel :: GiveDomainDataReader(int domainNum, int domainSerNum, ContextFileMode cmode)
//
//
// returns domain i/o file 
// returns nonzero on success
//
{
 
  char fname[MAX_FILENAME_LENGTH];
  sprintf (fname, "%s.domain.%d.%d.din", this->dataOutputFileName, domainNum, domainSerNum);
 DataReader* dr;
 
 if ((dr = new OOFEMTXTDataReader(fname)) == NULL)
  _error ("Creation of DataReader failed");
 
  return dr;
}

void 
EngngModel :: error (char* file, int line, char *format, ...) const
{
  char buffer[MAX_ERROR_MSG_LENGTH];
	va_list args;

	va_start(args, format);
	vsprintf(buffer, format, args);
	va_end(args);
#ifdef __PARALLEL_MODE
  __OOFEM_ERROR4 (file, line, "Class: %s, Rank: %d\n%s",giveClassName(),rank,buffer);
#else
  __OOFEM_ERROR3 (file, line, "Class: %s\n%s",giveClassName(),buffer);
#endif
}


void  EngngModel :: warning (char* file, int line, char *format, ...) const
//
// this function handles error reporting
// prints errorMsg enriched by ClasName and Number
// to the standart error stream
{
  char buffer[MAX_ERROR_MSG_LENGTH];
	va_list args;

	va_start(args, format);
	vsprintf(buffer, format, args);
	va_end(args);
#ifdef __PARALLEL_MODE
  __OOFEM_WARNING4 (file, line, "Class: %s, Rank: %d\n%s",giveClassName(),rank,buffer);
#else
  __OOFEM_WARNING3 (file, line, "Class: %s\n%s",giveClassName(),buffer);
#endif
}

Domain* 
EngngModel::giveDomain (int i)
{
 if ((i > 0) && (i <= this->ndomains)) {
  return this->domainList->at(i);
 } else {
  _error ("giveDomain: Undefined domain");
 }
 return NULL;
}


#ifdef __PETSC_MODULE
PetscContext* 
EngngModel::givePetscContext (int i, EquationID ut)
{
 if ((i > 0) && (i <= this->ndomains)) {
   if  (i> petscContextList->giveSize()) {
     _error ("givePetscContext: petsc context not initialized for this problem");
   }
  return this->petscContextList->at(i);
 } else {
  _error ("givePetscContext: Undefined domain or ");
 }
 return NULL;
}

void
EngngModel::initPetscContexts () {
}
#endif   



MetaStep* 
EngngModel::giveMetaStep (int i)
{
 if ((i > 0) && (i <= this->nMetaSteps)) {
  return this->metaStepList->at(i);
 } else {
   _error2 ("giveMetaStep: undefined metaStep (%d)", i);
 }
 return NULL;
}
/*
char*  
EngngModel :: giveInputDataFileName ()

{
 // Returns the name of the file containing the data of the problem.
 char s[MAX_FILENAME_LENGTH] ;
 
 if (! dataInputFileName) {
  printf ("please enter the name of the input data file : \n") ;
  gets (s) ;
#ifndef __PARALLEL_MODE
  dataInputFileName = new char[strlen(s)+1] ;
  strcpy (dataInputFileName,s) ;
#else
  dataInputFileName = new char[strlen(s)+10] ;
  sprintf (dataInputFileName, "%s.%d", s, rank);
#endif
 }
 
 return dataInputFileName ;
}


FILE*  
EngngModel :: giveInputStream ()
   // Returns an input stream on the data file of the receiver.
{
   if (inputStream)
   return inputStream ;
  
   else {
     char *fname = this->giveInputDataFileName ();
     if ((inputStream = fopen(fname,"r")) == NULL)
       {
     fprintf (stderr,"\nEngngModel->giveInputStream:Can't open input stream %s\a\n",
         fname) ;
     exit (1);}
   
     return inputStream ;}
 }
*/

FILE*
EngngModel :: giveOutputStream ()
   // Returns an output stream on the data file of the receiver.
{
   if (! outputStream) {
   char* tmp = tmpnam (NULL);
   _warning2 ("giveOutputStream: using default output stream %s",tmp);
   outputStream = fopen (tmp,"w") ;
  }

   return outputStream ;
}


/*
char*  
EngngModel :: giveLineFromInput (char* line)
//
// reads one line from inputStream - for private use only.
//
{
 char *ptr;
  do {
    fgets(line,OOFEM_MAX_LINE_LENGTH,this->giveInputStream());
    if (line == NULL) {
      fprintf (stderr,"\nEnd of file encountered \a\n");
      exit (1); 
    }
  } while (*line == '#');   // skip comments
 // convert line to lowercase
 for (ptr=line; (*ptr = tolower (*ptr)); ptr++);
  return line;
}
*/

/*
time_t  
EngngModel ::  getTime ()
{
  time_t t;
  t = time(NULL);
  return t;
}
  
clock_t  
EngngModel ::  getClock ()
{
  clock_t t;
  t = clock();
  return t;
}
*/


void
EngngModel :: terminateAnalysis ()
{
 long int tsec; 
 int nsec=0, nmin = 0, nhrs=0;
 FILE* out = this->giveOutputStream ();
 time_t endTime = ::getTime();
 this->timer.stopTimer (EngngModelTimer::EMTT_AnalysisTimer);
 

 fprintf (out, "\nFinishing analysis on: %s\n",ctime(&endTime));
 LOG_FORCED_MSG(oofem_logger, "\n\n____________________________________________________\n");
 // compute real time consumed
 tsec = this->timer.getWtime (EngngModelTimer::EMTT_AnalysisTimer);
 this->timer.convert2HMS (nhrs, nmin, nsec, tsec);
 fprintf (out, "Real time consumed: %03dh:%02dm:%02ds\n",nhrs,nmin,nsec);
 LOG_FORCED_MSG(oofem_logger, "ANALYSIS FINISHED (real time consumed: %03dh:%02dm:%02ds)\n",nhrs,nmin,nsec);
 // compute processor time used by the program
 // nsec = (endClock - startClock) / CLOCKS_PER_SEC;
 tsec = this->timer.getUtime (EngngModelTimer::EMTT_AnalysisTimer);
 this->timer.convert2HMS (nhrs, nmin, nsec, tsec);
 fprintf (out, "User time consumed: %03dh:%02dm:%02ds\n\n\n",nhrs,nmin,nsec);
 LOG_FORCED_MSG(oofem_logger, "                  (user time consumed: %03dh:%02dm:%02ds)\n",nhrs,nmin,nsec);
 
 LOG_FORCED_MSG(oofem_logger, "____________________________________________________\n\a");

 exportModuleManager->terminate();
} 

int
EngngModel :: checkProblemConsistency ()
{
 int result = 1;

 result &= this->checkConsistency();

 int ndomains = this->giveNumberOfDomains();
 for (int i = 1; i<= ndomains; i++)
  result&= this->giveDomain(i)->checkConsistency();

#  ifdef VERBOSE
  if (result)
  VERBOSE_PRINTS("Consistency check","ok")
  else {
  VERBOSE_PRINTS("Consistency check","failed")
  exit (1);
  }
#  endif

 return result;
}



#ifdef __PARALLEL_MODE
#ifdef __USE_MPI
void
EngngModel :: initParallel () 
{
 int len;
 
 MPI_Get_processor_name (processor_name, &len);
 MPI_Comm_rank (MPI_COMM_WORLD, &this->rank);
 MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
#ifdef __VERBOSE_PARALLEL
 OOFEM_LOG_RELEVANT("[%d/%d] Running on %s\n",rank, numProcs, processor_name);
#endif
}

#endif
#endif
 
#ifdef __OOFEG
void EngngModel :: drawYourself (oofegGraphicContext& context) {
 this->giveDomain(context.getActiveDomain())->drawYourself (context);
}

void EngngModel :: drawElements (oofegGraphicContext& context) {
 this->giveDomain(context.getActiveDomain())->drawElements (context);
}

void  EngngModel :: drawNodes (oofegGraphicContext& context) {
 this->giveDomain(context.getActiveDomain())->drawNodes (context);
}

#endif


#ifdef __PARALLEL_MODE
void 
EngngModel::ballanceLoad (TimeStep* atTime)
{
  LoadBallancerMonitor::LoadBallancerDecisionType _d;
  this->giveLoadBallancerMonitor();
  this->giveLoadBallancer();

  _d = lbm->decide();
  if ((_d == LoadBallancerMonitor::LBD_RECOVER) || 
      ((atTime->isTheFirstStep()) && force_load_rebalance_in_first_step)) {
    // determine nwe partitioning
    lb->calculateLoadTransfer();
    // pack e-model solution data into dof dictionaries
    this->packMigratingData(atTime);
    // migrate data 
    this->giveDomain(1)->migrateLoad (lb);
    // renumber itself
    this->forceEquationNumbering();
#if __VERBOSE_PARALLEL
    // debug print
    int i, j, nnodes=giveDomain(1)->giveNumberOfDofManagers();
    int myrank = this->giveRank();
    fprintf (stderr, "\n[%d] Nodal Table\n", myrank);
    for (i=1; i<=nnodes; i++) {
      if (giveDomain(1)->giveDofManager(i)->giveParallelMode()==DofManager_local) 
        fprintf (stderr, "[%d]: %5d[%d] local ", myrank, i, giveDomain(1)->giveDofManager(i)->giveGlobalNumber());
      else if (giveDomain(1)->giveDofManager(i)->giveParallelMode()==DofManager_shared) {
        fprintf (stderr, "[%d]: %5d[%d] shared ", myrank, i, giveDomain(1)->giveDofManager(i)->giveGlobalNumber());
      }
      for (j=1; j<=giveDomain(1)->giveDofManager(i)->giveNumberOfDofs(); j++) {
        fprintf (stderr, "(%d)", giveDomain(1)->giveDofManager(i)->giveDof(j)->giveEquationNumber());
      }
      fprintf (stderr, "\n");
    }
    
#endif
    // unpack (restore) e-model solution data from dof dictionaries
    this->unpackMigratingData(atTime);
    
    
    
    
  }
}
#endif
