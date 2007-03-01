/* $Header: /home/cvs/bp/oofem/tm/src/staggeredproblem.h,v 1.1.4.1 2004/04/05 15:19:53 bp Exp $ */
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


#ifndef staggeredproblem_h

#include "engngm.h"
#include "alist.h"
#include "cltypes.h"
#include "inputrecord.h"
#include "fieldmanager.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

/**
 Implementation of general sequence (staggered) problem. The problem consists in sequence of
 low level problems (slaves) which are executed sequentially and where the results 
 of particular slave depends on the results of previous slaves in sequence.
 Typical example is heat&mass transfer analysis followed by mechanical one, which
 takes into account the temperature field from the first analysis.

 The sequence problem is represented by this SerialProblem class. It maintains list
 of subsequent (slave) problems and it is executes the slave problems. It is responsible
 for solutionn step generation and synchronization between slave problems.
 The transfer of required state variables is done by mapping of corresponding variables 
 between problem domains. This allows to to transfer primary (nodal) values of one problem to 
 integration points of subsequent problem or to use completely different discretizations for
 slave problems.
 
 Since the master problem is responsible for synchronization, it is responsible for
 generation the solution steps. Therefore, the solution step specification, as weel as
 relevant meta step attributes are specified at master level. 
 //To avoid confusion, 
 //the slaves are treated in so-called maintained mode. In this mode, the attributes and 
 //meta step attributes are taken from the master. The local attributes, even if specified,
 //are ignored. 
*/
class StaggeredProblem : public EngngModel 
{
protected:
 int nModels;
 AList<EngngModel>*  emodelList;
 double deltaT;
 FieldManager fieldManager;
 char** inputStreamNames;
 /// Associated time function for time step increment
 int dtTimeFunction ;

public:
 /**
  Constructor. Creates Engng model with number i belonging to domain d.
  */
 StaggeredProblem (int i, EngngModel* _master = NULL) ;    // constructor
 /// Destructor.
 virtual ~StaggeredProblem ()  ;      // destructor

 /**
  Sets context output mode of receiver.
  @param contextMode domain context mode.
  */
 void               setContextOutputMode (ContextOutputMode contextMode) ;
 /**
  Sets user defined context output mode (it sets contextOutputMode to contextOutputMode), 
  setting contextOutputStep to given value.
  @param cStep new context output step
  */
 void               setUDContextOutputMode (int cStep);
 /**
  Sets domain mode to given mode.
  @param mode domain mode.
  */
 void               setProblemMode (problemMode mode);
  /// Sets the renumber flag to TRUE
 void                setRenumberFlag();
  
 // solving
 /**
  Solves problem for given time step. Should assemble characteristic matrices and vectors 
  if necessary and solve problem using appropriate numerical method. After finishing solution,
  this->updateYourself function for updating solution state and then this->terminate
  function (for updating nodal and element values) should be called.
  */
 virtual void               solveYourselfAt (TimeStep*) ;
 //virtual int                requiresNewLhs () {return 1;}
 /**
  Updates internal state after finishing time step. (for example total values may be 
  updated according to previously solved increments).  Then element values are also updated
  (together with related integration points and material statuses).
  */
 virtual void               updateYourself (TimeStep* stepN);
 /**
  Provides the oportunity to initialize state variables stored in element
    integration points acording to
   initial conditions using function initializeYourself() on element level.
  Should be called when curent time step is time step when IC will aply 
  (see EngngModel::giveNumberOfTimeStepWhenIcApply)
  somewhere from solveYourselfAt function). Implementation must be provided.
  Default implementation is empty.
  */
 virtual void               initializeYourself (TimeStep*) {}
  /**
  Initializes the newly generated discretization state acording to previous solution.
  This process should typically include restoring old solution, instanciating newly
  generated domain(s) and by mapping procedure. 
  */
  virtual int                initializeAdaptive (int stepNumber) {return 0;}
  /**
  Prints the ouput of the solution step (using virtual this->printOutputAtservice) 
  to the stream detemined using this->giveOutputStream() method
  and calls exportModuleManager to do output.
  */
  virtual void                       doStepOutput(TimeStep*);

 /**
  Initializes whole problem acording to its description stored in inputStream.
  Prints header, opens the outFileName, instanciate itself the receicer using 
  using virtual initializeFrom service and instancites all problem domains.
 */
  int instanciateYourself (DataReader* dr, InputRecord* ir, char* outFileName, char* desc) ;
  /**
  Initializes receiver acording to object description in input reader.
  InitString can be imagined as data record in component database
  belonging to receiver. Receiver may use value-name extracting functions 
  to extract particular field from record.*/
  virtual IRResultType initializeFrom (InputRecord* ir);
  /**
  Update receiver attributes according to step metaStep attributes.
  Allows the certain parameters or attributes to be updated for particular metastep.
  The metastep provides the attributes record, from which the corresponding attributes can
  be read. The service takes TimeStep as parameter, from which corresponding MetaStep is
  requested. It is recomended, to implement this service in such way, that multiple calls
  for steps belonging to same MetaStep does not change response. 
  The default implementation updates the numerical method attributes.
  @param TimeStep time step.
  */
 virtual void updateAttributes (TimeStep*);

 /** 
  Stores the  state of model to output stream. Stores not only the receiver state,
  but also same function is invoked for all DofManagers and Elements in associated
  domain. Note that by storing element  context also contexts of all associated
  integration points (and material statuses) are stored.
  Stored context is associated with current time step. One time step can have only 
  one associated context. Multiple call to saveContext within same time step
  owerride previously saved context for this step.
  By default the stream paprameter is used to store data and is not closed. 
  If stream is NULL, new file descriptor is created and this must be also closed at the end. 
  @param stream - context stream. If NULL then new file descriptor will be openned and closed
  at the end else the stream given as parameter will be used and not closed at the end.
  @param mode determines ammount of info required in stream (state, definition,...)
  @return contextIOResultType.
  @exception throws an ContextIOERR exception if error encountered
  */
  virtual contextIOResultType                saveContext (DataStream *stream, ContextMode mode, void *obj = NULL) ;
 /**
  Restores the  state of model from output stream. Restores not only the receiver state,
  but also same function is invoked for all DofManagers and Elements in associated
  domain. Note that by restoring element  context also contexts of all associated
  integration points (and material statuses) are restored.
  Each context is associated with unique time step. Only one context per time step is
  allowed. Restore context function will restore such contex, which is related 
  (through its step number) to time step number given in obj parameter. 
  Restoring context will change current time step in order to correspond to newly restored
  context.
  @param stream context file
  @param mode determines ammount of info required in stream (state, definition,...)
  @param obj context corresponding to  time step number (retyped as (void*)) will be restored.
  @return contextIOResultType.
  @exception throws an ContextIOERR exception if error encountered.
  */
  virtual contextIOResultType    restoreContext (DataStream* stream, ContextMode mode, void* obj = NULL) ;
   /**
   Updates domain links after the domains of receiver have changed. Used mainly after 
   restoring context - the domains may change and this service is then used
   to update domain variables in all components belonging to receiver
   like errorestimators, solvers, etc, having domains as attributes.
   */
   virtual void updateDomainLinks() ;
 /** 
  Prints output of receiver to ouput domain stream, for given time step.
  Corresponding function for element gauss points is invoked
  (gaussPoint::printOutputAt).
  */
 virtual void                  printOutputAt (FILE *, TimeStep*) ;


 // input / output
 /// Prints stete of receiver. Usefull for debugging.
 void printYourself () ;

  /** DOF printing routine. Called by DofManagers to print Dof specific part.
  Dof class provides component printing routines, but emodel is responsible
  for what will be printed at DOF level.
  @param stream output stream
  @param iDof dof to be processed
  @param atTime solution step
  */
 virtual void printDofOutputAt (FILE* stream, Dof* iDof, TimeStep* atTime) {}
  TimeStep* giveNextStep ();
 TimeStep*  giveSolutionStepWhenIcApply();


      // identification 
  /// Returns class name of the receiver.
 const char*  giveClassName () const { return "StaggeredProblem" ;}
 /// Returns classType id of receiver.
 classType giveClassID () const { return StaggeredProblemClass ;}
 /// Returns nonzero if receiver does incremental analysis.
 virtual int isIncremental () {return 0;}
  /// Returns nonzero if nonlocal stiffness option activated.
 virtual int useNonlocalStiffnessOption () {return 0;}
 /**
  Indicates type of non linear computation (total or updated formulation).
  This is used for example on Nodal level to update coordinates 
  if updated formulation 
  is done, or on element level, when non linear contributions are computed.
  */
 virtual fMode giveFormulation () {return UNKNOWN;}  // for non-linear computation
 /*
  Returns Load Response Mode of receiver.
  This value indicates, whether nodes and elements should assemble
  total or incremental load vectors.
  
 virtual  LoadResponseMode giveLoadResponseMode () {return TotalLoad;}
 */

 /**
    Returns time function for time step increment.
    Used time function should provide step lengths as function of step number. 
    Initial step with number 0 is considered as [ -dt(0), 0 ], first step is [ 0, dt(1) ], ...
 */
 LoadTimeFunction*  giveDtTimeFunction () ;
 
 /**
    Returns the timestep length for given step number n, initial step is number 0
 */
 double giveDeltaT (int n);



#ifdef __OOFEG   
 void               drawYourself (oofegGraphicContext& context);  
 void               drawElements (oofegGraphicContext& context);
 void               drawNodes (oofegGraphicContext& context);
  /**
  Shows the sparse structure of required matrix, type == 1 stiffness.
  */
  virtual void       showSparseMtrxStructure (int type, oofegGraphicContext& context, TimeStep* atTime) {}
#endif

 /** Allows programmer to test some receiver's internal data, before computation begins.
  @return nonzero if receiver check is o.k. */
 virtual int checkConsistency () ; 

  /**Returns i-th slave problem */
 virtual EngngModel* giveSlaveProblem (int i);
  /**Returns number of slave problems */
  virtual int giveNumberOfSlaveProblems() {return nModels;}

 /// Returns number of first time step used by receiver.
 virtual int        giveNumberOfFirstStep () {if (master) return master->giveNumberOfFirstStep(); else return 1;}
 /// Returns the time step number, when initial conditions should apply.
 virtual  int       giveNumberOfTimeStepWhenIcApply() {if (master) return master->giveNumberOfTimeStepWhenIcApply(); 
 else return 0;}


protected:
 int instanciateSlaveProblems ();
};

#define staggeredproblem_h
#endif









