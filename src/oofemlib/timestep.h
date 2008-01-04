/* $Header: /home/cvs/bp/oofem/oofemlib/src/timestep.h,v 1.11 2003/04/06 14:08:26 bp Exp $ */
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
/*
 The original idea for this class comes from 
  Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 PhD Thesis, EPFL, Lausanne, 1992.
*/

//   ***********************
//   *** CLASS TIME STEP ***
//   ***********************


#ifndef timestep_h

#include "femcmpnn.h"
#include "engngm.h"
#include "compiler.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif


/**
 Class representing solution step. The Solution step instance may represent either 
 time step, load increment, or load case depending on current Engineering model used.

 Solution step maintain the reference to correspoding Engineering model class instance.
 It maintain also its "intrinsic time" and corresponding time increment. The meaning of these 
 values is dependent on current Engineering model used. The time may represent either
 current time, load increment number or load case number. See corresponding 
 Engng model reference for details.
 
 Some components (typically integration points real stresses or integration points nonlocal values)
 are computationally wery demanding. Because in typical code, there are number of requests for same value 
 during the computation process, it may be efficient to store these values and compute them only once.
 The principal problem is to recognize, when is necessary to re-compute these stored values to reflect 
 newly reached state. This cannot be determined form solution step "time", because solution step may 
 represent for example load increment, inside which typically many iterations are needed to reach 
 convergence. For this purpose, a concept of solution state counters is provided.
 Whenever the solution state changes, the engineering model updates the solution state counter.
 The solution state counter is guaranteed to grow up smoothly (it newer decreases) during solution process.
 Other components of program (integration points) can then store their computationally expensive values
 but have to store also corresponding solution state counter value valid when these were computed.
 Then their can easily check for difference between freezed solution state counter for their value with 
 current solution state requested from solution step and recompute the values if necessary.
*/
class TimeStep
{
/*
   This class implements a step in a time history.
 DESCRIPTION :
   A time step n is characterized by its date 't', and by the time increment
   'dt' between step n-1 and step n.
   The time step knows the time history 'scheme' it belongs to. 

 TASKS :
   - returning its time 't' and its time increment 'dt'.

*/
protected:
 /// Engng model reference.
 EngngModel             *eModel ;
 /// Current intrinsic time. 
 double                  t ;
 /// Current intrinsic time increment.
 double                  deltaT ;
 /// Solution state counter.
 StateCounterType        solutionStateCounter;
 /// receiver's number
 int number;
 /** receiver's version, used for special applicatons; default set to 0.
  Typically, new version of same step is generated after adaptive restart, when
  the restarted step is equlibrated on new domain
  */
 int version;
  /// corresponding meta step number
 int mstepNumber;

public:
 /*
   Constructor. Creates a new solution step.
   @param n solution step number
   @param e reference to corresponding engng model
   @param mn meta step number
   @param tt intrinsic time
   @param dt intrinsic time increment
   @param counter solution state counter 
  */
 TimeStep (int n, EngngModel* e, int mn, double tt,double dt,StateCounterType counter) ;       // constructors
 TimeStep (const TimeStep&);  
 TimeStep (EngngModel* e);
 TimeStep&  operator=  (const TimeStep&);  // assignment: cleanup and copy

 /// Returns receiver's number
 int giveNumber () {return number;}
 /// Returns receiver's version
 int giveVersion () {return version;}
 /// Returns receiver's meta step number
 int giveMetaStepNumber() {return mstepNumber;}
 /**
  Returns class name of receiver.
  @param s buffer to store name
  @return pointer to s parameter filled with name
  */
 const char* giveClassName () const   { return "TimeStep" ;}
 /// Returns poiter to previous solution step.
 TimeStep*  givePreviousStep() ;
 /// Returns solution step associated time.
 double     giveTime ()                { return t ;}
 /// Returns solution step associated time increment.
 double     giveTimeIncrement ()       { return deltaT ;}
 /// Sets solution step time increment.
  void       setTimeIncrement (double newDt) {deltaT = newDt;}
 /// Sets solution step time.
  void       setTime (double newt) {t = newt;}

 
 /**
  Tests if solution step is not the last step.
  @return nonzero if not last step, zero otherwise.
  */
 int        isNotTheLastStep () ;
 /** 
  Tests if receiver is first step.
  @return nonzero if receiver is the first step, zero otherwise.
  */
 int        isTheFirstStep () ;
 /**
  Test if receiver is solution step when initial conditions should apply.
  @return nonzero if ic apply, zero otherwise.
  */
 int        isIcApply () ;
 /**
  Returns current solution state counter.
  */
 StateCounterType giveSolutionStateCounter () {return solutionStateCounter;}
 /// Updates solution state counter.
 void       incrementStateCounter() {solutionStateCounter ++;}
 /// Increments receiver's version
 void       incrementVersion () {version++;}
 // LoadResponseMode giveLoadResponseMode () {return eModel->giveLoadResponseMode();} 
 // indicate whether incremetal  or total load is used

 /**
  Tests if receiver is currentsolution step.
  @returns nonzero, if receiver is current step, zero otherwise.
  */
 int        isTheCurrentTimeStep () ;
 //int        requiresNewLhs () {return eModel->requiresNewLhs();}
 IRResultType initializeFrom (InputRecord* ir) {return IRRT_OK;}
 /** Stores receiver state to output stream. 
  Receiver should write class-id first in order to allow test
  whether correct data are then restored.
  @param stream output stream 
  @param mode determines ammount of info required in stream (state, definition,...)
  @param obj special parameter, used only to send particular integration
  point to material class version of this method. Except this 
  case, obj parameter is always NULL pointer.
  @exception throws an ContextIOERR exception if error encountered.
  */
  contextIOResultType    saveContext (DataStream* stream, ContextMode mode, void *obj = NULL) ;
 /** Restores the receiver state previously written in stream.
  @see saveContext member function.*/
 contextIOResultType    restoreContext(DataStream* stream, ContextMode mode, void *obj = NULL) ;

} ;

#define timestep_h
#endif
