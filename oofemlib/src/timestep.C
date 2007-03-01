/* $Header: /home/cvs/bp/oofem/oofemlib/src/timestep.C,v 1.10.4.1 2004/04/05 15:19:44 bp Exp $ */
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

//   file TIMESTEP.CC

#include "timestep.h"
#include "domain.h"
#include "engngm.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

TimeStep :: TimeStep (int n, EngngModel* e, int mn, double tt, double dt, StateCounterType counter)
   // Constructor. Creates a new time step, with number n, and belonging to
   // the time history of s. Used for the initial step (0 or 1).
{
 eModel = e ;
 deltaT = dt ;
 t      = tt ;
 solutionStateCounter = counter;
 number = n;
 version = 0;
 mstepNumber = mn;
}

TimeStep :: TimeStep (EngngModel* e)
{
 eModel = e;
 deltaT = 0.0 ;
 t      = 0.0 ;
 solutionStateCounter = 0;
 number = -1;
 version = 0;
 mstepNumber = 0;
}

TimeStep :: TimeStep (const TimeStep& src)
{
 eModel = src.eModel;
 t      = src.t;
 deltaT = src.deltaT;
 solutionStateCounter = src.solutionStateCounter;
 number = src.number;
 version = src.version;
 mstepNumber = src.mstepNumber;
}

TimeStep&  
TimeStep :: operator=  (const TimeStep& src)
{
 eModel = src.eModel;
 t      = src.t;
 deltaT = src.deltaT;
 solutionStateCounter = src.solutionStateCounter;
 number = src.number;
 version = src.version;
 mstepNumber = src.mstepNumber;

 return *this;
}




TimeStep*  TimeStep :: givePreviousStep ()
   // Not accepted in-line.
{
 if (isTheCurrentTimeStep())
  return eModel->givePreviousStep() ;
 else {
  OOFEM_ERROR ("TimeStep::givePreviousStep Could not return previous step of noncurrent step");
 }
 return NULL; // to make compiler happy
}


int  TimeStep :: isNotTheLastStep ()
   // Returns True if the time history contains steps after the receiver,
   // else returns False.
{
   return  (number != eModel->giveNumberOfSteps()) ;
}

int  TimeStep :: isTheFirstStep ()
{
 // Returns True if the receiver is the first time step, 
 // according to first step number
 // else returns False.

   return  (number == eModel->giveNumberOfFirstStep()) ;
}
   

int  TimeStep :: isIcApply ()
{
 // Returns True if the receiver is the  time step, 
 // when Initial conditions apply
 // else returns False.
 
 return  (number == eModel->giveNumberOfTimeStepWhenIcApply()) ;
}



int  TimeStep :: isTheCurrentTimeStep ()
   // Not accepted in-line.
{
   return this==eModel->giveCurrentStep() ;
}


contextIOResultType    
TimeStep :: saveContext (FILE* stream, void *obj) 
{
  int type_id = TimeStepClass;
  // write class header
  if (fwrite(&type_id,sizeof(int),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
 // write step number
  if (fwrite(&number,sizeof(int),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
 // write meta step number
  if (fwrite(&mstepNumber,sizeof(int),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
  // write time
  if (fwrite(&this->t,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
  // write deltaT
  if (fwrite(&this->deltaT,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
  // write solutionStateCounter
  if (fwrite(&this->solutionStateCounter,sizeof(StateCounterType),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
 
  // return result back
  return CIO_OK;
}

contextIOResultType   
TimeStep ::  restoreContext(FILE* stream, void *obj) 
{
  int class_id;
  // read class header
  if (fread(&class_id,sizeof(int),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
  if (class_id != TimeStepClass) THROW_CIOERR(CIO_BADVERSION);

 // read step number
  if (fread(&number,sizeof(int),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
 // read meta step number
  if (fread(&mstepNumber,sizeof(int),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
  // read time
  if (fread(&this->t,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
  // read deltaT
  if (fread(&this->deltaT,sizeof(double),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
  // read solutionStateCounter
  if (fread(&this->solutionStateCounter,sizeof(StateCounterType),1,stream) != 1) THROW_CIOERR(CIO_IOERR);
  // return result back
  return CIO_OK;

}

