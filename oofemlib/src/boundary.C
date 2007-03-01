/* $Header: /home/cvs/bp/oofem/oofemlib/src/boundary.C,v 1.8.4.1 2004/04/05 15:19:43 bp Exp $ */
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
//   file BOUNDARY.C
//

 
#include "boundary.h"
#include "timestep.h"
#include "loadtime.h"
#include "verbose.h"
#include "cltypes.h"

/*
double  BoundaryCondition :: give (char u, TimeStep* stepN)
   // Returns the value at stepN of the prescribed value of the kinematic 
   // unknown 'u'. Returns 0 if 'u' has no prescribed value.
{
   double value,factor ;
 
   if (! prescribedValueDictionary)
      _error ("give: prescribedValueDictionary is not defined");

   if (prescribedValueDictionary -> includes(u)) {
      value  = prescribedValueDictionary -> at(u) ;
      factor = this -> giveLoadTimeFunction() -> at(stepN->giveTime()) ;
   if ((stepN->giveLoadResponseMode()==IncrementOfLoad) && (!stepN->isTheFirstStep()))
    //factor -= this->giveLoadTimeFunction()->at(stepN->givePreviousStep()->giveTime()) ;
    factor -= this->giveLoadTimeFunction()->at(stepN->giveTime()-stepN->giveTimeIncrement()) ;

      return value*factor ;}
   else
      return 0. ;
}
*/

double  BoundaryCondition :: give (Dof* dof, ValueModeType mode, TimeStep* stepN)
   // Returns the value at stepN of the prescribed value of the kinematic 
   // unknown 'u'. Returns 0 if 'u' has no prescribed value.
{
   double factor ;

   /*
     factor = this -> giveLoadTimeFunction() -> at(stepN->giveTime()) ;
     if ((isUnknownTypeModeIncremental(mode)) && (!stepN->isTheFirstStep()))
     //factor -= this->giveLoadTimeFunction()->at(stepN->givePreviousStep()->giveTime()) ;
     factor -= this->giveLoadTimeFunction()->at(stepN->giveTime()-stepN->giveTimeIncrement());
   */
   factor = this -> giveLoadTimeFunction() -> evaluate(stepN, mode) ;
   return prescribedValue*factor ;

/* 
   double value,factor ;
   char u;

   if (! prescribedValueDictionary)
      _error ("give: prescribedValueDictionary is not defined");

   u = cltypesGiveUnknownTypeKey (type);
   if (prescribedValueDictionary -> includes(u)) {
      value  = prescribedValueDictionary -> at(u) ;
      factor = this -> giveLoadTimeFunction() -> at(stepN->giveTime()) ;
   if ((isUnknownTypeModeIncremental(mode)) && (!stepN->isTheFirstStep()))
    //factor -= this->giveLoadTimeFunction()->at(stepN->givePreviousStep()->giveTime()) ;
    factor -= this->giveLoadTimeFunction()->at(stepN->giveTime()-stepN->giveTimeIncrement());

      return value*factor ;}
   else
      return 0. ;
*/

}


int BoundaryCondition :: isImposed (TimeStep* tStep)
{
// returs a value of isImposedTimeFunction, indicating whether b.c. is imposed or not
// in given time (nonzero indicates imposed b.c.).

 if (isImposedTimeFunction) {
  int flag;
  flag = (domain -> giveLoadTimeFunction(isImposedTimeFunction)->evaluate(tStep, VM_Total) != 0.);
  return flag;
 } else {
  // zero value indicates default behaviour -> b.c. is imposed
  // anytime
  return 1;
 }
}

IRResultType 
BoundaryCondition::initializeFrom (InputRecord* ir)
// Sets up the dictionary where the receiver stores the conditions it
// imposes.
{
//   char   key [32] ; // 'key' is eventually of size 1, but some words that are
//   double value ;    // read in the meantime in the data file can be larger !
//   int    nCond,i,j ;
//
//   prescribedValueDictionary = new Dictionary() ;
//   i     = 0 ;
//   nCond = this->readInteger("conditions",++i) ;

//   for (j=1 ; j<=nCond ; j++) {
//      this -> readString("conditions",++i,key) ;
//      value = this->read("conditions",++i) ;
//      prescribedValueDictionary -> add(key[0],value) ;}

 const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
 IRResultType result;                   // Required by IR_GIVE_FIELD macro

 GeneralBoundaryCondition::initializeFrom(ir);

 isImposedTimeFunction = 0;
 IR_GIVE_OPTIONAL_FIELD (ir, isImposedTimeFunction, IFT_BoundaryCondition_IsImposedTimeFunct, "isimposedtimefunction"); // Macro

 if (ir->hasField(IFT_BoundaryCondition_PrescribedValue, "prescribedvalue")) {
  IR_GIVE_FIELD (ir, prescribedValue, IFT_BoundaryCondition_PrescribedValue, "prescribedvalue"); // Macro
 } else {
  IR_GIVE_FIELD (ir, prescribedValue, IFT_BoundaryCondition_PrescribedValue, "d"); // Macro
 }

  return IRRT_OK;

}



int
BoundaryCondition::giveInputRecordString(std::string &str, bool keyword)
{
 char buff[1024];

 GeneralBoundaryCondition::giveInputRecordString(str, keyword);
 sprintf(buff, " isimposedtimefunction %d prescribedvalue %e", 
     this -> isImposedTimeFunction, this -> prescribedValue);
 str += buff;

 return 1;
}

