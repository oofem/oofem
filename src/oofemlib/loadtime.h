/* $Header: /home/cvs/bp/oofem/oofemlib/src/loadtime.h,v 1.9 2003/04/06 14:08:25 bp Exp $ */
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


//   ********************************
//   *** CLASS LOAD-TIME FUNCTION ***
//   ********************************
 

#ifndef loadtime_h

#include "femcmpnn.h"
#include "domain.h"

/**
 Abstract base class representing load time function. Classes derived from Load class typically 
 describe load from spatial point of view. The purpose of introducing load time function is to express
 variation of some components in time. Load time function typically belongs to domain and is 
 attribute of one or more loads. Generally load time function is real function of time (\f$y=f(t)\f$).
*/
class LoadTimeFunction : public FEMComponent
{
/*
   This abstract class is the superclass of the classes that implement
   functions  y = f(t) , where t is the time. These function are used for
   weighing the loads (see TJR Hughes, "The Finite Element Method", p 677).
   A load-time function is attribute of the domain. It is usually also at-
   tribute of one or more loads.
 DESCRIPTION
   The parameters which describe the function are defined by the subclasses
   of this class.
 TASK
   Returning the value 'y' at any abscissa 't'.
  */
 
public:
 
 /**
  Constructor. Creates load time function with given number, belonging to given domain.
  @param n load time function number
  @param d domain to which new object will belongs.
  */
 LoadTimeFunction (int i,Domain* d) : FEMComponent(i,d) {} 
 /// Destructor
 virtual ~LoadTimeFunction ()  {}
 
 // computations 
 /**
  Returns the value of load time function at given time. Abstract service. 
  Must be implemented by derived classes.
  @param t time
  @return load time function value
  */
 double     evaluate (TimeStep* atTime, ValueModeType mode);
 
 // definition of a function
 
 /**
  Returns a newly allocated load time function, with type depending on parameter.
  Creates new object for following classes ConstantFunction
  otherwise calls directly CreateUsrDefLoadTimeFunctionOfType global function to allocate
  new instance of load time function of given type.
  @param aClass string with load time function ID  name
  @return newly allocated load time function of required type.
  @see CreateUsrDefLoadTimeFunctionOfType function.
 */
 LoadTimeFunction*  ofType (char*) ;
 /// Returns classType id of receiver.
 classType   giveClassID () const { return LoadTimeFunctionClass;}
 /// Returns class name of the receiver.
 const char*  giveClassName () const { return "LoadTimeFunction" ;}
 /**
  Initializes receiver acording to object description stored in input record.
  Must be implemented in derived classes
  */
 IRResultType initializeFrom (InputRecord* ir) {return IRRT_OK;}


 /**
    Returns the value of load time function at given time. Abstract service. 
    Must be implemented by derived classes.
    @param t time
    @return load time function value
 */
 virtual double     __at (double)            { return 0. ;}
 /**
    Returns the first time derivative of load time function at given time. Abstract service. 
    Must be implemented by derived classes.
    @param t time
    @return load time function value
 */
 virtual double    __derAt (double) {return 0.;}
 /**  
    Returns the second time derivative of load time function at given time. Abstract service. 
    Must be implemented by derived classes.
    @param t time
    @return load time function value
 */
 virtual double    __accelAt (double) {return 0.;}

} ;

#define loadtime_h
#endif









