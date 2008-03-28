/* $Header: /home/cvs/bp/oofem/sm/src/peak.h,v 1.4 2003/04/06 14:08:31 bp Exp $ */
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

//   ***************************
//   *** CLASS PEAK FUNCTION ***
//   ***************************
#ifndef peak_h
#define peak_h

#include "loadtime.h"

class PeakFunction : public LoadTimeFunction
{
/*
   This class implements a function that is 0 everywhere, except in a single
   point.
 DESCRIPTION
   't' is the only abscissa (time) where the function is not 0. 'value' is
   the value of the function 't'. Both 't' and 'value' are pointers, rather
   than numbers, so that their state (initialized or not) can be checked.
*/

private:
  double  t ;
  double  value ;

public:
  PeakFunction (int i,Domain* d) : LoadTimeFunction(i,d)
    { t=0.0 ; value=0.0 ;}
  ~PeakFunction ()                       {}
  
  //      void    getCoefficients () ;
  IRResultType initializeFrom (InputRecord* ir);
  classType   giveClassID () const { return PeakFunctionClass;}
  const char*  giveClassName ()  const { return "PeakFunction" ;}   

  double  __at (double) ;

} ;


#endif // peak_h





