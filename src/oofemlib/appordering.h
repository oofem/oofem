/* $Header: /home/cvs/bp/oofem/oofemlib/src/Attic/appordering.h,v 1.1.2.1 2004/04/05 15:19:43 bp Exp $ */
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

#ifndef appordering_h
#define appordering_h

#include "cltypes.h"
class EngngModel;
class IntArray;

/**
*/
class ApplicationOrdering {
public:
  enum EquationType {et_standard, et_prescribed};
protected:
public:
  ApplicationOrdering() {}
  virtual ~ApplicationOrdering() {}

  virtual void init (EngngModel* em, EquationID ut, int di, EquationType et = et_standard) = 0;

  /** 
    Returns number of local eqs; ie. those that belong to receiver processor;
    Note that some eqs may be owned by remote processors (some shared nodes,...).
    The summ of local eqs for all processors should give total number of eqs.
    */
  virtual int giveNumberOfLocalEqs () {return 0;}
  /** Returns the total number of eqs of the problem */
  virtual int giveNumberOfGlobalEqs () {return 0;} 

  virtual int giveNewEq (int leq) = 0;
  virtual int giveOldEq (int eq) = 0;
  
  virtual void map2New (IntArray& answer, const IntArray& src, int baseOffset = 0) = 0;
  virtual void map2Old (IntArray& answer, const IntArray& src, int baseOffset = 0) = 0;

  
};

#endif // appordering_h
