/* $Header: $ */
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

#ifndef materialinterface_h
#include "femcmpnn.h"


/**
   Abstract base class representing (moving) material interfaces.
   Its typical use to model moving interface (such as free surface)
   in a fixed-grid methods (as typically used in CFD).
   The basic tasks are representation of interface and its updating.
 */
class MaterialInterface : public FEMComponent
{
 protected:
 public:
  /** Constructor. Takes two two arguments. Creates 
      MaterialInterface instance with given number and belonging to given domain.
      @param n component number in particular domain. For instance, can represent  
      node number in particular domain.
      @param d domain to which component belongs to 
  */
  MaterialInterface (int n,Domain* d) : FEMComponent (n,d) {}

  /**
     Updates the position of interface according to state reached in given solution step.
   */
  virtual void updatePosition (TimeStep* atTime) = 0;
  /**
    Updates element state after equlibrium in time step has been reached.
    All temporary history variables,
    which now describe equlibrium state should be copied into equilibrium ones.
    The existing internal state is used for update.
   */
  virtual void updateYourself (TimeStep* tStep) = 0;
  /**
     Computes critical time step induced by receiver integration algorithm
  */
  virtual double computeCriticalTimeStep (TimeStep*) = 0;
  
};

#define materialinterface_h
#endif
