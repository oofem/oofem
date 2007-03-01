/* $Header: /home/cvs/bp/oofem/oofemlib/src/nummet.h,v 1.13 2003/04/06 14:08:25 bp Exp $ */
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


//   *****************************
//   *** CLASS NumericalMethod ***
//   *****************************

 
#ifndef nummet_h


#include "femcmpnn.h"
#include "alist.h"
#include "cltypes.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

class EngngModel ; 

/**
This base class is an abstraction for numerical
algorithm. Generally, the particular instances (instances of derived
classes) perform some sequence of numerical operations on given data
to obtain the solution. 
The derived class should declare the interface for specific problem type
(like solution of linear system). The interface usualy consist in
declaring virtual abstract function solve, with parameters corresponding
to problem under consideration. The solve method shoud return
value of NM_Status type. (Other parameters can be provided via instanciateFrom
service, which receives the init record of Engng method).
The data are specified using parameters passed to solve method (so called
mapping). Typically, each particular Engng model instance is responsible
for mapping of its governing equation components to corresponding
numerical components. Such mapping allows the numerical method implementation
to be independent of a particular physical problem by strictly dealing
with numerical components, which are mapped to corresponding physical
components of governing equation, that are  hidden to numerical method.
It should be pointed out, that all numerical methods solving the same
numerical problem use the same compulsory names for the corresponding numerical components - this is
enforced by using the same base problem-specific class. 
It is therefore possible to use any suitable
instance of the Numerical method class to the solve problem, and leave the  whole engineering model code,
including mapping, unchanged, because all instances of the Numerical
method class provide the common interface.
Similarly, a high-level numerical method instance may use services of another
low-level numerical method instance. Th Numerical method
instance may also represent interface to an existing  procedure
written in C or Fortran.
*/
class NumericalMethod : public FEMComponent
{
protected:
 /// Pointer to Engng Model. 
 EngngModel* engngModel ;

public:
 /// Constructor
 NumericalMethod (int i, Domain* d,EngngModel* m):FEMComponent(i,d)
  {engngModel = m ; }        // constructor
 /// Destructor
 ~NumericalMethod () {}                     // destructor

 EngngModel*  giveEngngModel() {return engngModel;}

 // identification 
 /// Returns class name of the receiver.
 const char*  giveClassName () const { return "NumericalMethod" ;}
 /** Returns classType id of receiver.
  @see FEMComponent::giveClassID 
  */
 classType giveClassID () const { return NumericalMethodClass ;}
 
 public:
 };

#define nummet_h
#endif









