/* $Header: /home/cvs/bp/oofem/sm/src/targe2interface.h,v 1.3 2003/04/06 14:08:31 bp Exp $ */
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

//   ******************************
//   *** CLASS TARGE2 INTERFACE ***
//   ******************************

#ifndef targe2interface_h 
#define targe2interface_h 

#include "mesherinterface.h"

class TimeStep;

/**
 This class represents the interface to Targe2 mesh generation package.
 This interface is primarly responsible for two main tasks:
 - to create input mesher file, containing all informations including the mesh density informations
 based on informations from remeshing criteria.
 - possibly to launch the mesher and transform its output to oofem input (using targe2oofem)
*/
class Targe2Interface : public MesherInterface {
 
public:
 /// Constructor
 Targe2Interface () : MesherInterface () {}
 /// Destructor
 virtual ~Targe2Interface() {}
 
 /// Runs the mesh generation, mesh will be written to corresponding domain din file
 virtual int createMesh (Domain* d, TimeStep* tStep, int domainNumber, int domainSerNum);

 
protected:
 /// Creates the mesher input, containing the required mesh density informations.
 int createInput (Domain* d, TimeStep* stepN);
};

#endif // targe2interface_h 
