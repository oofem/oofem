/* $Header: /home/cvs/bp/oofem/oofemlib/src/materialmapperinterface.h,v 1.6 2003/04/06 14:08:25 bp Exp $ */
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


//   **********************************************
//   *** CLASS MATERIAL MODEL MAPPING INTERFACE ***
//   **********************************************

#ifndef materilmapperinterface_h 

#include "compiler.h"
#include "cltypes.h"
#include "interface.h"

class Domain;
class Element;
class TimeStep;

/**
 The class representing the general material model adaptive mapping interface.
 The basic task is to define the algorithm for mapping of the internal material model  
 variables from one (old) mesh to the given IP of new mesh.
*/
class MaterialModelMapperInterface : public Interface {

protected:

public:
 /// Constructor
 MaterialModelMapperInterface () : Interface () {}
 /// Destructor
 virtual ~MaterialModelMapperInterface() {}
 /** Maps the required internal state variables from 
  old mesh oldd to given ip. The result is stored in gp status.
  @param gp Integration point belonging to new domain which values will be mapped
  @param oldd old mesh reference
  @param tStep time step
  @return nonzero if o.k.
  */
 virtual int MMI_map (GaussPoint* gp, Domain* oldd, TimeStep* tStep) = 0;
 /** Updates the required internal state variables from previously mapped values. 
  The result is stored in gp status. This map and update splitting is necessary,
  for example for nonlocal models tahe local quantity to be averaged must be mapped in all ips
  and then update can happen, because it may depend on nonlocal variable, which is computed
  from local values.
  @param gp Integration point belonging to new domain which values will be mapped
  @param oldd old mesh reference
  @param tStep time step
  @param elemGPVec vector passed to MMI_update at material level, probably computed from primary unknowns
   (for structural elements this represent strain vector).
  @return nonzero if o.k.
  */
 virtual int MMI_update (GaussPoint* gp, TimeStep* tStep, FloatArray* elemGPVec= NULL) = 0;
 /**
  Finishes the mapping for given time step. Used to perform cleanup. 
  Typically some mappers reguire to compute some global mesh data related to
  current step, which are valid for example to all IPs - so they are computed only once for
  all IPs, stored and they need to be dealocated. These mappers are typically class variables,
  but their finish is invoked by all members.
  */
 virtual int MMI_finish (TimeStep* tStep) = 0;
 
protected:
};

#define materilmapperinterface_h
#endif






