/* $Header: /home/cvs/bp/oofem/sm/src/microplane.h,v 1.4 2003/04/06 14:08:31 bp Exp $ */
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

//   ************************
//   *** CLASS MICROPLANE ***
//   ************************
 
#ifndef microplane_h

#include "gausspnt.h"
#include "microplanematerial.h"
#include "flotarry.h"
#include "element.h"
#include "matstatus.h"
#include "tdictionary.h"

class Element ; class Material ; class LayeredCrossSection;

/**
 Class representing microplane integration point in finite element program.
 Microplane always belongs to macro integration point and 
 represents micro material level. Microplanes are created by particular microplane
 material model.
 Because all corresponding microplanes in all macro integration points
 subjected to particular material share the same properties (same normals, 
 same integration weights), these are stored only once at material model level.
 Thus requests for microplane integration weights and normal are forwarded
 to corresponding material model.

 Generally, every integretion point must hold
 its own copy of history variables (which are related to corresponding
 material model used). These material type dependent history variables 
 are stored in material type related material status, which can be 
 managed by integration point. 
 Each material model class should introduce related material status class 
 (derived from material status class or from its children), where necesary
 history variables are kept and can be accesed by material.
 Material class then creates unique copy of related status in all necessary
 integration points. Because integration point is compulsory parameter of
 all member functions of material class, particular material then can 
 easily access its associated status from integration point and therefore its
 history variables for particular integration point.
*/
class Microplane : public GaussPoint
{

protected:

public:
 /**
  Constructor.  Creates microplane integration point belonging 
  to given element, with given number, integration weight, coordinates and material mode.
  @param e element to which integration point belongs to.
  @param n integration point number
  @param mode material mode
  */
 Microplane  (IntegrationRule* ir,int n,MaterialMode mode);
 /// Destructor
 ~Microplane  () ;                                  // destructor
 
 /// Returns  integration weight of receiver.
 double       giveWeight ()               
  { return ((MicroplaneMaterial*)this->giveMaterial())->giveMicroplaneIntegrationWeight (this) ;}
  /// Returns normal of microplane.
 void giveMicroplaneNormal (FloatArray& answer) 
  { ((MicroplaneMaterial*)this->giveMaterial())->giveMicroplaneNormal (answer, this);}
 /**
  Prints output of receiver to file. Corresponding printOutputAt  function for
  associated status is called. Overloaded in order to avoid report microplane state.
  */
 void         printOutputAt (FILE *, TimeStep*) ;

 /// Returns classType id of receiver.
 classType giveClassID () const {return MicroplaneClass;}
 /// Returns class name of the receiver.
 const char*  giveClassName () const {return "Microplane" ;}
 ///Initializes receiver acording to object description stored in input record.
 IRResultType initializeFrom (InputRecord* ir) {return IRRT_OK;}
 
} ;

#define microplane_h
#endif





