/* $Header: /home/cvs/bp/oofem/oofemlib/src/structuralms.h,v 1.9 2003/04/06 14:08:26 bp Exp $ */
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


#ifndef structuralms_h 
#define structuralms_h


#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "matstatus.h"
#include "gausspnt.h"

/*
   This class implements a structural material status information. It is atribute of
 gaussPoint. This is only an abstract class, for every instance of material class 
 there should be specialized derived class, which handles are history variables.

 DESCRIPTION
   This is a base class for all material statuses coresponding to materials derived from
 structural matarial class.
 It defines stress and strain vectors and their increments.
 Functions for accessing these components are defined.
 
 TASKS
   This is abstract class - only basic functionality is supported like:
 - maintaining and providing access to stress and strain vectors 
   (including their increments)
 - storing and restoring status on tape
 - printingYourself()
 - updating Yourself after a new equlibrium state has been reached.

 REMARK
   Materials statuses are atributes of GaussPoints, they are stored in 
  MatStatus variable
*/

class GaussPoint ; class Dictionary; 
class Domain; 

/**
 This class implements a structural material status information. It is atribute of
 gaussPoint. This is only an abstract class, for every instance of material class 
 there should be specialized derived class, which handles are history variables.
 It only adds attributes common to all "structural analysis" material models - the
 strain and stress vectors (both the temp and equlibrium). The corresponding services
 for accessing, setting, initializing and updating these attributes are provided.
 
 For general description of material status, and its role @see MaterialStatus class.
*/
class StructuralMaterialStatus : public MaterialStatus
{
protected:

 /// Equilibrated strain vector in reduced form
  FloatArray  strainVector ;  // reduced form
  /// Equilibrated stress vector in reduced form
  FloatArray  stressVector ;  // reduced form
  /// Temp stress vector in reduced form
  FloatArray  tempStressVector;// increments are used mainly in nonlinear analysis
  /// Temp strain vector in reduced form
  FloatArray  tempStrainVector;// to find balanced state. for update call gp->update()

public:
 /// Constructor - creates new StructuralMaterialStatus with number n, belonging to domain d and IntegrationPoint g.
  StructuralMaterialStatus (int n, Domain* d, GaussPoint* g) ;
  /// Destructor
  ~StructuralMaterialStatus () ;

 /// Print receiver's output to given stream.
  void   printOutputAt (FILE *, TimeStep*) ;

  /**
  Initializes the temporary internal variables (stresss and strains vectors), 
  describing the current state according to 
  previously reached equilibrium internal variables.
  */
  virtual void initTempStatus ()         ;
  /**
  Update equilibrium history variables according to temp-variables.
  Invoked, after new equilibrium state has been reached.
  */
  virtual void updateYourself(TimeStep*) ; // update after new equilibrium state reached


 /**
  Stores context of receiver into given stream (the equilibriun stress and strains vectors are stored). 
  Generally, only non-temp internal history variables should be stored.
  @param stream stream where to write data
  @param mode determines ammount of info required in stream (state, definition,...)
  @param obj pointer to integration point, which invokes this method
  @return contextIOResultType.
  @exception throws an ContextIOERR exception if error encountered.
  */
  contextIOResultType    saveContext (DataStream* stream, ContextMode mode, void *obj = NULL)   ;
 /**
  Restores context of receiver from given stream (the equilibriun stress and strains vectors are restored).
  @param stream stream where to read data
  @param mode determines ammount of info required in stream (state, definition,...)
  @param obj pointer to integration point, which invokes this method
  @return contextIOResultType.
  @exception throws an ContextIOERR exception if error encountered.
  */
  contextIOResultType    restoreContext(DataStream* stream, ContextMode mode, void *obj = NULL) ;
  // saves current context(state) into stream

  /// Returns the const pointer to receiver's strainVector attribute
  const FloatArray&  giveStrainVector ()         { return strainVector ;}
  /// Returns the const pointer to receiver's stressVector
  const FloatArray&  giveStressVector ()         { return stressVector ;}
  /// Returns the const pointer to receiver's tempStrainVector
  const FloatArray&  giveTempStrainVector()      { return tempStrainVector;}
  /// Returns the const pointer to receiver's tempStressVector
  const FloatArray&  giveTempStressVector()      { return tempStressVector;}
  /// Assigns strainVector to given vector v
  void         letStrainVectorBe (const FloatArray& v) 
  { strainVector = v ;}
  /// Assigns stressVector to given vector v
  void         letStressVectorBe (const FloatArray& v) 
  { stressVector = v ;}
  /// Assigns tempStressVector to given vector v
  void         letTempStressVectorBe (const FloatArray& v) 
  { tempStressVector = v ;}
  /// Assigns tempStrainVector to given vector v
  void         letTempStrainVectorBe (const FloatArray& v) 
  { tempStrainVector = v ;}
  
  // definition
  /// Returns "StructuralMaterialStatus" - class name of the receiver.
  const char* giveClassName () const { return "StructuralMaterialStatus" ;}
  /// Returns StructuralMaterialStatusClass - classType id of receiver.
  classType                giveClassID () const
  { return StructuralMaterialStatusClass; }
 
};

#endif // structuralms_h




