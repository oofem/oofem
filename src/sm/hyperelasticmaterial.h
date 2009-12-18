/* $Header: /home/cvs/bp/oofem/oofemlib/src/structuralmaterial.h,v 1.11.4.1 2004/04/05 15:19:44 bp Exp $ */
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

//   *********************************
//   *** CLASS STRUCTURAL MATERIAL ***
//   *********************************

#ifndef hyperelasticmaterial_h 
#define hyperelasticmaterial_h

#include "structuralmaterial.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"
#include "matconst.h"
#include "matstatus.h"


#include "linearelasticmaterial.h"
#include "dictionr.h"

#include "structuralms.h"

namespace oofem {

/**
 This class implements associated Material Status to MyMaterial.
 It is atribute of matStatusDictionary at every GaussPoint, for which this material 
 is active.
*/
class HyperElasticMaterialStatus : public StructuralMaterialStatus
{
 protected:
  /*
    Declare state variables here
  */
 public: 
  /// Constructor
 //MyMaterialStatus (int n, Domain*d, GaussPoint* g) ;
   HyperElasticMaterialStatus(int n, Domain*d, GaussPoint* g) ;
  /// Destructor
   ~HyperElasticMaterialStatus() ;
  /// Prints the receiver state to stream
  void   printOutputAt (FILE *file, TimeStep* tStep) ;

  
  /* 
     declare state variable access amd modification methods
  */

  // definition
  const char* giveClassName () const { return "HyperElasticMaterialStatus" ;}
  classType             giveClassID () const { return HyperElasticMaterialStatusClass; }
 
  /**
     Initializes the temporary internal variables, describing the current state according to 
     previously reached equilibrium internal variables.
  */
  virtual void initTempStatus ();
  /**
     Update equilibrium history variables according to temp-variables.
     Invoked, after new equilibrium state has been reached.
  */
  virtual void updateYourself(TimeStep*); // update after new equilibrium state reached
  
  // saves current context(state) into stream
  /**
     Stores context of receiver into given stream. 
     Only non-temp internal history variables are stored.
     @param stream stream where to write data
     @param obj pointer to integration point, which invokes this method
     @return contextIOResultType.
  */
  contextIOResultType    saveContext (FILE* stream, void *obj = NULL);
  /**
     Restores context of receiver from given stream. 
     @param stream stream where to read data
     @param obj pointer to integration point, which invokes this method
     @return contextIOResultType.
  */
  contextIOResultType    restoreContext(FILE* stream, void *obj = NULL);
 
}; 




class HyperElasticMaterial : public StructuralMaterial
{
 protected:
 double K,G; // declare material properties here

 public:
  HyperElasticMaterial (int n, Domain* d);


/**
 Method for computing 1d  stifness matrix of receiver.
 Default implementation computes 3d stifness matrix using give3dMaterialStiffnessMatrix and
 reduces it to 1d stiffness using reduce method described above.
 Howewer, this reduction is quite time consuming and if it is possible, 
 it is recomended to overload this method and provide direct method for computing
 particular stifness matrix.
 @param answer stifness matrix
 @param form material response form
 @param mode material response mode
 @param gp integration point, which load history is used
 @param atTime time step (most models are able to respond only when atTime is current time step)
*/
 virtual void give3dMaterialStiffnessMatrix (FloatMatrix& answer,
                   MatResponseForm,MatResponseMode,GaussPoint* gp,
                   TimeStep* atTime);
  /**
     Computes the real stress vector for given total strain and integration point. 
     The total strain is defined as strain computed directly from displacement field at given time.
     The stress independent parts (temperature, eigen strains) are substracted in constitutive
     driver.
     The service should use previously reached equilibrium history variables. Also
     it should update temporary history variables in status according to newly reached state.
     The temporary history variables are moved into equilibrium ones after global structure
     equlibrium has been reached by iteration process.
     @param answer contains result
     @param form material response form
     @param gp integration point 
     @param reducedStrain strain vector in reduced form
     @param tStep current time step (most models are able to respond only when atTime is current time step)
  */
  virtual void giveRealStressVector (FloatArray& answer, MatResponseForm, GaussPoint*, 
                                     const FloatArray& , TimeStep*);


MaterialStatus* CreateStatus (GaussPoint* gp) const ;
//MaterialStatus* CreateStatus (GaussPoint *gp) const {return new HyperElasticMaterialStatus (1,HyperElasticMaterial::domain,gp);}
  
/**
 Requests material mode capability.
 @param mode material mode requested
 @return nonzero if available
 */
 virtual int hasMaterialModeCapability (MaterialMode ) ;
  /// Returns class name of the receiver.
 const char* giveClassName () const { return "HyperElasticMaterial" ;}
 /// Returns classType id of receiver.
 classType giveClassID ()         const {return HyperElasticMaterialClass; }//MyMaterialClass;}
/**
 Initializes receiver acording to object description stored in input record.
 The density of material is read into property dictionary (keyword 'd')
 */
 IRResultType initializeFrom (InputRecord* ir);
 };

} // end namespace oofem
#endif
