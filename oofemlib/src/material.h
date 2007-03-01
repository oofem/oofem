/* $Header: /home/cvs/bp/oofem/oofemlib/src/material.h,v 1.14 2003/04/14 16:00:47 bp Exp $ */
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


//   **********************
//   *** CLASS MATERIAL ***
//   **********************

#ifndef material_h 

#include "femcmpnn.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"
#include "matconst.h"
#include "structuralelement.h"
#include "matstatus.h"
#include "nonlocalmaterialext.h"
#include "timestep.h"

#define STRAIN_STEPS 10.0

class GaussPoint ;


/**
 Abstract base class for all material models. Declares the basic common interface 
 to all material models. Derived classes should expand this interface, because they are
 assumed to be  base classes for analysis specific tasks (for example mechanical or
 thermal analysis).
 Instance of integration point class is assumed to be implicit argument to 
 all method, depending on internal state in point of consideration.
 To provide oportunity for storing arbitrary material model related history variables
 in integration points, associated material status class is introduced.
 Each new material model class should be declared together with its associted status class
 (derived from MaterialStatus class). This status can be seen as simple container,
 storing necessary history variables and providing some access and modification methods.
 Each integration point can contain material status. Material model should create
 unique copy of its associated status in each integration point. 
 Because integration point is parameter of all mesages to material model
 class, material model therefore can easily access  all history variables it needs.
 @see MaterialStaus class
 @see GaussPoint class
*/
class Material : public FEMComponent
{
/*
   This class implements a material in a finite element problem. A material
   is an attribute of a domain. It is usually also attribute of many elements.
 DESCRIPTION
 The attribute 'propertyDictionary' contains all the properties of a mate-
   rial, like its Young modulus, its mass density or poisson ratio.
 TASK
  - indicate, whether there required material mode is valid for receiver
   (method hasMaterialModeCapability). Note: for some material models and linear materials
  there need not exist support for assembling material char matrix at material level,
  all is handled properly at crossSection level (_2dBeam mode, 3dShellMode, ...).
  But this function must indicate, whwether mode is valid or not for real stress computation.
*/
protected:
 /**
  Property dictionary. 
  Can be used to store constant material parameters, which
  are same for all integration points. 
  Note: Try to avoid using, more preferably use separate variables to
  store marerial params, because of wery slow access to dictionary.
  */
 Dictionary*    propertyDictionary ;

 /** 
     Casting time. For solution time less than casting time the material 
     is assumed to have no stiffness etc. This attribute is declared here, 
     but support for this functionality must be incorporated by particular 
     material model
 */
 double castingTime;
public:
 
 /**
  Constructor. Creates material with given number, belonging to given domain.
  @param n material number
  @param d domain to which new material will belong
  */
 Material (int n,Domain* d) : FEMComponent(n,d)
        { propertyDictionary = new Dictionary() ;}
 /// Destructor.
 ~Material ()                { delete propertyDictionary ;}
 
 // standart material stiffness matrices
 /**
  Computes characteristic matrix of receiver in given integration point.
  @param answer contains result
  @param form material response form
  @param mode  material response mode
  @param gp integration point
  @param atTime time step 
  */
 virtual void  giveCharacteristicMatrix (FloatMatrix& answer, 
                                         MatResponseForm form,
                                         MatResponseMode mode,
                                         GaussPoint* gp,
                                         TimeStep* atTime);
 /**
    Computes the characteristic value of receiver in given integration point, respecting its history.
    The algorithm should use temporary or equlibrium  history variables stored in integration point status
    to compute and return required result. 
    @param mode material response mode
    @param gp integration point
    @param atTime time step (most models are able to respond only when atTime is current time step)
 */
 virtual double  giveCharacteristicValue (MatResponseMode mode,
                                          GaussPoint* gp,
                                          TimeStep* atTime) ;

 /** 
   Returns true if stiffness matrix of receiver is symmetric
   Default implementation returns true.
   */
 virtual bool isCharacteristicMtrxSymmetric (MatResponseMode rMode) {return true;}

 // uptates MatStatus to newly reched (equilibrium) state
 /**
  Updates internal state of material in integration point after finishing time step.
  Default implementation  extract material status from integration point,
  and invokes updateYourself function on it.
  @param gp integration point, where to update state
  @param atTime time step
  */
 virtual void updateYourself (GaussPoint* gp, TimeStep* atTime) ;
 
 // non-standart - returns time independent material constant
 /**
  Returns the value of material property 'aProperty'. Property must be identified 
  by unique int id.
  @aProperty id of peroperty requested
  @return property value
  */
 virtual double   give (int aProperty) ;
 /**
    Returns casting time of the receiver
 */
 double giveCastingTime () {return this->castingTime;}
 /**
    Returns true, if material activated for given solution step
 */
 bool   isActivated (TimeStep* atTime) {if (atTime) return (atTime->giveTime() >= this->castingTime); else return true;}

  // identification and auxiliary functions
 /**
  Returns nonzero if receiver is non linear
  */
  virtual int hasNonLinearBehaviour ()   { return 0 ;}
  //virtual int hasStructuralCapability () { return 0;}
  //virtual int hasHeatCapability       () { return 0;}
 /**
  Returns nonzero, if receiver implements required extension.
  @ext required extensison
  @return nonzero, if supported, zero otherwise
  */
  virtual int testMaterialExtension      (MaterialExtension ext) {return 0;}
 /**
  Tests, if material supports material mode.
  @param mode required material mode
  @return nonzero if supported, zero otherwise
  */
  virtual int hasMaterialModeCapability (MaterialMode mode) ;
 /**
  Returns the integration point corresponding value in Reduced form.
  @param answer contain corresponding ip value, zero sized if not available
  @param aGaussPoint integration point
  @param type determines the type of internal variable
  @param type determines the type of internal variable
  @returns nonzero if ok, zero if var not supported
  */
 virtual int giveIPValue (FloatArray& answer, GaussPoint* aGaussPoint, InternalStateType type, TimeStep* atTime) 
  {answer.resize(0); return 0;}
 /**
  Returns the corresponding integration point  value size in Reduced form.
  @param type determines the type of internal variable
  @returns var size, zero if var not supported
  */
 virtual int giveIPValueSize (InternalStateType type, GaussPoint* aGaussPoint) 
  {return 0;}
 /**
  Returns the mask of reduced indexes of Internal Variable component .
  @param answer mask of Full VectorSize, with components beeing the indexes to reduced form vectors.
  @param type determines the internal variable requested (physical meaning)
  @returns nonzero if ok or error is generated for unknown mat mode.
  */
 virtual int giveIntVarCompFullIndx (IntArray& answer, InternalStateType type, MaterialMode mmode)
  {answer.resize(0); return 0;}
 /**
  Returns the type of internal variable (scalar, vector, tensor,...).
  @param type determines the type of internal variable
  @returns type of internal variable
  */
 virtual  InternalStateValueType giveIPValueType (InternalStateType type) 
  {_error ("giveIPValueType: unsupported InternalStateType"); return ISVT_UNDEFINED;}
  
 /// Returns class name of the receiver.
  const char*    giveClassName () const { return "Material" ;}
 /// Returns classType id of receiver.
 classType giveClassID ()         const {return MaterialClass;}
 /**
  Initializes receiver acording to object description stored in input record.
  The density of material is read into property dictionary (keyword 'd')
  */
 IRResultType initializeFrom (InputRecord* ir);
 /** Setups the input record string of receiver
  @param str string to be filled by input record
  @param keyword print record keyword (default true)
  */
 virtual int giveInputRecordString(std::string &str, bool keyword = true);
 /// Prints receiver state on stdout. Useful for debugging.
  void     printYourself () ;

  // definition
 /**
  Returns a newly allocated material, with type depending on parameter.
  Creates new object for following classes: IsotropicLinearElasticMaterial and IsotropicLinearHeatMaterial
  otherwise calls directly CreateUsrDefMaterialOfType global function to allocate
  new instance of material of given type.
  @param aClass string with material model name
  @return newly allocated material model with required type.
  @see CreateUsrDefMaterialOfType function.
  */
  Material*      ofType (char*) ; 
  // store & restore context functions
 /**
  Stores context of receiver into given stream. This method is called from 
  integration point saveContext function, to store material related staus in 
  integration point. Integration point passes itself as obj parameter, when
  invokes this method. Default implementation simply retypes obj parameter 
  (which is void*) to integration point pointer, requests material status, and 
  invokes saveContext function on this status.
  @param stream stream where to write data
  @param mode determines ammount of info in stream (state, definition,...)
  @param obj pointer to integration point, which invokes this method
  @return contextIOResultType.
  @exception throws an ContextIOERR exception if error encountered.
  */
  contextIOResultType    saveContext (DataStream* stream, ContextMode mode, void *obj = NULL);
 /**
  Restores context of receiver from given stream. This method is called from 
  integration point restoeContext function, to restore material related staus in 
  integration point. Integration point passes itself as obj parameter, when
  invokes this method. Default implementation simply retypes obj parameter 
  (which is void*) to integration point pointer, requests material status, and 
  invokes restoreContext function on this status.
  @param stream stream where to read data
  @param mode determines ammount of info in stream (state, definition,...)
  @param obj pointer to integration point, which invokes this method
  @return contextIOResultType.
  @exception throws an ContextIOERR exception if error encountered.
  */
  contextIOResultType    restoreContext(DataStream* stream, ContextMode mode, void *obj = NULL);
  
  // initialize gp record at the begining of new load Increment
 /**
  Initialize integration point (must be inside  material region associated to receiver)
  when new time step begins. Simply calls initTempStatus function with integration point
  as parameter.
  @param gp integration point to initialize
  */
  virtual void initGpForNewStep (GaussPoint* gp);
 /**
  Returns material status of receiver in given integration point.
  If status doe not exist yet, it is created using CreateStatus  member function.
  @param gp Returns reference to material status belonging to integration 
  point gp.
  @return material status associated with given integration point.
  */
  virtual MaterialStatus* giveStatus (GaussPoint* gp) const;
 /*
  In the case of nonlocal constitutive models, 
  the use of multiple inheritance is assumed. Typically, the class representing nonlocal 
  constitutive model is derived both from class representing local model and from class 
  NonlocalMaterialExtension or from one of its derived classes 
  (which declare services and variables corresponding to specific analysis type).
  Or artenatively, material model can introduce stronger form of this relation,
  it can posses object of NonlocalMaterialExtension. 
  @return In both cases, this function returns pointer to this object, obtained by 
  returning adress of component or using pointer conversion from receiver to base class 
  NonlocalMaterialExtension. If no nonlocal extension exists, NULL pointer is returned.
  */
 //virtual NonlocalMaterialExtension* giveNonlocalMaterialExtensionPtr () {return NULL;}

#ifdef __PARALLEL_MODE
  /**
  Pack all necessary data of integration point (according to element parallel_mode) 
  into given communication buffer. The nature of packed data is material model dependent.
  Typically, for material of "local" response (response depeneds only on integration point local state)
  no data are exchanged. For "nonlocal" constitutive models the send/receive of local values which
  undergo averaging is performed between local and corressponding remote elements.
  @param buff communication buffer
  @param stepN solution step
  @param ip integration point
  */
  virtual int packUnknowns (CommunicationBuffer& buff, TimeStep* stepN, GaussPoint* ip)
  { return 1;}
  /**
  Unpack and updates all necessary data of given integration point (according to element parallel_mode) 
  into given communication buffer. 
  @see packUnknowns service.
  @param buff communication buffer
  @param stepN solution step.
  @param ip integration point
  */
  virtual int unpackAndUpdateUnknowns (CommunicationBuffer& buff, TimeStep* stepN, GaussPoint* ip)
  {return 1;}
 /**
  Estimates the necessary pack size to hold all packed data of receiver.
  */
 virtual int estimatePackSize (CommunicationBuffer& buff, GaussPoint* ip)
  {return 0;}
#endif


 /**
  Creates new copy of associated status and inserts it into given integration point.
  @param gp Integration point where newly created status will be stored.
  @return reference to new status.
  */
  virtual MaterialStatus* CreateStatus (GaussPoint* gp) const
  {return NULL;}

protected:
 /**
  Initializes temporary variables stored in integration point status
  at the begining of new time step.
  (Temporary history variables (they describe state of material during 
  solution of time step) are initialized according to history variables, which 
  describe state corresponding to previous equilibrium solution).
  Default implementation simply extracts status from integration point and
  calls its initTempStatus method.
  */
  virtual void initTempStatus (GaussPoint* gp);

} ;


#define material_h
#endif
