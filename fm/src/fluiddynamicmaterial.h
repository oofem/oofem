/* $Header: /home/cvs/bp/oofem/tm/src/transportmaterial.h,v 1.1 2003/04/14 16:01:40 bp Exp $ */
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

#ifndef fluiddynamicmaterial_h 

#include "material.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"
#include "matconst.h"
#include "structuralelement.h"
#include "matstatus.h"

class GaussPoint ;

/**
 This class implements a transport material status information. It is atribute of
 gaussPoint. This is only an abstract class, for every instance of material class 
 there should be specialized derived class, which handles are history variables.
 It only adds attributes common to all "transport problem" material models - the
 state value vectors (both the temp and equlibrium) containing the state values 
 in associated integration point. The corresponding services
 for accessing, setting, initializing and updating these attributes are provided.
 
 For general description of material status, and its role @see MaterialStatus class.
*/
class FluidDynamicMaterialStatus : public MaterialStatus
{

protected:

 /// Equilibrated state vector in reduced form
  FloatArray  deviatoricStressVector ;  // reduced form

public:
 /// Constructor - creates new TransportMaterialStatus with number n, belonging to domain d and IntegrationPoint g.
  FluidDynamicMaterialStatus (int n, Domain* d, GaussPoint* g) ;
  /// Destructor
  ~FluidDynamicMaterialStatus () {}

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
  @param obj pointer to integration point, which invokes this method
  @return contextIOResultType.
  @exception throws an ContextIOERR exception if error encountered.
  */
  contextIOResultType    saveContext (FILE* stream, void *obj = NULL)   ;
 /**
  Restores context of receiver from given stream (the equilibriun stress and strains vectors are restored).
  @param stream stream where to read data
  @param obj pointer to integration point, which invokes this method
  @return contextIOResultType.
  @exception throws an ContextIOERR exception if error encountered.
  */
  contextIOResultType    restoreContext(FILE* stream, void *obj = NULL) ;
  // saves current context(state) into stream

  /// Returns the const pointer to receiver's stateVector attribute
  const FloatArray&  giveDeviatoricStressVector ()         { return deviatoricStressVector ;}
  void         letTempDeviatoricStressVectorBe (const FloatArray& v) 
  { deviatoricStressVector = v ;}
  
  /// Returns "TransportMaterialStatus" - class name of the receiver.
  const char* giveClassName () const { return "FluidDynamicMaterialStatus" ;}
  /// Returns TransportMaterialStatusClass - classType id of receiver.
  classType                giveClassID () const
  { return FluidDynamicMaterialStatusClass; }
 
};


/**
 Abstract base class for all constitutive models for transport problems. It declares common  services provided 
 by all structural material models. The implementation of these services is partly left on derived classes,
 which will implement constitutive model dependent part.
 Some general purpose services are implemented on this level. For details, how to store 
 material model related history variables in integration points, see base class \ref Material documentation.
 
 The constitutive model can in general support several material modes (plane stress, plane strain ,... modes).
 Its capabilities can be examined using hasMaterialModeCapability  service. 
 It is generally assumed, that results obtained from constitutive model services are according to
 valid material mode. This mode is determined from integration point, which is compulsory parameter of all material
 services.
*/
class FluidDynamicMaterial : public Material
{

protected:
public:
 
 /**
  Constructor. Creates material with given number, belonging to given domain.
  @param n material number
  @param d domain to which new material will belong
  */
  FluidDynamicMaterial (int n,Domain* d) : Material (n,d) {}
 /// Destructor.
  ~FluidDynamicMaterial ()                {}
 
 /**
    Computes devatoric stress vector from given strain 
 */
 virtual void computeDeviatoricStressVector (FloatArray& answer, GaussPoint* gp, const FloatArray& eps, TimeStep* tStep) =0;

 /**
    Computes Deviatoric stiffness (derivative of deviatoric stress tensor with respect to strain)
 */
 virtual void giveDeviatoricStiffnessMatrix (FloatMatrix& answer, MatResponseMode,GaussPoint* gp,
					     TimeStep* atTime) = 0;

 /**
  Updates internal state of material according to new state vector.
  @param vec new state vector
  @param gp integration point
  @param tStep solution step
  */
 virtual void updateInternalState (const FloatArray& vec, GaussPoint* gp, TimeStep*);
 /**
  Request material extension. 
  @param ext material extension tested
  @return nonzero if implemented
  */
 virtual int testMaterialExtension (MaterialExtension ext) {return ((ext==Material_FluidDynamicsCapability)?1:0);}
 /// Returns class name of the receiver.
  const char* giveClassName () const { return "FluidDynamicMaterial" ;}
 /// Returns classType id of receiver.
  classType giveClassID ()         const {return FluidDynamicMaterialClass;}
 
#ifdef __OOFEG
#endif
} ;


#define fluiddynamicmaterial_h
#endif





