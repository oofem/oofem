/* $Header: /home/cvs/bp/oofem/oofemlib/src/crosssection.h,v 1.11 2003/04/06 14:08:23 bp Exp $ */
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


//   ***************************
//   *** CLASS CROSSSECTOION ***
//   ***************************

#ifndef crosssection_h 

#include "femcmpnn.h"
#include "material.h"
//#include "perfectlyplasticmaterial.h"
#include "gausspnt.h"
//#include "element.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"

#define THICKNESS 400
#define WIDTH     401
#define BEAM_SHEAR_COEFF 402
#define AREA      403
#define INERTIA_MOMENT_Y 404
#define INERTIA_MOMENT_Z 405
#define TORSION_MOMENT_X 406

// id's for layered model
#define TOPZCOORD 498
#define BOTTOMZCOORD 499




enum updateGpMode {
 doNotUpdate,
 Update
};

/**
 Base abstract class representing cross section in finite element mesh.
 
 The main idea, why cross section has been introduced, is to hide all details 
 of cross section description from particular element. Generally elements 
 do not comunicate directly with material but comunicate through cross section interface, 
 which therefore can perform necessary integration (for example over layers of fibers).
 
 The derived classes are supposed to be base cross section classes for particular
 type of analyses. They should declare general interface methods necessary.
 
 In particular cross section implementation, where is necessary to perform integration 
 over cross section volume (over layers, fibers, ...) and therefore generally one must keep 
 complete load history in these integration points, the concept of master-slave integration 
 points should be used.
 Integration point generally can contain list of slave integration points
 therefore is called as master point. Slaves are used for example to implement
 layered or fibred cross sections by cross section class. Then in one
 "macro" master gauss point, cross section creates few slaves (one per layer)
 and puts them into master list. When cross sections completes requests for
 particular master integration point, it performs integration over layers.
 It therefore calls material class for each layer, sending corresponding
 slave as parameter and integrates results.
 @see GaussPoint class for more detail.
 */
class CrossSection : public FEMComponent
{
/*
   This abstract class implements a cross section in a finite element problem. A cross
 section  is an attribute of a domain. It is usually also attribute of many 
 elements. This class is base class for SimpleCrossSection, LayeredCrossSection,
 FibredCrossSection, ....

 DESCRIPTION
   The attribute 'propertyDictionary' contains all the properties of a 
 cross section, like its area or thickness.

 TASK
 - Returning a properties of cross section like thickness or area.
*/
protected:
 /**
  Dictionary for  storing cross section parameters (like dimensions).
  More preferably, (due to slow access into dictionary values) one should use
  corresponding varibles declared inside class
  */
 Dictionary*    propertyDictionary ;

public:
 /**
  Constructor. Creates cross section with number n belonging to domain d.
  @param n cross section number
  @param d domain
  */
 CrossSection (int n,Domain* d) : FEMComponent(n,d)
  { propertyDictionary = new Dictionary() ;}
 /// Destructor.
 ~CrossSection ()                { delete propertyDictionary ;}

 /**
  Returns the value of cross section property 'aProperty'. Property must be identified 
  by unique int id.
  @aProperty id of peroperty requested
  @return property value
  */
 virtual double   give (int) ;
 /** 
   Returns true if stiffness matrix of receiver is symmetric
   Default implementation returns true.
   It can be moved to base Cross section class in the future.
   */
 virtual bool isCharacteristicMtrxSymmetric (MatResponseMode rMode, int mat);
 
 // identification and auxiliary functions
 /// Returns class name of the receiver.
 const char*    giveClassName () const      {return "CrossSection" ;}
 /// Returns classType id of receiver.
 classType giveClassID ()         const      {return CrossSectionClass;}
 /// Initializes receiver acording to object description stored in input record.
 IRResultType initializeFrom (InputRecord* ir);
 /// Prints receiver state on stdout. Useful for debugging.
 void     printYourself () ;
   
 // definition
 /**
  Returns a newly allocated cross section, with type depending on parameter.
  Creates new object for following classes SimpleCrossSection and HeatCrossSection
  otherwise calls directly CreateUsrDefCrossSectionOfType global function to allocate
  new instance of cross section of given type.
  @param aClass string with cross section model name
  @return newly allocated cross section model with required type.
  @see CreateUsrDefCrossSectionOfType function.
  */
 CrossSection*      ofType (char*) ; 

 /**
  Returns nonzero, if receiver implements required extension.
  @ext required extensison
  @return nonzero, if supported, zero otherwise
  */
  virtual int testCrossSectionExtension (CrossSectExtension ext) {return 0;}
  //virtual int hasStructuralCapability () {return 0;}
  //virtual int hasHeatCapability       () {return 0;}

 // store & restore context functions
 /**
  Stores context of receiver into given stream. 
  Default implementation simply retypes obj parameter 
  (which is void*) to integration point pointer and 
  invokes saveContext function on this integration point.
  Derived classes, if they use slave concept, must invoke saveContext function also on all
  defined slaves of master.
  @param stream stream where to write data
  @param mode determines ammount of info in stream
  @param obj pointer to integration point, which related context should be stored
  @return contextIOResultType value.
  @exception throws an ContextIOERR exception if error encountered
  */
 contextIOResultType    saveContext (DataStream* stream, ContextMode mode, void *obj = NULL);
 /**
  Restores context of receiver from given stream. 
  Default implementation simply retypes obj parameter 
  (which is void*) to integration point pointer and 
  invokes restoreContext function on this integration point.
  Derived classes, if they use slave concept, must invoke restoreContext function also on all
  defined slaves of master.
  @param stream stream where to read data
  @param mode determines ammount of info in stream
  @param obj pointer to integration point, which invokes this method
  @return contextIOResultType value.
  @exception throws an ContextIOERR exception if error encountered
  */
 contextIOResultType    restoreContext(DataStream* stream, ContextMode mode, void *obj = NULL);

 /**
  Returns the integration point corresponding value in Reduced form.
  @param answer contain corresponding ip value, zero sized if not available
  @param aGaussPoint integration point
  @param type determines the type of internal variable
  @param atTime time step
  @return nonzero if o.k, zero otherwise
 */
 virtual int giveIPValue (FloatArray& answer, GaussPoint* aGaussPoint, InternalStateType type, TimeStep* atTime) 
  {return aGaussPoint->giveMaterial()->giveIPValue (answer, aGaussPoint, type, atTime);}
 /**
  Returns the corresponding integration point value size in Reduced form.
  @param type determines the type of internal variable
  @param mat corresponding material
  @return nonzero if o.k, zero otherwise
 */
 virtual int giveIPValueSize (InternalStateType type, GaussPoint* aGaussPoint) 
  {return aGaussPoint->giveMaterial()->giveIPValueSize (type, aGaussPoint);}
 /**
  Returns the mask of reduced indexes of Internal Variable component .
  @param answer mask of Full VectorSize, with components beeing the indexes to reduced form vectors.
  @param type determines the internal variable requested (physical meaning)
  @returns nonzero if ok or error is generated for unknown mat mode.
  */
 virtual int giveIntVarCompFullIndx (IntArray& answer, InternalStateType type, MaterialMode mmode, Material* mat)
  {return mat->giveIntVarCompFullIndx (answer, type, mmode);}
 /**
  Returns the type of internal variable (scalar, vector, tensor,...).
  @param type determines the type of internal variable
  @returns type of internal variable
  */
 virtual  InternalStateValueType giveIPValueType (InternalStateType type, Material* mat) 
  {return mat->giveIPValueType (type);}
 
#ifdef __PARALLEL_MODE
  /**
  Pack all necessary data of integration point (according to element parallel_mode) 
  into given communication buffer. The corresponding material model service for particular integration point
  is invoked. The nature of packed data is material model dependent.
  Typically, for material of "local" response (response depeneds only on integration point local state)
  no data are exchanged. For "nonlocal" constitutive models the send/receive of local values which
  undergo averaging is performed between local and corressponding remote elements.
  @param buff communication buffer
  @param stepN solution step
  @param ip integration point
  */
  virtual int packUnknowns (CommunicationBuffer& buff, TimeStep* stepN, GaussPoint* ip)
  { return ip->giveMaterial()->packUnknowns (buff, stepN, ip);}
  /**
  Unpack and updates all necessary data of given integration point (according to element parallel_mode) 
  into given communication buffer. 
  @see packUnknowns service.
  @param buff communication buffer
  @param stepN solution step.
  @param ip integration point
  */
  virtual int unpackAndUpdateUnknowns (CommunicationBuffer& buff, TimeStep* stepN, GaussPoint* ip)
  {return ip->giveMaterial()->unpackAndUpdateUnknowns (buff, stepN, ip);}
 /**
  Estimates the necessary pack size to hold all packed data of receiver.
  The corresponding material model  service is invoked. The 
  nature of packed data is typically material model dependent.  
  */
 virtual int estimatePackSize (CommunicationBuffer& buff, GaussPoint* ip)
  {return ip->giveMaterial()->estimatePackSize (buff, ip);}
#endif

 friend class Material;
} ;


#define crosssection_h
#endif

