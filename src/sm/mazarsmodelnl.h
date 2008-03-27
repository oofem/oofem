/* $Header: /home/cvs/bp/oofem/sm/src/mazarsmodelnl.h,v 1.8 2003/04/06 14:08:31 bp Exp $ */
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

  //   *********************************************************
  //   *** CLASS NONLOCAL MAZARS MODEL FOR CONCRETE ************
  //   *********************************************************

#ifndef mazarsmodelnl_h 
#define mazarsmodelnl_h 
 
#include "mazarsmodel.h"
#include "structuralnonlocalmaterialext.h"
#include "cltypes.h"

  class GaussPoint ;

/**
   This class implements associated Material Status to MazarsNLModel.
*/
class MazarsNLMaterialStatus : public MazarsMaterialStatus, public StructuralNonlocalMaterialStatusExtensionInterface
{

 protected:
 
  /// Equivalent strain for avaraging 
  double localEquivalentStrainForAverage;

 public:
 
  /// constructor
  MazarsNLMaterialStatus (int n, Domain*d, GaussPoint* g);
  /// Destructor
  ~MazarsNLMaterialStatus ();

  /// Prints the receiver state to given stream
  void   printOutputAt (FILE *file, TimeStep* tStep) ;

  /// Returns the local  equivalent strain to be averaged
  double giveLocalEquivalentStrainForAverage ()     {return localEquivalentStrainForAverage;}
  /// Sets the localEquivalentStrainForAverage to given value
  void   setLocalEquivalentStrainForAverage (double ls) {localEquivalentStrainForAverage = ls;}

  // definition
  const char* giveClassName () const { return "MazarsNLMaterialStatus" ;}
  classType             giveClassID () const { return IsotropicDamageMaterialStatusClass; }
 

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
     Corresponding parent method invoked.
     @param stream stream where to write data
     @param mode determines ammount of info required in stream (state, definition,...)
     @param obj pointer to integration point, which invokes this method
     @return contextIOResultType.
  */
  contextIOResultType    saveContext (DataStream* stream, ContextMode mode, void *obj = NULL);
  /**
     Restores context of receiver from given stream. 
     Corresponding parent method invoked.
     @param stream stream where to read data
     @param mode determines ammount of info required in stream (state, definition,...)
     @param obj pointer to integration point, which invokes this method
     @return contextIOResultType.
  */
  contextIOResultType    restoreContext(DataStream* stream, ContextMode mode, void *obj = NULL);
  /**
     Interface requesting service.
     In the case of nonlocal constitutive models, 
     the use of multiple inheritance is assumed. Typically, the class representing nonlocal 
     constitutive model status is derived both from class representing local status and from class 
     NonlocalMaterialStatusExtension or from one of its derived classes 
     (which declare services and variables corresponding to specific analysis type).
     @return In both cases, this function returns pointer to this object, obtained by 
     returning adress of component or using pointer conversion from receiver to base class 
     NonlocalMaterialStatusExtension. If no nonlocal extension exists, NULL pointer is returned.
  */
  virtual Interface* giveInterface (InterfaceType) ;
 
}; 


/**
   This class implements a Nonlocal Mazars Damage Model for Concrete 
   Model based on nonlocal averaging of equivalent strain.
*/ 
class MazarsNLMaterial : public MazarsMaterial, public StructuralNonlocalMaterialExtensionInterface
{

 protected:
  /// Interaction radius, related to the nonlocal characteristic length of material.
  double R;

 public:
 
  /// Constructor
  MazarsNLMaterial (int n,Domain* d) ;
  /// Destructor
  ~MazarsNLMaterial () ;
 
  // identification and auxiliary functions
  const char*    giveClassName () const { return "MazarsNLMaterial" ;}
  classType giveClassID ()         const { return MazarsMaterialClass;}
 
  /// Initializes the receiver from given record
  IRResultType initializeFrom (InputRecord* ir);
  /** Interface requesting service */
  virtual Interface* giveInterface (InterfaceType);


  /**
     Computes the equivalent nonlocal strain measure from given strain vector (full form).
     @param kappa return param, comtaining the corresponding equivalent strain
     @param strain total strain vector in full form
     @param gp integration point
     @param atTime time step
  */
  virtual void computeEquivalentStrain (double& kappa, const FloatArray& strain, GaussPoint* gp, TimeStep* atTime) ;
  /**
     Computes the equivalent local strain measure from given strain vector (full form).
     @param kappa return param, comtaining the corresponding equivalent strain
     @param strain total strain vector in full form
     @param gp integration point
     @param atTime time step
  */
  void computeLocalEquivalentStrain (double& kappa, const FloatArray& strain, GaussPoint* gp, TimeStep* atTime) 
    {MazarsMaterial::computeEquivalentStrain (kappa, strain, gp, atTime);}

  /**
     Implements the service updating local variables in given integration points, 
     which take part in nonlocal average process. Actually, no update is necessary,
     because the value used for nonlocal averaging is strain vector used for nonlocal secant stiffness
     computation. It is therefore necessary only to store local strain in corresponding status.
     This service is declared at StructuralNonlocalMaterial level.
     @param equivalentStrain equivalent strain vector in given integration point.
     @param gp integration point to update.
     @param atTime solution step indicating time of update.
  */
  virtual void updateBeforeNonlocAverage (const FloatArray& strainVector, GaussPoint* gp, TimeStep* atTime) ;

  /**
     Computes the value of nonlocal weight function in given point. 
     @param src coordinates of source point.
     @param coord coordinates of point, where nonlocal weight function is evaluated.
     @return value of weight function.
  */
  virtual double computeWeightFunction (const FloatArray& src, const FloatArray& coord) ;
  /**
     Determines, whether receiver has bounded weighting function (limited support)
     @return true if weighting function bounded, zero otherwise
  */
  virtual int hasBoundedSupport () {return 1;}
  /**
     Determines the width (radius) of limited support of weighting function
  */
  virtual void giveSupportRadius (double& radius) {radius = this->R;}

#ifdef __PARALLEL_MODE
  /**
     Updates domain before nonloc average (using updateDomainBeforeNonlocAverage service)
     to ensure, that the localStrainVectorForAverage variable is correctly updated and
     pack this localStrainVectorForAverage into given buffer.
     @see Material::packUnknowns for description.
     @param buff communication buffer
     @param stepN solution step
     @param ip integration point
  */
  int packUnknowns (CommunicationBuffer& buff, TimeStep* stepN, GaussPoint* ip);
  /**
     Unpack localStrainVectorForAverage value from given buffer.
     @see Material::unpackAndUpdateUnknowns service.
     @param buff communication buffer
     @param stepN solution step.
     @param ip integration point
  */
  int unpackAndUpdateUnknowns (CommunicationBuffer& buff, TimeStep* stepN, GaussPoint* ip);
  /**
     Estimates the necessary pack size to hold all packed data of receiver.
  */
  int estimatePackSize (CommunicationBuffer& buff, GaussPoint* ip);
#endif
 
  /// Creates the corresponding material status
  MaterialStatus* CreateStatus (GaussPoint *gp) const {return new MazarsNLMaterialStatus (1,MazarsMaterial::domain,gp);}

 protected: 
  /** 
      Perfoms initialization, when damage first appear. The Le characteristic length is
      set equal to 1.0, it doesnot matter - nonlocal approach is used. The only 
      same value with reference length (which is used in local model, which 
      computeDmaga function is reused).
      @param kappa scalar measure of strain level
      @param totalStrainVector current total strain vector
      @param gp integration point
  */
  void initDamaged (double kappa, FloatArray& totalStrainVector, GaussPoint* gp);

} ;


#endif




