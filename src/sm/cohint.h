/*
CohesiveInterface added by Milan Jirasek on 1 Feb 2010

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2010   Borek Patzak                                       



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

//   *************************************************************************
//   *** CLASS COHESIVE INTERFACE FOR THE COHESIVE PARTICLE MODEL ************
//   *************************************************************************

#ifndef cohint_h 
#define cohint_h

#include "material.h"
#include "linearelasticmaterial.h"
#include "structuralmaterial.h"
#include "structuralms.h"
#include "cltypes.h"

// material contant's keys for give()
class GaussPoint ;

namespace oofem {

//---------------------------------------------------------------------------------------------------
// c l a s s   CohesiveInterfaceMaterialStatus
//---------------------------------------------------------------------------------------------------

/**
 This class implements associated Material Status to CohesiveInterfaceMaterial.
 It is atribute of matStatusDictionary at every GaussPoint, for which this material 
 is active.
*/
class CohesiveInterfaceMaterialStatus : public StructuralMaterialStatus
{
/*
 This class implements associated Material Status to CohesiveInterfaceMaterial.
 It is atribute of matStatusDictionary at every GaussPoint, for which this material 
 is active.
 DESCRIPTION:
   Idea used there is that we have variables
   describing:
   1) state at previous equilibrium state (variables without temp)
   2) state during searching new equilibrium (variables with temp)
   when we start search new state from previous equilibrium one we copy 
   non-tem variables into temp ones. And after we reach new equilibrium 
   (now decribed by temp variables) we copy tem-var into non-tepm ones
   (see function updateYourself).
   
   variables description:

  kappa      - scalar measure of the largest strain level ever reached in material
  tempKappa

*/

protected:
 /// scalar measure of the largest equivalent strain ever reached in material
 double kappa;
  /// non-equilibrated scalar measure of the largest equivalent strain
  double tempKappa;
  /// damage variable
 double damage;
  /// non-equilibrated damage variable
 double tempDamage;
 /// plastic shear strains
 double gam1p, gam2p;
 /// non-equilibrated plastic shear strains
 double tempGam1p, tempGam2p;
public: 
 /// Constructor
 CohesiveInterfaceMaterialStatus (int n, Domain*d, GaussPoint* g) ;
  /// Destructor
 ~CohesiveInterfaceMaterialStatus ();

  /// Prints the receiver state to stream
 void   printOutputAt (FILE *file, TimeStep* tStep) ;

  /// Returns the last equilibrated scalar measure of the largest strain level
  double giveKappa () {return kappa;}
  /// Returns the temp. scalar measure of the largest strain level
  double giveTempKappa () {return tempKappa;}
  /// Sets the temp scalar measure of the largest strain level to given value
  void   setTempKappa (double newKappa) { tempKappa = newKappa;}
  /// Returns the last equilibrated damage level
  double giveDamage () {return damage;}
  /// Returns the temp. damage level
  double giveTempDamage () {return tempDamage;}
  /// Sets the temp damage level to given value
  void   setTempDamage (double newDamage) { tempDamage = newDamage;}
  /// Returns the last equilibrated plastic shear strain component 1
  double giveGam1p () {return gam1p;}
  /// Returns the last equilibrated plastic shear strain component 2
  double giveGam2p () {return gam2p;}
  /// Sets the temp plastic shear strain component 1 to given value
  void   setTempGam1p (double gam) {tempGam1p = gam;}
  /// Sets the temp plastic shear strain component 2 to given value
  void   setTempGam2p (double gam) {tempGam2p = gam;}

 // definition
  const char* giveClassName () const { return "CohesiveInterfaceMaterialStatus" ;}
  classType giveClassID () const { return CohesiveInterfaceMaterialStatusClass; }
 
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
  @param mode determines ammount of info required in stream (state, definition,...)
  @param obj pointer to integration point, which invokes this method
  @return contextIOResultType.
  */
 contextIOResultType    saveContext (DataStream* stream, ContextMode mode, void *obj = NULL);
 /**
  Restores context of receiver from given stream. 
  @param stream stream where to read data
  @param mode determines ammount of info required in stream (state, definition,...)
  @param obj pointer to integration point, which invokes this method
  @return contextIOResultType.
  */
 contextIOResultType    restoreContext(DataStream* stream, ContextMode mode, void *obj = NULL);
 
}; 


//---------------------------------------------------------------------------------------------------
// c l a s s   CohesiveInterfaceMaterial
//---------------------------------------------------------------------------------------------------

/**
 Class representing cohesive interface material model.
*/
class CohesiveInterfaceMaterial : public StructuralMaterial
{
 protected:
  /// Coefficient of thermal dilatation
  double tempDillatCoeff;
  /// elastic properties (normal and shear moduli) 
  double kn, ks;
  /// parameters controling tensile softening
  double e0, ef, ksi;
  /// parameters controling shear (cohesion, friction coefficient)
  double coh, tanphi;
  /// parameters for the rate-dependent version
  double damchartime, damrateexp, plchartime, plrateexp;

 public:
  /// Constructor
  CohesiveInterfaceMaterial (int n,Domain* d);
  /// Destructor
  ~CohesiveInterfaceMaterial (){};
 
 /// Returns nonzero indicating that receiver is nonlinear
 int hasNonLinearBehaviour ()   {return 1;}
 /**
  Tests, if material supports material mode.
  @param mode required material mode
  @return nonzero if supported, zero otherwise
  */
 int hasMaterialModeCapability (MaterialMode mode){return(mode==_3dInterface);}
 const char* giveClassName () const {return "CohesiveInterfaceMaterial";}
 classType giveClassID ()     const {return CohesiveInterfaceMaterialClass;}

/**
 Computes full 3d material stiffness matrix at given integration point, time, respecting load history 
 in integration point.
 @param answer computed results
 @param form material response form
 @param mode material response mode
 @param gp integration point
 @param atTime time step (most models are able to respond only when atTime is current time step)
 */

 virtual void  give3dMaterialStiffnessMatrix (FloatMatrix& answer,
                        MatResponseForm, MatResponseMode, 
                        GaussPoint* gp,
                        TimeStep* atTime);

 double computeVolumetricStrain(GaussPoint* gp, TimeStep* atTime);

/**
 Computes the real stress vector for given strain increment and integration point. 
 Temp vars are updated accordingly
 @param answer contains result
 @param form material response form
 @param gp integration point 
 @param reducedStrain strain  vector in reduced form
 @param tStep current time step (most models are able to respond only when atTime is current time step)
 */
void giveRealStressVector (FloatArray& answer,  MatResponseForm, GaussPoint*, 
               const FloatArray&,TimeStep* );
  /**
     Computes the equivalent strain measure from given strain vector (full form).
     @param strain total strain vector in full form
  */
  virtual double computeEquivalentStrain (const FloatArray& strain);


  /**
  computes the value of damage parameter omega, based on given value of equivalent strain
  @param kappa equivalent strain measure
  */
  virtual double computeDamage(double kappa);
 
  virtual void  giveCharacteristicMatrix (FloatMatrix& answer,
                                          MatResponseForm form,
                                          MatResponseMode mode,
                                          GaussPoint* gp,
                                          TimeStep* atTime);

  virtual int     giveStressStrainComponentIndOf (MatResponseForm,MaterialMode mmode, int);
  virtual void    giveStressStrainMask (IntArray& answer, MatResponseForm, MaterialMode mmode) const;
  virtual int  giveSizeOfReducedStressStrainVector (MaterialMode);
  void giveReducedCharacteristicVector (FloatArray& answer, GaussPoint*, 
                                        const FloatArray& charVector3d);
  void giveFullCharacteristicVector (FloatArray& answer,  GaussPoint*, 
                                     const FloatArray&) ;
 
#ifdef __OOFEG
#endif
 /**
  Returns the integration point corresponding value in Reduced form.
  @param answer contain corresponding ip value, zero sized if not available
  @param aGaussPoint integration point
  @param type determines the type of internal variable
  @param type determines the type of internal variable
  @returns nonzero if ok, zero if var not supported
  */
 virtual int giveIPValue (FloatArray& answer, GaussPoint* aGaussPoint, InternalStateType type, TimeStep* atTime) ;
 /**
  Returns the mask of reduced indexes of Internal Variable component .
  @param answer mask of Full VectorSize, with components beeing the indexes to reduced form vectors.
  @param type determines the internal variable requested (physical meaning)
  @returns nonzero if ok or error is generated for unknown mat mode.
  */
 virtual int giveIntVarCompFullIndx (IntArray& answer, InternalStateType type, MaterialMode mmode);


 /**
  Returns the type of internal variable (scalar, vector, tensor,...).
  @param type determines the type of internal variable
  @returns type of internal variable
  */
 virtual  InternalStateValueType giveIPValueType (InternalStateType type) ;
 /**
  Returns the corresponding integration point  value size in Reduced form.
  @param type determines the type of internal variable
  @returns var size, zero if var not supported
  */
 virtual int giveIPValueSize (InternalStateType type, GaussPoint* aGaussPoint) ;
/**
 Returns a vector of coefficients of thermal dilatation in direction
 of each material principal (local) axis.
 @param answer vector of thermal dilatation coefficients
 @param gp integration point
 @param tStep time step (most models are able to respond only when atTime is current time step)
 */
  virtual void giveThermalDilatationVector (FloatArray& answer, GaussPoint*, TimeStep*) ;

 /**
 Instanciates the receiver from input record
 */
  IRResultType initializeFrom (InputRecord* ir);
  /** Setups the input record string of receiver
  @param str string to be filled by input record
  @param keyword print record keyword (default true)
  */
  virtual int giveInputRecordString(std::string &str, bool keyword = true);

  /// Creates corresponding material status 
  MaterialStatus* CreateStatus (GaussPoint *gp) const {return new CohesiveInterfaceMaterialStatus (1,FEMComponent::domain,gp);}

protected: 

  // Overloaded to use specialized versions of these services possibly implemented by linearElastic member


 void give3dInterfaceMaterialStiffnessMatrix (FloatMatrix& answer, MatResponseForm form, MatResponseMode rMode,
                                              GaussPoint* gp, TimeStep* atTime);

 // Auxiliary functions used by the rate-dependent version
 double computeDamageOverstress(double eps, double& damstrain, double omega, double dt);
 double solveBeta(double c, double N);
 double computeViscoplasticScalingFactor(double tauTrial, double tauYield, double dt);
} ;

} // namespace oofem
#endif

