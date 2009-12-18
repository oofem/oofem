/*
                                                                           
                                                                       
                         Author: Peter Grassl                    
                                  (C) 2003-2009                                   
                                                                             
                                   NOTICE                                    
                                                                             
   Permission to use, copy, modify, and distribute this software and         
   its documentation for any purpose and without fee is hereby granted       
   provided that the above copyright notice appear in all copies and         
   that both the copyright notice and this permission notice appear in       
   supporting documentation.                                                 
                                                                             
   Neither the Institution (Czech Technical University) nor the Authors 
	make any representations about the suitability of this software for
   any purpose.  This software is provided ``as is''without expressed 
	or implied warranty.                  
                                                                             

*/
// file : concretedpm.h

//   ************************************************************
//   *** CLASS CONCRETE DAMAGE-PLASTIC MATERIAL MODEL STATUS  ***
//   ************************************************************

#ifndef ConcreteDPM_h
#define ConcreteDPM_h

#include "structuralmaterial.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"
#include "structuralms.h"
#include "strainvector.h"
#include "stressvector.h"
#include "isolinearelasticmaterial.h"


namespace oofem {

class ConcreteDPMStatus : public StructuralMaterialStatus 
{

  /* 
     This class implements the material status associated to ConcreteDPM.

     At every GaussPoint for which this material is active, it is an attribute 
     of matStatusDictionary.
	 
     DESCRIPTION:
     Idea used there is that we have variables describing:
	 
     temp-variables: 
     at the beginning of an iteration for material equilibrium,
     these variables carry the previous equilibrium state. 
     When material equilibrium is achieved, 
     they carry the values associated to this state. If it corresponds
     to a new equilibrium, i.e. if convergence in the force balance is satisfied, 
     they are copied to the

     non-temp-variables: 
     they always carry the state corresponding to the latest global equilibrium.
  */



public:


  /// values of history variable state_flag
  enum state_flag_values {ConcreteDPM_Elastic, ConcreteDPM_Unloading, ConcreteDPM_Plastic, ConcreteDPM_Damage, ConcreteDPM_PlasticDamage,  ConcreteDPM_VertexCompression, ConcreteDPM_VertexTension, ConcreteDPM_VertexCompressionDamage, ConcreteDPM_VertexTensionDamage};
protected:
  /// history variables of the plasticity model
  //plastic strain
  StrainVector plasticStrain;
  StrainVector tempPlasticStrain;
  
  double tempVolumetricPlasticStrain;

  double dFDKappa;
  double deltaLambda;

  /// hardening variable
  double kappaP ;
  double tempKappaP ;


  double le;

  /// history variables of the damage model
  
  double equivStrain;
  double tempEquivStrain;

  double kappaD ;
  double tempKappaD ;
 
  double damage ;
  double tempDamage ;

  double deltaEquivStrain;

  /// Indicates the state (i.e. elastic, unloading, plastic, damage, vertex) of the Gauss point
  int state_flag ;
  int temp_state_flag ;

public:
  /// constructor
  ConcreteDPMStatus (int, Domain*, GaussPoint*) ;

  /// destructor
  ~ConcreteDPMStatus () ;

  void  initTempStatus () ;
  void  updateYourself (TimeStep*) ;
  void  printOutputAt (FILE *file, TimeStep* tStep) ;
  contextIOResultType  saveContext (DataStream*, ContextMode mode, void *obj = NULL) ;
  contextIOResultType  restoreContext (DataStream*, ContextMode mode, void *obj = NULL) ;
  const char*  giveClassName () const 
  { return "ConcreteDPMStatus" ; }
  classType  giveClassID () const
  { return ConcreteDPMStatusClass ; }

  // Inline functions for access to state variables 
  // give:
  // Functions used to access the value of a internal variable.
  /**
     Get the full plastic strain vector from the material status.
     @param answer plastic strain vector.
  */
  void  giveFullPlasticStrainVector (StrainVector& answer) const 
  { 
    StrainVector plasticStrainVector(_Unknown);
    givePlasticStrain(plasticStrainVector);
    plasticStrainVector.convertToFullForm(answer) ;
  }
  /**
     Get the plastic strain deviator from the material status.
     @param answer plastic strain deviator.
  */
  void  givePlasticStrain(StrainVector& answer) const 
  { answer = plasticStrain; }


  /** 
      Get the deviatoric plastic strain norm from the material status.
      @returns deiatoric plasticStrainNorm
  */
  double giveDeviatoricPlasticStrainNorm()
  {
    StrainVector deviatoricPlasticStrain(gp->giveMaterialMode());
    double volumetricPlasticStrain;
    plasticStrain.computeDeviatoricVolumetricSplit(deviatoricPlasticStrain, 
                                                   volumetricPlasticStrain);
    return deviatoricPlasticStrain.computeStrainNorm(); 
  }
  

  double giveVolumetricPlasticStrain() const
  {
    return 1./3.*(plasticStrain(0) + plasticStrain(1) + plasticStrain(2));
  }
  

  /**
     Get the hardening variable of the plasticity model.
     @returns the hardening variable of the plasticity model
  */
  double  giveKappaP () const 
  { return kappaP ; }
  
  /**
     Get the hardening variable of the damage model from the 
     material status.
     @returns hardening variable kappaD
  */
  double  giveKappaD () const 
  { return kappaD ; }

  /**
     Get the hardening variable of the damage model from the 
     material status.
     @returns hardening variable kappaD
  */
  double  giveEquivStrain () const 
  { return equivStrain; }
  

  /**
     Get the damage variable of the damage model from the 
     material status.
     @returns damage variable damage
  */
  double  giveDamage () const 
  { return damage; }


  /**
     Get the state flag from the material status.
     @returns state flag (i.e. elastic, unloading, yielding, vertex case yielding)
  */
  int giveStateFlag () const 
  { return state_flag ; }


  // giveTemp:
  // Functions used to access the temp variables.
  /**
     Get the temp value of the full plastic strain vector from the material status.
     @param answer temp value of plastic strain vector.
  */
  void  giveTempPlasticStrain (StrainVector& answer) const 
  { answer = tempPlasticStrain;}



  /** 
      Get the temp value of the volumetric plastic strain in plane stress 
  */
  double  giveTempVolumetricPlasticStrain() const
  {return tempVolumetricPlasticStrain;}


  /**
     Get the temp value of the hardening variable of the plasticity model
     from the material status.
     @returns temp value ofhardening variable kappaP
  */
  double  giveTempKappaP () const 
  { return tempKappaP ; }

 /**
     Get the temp value of the hardening variable of the damage model 
     from the material status.
     @returns temp value ofhardening variable kappaD
  */
  double  giveTempKappaD () const 
  { return tempKappaD ; }

  /**
     Get the temp value of the hardening variable of the damage model 
     from the material status.
     @returns temp value of the damage variable damage
  */
  double  giveTempDamage () const 
  { return tempDamage; }



  /**
     Get the temp value of the hardening variable of the damage model 
     from the material status.
     @returns temp value of the damage variable damage
  */
  double  giveDeltaEquivStrain () const 
  { return deltaEquivStrain; }


  /**
     Get the temp value of the state flag from the material status.
     @returns the temp value of the state flag (i.e. elastic, unloading, 
     yielding, vertex case yielding)
  */
  int  giveTempStateFlag () const 
  { return temp_state_flag ; }

  // letTemp...be :
  // Functions used by the material to assign a new value to a temp variable.
  /**
     Assign the temp value of deviatoric plastic strain.
     @v new temp value of deviatoric plastic strain
  */
  void  letTempPlasticStrainBe (const StrainVector& v) 
  { tempPlasticStrain = v ; }



  void  letDeltaLambdaBe(const double v)
  {deltaLambda = v;}

  /** 
      Assign the temp value of the volumetric
      plastic strain in plane stress 
  */
  void  letTempVolumetricPlasticStrainBe(const double v)
  {tempVolumetricPlasticStrain = v;}

  /**
     Assign the temp value of the hardening variable of the plasticity model.
     @v new temp value of the hardening variable
  */
  void  letTempKappaPBe (const double v) 
  { tempKappaP = v ; }


  /**
     Assign the temp value of the hardening variable of the damage model.
     @v new temp value of the hardening variable
  */
  void  letTempKappaDBe (const double v) 
  { tempKappaD = v ; }

  /**
     Assign the temp value of the hardening variable of the damage model.
     @v new temp value of the hardening variable
  */
  void  letTempEquivStrainBe (const double v) 
  { tempEquivStrain = v ; }


  /**
     Assign the temp value of the damage variable of the damage model.
     @v new temp value of the damage variable
  */
  void  letTempDamageBe (const double v) 
  { tempDamage = v ; }

  /**
     Assign the temp value of the damage variable of the damage model.
     @v new temp value of the damage variable
  */
  void  letDeltaEquivStrainBe (const double v) 
  { deltaEquivStrain = v ; }
  
  /** 
      Gives the characteristic length
  */
  double giveLe ()  {return le;}
  
  /** 
      Sets the characteristic length
  */

  void   setLe (double ls) 
  {le = ls;}
  /**
     Assign the temp value of the state flag.
     @v new temp value of the state flag (i.e. elastic, unloading, yielding, 
     vertex case yielding)
  */
  void  letTempStateFlagBe (const int v) 
  { temp_state_flag = v ; }

};


//   ********************************************************************
//   *** CLASS CONCRETE PLASTICITY ISOTROPIC DAMAGE MATERIAL STATUS   ***
//   ********************************************************************

/**
This class contains the combination of a local plasticity model for concrete with a local isotropic damage model. The yield surface of the plasticity model is based on the extension of the Menetrey and Willam yield criterion. The flow rule is nonassociated. The evolution laws of the hardening variables depend on the stress state. The plasticity model describes only hardening and perfect plasticity. It is based on h eeffective stress. The damage parameter of the isotropic damage model is based on the total volumetric strain. An exponential softening law is implemented.
*/
class ConcreteDPM : public StructuralMaterial
{

public:
  /// values of history variable state_flag


protected:
  enum Concrete_VertexType{VT_Regular, VT_Tension, VT_Compression};
  Concrete_VertexType vertexType;
  
  /**parameters of the yield surface of the plasticity model. fc is the uniaxial compressive strength, ft the uniaxial tensile strength and ecc controls the out of roundness of the deviatoric section.*/
  double fc, ft, ecc;
  
  ///parameter of the ductilityMeasure of the plasticity model.
  double AHard;
  ///parameter of the ductilityMeasure of the plasticity model.
  double BHard;
  ///parameter of the ductilityMeasure of the plasticity model.
  double CHard;
 ///parameter of the ductilityMeasure of the plasticity model.
  double DHard;
  
  ///parameter of the ductilityMeasure of the damage model.
  double ASoft; 
  ///parameter of the ductilityMeasure of the damage model.
  double BSoft;

  ///parameter of the hardening law of the plasticity model.
  double yieldHardInitial;

  ///Control parameter for te volumetric plastic flow of the plastic potential.
  double dilationConst;
  
  /// plastic multiplier of the plasticity model
  double deltaLambda;
  
  /// the volumetric stress.
  double sig;
  /// the length of the deviatoric stress.*/
  double rho; 

  /// the lode angle of the trial stress.
  double thetaTrial;

  ///the friction parameter of the yield surface
  double m;

  /// the dilation parameter of the plastic potential
  double mQ;


  /// pointer for linear elastic material
  //  IsotropicLinearElasticMaterial *ILEMaterial;
  LinearElasticMaterial *linearElasticMaterial;  

  /// elastic Young's modulus
  double eM;
  /// elastic shear modulus
  double gM; 
  /// elastic bulk modulus
  double kM;
  /// elastic poisson's ration
  double nu;

  /// hardening variable of plasticity model
  double kappaP ;
  double tempKappaP ;

  /// hardening variable of damage model
  double kappaD ;
  double tempKappaD ;

  ///damage variable of damage model
  double damage ;
  double tempDamage ;

  /// controll parameter for the exponential softening law.
  double ef; 

  /// yield tolerance for the plasticity model.
  double yieldTol ;

  /// Maximum number of iterations for stress return.
  int newtonIter ;

  // Material mode for convinient access
  MaterialMode matMode ;

  ///  stress and deviatoric part of it 
  StressVector effectiveStress;

public:
  /// constructor
  ConcreteDPM(int n,Domain* d) ;
  /// destructor
  ~ConcreteDPM () ;
 IRResultType initializeFrom (InputRecord* ir);


  /**
     Tests, if material supports material mode.
     @param mode required material mode
     @return nonzero if supported, zero otherwise
  */
  int  hasMaterialModeCapability (MaterialMode mMode) ;

  const char* giveClassName () const 
  { return "ConcreteDPM" ;}
  classType giveClassID () const 
  { return ConcreteDPMClass;}

  /// overloaded from Material
  virtual ConcreteDPMStatus* giveStatus (GaussPoint* gp) const
  { return static_cast <ConcreteDPMStatus*> (this -> Material::giveStatus(gp)); }

  LinearElasticMaterial* giveLinearElasticMaterial () {return linearElasticMaterial;}
  

  /**
     Computes real macro stress in corresponding macro integration point for
     given macroscopic strain vector.
  */
  void  giveRealStressVector (FloatArray& answer,  
                              MatResponseForm form, 
                              GaussPoint* gp, 
                              const FloatArray& strainVector, 
                              TimeStep* atTime) ;


  /// Perform stress return of the plasticity model and compute history variables.
  /**
     @param gp Gauss point
     @param strain strain vector of this Gauss point
  */

  void  performPlasticityReturn(GaussPoint* gp, 
                                StrainVector& strain) ;

  /// Check for vertex case.
  /**
     Check if the trial stress state falls within the vertex region of the plasticity model at the apex of triaxial extension or triaxial compression.
     @returns true for vertex case and false if regular stress return can be used.
     @param answer is the volumetric apex stress.
     @param sig is the volumetric stress.
     @param tempKappa is the hardening variable
  */

  bool checkForVertexCase(double answer,
                          const double sig,
                          const double tempKappa) ;
  
  /// Perform regular stress return for plasticity case.
  /**
     Perform regular stress return for the plasticity model, i.e. if the trial stress state does not lie in the vertex region.
     @param stress stress vector which is computed
     @param gp Gauss point
  */
  void performRegularReturn(StressVector& stress,
                            GaussPoint* gp);
  

  /// Perform stress return for vertex case of the plasticity model.
  /**
     Perform stress return for vertex case of the plasticity model, i.e. if the trial stress state lies within the vertex region.
     @param strain strain vector of this Gauss point
     @param apexStress volumetric stress at the apex of the yield surface.
     @param gp Gauss point
  */
  void performVertexReturn(StressVector& stress,
                           double apexStress,
                           GaussPoint* gp) ;

  /// Compute yield value of the plasticity surface.
  /**
     Compute the yield value based on stress and hardening variable.
     @param sig volumetric stress
     @param rho length of the deviatoric stress
     @param theta lode angle of the stress state
     @param tempKappa hardening variable
     @returns yield value
  */
  double computeYieldValue(const double sig,
                           const double rho,
                           const double theta,
                           const double tempKappa) const ;

  /// compute the value of hardening function.
  /**
     Compute the value of the hardening function based on the hardening 
     variable.
     @param tempKappa hardening variable
     @returns value of the hardening function
  */
  double computeHardeningOne(const double tempKappa) const ;
  
  /// compute the derivative of the hardening function
  /** Compute the derivative of the hardening function based on the 
      hardening parameter.
      @param tempKappa hardening variable
      @returns the derivative of the hardening function */
  double computeHardeningOnePrime(const double tempKappa) const ;
  
  
  /// compute the derivative of the yield surface with respect 
  /// to the hardening variable 
  /** compute the derivative of the yield surface with respect to the hardening
      variable based on the stress state and the hardening variable
      @param sig volumetric stress
      @param rho deviatoric length
      @param tempKappa hardening variable
      @returns the derivative of the yield surface*/
  double computeDFDKappa(const double sig,
                         const double rho,
                         const double tempKappa);
  

  ///compute the derivative of kappa with respect to delta lambda is computed
  /** compute the derivative of kappa with respect of delta lambda based on the stress state and the hardening variable.
      @param sig volumetric stress
      @param rho length of the deviatoric stress
      @param tempKappa the hardening variable
      @returns derivative of kappa with respect to delta lambda
  */
  double computeDKappaDDeltaLambda(const double sig,
                                   const double rho,
                                   const double tempKappa);
  



  /// compute the ductility measure is computed
  /** compute the ductility measure based on the stress state
      @param sig volumetric stress
      @param rho length of the deviatoric strength
      @returns ductility measure*/
  virtual double computeDuctilityMeasure(const double sig,
                                         const double rho);


  /// compute the first derivative of the ductility measure with respect to the invariants sig and rho
  /** compute the first derivative of the ductility measure with respect to the invariants sig and rho based on the stress state and the hardening parameter*/
  void computeDDuctilityMeasureDInv(FloatArray& answer,
                                    const double sig,
                                    const double rho,
                                    const double tempKappa);
  
  /* This matrix is the core of the closest point projection and collects 
     the derivative of the flow rule and the hardening parameters */
  void computeAMatrix(FloatMatrix& answer, 
                      const double sig,
                      const double rho,
                      const double tempKappa);

  /* Here, the first derivative of the plastic potential with respect 
     to the invariants sig and rho are computed */
  void computeDGDInv(FloatArray& answer, 
                     const double sig,
                     const double rho,
                     const double tempKappa);

  /* This function computes the ratio of the volumtric and deviatoric component
     of the flow direction. It is used within the vertex return to check, 
     if the vertex return is admissible. */
  double computeRatioPotential(const double sig,
                               const double tempKappa);

  /* Here, the second derivative of the plastic potential with respect to the 
     invariants sig and rho are computed */
  void computeDDGDDInv(FloatMatrix& answer, 
                       const double sig,
                       const double rho,
                       const double tempKappa);
  
  /* Here, the mixed derivative of the plastic potential with respect 
     to the invariants and the hardening parameter are determined */
  void computeDDGDInvDKappa(FloatArray& answer, 
                            const double sig,
                            const double rho,
                            const double tempKappa);

  /* Computes the mixed derivative of the hardening parameter kappa with 
     respect to the plastic multiplier delta Lambda and the invariants sig 
     and rho. */
  void computeDDKappaDDeltaLambdaDInv(FloatArray& answer,
                                      const double sig,
                                      const double rho,
                                      const double tempKappa);

  /* Computes the derivative of the evolution law of the hardening parameter kappa with respect to the hardening variable kappa */
  double computeDDKappaDDeltaLambdaDKappa(const double sig,
                                          const double rho,
                                          const double tempKappa);
  


  /** 
     Computes the derivative of the yield surface with respect to the 
     invariants sig and rho. */
  void computeDFDInv(FloatArray& answer, 
                     const double sig,
                     const double rho,
                     const double tempKappa) const;

  /**
     compute tempKappa
  */
  double computeTempKappa(const double kappaInitial,
                          const double sigTrial,
                          const double rhoTrial,
                          const double sig);



  /// Perform stress return for damage case.
  /**
     Perform stress return for the damage model, i.e. if the trial stress state does not violate the plasticity surface.
     @param elastic strain
     @param equivStrain
     @return tempKappaD
     @param gp GaussPoint
  */
  double  computeDamage(const StrainVector& strain, GaussPoint* gp, TimeStep* atTime);



  /// Compute damage parameter
virtual  double computeDamageParam (double kappa, GaussPoint* gp);
    
  /// Compute equivalent strain value
virtual  void computeEquivalentStrain(double& kappaD, const StrainVector& elasticStrain, GaussPoint* gp, TimeStep* atTime);


  /// Compute the ductility measure for the damage model
  double computeDuctilityMeasureDamage(const StrainVector& strain, GaussPoint* gp);

  /**
     Initialize the characteristic length, if damage is not yet activated
  */
  void initDamaged (double kappa, 
                    const StrainVector& elasticStrain, 
                    GaussPoint* gp);


  /// Compute the trial coordinates
  void computeTrialCoordinates (const StressVector& stress);


  /// Assign state flag
  void assignStateFlag(GaussPoint* gp);

  ///Computes the derivative of rho with respect to the stress
  void computeDRhoDStress(FloatArray& answer,
                          const StressVector& stress) const;
  
  ///Computes the derivative of sig with respect to the stress
  void computeDSigDStress(FloatArray& answer) const;
  
  ///Computes the seconfd derivative of rho with the respect to the stress
  void computeDDRhoDDStress(FloatMatrix& answer,
                                const StressVector& stress) const;
  
  ///Computes the derivative of costheta with respect to the stress
  void computeDCosThetaDStress(FloatArray& answer,
                               const StressVector& stress) const;

  ///Compute the derivative of R with respect to costheta
  double computeDRDCosTheta(const double theta, const double ecc) const;

  ///give the the 3d material stiffness
  void give3dMaterialStiffnessMatrix(FloatMatrix&, MatResponseForm, 
                                     MatResponseMode, GaussPoint*, TimeStep*);


  virtual bool isCharacteristicMtrxSymmetric (MatResponseMode rMode) {return false;}
  
    /**
       Returns the value of a state variable at the specified integration point in reduced form.
       @param answer contains corresponding ip value, 
       zero sized if not available
       @param gp integration point
       @param type determines the type of internal variable
     @param atTime time step
     @returns nonzero if ok, zero if var not supported
  */
  int giveIPValue (FloatArray& answer, 
                   GaussPoint* gp, 
                   InternalStateType type, 
                   TimeStep* atTime) ;

  /**
     Returns the size of a state variable at the specified integration point in reduced form.
     @param type determines the type of internal variable
     @param gp Gauss point
     @returns var size, zero if var not supported
  */
  int giveIPValueSize (InternalStateType type, 
                       GaussPoint* gp) ;

  /**
     Returns the mask of reduced indexes of Internal Variable component.
     @param answer mask of Full VectorSize, with components being the indexes to reduced form vectors.
     @param type determines the internal variable requested (physical meaning)
     @param mmode material mode
     @returns nonzero if ok or error is generated for unknown mat mode.
  */
  int giveIntVarCompFullIndx (IntArray& answer, 
                              InternalStateType type, 
                              MaterialMode mmode) ;

  /**
     Returns the type of internal variable (scalar, vector, tensor,...).
     @param type determines the type of internal variable
     @returns type of internal variable
  */
  InternalStateValueType giveIPValueType (InternalStateType type) ;

protected:   
  MaterialStatus*  CreateStatus (GaussPoint* gp) const ;
} ;

} // end namespace oofem
#endif
