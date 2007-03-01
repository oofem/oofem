/* $Header: /home/cvs/bp/oofem/sm/src/Attic/mplasticmaterial2.h,v 1.1.2.1 2004/04/05 15:19:47 bp Exp $ */
/*

                   *****    *****   ******  ******  ***   ***                            
                 **   **  **   **  **      **      ** *** **                             
                **   **  **   **  ****    ****    **  *  **                              
               **   **  **   **  **      **      **     **                               
              **   **  **   **  **      **      **     **                                
              *****    *****   **      ******  **     **         
            
                                                                   
               OOFEM : Object Oriented Finite Element Code                 
                    
                 Copyright (C) 1993 - 2003   Borek Patzak                                       



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

#ifndef plasticmaterial2_h 

#include "structuralmaterial.h"
#include "linearelasticmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"
#include "structuralms.h"
#ifndef __MAKEDEPEND
#include <vector>
#include <set>
#endif

class GaussPoint ;

/**
   This class implements associated Material Status to MPlasticMaterial.
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
   
*/
class MPlasticMaterial2Status : public StructuralMaterialStatus
{
 public:
  enum state_flag_values {PM_Elastic, PM_Yielding, PM_Unloading};
 
 protected:
  
  /// plastic strain vector
  FloatArray plasticStrainVector; 
  FloatArray tempPlasticStrainVector;
  
  /// strain space hardening variables
  FloatArray strainSpaceHardeningVarsVector;
  FloatArray tempStrainSpaceHardeningVarsVector;
  
  /// yield function status indicator
  int state_flag;
  int temp_state_flag;

  /// consistency parameter values (needed for algorithmic stiffness)
  FloatArray gamma, tempGamma;
  /// active set of yield functions (needed for algorithmic stiffness)
  IntArray activeConditionMap, tempActiveConditionMap;

 public:
  MPlasticMaterial2Status (int n, Domain*d, GaussPoint* g) ;
  ~MPlasticMaterial2Status ();
  
  void   printOutputAt (FILE *file, TimeStep* tStep) ;
  
  virtual void initTempStatus ();
  /// update state after new equilibrium state reached
  virtual void updateYourself(TimeStep*);
  
  /// saves current context(state) into stream
  contextIOResultType    saveContext (FILE* stream, void *obj = NULL);
  /// restores state from stream
  contextIOResultType    restoreContext(FILE* stream, void *obj = NULL);
  
  /// returns the equilibrated strain vector
  void givePlasticStrainVector (FloatArray& answer) const {answer = plasticStrainVector;}
  /// returns the actual (temp) strain vector
  void giveTempPlasticStrainVector (FloatArray& answer) const {answer = tempPlasticStrainVector;}
  /// returns the equilibrated hardening variable vector
  void giveStrainSpaceHardeningVars (FloatArray& answer) const 
    {answer = strainSpaceHardeningVarsVector;}
  /// returns the actual (temp) hardening variable vector
  void giveTempStrainSpaceHardeningVarsVector (FloatArray& answer) const 
    {answer = tempStrainSpaceHardeningVarsVector;}
  
  void         letPlasticStrainVectorBe (const FloatArray& v) 
    { plasticStrainVector = v ;}
  void         letTempPlasticStrainVectorBe (const FloatArray& v) 
    { tempPlasticStrainVector = v ;}
  void         letTempStrainSpaceHardeningVarsVectorBe (const FloatArray& v) 
    { tempStrainSpaceHardeningVarsVector = v ;}
  void         letStrainSpaceHardeningVarsVectorBe (const FloatArray& v) 
    { strainSpaceHardeningVarsVector = v ;}
  
  int giveStateFlag () {return state_flag;}
  int giveTempStateFlag () {return temp_state_flag;}
  void letTempStateFlagBe (int v) {temp_state_flag = v;}

  void  giveTempActiveConditionMap (IntArray& answer) {answer=tempActiveConditionMap;}
  void  setTempActiveConditionMap (const IntArray& v) {tempActiveConditionMap=v;}
  void  giveTempGamma (FloatArray& answer) {answer= tempGamma;}
  void  setTempGamma  (const FloatArray& v) {tempGamma=v;}

  // definition
  const char* giveClassName () const { return "MPlasticMaterial2Status" ;}
  classType             giveClassID () const
    { return MPlasticMaterialStatusClass; }
  
}; 

/**
   This class represents a base class for non-associated multisurface plasticity.
   The Multisurface plasticity is characterized by the following:
   Let \f$\sigma, \varepsilon\f$, and \f$\varepsilon^p\f$ be the stress, total strain, and plastic strain vectors, respectively.
   It is assumed that the total strain is decomposed into reversible elastic and irreversible plastic parts
   $\varepsilon=\varepsilon^e+\varepsilon^p$.
   The elastic response is characterized in terms of elastic constitutive matrix $D^e$ as $\sigma=D^e(\varepsilon-\varepsilon^e)$
   As long as the stress remains inside the elastic domain, the deformation process is purely elastic and the 
   plastic strain does not change.

   It is assumed that the elastic domain, denoted as $IE$ is bounded by a composite yield surface. It is defined as
   \f[
   IE=\{(\sigma,\kappa)|f_i(\sigma,\kappa)<0, \rm{for\ all\ }i\in\{1,\cdots,m\}\}
   \f]
   where \f$f_i(\sigma,\kappa)\f$ are \f$m\ge1\f$ yield functions intersecting in a possibly non-smooth fashion. The 
   vector \f$\kappa\f$ contains internal variables controlling the evolution of yield surfaces (amount of hardening or softening). 
   The evolution of plastic strain $\ep$ is expressed in Koiter's form. Assuming the non-associated plasticity, this reads
   \f[
   \label{epe}
   \varepsilon^p=\sum^{m}_{i=1} \lambda^i \partial_{\sigma}g_i(\sigma,\kappa)
   \f]
   where \f$g_i\f$ are plastic potential functions. The \f$\lambda^i\f$ are referred as plastic consistency parameters, which satisfy the following Kuhn-Tucker conditions
   \f[
   \label{ktc}
   \lambda^i\ge0,\;f_i\le0,\;{\rm and}\ \lambda^i f_i=0
   \f]
   These conditions imply  that in the elastic regime the yield function must remain negative and the rate of the plastic multiplier is zero (plastic strain remains constant) while in the plastic regime the yield function must be equal to zero (stress remains on the surface) and the rate of the plastic multiplier is positive.
   The evolution of vector of internal hardening/softening variables \f$\kappa\f$  is expressed in terms of a general 
   hardening/softening law of the form
   \f[
   \dot{\kappa} = \dot{\kappa}(\sigma, \lambda)
   \f]
   where \f$\lambda\f$ is the vector of plastic consistency parameters \f$\lambda_i\f$.
   
*/
class MPlasticMaterial2 : public StructuralMaterial
{
protected:
  // add common (same for every gauss point) material parameters here
  
  /// Reference to bulk (undamaged) material
  LinearElasticMaterial *linearElasticMaterial; 
  /// number of yield surfaces
  int nsurf;
  /// protected type to determine the return mapping algorithm
  enum ReturnMappingAlgoType {mpm_ClosestPoint, mpm_CuttingPlane} rmType;
  /// type that allows to distinguish between yield function and loading function
  enum functType {yieldFunction, loadFunction};
  enum plastType {associatedPT, nonassociatedPT} plType;
  /// flag indicating whether iterative update of a set of active yield conditions takes place
  bool iterativeUpdateOfActiveConds;
  /// set for keeping record of generated populations of active yield condiditons during return
  std::set<long> populationSet;

 public:
 
  MPlasticMaterial2 (int n,Domain* d) ;
  ~MPlasticMaterial2 () ;
 
  // identification and auxiliary functions
  int hasNonLinearBehaviour ()   { return 1 ;}
  int hasMaterialModeCapability (MaterialMode mode);
  const char*    giveClassName () const { return "MPlasticMaterial2" ;}
  classType giveClassID ()         const { return MPlasticMaterialClass;}

  /// Returns reference to undamaged (bulk) material
  LinearElasticMaterial* giveLinearElasticMaterial () {return linearElasticMaterial;}

  /** 
      Returns true if stiffness matrix of receiver is symmetric
      Default implementation returns true.
  */
  bool isCharacteristicMtrxSymmetric (MatResponseMode rMode) {return true;}

  // non-standart - returns time independent material constant
  // double   give (int) ; 

  virtual void give3dMaterialStiffnessMatrix (FloatMatrix& answer, 
                                              MatResponseForm, MatResponseMode,
                                              GaussPoint* gp,
                                              TimeStep* atTime);

  
  void giveRealStressVector (FloatArray& answer,  MatResponseForm, GaussPoint*, 
                             const FloatArray&, TimeStep* );

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

#ifdef __OOFEG
#endif

  // auxiliary functions
 virtual int         giveSizeOfFullHardeningVarsVector()  {return 0;}
 virtual int         giveSizeOfReducedHardeningVarsVector(GaussPoint*)  {return 0;}

 MaterialStatus* CreateStatus (GaussPoint* gp) const ;

protected: 
  
 virtual int giveMaxNumberOfActiveYieldConds (GaussPoint* gp)=0;
 void closestPointReturn (FloatArray& answer, IntArray& activeConditionMap, FloatArray& gamma,
                          MatResponseForm form, GaussPoint* gp, 
                          const FloatArray& totalStrain, FloatArray& plasticStrainR, 
                          FloatArray& strainSpaceHardeningVariables, TimeStep* atTime);
 
 void cuttingPlaneReturn (FloatArray& answer, IntArray& activeConditionMap, FloatArray& gamma,
                          MatResponseForm form, GaussPoint* gp, 
                          const FloatArray& totalStrain, FloatArray& plasticStrainR, 
                          FloatArray& strainSpaceHardeningVariables, TimeStep* atTime);
 
 // add here some auxiliary functions if needed
 /* void  computeGradientVector (FloatArray& answer, functType ftype, int isurf, GaussPoint* gp, const FloatArray& fullStressVector,
    const FloatArray& fullStressSpaceHardeningVars);*/
 void computeResidualVector(FloatArray& answer, GaussPoint *gp, const FloatArray& gamma, 
                            const IntArray& activeConditionMap, const FloatArray& plasticStrainVectorR, 
                            const FloatArray& strainSpaceHardeningVariables, std::vector<FloatArray>& gradVec);
 virtual void giveConsistentStiffnessMatrix (FloatMatrix& answer, 
                                             MatResponseForm, 
                                             MatResponseMode, 
                                             GaussPoint* gp,
                                             TimeStep* atTime);
 
 virtual void giveElastoPlasticStiffnessMatrix (FloatMatrix &answer,
                                                MatResponseForm form, 
                                                MatResponseMode mode, 
                                                GaussPoint* gp,
                                                TimeStep* atTime);

 void  computeAlgorithmicModuli (FloatMatrix& answer,
                                 GaussPoint *gp, const FloatMatrix &elasticModuliInverse,
                                 const FloatArray& gamma, const IntArray& activeConditionMap, 
                                 const FloatArray& fullStressVector,
                                 const FloatArray& strainSpaceHardeningVariables);

 /*void  computeDiagModuli(FloatMatrix& answer,
                         GaussPoint *gp, FloatMatrix &elasticModuliInverse,
                         FloatMatrix &hardeningModuliInverse);*/


  /// Computes the value of yield function 
 virtual double computeYieldValueAt (GaussPoint *gp, int isurf, const FloatArray& stressVector, 
                                     const FloatArray& strainSpaceHardeningVariables) = 0;

 /// Computes the stress gradient of yield/loading function (df/d_sigma)
 virtual void computeStressGradientVector (FloatArray& answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray& stressVector,
                                           const FloatArray& strainSpaceHardeningVariables) = 0;
 void computeReducedStressGradientVector(FloatArray& answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray& stressVector,
                                         const FloatArray& strainSpaceHardeningVariables);

 /**
    Computes the increment of strain-space hardening variables
    @param answer result
    @param stress updated stress (corresponbds to newly reached state)
    @param dlambda increment of consistency parameters
    @param dplasticStrain actual plastic strain increment
    @param activeConditionMas array of active yield conditions
  */
 virtual void computeStrainHardeningVarsIncrement (FloatArray& answer, GaussPoint* gp,
                                                   const FloatArray& stress, const FloatArray& dlambda,
                                                   const FloatArray& dplasticStrain, const IntArray& activeConditionMap) = 0;
  /**
     Computes the derivative of yield/loading function with respect to \f$\kappa\f$ vector
  */
 virtual void computeKGradientVector (FloatArray& answer, functType ftype, int isurf, GaussPoint *gp, FloatArray& fullStressVector,
                                      const FloatArray& strainSpaceHardeningVariables)=0;
 
  /**
     Computes derivative of $\kappa$ vector with respect to stress
   */
 virtual void computeReducedHardeningVarsSigmaGradient(FloatMatrix& answer, GaussPoint* gp, const IntArray& activeConditionMap,
                                                       const FloatArray&fullStressVector, 
                                                       const FloatArray& strainSpaceHardeningVars,
                                                       const FloatArray& gamma) =0;
 /// computes derivative of $\kappa$ vector with respect to lambda vector
 virtual void computeReducedHardeningVarsLamGradient(FloatMatrix& answer, GaussPoint* gp, int actSurf, 
                                                     const IntArray& activeConditionMap, 
                                                     const FloatArray&fullStressVector, 
                                                     const FloatArray& strainSpaceHardeningVars,
                                                     const FloatArray& gamma) =0;
 /**
    Indicates, whether receiver model has hardening/softening behaviour or behaves according to perfect plasticity theory
  */
 virtual int hasHardening () = 0;
 /* virtual void  computeReducedGradientMatrix (FloatMatrix& answer, int isurf,
                                             GaussPoint *gp,
                                             const FloatArray& stressVector, 
                                             const FloatArray& stressSpaceHardeningVars) = 0;*/
 /// Computes second derivative of loading finction with respect to stress.
 virtual void  computeReducedSSGradientMatrix (FloatMatrix& gradientMatrix,  int i, GaussPoint* gp , const FloatArray& fullStressVector, 
                                               const FloatArray& strainSpaceHardeningVariables)=0;
 /// Computes second derivative of loading finction with respect to stress and hardening vars 
 virtual void  computeReducedSKGradientMatrix (FloatMatrix& gradientMatrix,  int i, GaussPoint* gp , const FloatArray& fullStressVector, 
                                               const FloatArray& strainSpaceHardeningVariables)=0;

 /**
    Computes full-space trial stress increment (elastic)
    @param answer contains result in full space
    @param gp integratioan point
    @param strainIncrement strain increment vector
    @param atTime solution step
  */
 virtual void computeTrialStressIncrement (FloatArray& answer, GaussPoint *gp, 
                                           const FloatArray& strainIncrement, TimeStep* atTime);
 virtual void computeReducedElasticModuli(FloatMatrix& answer, GaussPoint *gp, 
                                          TimeStep *atTime);
 //virtual void compute3dElasticModuli(FloatMatrix& answer, GaussPoint *gp,
 //                                    TimeStep *atTime) = 0;
 
 // next functions overloaded rom structural material level
 virtual void givePlaneStressStiffMtrx (FloatMatrix& answer,
                                        MatResponseForm,MatResponseMode,
                                        GaussPoint* gp,
                                        TimeStep* atTime);
 virtual void givePlaneStrainStiffMtrx (FloatMatrix& answer,
                                        MatResponseForm,MatResponseMode,
                                        GaussPoint* gp, 
                                        TimeStep* atTime);
  virtual void give1dStressStiffMtrx (FloatMatrix& answer,
                                      MatResponseForm,MatResponseMode,
                                      GaussPoint* gp,
                                      TimeStep* atTime);
  virtual void give2dBeamLayerStiffMtrx (FloatMatrix& answer,
                                         MatResponseForm,MatResponseMode,
                                         GaussPoint* gp,
                                         TimeStep* atTime);
  virtual void give2dPlateLayerStiffMtrx( FloatMatrix& answer,
                                          MatResponseForm,MatResponseMode,
                                          GaussPoint* gp,
                                          TimeStep* atTime);
  
  virtual void give1dFiberStiffMtrx(FloatMatrix& answer,
                                    MatResponseForm,MatResponseMode,GaussPoint* gp,
                                    TimeStep* atTime);
  
  virtual void give3dShellLayerStiffMtrx(FloatMatrix& answer,
                                         MatResponseForm,MatResponseMode,
                                         GaussPoint* gp,
                                         TimeStep* atTime);

 protected:
    long getPopulationSignature (IntArray &mask);
    int testPopulation (long pop);
    void clearPopulationSet ();
    void addNewPopulation (IntArray &mask);
    int getNewPopulation (IntArray& result, IntArray &candidateMask, int degree, int size);
};


#define plasticmaterial2_h
#endif

