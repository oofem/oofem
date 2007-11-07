/* $Header: /home/cvs/bp/oofem/sm/src/mdm.h,v 1.12 2003/05/19 13:04:00 bp Exp $ */
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

//   ****************************************
//   *** CLASS MICROPLANE DAMAGE MATERIAL ***
//   ****************************************

#ifndef mdm_h 

/**
   Select the mapping algorithm. The IDM_USE_MMAShapeFunctProjection does not work, since
   this mapper does not preserve the max. property of damage and equivalent strain.
*/
#define MDM_USE_MMAClosestIPTransfer
//#define MDM_USE_MMAShapeFunctProjection
//#define MDM_USE_MMALeastSquareProjection


// Allows to select different mapping algorithms
#define MDM_MAPPING_DEBUG 1

#include "microplanematerial.h"
#include "structuralnonlocalmaterialext.h"
#include "gausspnt.h"
#include "matconst.h"
#include "structuralms.h"
#include "materialmapperinterface.h"
#include "mmaclosestiptransfer.h"

#ifdef MDM_MAPPING_DEBUG
#include "mmashapefunctprojection.h"
#include "mmaleastsquareprojection.h"

#else

#ifdef MDM_USE_MMAShapeFunctProjection
#include "mmashapefunctprojection.h"
#endif
#ifdef MDM_USE_MMALeastSquareProjection
#include "mmaleastsquareprojection.h"
#endif
#endif

#include "mmashapefunctprojection.h"
#include "mmaleastsquareprojection.h"

class Microplane;

/**
   Material status class MDMStatus associated to MDM matarial
*/
class MDMStatus : public StructuralMaterialStatus, public StructuralNonlocalMaterialStatusExtensionInterface
{

 protected:
  // add history variables declaration here
  FloatArray Psi, PsiTemp;  // damage values on individual microplanes
  FloatMatrix DamageTensor, DamageTensorTemp, localDamageTensor;  // macroscopic damage tensor
  //EigenSystem* DamageSystem;  // principal damage directions
  FloatArray tempDamageTensorEigenValues, damageTensorEigenValues;
  FloatMatrix tempDamageTensorEigenVectors, damageTensorEigenVectors;

 public:
  MDMStatus (int n, int nsd, int nmplanes, Domain*d, GaussPoint* g) ;
  ~MDMStatus ();

  // add declaration of access and update functions for history variables
  void setTempDamageTensorEigenVals (const FloatArray& src) {tempDamageTensorEigenValues = src;}
  void setTempDamageTensorEigenVec  (const FloatMatrix& src) {tempDamageTensorEigenVectors = src;}
  const FloatArray* giveTempDamageTensorEigenVals () {return &tempDamageTensorEigenValues;}
  const FloatArray* giveDamageTensorEigenVals () {return &damageTensorEigenValues;}
  const FloatMatrix* giveTempDamageTensorEigenVec () {return &tempDamageTensorEigenVectors;}
  const FloatMatrix* giveDamageTensorEigenVec () {return &damageTensorEigenVectors;}



  double giveMicroplaneTempDamage (int m) {return PsiTemp.at(m);}
  double giveMicroplaneDamage (int m) {return Psi.at(m);}
  void   setMicroplaneTempDamage (int m, double val) {PsiTemp.at(m) = val;}
  void   giveMicroplaneDamageValues (FloatArray& answer) {answer=Psi;}
  void   setMicroplaneTempDamageValues  (FloatArray& src) {PsiTemp = src;}

  void   giveTempDamageTensor (FloatMatrix& answer) {answer = DamageTensorTemp;}
  void   giveDamageTensor (FloatMatrix& answer) {answer = DamageTensor;}
  void   setTempDamageTensor (FloatMatrix& src) {DamageTensorTemp = src;}
  void   setLocalDamageTensorForAverage (FloatMatrix& src) {localDamageTensor = src;}
  void   giveLocalDamageTensorForAverage (FloatMatrix& answer) {answer = localDamageTensor;}
  const FloatMatrix* giveLocalDamageTensorForAveragePtr () {return &localDamageTensor;}


  // definition
  /// Returns class name of receiver
  const char* giveClassName () const { return "MDMStatus" ;}
  /// Returns class ID of receiver
  classType             giveClassID () const { return MicroplaneDamageMaterialStatusClass; }

  void   printOutputAt (FILE *file, TimeStep* tStep) ;
  /**
     Initializes temporary history variables of receiver to previous equilibrium values.
  */
  virtual void initTempStatus ();
  /**
     Updates history variables after new equilibrium state has been reached.
     (e.g. temporary variables are copied into equilibrium values
  */
  virtual void updateYourself(TimeStep*);

  /// Stores receiver's state to stream
  contextIOResultType    saveContext (DataStream* stream, ContextMode mode, void *obj = NULL);
  /// Restores receiver's state from stream
  contextIOResultType    restoreContext(DataStream* stream, ContextMode mode, void *obj = NULL);
  /**
     Interface requesting service.  In the case of nonlocal constitutive models, 
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
   Impolementation of microplane damage material (According to M.Jirasek). 
*/
class MDM : public MicroplaneMaterial, public StructuralNonlocalMaterialExtensionInterface, public MaterialModelMapperInterface
{

 public:
  enum MDMFormulatrionType {COMPLIANCE_DAMAGE, STIFFNESS_DAMAGE};
  enum MDMModeType {mdm_3d, mdm_2d};

 protected:
  int ndc; // number of damage components
  int nsd; // number of spatial dimensions
  int type_dam;
  int type_soft;

  /// parameter controlling the elastic limit
  double mdm_Ep;
  /// prescribed value of ef-ep 
  double mdm_Efp;  
  double ParMd; // (m/E*Ep)
  double tempDillatCoeff; // temperature dillatation coeff

  /// Fracture energy (necessary to determine Ep and Efp if not not given)
  double Gf;  
  /// Macroscopic tensile strength (necessary to determine Ep and Efp if not not given)
  double Ft; 

  MDMFormulatrionType formulation;
  MDMModeType mdmMode;

  /// Reference to bulk (undamaged) material
  StructuralMaterial *linearElasticMaterial; 

  /// flag indicating local ro nonlocal mode
  int nonlocal;
  /// Interaction radius, related to the nonlocal characteristic length of material.
  double R;

#ifdef MDM_MAPPING_DEBUG
  /// Mapper used to map internal variables in adaptivity
  static MMAShapeFunctProjection mapperSFT;
  static MMALeastSquareProjection mapperLST;

  // determines the mapping algorithm to be used
  enum MDMMapperType {mdm_cpt=0, mdm_sft=1, mdm_lst=2};
  MDMMapperType mapperType;
#else

#ifdef MDM_USE_MMAClosestIPTransfer
  /// Mapper used to map internal variables in adaptivity
  static MMAClosestIPTransfer mapper;
#endif
#ifdef MDM_USE_MMAShapeFunctProjection
  /// Mapper used to map internal variables in adaptivity
  static MMAShapeFunctProjection mapper;
#endif
#ifdef MDM_USE_MMALeastSquareProjection
  /// Mapper used to map internal variables in adaptivity
  static MMALeastSquareProjection mapper;
#endif
#endif

  /// Mapper used to map stresses in adaptivity
  static MMAClosestIPTransfer mapper2;


 public: 
 
  /**
     Constructor. Creates Microplane Material belonging to domain d, with number n.
     @param n material number
     @param d domain to which newly created material belongs
  */
  MDM (int n,Domain* d) : MicroplaneMaterial (n,d), StructuralNonlocalMaterialExtensionInterface(d), MaterialModelMapperInterface()
    {linearElasticMaterial = NULL; nonlocal = 0; type_dam = 0; type_soft = 0; mdm_Ep = mdm_Efp = -1.0;}
  /// Destructor.
  ~MDM ()                {if (linearElasticMaterial) delete linearElasticMaterial;}

  /**
     Tests, if material supports material mode.
     @param mode required material mode
     @return nonzero if supported, zero otherwise
  */
  int hasMaterialModeCapability (MaterialMode mode);
  /**
     Computes real macro stress in corresponding macro integration point for
     given macroscopic strain  vector.
  */
  virtual void giveRealStressVector (FloatArray& answer, MatResponseForm, GaussPoint*, 
				     const FloatArray& , TimeStep*) ;

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

  /// Instanciates receiver from input record.
  IRResultType initializeFrom (InputRecord* ir);

  // identification and auxiliary functions
  /// Returns class name of the receiver.
  const char*    giveClassName () const { return "MDM (Microplane Damaga Material)" ;}
  /// Returns classType id of receiver.
  classType giveClassID ()         const {return MicroplaneDamageMaterialClass;}

  /**
     Computes real stress vector on given microplane
     (the meaning of  values depends on particular implementaion,
     e.g, can contain volumetric, devatoric normal srtresses and shear streses on microplane)
     for given increment of microplane strains.
     @param answer computed result
     @param mplane pointer to microplane object, for which response is computed
     @param strain strain vector 
     @param tStep time step
  */
  virtual void giveRealMicroplaneStressVector (FloatArray& answer, Microplane* mplane, 
					       const FloatArray& strain, TimeStep* tStep) {};

  /** 
      Initializes internal data (integration weights, 
      microplane normals and computes projection tensors)
      @param numberOfMicroplanes number of required microplanes
  */
  virtual void initializeData (int numberOfMicroplanes);
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
  virtual void giveSupportRadius (double& radius) {radius = this-> R;}

  /** Interface requesting service */
  virtual Interface* giveInterface (InterfaceType);

  /**
     @name The interface required by MaterialModelMapperInterface
  */
  //@{
  /** Maps the required internal state variables from 
      old mesh oldd to given ip. The result is stored in gp status.
      @param gp Integration point belonging to new domain which values will be mapped
      @param oldd old mesh reference
      @param tStep time step
      @return nonzero if o.k.
  */
  virtual int MMI_map (GaussPoint* gp, Domain* oldd, TimeStep* tStep);
  /** Updates the required internal state variables from previously mapped values. 
      The result is stored in gp status. This map and update splitting is necessary,
      for example for nonlocal models tahe local quantity to be averaged must be mapped in all ips
      and then update can happen, because it may depend on nonlocal variable, which is computed
      from local values.
      @param gp Integration point belonging to new domain which values will be mapped
      @param tStep time step
      @return nonzero if o.k.
  */
  virtual int MMI_update (GaussPoint* gp, TimeStep* tStep, FloatArray* estrain = NULL);
  /**
     Finishes the mapping for given time step. Used to perform cleanup. 
     Typically some mappers reguire to compute some global mesh data related to
     current step, which are valid for example to all IPs - so they are computed only once for
     all IPs, stored and they need to be dealocated. These mappers are typically class variables,
     but their finish is invoked by all members.
     @return nonzero if ok.
  */
  virtual int MMI_finish (TimeStep* tStep) ;
  //@}

#ifdef __PARALLEL_MODE
  /**
  Updates domain before nonloc average (using updateDomainBeforeNonlocAverage service)
  to ensure, that the tempDamageTensor status variable is correctly updated and
  pack this into given buffer to be send to remote partition and receiver's remote conterpart.
  Packs only when in nonlocal mode.
  @see Material::packUnknowns for description.
  @param buff communication buffer
  @param stepN solution step
  @param ip integration point
  */
  int packUnknowns (CommunicationBuffer& buff, TimeStep* stepN, GaussPoint* ip);
  /**
  Unpack tempDamageTensor value from given buffer.  Unpacks only when in nonlocal mode.
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
 /**
     Returns the weight representing relative computational cost of receiver
     The reference material model  is linear isotropic material - its weight is set to 100. 
     The other material models should compare to this reference model.
  */
 virtual int predictRelativeComputationalCost (GaussPoint* gp);
  /**
     Returns the relative redistribution cost of the receiver
  */
  virtual int predictRelativeRedistributionCost (GaussPoint* gp) {return 100;}

#endif

  virtual MaterialStatus* CreateStatus (GaussPoint* gp) const;

 protected:
  /// Returns reference to undamaged (bulk) material
  StructuralMaterial* giveLinearElasticMaterial () {return linearElasticMaterial;}

  virtual MaterialStatus* CreateMicroplaneStatus (GaussPoint* gp) 
    {return NULL;}
  void computeDamageTensor(FloatMatrix& tempDamageTensor, const FloatArray& totalStrain, 
			   GaussPoint* gp, TimeStep* atTime);
  void computeLocalDamageTensor(FloatMatrix& tempDamageTensor, const FloatArray& totalStrain, 
				GaussPoint* gp, TimeStep* atTime);
  double computeDamageOnPlane (GaussPoint* gp, Microplane* mplane, const FloatArray& strain);
  void   computePDC (FloatMatrix& tempDamageTensor, FloatArray& tempDamageTensorEigenVals, 
		     FloatMatrix& tempDamageTensorEigenVec);
  void  transformStrainToPDC(FloatArray& answer, FloatArray& strain, 
			     FloatMatrix& t, GaussPoint* gp);
  void applyDamageTranformation(FloatArray& strainPDC, const FloatArray& tempDamageTensorEigenVals);
  void transformStressFromPDC(FloatArray& answer, const FloatArray& stressPDC, const FloatMatrix& t, GaussPoint* gp);
  void computeEffectiveStress (FloatArray& stressPDC, const FloatArray& strainPDC, 
			       GaussPoint* gp, TimeStep* atTime);
  void giveMaterialStiffnessMatrix (FloatMatrix& answer, MatResponseForm form,
				    MatResponseMode mode, GaussPoint* gp,
				    TimeStep* atTime);
  void applyDamageToStiffness(FloatMatrix& d, GaussPoint* gp);
  void transformStiffnessfromPDC (FloatMatrix& de, const FloatMatrix& t);
  /**
     Method for computing plane stress stifness matrix of receiver.
     Default implementation overloaded to use direct implementation of 
     corresponding service at bulk material model level.
     @param answer stifness matrix
     @param form material response form
     @param mode material response mode
     @param gp integration point, which load history is used
     @param atTime time step (most models are able to respond only when atTime is current time step)
  */
  void givePlaneStressStiffMtrx (FloatMatrix& answer, MatResponseForm form,MatResponseMode,
				 GaussPoint* gp,
				 TimeStep* atTime);
  /**
     Method for computing plane stress stifness matrix of receiver.
     Default implementation overloaded to use direct implementation of 
     corresponding service at bulk material model level.
     @param answer stifness matrix
     @param form material response form
     @param mode material response mode
     @param gp integration point, which load history is used
     @param atTime time step (most models are able to respond only when atTime is current time step)
  */
  void givePlaneStrainStiffMtrx (FloatMatrix& answer,
				 MatResponseForm,MatResponseMode,GaussPoint* gp,
				 TimeStep* atTime);

  void  rotateTensor4(FloatMatrix& Dlocal, const FloatMatrix& t);
  void  formTransformationMatrix(FloatMatrix& answer, const FloatMatrix& t, int n);


  //double giveParameterEfp(const FloatArray& reducedStrain, GaussPoint* gp);
  void   giveRawMDMParameters (double& Efp, double& Ep, const FloatArray& reducedStrain,GaussPoint* gp);
} ;


#define mdm_h
#endif
