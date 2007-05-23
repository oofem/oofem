/* $Header: /home/cvs/bp/oofem/oofemlib/src/cltypes.h,v 1.45.4.1 2004/04/05 15:19:43 bp Exp $ */
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

//
// FILE: cltypes.h
//

#ifndef cltypes_h

/// Maximum input line length read.
#define OOFEM_MAX_LINE_LENGTH 32768
// size of token buffer, which holds tokens
#define MAX_TOKENS_LENGTH 32768
/// Maximum keyword name string length.
#define MAX_NAME_LENGTH 30
/// Maximum file name path string length. 
#define MAX_FILENAME_LENGTH 120
/// Maximum class name string length
#define MAX_CLASSNAME_LENGTH 50
/// max number of tokens 
#define MAX_TOKENS 8000
/// max legth of error message
#define MAX_ERROR_MSG_LENGTH 2048



#define ERR_HEADER "_______________________________________________________\nError : \n"
#define WRN_HEADER "_______________________________________________________\nWarning : \n"

#define ERR_TAIL   "_______________________________________________________\a\n"

/**
 Type introduced to distinguish between classes. Intended for run-time
 type checking. Every derived class from FemComponent base class should 
 owerload the abstract giveClassID service. This service can be used
 to distinguish between multiple derived classes. 
 */
enum classType {
  
  FEMComponentClass,

 GeneralBoundaryConditionClass,
  LoadClass,
  BoundaryConditionClass,
  InitialConditionClass, 
  NodalLoadClass,
  DeadWeightClass,
  BodyLoadClass,
  StructuralTemperatureLoadClass,
  BoundaryLoadClass,
  LinearEdgeLoadClass,
  ConstantEdgeLoadClass,
  ConstantSurfaceLoadClass,
  PointLoadClass,
  UserDefinedTemperatureFieldClass,
  TransmissionBCClass,
  ConvectionBCClass,

  LoadTimeFunctionClass,
  ConstantFunctionClass,
  PeakFunctionClass,
  PiecewiceClass,
  PeriodicPiecewiseClass,
  HeavisideLTFClass,
 UserDefinedLoadTimeFunctionClass,

  DofClass,
  MasterDofClass,
  SimpleSlaveDofClass,
  SlaveDofClass,
  DofManagerClass,
  NodeClass,
  RigidArmNodeClass,
  HangingNodeClass,
  ElementSideClass,
  RemoteMasterDofClass,
  SharedMasterDofClass,
  NullDofClass,
  RigidArmSlaveDofClass,
  HangingDofClass,

  ElementClass,
  StructuralElementClass,
  TransportElementClass,
  FMElementClass,
  CBSElementClass,
  SUPGElementClass,

  TrStructuralElementClass,
  NLStructuralElementClass,
  PlaneStress2dClass,
  Quad1PlaneStrainClass,
  QPlaneStress2dClass,
  TrPlaneStress2dClass,
  CCTPlateClass,
  LSpaceClass,
  QSpaceClass,
  LTRSpaceClass,
  Truss2dClass,
  LIBeam2dClass,
  TrPlaneStrRotClass,
  Axisymm3dClass,
  QTrPlaneStress2dClass,
  Q4AxisymmClass,
  L4AxisymmClass,
  RerShellClass,
  PPdeElementClass,
  LTrElementPPDEClass,
  Beam2dClass,   
  Beam3dClass,
  LIBeam3dClass,
  Truss3dClass,
  Truss1dClass,
  TrPlaneStrainClass,
  LIBeam3dNLClass,
  Quad1_htClass,
  Tr1_htClass,
  QuadAxisym1_htClass,
  TrAxisym1_htClass,
  Brick1_htClass,
  Tetrah1_htClass,
  
  InterfaceElem1dClass,

  DomainClass,

  EngngModelClass,
  EigenValueDynamicClass,
  LinearStaticClass,
  NonLinearStaticClass,
  DEIDynamicClass,
  NlDEIDynamicClass,
  DIIDynamicClass,
  IncrementalLinearStaticClass,
  StationaryFlowClass,
  LinearStabilityClass,
  StructuralEngngModelClass,
  PNlDEIDynamicClass,
  AdaptiveLinearStaticClass,
  AdaptiveNonLinearStaticClass,
  StationaryTransportProblemClass,
  NonStationaryTransportProblemClass,
  NLTransientTransportProblemClass,
  StaggeredProblemClass,
  CBSClass,
  SUPGClass,

  NumericalMethodClass,
  SparseLinearSystemNMClass,
  SparseNonLinearSystemNMClass,
  SparseGeneralEigenValueSystemNMClass,
  LDLTFactorizationClass,
  GeneralizedJacobiSolverClass,
  SubspaceIterationSolverClass,
  InverseIterationSolverClass,
  CylindricalALMSolverClass,
  SolutionOfLinearEquationsClass,
  NRSolverClass,
  IMLSolverClass,
  SpoolesSolverClass,
  PetscSolverClass,
  
  MaterialClass,
  StructuralMaterialClass,
  LinearElasticMaterialClass,
  IsotropicLinearElasticMaterialClass,
  OrthotropicLinearElasticMaterialClass,
  PerfectlyPlasticMaterialClass,
  Steel1MaterialClass,
  SmearedCrackingMaterialClass,
  PlasticSmearedCrackingMaterialClass,
  RCMMaterialClass,
  Concrete1Class,
  MaxwelllChainMaterialClass,
  CebFip78MaterialClass,
  DoublePowerLawMaterialClass,
  B3MaterialClass,
  Concrete2Class,
  Concrete3Class,
  DeformationTheoryMaterialClass,
  HeatMaterialClass,
  IsotropicLinearHeatMaterialClass,
  RCSDMaterialStatusClass,
  RCSDEMaterialStatusClass,
  //NonlocalMaterialExtensionClass,
  //StructuralNonlocalMaterialExtensionClass,
  RCSDNLMaterialStatusClass,
  MicroplaneMaterialClass,
  MicroplaneMaterial_BazantClass,
  M4MaterialClass,
  IsotropicDamageMaterialClass,
  IsotropicDamageMaterial1Class,
  MazarsMaterialClass,
  MicroplaneDamageMaterialClass,
  DruckerPragerPlasticitySMClass,
  MPlasticMaterialClass,
  J2MatClass,

  TransportMaterialClass,
  IsotropicHeatTransferMaterialClass,
  HeMoTKMaterialClass,
  
  HydratingIsoHeatMaterialClass,
  HydratingTransportMaterialStatusClass,
  HydratingHeMoMaterialClass,
  
  HellmichMaterialClass,
  HellmichMaterialStatusClass,
  HydrationModelClass,
  HydrationModelStatusClass,


  MaterialStatusClass,
  StructuralMaterialStatusClass,
  TransportMaterialStatusClass,
  PerfectlyPlasticMaterialStatusClass,
  SmearedCrackingMaterialStatusClass,
  PlasticSmearedCrackingMaterialStatusClass,
  Concrete2MaterialStatusClass,
  MaxwellChainMaterialStatusClass,
  RCMMaterialStatusClass,
  RCSDMaterialClass,
  RCSDEMaterialClass,
  //NonlocalMaterialStatusExtensionClass,
  //StructuralNonlocalMaterialStatusExtensionClass, 
  RCSDNLMaterialClass,
  M4MaterialStatusClass,
  IsotropicDamageMaterialStatusClass,
  MicroplaneDamageMaterialStatusClass,
  DruckerPragerPlasticitySMStatusClass,
  MPlasticMaterialStatusClass,

  FluidDynamicMaterialClass,
  FluidDynamicMaterialStatusClass,
  NewtonianFluidMaterialClass,
  BinghamFluidMaterialClass,
  BinghamFluidMaterialStatusClass,
  TwoFluidMaterialClass,

  CrossSectionClass,
  StructuralCrossSectionClass,
  SimpleCrossSectionClass,
  LayeredCrossSectionClass,
  FiberedCrossSectionClass,
  HeatCrossSectionClass,
  EmptyCrossSectionClass,

  GaussPointClass,
  IntegrationRuleClass,
  GaussIntegrationRuleClass,
  LobattoIntegrationRuleClass,
  MicroplaneClass,

  FloatArrayClass,
  FloatMatrixClass,
  IntArrayClass,
  DictionaryClass,
  TDictionaryClass,
  TimeStepClass,
  MetaStepClass,
  
  ErrorEstimatorClass,
  ScalarErrorIndicatorClass,
  ZZErrorEstimatorClass,
  CombinedZZSIErrorEstimatorClass,
  HuertaErrorEstimatorClass,
  
  RemeshingCriteriaClass,
  ZZRemeshingCriteriaClass,
  CombinedZZSIRemeshingCriteriaClass,
  HuertaRemeshingCriteriaClass,
  
  SpatialLocalizerClass,
  DummySpatialLocalizerClass,
  OctreeSpatialLocalizerClass,
  
  NonlocalBarrierClass,
  PolylineNonlocalBarrierClass,
  
  PrimaryFieldClass,
  
  LEPlicClass,
  LevelSetPCSClass,
  FastMarchingMethodClass,


  WallClockLoadBallancerMonitorClass,
  ParmetisLoadBallancerClass
    
};

/**
 Type representing kind of characteristic value (of scalar, verot or tensorial character) or 
 unknown, which is required, requested, returned, or passed to/from various general services. 
 It typically describes the physical meaning of corresponding component.
 Typically, many top base classes declare general services for requesting or computing 
 some "characteristic" values of given type. Then only one service for all values of sane type 
 (like vector, scalar) is declared, passing the type of required value (of CharType type) as parameter.
 Particular implementation based on passed CharType value, usually invokes corresponding specialized servises
 and returns result. If passed CharType value is of unsupported value, error is generated.
 @see see also ValueModeType type.
 */
enum CharType {
 // if modified , modify also isCharMtrxIncrementalValue function in cltypes.C

 // characteristic matrices 

  UnknownCharType,
  StiffnessMatrix,
  TangentStiffnessMatrix,
  SecantStiffnessMatrix,
  ElasticStiffnessMatrix,
  MassMatrix,
  LumpedMassMatrix,
  DIIModifiedStiffnessMatrix,
  ConductivityMatrix,
  CapacityMatrix,
  InitialStressMatrix,
  HeatAndMoistureCharMatrix,

  // CBS
  IntermediateConvectionTerm,
  IntermediateDiffusionTerm,
  DensityRhsVelocityTerms,
  DensityRhsPressureTerms,
  DensityPrescribedTractionPressure,
  NumberOfNodalPrescribedTractionPressureContributions,
  PressureLhs,
  CorrectionRhs,
  CriticalTimeStep,
  PrescribedVelocityRhsVector,
  PrescribedDensityRhsVector,
  // SUPG/PSPG
  AccelerationTerm_MB,
  AdvectionDerivativeTerm_MB,
  DiffusionDerivativeTerm_MB,
  SecantDiffusionDerivativeTerm_MB,
  TangentDiffusionDerivativeTerm_MB,
  InitialDiffusionDerivativeTerm_MB,
  PressureTerm_MB,
  LSICStabilizationTerm_MB,
  LinearAdvectionTerm_MC,
  AdvectionTerm_MC,
  AdvectionDerivativeTerm_MC,
  AccelerationTerm_MC,
  DiffusionDerivativeTerm_MC,
  DiffusionTerm_MC,
  PressureTerm_MC,
  BCRhsTerm_MB,
  BCRhsTerm_MC,
  AlgorithmicRhsTerm_MB,
  AlgorithmicRhsTerm_MC,
  AdvectionTerm_MB,
  DiffusionTerm_MB,


  // characteristic vectors
  
  LoadVector,
  NodalInternalForcesVector,
  LastEquilibratedNodalInternalForcesVector,
  
  ElementPPDELoadVector,
  ElementForceLoadVector,
  ElementNonForceLoadVector,
  NodalLoadVector,
  BcLhsDueToConvection,
  ElementHEMOLoadVector,
  ElementBCTransportVector,
  ElementInternalSourceVector,
  LHSBCMatrix,  // LHS due to Boundary Conditions (Transport problems)
  
  NSTP_MidpointLhs, // NonStationaryTransportProblem - LHS for midpoint disretization alg.
  NSTP_MidpointRhs, // NonStationaryTransportProblem - RHS for midpoint disretization alg.
  IntSourceLHSMatrix,  // LHS due to material internal source (Transport problems)
  ElementForceLoadVectorOfPrescribed,   // Prescribed here means corresponding to prescribed dofs
  NodalLoadVectorOfPrescribed,
  PrescribedRhsVector
  
/*
  LoadVector,
  NodalInternalForcesVector_Total,
  NodalInternalForcesVector_Incremental,
 ElementPPDELoadVector_Total,
 ElementPPDELoadVector_Incremental,
  ElementForceLoadVector_Total,
  ElementForceLoadVector_Incremental,
 ElementNonForceLoadVector_Total,
 ElementNonForceLoadVector_Incremental,
  NodalLoadVector_Total,
  NodalLoadVector_Incremental,
  BcLhsDueToConvection,
 ElementHEMOLoadVector_Total,
 ElementBCTransportVector_Total,
 ElementBCTransportVector_Incremental,
 ElementInternalSourceVector_Total,
 ElementInternalSourceVector_Incremental,
 LHSBCMatrix_Total,  // LHS due to Boundary Conditions (Transport problems)
 LHSBCMatrix_Incremental,

 ElementForceLoadVectorOfPrescribed_Total,   // Prescribed here means corresponding to prescribed dofs
 ElementForceLoadVectorOfPrescribed_Incremental,
 NodalLoadVectorOfPrescribed_Total,
 NodalLoadVectorOfPrescribed_Incremental
*/
};

//
// following type determine the mode of some value.
// which can be requested from various specialized methods.
// particular specialized methods (for example method for computing the load vector)
// are general, i.e., they are able to compute response for
// both totalLoadVector and incrementalLoadVector charTypes. 
// The particular type of response is then requested using parameter of ValueModeType type.
//

/**
 Type  representing the physical meaning of element or constitutive model internal variable.
 Values of this type are used, when these internal variables are requested.
*/
enum InternalStateType {
 IST_Undedfined                    = 0,                   
 IST_StressTensor                  = 1,
 IST_PrincipalStressTensor         = 2,
 IST_PrincipalStressTempTensor     = 3,
 IST_StrainTensor                  = 4,
 IST_PrincipalStrainTensor         = 5,
 IST_PrincipalStrainTempTensor     = 6,
 IST_BeamForceMomentumTensor       = 7,
 IST_BeamStrainCurvatureTensor     = 8,
 
 IST_ShellForceMomentumTensor      = 9,
 IST_ShellStrainCurvatureTensor    = 10,
 //IST_ForceTensor,
 //IST_MomentumTensor,
 IST_CurvatureTensor               = 11,
 IST_DisplacementVector            = 12,
 IST_DamageTensor                  = 13,
 IST_DamageInvTensor               = 14,
 IST_PrincipalDamageTensor         = 15,
 IST_PrincipalDamageTempTensor     = 16,
 IST_CrackState                    = 17,

 IST_StressTensorTemp              = 18,
 IST_StrainTensorTemp              = 19,
 IST_ForceTensorTemp               = 20,
 IST_MomentumTensorTemp            = 21,
 IST_CurvatureTensorTemp           = 22,
 IST_DisplacementVectorTemp        = 23,
 IST_DamageTensorTemp              = 24,
 IST_DamageInvTensorTemp           = 25,
 IST_CrackStateTemp                = 26,
 IST_PlasticStrainTensor           = 27,
 IST_PrincipalPlasticStrainTensor  = 28,

 IST_CylindricalStressTensor       = 29,
 IST_CylindricalStrainTensor       = 30,

 IST_MaxEquivalentStrainLevel      = 31,
 IST_ErrorIndicatorLevel           = 32,
 IST_InternalStressError           = 33,
 IST_PrimaryUnknownError           = 34,
 IST_RelMeshDensity                = 35,

 IST_MicroplaneDamageValues        = 36,

 IST_Temperature                   = 37,
 IST_MassConcentration_1           = 38,

 IST_HydrationDegree               = 39,
 IST_Humidity                      = 40,

 IST_Velocity                      = 41,
 IST_Pressure                      = 42,

 IST_VOFFraction                   = 43,
 IST_Density                       = 44,
 
 IST_MaterialInterfaceVal          = 45,


  //IST_StressVector,
  //IST_StrainVector,
  PlasticStrainVector,
  CrackStatuses,
 CrackedFlag,
  PositiveEffStrains,
  CrackStrains,
  MaxCrackStrains,
  MaxEffTotalStrains,
  ReachedSofteningStress,
  CrackDirs,
  MinEffStrainsForFullyOpenCrack,
  CharLengths,
  OldPrincipalStrain,
  PrincipalStrain,
  OldPrincipalStress,
  PrincipalStress,
  CrackMap,
  HiddenStress,
 NonlocalStrainVector

};


/// type determining the scale corresponding to particular variable
enum VarScaleType {
  VST_Length,
  VST_Velocity,
  VST_Time,
  VST_Density,
  VST_Pressure,
  VST_Force,
  VST_Viscosity,
  VST_ReynoldsNumber
};

/// enum determining the type of internal variable
enum InternalStateValueType {
 ISVT_UNDEFINED,
 ISVT_SCALAR,
 ISVT_VECTOR,
 ISVT_TENSOR_S3, // symmetric 3x3 tensor
 ISVT_TENSOR_S3E,// symmetric 3x3 tensor, packed with off diagonal components multiplied by 2 
                  // (engineering strain vector, for example)
 ISVT_TENSOR_G   // general tensor
};

/// enum determining the mode of internal variable
enum InternalStateMode {
 ISM_local,
 ISM_recovered
};


/**
 Type representing particular unknown (its physical meaning). 
*/
enum UnknownType {
  UnknownType_Unknown,
  DisplacementVector,
  GeneralizedDisplacementVector,
  //HeMaCVector, // heat and mass concetration vector
  FluxVector,
  VelocityVector,
  PressureVector,

  EigenValue,
  EigenVector,
  TotalLoadLevel,
  
  ReynoldsNumber,
  Theta_1, // CBS integration constant
  Theta_2, // CBS integration constant
  PrescribedTractionPressure // CBS prescribed pressure due to applied traction

};


/**
   This type identifies the governing equation
*/
enum EquationID {
  EID_MomentumBalance,
  EID_AuxMomentumBalance, // Auxiliary equation for characteristic-split based methods
  EID_ConservationEquation,
  EID_MomentumBalance_ConservationEquation, // Coupled system
};

//
// following mode determines the mode of particular unknown
// which can be requested on DOF.
// particular DOF contain for example displacement type unknown,
// but we can request total value, increment of value or velocity or acceleration of
// this unknown. This has been done mainly in  order to  improve runtime checking 
// ad Dof level.
//
// see also isUnknownModeIncrementalMode() function (cltypes.C)
// when adding new ValueModeType mode.
/**
 Type representing the mode of UnknownType or CharType, or similar types. 
 Afore mentioned types usually describes the physical meaning of 
 value and ValueModeType provides the further necessary classification. For example "DisplacementVector" 
  value can be futher classified to be total displacement (TotalMode) or  velocity of 
 displacement (VelocityMode) an so on.
*/
enum ValueModeType {
 VM_Unknown      = 0, 
 VM_Total        = 1,
 VM_Velocity     = 2,
 VM_Acceleration = 3,
 VM_Incremental

};

/**
 Type representing numerical component. The components of characteristic equations are mapped 
 to their corresponding numrical counterparts using these common component types.
 All numerical methods solving the same problem have to use the same and compulsory
 NumericalCmpn values. This allows to use generally any numerical method instance (even added in future) 
 without changing any code.
 */ 
enum NumericalCmpn {
  
  LinearEquationLhs,
  LinearEquationRhs,     /* for linear equation problem in the form Ax=b */
  LinearEquationSolution,
  EigenValues,
  EigenVectors,
  AEigvMtrx,
  BEigvMtrx,
  NumberOfEigenValues,
  PrescribedTolerancy,
  ReachedLevel,
  IncrementOfSolution,
  RequiredIterations,
  NonLinearLhs,
  NonLinearRhs_Total,
  NonLinearRhs_Incremental,
  InitialNonLinearRhs,
  TotalNonLinearSolution,
  StepLength,
  CurrentLevel,
  InternalRhs,
  IncrementOfNonlinearSolution,
  InternalState

};

/**
 Type representing material mode of integration point.
 */
enum MaterialMode {  // characteristic (material) mode for gaussian
                   // point (material point)
  _Unknown,
  _3dMat,
  _PlaneStress,
  _PlaneStrain,
  _2dPlate,
  _1dMat,
  _2dBeam,
  _3dBeam,
  _3dShell,
  _3dRotContinuum, // axisymmetry
  
  _2dPlateLayer,
  _2dBeamLayer,
  _3dShellLayer,
  _PlaneStressRot,

  _1dFiber,
  _3dMicroplane,
  _2dInterface,
  _1dInterface,

  _2dHeat,// 2d heat
  _2dHeMo,// 2d heat and mass (one component) transfer
  _3dHeat,
  _3dHeMo,

  _2dFlow,
  _2dAxiFlow

} ;


/**
 Describes the character of characteristic material matrix.
*/
enum MatResponseMode {
 TangentStiffness,
 SecantStiffness,
 ElasticStiffness,
 Conductivity,    // element level conductivity matrix
 Conductivity_ww, // material level conductivity submatrix
 Conductivity_hh, // material level conductivity submatrix
 Conductivity_hw, // material level conductivity submatrix
 Conductivity_wh, // material level conductivity submatrix
 Capacity,
 Capacity_ww, // material level capacity submatrix
 Capacity_hh, // material level capacity submatrix
 Capacity_hw, // material level capacity submatrix
 Capacity_wh, // material level capacity submatrix

 IntSource,
 IntSource_ww, // material level internal source submatrix - water source
 IntSource_hh, //  - heat source
 IntSource_hw, //  - heat source dependency on water content change
 IntSource_wh, //  - water source dependency on temperature change

 MRM_Density,      // material density
 MRM_Viscosity
};

/**
 Type representing the form of returned characteristic value (for cross section and material models).
 The response can be returned in so called full or reduced form.
 Generally, the full form contain all components, even if they are generally always zero (based on MaterialMode of
 given integration point). On the other hand, the reduced form contain only generally nonzero components.
 For example the "full-like" strain vector contains six components. For integration point in plane stress mode 
 the "reduced-like" strain vector contains only 3 generally nonzero components. The contens of full and reduced forms
 is defined by corresponding (material or cross section level) base classes.
 */
enum MatResponseForm {
 // identifies return form for material stiffness matrix
 ReducedForm, // only stiffness for necesarry stresses are given
 FullForm     // all component of 3d stresses are available, even if they eq 0.
 };

/**
 Type determining the type of general boundary condition.
 */
enum bcValType {
  // to identify load type
  UnknownBVT,
  TemperatureBVT,
  ForceLoadBVT,
  PressureBVT,
  HumidityBVT,
  VelocityBVT,
  DisplacementBVT
  
};
  
/**
 Type representing the geometric character of loading.
 */
enum bcGeomType {
  // to indentify load gemetry type
  UnknownBGT,
  NodalLoadBGT,
  BodyLoadBGT,
  EdgeLoadBGT,
  SurfaceLoadBGT,
  PointLoadBGT
};


/**Type representing the type of bc.*/
enum bcType {
  UnknownBT,
  DirichletBT,
  TransmissionBC, // Neumann type (prescribed flux)
  ConvectionBC    // Newton type
};

/**
 Type representing the required character of load vector. 
 */
enum LoadResponseMode {
 TotalLoad,
 IncrementOfLoad
};

/**
 Type representing material model extension.
 */
enum MaterialExtension {
  // id indicating material class extensions for run time testing
  Material_StructuralCapability,
  Material_HeatCapability,
  Material_TransportCapability,
  Material_FluidDynamicsCapability
};

/**
 Type representing cross section extension.
 */
enum CrossSectExtension {
  // id indicating material class extensions for run time testing
  CS_StructuralCapability,
  CS_HeatCapability
};

/**
Type representing element extension.
*/
enum ElementExtension {
  // id indicating element class extensions for run time testing
  Element_SurfaceLoadSupport,
  Element_EdgeLoadSupport

};

#ifdef __PARALLEL_MODE
/**
   The communicator mode determines the communication:
   
   (Static) The mode can be static, meaning that each node can assemble its communication maps
   independently (or by independent communication). This implies that the size of
   communication buffers is known in advance. Also if no data are planned to send to remote node, there
   is no communication with this node (both sender and receiver know that there will be no data to send).
   
   (Dynamic) In this case the communication pattern and the ammount of data sent between nodes is 
   not known in advance. This requires to use dynamic (packeted) buffering.
*/
enum CommunicatorMode {
  CommMode_Static,
  CommMode_Dynamic
};

enum CommBuffType {CBT_static, CBT_dynamic};

#endif

#ifdef __OOFEG

enum DrawMode {
  unknown,
  
  rawGeometry,
  deformedGeometry,
  eigenVectorGeometry,
  nodeAnnotation,
  appliedPrimaryBc,
  
  internalStateBegin,
  mxForce,
  myForce,
  mzForce,
  myzForce,
  mzxForce,
  mxyForce,
  //qxzForce,
  //qyzForce,
  sxForce,
  syForce,
  szForce,
  syzForce,
  szxForce,
  sxyForce,
  yieldState,
  crackedState,
  stressErrorState,
  requiredAdaptiveMeshSizeState,
  damageLevel,
  errorIndicatorLevel,
  relativeMeshSizeDensity,
  temperatureField,
  massConcentration1Field,
  velocityField,
  pressureField,
  vofField,
  densityField,

  hydrationDegreeState,
  humidityState,

  internalStateEnd
  
} ;

enum MatStatusVar {
  ms_unknown,
  ms_yield_flag,
  ms_isCracked_flag,
  ms_crackDirection_flag,
  ms_crackStatus_flag,
  ms_damage_flag

} ;

#endif

enum ContextOutputMode {
  NOCONTEXT,
  ALWAYS,         // enable for post-processing
  REQUIRED,       // if required (for backtracing computation)
  USERDEFINED     // input attribute of domain (each n-th step)
};

/**
 Determines the input/output mode of context file.
 contextMode_read - context file is openned for reading.
 contextMode_write- context mode is openned for writting,
                    if not exist then is created, otherwise 
           it will be truncated.
*/
enum ContextFileMode {
 contextMode_read,
 contextMode_write
};
  
/**
 Type representing the type of formulation (total or updated) of non-linear computation.
 */
enum fMode {    
  // type of formulation (non-linear computation only)
  UNKNOWN,
  TL,        // Total Lagrange
  AL         // Updated Lagrange
};
    
/** 
 Type representing type of domain.
 Domain type (the member value of Domain class) is used to determine the default 
 number of DOFs per node and side and to determine their corresponding physical meaning.
 */
enum domainType {
  _unknownMode,
  _2dPlaneStressMode,
  _PlaneStrainMode,
  _2dPlaneStressRotMode,
  _3dMode,
  _3dAxisymmMode,
  _2dMindlinPlateMode,
  _3dShellMode,
  _2dTrussMode,
  _1dTrussMode,
  _2dBeamMode,
  _HeatTransferMode,
  _HeatMass1Mode,   // Coupled heat and mass (1 matter) transfer
  _2dIncompressibleFlow // 2d Incompressible flow, no energy eq
};

enum problemMode {
  _processor,
  _postProcessor
};


enum stressStrainPrincMode {
/*
   we have only one algorithm for computing eigen values and vectors
   in order to be able to distinguish between some different modes we define
   this new type
*/
  principal_strain, // for computing principal_strains - from engineering strains
  principal_stress, // for computing principal_stress
  principal_deviatoricstress
};

enum integrationDomain {
/*
Used by integrator class to suply 
integration points for proper domain to be 
integrated (Area,Volume and its shape)

*/
  _Unknown_integrationDomain,
  _Line,
  _Triangle,
  _Square,
  _Cube,
  _Tetrahedra,

  _Embedded2dLine
};

/* mask defing NumMetod Status; which can be asked after
 finishing computation by Numerical Method.
 this mask should report some sitiuation, for
 exaple:
 None,
 Success, NoSuccess  -> at now inform about succes,
 if no succes, currently NoSuccess is not reported 
 but NoSuccess causes exit();
 
 some more detailed messages can be incorporated,
 like:
 
 KeepTangent -> used by some non-linear solvers telling
                don't assemble new tangent, but use previous.

 note: every mask flag begins with NM_ to avoid possible multiple
 definition
*/
 
typedef unsigned long  NM_Status;
/* Mask selecting status */

#define NM_None         0
#define NM_Success      (1L << 1) 
#define NM_NoSuccess    (1L << 2)
#define NM_KeepTangent  (1L << 3)
#define NM_ForceRestart (1L << 4)

      
/* Dof Type, determines the type of DOF created 
 */
enum dofType {
  DT_master = 0,
  DT_simpleSlave = 1,
  DT_slave = 2
};


/* mask definning the physical meaning of particular DOF in node.
  mask array are also used in elements, where these arrays
  are determining required DOFs needed by element and which are then 
  requsted on particular nodes. They are used in function 
 Node::giveLOcationArray which returns equations numbers
 corresponding to selected dofs. 
*/

typedef char DofID;
/** 
 Type representing particular dof type. Values of this type describe the physical meaning of
 available DOFs. 
 Note: implementation of Node::computeGNTransformation rely on D_u, D_v and D_w (R_u, R_v, R_w) order.
 Do not change their order and do not insert any values between these values. 
 */
enum DofIDItem {
 Undef, // Erorr value
 D_u = 1,   // u-displacement (in direction of x-axis)
 D_v = 2,   // v-displacement (in direction of y-axis)
 D_w = 3,   // w-displacement (in direction of z-axis)
 R_u = 4,   // rotation around x-axis (right hand rule assumed)
 R_v = 5,   // rotation around y-axis
 R_w = 6,   // rotation around z-axis 

 V_u = 7,   // u-velocity (in direction of x-axis)
 V_v = 8,   // v-velocity (in direction of y-axis)
 V_w = 9,   // w-velocity (in direction of z-axis)

 T_f =10,   // temperature field
 P_f =11,   // pressure field
 G_0 =12,   // DOF for gradient formulation no. 0
 G_1 =13,   // DOF for gradient formulation no. 1
 C_1 =14    // mass concentration of the first constituent
  
};
// max length of text string with DofIdName + 1
// see Dof::giveDofIDName function
#define DofIdNameMaxLength 5

/// StateCounterType type used to indicate solution state.
typedef unsigned long StateCounterType;


/**
 Enumerative type, used to specify type of trasformation required from dofManager (node).
 The _toGlobalCS value requires transformation from node-depenedent coordinate system
 to gbal coordinte system in node to be assebled. Then global vector fg can be obtained by
 followwing operation fg = T fn, where T is transformation matrix and fn is vector expressed in 
 nodal coordinate system).
 The _toNodalCS value represent transformation from global c.s in node to node-dependent
 coordinate system. 
*/
enum DofManTrasfType {
 _toGlobalCS,
 _toNodalCS
};

/**
 Enumerative type, used to identify interface type.
 @see Interface class for details
*/
enum InterfaceType {
 UnknownInterfaceType,

 LayeredCrossSectionInterfaceType,
 FiberedCrossSectionInterfaceType,

 ZZNodalRecoveryModelInterfaceType,
 NodalAveragingRecoveryModelInterfaceType,
 SPRNodalRecoveryModelInterfaceType,

 ZZErrorEstimatorInterfaceType,
 HuertaErrorEstimatorInterfaceType,
 Huerta1dErrorEstimatorInterfaceType, // experimental

 DirectErrorIndicatorRCInterfaceType,
 ZZRemeshingCriteriaInterfaceType,
 HuertaRemeshingCriteriaInterfaceType,

 SpatialLocalizerInterfaceType,

 EIPrimaryUnknownMapperInterfaceType,
 MMAShapeFunctProjectionInterfaceType,
 EIPrimaryFieldInterfaceType,

 NonlocalMaterialStatusExtensionInterfaceType,

 NonlocalMaterialExtensionInterfaceType,
 NonlocalMaterialStiffnessInterfaceType,
 MaterialModelMapperInterfaceType,

 HydrationModelInterfaceType,
 HydrationModelStatusInterfaceType,

 LEPlicElementInterfaceType,
 LevelSetPCSElementInterfaceType,
};


/**
 Enumerative type used to identify the sparse matrix type
*/
enum SparseMtrxType {
 SMT_Skyline,          // symmetric skyline
 SMT_SkylineU,         // unsymmetric skyline
 SMT_CompCol,          // compressed column
 SMT_DynCompCol,       // dynamically growing compressed column
 SMT_SymCompCol,       // symmetric compressed column
 SMT_DynCompRow,       // dynamically growing compressed row
 SMT_SpoolesMtrx,      // spooles sparse mtrx representation
 SMT_PetscMtrx,        // PETSc library mtrx representation
 SMT_DSS_sym_LDL,  // Richard Vondracek's sparse direct solver 
 SMT_DSS_sym_LL,   // Richard Vondracek's sparse direct solver 
 SMT_DSS_unsym_LU // Richard Vondracek's sparse direct solver 
 
};


enum LinSystSolverType {
 ST_Direct, 
 ST_CG, 
 ST_GMRES,
 ST_Spooles,
 ST_Petsc,
 ST_DSS,
 ST_Feti
};

enum GenEigvalSolverType {
  GES_SubspaceIt,
  GES_InverseIt
};


enum ErrorEstimatorType {
 EET_SEI,      // Scalar Error Indicator 
 EET_ZZEE,     // Zienkiewicz-Zhu EE
 EET_CZZSI,    // Combined ZZ and ScalarIndicator EE
 EET_HEE       // Huerta EE
};

/**
 Enumerative type used to classify element geometry
 Poosible values are:
 EGT_line_1 - line elements with two nodes  1-------2
 EGT_line_2 - line element with three nodes 1---2---3
 EGT_triangle_1 - triangle element with three nodes
 EGT_triangle_2 - triangle element with 6 nodes
                       3
                    6     5
                 1     4     2
 
 EGT_quad_1 - quadrialateral with 4 nodes
  EGT_tetra_1 - tetrahedron with 4 nodes
 EGT_hexa_1  - hexahedron with 8 nodes
*/
enum Element_Geometry_Type {
 EGT_line_1,   // line elements with two nodes  1-------2
 EGT_line_2,   // line element with three nodes 1---2---3
 EGT_triangle_1, // triangle element with three nodes
 EGT_triangle_2, // triangle element with 6 nodes
 EGT_quad_1,     // quadrialateral with 4 nodes
 EGT_tetra_1,    // tetrahedron with 4 nodes
 EGT_hexa_1,     // hexahedron with 8 nodes
 
 EGT_unknown     // unknown element geometry type
};

/**
   Type allowing to specify the required renumbering scheme;
   One can have a renumbering scheme for dof managers
   and another one for elements;
 */
enum EntityRenumberingScheme {
  ERS_DofManager,
  ERS_Element
};

enum contextIOResultType {
 CIO_OK = 0,        // ok
 CIO_BADVERSION,    // incompatible context file
 CIO_BADOBJ,        // bad object passed 
 CIO_IOERR          // general io error
};

/// context io exception class
class ContextIOERR 
{
 contextIOResultType error;
 char *msg, *file;
 int line;

public:

 ContextIOERR (contextIOResultType e, char* file, int line);
 ContextIOERR (contextIOResultType e, char* msg , char* file, int line);
 ~ContextIOERR ();

 void print ();

};

#define THROW_CIOERR(e) throw ContextIOERR (e,__FILE__, __LINE__);
#define THROW_CIOERRM(e,m) throw ContextIOERR (e,m,__FILE__, __LINE__);

/**
   Context mode (mask), defining the type of information written/read to/from context
*/
typedef unsigned long  ContextMode;
/* Mask selecting status */
#define CM_None         0
#define CM_State        (1L << 1) 
#define CM_Definition   (1L << 2)


/// oofem terminate exception class


class OOFEM_Terminate
{
public:
  enum OOFEM_exit_status {
    ES_OK,
    ES_UnknownErr
  };

  OOFEM_exit_status status;
  OOFEM_Terminate (OOFEM_exit_status s = ES_OK) {status = s;}
};


//char   cltypesGiveUnknownTypeKey (UnknownType type);
char   cltypesGiveUnknownTypeModeKey (ValueModeType mode);
int    isUnknownTypeModeIncremental (ValueModeType) ;
/// Returns the value type of corresponding InternalStateType
InternalStateValueType giveInternalStateValueType (InternalStateType type);

#define cltypes_h
#endif












