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


/*
  The enum declaration and to_string conversion
  inspired by X-Macros technique, as described in Wikipedia entry on C preprocessor 
  (http://en.wikipedia.org/wiki/C_preprocessor)
*/


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
  //RemoteMasterDofClass,
  //SharedMasterDofClass,
  //NullDofClass,
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
 // if modified , modify also isCharMtrxIncrementalValue function in cltypes.C

#define CharType_DEF				\
  ENUM_ITEM(UnknownCharType)			\
  ENUM_ITEM(StiffnessMatrix)			\
  ENUM_ITEM(TangentStiffnessMatrix)		\
  ENUM_ITEM(SecantStiffnessMatrix)		\
  ENUM_ITEM(ElasticStiffnessMatrix)		\
  ENUM_ITEM(MassMatrix)				\
  ENUM_ITEM(LumpedMassMatrix)			\
  ENUM_ITEM(DIIModifiedStiffnessMatrix)		\
  ENUM_ITEM(ConductivityMatrix)			\
  ENUM_ITEM(CapacityMatrix)			\
  ENUM_ITEM(InitialStressMatrix)		\
  ENUM_ITEM(HeatAndMoistureCharMatrix)		\
  /* CBS */					\
  ENUM_ITEM(IntermediateConvectionTerm)		\
  ENUM_ITEM(IntermediateDiffusionTerm)		\
  ENUM_ITEM(DensityRhsVelocityTerms)		\
  ENUM_ITEM(DensityRhsPressureTerms)		\
  ENUM_ITEM(DensityPrescribedTractionPressure)		\
  ENUM_ITEM(NumberOfNodalPrescribedTractionPressureContributions)	\
  ENUM_ITEM(PressureLhs)						\
  ENUM_ITEM(CorrectionRhs)						\
  ENUM_ITEM(CriticalTimeStep)						\
  ENUM_ITEM(PrescribedVelocityRhsVector)				\
  ENUM_ITEM(PrescribedDensityRhsVector)					\
  /* SUPG/PSPG */							\
  ENUM_ITEM(AccelerationTerm_MB)					\
  ENUM_ITEM(AdvectionDerivativeTerm_MB)					\
  ENUM_ITEM(DiffusionDerivativeTerm_MB)					\
  ENUM_ITEM(SecantDiffusionDerivativeTerm_MB)				\
  ENUM_ITEM(TangentDiffusionDerivativeTerm_MB)				\
  ENUM_ITEM(InitialDiffusionDerivativeTerm_MB)				\
  ENUM_ITEM(PressureTerm_MB)						\
  ENUM_ITEM(LSICStabilizationTerm_MB)					\
  ENUM_ITEM(LinearAdvectionTerm_MC)					\
  ENUM_ITEM(AdvectionTerm_MC)						\
  ENUM_ITEM(AdvectionDerivativeTerm_MC)					\
  ENUM_ITEM(AccelerationTerm_MC)					\
  ENUM_ITEM(DiffusionDerivativeTerm_MC)					\
  ENUM_ITEM(DiffusionTerm_MC)						\
  ENUM_ITEM(PressureTerm_MC)						\
  ENUM_ITEM(BCRhsTerm_MB)						\
  ENUM_ITEM(BCRhsTerm_MC)						\
  ENUM_ITEM(AlgorithmicRhsTerm_MB)					\
  ENUM_ITEM(AlgorithmicRhsTerm_MC)					\
  ENUM_ITEM(AdvectionTerm_MB)						\
  ENUM_ITEM(DiffusionTerm_MB)						\
  /* characteristic vectors */						\
  ENUM_ITEM(LoadVector)							\
  ENUM_ITEM(NodalInternalForcesVector)					\
  ENUM_ITEM(LastEquilibratedNodalInternalForcesVector)			\
  ENUM_ITEM(ElementPPDELoadVector)					\
  ENUM_ITEM(ElementForceLoadVector)					\
  ENUM_ITEM(ElementNonForceLoadVector)					\
  ENUM_ITEM(NodalLoadVector)						\
  ENUM_ITEM(BcLhsDueToConvection)					\
  ENUM_ITEM(ElementHEMOLoadVector)					\
  ENUM_ITEM(ElementBCTransportVector)					\
  ENUM_ITEM(ElementInternalSourceVector)				\
  ENUM_ITEM(LHSBCMatrix)  /* LHS due to Boundary Conditions (Transport problems) */ \
  ENUM_ITEM(NSTP_MidpointLhs) /* NonStationaryTransportProblem - LHS for midpoint disretization alg. */	\
  ENUM_ITEM(NSTP_MidpointRhs) /* NonStationaryTransportProblem - RHS for midpoint disretization alg. */	\
  ENUM_ITEM(IntSourceLHSMatrix)  /* LHS due to material internal source (Transport problems) */	\
  ENUM_ITEM(ElementForceLoadVectorOfPrescribed)  /* Prescribed here means corresponding to prescribed dofs*/ \
  ENUM_ITEM(NodalLoadVectorOfPrescribed)				\
  ENUM_ITEM(PrescribedRhsVector)


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
#define InternalStateType_DEF\
  ENUM_ITEM_WITH_VALUE(IST_Undedfined,0)\
  ENUM_ITEM_WITH_VALUE(IST_StressTensor,1)\
  ENUM_ITEM_WITH_VALUE(IST_PrincipalStressTensor,2)\
  ENUM_ITEM_WITH_VALUE(IST_PrincipalStressTempTensor,3)\
  ENUM_ITEM_WITH_VALUE(IST_StrainTensor,4)\
  ENUM_ITEM_WITH_VALUE(IST_PrincipalStrainTensor,5)\
  ENUM_ITEM_WITH_VALUE(IST_PrincipalStrainTempTensor,6)\
  ENUM_ITEM_WITH_VALUE(IST_BeamForceMomentumTensor,7)\
  ENUM_ITEM_WITH_VALUE(IST_BeamStrainCurvatureTensor,8)\
  ENUM_ITEM_WITH_VALUE(IST_ShellForceMomentumTensor,9)\
  ENUM_ITEM_WITH_VALUE(IST_ShellStrainCurvatureTensor,10)\
  ENUM_ITEM_WITH_VALUE(IST_CurvatureTensor,11)\
  ENUM_ITEM_WITH_VALUE(IST_DisplacementVector,12)\
  ENUM_ITEM_WITH_VALUE(IST_DamageTensor,13)\
  ENUM_ITEM_WITH_VALUE(IST_DamageInvTensor,14)\
  ENUM_ITEM_WITH_VALUE(IST_PrincipalDamageTensor,15)\
  ENUM_ITEM_WITH_VALUE(IST_PrincipalDamageTempTensor,16)\
  ENUM_ITEM_WITH_VALUE(IST_CrackState,17)\
  ENUM_ITEM_WITH_VALUE(IST_StressTensorTemp,18)\
  ENUM_ITEM_WITH_VALUE(IST_StrainTensorTemp,19)\
  ENUM_ITEM_WITH_VALUE(IST_ForceTensorTemp,20)\
  ENUM_ITEM_WITH_VALUE(IST_MomentumTensorTemp,21)\
  ENUM_ITEM_WITH_VALUE(IST_CurvatureTensorTemp,22)\
  ENUM_ITEM_WITH_VALUE(IST_DisplacementVectorTemp,23)\
  ENUM_ITEM_WITH_VALUE(IST_DamageTensorTemp,24)\
  ENUM_ITEM_WITH_VALUE(IST_DamageInvTensorTemp,25)\
  ENUM_ITEM_WITH_VALUE(IST_CrackStateTemp,26)\
  ENUM_ITEM_WITH_VALUE(IST_PlasticStrainTensor,27)\
  ENUM_ITEM_WITH_VALUE(IST_PrincipalPlasticStrainTensor,28)\
  ENUM_ITEM_WITH_VALUE(IST_CylindricalStressTensor,29)\
  ENUM_ITEM_WITH_VALUE(IST_CylindricalStrainTensor,30)\
  ENUM_ITEM_WITH_VALUE(IST_MaxEquivalentStrainLevel,31)\
  ENUM_ITEM_WITH_VALUE(IST_ErrorIndicatorLevel,32)\
  ENUM_ITEM_WITH_VALUE(IST_InternalStressError,33)\
  ENUM_ITEM_WITH_VALUE(IST_PrimaryUnknownError,34)\
  ENUM_ITEM_WITH_VALUE(IST_RelMeshDensity,35)\
  ENUM_ITEM_WITH_VALUE(IST_MicroplaneDamageValues,36)\
  ENUM_ITEM_WITH_VALUE(IST_Temperature,37)\
  ENUM_ITEM_WITH_VALUE(IST_MassConcentration_1,38)\
  ENUM_ITEM_WITH_VALUE(IST_HydrationDegree,39)\
  ENUM_ITEM_WITH_VALUE(IST_Humidity,40)\
  ENUM_ITEM_WITH_VALUE(IST_Velocity,41)\
  ENUM_ITEM_WITH_VALUE(IST_Pressure,42)\
  ENUM_ITEM_WITH_VALUE(IST_VOFFraction,43)\
  ENUM_ITEM_WITH_VALUE(IST_Density,44)\
  ENUM_ITEM_WITH_VALUE(IST_MaterialInterfaceVal,45)\
  \
  ENUM_ITEM(CrackStatuses)\
  ENUM_ITEM(CrackedFlag)\
  ENUM_ITEM(CrackDirs)\


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

#define UnknownType_DEF\
  ENUM_ITEM(UnknownType_Unknown)\
  ENUM_ITEM(DisplacementVector)\
  ENUM_ITEM(GeneralizedDisplacementVector)\
  ENUM_ITEM(FluxVector)\
  ENUM_ITEM(VelocityVector)			\
  ENUM_ITEM(PressureVector)			\
  ENUM_ITEM(TemperatureVector)			\
  ENUM_ITEM(EigenValue)				\
  ENUM_ITEM(EigenVector)			\
  ENUM_ITEM(TotalLoadLevel)			\
  ENUM_ITEM(ReynoldsNumber)						\
  ENUM_ITEM(Theta_1) /* CBS integration constan)*/			\
  ENUM_ITEM(Theta_2) /* CBS integration constan)*/			\
  ENUM_ITEM(PrescribedTractionPressure) /* CBS prescribed pressure due to applied tractio)*/


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
#define ValueModeType_DEF\
  ENUM_ITEM_WITH_VALUE(VM_Unknown,0)\
  ENUM_ITEM_WITH_VALUE(VM_Total,1)		\
  ENUM_ITEM_WITH_VALUE(VM_Velocity,2)		\
  ENUM_ITEM_WITH_VALUE(VM_Acceleration,3)	\
  ENUM_ITEM_WITH_VALUE(VM_Incremental,4)	\
  ENUM_ITEM_WITH_VALUE(VM_RhsTotal,5)		\
  ENUM_ITEM_WITH_VALUE(VM_RhsIncremental,6)	\
  ENUM_ITEM_WITH_VALUE(VM_RhsInitial,7)


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
#define MaterialMode_DEF\
  ENUM_ITEM(_Unknown)	\
  ENUM_ITEM(_3dMat)\
  ENUM_ITEM(_PlaneStress)\
  ENUM_ITEM(_PlaneStrain)\
  ENUM_ITEM(_2dPlate)\
  ENUM_ITEM(_1dMat)\
  ENUM_ITEM(_2dBeam)\
  ENUM_ITEM(_3dBeam)\
  ENUM_ITEM(_3dShell)\
  ENUM_ITEM(_3dRotContinuum) /* axisymmetry */\
  \
  ENUM_ITEM(_2dPlateLayer)\
  ENUM_ITEM(_2dBeamLayer)\
  ENUM_ITEM(_3dShellLayer)\
  ENUM_ITEM(_PlaneStressRot)\
  \
  ENUM_ITEM(_1dFiber)\
  ENUM_ITEM(_3dMicroplane)\
  ENUM_ITEM(_2dInterface)\
  ENUM_ITEM(_1dInterface)\
  \
  ENUM_ITEM(_2dHeat) /* 2d heat */ \
  ENUM_ITEM(_2dHeMo) /* 2d heat and mass (one component) transfer */\
  ENUM_ITEM(_3dHeat)\
  ENUM_ITEM(_3dHeMo)\
  \
  ENUM_ITEM(_2dFlow)\
  ENUM_ITEM(_2dAxiFlow)\
  ENUM_ITEM(_3dFlow)\


/**
 Describes the character of characteristic material matrix.
*/
#define MatResponseMode_DEF\
  ENUM_ITEM(TangentStiffness)\
  ENUM_ITEM(SecantStiffness) \
  ENUM_ITEM(ElasticStiffness)						\
  ENUM_ITEM(Conductivity)    /* element level conductivity matrix */	\
  ENUM_ITEM(Conductivity_ww) /* material level conductivity submatrix */ \
  ENUM_ITEM(Conductivity_hh) /* material level conductivity submatrix */ \
  ENUM_ITEM(Conductivity_hw) /* material level conductivity submatrix */ \
  ENUM_ITEM(Conductivity_wh) /* material level conductivity submatrix */ \
  ENUM_ITEM(Capacity)							\
  ENUM_ITEM(Capacity_ww) /* material level capacity submatrix */	\
  ENUM_ITEM(Capacity_hh) /* material level capacity submatrix */	\
  ENUM_ITEM(Capacity_hw) /* material level capacity submatrix */	\
  ENUM_ITEM(Capacity_wh) /* material level capacity submatrix */	\
  ENUM_ITEM(IntSource)							\
  ENUM_ITEM(IntSource_ww) /* material level internal source submatrix - water source */ \
  ENUM_ITEM(IntSource_hh) /*  - heat source */				\
  ENUM_ITEM(IntSource_hw) /*  - heat source dependency on water content change */ \
  ENUM_ITEM(IntSource_wh) /*  - water source dependency on temperature change */ \
  ENUM_ITEM(MRM_Density)  /* material density */			\
  ENUM_ITEM(MRM_Viscosity)


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
#define domainType_DEF\
  ENUM_ITEM(_unknownMode)\
  ENUM_ITEM(_2dPlaneStressMode)\
  ENUM_ITEM(_PlaneStrainMode)\
  ENUM_ITEM(_2dPlaneStressRotMode)\
  ENUM_ITEM(_3dMode)\
  ENUM_ITEM(_3dAxisymmMode)\
  ENUM_ITEM(_2dMindlinPlateMode)\
  ENUM_ITEM(_3dShellMode)\
  ENUM_ITEM(_2dTrussMode)\
  ENUM_ITEM(_1dTrussMode)\
  ENUM_ITEM(_2dBeamMode)\
  ENUM_ITEM(_HeatTransferMode)\
  ENUM_ITEM(_HeatMass1Mode)   /* Coupled heat and mass (1 matter) transfer */ \
  ENUM_ITEM(_2dIncompressibleFlow) /* 2d Incompressible flow, no energy eq */ \
  ENUM_ITEM(_3dIncompressibleFlow) /* 3d Incompressible flow, no energy eq */ 


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
#define dofType_DEF\
  ENUM_ITEM_WITH_VALUE(DT_master,0)\
  ENUM_ITEM_WITH_VALUE(DT_simpleSlave,1)\
  ENUM_ITEM_WITH_VALUE(DT_slave,2)


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
#define DofIDItem_DEF\
  ENUM_ITEM_WITH_VALUE(Undef,0) /* Erorr value */	\
  ENUM_ITEM_WITH_VALUE(D_u,1)   /* u-displacement (in direction of x-axis) */ \
  ENUM_ITEM_WITH_VALUE(D_v,2)   /* v-displacement (in direction of y-axis) */ \
  ENUM_ITEM_WITH_VALUE(D_w,3)   /* w-displacement (in direction of z-axis) */ \
  ENUM_ITEM_WITH_VALUE(R_u,4)   /* rotation around x-axis (right hand rule assumed) */ \
  ENUM_ITEM_WITH_VALUE(R_v,5)   /* rotation around y-axis */		\
  ENUM_ITEM_WITH_VALUE(R_w,6)   /* rotation around z-axis */		\
  \
  ENUM_ITEM_WITH_VALUE(V_u,7)   /* u-velocity (in direction of x-axis) */ \
  ENUM_ITEM_WITH_VALUE(V_v,8)   /* v-velocity (in direction of y-axis) */ \
  ENUM_ITEM_WITH_VALUE(V_w,9)   /* w-velocity (in direction of z-axis) */ \
  \
  ENUM_ITEM_WITH_VALUE(T_f,10)  /* temperature field */	\
  ENUM_ITEM_WITH_VALUE(P_f,11)  /* pressure field */			\
  ENUM_ITEM_WITH_VALUE(G_0,12)  /* DOF for gradient formulation no. 0 */ \
  ENUM_ITEM_WITH_VALUE(G_1,13)  /* DOF for gradient formulation no. 1 */ \
  ENUM_ITEM_WITH_VALUE(C_1,14)  /* mass concentration of the first constituent */ 
  

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

/*
  The values of this type should be related not to specific solvers, 
  but more to specific packages that provide linear solver interface
  (possibly with many solver types) and are represented by a class
  derived from SparseLinearSystemNM. 
  The selection of particular solver from package should be done using keywords,
  related to particular package.
*/
enum LinSystSolverType {
 ST_Direct = 0, 
 ST_IML    = 1,
 ST_Spooles= 2,
 ST_Petsc  = 3,
 ST_DSS    = 4,
 ST_Feti   = 5
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
#define Element_Geometry_Type_DEF\
  ENUM_ITEM(EGT_line_1)   /* line elements with two nodes  1-------2 */\
  ENUM_ITEM(EGT_line_2)   /* line element with three nodes 1---2---3 */	\
  ENUM_ITEM(EGT_triangle_1) /* triangle element with three nodes */\
  ENUM_ITEM(EGT_triangle_2) /* triangle element with 6 nodes */\
  ENUM_ITEM(EGT_quad_1)     /* quadrialateral with 4 nodes */\
  ENUM_ITEM(EGT_tetra_1)    /* tetrahedron with 4 nodes */\
  ENUM_ITEM(EGT_hexa_1)     /* hexahedron with 8 nodes */ \
  ENUM_ITEM(EGT_unknown)    /* unknown element geometry type */


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
 const char *msg, *file;
 int line;

public:

 ContextIOERR (contextIOResultType e, const char* file, int line);
 ContextIOERR (contextIOResultType e, const char* msg , const char* file, int line);
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
#define CM_DefinitionGlobal (1L << 3)
#define CM_UnknownDictState (1L << 4)


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

/*
  The enum declaration and to_string conversion
  inspired by X-Macros technique, as described in Wikipedia entry on C preprocessor 
  (http://en.wikipedia.org/wiki/C_preprocessor)
*/


#define ENUM_ITEM(element) element,
#define ENUM_ITEM_WITH_VALUE(element,val) element=val,
enum InternalStateType {
  InternalStateType_DEF
};

enum UnknownType {
  UnknownType_DEF
};


enum dofType {
  dofType_DEF
};

enum domainType {
  domainType_DEF
};

enum MaterialMode {
  MaterialMode_DEF
};

enum Element_Geometry_Type {
  Element_Geometry_Type_DEF
};

enum ValueModeType {
  ValueModeType_DEF
};

enum MatResponseMode {
  MatResponseMode_DEF
};

enum DofIDItem {
  DofIDItem_DEF
};

enum CharType {
  CharType_DEF
};

char* __InternalStateTypeToString (InternalStateType _value);
char* __UnknownTypeToString (UnknownType _value);
char* __dofTypeToString (dofType _value);
char* __domainTypeToString (domainType _value);
char* __MaterialModeToString (MaterialMode _value);
char* __Element_Geometry_TypeToString(Element_Geometry_Type _value);
char* __ValueModeTypeToString(ValueModeType _value);
char* __MatResponseModeToString(MatResponseMode _value);
char* __DofIDItemToString(DofIDItem _value);
char* __CharTypeToString(CharType _value);

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE


















//char   cltypesGiveUnknownTypeKey (UnknownType type);
char   cltypesGiveUnknownTypeModeKey (ValueModeType mode);
int    isUnknownTypeModeIncremental (ValueModeType) ;
/// Returns the value type of corresponding InternalStateType
InternalStateValueType giveInternalStateValueType (InternalStateType type);

#define cltypes_h
#endif

