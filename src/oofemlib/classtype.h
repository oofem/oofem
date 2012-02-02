/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2010   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

//
// FILE: classtype.h
//

#ifndef classtype_h
#define classtype_h

namespace oofem {
/**
 * Type introduced to distinguish between classes. Intended for run-time
 * type checking. Every derived class from FemComponent base class should
 * overload the abstract giveClassID service. This service can be used
 * to distinguish between multiple derived classes.
 */
enum classType {
    FEMComponentClass,

    GeneralBoundaryConditionClass,
    ActiveBoundaryConditionClass,
    SurfaceTensionBoundaryConditionClass,
    LoadClass,
    BoundaryConditionClass,
    PrescribedGradientClass,
    MixedGradientClass,
    InitialConditionClass,
    NodalLoadClass,
    DeadWeightClass,
    BodyLoadClass,
    StructuralTemperatureLoadClass,
    StructuralEigenstrainLoadClass,
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
    ActiveDofClass,
    DofManagerClass,
    NodeClass,
    SlaveNodeClass,
    RigidArmNodeClass,
    HangingNodeClass,
    ElementSideClass,
    ElementDofManagerClass,
    RigidArmSlaveDofClass,
    HangingDofClass,
    ParticleClass,
    CohesiveSurface3dClass,
    CohesiveInterfaceMaterialClass,
    CohesiveInterfaceMaterialStatusClass,

    ElementClass,
    StructuralElementClass,
    TransportElementClass,
    FMElementClass,
    CBSElementClass,
    SUPGElementClass,
    Tr21StokesElementClass,
    Tet21StokesElementClass,
    LineSurfaceTensionElementClass,
    Line2SurfaceTensionElementClass,

    BoundaryElementClass,
    Line2BoundaryElementClass,

    GradDpElementClass,
    TrStructuralElementClass,
    NLStructuralElementClass,
    PlaneStress2dClass,
    Quad1PlaneStrainClass,
    QPlaneStrainClass,
    QPlaneStrainGradClass,
    QTrPlaneStrainClass,
    QTrPlaneStrainGradClass,
    QPlaneStress2dClass,
    QPlaneStressGradClass,
    TrPlaneStress2dClass,
    QTrPlaneStressGradClass,
    CCTPlateClass,
    CCTPlate3dClass,
    LSpaceClass,
    LSpaceBBClass,
    QSpaceClass,
    QSpaceGradClass,
    LTRSpaceClass,
    Truss2dClass,
    LIBeam2dClass,
    TrPlaneStrRotClass,
    TrPlaneStrRot3dClass,
    Axisymm3dClass,
    QTrPlaneStress2dClass,
    Q4AxisymmClass,
    L4AxisymmClass,
    RerShellClass,
    TR_SHELL01Class,
    PPdeElementClass,
    LTrElementPPDEClass,
    Beam2dClass,
    Beam3dClass,
    LIBeam3dClass,
    Truss3dClass,
    Truss1dClass,
    QTruss1dClass,
    QTruss1dGradClass,
    TrPlaneStrainClass,
    LIBeam3dNLClass,
    Quad1_htClass,
    Quad1_hmtClass,
    Tr1_htClass,
    QuadAxisym1_htClass,
    QuadAxisym1_hmtClass,
    TrAxisym1_htClass,
    Brick1_htClass,
    Brick1_hmtClass,
    Tetrah1_htClass,
    Tetrah1_hmtClass,
    PlaneStress2dXfemClass,
    LumpedMassElementClass,
    SpringElementClass,

    MacroLSpaceClass,
    MicroMaterialClass,

    InterfaceElem1dClass,
    InterfaceElem2dQuadClass,
    InterfaceElement3dTrLinClass,

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
    StokesFlowClass,
    StokesFlowStressHomogenizationClass,

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
    SlepcSolverClass,
    PETScSNESNMClass,

    MaterialClass,
    DummyMaterialClass,
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
    MaxwellChainMaterialClass,
    KelvinChainMaterialClass,
    KelvinChainSolidMaterialClass,
    CebFip78MaterialClass,
    DoublePowerLawMaterialClass,
    B3MaterialClass,
    B3SolidMaterialClass,
    MPSMaterialClass,
    Concrete2Class,
    Concrete3Class,
    DeformationTheoryMaterialClass,
    HeatMaterialClass,
    IsotropicLinearHeatMaterialClass,
    RCSDMaterialStatusClass,
    RCSDEMaterialStatusClass,
    RCSDNLMaterialStatusClass,
    MicroplaneMaterialClass,
    MicroplaneMaterial_BazantClass,
    M4MaterialClass,
    IsotropicDamageMaterialClass,
    IsotropicDamageMaterial1Class,
    MazarsMaterialClass,
    CompositeDamageMaterialClass,
    MicroplaneDamageMaterialClass,
    DruckerPragerPlasticitySMClass,
    MPlasticMaterialClass,
    J2MatClass,
    HyperElasticMaterialClass,
    MisesMatClass,
    MisesMatGradClass,
    MisesMatNlClass,
    RankineMatClass,
    RheoChainMaterialClass,
    TrabBoneMaterialClass,
    TrabBone3DClass,
    TrabBoneEmbedClass,
    ConcreteDPMClass,
    FE2SinteringMaterialClass,
    AbaqusUserMaterialClass,
    AbaqusUserMaterialStatusClass,
    ConcreteDPM2Class,

    TransportMaterialClass,
    IsotropicHeatTransferMaterialClass,
    HeMoTKMaterialClass,

    HydratingIsoHeatMaterialClass,
    HydratingTransportMaterialStatusClass,
    HydratingHeMoMaterialClass,
    HydratingConcreteMatClass,

    HellmichMaterialClass,
    HellmichMaterialStatusClass,
    HydrationModelClass,
    HydrationModelStatusClass,

    CemhydMatClass,
    CemhydMatStatusClass,

    MaterialStatusClass,
    StructuralMaterialStatusClass,
    TransportMaterialStatusClass,
    PerfectlyPlasticMaterialStatusClass,
    SmearedCrackingMaterialStatusClass,
    PlasticSmearedCrackingMaterialStatusClass,
    Concrete2MaterialStatusClass,
    MaxwellChainMaterialStatusClass,
    KelvinChainMaterialStatusClass,
    KelvinChainSolidMaterialStatusClass,
    B3SolidMaterialStatusClass,
    MPSMaterialStatusClass,
    RCMMaterialStatusClass,
    RCSDMaterialClass,
    RCSDEMaterialClass,
    RCSDNLMaterialClass,
    M4MaterialStatusClass,
    IsotropicDamageMaterialStatusClass,
    CompoDamageMatStatusClass,
    MicroplaneDamageMaterialStatusClass,
    DruckerPragerPlasticitySMStatusClass,
    MPlasticMaterialStatusClass,
    MicroMaterialStatusClass,
    HyperElasticMaterialStatusClass,
    MisesMatStatusClass,
    RankineMatStatusClass,
    RheoChainMaterialStatusClass,
    TrabBoneMaterialStatusClass,
    TrabBone3DStatusClass,
    TrabBoneEmbedStatusClass,
    ConcreteDPMStatusClass,
    FE2SinteringMaterialStatusClass,
    ConcreteDPM2StatusClass,

    FluidDynamicMaterialClass,
    FluidDynamicMaterialStatusClass,
    NewtonianFluidMaterialClass,
    BinghamFluidMaterialClass,
    BinghamFluidMaterialStatusClass,
    TwoFluidMaterialClass,

    SurfaceTensionMaterialClass,
    SurfaceTensionMaterialStatusClass,

    SimpleInterfaceMaterialClass,
    SimpleInterfaceMaterialStatusClass,

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
    PatchIntegrationRuleClass,
    MicroplaneClass,

    FloatArrayClass,
    FloatMatrixClass,
    IntArrayClass,
    DictionaryClass,
    TDictionaryClass,
    TimeStepClass,
    MetaStepClass,

    TopologyDescriptionClass,
    MeshTopologyDescriptionClass,
    ParticleTopologyDescriptionClass,

    DofManExportModuleClass,

    ErrorEstimatorClass,
    ScalarErrorIndicatorClass,
    ZZErrorEstimatorClass,
    CombinedZZSIErrorEstimatorClass,
    HuertaErrorEstimatorClass,
    MeshQualityErrorEstimatorClass,

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
    InternalVariableFieldClass,

    LEPlicClass,
    LevelSetPCSClass,
    FastMarchingMethodClass,

    WallClockLoadBalancerMonitorClass,
    ParmetisLoadBalancerClass
    
    
};
} // end namespace oofem
#endif // classtype_h
