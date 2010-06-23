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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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
 * owerload the abstract giveClassID service. This service can be used
 * to distinguish between multiple derived classes.
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
    DofManagerClass,
    NodeClass,
    RigidArmNodeClass,
    HangingNodeClass,
    ElementSideClass,
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

    TrStructuralElementClass,
    NLStructuralElementClass,
    PlaneStress2dClass,
    Quad1PlaneStrainClass,
    QPlaneStress2dClass,
    TrPlaneStress2dClass,
    CCTPlateClass,
    CCTPlate3dClass,
    LSpaceClass,
    LSpaceBBClass,
    QSpaceClass,
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
    TrPlaneStrainClass,
    LIBeam3dNLClass,
    Quad1_htClass,
    Tr1_htClass,
    QuadAxisym1_htClass,
    TrAxisym1_htClass,
    Brick1_htClass,
    Tetrah1_htClass,
    PlaneStress2dXfemClass,

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
    MaxwellChainMaterialClass,
    KelvinChainMaterialClass,
    CebFip78MaterialClass,
    DoublePowerLawMaterialClass,
    B3MaterialClass,
    B3SolidMaterialClass,
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
    IsotropicDamageMaterial2Class,
    MazarsMaterialClass,
    CompositeDamageMaterialClass,
    MicroplaneDamageMaterialClass,
    DruckerPragerPlasticitySMClass,
    MPlasticMaterialClass,
    J2MatClass,
    HyperElasticMaterialClass,
    MisesMatClass,
    RheoChainMaterialClass,
    TrabBoneMaterialClass,
    TrabBone3DClass,
    TrabBoneEmbedClass,
    ConcreteDPMClass,

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
    KelvinChainMaterialStatusClass,
    B3SolidMaterialStatusClass,
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
    RheoChainMaterialStatusClass,
    TrabBoneMaterialStatusClass,
    TrabBone3DStatusClass,
    TrabBoneEmbedStatusClass,
    ConcreteDPMStatusClass,

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
    PatchIntegrationRuleClass,
    MicroplaneClass,

    FloatArrayClass,
    FloatMatrixClass,
    IntArrayClass,
    DictionaryClass,
    TDictionaryClass,
    TimeStepClass,
    MetaStepClass,

    DofManExportModuleClass,

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
    InternalVariableFieldClass,

    LEPlicClass,
    LevelSetPCSClass,
    FastMarchingMethodClass,

    WallClockLoadBalancerMonitorClass,
    ParmetisLoadBalancerClass
};

} // end namespace oofem
#endif // classtype_h
