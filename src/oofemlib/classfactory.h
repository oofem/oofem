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
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef classfactory_h
#define classfactory_h

#include "oofemcfg.h"
#include "sparsemtrxtype.h"
#include "errorestimatortype.h"
#include "doftype.h"
#include "linsystsolvertype.h"
//#include "patch.h" // for PatchType
#include "nodalrecoverymodel.h" // for NodalRecoveryModelType
#include "integrationrule.h" // for IntegrationRuleType
#include "geneigvalsolvertype.h"
#include "materialmappingalgorithmtype.h"
#include "meshpackagetype.h"
#include "dofiditem.h"

#include <map>
#include <string>
#include <cstring>
#include <memory>

namespace oofem {
class Domain;
class EngngModel;
class Element;
class DofManager;
class GeneralBoundaryCondition;
class CrossSection;
class Material;
class Function;
class NonlocalBarrier;
class ExportModule;
class SparseNonLinearSystemNM;
class InitModule;
class TopologyDescription;
class Monitor;

class Dof;
class SparseMtrx;
class SparseLinearSystemNM;
class ErrorEstimator;
class InitialCondition;
class Patch;
class NodalRecoveryModel;
class SparseGeneralEigenValueSystemNM;
class IntegrationRule;
class MaterialMappingAlgorithm;
class MesherInterface;

class LoadBalancerMonitor;
class LoadBalancer;

class XfemManager;
class EnrichmentItem;
class NucleationCriterion;
class EnrichmentFunction;
class EnrichmentDomain;
class BasicGeometry;
class EnrichmentFront;
class PropagationLaw;

class FractureManager;
class FailureCriteriaStatus;
class FailureCriteria;

class ContactManager;
class ContactDefinition;

#ifdef _GNUC
#define OOFEM_ATTR_UNUSED __attribute__((unused))
#else
#define OOFEM_ATTR_UNUSED 
#endif

// Templates to wrap constructors into functions

// Wraps any ctor into 

/**
 * Wrapper for ctors, creates a function that returns a unique_ptr to the object. 
 * The base class (B) template argument must be specified explicitly, becuase we can't "downcast" the function pointer.
 * E.g.
 * CTOR<SparseLinearSystemNM, LDLTFactorization, Domain*, EngngModel*>
 * CTOR<SparseMtrx, PetscSparseMtrx>
 */
template< typename B, typename T, typename... V > std::unique_ptr<B> CTOR(V... x) { return std::make_unique<T>(x...); }

template< typename T > Dof *dofCreator(DofIDItem dofid, DofManager *dman) { return new T(dman, dofid); }

///@name Macros for registering new components. Unique dummy variables must be created as a result (design flaw in C++).
//@{
#define REGISTER_Element(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerElement(_IFT_ ## class ## _Name, CTOR< Element, class, int, Domain* > );
#define REGISTER_DofManager(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerDofManager(_IFT_ ## class ## _Name, CTOR< DofManager, class, int, Domain* > );
#define REGISTER_BoundaryCondition(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerBoundaryCondition(_IFT_ ## class ## _Name, CTOR< GeneralBoundaryCondition, class, int, Domain* > );
#define REGISTER_CrossSection(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerCrossSection(_IFT_ ## class ## _Name, CTOR< CrossSection, class, int, Domain* > );
#define REGISTER_Material(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerMaterial(_IFT_ ## class ## _Name, CTOR< Material, class, int, Domain* > );
#define REGISTER_Material_Alt(class, altname) static bool __dummy_ ## class ## altname OOFEM_ATTR_UNUSED = GiveClassFactory().registerMaterial(#altname, CTOR< Material, class, int, Domain* > );
#define REGISTER_EngngModel(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerEngngModel(_IFT_ ## class ## _Name, CTOR< EngngModel, class, int, EngngModel* > );
#define REGISTER_Function(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerFunction(_IFT_ ## class ## _Name, CTOR< Function, class, int, Domain* > );
#define REGISTER_NonlocalBarrier(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerNonlocalBarrier(_IFT_ ## class ## _Name, CTOR< NonlocalBarrier, class, int, Domain* > );
#define REGISTER_ExportModule(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerExportModule(_IFT_ ## class ## _Name, CTOR< ExportModule, class, int, EngngModel* > );
#define REGISTER_SparseNonLinearSystemNM(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerSparseNonLinearSystemNM(_IFT_ ## class ## _Name, CTOR< SparseNonLinearSystemNM, class, Domain*, EngngModel* > );
#define REGISTER_InitModule(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerInitModule(_IFT_ ## class ## _Name, CTOR< InitModule, class, int, EngngModel* > );
#define REGISTER_TopologyDescription(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerTopologyDescription(_IFT_ ## class ## _Name, CTOR< TopologyDescription, class, Domain* > );
#define REGISTER_LoadBalancerMonitor(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerLoadBalancerMonitor(_IFT_ ## class ## _Name, CTOR< LoadBalancerMonitor, class, EngngModel* > );
#define REGISTER_LoadBalancer(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerLoadBalancer(_IFT_ ## class ## _Name, CTOR< LoadBalancer, class, EngngModel* > );
#define REGISTER_Monitor(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerMonitor(_IFT_ ## class ## _Name, CTOR< Monitor, class, int > );

// These should be converted to use strings.
#define REGISTER_SparseMtrx(class, type) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerSparseMtrx(type, CTOR< SparseMtrx, class > );
#define REGISTER_SparseLinSolver(class, type) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerSparseLinSolver(type, CTOR< SparseLinearSystemNM, class, Domain*, EngngModel* > );
#define REGISTER_ErrorEstimator(class, type) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerErrorEstimator(type, CTOR< ErrorEstimator, class, int, Domain* > );
#define REGISTER_NodalRecoveryModel(class, type) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerNodalRecoveryModel(type, CTOR< NodalRecoveryModel, class, Domain* > );
#define REGISTER_GeneralizedEigenValueSolver(class, type) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerGeneralizedEigenValueSolver(type, CTOR< SparseGeneralEigenValueSystemNM, class, Domain*, EngngModel* > );
#define REGISTER_Mesher(class, type) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerMesherInterface(type, CTOR< MesherInterface, class, Domain* > );
#define REGISTER_MaterialMappingAlgorithm(class, type) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerMaterialMappingAlgorithm(type, CTOR< MaterialMappingAlgorithm, class > );

#define REGISTER_XfemManager(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerXfemManager(_IFT_ ## class ## _Name, CTOR< XfemManager, class, Domain* > );
#define REGISTER_EnrichmentItem(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerEnrichmentItem(_IFT_ ## class ## _Name, CTOR< EnrichmentItem, class, int, XfemManager*, Domain* > );
#define REGISTER_NucleationCriterion(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerNucleationCriterion(_IFT_ ## class ## _Name, CTOR< NucleationCriterion, class, Domain* > );
#define REGISTER_EnrichmentFunction(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerEnrichmentFunction(_IFT_ ## class ## _Name, CTOR< EnrichmentFunction, class, int, Domain* > );
//#define REGISTER_EnrichmentDomain(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerEnrichmentDomain(_IFT_ ## class ## _Name, CTOR< EnrichmentDomain, class > );
#define REGISTER_EnrichmentFront(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerEnrichmentFront(_IFT_ ## class ## _Name, CTOR< EnrichmentFront, class > );
#define REGISTER_PropagationLaw(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerPropagationLaw(_IFT_ ## class ## _Name, CTOR< PropagationLaw, class > );
#define REGISTER_Geometry(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerGeometry(_IFT_ ## class ## _Name, CTOR< BasicGeometry, class > );

#define REGISTER_FailureCriteria(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerFailureCriteria(_IFT_ ## class ## _Name, CTOR< FailureCriteria, class, int, FractureManager* > );
#define REGISTER_FailureCriteriaStatus(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerFailureCriteriaStatus(_IFT_ ## class ## _Name, CTOR< FailureCriteriaStatus, class, FailureCriteria* > );

#define REGISTER_ContactManager(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerContactManager(_IFT_ ## class ## _Name, CTOR< ContactManager, class, Domain* > );
#define REGISTER_ContactDefinition(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerContactDefinition(_IFT_ ## class ## _Name, CTOR< ContactDefinition, class, ContactManager* > );
///@todo What is this? Doesn't seem needed / Mikael
#define REGISTER_Quasicontinuum(class) static bool __dummy_ ## class OOFEM_ATTR_UNUSED = GiveClassFactory().registerQuasicontinuum(_IFT_ ## class ## _Name, < QuasiContinuum, class, ????? > );
//@}

/**
 * Class Factory allows to register terminal oofem classes, based on their membership
 * (classes representing elements, dof managers, material models, etc)
 * and create them on demand according to their name or id.
 * Global instance of ClassFactory named classFactory is created on startup.
 *
 * @note To register new elements on startup, you must call GiveClassFactory to ensure that the global class factory is created first. This is ensured if you use the corresponding macro.
 */
class OOFEM_EXPORT ClassFactory
{
private:
    /// Associative container containing element creators with element name as key.
    std :: map < std :: string, std::unique_ptr<Element> ( * )(int, Domain *) > elemList;
    /// Associative container containing dofmanager creators with dofmanager  name as key.
    std :: map < std :: string, std::unique_ptr<DofManager> ( * )(int, Domain *) > dofmanList;
    /// Associative container containing boundary condition creators with bc  name as key.
    std :: map < std :: string, std::unique_ptr<GeneralBoundaryCondition> ( * )(int, Domain *) > bcList;
    /// Associative container containing cross section creators with cross section name as key.
    std :: map < std :: string, std::unique_ptr<CrossSection> ( * )(int, Domain *) > csList;
    /// Associative container containing material creators with material name as key.
    std :: map < std :: string, std::unique_ptr<Material> ( * )(int, Domain *) > matList;
    /// Associative container containing engng model creators with engng model name as key.
    std :: map < std :: string, std::unique_ptr<EngngModel> ( * )(int, EngngModel *) > engngList;
    /// Associative container containing load time function creators with function name as key.
    std :: map < std :: string, std::unique_ptr<Function> ( * )(int, Domain *) > funcList;
    /// Associative container containing nonlocal barriers creators with barrier name as key.
    std :: map < std :: string, std::unique_ptr<NonlocalBarrier> ( * )(int, Domain *) > nlbList;
    /// Associative container containing export module creators.
    std :: map < std :: string, std::unique_ptr<ExportModule> ( * )(int, EngngModel *) > exportList;
    /// Associative container containing monitor creators.
    std :: map < std :: string, std::unique_ptr<Monitor> ( * )(int) > monitorList;
    /// Associative container containing nonlinear solver creators.
    std :: map < std :: string, std::unique_ptr<SparseNonLinearSystemNM> ( * )(Domain *, EngngModel *) > nonlinList;
    /// Associative container containing init module creators.
    std :: map < std :: string, std::unique_ptr<InitModule> ( * )(int, EngngModel *) > initList;
    /// Associative container containing topology description creators.
    std :: map < std :: string, std::unique_ptr<TopologyDescription> ( * )(Domain *) > topologyList;
    // Internal structures (accessed by hard-coded enum values)
    /// Associative container containing load balancer creators.
    std :: map < std :: string, std::unique_ptr<LoadBalancer> ( * )(Domain *) > loadBalancerList;
    /// Associative container containing load balancer monitor creators.
    std :: map < std :: string, std::unique_ptr<LoadBalancerMonitor> ( * )(EngngModel *) > loadMonitorList;
    /// Associative container containing sparse matrix creators.
    std :: map < SparseMtrxType, std::unique_ptr<SparseMtrx> ( * )() > sparseMtrxList;
    /// Associative container containing dof creators.
    std :: map < dofType, Dof * ( * )(DofIDItem, DofManager *) > dofList;
    /// Associative container containing error estimator creators.
    std :: map < ErrorEstimatorType, std::unique_ptr<ErrorEstimator> ( * )(int, Domain *) > errEstList;
    /// Associative container containing sparse linear solver creators
    std :: map < LinSystSolverType, std::unique_ptr<SparseLinearSystemNM> ( * )(Domain *, EngngModel *) > sparseLinSolList;
    /// Associative container containing nodal recovery model creators
    std :: map < NodalRecoveryModel :: NodalRecoveryModelType, std::unique_ptr<NodalRecoveryModel> ( * )(Domain *) > nodalRecoveryModelList;
    /// Associative container containing sparse generalized eigenvalue creators
    std :: map < GenEigvalSolverType, std::unique_ptr<SparseGeneralEigenValueSystemNM> ( * )(Domain *, EngngModel *) > generalizedEigenValueSolverList;
    /// Associative container containing material mapping algorithm creators
    std :: map < MaterialMappingAlgorithmType, std::unique_ptr<MaterialMappingAlgorithm> ( * )() > materialMappingList;
    /// Associative container containing mesher interface creators
    std :: map < MeshPackageType, std::unique_ptr<MesherInterface> ( * )(Domain *) > mesherInterfaceList;

    // XFEM:
    /// Associative container containing enrichment item creators
    std :: map < std :: string, std::unique_ptr<EnrichmentItem> ( * )(int, XfemManager *, Domain *) > enrichItemList;
    /// Associative container containing nucleation criterion creators
    std :: map < std :: string, std::unique_ptr<NucleationCriterion> ( * )(Domain *) > nucleationCritList;
    /// Associative container containing enrichment function creators
    std :: map < std :: string, std::unique_ptr<EnrichmentFunction> ( * )(int, Domain *) > enrichFuncList;
    /// Associative container containing geometry creators
    std :: map < std :: string, std::unique_ptr<BasicGeometry> ( * )() > geometryList;
    /// Associative container containing enrichment-domain creators
    std :: map < std :: string, std::unique_ptr<EnrichmentDomain> ( * )() > enrichmentDomainList;
    /// Associative container containing enrichment front creators
    std :: map < std :: string, std::unique_ptr<EnrichmentFront> ( * )() > enrichmentFrontList;
    /// Associative container containing propagation law creators
    std :: map < std :: string, std::unique_ptr<PropagationLaw> ( * )() > propagationLawList;
    /// Associative container containing XfemManager creators
    std :: map < std :: string, std::unique_ptr<XfemManager> ( * )(Domain *) > xManList;


    /// Associative container containing failure criteria creators
    std :: map < std :: string, std::unique_ptr<FailureCriteria> ( * )(int, FractureManager *) > failureCriteriaList;
    std :: map < std :: string, std::unique_ptr<FailureCriteriaStatus> ( * )(int, FailureCriteria *) > failureCriteriaStatusList;

    /// Associative container containing ContactManager creators
    std :: map < std :: string, std::unique_ptr<ContactManager> ( * )(Domain *) > contactManList;
    std :: map < std :: string, std::unique_ptr<ContactDefinition> ( * )(ContactManager *) > contactDefList;

public:
    /// Creates empty factory
    ClassFactory();

    /**
     * Creates new instance of element corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  Element number.
     * @param domain Domain assigned to new element.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    std::unique_ptr<Element> createElement(const char *name, int num, Domain *domain);
    /**
     * Registers a new element in the class factory.
     * @param name Keyword string.
     */
    bool registerElement( const char *name, std::unique_ptr<Element> ( *creator )( int, Domain * ) );
    /**
     * Creates new instance of Dof manager corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  dofmanager number.
     * @param domain Domain assigned to new instance.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    std::unique_ptr<DofManager> createDofManager(const char *name, int num, Domain *domain);
    /**
     * Registers a new dof manager in the class factory.
     * @param name Keyword string.
     */
    bool registerDofManager( const char *name, std::unique_ptr<DofManager> ( *creator )( int, Domain * ) );
    /**
     * Creates new instance of boundary condition corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  number of new object.
     * @param domain Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    std::unique_ptr<GeneralBoundaryCondition> createBoundaryCondition(const char *name, int num, Domain *domain);
    /**
     * Registers a new boundary condition in the class factory.
     * @param name Keyword string.
     */
    bool registerBoundaryCondition( const char *name, std::unique_ptr<GeneralBoundaryCondition> ( *creator )( int, Domain * ) );
    /**
     * Creates new instance of cross section corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  object's number.
     * @param domain Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    std::unique_ptr<CrossSection> createCrossSection(const char *name, int num, Domain *domain);
    /**
     * Registers a new cross section in the class factory.
     * @param name Keyword string.
     */
    bool registerCrossSection( const char *name, std::unique_ptr<CrossSection> ( *creator )( int, Domain * ) );
    /**
     * Creates new instance of material corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  material number.
     * @param domain Domain assigned to new material.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    std::unique_ptr<Material> createMaterial(const char *name, int num, Domain *domain);
    /**
     * Registers a new material in the class factory.
     * @param name Keyword string.
     */
    bool registerMaterial( const char *name, std::unique_ptr<Material> ( *creator )( int, Domain * ) );
    /**
     * Creates new instance of engng model corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  Engng model number.
     * @param master Master engineering model (used for multiscale modeling).
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    std::unique_ptr<EngngModel> createEngngModel(const char *name, int num, EngngModel *master);
    /**
     * Registers a new engineering model in the class factory.
     * @param name Keyword string.
     */
    bool registerEngngModel( const char *name, std::unique_ptr<EngngModel> ( *creator )( int, EngngModel * ) );
    /**
     * Creates new instance of load time function corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  object's number.
     * @param domain Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    std::unique_ptr<Function> createFunction(const char *name, int num, Domain *domain);
    /**
     * Registers a new load time function in the class factory.
     * @param name Keyword string.
     */
    bool registerFunction( const char *name, std::unique_ptr<Function> ( *creator )( int, Domain * ) );
    /**
     * Creates new instance of nonlocal barrier corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  object's number.
     * @param domain Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    std::unique_ptr<NonlocalBarrier> createNonlocalBarrier(const char *name, int num, Domain *domain);
    /**
     * Registers a new nonlocal barrier in the class factory.
     * @param name Keyword string.
     */
    bool registerNonlocalBarrier( const char *name, std::unique_ptr<NonlocalBarrier> ( *creator )( int, Domain * ) );
    /**
     * Creates new instance of export module corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  Export module number.
     * @param emodel Engineering model that object belongs to.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    std::unique_ptr<ExportModule> createExportModule(const char *name, int num, EngngModel *emodel);
    /**
     * Registers a new export module in the class factory.
     * @param name Keyword string.
     */
    bool registerExportModule( const char *name, std::unique_ptr<ExportModule> ( *creator )( int, EngngModel * ) );
    /**
     * Creates new instance of monitor corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  Monitor number.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    std::unique_ptr<Monitor> createMonitor(const char *name, int num);
    /**
     * Registers a new monitor in the class factory.
     * @param name Keyword string.
     */
    bool registerMonitor( const char *name, std::unique_ptr<Monitor> ( *creator )( int ) );

    /**
     * Creates new instance of nonlinear solver corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param domain Domain assigned to new object.
     * @param emodel Engineering model that object belongs to.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    std::unique_ptr<SparseNonLinearSystemNM> createNonLinearSolver(const char *name, Domain *domain, EngngModel *emodel);
    /**
     * Registers a new nonlinear solver in the class factory.
     * @param name Keyword string.
     */
    bool registerSparseNonLinearSystemNM( const char *name, std::unique_ptr<SparseNonLinearSystemNM> ( *creator )( Domain *, EngngModel * ) );
    /**
     * Creates new instance of init module corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  Init module number.
     * @param emodel Engineering model that object belongs to.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    std::unique_ptr<InitModule> createInitModule(const char *name, int num, EngngModel *emodel);
    /**
     * Registers a new init module in the class factory.
     * @param name Keyword string.
     */
    bool registerInitModule( const char *name, std::unique_ptr<InitModule> ( *creator )( int, EngngModel * ) );
    /**
     * Creates new instance of Initial Condition corresponding to given type.
     * @param name Keyword string determining the type of new instance.
     * @param num  object's number.
     * @param d Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    std::unique_ptr<InitialCondition> createInitialCondition(const char *name, int num, Domain *d);
    /**
     * Creates new instance of topology description corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param domain Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    std::unique_ptr<TopologyDescription> createTopology(const char *name, Domain *domain);
    /**
     * Registers a new topology description in the class factory.
     * @param name Keyword string.
     */
    bool registerTopologyDescription( const char *name, std::unique_ptr<TopologyDescription> ( *creator )( Domain * ) );

    /**
     * Creates new instance of sparse matrix corresponding to given keyword.
     * @param type SparseMtrxType id determining the type of new instance.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    std::unique_ptr<SparseMtrx> createSparseMtrx(SparseMtrxType type);
    /**
     * Registers a sparse matrix type.
     * @param type SparseMtrxType id determining the type of new instance.
     */
    bool registerSparseMtrx( SparseMtrxType type, std::unique_ptr<SparseMtrx> ( *creator )(void) );
    /**
     * Creates new instance of DOF corresponding to given keyword.
     * @param type ID determining the type of new instance.
     * @param dofid The dof ID.
     * @param dman Dof manager assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    Dof *createDof(dofType type, DofIDItem dofid, DofManager *dman);
    /**
     * Creates new instance of SparseLinearSystemNM corresponding
     * to given type.
     * @param st LinSystSolverType id determining the type of new instance.
     * @param d Domain assigned to new object.
     * @param m EngngModel  assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    std::unique_ptr<SparseLinearSystemNM> createSparseLinSolver(LinSystSolverType st, Domain *d, EngngModel *m);
    /**
     * Registers a sparse linear system solver.
     * @param type LinSystSolverType id determining the type of new instance.
     */
    bool registerSparseLinSolver( LinSystSolverType type, std::unique_ptr<SparseLinearSystemNM> ( *creator )(Domain *, EngngModel *) );
    /**
     * Creates new instance of ErrorEstimator corresponding
     * to given type.
     * @param type ErrorEstimatorType id determining the type of new instance.
     * @param num  object's number.
     * @param d Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    std::unique_ptr<ErrorEstimator> createErrorEstimator(ErrorEstimatorType type, int num, Domain *d);
    /**
     * Registers a new  error estimator.
     * @param type ErrorEstimatorType id determining the type of new instance.
     */
    bool registerErrorEstimator( ErrorEstimatorType type, std::unique_ptr<ErrorEstimator> ( *creator )(int, Domain *) );
    /**
     * Registers a new nodal recovery model.
     * @param name Indentifier.
     */
    bool registerNodalRecoveryModel( NodalRecoveryModel :: NodalRecoveryModelType name, std::unique_ptr<NodalRecoveryModel> ( *creator )(Domain *) );
    /**
     * Creates new instance of nodal recovery model corresponding to given type.
     * @param type ID determining the type of new instance.
     * @param d Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    std::unique_ptr<NodalRecoveryModel> createNodalRecoveryModel(NodalRecoveryModel :: NodalRecoveryModelType type, Domain *d);

    // XFEM:
    std::unique_ptr<EnrichmentItem> createEnrichmentItem(const char *name, int num, XfemManager *xm, Domain *domain);
    bool registerEnrichmentItem( const char *name, std::unique_ptr<EnrichmentItem> ( *creator )( int, XfemManager *, Domain * ) );

    std::unique_ptr<NucleationCriterion> createNucleationCriterion(const char *name, Domain *domain);
    bool registerNucleationCriterion( const char *name, std::unique_ptr<NucleationCriterion> ( *creator )( Domain * ) );

    std::unique_ptr<EnrichmentFunction> createEnrichmentFunction(const char *name, int num, Domain *domain);
    bool registerEnrichmentFunction( const char *name, std::unique_ptr<EnrichmentFunction> ( *creator )( int, Domain * ) );

#if 0
    std::unique_ptr<EnrichmentDomain> createEnrichmentDomain(const char *name);
    bool registerEnrichmentDomain( const char *name, std::unique_ptr<EnrichmentDomain> ( *creator )( ) );
#endif

    std::unique_ptr<EnrichmentFront> createEnrichmentFront(const char *name);
    bool registerEnrichmentFront( const char *name, std::unique_ptr<EnrichmentFront> ( *creator )( ) );

    std::unique_ptr<PropagationLaw> createPropagationLaw(const char *name);
    bool registerPropagationLaw( const char *name, std::unique_ptr<PropagationLaw> ( *creator )( ) );

    std::unique_ptr<BasicGeometry> createGeometry(const char *name);
    bool registerGeometry( const char *name, std::unique_ptr<BasicGeometry> ( *creator )( ) );

    std::unique_ptr<XfemManager> createXfemManager(const char *name, Domain *domain);
    bool registerXfemManager( const char *name, std::unique_ptr<XfemManager> ( *creator )( Domain * ) );


    std::unique_ptr<ContactManager> createContactManager(const char *name, Domain *domain);
    bool registerContactManager( const char *name, std::unique_ptr<ContactManager> ( *creator )( Domain * ) );

    std::unique_ptr<ContactDefinition> createContactDefinition(const char *name, ContactManager *cMan);
    bool registerContactDefinition( const char *name, std::unique_ptr<ContactDefinition> ( *creator )( ContactManager * ) );

    // Failure module (in development!)
    std::unique_ptr<FailureCriteria> createFailureCriteria(const char *name, int num, FractureManager *fracManager);
    bool registerFailureCriteria( const char *name, std::unique_ptr<FailureCriteria> ( *creator )( int, FractureManager * ) );

    std::unique_ptr<FailureCriteriaStatus> createFailureCriteriaStatus(const char *name, int num, FailureCriteria *critManager);
    bool registerFailureCriteriaStatus( const char *name, std::unique_ptr<FailureCriteriaStatus> ( *creator )( int, FailureCriteria * ) );

    std::unique_ptr<SparseGeneralEigenValueSystemNM> createGeneralizedEigenValueSolver(GenEigvalSolverType name, Domain *d, EngngModel *m);
    bool registerGeneralizedEigenValueSolver( GenEigvalSolverType name, std::unique_ptr<SparseGeneralEigenValueSystemNM> ( *creator )(Domain *, EngngModel *) );

    std::unique_ptr<IntegrationRule> createIRule(IntegrationRuleType name, int number, Element *e);

    std::unique_ptr<MaterialMappingAlgorithm> createMaterialMappingAlgorithm(MaterialMappingAlgorithmType name);
    bool registerMaterialMappingAlgorithm( MaterialMappingAlgorithmType name, std::unique_ptr<MaterialMappingAlgorithm> ( *creator )( ) );

    std::unique_ptr<MesherInterface> createMesherInterface(MeshPackageType name, Domain *d);
    bool registerMesherInterface( MeshPackageType name, std::unique_ptr<MesherInterface> ( *creator )( Domain * ) );

    std::unique_ptr<LoadBalancerMonitor> createLoadBalancerMonitor(const char *name, EngngModel *e);
    bool registerLoadBalancerMonitor( const char *name, std::unique_ptr<LoadBalancerMonitor> ( *creator )( EngngModel * ) );

    std::unique_ptr<LoadBalancer> createLoadBalancer(const char *name, Domain *d);
    bool registerLoadBalancer( const char *name, std::unique_ptr<LoadBalancer> ( *creator )( Domain * ) );
};

extern ClassFactory &classFactory;

/**
 * This function must be used by all code that run at link time to ensure that the classFactory is constructed first.
 * See "static initialization order fiasco" for explanation.
 */
ClassFactory &GiveClassFactory();
} // end namespace oofem
#endif // clasfactort_h
