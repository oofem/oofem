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

#include "classfactory.h"

#include <string>
#include <algorithm>
#include <cctype>

#include "masterdof.h"
#include "slavedof.h"
#include "simpleslavedof.h"
#include "activedof.h"

#include "gaussintegrationrule.h"
#include "lobattoir.h"

#include "initialcondition.h"


namespace oofem {
ClassFactory &GiveClassFactory()
{
    static ClassFactory ans;
    return ans;
}

ClassFactory &classFactory = GiveClassFactory();

std :: string conv2lower(std :: string input)
{
    std::transform(input.begin(), input.end(), input.begin(), ::tolower);
    return input;
}

// Non-string names (should be changed eventually):
template< typename C, typename T, typename V, typename ... As> C *cf_create2(const T &list, V name, As ... args)
{
    auto creator = list.find(name);
    return creator != list.end() ? creator->second(args...) : nullptr;
}

// Helper for storing creators
template< typename T, typename V, typename C> bool cf_store2(T &list, V name, C &creator)
{
    list[ name ] = creator;
    return true;
}

// Helper for creating objects
template< typename C, typename T, typename ... As> C *cf_create(const T &list, const char *name, As ... args)
{
    auto creator = list.find(conv2lower(name));
    return creator != list.end() ? creator->second(args...) : nullptr;
}

// Helper for storing creators
template< typename T, typename C> bool cf_store(T &list, const char *name, C &creator)
{
    list[ conv2lower(name) ] = creator;
    return true;
}


ClassFactory :: ClassFactory()
{
    // Fixed list for DOF types. No new components can register for these since these are part of the internal structure in OOFEM.
    dofList [ DT_master ] = dofCreator< MasterDof >;
    dofList [ DT_simpleSlave ] = dofCreator< SimpleSlaveDof >;
    dofList [ DT_slave ] = dofCreator< SlaveDof >;
    dofList [ DT_active ] = dofCreator< ActiveDof >;
}

SparseMtrx *ClassFactory :: createSparseMtrx(SparseMtrxType name)
{
    return cf_create2<SparseMtrx>(sparseMtrxList, name);
}

bool ClassFactory :: registerSparseMtrx( SparseMtrxType name, SparseMtrx * ( *creator )( ) )
{
    return cf_store2(sparseMtrxList, name, creator);
}

Dof *ClassFactory :: createDof(dofType name, DofIDItem dofid, DofManager *dman)
{
    return cf_create2<Dof>(dofList, name, dofid, dman);
}

SparseLinearSystemNM *ClassFactory :: createSparseLinSolver(LinSystSolverType name, Domain *domain, EngngModel *emodel)
{
    return cf_create2<SparseLinearSystemNM>(sparseLinSolList, name, domain, emodel);
}

bool ClassFactory :: registerSparseLinSolver( LinSystSolverType name, SparseLinearSystemNM * ( *creator )( Domain *, EngngModel * ) )
{
    return cf_store2(sparseLinSolList, name, creator);
}

ErrorEstimator *ClassFactory :: createErrorEstimator(ErrorEstimatorType name, int number, Domain *domain)
{
    return cf_create2<ErrorEstimator>(errEstList, name, number, domain);
}

bool ClassFactory :: registerErrorEstimator( ErrorEstimatorType name, ErrorEstimator * ( *creator )( int, Domain * ) )
{
    return cf_store2(errEstList, name, creator);
}

InitialCondition *ClassFactory :: createInitialCondition(const char *name, int number, Domain *domain)
{
    if ( conv2lower(name).compare("initialcondition") == 0 ) {
        return new InitialCondition(number, domain);
    }
    return nullptr;
}

bool ClassFactory :: registerNodalRecoveryModel( NodalRecoveryModel :: NodalRecoveryModelType name, NodalRecoveryModel * ( *creator )(Domain *) )
{
    return cf_store2(nodalRecoveryModelList, name, creator);
}

NodalRecoveryModel *ClassFactory :: createNodalRecoveryModel(NodalRecoveryModel :: NodalRecoveryModelType name, Domain *domain)
{
    return cf_create2<NodalRecoveryModel>(nodalRecoveryModelList, name, domain);
}


Element *ClassFactory :: createElement(const char *name, int number, Domain *domain)
{
    return cf_create<Element>(elemList, name,number, domain);
}

bool ClassFactory :: registerElement( const char *name, Element * ( *creator )( int, Domain * ) )
{
    return cf_store(elemList, name, creator);
}

DofManager *ClassFactory :: createDofManager(const char *name, int number, Domain *domain)
{
    return cf_create<DofManager>(dofmanList, name, number, domain);
}

bool ClassFactory :: registerDofManager( const char *name, DofManager * ( *creator )( int, Domain * ) )
{
    return cf_store(dofmanList, name, creator);
}

GeneralBoundaryCondition *ClassFactory :: createBoundaryCondition(const char *name, int number, Domain *domain)
{
    return cf_create<GeneralBoundaryCondition>(bcList, name, number, domain);
}

bool ClassFactory :: registerBoundaryCondition( const char *name, GeneralBoundaryCondition * ( *creator )( int, Domain * ) )
{
    return cf_store(bcList, name, creator);
}

CrossSection *ClassFactory :: createCrossSection(const char *name, int number, Domain *domain)
{
    return cf_create<CrossSection>(csList, name, number, domain);
}

bool ClassFactory :: registerCrossSection( const char *name, CrossSection * ( *creator )( int, Domain * ) )
{
    return cf_store(csList, name, creator);
}

Material *ClassFactory :: createMaterial(const char *name, int number, Domain *domain)
{
    return cf_create<Material>(matList, name, number, domain);
}

bool ClassFactory :: registerMaterial( const char *name, Material * ( *creator )( int, Domain * ) )
{
    return cf_store(matList, name, creator);
}

EngngModel *ClassFactory :: createEngngModel(const char *name, int number, EngngModel *master)
{
    return cf_create<EngngModel>(engngList, name, number, master);
}

bool ClassFactory :: registerEngngModel( const char *name, EngngModel * ( *creator )( int, EngngModel * ) )
{
    return cf_store(engngList, name, creator);
}

Function *ClassFactory :: createFunction(const char *name, int number, Domain *domain)
{
    return cf_create<Function>(funcList, name, number, domain);
}

bool ClassFactory :: registerFunction( const char *name, Function * ( *creator )( int, Domain * ) )
{
    return cf_store(funcList, name, creator);
}

NonlocalBarrier *ClassFactory :: createNonlocalBarrier(const char *name, int number, Domain *domain)
{
    return cf_create<NonlocalBarrier>(nlbList, name, number, domain);
}

bool ClassFactory :: registerNonlocalBarrier( const char *name, NonlocalBarrier * ( *creator )( int, Domain * ) )
{
    return cf_store(nlbList, name, creator);
}

ExportModule *ClassFactory :: createExportModule(const char *name, int number, EngngModel *emodel)
{
    return cf_create<ExportModule>(exportList, name, number, emodel);
}

bool ClassFactory :: registerExportModule( const char *name, ExportModule * ( *creator )( int, EngngModel * ) )
{
    return cf_store(exportList, name, creator);
}

SparseNonLinearSystemNM *ClassFactory :: createNonLinearSolver(const char *name, Domain *domain, EngngModel *emodel)
{
    return cf_create<SparseNonLinearSystemNM>(nonlinList, name, domain, emodel);
}

bool ClassFactory :: registerSparseNonLinearSystemNM(const char *name, SparseNonLinearSystemNM * ( *creator )( Domain *, EngngModel * ) )
{
    return cf_store(nonlinList, name, creator);
}

InitModule *ClassFactory :: createInitModule(const char *name, int number, EngngModel *emodel)
{
    return cf_create<InitModule>(initList, name, number, emodel);
}

bool ClassFactory :: registerInitModule(const char *name, InitModule * ( *creator )( int, EngngModel * ) )
{
    return cf_store(initList, name, creator);
}

TopologyDescription *ClassFactory :: createTopology(const char *name, Domain *domain)
{
    return cf_create<TopologyDescription>(topologyList, name, domain);
}

bool ClassFactory :: registerTopologyDescription(const char *name, TopologyDescription * ( *creator )( Domain * ) )
{
    return cf_store(topologyList, name, creator);
}


// XFEM:
EnrichmentItem *ClassFactory :: createEnrichmentItem(const char *name, int number, XfemManager *xm, Domain *domain)
{
    return cf_create<EnrichmentItem>(enrichItemList, name, number, xm, domain);
}

bool ClassFactory :: registerEnrichmentItem(const char *name, EnrichmentItem * ( *creator )( int, XfemManager *, Domain * ) )
{
    return cf_store(enrichItemList, name, creator);
}

NucleationCriterion *ClassFactory :: createNucleationCriterion(const char *name, Domain *domain)
{
    return cf_create<NucleationCriterion>(nucleationCritList, name, domain);
}

bool ClassFactory :: registerNucleationCriterion(const char *name, NucleationCriterion * ( *creator )( Domain * ) )
{
    return cf_store(nucleationCritList, name, creator);
}

EnrichmentFunction *ClassFactory :: createEnrichmentFunction(const char *name, int number, Domain *domain)
{
    return cf_create<EnrichmentFunction>(enrichFuncList, name, number, domain);
}

bool ClassFactory :: registerEnrichmentFunction(const char *name, EnrichmentFunction * ( *creator )( int, Domain * ) )
{
    return cf_store(enrichFuncList, name, creator);
}

EnrichmentDomain *ClassFactory :: createEnrichmentDomain(const char *name)
{
    return cf_create<EnrichmentDomain>(enrichmentDomainList, name);
}

bool ClassFactory :: registerEnrichmentDomain(const char *name, EnrichmentDomain * ( *creator )( ) )
{
    return cf_store(enrichmentDomainList, name, creator);
}

EnrichmentFront *ClassFactory :: createEnrichmentFront(const char *name)
{
    return cf_create<EnrichmentFront>(enrichmentFrontList, name);
}

bool ClassFactory :: registerEnrichmentFront(const char *name, EnrichmentFront * ( *creator )( ) )
{
    return cf_store(enrichmentFrontList, name, creator);
}

PropagationLaw *ClassFactory :: createPropagationLaw(const char *name)
{
    return cf_create<PropagationLaw>(propagationLawList, name);
}

bool ClassFactory :: registerPropagationLaw( const char *name, PropagationLaw * ( *creator )( ) )
{
    return cf_store(propagationLawList, name, creator);
}

BasicGeometry *ClassFactory :: createGeometry(const char *name)
{
    return cf_create<BasicGeometry>(geometryList, name);
}

bool ClassFactory :: registerGeometry(const char *name, BasicGeometry * ( *creator )( ) )
{
    return cf_store(geometryList, name, creator);
}

XfemManager *ClassFactory :: createXfemManager(const char *name, Domain *domain)
{
    return cf_create<XfemManager>(xManList, name, domain);
}

bool ClassFactory :: registerXfemManager(const char *name, XfemManager * ( *creator )( Domain * ) )
{
    return cf_store(xManList, name, creator);
}


// Failure module:

FailureCriteria *ClassFactory :: createFailureCriteria(const char *name, int number, FractureManager *fracManager)
{
    return cf_create<FailureCriteria>(failureCriteriaList, name, number, fracManager);
}

bool ClassFactory :: registerFailureCriteria( const char *name, FailureCriteria * ( *creator )( int, FractureManager * ) )
{
    return cf_store(failureCriteriaList, name, creator);
}

FailureCriteriaStatus *ClassFactory :: createFailureCriteriaStatus(const char *name, int number, FailureCriteria *fc)
{
    return cf_create<FailureCriteriaStatus>(failureCriteriaStatusList, name, number, fc);
}

bool ClassFactory :: registerFailureCriteriaStatus(const char *name, FailureCriteriaStatus * ( *creator )( int, FailureCriteria * ) )
{
    return cf_store(failureCriteriaStatusList, name, creator);
}


ContactManager *ClassFactory :: createContactManager(const char *name, Domain *domain)
{
    return cf_create<ContactManager>(contactManList, name, domain);
}

bool ClassFactory :: registerContactManager(const char *name, ContactManager * ( *creator )( Domain * ) )
{
    return cf_store(contactManList, name, creator);
}


ContactDefinition *ClassFactory :: createContactDefinition(const char *name, ContactManager *cMan)
{
    return cf_create<ContactDefinition>(contactDefList, name, cMan);
}

bool ClassFactory :: registerContactDefinition( const char *name, ContactDefinition * ( *creator )( ContactManager * ) )
{
    return cf_store(contactDefList, name, creator);
}

bool ClassFactory :: registerGeneralizedEigenValueSolver( GenEigvalSolverType name, SparseGeneralEigenValueSystemNM * ( *creator )(Domain *, EngngModel *) )
{
    return cf_store2(generalizedEigenValueSolverList, name, creator);
}

SparseGeneralEigenValueSystemNM *ClassFactory :: createGeneralizedEigenValueSolver(GenEigvalSolverType name, Domain *domain, EngngModel *emodel)
{
    return cf_create2<SparseGeneralEigenValueSystemNM>(generalizedEigenValueSolverList, name, domain, emodel);
}

IntegrationRule *ClassFactory :: createIRule(IntegrationRuleType type, int number, Element *e)
{
    if ( type == IRT_Gauss ) {
        return new GaussIntegrationRule(number, e);
    } else if ( type == IRT_Lobatto ) {
        return new LobattoIntegrationRule(number, e);
    }
    return nullptr;
}

bool ClassFactory :: registerMaterialMappingAlgorithm( MaterialMappingAlgorithmType name, MaterialMappingAlgorithm * ( *creator )( ) )
{
    return cf_store2(materialMappingList, name, creator);
}

MaterialMappingAlgorithm *ClassFactory :: createMaterialMappingAlgorithm(MaterialMappingAlgorithmType name)
{
    return cf_create2<MaterialMappingAlgorithm>(materialMappingList, name);
}

bool ClassFactory :: registerMesherInterface( MeshPackageType name, MesherInterface * ( *creator )( Domain * ) )
{
    return cf_store2(mesherInterfaceList, name, creator);
}

MesherInterface *ClassFactory :: createMesherInterface(MeshPackageType name, Domain *domain)
{
    return cf_create2<MesherInterface>(mesherInterfaceList, name, domain);
}


LoadBalancerMonitor *ClassFactory :: createLoadBalancerMonitor(const char *name, EngngModel *emodel)
{
    return cf_create<LoadBalancerMonitor>(loadMonitorList, name, emodel);
}

bool ClassFactory :: registerLoadBalancerMonitor(const char *name, LoadBalancerMonitor * ( *creator )( EngngModel * ) )
{
    return cf_store(loadMonitorList, name, creator);
}

LoadBalancer *ClassFactory :: createLoadBalancer(const char *name, Domain *domain)
{
    return cf_create<LoadBalancer>(loadBalancerList, name, domain);
}

bool ClassFactory :: registerLoadBalancer(const char *name, LoadBalancer * ( *creator )( Domain * ) )
{
    return cf_store(loadBalancerList, name, creator);
}

} // End namespace oofem
