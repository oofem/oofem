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
#include <cctype>

#include "masterdof.h"
#include "slavedof.h"
#include "simpleslavedof.h"
#include "activedof.h"

#include "gaussintegrationrule.h"
#include "lobattoir.h"

#include "initialcondition.h"

#include "subspaceit.h"
#include "inverseit.h"
#ifdef __SLEPC_MODULE
 #include "slepcsolver.h"
#endif

#include "nodalaveragingrecoverymodel.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"

// mesher interfaces
#include "t3dinterface.h"
#include "targe2interface.h"
#include "subdivision.h"
#include "freeminterface.h"

#include "mmaclosestiptransfer.h"
#include "mmaleastsquareprojection.h"
#include "mmashapefunctprojection.h"
#include "mmacontainingelementprojection.h" ///@todo This doesn't seem to be included? It is broken?

// Helper macro for creating objects
#define CF_CREATE(list, ...) \
    auto creator = list.find(conv2lower(name));\
    return creator != list.end() ? creator->second(__VA_ARGS__) : nullptr;

// Helper macro for storing objects
#define CF_STORE(list) \
    list[ conv2lower(name) ] = creator;\
    return true;

namespace oofem {
ClassFactory &GiveClassFactory()
{
    static ClassFactory ans;
    return ans;
}

ClassFactory &classFactory = GiveClassFactory();

std :: string conv2lower(std :: string input)
{
    for ( std :: size_t i = 0; i < input.size(); i++ ) {
        input [ i ] = (char)std :: tolower(input [ i ]);
    }
    return input;
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
    //CF_CREATE(sparseMtrxList)
    return ( sparseMtrxList.count(name) == 1 ) ? sparseMtrxList [ name ]() : NULL;
}

bool ClassFactory :: registerSparseMtrx( SparseMtrxType name, SparseMtrx * ( *creator )( ) )
{
    //CF_STORE(sparseMtrxList);
    sparseMtrxList[ name ] = creator;
    return true;
}

Dof *ClassFactory :: createDof(dofType name, DofIDItem dofid, DofManager *dman)
{
    //CF_CREATE(dofList, dofid, dman)
    return ( dofList.count(name) == 1 ) ? dofList [ name ](dofid, dman) : NULL;
}

SparseLinearSystemNM *ClassFactory :: createSparseLinSolver(LinSystSolverType name, Domain *domain, EngngModel *emodel)
{
    //CF_CREATE(sparseLinSolList, domain, emodel)
    return ( sparseLinSolList.count(name) == 1 ) ? sparseLinSolList [ name ](domain, emodel) : NULL;
}

bool ClassFactory :: registerSparseLinSolver( LinSystSolverType name, SparseLinearSystemNM * ( *creator )( Domain *, EngngModel * ) )
{
    //CF_STORE(sparseLinSolList)
    sparseLinSolList[ name ] = creator;
    return true;
}

ErrorEstimator *ClassFactory :: createErrorEstimator(ErrorEstimatorType name, int number, Domain *domain)
{
    //CF_CREATE(errEstList, number, domain)
    return ( errEstList.count(name) == 1 ) ? errEstList [ name ](number, domain) : NULL;
}

bool ClassFactory :: registerErrorEstimator( ErrorEstimatorType name, ErrorEstimator * ( *creator )( int, Domain * ) )
{
    //CF_STORE(errEstList)
    errEstList[ name ] = creator;
    return true;
}

InitialCondition *ClassFactory :: createInitialCondition(const char *name, int number, Domain *domain)
{
    if ( conv2lower(name).compare("initialcondition") == 0 ) {
        return new InitialCondition(number, domain);
    }
    return NULL;
}

NodalRecoveryModel *ClassFactory :: createNodalRecoveryModel(NodalRecoveryModel :: NodalRecoveryModelType type, Domain *domain)
{
    //CF_CREATE(nodalRecoveryModelList, domain)
    if ( type == NodalRecoveryModel :: NRM_NodalAveraging ) {
        return new NodalAveragingRecoveryModel(domain);
    } else if ( type == NodalRecoveryModel :: NRM_ZienkiewiczZhu ) {
        return new ZZNodalRecoveryModel(domain);
    } else if ( type == NodalRecoveryModel :: NRM_SPR ) {
        return new SPRNodalRecoveryModel(domain);
    }
    return NULL;
}


Element *ClassFactory :: createElement(const char *name, int number, Domain *domain)
{
    CF_CREATE(elemList, number, domain)
}

bool ClassFactory :: registerElement( const char *name, Element * ( *creator )( int, Domain * ) )
{
    CF_STORE(elemList)
}

DofManager *ClassFactory :: createDofManager(const char *name, int number, Domain *domain)
{
    CF_CREATE(dofmanList, number, domain)
}

bool ClassFactory :: registerDofManager( const char *name, DofManager * ( *creator )( int, Domain * ) )
{
    CF_STORE(dofmanList)
}

GeneralBoundaryCondition *ClassFactory :: createBoundaryCondition(const char *name, int number, Domain *domain)
{
    CF_CREATE(bcList, number, domain)
}

bool ClassFactory :: registerBoundaryCondition( const char *name, GeneralBoundaryCondition * ( *creator )( int, Domain * ) )
{
    CF_STORE(bcList)
}

CrossSection *ClassFactory :: createCrossSection(const char *name, int number, Domain *domain)
{
    CF_CREATE(csList, number, domain)
}

bool ClassFactory :: registerCrossSection( const char *name, CrossSection * ( *creator )( int, Domain * ) )
{
    CF_STORE(csList)
}

Material *ClassFactory :: createMaterial(const char *name, int number, Domain *domain)
{
    CF_CREATE(matList, number, domain)
}

bool ClassFactory :: registerMaterial( const char *name, Material * ( *creator )( int, Domain * ) )
{
    CF_STORE(matList)
}

EngngModel *ClassFactory :: createEngngModel(const char *name, int number, EngngModel *master)
{
    CF_CREATE(engngList, number, master)
}

bool ClassFactory :: registerEngngModel( const char *name, EngngModel * ( *creator )( int, EngngModel * ) )
{
    CF_STORE(engngList)
}

Function *ClassFactory :: createFunction(const char *name, int number, Domain *domain)
{
    CF_CREATE(funcList, number, domain)
}

bool ClassFactory :: registerFunction( const char *name, Function * ( *creator )( int, Domain * ) )
{
    CF_STORE(funcList)
}

NonlocalBarrier *ClassFactory :: createNonlocalBarrier(const char *name, int number, Domain *domain)
{
    CF_CREATE(nlbList, number, domain)
}

bool ClassFactory :: registerNonlocalBarrier( const char *name, NonlocalBarrier * ( *creator )( int, Domain * ) )
{
    CF_STORE(nlbList)
}

ExportModule *ClassFactory :: createExportModule(const char *name, int number, EngngModel *emodel)
{
    CF_CREATE(exportList, number, emodel)
}

bool ClassFactory :: registerExportModule( const char *name, ExportModule * ( *creator )( int, EngngModel * ) )
{
    CF_STORE(exportList)
}

SparseNonLinearSystemNM *ClassFactory :: createNonLinearSolver(const char *name, Domain *domain, EngngModel *emodel)
{
    CF_CREATE(nonlinList, domain, emodel)
}

bool ClassFactory :: registerSparseNonLinearSystemNM( const char *name, SparseNonLinearSystemNM * ( *creator )( Domain *, EngngModel * ) )
{
    CF_STORE(nonlinList)
}

InitModule *ClassFactory :: createInitModule(const char *name, int number, EngngModel *emodel)
{
    CF_CREATE(initList, number, emodel)
}

bool ClassFactory :: registerInitModule( const char *name, InitModule * ( *creator )( int, EngngModel * ) )
{
    CF_STORE(initList)
}

TopologyDescription *ClassFactory :: createTopology(const char *name, Domain *domain)
{
    CF_CREATE(topologyList, domain)
}

bool ClassFactory :: registerTopologyDescription( const char *name, TopologyDescription * ( *creator )( Domain * ) )
{
    CF_STORE(topologyList)
}


// XFEM:
EnrichmentItem *ClassFactory :: createEnrichmentItem(const char *name, int number, XfemManager *xm, Domain *domain)
{
    CF_CREATE(enrichItemList, number, xm, domain)
}

bool ClassFactory :: registerEnrichmentItem( const char *name, EnrichmentItem * ( *creator )( int, XfemManager *, Domain * ) )
{
    CF_STORE(enrichItemList)
}

EnrichmentFunction *ClassFactory :: createEnrichmentFunction(const char *name, int number, Domain *domain)
{
    CF_CREATE(enrichFuncList, number, domain)
}

bool ClassFactory :: registerEnrichmentFunction( const char *name, EnrichmentFunction * ( *creator )( int, Domain * ) )
{
    CF_STORE(enrichFuncList)
}

EnrichmentDomain *ClassFactory :: createEnrichmentDomain(const char *name)
{
    CF_CREATE(enrichmentDomainList)
}

bool ClassFactory :: registerEnrichmentDomain( const char *name, EnrichmentDomain * ( *creator )( ) )
{
    CF_STORE(enrichmentDomainList)
}

EnrichmentFront *ClassFactory :: createEnrichmentFront(const char *name)
{
    CF_CREATE(enrichmentFrontList)
}

bool ClassFactory :: registerEnrichmentFront( const char *name, EnrichmentFront * ( *creator )( ) )
{
    CF_STORE(enrichmentFrontList)
}

PropagationLaw *ClassFactory :: createPropagationLaw(const char *name)
{
    CF_CREATE(propagationLawList)
}

bool ClassFactory :: registerPropagationLaw( const char *name, PropagationLaw * ( *creator )( ) )
{
    CF_STORE(propagationLawList)
}

BasicGeometry *ClassFactory :: createGeometry(const char *name)
{
    CF_CREATE(geometryList)
}

bool ClassFactory :: registerGeometry( const char *name, BasicGeometry * ( *creator )( ) )
{
    CF_STORE(geometryList)
}

XfemManager *ClassFactory :: createXfemManager(const char *name, Domain *domain)
{
    CF_CREATE(xManList, domain)
}

bool ClassFactory :: registerXfemManager( const char *name, XfemManager * ( *creator )( Domain * ) )
{
    CF_STORE(xManList)
}


// Failure module:

FailureCriteria *ClassFactory :: createFailureCriteria(const char *name, int number, FractureManager *fracManager)
{
    CF_CREATE(failureCriteriaList, number, fracManager)
}

bool ClassFactory :: registerFailureCriteria( const char *name, FailureCriteria * ( *creator )( int, FractureManager * ) )
{
    CF_STORE(failureCriteriaList)
}

FailureCriteriaStatus *ClassFactory :: createFailureCriteriaStatus(const char *name, int number, FailureCriteria *fc)
{
    CF_CREATE(failureCriteriaStatusList, number, fc)
}

bool ClassFactory :: registerFailureCriteriaStatus( const char *name, FailureCriteriaStatus * ( *creator )( int, FailureCriteria * ) )
{
    CF_STORE(failureCriteriaStatusList)
}


ContactManager *ClassFactory :: createContactManager(const char *name, Domain *domain)
{
    CF_CREATE(contactManList, domain)
}

bool ClassFactory :: registerContactManager( const char *name, ContactManager * ( *creator )( Domain * ) )
{
    CF_STORE(contactManList)
}


ContactDefinition *ClassFactory :: createContactDefinition(const char *name, ContactManager *cMan)
{
    CF_CREATE(contactDefList, cMan)
}

bool ClassFactory :: registerContactDefinition( const char *name, ContactDefinition * ( *creator )( ContactManager * ) )
{
    CF_STORE(contactDefList)
}



SparseGeneralEigenValueSystemNM *ClassFactory :: createGeneralizedEigenValueSolver(GenEigvalSolverType st, Domain *domain, EngngModel *emodel)
{
    //CF_CREATE(mesherInterfaceList, domain)
    if ( st == GES_SubspaceIt ) {
        return new SubspaceIteration(domain, emodel);
    } else if ( st == GES_InverseIt ) {
        return new InverseIteration(domain, emodel);
    }
#ifdef __SLEPC_MODULE
    else if ( st == GES_SLEPc ) {
        return new SLEPcSolver(domain, emodel);
    }
#endif
    return NULL;
}

IntegrationRule *ClassFactory :: createIRule(IntegrationRuleType type, int number, Element *e)
{
    if ( type == IRT_Gauss ) {
        return new GaussIntegrationRule(number, e);
    } else if ( type == IRT_Lobatto ) {
        return new LobattoIntegrationRule(number, e);
    }
    return NULL;
}

MaterialMappingAlgorithm *ClassFactory :: createMaterialMappingAlgorithm(MaterialMappingAlgorithmType type)
{
    //CF_CREATE(materialMappingList)
    if ( type == MMA_ClosestPoint ) {
        return new MMAClosestIPTransfer();
    } else if ( type == MMA_LeastSquareProjection ) {
        return new MMALeastSquareProjection();
    } else if ( type == MMA_ShapeFunctionProjection ) {
        return new MMAShapeFunctProjection();
    }
    return NULL;
}

MesherInterface *ClassFactory :: createMesherInterface(MeshPackageType type, Domain *domain)
{
    //CF_CREATE(mesherInterfaceList, domain)
    if ( type == MPT_T3D ) {
        return new T3DInterface(domain);
    } else if ( type == MPT_TARGE2 ) {
        return new Targe2Interface(domain);
    } else if ( type == MPT_FREEM ) {
        return new FreemInterface(domain);
    } else if ( type == MPT_SUBDIVISION ) {
        return new Subdivision(domain);
    }
    return NULL;
}


LoadBalancerMonitor *ClassFactory :: createLoadBalancerMonitor(const char *name, EngngModel *emodel)
{
    CF_CREATE(loadMonitorList, emodel)
}

bool ClassFactory :: registerLoadBalancerMonitor( const char *name, LoadBalancerMonitor * ( *creator )( EngngModel * ) )
{
    CF_STORE(loadMonitorList)
}

LoadBalancer *ClassFactory :: createLoadBalancer(const char *name, Domain *domain)
{
    CF_CREATE(loadBalancerList, domain)
}

bool ClassFactory :: registerLoadBalancer( const char *name, LoadBalancer * ( *creator )( Domain * ) )
{
    CF_STORE(loadBalancerList)
}

} // End namespace oofem
