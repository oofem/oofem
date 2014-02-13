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

#include <cstring>
#ifdef HAVE_STRINGS_H
 #include <strings.h>
#endif
#include "compiler.h"

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

#ifdef __PARALLEL_MODE
 #include "loadbalancer.h"
 #include "parmetisloadbalancer.h"
#endif

namespace oofem {
ClassFactory &GiveClassFactory()
{
    static ClassFactory ans;
    return ans;
}

ClassFactory &classFactory = GiveClassFactory();

int ClassFactory :: CaseComp :: operator() (const std :: string & a, const std :: string & b) const
{
    return strncasecmp( a.c_str(), b.c_str(), ( int ) b.length() ) < 0;
}


ClassFactory :: ClassFactory()
{
    // Fixed list for DOF types. No new components can register for these since these are part of the internal structure in OOFEM.
    dofList [ DT_master ] = dofCreator< MasterDof >;
    dofList [ DT_simpleSlave ] = dofCreator< SimpleSlaveDof >;
    dofList [ DT_slave ] = dofCreator< SlaveDof >;
    dofList [ DT_active ] = dofCreator< ActiveDof >;

#ifdef __PARALLEL_MODE
    loadBalancerList [ "parmetis" ] = loadBalancerCreator< ParmetisLoadBalancer >;
    loadMonitorList [ "wallclock" ] = loadMonitorCreator< WallClockLoadBalancerMonitor >;
#endif
}

SparseMtrx *ClassFactory :: createSparseMtrx(SparseMtrxType name)
{
    return ( sparseMtrxList.count(name) == 1 ) ? sparseMtrxList [ name ]() : NULL;
}

bool ClassFactory :: registerSparseMtrx( SparseMtrxType name, SparseMtrx * ( *creator )( ) )
{
    sparseMtrxList [ name ] = creator;
    return true;
}

Dof *ClassFactory :: createDof(dofType name, int num, DofManager *dman)
{
    return ( dofList.count(name) == 1 ) ? dofList [ name ](num, dman) : NULL;
}

SparseLinearSystemNM *ClassFactory :: createSparseLinSolver(LinSystSolverType name, Domain *d, EngngModel *m)
{
    return ( sparseLinSolList.count(name) == 1 ) ? sparseLinSolList [ name ](d, m) : NULL;
}

bool ClassFactory :: registerSparseLinSolver( LinSystSolverType name, SparseLinearSystemNM * ( *creator )( Domain *, EngngModel * ) )
{
    sparseLinSolList [ name ] = creator;
    return true;
}

ErrorEstimator *ClassFactory :: createErrorEstimator(ErrorEstimatorType name, int num, Domain *d)
{
    return ( errEstList.count(name) == 1 ) ? errEstList [ name ](num, d) : NULL;
}

bool ClassFactory :: registerErrorEstimator( ErrorEstimatorType name, ErrorEstimator * ( *creator )( int, Domain * ) )
{
    errEstList [ name ] = creator;
    return true;
}

InitialCondition *ClassFactory :: createInitialCondition(const char *name, int num, Domain *d)
{
    CaseComp c;
    if ( c(name, "initialcondition") ) {
        return new InitialCondition(num, d);
    }
    return NULL;
}

NodalRecoveryModel *ClassFactory :: createNodalRecoveryModel(NodalRecoveryModel :: NodalRecoveryModelType type, Domain *d)
{
    if ( type == NodalRecoveryModel :: NRM_NodalAveraging ) {
        return new NodalAveragingRecoveryModel(d);
    } else if ( type == NodalRecoveryModel :: NRM_ZienkiewiczZhu ) {
        return new ZZNodalRecoveryModel(d);
    } else if ( type == NodalRecoveryModel :: NRM_SPR ) {
        return new SPRNodalRecoveryModel(d);
    }
    return NULL;
}


Element *ClassFactory :: createElement(const char *name, int number, Domain *domain)
{
    return ( elemList.count(name) == 1 ) ? elemList [ name ](number, domain) : NULL;
}

bool ClassFactory :: registerElement( const char *name, Element * ( *creator )( int, Domain * ) )
{
    elemList [ name ] = creator;
    return true;
}

DofManager *ClassFactory :: createDofManager(const char *name, int number, Domain *domain)
{
    return ( dofmanList.count(name) == 1 ) ? dofmanList [ name ](number, domain) : NULL;
}

bool ClassFactory :: registerDofManager( const char *name, DofManager * ( *creator )( int, Domain * ) )
{
    dofmanList [ name ] = creator;
    return true;
}

GeneralBoundaryCondition *ClassFactory :: createBoundaryCondition(const char *name, int number, Domain *domain)
{
    return ( bcList.count(name) == 1 ) ? bcList [ name ](number, domain) : NULL;
}

bool ClassFactory :: registerBoundaryCondition( const char *name, GeneralBoundaryCondition * ( *creator )( int, Domain * ) )
{
    bcList [ name ] = creator;
    return true;
}

CrossSection *ClassFactory :: createCrossSection(const char *name, int number, Domain *domain)
{
    return ( csList.count(name) == 1 ) ? csList [ name ](number, domain) : NULL;
}

bool ClassFactory :: registerCrossSection( const char *name, CrossSection * ( *creator )( int, Domain * ) )
{
    csList [ name ] = creator;
    return true;
}

Material *ClassFactory :: createMaterial(const char *name, int number, Domain *domain)
{
    return ( matList.count(name) == 1 ) ? matList [ name ](number, domain) : NULL;
}

bool ClassFactory :: registerMaterial( const char *name, Material * ( *creator )( int, Domain * ) )
{
    matList [ name ] = creator;
    return true;
}

EngngModel *ClassFactory :: createEngngModel(const char *name, int number, EngngModel *master)
{
    return ( engngList.count(name) == 1 ) ? engngList [ name ](number, master) : NULL;
}

bool ClassFactory :: registerEngngModel( const char *name, EngngModel * ( *creator )( int, EngngModel * ) )
{
    engngList [ name ] = creator;
    return true;
}

Function *ClassFactory :: createFunction(const char *name, int number, Domain *domain)
{
    return ( funcList.count(name) == 1 ) ? funcList [ name ](number, domain) : NULL;
}

bool ClassFactory :: registerFunction( const char *name, Function * ( * creator )( int, Domain * ) )
{
    funcList [ name ] = creator;
    return true;
}

NonlocalBarrier *ClassFactory :: createNonlocalBarrier(const char *name, int number, Domain *domain)
{
    return ( nlbList.count(name) == 1 ) ? nlbList [ name ](number, domain) : NULL;
}

bool ClassFactory :: registerNonlocalBarrier( const char *name, NonlocalBarrier * ( *creator )( int, Domain * ) )
{
    nlbList [ name ] = creator;
    return true;
}

RandomFieldGenerator *ClassFactory :: createRandomFieldGenerator(const char *name, int number, Domain *domain)
{
    return ( rfgList.count(name) == 1 ) ? rfgList [ name ](number, domain) : NULL;
}

bool ClassFactory :: registerRandomFieldGenerator( const char *name, RandomFieldGenerator * ( *creator )( int, Domain * ) )
{
    rfgList [ name ] = creator;
    return true;
}

ExportModule *ClassFactory :: createExportModule(const char *name, int number, EngngModel *emodel)
{
    return ( exportList.count(name) == 1 ) ? exportList [ name ](number, emodel) : NULL;
}

bool ClassFactory :: registerExportModule( const char *name, ExportModule * ( *creator )( int, EngngModel * ) )
{
    exportList [ name ] = creator;
    return true;
}

SparseNonLinearSystemNM *ClassFactory :: createNonLinearSolver(const char *name, Domain *d, EngngModel *emodel)
{
    return ( nonlinList.count(name) == 1 ) ? nonlinList [ name ](d, emodel) : NULL;
}

bool ClassFactory :: registerSparseNonLinearSystemNM( const char *name, SparseNonLinearSystemNM * ( *creator )( Domain *, EngngModel * ) )
{
    nonlinList [ name ] = creator;
    return true;
}

InitModule *ClassFactory :: createInitModule(const char *name, int number, EngngModel *emodel)
{
    return ( initList.count(name) == 1 ) ? initList [ name ](number, emodel) : NULL;
}

bool ClassFactory :: registerInitModule( const char *name, InitModule * ( *creator )( int, EngngModel * ) )
{
    initList [ name ] = creator;
    return true;
}

TopologyDescription *ClassFactory :: createTopology(const char *name, Domain *domain)
{
    return ( topologyList.count(name) == 1 ) ? topologyList [ name ](domain) : NULL;
}

bool ClassFactory :: registerTopologyDescription( const char *name, TopologyDescription * ( *creator )( Domain * ) )
{
    topologyList [ name ] = creator;
    return true;
}


// XFEM:
EnrichmentItem *ClassFactory :: createEnrichmentItem(const char *name, int number, XfemManager *xm, Domain *domain)
{
    return ( enrichItemList.count(name) == 1 ) ? enrichItemList [ name ](number, xm, domain) : NULL;
}

bool ClassFactory :: registerEnrichmentItem( const char *name, EnrichmentItem * ( *creator )( int, XfemManager *, Domain * ) )
{
    enrichItemList [ name ] = creator;
    return true;
}

EnrichmentFunction *ClassFactory :: createEnrichmentFunction(const char *name, int number, Domain *domain)
{
    return ( enrichFuncList.count(name) == 1 ) ? enrichFuncList [ name ](number, domain) : NULL;
}

bool ClassFactory :: registerEnrichmentFunction( const char *name, EnrichmentFunction * ( *creator )( int, Domain * ) )
{
    enrichFuncList [ name ] = creator;
    return true;
}

EnrichmentDomain *ClassFactory :: createEnrichmentDomain(const char *name)
{
    return ( enrichmentDomainList.count(name) == 1 ) ? enrichmentDomainList [ name ]() : NULL;
}

bool ClassFactory :: registerEnrichmentDomain( const char *name, EnrichmentDomain * ( *creator )( ) )
{
    enrichmentDomainList [ name ] = creator;
    return true;
}

EnrichmentFront *ClassFactory :: createEnrichmentFront(const char *name)
{
    return ( enrichmentFrontList.count(name) == 1 ) ? enrichmentFrontList [ name ]() : NULL;
}

bool ClassFactory :: registerEnrichmentFront( const char *name, EnrichmentFront * ( *creator )( ) )
{
    enrichmentFrontList [ name ] = creator;
    return true;
}

PropagationLaw *ClassFactory :: createPropagationLaw(const char *name)
{
    return ( propagationLawList.count(name) == 1 ) ? propagationLawList [ name ]() : NULL;
}

bool ClassFactory :: registerPropagationLaw( const char *name, PropagationLaw * ( *creator )( ) )
{
    propagationLawList [ name ] = creator;
    return true;
}

BasicGeometry *ClassFactory :: createGeometry(const char *name)
{
    return ( geometryList.count(name) == 1 ) ? geometryList [ name ]() : NULL;
}

bool ClassFactory :: registerGeometry( const char *name, BasicGeometry * ( *creator )( ) )
{
    geometryList [ name ] = creator;
    return true;
}


// Failure module:

FailureCriteria *ClassFactory :: createFailureCriteria(const char *name, int number, FractureManager *fracManager)
{
    return ( failureCriteriaList.count(name) == 1 ) ? failureCriteriaList [ name ](number, fracManager) : NULL;
}

bool ClassFactory :: registerFailureCriteria( const char *name, FailureCriteria * ( *creator )( int, FractureManager * ) )
{
    failureCriteriaList [ name ] = creator;
    return true;
}

FailureCriteriaStatus *ClassFactory :: createFailureCriteriaStatus(const char *name, int number, FailureCriteria *fc)
{
    return ( failureCriteriaList.count(name) == 1 ) ? failureCriteriaStatusList [ name ](number, fc) : NULL;
}

bool ClassFactory :: registerFailureCriteriaStatus( const char *name, FailureCriteriaStatus * ( *creator )( int, FailureCriteria * ) )
{
    failureCriteriaStatusList [ name ] = creator;
    return true;
}




SparseGeneralEigenValueSystemNM *ClassFactory :: createGeneralizedEigenValueSolver(GenEigvalSolverType st, Domain *d, EngngModel *m)
{
    if ( st == GES_SubspaceIt ) {
        return new SubspaceIteration(d, m);
    } else if ( st == GES_InverseIt ) {
        return new InverseIteration(d, m);
    }
#ifdef __SLEPC_MODULE
    else if ( st == GES_SLEPc ) {
        return new SLEPcSolver(d, m);
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
    if ( type == MMA_ClosestPoint ) {
        return new MMAClosestIPTransfer();
    } else if ( type == MMA_LeastSquareProjection ) {
        return new MMALeastSquareProjection();
    } else if ( type == MMA_ShapeFunctionProjection ) {
        return new MMAShapeFunctProjection();
    }
    return NULL;
}

MesherInterface *ClassFactory :: createMesherInterface(MeshPackageType type, Domain *d)
{
    if ( type == MPT_T3D ) {
        return new T3DInterface(d);
    } else if ( type == MPT_TARGE2 ) {
        return new Targe2Interface(d);
    } else if ( type == MPT_FREEM ) {
        return new FreemInterface(d);
    } else if ( type == MPT_SUBDIVISION ) {
        return new Subdivision(d);
    }
    return NULL;
}


#ifdef __PARALLEL_MODE
LoadBalancerMonitor *ClassFactory :: createLoadBalancerMonitor(const char *name, EngngModel *e)
{
    return ( loadMonitorList.count(name) == 1 ) ? loadMonitorList [ name ](e) : NULL;
}

LoadBalancer *ClassFactory :: createLoadBalancer(const char *name, Domain *d)
{
    return ( loadBalancerList.count(name) == 1 ) ? loadBalancerList [ name ](d) : NULL;
}
#endif
} // End namespace oofem
