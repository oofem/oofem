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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "usrdefsub.h"
#include "classfactory.h"
#include "compiler.h" // to supply missing strncasecmp on some platforms
#ifndef __MAKEDEPEND
 #include <string.h>
 #ifdef HAVE_STRINGS_H
  #include <strings.h>
 #endif
#endif


// __OOFEMLIB_MODULE
#include "node.h"
#include "element.h"
#include "engngm.h"
#include "xfemmanager.h"
#include "load.h"
#include "loadtime.h"
#include "material.h"
#include "gaussintegrationrule.h"

#include "subspaceit.h"
#include "inverseit.h"
#include "slepcsolver.h"

// Nonlinear solvers
#include "nrsolver.h"
#include "nrsolver2.h"
#include "calmls.h"

// export modules
#include "vtkexportmodule.h"
#include "vtkxmlexportmodule.h"

// nodal recovery models
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"

#if 0 // Soon
 #include "particletopologydescription.h" // Soon
#endif
// end __OOFEMLIB_MODULE


#ifdef __SM_MODULE

// export modules
 #include "poiexportmodule.h"
 #include "homexportmodule.h"
 #include "dmexportmodule.h"
 #include "gpexportmodule.h"

// init modules
 #include "gpinitmodule.h"

// mesher interfaces
 #include "t3dinterface.h"
 #include "targe2interface.h"
 #include "freeminterface.h"
 #include "subdivision.h"

 #include "dss.h"

 #include "mmaclosestiptransfer.h"
 #include "mmaleastsquareprojection.h"
#endif //__SM_MODULE



#ifdef __TM_MODULE
#endif //__TM_MODULE


#ifdef __FM_MODULE
#endif // __FM_Module

// GENERAL
#include "masterdof.h"
#include "slavedof.h"
#include "simpleslavedof.h"
#include "activedof.h"

// XFEM
#include "enrichmentfunction.h"
#include "enrichmentitem.h"

#ifdef __PARALLEL_MODE
 #include "loadbalancer.h"
 #include "parmetisloadbalancer.h"
#endif

#include <string>

/* ======================== Note =============================================
 * The system of new class registration has been changed to more flexible one.
 * New classes are registered in _baseclassname_classfactory.h files
 * where _baseclassname_ is the name of corresponding base class name.
 * Following factories are supported:
 *   Elements            - elementclassfactory.h
 *   Dof managers        - dofmanclassfactory.h
 *   DOFs                - dofclassfactory.h
 *   Boundary conditions - boundaryconditionclassfactory.h
 *   Materials           - materialclassfactory.h
 *   Cross sections      - crosssectionclassfactory.h
 *   Engng models        - engngmodelclassfactory.h
 *   Load time functions - ltfclassfactory.h
 *   Sparse matrices     - sparsemtrxclassfactory.h
 *   Sparse linear system solvers - sparselinearsystemsolverclassfactory.h
 *   Error estimators    - errorestimatorclassfactory.h
 *
 * Refer to these files when registering new classes.
 * ===========================================================================
 */

namespace oofem {
// Comparison operator for strings. Just strcasecmp here?
struct CaseComp
{
    int operator() (std :: string a, std :: string b) const { return strncasecmp( a.c_str(), b.c_str(), b.length() ) < 0; }
};
struct CaseCompId
{
    int operator() (classType a, classType b) const { return ( ( ( int ) a ) < ( ( int ) b ) ); }
};


// ================ELEMENT CLASS FACTORY==================

Element *CreateUsrDefElementOfType(const char *aClass, int number, Domain *domain)
{
    return classFactory.createElement(aClass, number, domain);
}

Element *CreateUsrDefElementOfType(classType type, int number, Domain *domain)
{
    return classFactory.createElement(type, number, domain);
}

// ================DOFMAN CLASS FACTORY==================

DofManager *CreateUsrDefDofManagerOfType(const char *aClass, int number, Domain *domain)
{
    return classFactory.createDofManager(aClass, number, domain);
}

DofManager *CreateUsrDefDofManagerOfType(classType type, int number, Domain *domain)
{
    return classFactory.createDofManager(type, number, domain);
}

// ================ BC CLASS FACTORY==================


GeneralBoundaryCondition *CreateUsrDefBoundaryConditionOfType(const char *aClass, int number, Domain *domain)
{
    return classFactory.createBoundaryCondition(aClass, number, domain);
}

GeneralBoundaryCondition *CreateUsrDefBoundaryConditionOfType(classType type, int number, Domain *domain)
{
    return classFactory.createBoundaryCondition(type, number, domain);
}
// ================ IC CLASS FACTORY==================
InitialCondition *CreateUsrDefInitialConditionOfType(classType type, int number, Domain *domain)
{
    return classFactory.createInitialCondition(type, number, domain);
}

// ================ CROSSSECTION CLASS FACTORY==================

CrossSection *CreateUsrDefCrossSectionOfType(const char *aClass, int number, Domain *domain)
{
    return classFactory.createCrossSection(aClass, number, domain);
}

CrossSection *CreateUsrDefCrossSectionOfType(classType type, int number, Domain *domain)
{
    return classFactory.createCrossSection(type, number, domain);
}

// ================ MATERIAL CLASS FACTORY==================

Material *CreateUsrDefMaterialOfType(const char *aClass, int number, Domain *domain)
{
    return classFactory.createMaterial(aClass, number, domain);
}

Material *CreateUsrDefMaterialOfType(classType type, int number, Domain *domain)
{
    return classFactory.createMaterial(type, number, domain);
}

// ================ ENGNG MODEL CLASS FACTORY==================
EngngModel *CreateUsrDefEngngModelOfType(const char *aClass, int number, EngngModel *master)
{
    return classFactory.createEngngModel(aClass, number, master);
}

EngngModel *CreateUsrDefEngngModelOfType(classType type, int number, EngngModel *master)
{
    return classFactory.createEngngModel(type, number, master);
}

// ================ LOADTIMEFUNCTION CLASS FACTORY==================
LoadTimeFunction *CreateUsrDefLoadTimeFunctionOfType(const char *aClass, int number, Domain *domain)
{
    return classFactory.createLoadTimeFunction(aClass, number, domain);
}

LoadTimeFunction *CreateUsrDefLoadTimeFunctionOfType(classType type, int number, Domain *domain)
{
    return classFactory.createLoadTimeFunction(type, number, domain);
}

// ================ SparseMtrx CLASS FACTORY==================
SparseMtrx *CreateUsrDefSparseMtrx(SparseMtrxType type)
{
    return classFactory.createSparseMtrx(type);
}

// ================ DOF CLASS FACTORY==================
Dof *CreateUsrDefDofOfType(classType type, int number, DofManager *dman)
{
  return classFactory.createDof (type, number, dman);
}

// ================ SparseLinearSystemNM CLASS FACTORY==================
SparseLinearSystemNM *CreateUsrDefSparseLinSolver(LinSystSolverType st, int i, Domain *d, EngngModel *m)
{
  return classFactory.createSparseLinSolver(st,i,d,m);
}

// ================ ErrorEstimator CLASS FACTORY==================
ErrorEstimator *CreateUsrDefErrorEstimator(ErrorEstimatorType type, int number, Domain *d)
{
  return classFactory.createErrorEstimator(type,number,d);
}

// ================ Nonlocal Barrier CLASS FACTORY==================
NonlocalBarrier *CreateUsrDefNonlocalBarrierOfType(const char *aClass, int number, Domain *domain)
{
  return classFactory.createNonlocalBarrier(aClass,number,domain);
}

NonlocalBarrier *CreateUsrDefNonlocalBarrierOfType(classType type, int number, Domain *domain)
{
  return classFactory.createNonlocalBarrier(type,number,domain);
}

// ================ Random Field Generator CLASS FACTORY==================
RandomFieldGenerator *CreateUsrDefRandomFieldGenerator(const char *aClass, int number, Domain *domain)
{
  return classFactory.createRandomFieldGenerator(aClass,number,domain);
}

RandomFieldGenerator *CreateUsrDefRandomFieldGenerator(classType type, int number, Domain *domain)
{
  return classFactory.createRandomFieldGenerator(type,number,domain);
}


//------------------OLD Style CreateUsrDef Functions-------------------------------------------

// ================ Topology CLASS FACTORY==================
// Template to wrap constructors into functions
template< typename T > TopologyDescription *topologyCreator(Domain *d) { return new T(d); }
std :: map < std :: string, TopologyDescription * ( * )(Domain *), CaseComp > topologyNameList;

TopologyDescription *CreateUsrDefTopologyOfType(const char *aClass, Domain *domain)
{
#if 0
    if ( topologyNameList.size() == 0 ) { topologyNameList["particletopology"] = topologyCreator< ParticleTopologyDescription >; }
#endif

    return ( topologyNameList.count(aClass) == 1 ) ? topologyNameList [ aClass ](domain) : NULL;
}

SparseGeneralEigenValueSystemNM *CreateUsrDefGeneralizedEigenValueSolver(GenEigvalSolverType st, int i, Domain *d, EngngModel *m)
{
    SparseGeneralEigenValueSystemNM *nm = NULL;
    if ( st == GES_SubspaceIt ) {
        nm = ( SparseGeneralEigenValueSystemNM * ) new SubspaceIteration(i, d, m);
        return nm;
    } else if ( st == GES_InverseIt ) {
        nm = ( SparseGeneralEigenValueSystemNM * ) new InverseIteration(i, d, m);
        return nm;
    } else if ( st == GES_SLEPc ) {
        nm = ( SparseGeneralEigenValueSystemNM * ) new SLEPcSolver(i, d, m);
        return nm;
    } else {
        OOFEM_ERROR("CreateUsrDefGeneralizedEigenValueSolver: Unknown solver type\n");
    }

    return nm;
}


template< typename T > SparseNonLinearSystemNM *nonlinCreator(int n, Domain *d, EngngModel *m, EquationID eid) { return ( new T(n, d, m, eid) ); }
std :: map < std :: string, SparseNonLinearSystemNM * ( * )(int, Domain *, EngngModel *, EquationID), CaseComp > nonlinList;

SparseNonLinearSystemNM *CreateUsrDefNonLinearSolver(const char *aClass, int number, Domain *d, EngngModel *emodel, EquationID eid)
{
    if ( nonlinList.size() == 0 ) {
        //nonlinList["snes"]       = nonlinCreator< PETScSNES >;
        nonlinList [ "nrsolver" ]   = nonlinCreator< NRSolver >;
        nonlinList [ "nrsolver2" ]  = nonlinCreator< NRSolver2 >;
        nonlinList [ "calm" ]       = nonlinCreator< CylindricalALM >;
    }

    return ( nonlinList.count(aClass) == 1 ) ? nonlinList [ aClass ](number, d, emodel, eid) : NULL;
}

template< typename T > ExportModule *exportCreator(int n, EngngModel *e) { return ( new T(n, e) ); }
std :: map < std :: string, ExportModule * ( * )(int, EngngModel *), CaseComp > exportList;

ExportModule *CreateUsrDefExportModuleOfType(const char *aClass, int number, EngngModel *emodel)
{
    if ( exportList.size() == 0 ) {
        exportList [ "vtkxml" ]    = exportCreator< VTKXMLExportModule >;
        exportList [ "vtk" ]       = exportCreator< VTKExportModule >;
#ifdef __SM_MODULE
        exportList [ "poi" ]       = exportCreator< POIExportModule >;
        exportList [ "hom" ]       = exportCreator< HOMExportModule >;
        exportList [ "dm" ]        = exportCreator< DofManExportModule >;
        exportList [ "gp" ]        = exportCreator< GPExportModule >;
#endif //__SM_MODULE
    }

    return ( exportList.count(aClass) == 1 ) ? exportList [ aClass ](number, emodel) : NULL;
}


template< typename T > InitModule *initCreator(int n, EngngModel *e) { return ( new T(n, e) ); }
std :: map < std :: string, InitModule * ( * )(int, EngngModel *), CaseComp > initList;

InitModule *CreateUsrDefInitModuleOfType(const char *aClass, int number, EngngModel *emodel)
{
    if ( initList.size() == 0 ) {
#ifdef __SM_MODULE
        initList [ "gpinitmodule" ] = initCreator< GPInitModule >;
#endif //__SM_MODULE
    }

    return ( initList.count(aClass) == 1 ) ? initList [ aClass ](number, emodel) : NULL;
}


IntegrationRule *CreateUsrDefIRuleOfType(classType type, int number, Element *e)
{
    IntegrationRule *answer = NULL;
    if ( type == GaussIntegrationRuleClass ) {
        answer = new GaussIntegrationRule(number, e);
    }

    if ( answer == NULL ) {
        OOFEM_ERROR2("CreateUsrDefIRuleOfType: Unknown integration rule type [%d]", type);
    }

    return answer;
}




MaterialMappingAlgorithm *CreateUsrDefMaterialMappingAlgorithm(MaterialMappingAlgorithmType type)
{
    MaterialMappingAlgorithm *answer = NULL;

#ifdef __SM_MODULE
    if ( type == MMA_ClosestPoint ) {
        answer = new MMAClosestIPTransfer();
    } else if ( type == MMA_LeastSquareProjection ) {
        answer = new MMALeastSquareProjection();
    } else if ( type == MMA_ShapeFunctionProjection ) {
        answer = new MMAShapeFunctProjection();
    }

#endif
    if ( answer == NULL ) {
        OOFEM_ERROR2("CreateUsrDefMaterialMappingAlgorithm: Unknown mma type [%d]", type);
    }

    return answer;
}

MesherInterface *CreateUsrDefMesherInterface(MeshPackageType type, Domain *d)
{
    MesherInterface *answer = NULL;
#ifdef __SM_MODULE
    if ( type == MPT_T3D ) {
        answer = new T3DInterface(d);
    } else if ( type == MPT_TARGE2 ) {
        answer = new Targe2Interface(d);
    } else if ( type == MPT_FREEM ) {
        answer = new FreemInterface(d);
    } else if ( type == MPT_SUBDIVISION ) {
        answer = new Subdivision(d);
    }

#endif
    if ( answer == NULL ) {
        OOFEM_ERROR2("CreateUsrDefMesherInterface: Unknown MI type [%d]", type);
    }

    return answer;
}


template< typename T > EnrichmentItem *enrichItemCreator(int n, XfemManager *x, Domain *d) { return new T(n, x, d); }
std :: map < std :: string, EnrichmentItem * ( * )(int, XfemManager *, Domain *), CaseComp > enrichItemList;

EnrichmentItem *CreateUsrDefEnrichmentItem(const char *aClass, int number, XfemManager *xm, Domain *domain)
{
    if ( enrichItemList.size() == 0 ) {
        enrichItemList [ "cracktip" ]      = enrichItemCreator< CrackTip >;
        enrichItemList [ "crackinterior" ] = enrichItemCreator< CrackInterior >;
        enrichItemList [ "inclusion" ]     = enrichItemCreator< Inclusion >;
    }

    return ( enrichItemList.count(aClass) == 1 ) ? enrichItemList [ aClass ](number, xm, domain) : NULL;
}

EnrichmentItem *CreateUsrDefEnrichmentItem(classType type, int number, XfemManager *xm, Domain *domain)
{
  OOFEM_ERROR("CreateUsrDefEnrichmentItem: Unknown type\n");
  return NULL;
}


template< typename T > EnrichmentFunction *enrichFuncCreator(int n, Domain *d) { return new T(n, d); }
std :: map < std :: string, EnrichmentFunction * ( * )(int, Domain *), CaseComp > enrichFuncList;

EnrichmentFunction *CreateUsrDefEnrichmentFunction(const char *aClass, int number, Domain *domain)
{
    if ( enrichFuncList.size() == 0 ) {
        enrichFuncList [ "discontinuousfunction" ] = enrichFuncCreator< DiscontinuousFunction >;
        enrichFuncList [ "branchfunction" ] = enrichFuncCreator< BranchFunction >;
        enrichFuncList [ "rampfunction" ] = enrichFuncCreator< RampFunction >;
    }

    return ( enrichFuncList.count(aClass) == 1 ) ? enrichFuncList [ aClass ](number, domain) : NULL;
}

EnrichmentFunction *CreateUsrDefEnrichmentFunction(classType type, int number, Domain *domain)
{
  OOFEM_ERROR("CreateUsrDefEnrichmentFunction: Unknown type\n");
  return NULL;
}

template< typename T > BasicGeometry *geometryCreator() { return new T(); }
std :: map < std :: string, BasicGeometry * ( * )(), CaseComp > geometryList;

BasicGeometry *CreateUsrDefGeometry(const char *aClass)
{
    if ( geometryList.size() == 0 ) {
        geometryList [ "line" ] = geometryCreator< Line >;
        geometryList [ "circle" ] = geometryCreator< Circle >;
    }

    return ( geometryList.count(aClass) == 1 ) ? geometryList [ aClass ]() : NULL;
}
BasicGeometry *CreateUsrDefGeometry(classType type)
{
  OOFEM_ERROR("CreateUsrDefGeometry: Unknown type\n");
  return NULL;
}

Patch *CreateUsrDefPatch(Patch :: PatchType ptype, Element *e)
{
    Patch *answer = NULL;
    if ( ptype == Patch :: PT_TrianglePatch ) {
        answer = new TrianglePatch(e);
    } else {
        OOFEM_ERROR2("CreateUsrDefPatch: Unknown PatchType [%d]", ptype);
    }

    return answer;
}

NodalRecoveryModel *
CreateUsrDefNodalRecoveryModel(NodalRecoveryModel :: NodalRecoveryModelType type, Domain *d)
{
    NodalRecoveryModel *answer = NULL;
    if ( type == NodalRecoveryModel :: NRM_NodalAveraging ) {
        answer = new NodalAveragingRecoveryModel(d);
    } else if ( type == NodalRecoveryModel :: NRM_ZienkiewiczZhu ) {
        answer = new ZZNodalRecoveryModel(d);
    } else if ( type == NodalRecoveryModel :: NRM_SPR ) {
        answer = new SPRNodalRecoveryModel(d);
    } else {
        OOFEM_ERROR2("CreateUsrDefNodalRecoveryModel: unsupported NodalRecoveryModelType [%d]", type);
    }

    return answer;
}


#ifdef __PARALLEL_MODE
LoadBalancerMonitor *CreateUsrDefLoadBalancerMonitorOfType(classType type, EngngModel *e)
{
    LoadBalancerMonitor *answer = NULL;
    if ( type == WallClockLoadBalancerMonitorClass ) {
        answer = new WallClockLoadBalancerMonitor(e);
    }

    if ( answer == NULL ) {
        OOFEM_ERROR2("CreateUsrDefLoadBalancerMonitorOfType: Unknown type [%d]", type);
    }

    return answer;
}

LoadBalancer *CreateUsrDefLoadBalancerOfType(classType type, Domain *d)
{
    LoadBalancer *answer = NULL;
    if ( type == ParmetisLoadBalancerClass ) {
        answer = new ParmetisLoadBalancer(d);
    }

    if ( answer == NULL ) {
        OOFEM_ERROR2("CreateUsrDefLoadBalancerOfType: Unknown type [%d]", type);
    }

    return answer;
}
#endif
} // end namespace oofem
