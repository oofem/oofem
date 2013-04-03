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

#include "classfactory.h"

#include <cstring>
#ifdef HAVE_STRINGS_H
 #include <strings.h>
#endif
#include "compiler.h"

#define REGISTER_CLASS(_class, name, id)
#include "elementclassfactory.h"
#include "dofmanclassfactory.h"
#include "boundaryconditionclassfactory.h"
#include "crosssectionclassfactory.h"
#include "materialclassfactory.h"
#include "engngmodelclassfactory.h"
#include "ltfclassfactory.h"
#include "nonlocalbarrierclassfactory.h"
#include "randomfieldgeneratorclassfactory.h"

// No class factory files for these;
#include "subspaceit.h"
#include "inverseit.h"
#include "slepcsolver.h"

#include "vtkexportmodule.h"
#include "vtkxmlexportmodule.h"
#include "matlabexportmodule.h"
#include "gpexportmodule.h"
#ifdef __SM_MODULE
 #include "poiexportmodule.h"
 #include "homexportmodule.h"
 #include "dmexportmodule.h"
 
 // mesher interfaces
 #include "t3dinterface.h"
 #include "targe2interface.h"
 #include "freeminterface.h"
 #include "subdivision.h"

 #include "mmaclosestiptransfer.h"
 #include "mmaleastsquareprojection.h"
 #include "mmashapefunctprojection.h"
#endif

// or these
#include "masterdof.h"
#include "slavedof.h"
#include "simpleslavedof.h"
#include "activedof.h"
#include "nrsolver.h"
#include "calmls.h"
#ifdef __SM_MODULE
 #include "gpinitmodule.h"
#endif
#if 0 // Soon
#include "particletopologydescription.h"
#endif

#ifdef __PARALLEL_MODE
 #include "loadbalancer.h"
 #include "parmetisloadbalancer.h"
#endif

#undef  REGISTER_CLASS
#define REGISTER_CLASS(_class, id)
#include "sparsemtrxclassfactory.h"
#include "dofclassfactory.h"
#include "sparselinearsystemsolverclassfactory.h"
#include "errorestimatorclassfactory.h"
#include "initialcondition.h"

namespace oofem {

ClassFactory classFactory;


int ClassFactory::CaseComp::operator()(const std::string &a, const std::string &b) const
{
    return strncasecmp ( a.c_str(), b.c_str(), ( int ) b.length() ) < 0;
}


ClassFactory :: ClassFactory() {
    // register elements
#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, name, id) \
    elemNameList [ name ]  = elemCreator< _class >;
#include "elementclassfactory.h"

#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, name, id) \
    elemIdList [ id ] = elemCreator< _class >;
#include "elementclassfactory.h"

    // register dof managers
#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, name, id) \
    dofmanNameList [ name ]  = dofmanCreator< _class >;
#include "dofmanclassfactory.h"

#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, name, id) \
    dofmanIdList [ id ] = dofmanCreator< _class >;
#include "dofmanclassfactory.h"

    // register boundary conditions
#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, name, id) \
    bcNameList [ name ]  = bcCreator< _class >;
#include "boundaryconditionclassfactory.h"

#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, name, id) \
    bcIdList [ id ] = bcCreator< _class >;
#include "boundaryconditionclassfactory.h"

    // register cross sections
#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, name, id) \
    csNameList [ name ]  = csCreator< _class >;
#include "crosssectionclassfactory.h"

#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, name, id) \
    csIdList [ id ] = csCreator< _class >;
#include "crosssectionclassfactory.h"

    // register materials
#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, name, id) \
    matNameList [ name ]  = matCreator< _class >;
#include "materialclassfactory.h"

#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, name, id) \
    matIdList [ id ] = matCreator< _class >;
#include "materialclassfactory.h"

    // register engng models
#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, name, id) \
    engngNameList [ name ]  = engngCreator< _class >;
#include "engngmodelclassfactory.h"

#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, name, id) \
    engngIdList [ id ] = engngCreator< _class >;
#include "engngmodelclassfactory.h"

    // register load time functions
#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, name, id) \
    ltfNameList [ name ]  = ltfCreator< _class >;
#include "ltfclassfactory.h"

#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, name, id) \
    ltfIdList [ id ] = ltfCreator< _class >;
#include "ltfclassfactory.h"
    // register nonlocal barriers
#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, name, id) \
    nlbNameList [ name ]  = nlbCreator< _class >;
#include "nonlocalbarrierclassfactory.h"

#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, name, id) \
    nlbIdList [ id ] = nlbCreator< _class >;
#include "nonlocalbarrierclassfactory.h"
    // register random field generators
#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, name, id) \
    rfgNameList [ name ]  = rfgCreator< _class >;
#include "randomfieldgeneratorclassfactory.h"

#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, name, id) \
    rfgIdList [ id ] = rfgCreator< _class >;
#include "randomfieldgeneratorclassfactory.h"

#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, id) \
    sparseMtrxList [ id ] = sparseMtrxCreator< _class >;
#include "sparsemtrxclassfactory.h"

#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, id) \
    dofList [ id ] = dofCreator< _class >;
#include "dofclassfactory.h"

#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, id) \
    errEstList [ id ] = errEstCreator< _class >;
#include "errorestimatorclassfactory.h"

#undef REGISTER_CLASS
#define REGISTER_CLASS(_class, id) \
    sparseLinSolList [ id ] = sparseLinSolCreator< _class >;
#include "sparselinearsystemsolverclassfactory.h"

    // No separate files for export modules (yet)
    exportList [ "vtkxml" ]    = exportCreator< VTKXMLExportModule >;
    exportList [ "vtk" ]       = exportCreator< VTKExportModule >;
    exportList [ "matlab" ]    = exportCreator< MatlabExportModule >;
    exportList [ "gp" ]        = exportCreator< GPExportModule >;
#ifdef __SM_MODULE
    exportList [ "poi" ]       = exportCreator< POIExportModule >;
    exportList [ "dm" ]        = exportCreator< DofManExportModule >;
    exportList [ "hom" ]       = exportCreator< HOMExportModule >;
#endif //__SM_MODULE

    nonlinList [ "nrsolver" ]   = nonlinCreator< NRSolver >;
    nonlinList [ "calm" ]       = nonlinCreator< CylindricalALM >;

#ifdef __SM_MODULE
    initList [ "gpinitmodule" ] = initCreator< GPInitModule >;
#endif //__SM_MODULE

#if 0
    topologyList["particletopology"] = topologyCreator< ParticleTopologyDescription >;
#endif

    // No separate files for xfem (yet)
    enrichItemList [ "cracktip" ]      = enrichItemCreator< CrackTip >;
    enrichItemList [ "crackinterior" ] = enrichItemCreator< CrackInterior >;
    enrichItemList [ "inclusion" ]     = enrichItemCreator< Inclusion >;
    enrichItemList [ "delamination" ]  = enrichItemCreator< Delamination >;
    enrichItemList [ "multipledelamination" ]  = enrichItemCreator< Delamination >;

    enrichFuncList [ "discontinuousfunction" ] = enrichFuncCreator< DiscontinuousFunction >;
    enrichFuncList [ "branchfunction" ] = enrichFuncCreator< BranchFunction >;
    enrichFuncList [ "rampfunction" ] = enrichFuncCreator< RampFunction >;

    enrichmentDomainList [ "dofmanlist" ]  = enrichmentDomainCreator< DofManList >;
    enrichmentDomainList [ "wholedomain" ] = enrichmentDomainCreator< WholeDomain >;
    enrichmentDomainList [ "circle" ]      = enrichmentDomainCreator< EDBGCircle >;
    //enrichmentDomainList [ "line" ]        = enrichmentDomainCreator< BasicGeometryDomain<Line> >;

    geometryList [ "line" ] = geometryCreator< Line >;
    geometryList [ "circle" ] = geometryCreator< Circle >;
    geometryList [ "pointswarm" ] = geometryCreator< PointSwarm >; // just temporary
}

SparseMtrx *ClassFactory :: createSparseMtrx(SparseMtrxType type)
{
    return ( sparseMtrxList.count ( type ) == 1 ) ? sparseMtrxList [ type ] () : NULL;
}

Dof *ClassFactory :: createDof(dofType type, int num, DofManager* dman)
{
    return ( dofList.count ( type ) == 1 ) ? dofList [ type ] (num, dman) : NULL;
}

SparseLinearSystemNM *ClassFactory :: createSparseLinSolver(LinSystSolverType type, Domain *d, EngngModel *m)
{
    return ( sparseLinSolList.count ( type ) == 1 ) ? sparseLinSolList [ type ] (d, m) : NULL;
}

ErrorEstimator * ClassFactory :: createErrorEstimator(ErrorEstimatorType type, int num, Domain *d)
{
    return ( errEstList.count ( type ) == 1 ) ? errEstList [ type ] (num, d) : NULL;
}

InitialCondition * ClassFactory::createInitialCondition(classType type, int num, Domain *d)
{
    if (type == InitialConditionClass) {
        return new InitialCondition(num, d);
    }
    return NULL;
}

Patch * ClassFactory::createPatch(Patch :: PatchType type, Element *e)
{
    if ( type == Patch :: PT_TrianglePatch ) {
        return new TrianglePatch(e);
    }
    return NULL;
}

NodalRecoveryModel * ClassFactory::createNodalRecoveryModel(NodalRecoveryModel :: NodalRecoveryModelType type, Domain *d)
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


Element* ClassFactory::createElement ( const char* name, int number, Domain* domain )
{
    return ( elemNameList.count ( name ) == 1 ) ? elemNameList [ name ] ( number, domain ) : NULL;
}

Element* ClassFactory::createElement ( classType type, int number, Domain* domain )
{
    return ( elemIdList.count ( type ) == 1 ) ? elemIdList [ type ] ( number, domain ) : NULL;
}

DofManager* ClassFactory::createDofManager ( const char* name, int number, Domain* domain )
{
    return ( dofmanNameList.count ( name ) == 1 ) ? dofmanNameList [ name ] ( number, domain ) : NULL;
}

DofManager* ClassFactory::createDofManager ( classType type, int number, Domain* domain )
{
    return ( dofmanIdList.count ( type ) == 1 ) ? dofmanIdList [ type ] ( number, domain ) : NULL;
}

GeneralBoundaryCondition* ClassFactory::createBoundaryCondition ( const char* name, int number, Domain* domain )
{
    return ( bcNameList.count ( name ) == 1 ) ? bcNameList [ name ] ( number, domain ) : NULL;
}

GeneralBoundaryCondition* ClassFactory::createBoundaryCondition ( classType type, int number, Domain* domain )
{
    return ( bcIdList.count ( type ) == 1 ) ? bcIdList [ type ] ( number, domain ) : NULL;
}

CrossSection* ClassFactory::createCrossSection ( const char* name, int number, Domain* domain )
{
    return ( csNameList.count ( name ) == 1 ) ? csNameList [ name ] ( number, domain ) : NULL;
}

CrossSection* ClassFactory::createCrossSection ( classType type, int number, Domain* domain )
{
    return ( csIdList.count ( type ) == 1 ) ? csIdList [ type ] ( number, domain ) : NULL;
}

Material* ClassFactory::createMaterial ( const char* name, int number, Domain* domain )
{
    return ( matNameList.count ( name ) == 1 ) ? matNameList [ name ] ( number, domain ) : NULL;
}

Material* ClassFactory::createMaterial ( classType type, int number, Domain* domain )
{
    return ( matIdList.count ( type ) == 1 ) ? matIdList [ type ] ( number, domain ) : NULL;
}

EngngModel* ClassFactory::createEngngModel ( const char* name, int number, EngngModel* master )
{
    return ( engngNameList.count ( name ) == 1 ) ? engngNameList [ name ] ( number, master ) : NULL;
}

EngngModel* ClassFactory::createEngngModel ( classType type, int number, EngngModel* master )
{
    return ( engngIdList.count ( type ) == 1 ) ? engngIdList [ type ] ( number, master ) : NULL;
}

LoadTimeFunction* ClassFactory::createLoadTimeFunction ( const char* name, int number, Domain* domain )
{
    return ( ltfNameList.count ( name ) == 1 ) ? ltfNameList [ name ] ( number, domain ) : NULL;
}

LoadTimeFunction* ClassFactory::createLoadTimeFunction ( classType type, int number, Domain* domain )
{
    return ( ltfIdList.count ( type ) == 1 ) ? ltfIdList [ type ] ( number, domain ) : NULL;
}

NonlocalBarrier* ClassFactory::createNonlocalBarrier ( const char* name, int number, Domain* domain )
{
    return ( nlbNameList.count ( name ) == 1 ) ? nlbNameList [ name ] ( number, domain ) : NULL;
}

NonlocalBarrier* ClassFactory::createNonlocalBarrier ( classType type, int number, Domain* domain )
{
    return ( nlbIdList.count ( type ) == 1 ) ? nlbIdList [ type ] ( number, domain ) : NULL;
}

RandomFieldGenerator* ClassFactory::createRandomFieldGenerator ( const char* name, int number, Domain* domain )
{
    return ( rfgNameList.count ( name ) == 1 ) ? rfgNameList [ name ] ( number, domain ) : NULL;
}

RandomFieldGenerator* ClassFactory::createRandomFieldGenerator ( classType type, int number, Domain* domain )
{
    return ( rfgIdList.count ( type ) == 1 ) ? rfgIdList [ type ] ( number, domain ) : NULL;
}

ExportModule* ClassFactory::createExportModule(const char *name, int number, EngngModel *emodel)
{
    return ( exportList.count(name) == 1 ) ? exportList [ name ](number, emodel) : NULL;
}

SparseNonLinearSystemNM* ClassFactory::createNonLinearSolver(const char *name, Domain *d, EngngModel *emodel, EquationID eid)
{
    return ( nonlinList.count(name) == 1 ) ? nonlinList [ name ](d, emodel, eid) : NULL;
}

InitModule* ClassFactory::createInitModule(const char *name, int number, EngngModel *emodel)
{
    return ( initList.count(name) == 1 ) ? initList [ name ](number, emodel) : NULL;
}

TopologyDescription* ClassFactory::createTopology(const char *name, Domain *domain)
{
    return ( topologyList.count(name) == 1 ) ? topologyList [ name ](domain) : NULL;
}


// XFEM:
EnrichmentItem* ClassFactory::createEnrichmentItem(const char *name, int number, XfemManager *xm, Domain *domain)
{
    return ( enrichItemList.count(name) == 1 ) ? enrichItemList [ name ](number, xm, domain) : NULL;
}

EnrichmentFunction* ClassFactory::createEnrichmentFunction(const char *name, int number, Domain *domain)
{
    return ( enrichFuncList.count(name) == 1 ) ? enrichFuncList [ name ](number, domain) : NULL;
}

EnrichmentDomain* ClassFactory::createEnrichmentDomain(const char *name)
{
    return ( enrichmentDomainList.count(name) == 1 ) ? enrichmentDomainList [ name ]() : NULL;
}

BasicGeometry* ClassFactory::createGeometry(const char *name)
{
    return ( geometryList.count(name) == 1 ) ? geometryList [ name ]() : NULL;
}


SparseGeneralEigenValueSystemNM* ClassFactory::createGeneralizedEigenValueSolver(GenEigvalSolverType st, Domain *d, EngngModel *m)
{
    if ( st == GES_SubspaceIt ) {
        return new SubspaceIteration(d, m);
    } else if ( st == GES_InverseIt ) {
        return new InverseIteration(d, m);
    } else if ( st == GES_SLEPc ) {
        return new SLEPcSolver(d, m);
    }
    return NULL;
}

IntegrationRule* ClassFactory::createIRule(classType type, int number, Element *e)
{
    if ( type == GaussIntegrationRuleClass ) {
        return new GaussIntegrationRule(number, e);
    }
    return NULL;
}

MaterialMappingAlgorithm* ClassFactory::createMaterialMappingAlgorithm(MaterialMappingAlgorithmType type)
{
#ifdef __SM_MODULE
    if ( type == MMA_ClosestPoint ) {
        return new MMAClosestIPTransfer();
    } else if ( type == MMA_LeastSquareProjection ) {
        return new MMALeastSquareProjection();
    } else if ( type == MMA_ShapeFunctionProjection ) {
        return new MMAShapeFunctProjection();
    }
#endif
    return NULL;
}

MesherInterface* ClassFactory::createMesherInterface(MeshPackageType type, Domain *d)
{
#ifdef __SM_MODULE
    if ( type == MPT_T3D ) {
        return new T3DInterface(d);
    } else if ( type == MPT_TARGE2 ) {
        return new Targe2Interface(d);
    } else if ( type == MPT_FREEM ) {
        return new FreemInterface(d);
    } else if ( type == MPT_SUBDIVISION ) {
        return new Subdivision(d);
    }
#endif
    return NULL;
}


#ifdef __PARALLEL_MODE
LoadBalancerMonitor* ClassFactory::createLoadBalancerMonitor(classType type, EngngModel *e)
{
    if ( type == WallClockLoadBalancerMonitorClass ) {
        return new WallClockLoadBalancerMonitor(e);
    }
    return NULL;
}

LoadBalancer* ClassFactory::createLoadBalancer(classType type, Domain *d)
{
    if ( type == ParmetisLoadBalancerClass ) {
        return new ParmetisLoadBalancer(d);
    }
    return NULL;
}
#endif
} // End namespace oofem
