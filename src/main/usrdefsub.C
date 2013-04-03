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

#include "usrdefsub.h"
#include "classfactory.h"

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
Dof *CreateUsrDefDofOfType(dofType type, int number, DofManager *dman)
{
    return classFactory.createDof(type, number, dman);
}

// ================ SparseLinearSystemNM CLASS FACTORY==================
SparseLinearSystemNM *CreateUsrDefSparseLinSolver(LinSystSolverType st, Domain *d, EngngModel *m)
{
    return classFactory.createSparseLinSolver(st,d,m);
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

SparseNonLinearSystemNM *CreateUsrDefNonLinearSolver(const char *aClass, Domain *d, EngngModel *emodel)
{
    return classFactory.createNonLinearSolver(aClass,d,emodel);
}

InitModule *CreateUsrDefInitModuleOfType(const char *aClass, int number, EngngModel *emodel)
{
    return classFactory.createInitModule(aClass, number, emodel);
}

ExportModule *CreateUsrDefExportModuleOfType(const char *aClass, int number, EngngModel *emodel)
{
    return classFactory.createExportModule(aClass,number,emodel);
}

TopologyDescription *CreateUsrDefTopologyOfType(const char *aClass, Domain *domain)
{
    return classFactory.createTopology(aClass, domain);
}

Patch *CreateUsrDefPatch(Patch :: PatchType ptype, Element *e)
{
    return classFactory.createPatch(ptype, e);
}

NodalRecoveryModel *CreateUsrDefNodalRecoveryModel(NodalRecoveryModel :: NodalRecoveryModelType type, Domain *d)
{
    return classFactory.createNodalRecoveryModel(type, d);
}

EnrichmentItem *CreateUsrDefEnrichmentItem(const char *aClass, int number, XfemManager *xm, Domain *domain)
{
    return classFactory.createEnrichmentItem(aClass, number, xm, domain);
}

EnrichmentFunction *CreateUsrDefEnrichmentFunction(const char *aClass, int number, Domain *domain)
{
    return classFactory.createEnrichmentFunction(aClass, number, domain);
}

EnrichmentDomain *CreateUsrDefEnrichmentDomain(const char *aClass)
{
    return classFactory.createEnrichmentDomain(aClass);
}

BasicGeometry *CreateUsrDefGeometry(const char *aClass)
{
    return classFactory.createGeometry(aClass);
}

SparseGeneralEigenValueSystemNM *CreateUsrDefGeneralizedEigenValueSolver(GenEigvalSolverType st, Domain *d, EngngModel *m)
{
    return classFactory.createGeneralizedEigenValueSolver(st, d, m);
}

IntegrationRule *CreateUsrDefIRuleOfType(classType type, int number, Element *e)
{
    return classFactory.createIRule(type, number, e);
}

MaterialMappingAlgorithm *CreateUsrDefMaterialMappingAlgorithm(MaterialMappingAlgorithmType type)
{
    return classFactory.createMaterialMappingAlgorithm(type);
}

MesherInterface *CreateUsrDefMesherInterface(MeshPackageType type, Domain *d)
{
    return classFactory.createMesherInterface(type, d);
}

#ifdef __PARALLEL_MODE
LoadBalancerMonitor *CreateUsrDefLoadBalancerMonitorOfType(classType type, EngngModel *e)
{
    return classFactory.createLoadBalancerMonitor(type, e);
}

LoadBalancer *CreateUsrDefLoadBalancerOfType(classType type, Domain *d)
{
    return classFactory.createLoadBalancer(type, d);
}
#endif
} // end namespace oofem
