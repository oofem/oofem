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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

/**
 * @file: usrdefsub.h
 * User defined Subroutines
 */

#ifndef usrdefsub_h
#define usrdefsub_h

#include "classtype.h"
#include "sparselinsystemnm.h" // for LinSystSolverType
#include "patch.h" // for PatchType
#include "nodalrecoverymodel.h" // for NodalRecoveryModelType
#include "materialmappingalgorithmtype.h"
#include "meshpackagetype.h"
#include "sparsemtrxtype.h"
#include "geneigvalsolvertype.h"
#include "errorestimatortype.h"

namespace oofem {

class Element;
class DofManager;
class CrossSection;
class EngngModel;
class GeneralBoundaryCondition;
class LoadTimeFunction;
class Material;
class TopologyDescription;
class SparseMtrx;
class SparseLinearSystemNM;
class SparseGeneralEigenValueSystemNM;
class SparseNonLinearSystemNM;
class ErrorEstimator;
class ExportModule;
class InitModule;
class NonlocalBarrier;
class RandomFieldGenerator;
class IntegrationRule;
class Dof;
class MaterialMappingAlgorithm;
class MesherInterface;
class EnrichmentItem;
class EnrichmentFunction;
class BasicGeometry;
class Patch;
class NodalRecoveryModel;
class LoadBalancerMonitor;
class LoadBalancer;

/**
 * Creates new instance of element corresponding to given element keyword.
 * @param name Keyword string determining the type of new instance.
 * @param num  Element number.
 * @param d    Domain assigned to new element.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
Element *CreateUsrDefElementOfType(const char *name, int num, Domain *d);
/**
 * Creates new instance of user defined dof manager corresponding to given keyword.
 * @param name Keyword string determining the type of new instance
 * @param num  Component number.
 * @param d    Domain assigned to new object.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
DofManager *CreateUsrDefDofManagerOfType(const char *name, int num, Domain *d);
/**
 * Creates new instance of user defined cross section model corresponding to given keyword.
 * @param name Keyword string determining the type of new instance
 * @param num  Component number
 * @param d    Domain assigned to new object
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
CrossSection *CreateUsrDefCrossSectionOfType(const char *name, int num, Domain *d);
/**
 * Creates new instance of user defined engineering model corresponding to given keyword.
 * @param name Keyword string determining the type of new instance.
 * @param num  Component number.
 * @param master Master engineering model.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
EngngModel *CreateUsrDefEngngModelOfType(const char *name, int num, EngngModel *master = NULL);
/**
 * Creates new instance of user defined load  corresponding to given keyword.
 * @param name Keyword string determining the type of new instance.
 * @param num  Component number.
 * @param d    Domain assigned to new object.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
GeneralBoundaryCondition *CreateUsrDefBoundaryConditionOfType(const char *name, int num, Domain *d);
/**
 * Creates new instance of user defined load time function corresponding to given keyword.
 * @param name Keyword string determining the type of new instance.
 * @param num  Component number.
 * @param d    Domain assigned to new object.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
LoadTimeFunction *CreateUsrDefLoadTimeFunctionOfType(const char *name, int num, Domain *d);
/**
 * Creates new instance of user defined material model corresponding to given keyword.
 * @param name Keyword string determining the type of new instance.
 * @param num  Component number.
 * @param d    Domain assigned to new object.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
Material *CreateUsrDefMaterialOfType(const char *name, int num, Domain *d);
/**
 * Creates new instance of user defined engineering model corresponding to given keyword.
 * @param name Keyword string determining the type of new instance.
 * @param d    Domain assigned to new object.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
TopologyDescription *CreateUsrDefTopologyOfType(const char *name, Domain *d);
/**
 * Creates new empty instance of sparse matrix of given type.
 * @param type Determines sparseMtrx type.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
SparseMtrx *CreateUsrDefSparseMtrx(SparseMtrxType type);
/**
 * Creates new Instance of Linear Sparse Solver of given type.
 * @param st Solver type.
 * @param i Component number.
 * @param d Domain.
 * @param m Engineering model.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
SparseLinearSystemNM *CreateUsrDefSparseLinSolver(LinSystSolverType st, int i, Domain *d, EngngModel *m);
/**
 * Creates new instance of Generalized Eigenvalue Solver of given type.
 * @param st Solver type.
 * @param i Component number.
 * @param d Domain.
 * @param m Engineering model.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
SparseGeneralEigenValueSystemNM *CreateUsrDefGeneralizedEigenValueSolver(GenEigvalSolverType st, int i, Domain *d, EngngModel *m);
/**
 * Creates new instance of a nonlinear solver of given name.
 * @param name Keyword string determining the type of new instance.
 * @param i Component number.
 * @param d Domain.
 * @param m Engineering model.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
SparseNonLinearSystemNM *CreateUsrDefNonlinearSolver(const char *name, int i, Domain *d, EngngModel *m);
/**
 * Creates new instance of Error Estimator of given type.
 * @param type Determines Error Estimator type.
 * @param number Component number.
 * @param d Domain associated to ee.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
ErrorEstimator *CreateUsrDefErrorEstimator(ErrorEstimatorType type, int number, Domain *d);

/**
 * Creates new instance of Export Module of given name.
 * @param name Export module keyword.
 * @param number Component number.
 * @param emodel Engineering model associated to new export module.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
ExportModule *CreateUsrDefExportModuleOfType(const char *name, int number, EngngModel *emodel);

/**
 * Creates new Instance of Initialization Module of given name.
 * @param name Initialization module keyword.
 * @param number Component number.
 * @param emodel Engineering model associated to new export module.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
InitModule *CreateUsrDefInitModuleOfType(const char *name, int number, EngngModel *emodel);

/**
 * Creates new Instance of Nonlocal Barrier class corresponding to given name.
 * @param name Nonlocal barrier keyword.
 * @param num  Component number.
 * @param d    Domain assigned to new object.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
NonlocalBarrier *CreateUsrDefNonlocalBarrierOfType(const char *name, int num, Domain *d);
/**
 * Creates new Instance of Random generator class corresponding to given name.
 * @param name Random generator keyword.
 * @param num  Component number.
 * @param d    Domain assigned to new object.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
RandomFieldGenerator *CreateUsrDefRandomFieldGenerator(const char *name, int num, Domain *d);

/**
 * Creates new instance of user defined integration rule corresponding to given keyword.
 * @param type Id determining the type of new instance.
 * @param num  Component number.
 * @param e    Element to which new IR belongs.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
IntegrationRule *CreateUsrDefIRuleOfType(classType type, int num, Element *e);

/**
 * Creates new instance of element corresponding to given element keyword.
 * @param type Element id determining the type of new instance.
 * @param num  Element number.
 * @param d    Domain assigned to new element.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
Element *CreateUsrDefElementOfType(classType type, int num, Domain *d);
/**
 * Creates new instance of user defined dof manager corresponding to given keyword.
 * @param type Id determining the type of new instance.
 * @param num  Component number.
 * @param d    Domain assigned to new object.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
DofManager *CreateUsrDefDofManagerOfType(classType type, int num, Domain *d);
/**
 * Creates new instance of user defined dof corresponding to given keyword.
 * @param type Id determining the type of new instance.
 * @param num  Component number.
 * @param d    Dofmanager to which new dof belongs.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
Dof *CreateUsrDefDofOfType(classType type, int num, DofManager *d);
/**
 * Creates new instance of material mapping algorithm, corresponding to given MaterialMappingAlgorithmType.
 * @param type Id determining the type of new instance.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
MaterialMappingAlgorithm *CreateUsrDefMaterialMappingAlgorithm(MaterialMappingAlgorithmType type);
/**
 * Creates new instance of mesher interface, corresponding to given MeshPackageType.
 * @param type id determining the type of new instance.
 * @param d Domain assigned to new object.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
MesherInterface *CreateUsrDefMesherInterface(MeshPackageType type, Domain *d);
/**
 * Creates new instance of enrichment item.
 * @param name XFEM keyword.
 * @param num Component number.
 * @param xm XFEM manager which item belongs to.
 * @param d Domain assigned to new object.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
EnrichmentItem *CreateUsrDefEnrichmentItem(const char *name, int num, XfemManager *xm, Domain *d);
/**
 * Creates new instance of enrichment function.
 * @param name Enrichment function keyword.
 * @param num Component number.
 * @param d Domain assigned to new object.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
EnrichmentFunction *CreateUsrDefEnrichmentFunction(const char *name, int num, Domain *d);
/**
 * Creates new instance of geometry.
 * @param name Geometry keyword.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
BasicGeometry *CreateUsrDefGeometry(const char *name);
/**
 * Creates new instance of patch.
 * @param ptype Id determining the type of new instance.
 * @param e Parent element.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
Patch *CreateUsrDefPatch(Patch :: PatchType ptype, Element *e);
/**
 * Creates new instance nodal recovery model of given type.
 * @param type Id determining the type of new instance.
 * @param d domain to which newly created recovery model is associated to.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
NodalRecoveryModel *CreateUsrDefNodalRecoveryModel(NodalRecoveryModel :: NodalRecoveryModelType type, Domain *d);

#ifdef __PARALLEL_MODE
/**
 * Creates new instance of load balance monitor corresponding to given keyword.
 * @param type Id determining the type of new instance.
 * @param e    Engineering model to which new monitor belongs.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
LoadBalancerMonitor *CreateUsrDefLoadBalancerMonitorOfType(classType type, EngngModel *e);
/**
 * Creates new instance of load balancer corresponding to given keyword.
 * @param type Id determining the type of new instance.
 * @param d    Domain to which new balancer is attached.
 * @return Newly allocated object of requested type, null if keyword not supported.
 */
LoadBalancer *CreateUsrDefLoadBalancerOfType(classType type, Domain *d);
#endif
} // end namespace oofem
#endif // usrdefsub_h
