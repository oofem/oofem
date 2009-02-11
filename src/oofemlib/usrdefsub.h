/* $Header: /home/cvs/bp/oofem/oofemlib/src/usrdefsub.h,v 1.14 2003/05/19 13:03:58 bp Exp $ */
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

//   ********************************************
//   file: usrdefsub.h - User defined Subroutines
//   ********************************************


#ifndef usrdefsub_h
#define usrdefsub_h

#include "domain.h"
#include "element.h"
#include "crosssection.h"
#include "engngm.h"
#include "load.h"
#include "loadtime.h"
#include "material.h"
#include "sparselinsystemnm.h"
#include "sparsegeneigenvalsystemnm.h"
#include "materialmappingalgorithm.h"
#include "materialmappingalgorithmtype.h"
#include "meshpackagetype.h"
#include "mesherinterface.h"
//#include "nonlocalbarrier.h"
#include "randomfieldgenerator.h"
#include "classtype.h"
#include "sparsemtrxtype.h"
#include "geneigvalsolvertype.h"
#include "errorestimatortype.h"
#include "enrichmentitem.h"

#include "errorestimator.h"
#ifdef __PARALLEL_MODE
 #include "loadbalancer.h"
#endif

/**
 * Creates new instance of element corresponding to given element keyword.
 * @param name element keyword string determining the type of new instance
 * @param num  element number
 * @param d    domain assigned to new element
 * @return newly allocated object of requested type, null if keyword not suppported
 */
Element *CreateUsrDefElementOfType(char *, int, Domain *);
/**
 * Creates new instance of user defined dof manager corresponding to given keyword.
 * @param name keyword string determining the type of new instance
 * @param num  component number
 * @param d    domain assigned to new object
 * @return newly allocated object of requested type, null if keyword not suppported
 */
DofManager *CreateUsrDefDofManagerOfType(char *, int, Domain *);
/**
 * Creates new instance of user defined cross section model corresponding to given keyword.
 * @param name keyword string determining the type of new instance
 * @param num  component number
 * @param d    domain assigned to new object
 * @return newly allocated object of requested type, null if keyword not suppported
 */
CrossSection *CreateUsrDefCrossSectionOfType(char *, int, Domain *);
/**
 * Creates new instance of user defined engng model corresponding to given keyword.
 * @param name keyword string determining the type of new instance
 * @param number  component number
 * @return newly allocated object of requested type, null if keyword not suppported
 */
EngngModel *CreateUsrDefEngngModelOfType(char *, int, EngngModel *_master = NULL);
/**
 * Creates new instance of user defined load  corresponding to given keyword.
 * @param name keyword string determining the type of new instance
 * @param num  component number
 * @param d    domain assigned to new object
 * @return newly allocated object of requested type, null if keyword not suppported
 */
GeneralBoundaryCondition *CreateUsrDefBoundaryConditionOfType(char *, int, Domain *);
/**
 * Creates new instance of user defined load time function corresponding to given keyword.
 * @param name keyword string determining the type of new instance
 * @param num  component number
 * @param d    domain assigned to new object
 * @return newly allocated object of requested type, null if keyword not suppported
 */
LoadTimeFunction *CreateUsrDefLoadTimeFunctionOfType(char *, int, Domain *);
/**
 * Creates new instance of user defined material model corresponding to given keyword.
 * @param name keyword string determining the type of new instance
 * @param num  component number
 * @param d    domain assigned to new object
 * @return newly allocated object of requested type, null if keyword not suppported
 */
Material *CreateUsrDefMaterialOfType(char *, int, Domain *);

/**
 * Creates new Empty Instance of Sparse matrix of given type (the constructor with no parameters is called).
 * @param type determines sparseMtrx type
 * @return newly allocated object of requested type, null if keyword not suppported
 */
SparseMtrx *CreateUsrDefSparseMtrx(SparseMtrxType type);
/**
 * Creates new Instance of Linear Sparse Solver of given type (the constructor with no parameters is called).
 * @param st solver type
 * @param i cpomponent number
 * @param d domain
 * @param m emodel
 * @return newly allocated object of requested type, null if keyword not suppported
 */
SparseLinearSystemNM *CreateUsrDefSparseLinSolver(LinSystSolverType st, int i, Domain *d, EngngModel *m);
/**
 * Creates new Instance of Generalized Eigen Value Solver of given type (the constructor with no parameters is called).
 * @param st solver type
 * @param i cpomponent number
 * @param d domain
 * @param m emodel
 * @return newly allocated object of requested type, null if keyword not suppported
 */
SparseGeneralEigenValueSystemNM *CreateUsrDefGeneralizedEigenValueSolver(GenEigvalSolverType st, int i, Domain *d, EngngModel *m);

/**
 * Creates new Instance of Error Estimator of given type.
 * @param type determines Error Estimator type
 * @param number EE number
 * @param d domain associated to ee
 * @return newly allocated object of requested type, null if keyword not suppported
 */
ErrorEstimator *CreateUsrDefErrorEstimator(ErrorEstimatorType type, int number, Domain *d);

/**
 * Creates new Instance of Export Module of given name.
 * @param name determines Export module type
 * @param emodel engn model associated to new export module
 * @return newly allocated object of requested type, null if keyword not suppported
 */
ExportModule *CreateUsrDefExportModuleOfType(char *name, EngngModel *emodel);


/**
 * Creates new Instance of Nonlocal Barrier class corresponding to given name.
 * @param name determines Nonlocal Barrier type
 * @param num  component number
 * @param d    domain assigned to new object
 * @return newly allocated object of requested type, null if keyword not suppported
 */
NonlocalBarrier *CreateUsrDefNonlocalBarrierOfType(char *name, int num, Domain *d);
/**
 * Creates new Instance of Random generator class corresponding to given name.
 * @param name determines random generator type
 * @param num  component number
 * @param d    domain assigned to new object
 * @return newly allocated object of requested type, null if keyword not suppported
 */
RandomFieldGenerator *CreateUsrDefRandomFieldGenerator(char *name, int num, Domain *d);

/**
 * Creates new instance of user defined integration rule corresponding to given keyword.
 * @param type id determining the type of new instance
 * @param num  component number
 * @param d    dofmanager to which new dof belongs
 * @return newly allocated object of requested type, null if keyword not suppported
 */
IntegrationRule *CreateUsrDefIRuleOfType(classType type, int, Element *);

/**
 * Creates new instance of element corresponding to given element keyword.
 * @param type element id determining the type of new instance
 * @param num  element number
 * @param d    domain assigned to new element
 * @return newly allocated object of requested type, null if keyword not suppported
 */
Element *CreateUsrDefElementOfType(classType type, int, Domain *);
/**
 * Creates new instance of user defined dof manager corresponding to given keyword.
 * @param type id determining the type of new instance
 * @param num  component number
 * @param d    domain assigned to new object
 * @return newly allocated object of requested type, null if keyword not suppported
 */
DofManager *CreateUsrDefDofManagerOfType(classType type, int, Domain *);
/**
 * Creates new instance of user defined dof corresponding to given keyword.
 * @param type id determining the type of new instance
 * @param num  component number
 * @param d    dofmanager to which new dof belongs
 * @return newly allocated object of requested type, null if keyword not suppported
 */
Dof *CreateUsrDefDofOfType(classType type, int, DofManager *);
/**
 * Creates new instance of matrial mapping algorithm, corresponding to given MaterialMappingAlgorithmType.
 * @param type id determining the type of new instance
 * @return newly allocated object of requested type, null if keyword not suppported
 */
MaterialMappingAlgorithm *CreateUsrDefMaterialMappingAlgorithm(MaterialMappingAlgorithmType type);
/**
 * Creates new instance of mesher interface, corresponding to given MeshPackageType.
 * @param type id determining the type of new instance
 * @return newly allocated object of requested type, null if keyword not suppported
 */
MesherInterface *CreateUsrDefMesherInterface(MeshPackageType type, Domain *d);
/**
 * Creates new instance of enrichment item.
 * @param type id determining the type of new instance
 * @return newly allocated object of requested type, null if keyword not suppported
 */
EnrichmentItem *CreateUsrDefEnrichmentItem(char *aClass, int num, Domain *d);
/**
 * Creates new instance of enrichment function.
 * @param type id determining the type of new instance
 * @return newly allocated object of requested type, null if keyword not suppported
 */
EnrichmentFunction *CreateUsrDefEnrichmentFunction(char *aClass, int num, Domain *d);
/**
 * Creates new instance of geometry.
 * @param type id determining the type of new instance
 * @return newly allocated object of requested type, null if keyword not suppported
 */
BasicGeometry *CreateUsrDefGeometry(char *aClass);
#ifdef __PARALLEL_MODE
/**
 * Creates new instance of load balance monitor corresponding to given keyword.
 * @param type id determining the type of new instance
 * @param e    engng model to which new monitor belongs
 * @return newly allocated object of requested type, null if keyword not suppported
 */
LoadBalancerMonitor *CreateUsrDefLoadBalancerMonitorOfType(classType type, EngngModel *);
/**
 * Creates new instance of load balancer corresponding to given keyword.
 * @param type id determining the type of new instance
 * @param d    domain to which new balancer is attached
 * @return newly allocated object of requested type, null if keyword not suppported
 */
LoadBalancer *CreateUsrDefLoadBalancerOfType(classType type, Domain *);
#endif
#endif // usrdefsub_h
