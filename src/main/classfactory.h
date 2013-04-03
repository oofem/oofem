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

#ifndef classfactory_h
#define classfactory_h

#include "classtype.h"
#include "sparsemtrxtype.h"
#include "errorestimatortype.h"
#include "doftype.h"
#include "linsystsolvertype.h"

#include <map>
#include <string>
#include <cstring>

namespace oofem {

class Domain;
class EngngModel;
class Element;
class DofManager;
class GeneralBoundaryCondition;
class CrossSection;
class Material;
class LoadTimeFunction;
class NonlocalBarrier;
class RandomFieldGenerator;

class Dof;
class SparseMtrx;
class SparseLinearSystemNM;
class ErrorEstimator;
class InitialCondition;

// Templates to wrap constructors into functions
template< typename T > Element *elemCreator(int n, Domain *d) { return new T(n, d); }
template< typename T > DofManager *dofmanCreator(int n, Domain *d) { return new T(n, d); }
template< typename T > GeneralBoundaryCondition *bcCreator(int n, Domain *d) { return new T(n, d); }
template< typename T > CrossSection *csCreator(int n, Domain *d) { return new T(n, d); }
template< typename T > Material *matCreator(int n, Domain *d) { return new T(n, d); }
template< typename T > EngngModel *engngCreator(int n, EngngModel *m) { return ( new T(n, m) ); }
template< typename T > LoadTimeFunction *ltfCreator(int n, Domain *d) { return new T(n, d); }
template< typename T > NonlocalBarrier *nlbCreator(int n, Domain *d) { return new T(n, d); }
template< typename T > RandomFieldGenerator *rfgCreator(int n, Domain *d) { return new T(n, d); }
//template < typename T > SparseMatrix* sparseMtrxCreator() { return new T(); }
//template < typename T > SparseMatrix* sparseMtrxCreator(DSSMatrix::dssType t) { return new T(t); }

/**
 * Class Factory allows to register terminal oofem classes, based on their membership
 * (classes representing elements, dof managers, material models, etc)
 * and create them on demand according to their name or id.
 * Global instance of ClassFactory named classFactory is created on startup.
 */
class ClassFactory
{
    struct CaseComp
    {
        int operator()(const std :: string &a, const std :: string &b) const;
    };

protected:
    /// Associative container containing element creators with element name as key.
    std :: map < std :: string, Element * ( * )(int, Domain *), CaseComp > elemNameList;
    /// Associative container containing element creators with element id as key.
    std :: map < classType, Element * ( * )(int, Domain *) > elemIdList;
    /// Associative container containing dofmanager creators with dofmanager  name as key.
    std :: map < std :: string, DofManager * ( * )(int, Domain *), CaseComp > dofmanNameList;
    /// Associative container containing dofmanager creators with dofmanager id as key.
    std :: map < classType, DofManager * ( * )(int, Domain *) > dofmanIdList;
    /// Associative container containing boundary condition creators with bc  name as key.
    std :: map < std :: string, GeneralBoundaryCondition * ( * )(int, Domain *), CaseComp > bcNameList;
    /// Associative container containing boundary condition creators with bc id as key.
    std :: map < classType, GeneralBoundaryCondition * ( * )(int, Domain *) > bcIdList;
    /// Associative container containing cross section creators with cross section name as key.
    std :: map < std :: string, CrossSection * ( * )(int, Domain *), CaseComp > csNameList;
    /// Associative container containing cross section creators with cross section id as key.
    std :: map < classType, CrossSection * ( * )(int, Domain *) > csIdList;
    /// Associative container containing material creators with material name as key.
    std :: map < std :: string, Material * ( * )(int, Domain *), CaseComp > matNameList;
    /// Associative container containing material creators with material id as key.
    std :: map < classType, Material * ( * )(int, Domain *) > matIdList;
    /// Associative container containing engng model creators with engng model name as key.
    std :: map < std :: string, EngngModel * ( * )(int, EngngModel *), CaseComp > engngNameList;
    /// Associative container containing engng model creators with engng model id as key.
    std :: map < classType, EngngModel * ( * )(int, EngngModel *) > engngIdList;
    /// Associative container containing load time function creators with ltf name as key.
    std :: map < std :: string, LoadTimeFunction * ( * )(int, Domain *), CaseComp > ltfNameList;
    /// Associative container containing load time function creators with ltf id as key.
    std :: map < classType, LoadTimeFunction * ( * )(int, Domain *) > ltfIdList;
    /// Associative container containing nonlocal barriers creators with barrier name as key.
    std :: map < std :: string, NonlocalBarrier * ( * )(int, Domain *), CaseComp > nlbNameList;
    /// Associative container containing nonlocal barriers creators with barrier id as key.
    std :: map < classType, NonlocalBarrier * ( * )(int, Domain *) > nlbIdList;
    /// Associative container containing random field generator creators with names as key.
    std :: map < std :: string, RandomFieldGenerator * ( * )(int, Domain *), CaseComp > rfgNameList;
    /// Associative container containing random field generator creators with class id as key.
    std :: map < classType, RandomFieldGenerator * ( * )(int, Domain *) > rfgIdList;

public:
    /// Constructor, registers all classes
    ClassFactory();

    /**
     * Creates new instance of element corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  Element number.
     * @param d    Domain assigned to new element.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    Element *createElement(const char *aClass, int number, Domain *domain);
    /**
     * Creates new instance of element corresponding to given element id.
     * @param type ClassId determining the type of new instance.
     * @param num  Element number.
     * @param d    Domain assigned to new element.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    Element *createElement(classType type, int number, Domain *domain);
    /**
     * Creates new instance of Dof manager corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  dofmanager number.
     * @param d    Domain assigned to new instance.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    DofManager *createDofManager(const char *aClass, int number, Domain *domain);
    /**
     * Creates new instance of dof manager corresponding to given id.
     * @param type ClassId determining the type of new instance.
     * @param num  object's number.
     * @param d    Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    DofManager *createDofManager(classType type, int number, Domain *domain);
    /**
     * Creates new instance of boundary condition corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  number of new object.
     * @param d    Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    GeneralBoundaryCondition *createBoundaryCondition(const char *aClass, int number, Domain *domain);
    /**
     * Creates new instance of boundary condition corresponding to given id.
     * @param type ClassId determining the type of new instance.
     * @param num  object's number.
     * @param d    Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    GeneralBoundaryCondition *createBoundaryCondition(classType type, int number, Domain *domain);
    /**
     * Creates new instance of cross section corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  object's number.
     * @param d    Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    CrossSection *createCrossSection(const char *aClass, int number, Domain *domain);
    /**
     * Creates new instance of cross section corresponding to given id.
     * @param type ClassId determining the type of new instance.
     * @param num  object's number.
     * @param d    Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    CrossSection *createCrossSection(classType type, int number, Domain *domain);
    /**
     * Creates new instance of material corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  material number.
     * @param d    Domain assigned to new material.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    Material *createMaterial(const char *aClass, int number, Domain *domain);
    /**
     * Creates new instance of material corresponding to given keyword.
     * @param type ClassId determining the type of new instance.
     * @param num  object's number.
     * @param d    Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    Material *createMaterial(classType type, int number, Domain *domain);
    /**
     * Creates new instance of engng model corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  Engng model number.
     * @param d    Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    EngngModel *createEngngModel(const char *aClass, int number, EngngModel *master);
    /**
     * Creates new instance of engng model corresponding to given id.
     * @param type ClassId determining the type of new instance.
     * @param num  object's number.
     * @param d    Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    EngngModel *createEngngModel(classType type, int number, EngngModel *master);
    /**
     * Creates new instance of load time function corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  object's number.
     * @param d    Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    LoadTimeFunction *createLoadTimeFunction(const char *aClass, int number, Domain *domain);
    /**
     * Creates new instance of load time function corresponding to given id.
     * @param type ClassId determining the type of new instance.
     * @param name Keyword string determining the type of new instance.
     * @param num  object's number.
     * @param d    Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    LoadTimeFunction *createLoadTimeFunction(classType type, int number, Domain *domain);
    /**
     * Creates new instance of nonlocal barrier corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  object's number.
     * @param d    Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    NonlocalBarrier *createNonlocalBarrier(const char *aClass, int number, Domain *domain);
    /**
     * Creates new instance of nonlocal barrier corresponding to given id.
     * @param type ClassId determining the type of new instance.
     * @param name Keyword string determining the type of new instance.
     * @param num  object's number.
     * @param d    Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    NonlocalBarrier *createNonlocalBarrier(classType type, int number, Domain *domain);
    /**
     * Creates new instance of random field generator corresponding to given keyword.
     * @param name Keyword string determining the type of new instance.
     * @param num  object's number.
     * @param d    Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    RandomFieldGenerator *createRandomFieldGenerator(const char *aClass, int number, Domain *domain);
    /**
     * Creates new instance of random field generator corresponding to given id.
     * @param type ClassId determining the type of new instance.
     * @param name Keyword string determining the type of new instance.
     * @param num  object's number.
     * @param d    Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    RandomFieldGenerator *createRandomFieldGenerator(classType type, int number, Domain *domain);
    /**
     * Creates new instance of sparse mtrx corresponding to given keyword.
     * @param type SparseMtrxType id determining the type of new instance.
     * @param num  object's number.
     * @param d    Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    SparseMtrx *createSparseMtrx(SparseMtrxType type);
    /**
     * Creates new instance of DOF corresponding to given keyword.
     * @param type classType id determining the type of new instance.
     * @param num  object's number.
     * @param dman Dof manager assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    Dof *createDof(dofType type, int num, DofManager *dman);
    /**
     * Creates new instance of SparseLinearSystemNM corresponding
     * to given type.
     * @param type LinSystSolverType id determining the type of new instance.
     * @param d Domain assigned to new object.
     * @param m EngngModel  assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    SparseLinearSystemNM *createSparseLinSolver(LinSystSolverType st, Domain *d, EngngModel *m);
    /**
     * Creates new instance of ErrorEstimator corresponding
     * to given type.
     * @param type ErrorEstimatorType id determining the type of new instance.
     * @param number  object's number.
     * @param d Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    ErrorEstimator *createErrorEstimator(ErrorEstimatorType type, int num, Domain *d);
    /**
     * Creates new instance of Initial Condition corresponding
     * to given type.
     * @param type classType id determining the type of new instance.
     * @param number  object's number.
     * @param d Domain assigned to new object.
     * @return Newly allocated object of requested type, null if keyword not supported.
     */
    InitialCondition *createInitialCondition(classType type, int num, Domain *d) ;

};

extern ClassFactory classFactory;
} // end namespace oofem
#endif // clasfactort_h
