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

#undef  REGISTER_CLASS
#define REGISTER_CLASS(_class, id)
#define REGISTER_CLASS_1(_class, id, type)
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
    ///@todo Is this even relevant to have here in the class factory? There is only one type available / Mikael
    if (type == InitialConditionClass) {
        return new InitialCondition(num, d);
    } else {
        return NULL;
    }
}

Element* ClassFactory::createElement ( const char* aClass, int number, Domain* domain )
{
    return ( elemNameList.count ( aClass ) == 1 ) ? elemNameList [ aClass ] ( number, domain ) : NULL;
}

Element* ClassFactory::createElement ( classType type, int number, Domain* domain )
{
    return ( elemIdList.count ( type ) == 1 ) ? elemIdList [ type ] ( number, domain ) : NULL;
}

DofManager* ClassFactory::createDofManager ( const char* aClass, int number, Domain* domain )
{
    return ( dofmanNameList.count ( aClass ) == 1 ) ? dofmanNameList [ aClass ] ( number, domain ) : NULL;
}

DofManager* ClassFactory::createDofManager ( classType type, int number, Domain* domain )
{
    return ( dofmanIdList.count ( type ) == 1 ) ? dofmanIdList [ type ] ( number, domain ) : NULL;
}

GeneralBoundaryCondition* ClassFactory::createBoundaryCondition ( const char* aClass, int number, Domain* domain )
{
    return ( bcNameList.count ( aClass ) == 1 ) ? bcNameList [ aClass ] ( number, domain ) : NULL;
}

GeneralBoundaryCondition* ClassFactory::createBoundaryCondition ( classType type, int number, Domain* domain )
{
    return ( bcIdList.count ( type ) == 1 ) ? bcIdList [ type ] ( number, domain ) : NULL;
}

CrossSection* ClassFactory::createCrossSection ( const char* aClass, int number, Domain* domain )
{
    return ( csNameList.count ( aClass ) == 1 ) ? csNameList [ aClass ] ( number, domain ) : NULL;
}

CrossSection* ClassFactory::createCrossSection ( classType type, int number, Domain* domain )
{
    return ( csIdList.count ( type ) == 1 ) ? csIdList [ type ] ( number, domain ) : NULL;
}

Material* ClassFactory::createMaterial ( const char* aClass, int number, Domain* domain )
{
    return ( matNameList.count ( aClass ) == 1 ) ? matNameList [ aClass ] ( number, domain ) : NULL;
}

Material* ClassFactory::createMaterial ( classType type, int number, Domain* domain )
{
    return ( matIdList.count ( type ) == 1 ) ? matIdList [ type ] ( number, domain ) : NULL;
}

EngngModel* ClassFactory::createEngngModel ( const char* aClass, int number, EngngModel* master )
{
    return ( engngNameList.count ( aClass ) == 1 ) ? engngNameList [ aClass ] ( number, master ) : NULL;
}

EngngModel* ClassFactory::createEngngModel ( classType type, int number, EngngModel* master )
{
    return ( engngIdList.count ( type ) == 1 ) ? engngIdList [ type ] ( number, master ) : NULL;
}

LoadTimeFunction* ClassFactory::createLoadTimeFunction ( const char* aClass, int number, Domain* domain )
{
    return ( ltfNameList.count ( aClass ) == 1 ) ? ltfNameList [ aClass ] ( number, domain ) : NULL;
}

LoadTimeFunction* ClassFactory::createLoadTimeFunction ( classType type, int number, Domain* domain )
{
    return ( ltfIdList.count ( type ) == 1 ) ? ltfIdList [ type ] ( number, domain ) : NULL;
}

NonlocalBarrier* ClassFactory::createNonlocalBarrier ( const char* aClass, int number, Domain* domain )
{
    return ( nlbNameList.count ( aClass ) == 1 ) ? nlbNameList [ aClass ] ( number, domain ) : NULL;
}

NonlocalBarrier* ClassFactory::createNonlocalBarrier ( classType type, int number, Domain* domain )
{
    return ( nlbIdList.count ( type ) == 1 ) ? nlbIdList [ type ] ( number, domain ) : NULL;
}

RandomFieldGenerator* ClassFactory::createRandomFieldGenerator ( const char* aClass, int number, Domain* domain )
{
    return ( rfgNameList.count ( aClass ) == 1 ) ? rfgNameList [ aClass ] ( number, domain ) : NULL;
}

RandomFieldGenerator* ClassFactory::createRandomFieldGenerator ( classType type, int number, Domain* domain )
{
    return ( rfgIdList.count ( type ) == 1 ) ? rfgIdList [ type ] ( number, domain ) : NULL;
}

} // End namespace oofem
