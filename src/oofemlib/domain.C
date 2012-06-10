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

#include "domain.h"
#include "element.h"
#include "timestep.h"
#include "node.h"
#include "elementside.h"
#include "material.h"
#include "crosssection.h"
#include "load.h"
#include "initial.h"
#include "loadtime.h"
#include "engngm.h"
#include "oofem_limits.h"
#include "entityrenumberingscheme.h"
#include "datastream.h"
#include "contextioerr.h"
#include "verbose.h"
#include "conTable.h"
#include "outputmanager.h"
#include "octreelocalizer.h"
#include "datareader.h"
#include "nodalrecoverymodel.h"
#include "nonlocalbarrier.h"
#include "usrdefsub.h"
#include "logger.h"
#include "xfemmanager.h"
#include "topologydescription.h"
#include "randomfieldgenerator.h"
#include "errorestimator.h"
#include "compiler.h"

#ifdef __PARALLEL_MODE
 #include "parallel.h"
 #include "processcomm.h"
 #include "datastream.h"
 #include "communicator.h"
 #include "domaintransactionmanager.h"
#endif

#ifndef __MAKEDEPEND
 #include <cstdarg>
 #include <cstring>
#endif

namespace oofem {
Domain :: Domain(int n, int serNum, EngngModel *e) : defaultNodeDofIDArry()
    // Constructor. Creates a new domain.
{
    this->engineeringModel = e;
    this->number = n;
    this->serialNumber = serNum;

    elementList           = new AList< Element >(0);
    dofManagerList        = new AList< DofManager >(0);
    materialList          = new AList< Material >(0);
    bcList                = new AList< GeneralBoundaryCondition >(0);
    icList                = new AList< InitialCondition >(0);
    loadTimeFunctionList  = new AList< LoadTimeFunction >(0);
    crossSectionList      = new AList< CrossSection >(0);
    nonlocalBarierList    = new AList< NonlocalBarrier >(0);
    randomFieldGeneratorList = new AList< RandomFieldGenerator >(0);
    xfemManager           = NULL;

    dType                 = _unknownMode;

    connectivityTable     = NULL;
    spatialLocalizer      = NULL;
    outputManager         = new OutputManager(this);
    smoother              = NULL;
    topology              = NULL;

    nonlocalUpdateStateCounter = 0;

#ifdef __PARALLEL_MODE
    dmanMapInitialized = elementMapInitialized = false;
    transactionManager = NULL;
#endif
}


Domain :: ~Domain()
// Destructor.
{
    delete elementList;
    delete dofManagerList;
    delete materialList;
    delete bcList;
    delete icList;
    delete loadTimeFunctionList;
    delete crossSectionList;
    delete nonlocalBarierList;
    delete randomFieldGeneratorList;
    delete xfemManager;
    delete connectivityTable;
    delete spatialLocalizer;
    delete outputManager;
    if ( smoother ) {
        delete smoother;
    }

    if ( topology ) {
        delete topology;
    }

#ifdef __PARALLEL_MODE
    if ( transactionManager ) {
        delete transactionManager;
    }

#endif
}

void
Domain :: clear()
// Clear receiver
{
    elementList->clear();
    dofManagerList->clear();
    materialList->clear();
    bcList->clear();
    icList->clear();
    loadTimeFunctionList->clear();
    crossSectionList->clear();
    nonlocalBarierList->clear();
    randomFieldGeneratorList->clear();

    if ( xfemManager ) {
        xfemManager->clear();
    }

    if ( connectivityTable ) {
        connectivityTable->reset();
    }

    if ( spatialLocalizer ) {
        delete spatialLocalizer;
        spatialLocalizer = NULL;
    }

    if ( smoother ) {
        smoother->clear();
    }

    // bp: how to clear/reset topology data?
    if ( topology ) {
        delete topology;
        topology = NULL;
    }

#ifdef __PARALLEL_MODE
    if ( transactionManager ) {
        delete transactionManager;
        transactionManager = NULL;
    }

#endif
}

Element *
Domain :: giveElement(int n)
// Returns the n-th element. Generates error if it is not defined yet.
{
    if ( elementList->includes(n) ) {
        return elementList->at(n);
    } else {
        _error2("giveElement: undefined element (%d)", n);
    }

    return NULL;
}


Load *
Domain :: giveLoad(int n)
// Returns the n-th load. Generates the error if not defined.
{
    Load *answer;

    if ( bcList->includes(n) ) {
        answer = dynamic_cast< Load * >( bcList->at(n) );
        if ( answer ) {
            return answer;
        } else {
            _error2("giveLoad: cannot cast boundary condition %d to Load class", n);
        }
    } else {
        _error2("giveLoad: undefined load (%d)", n);
    }

    return NULL;
}


GeneralBoundaryCondition *
Domain :: giveBc(int n)
// Returns the n-th bc. Generates the error if not defined.
{
    if ( bcList->includes(n) ) {
        return bcList->at(n);
    } else {
        _error2("giveBc: undefined bc (%d)", n);
    }

    return NULL;
}


InitialCondition *
Domain :: giveIc(int n)
// Returns the n-th ic. Generates the error if not defined.
{
    if ( icList->includes(n) ) {
        return icList->at(n);
    } else {
        _error2("giveIc: undefined ic (%d)", n);
    }

    return NULL;
}


LoadTimeFunction *
Domain :: giveLoadTimeFunction(int n)
// Returns the n-th load-time function. Creates this fuction if it does
// not exist yet.
{
    if ( loadTimeFunctionList->includes(n) ) {
        return loadTimeFunctionList->at(n);
    } else {
        _error2("giveLoadTimeFunction: undefined load-time function (%d)", n);
        //      ltf = (LoadTimeFunction*) LoadTimeFunction(n,this).typed() ;
        //      loadTimeFunctionList -> put(n,ltf) ;}
    }

    return NULL;
}


Material *
Domain :: giveMaterial(int n)
// Returns the n-th material. Creates this material if it does not exist
// yet.
{
    if ( materialList->includes(n) ) {
        return materialList->at(n);
    } else {
        _error2("giveMaterial: undefined material (%d)", n);
        //      mat = new Material(n,this) ;
        //      materialList  -> put(n,mat) ;}
    }

    return NULL;
}


Node *
Domain :: giveNode(int n)
// Returns the n-th node if it exists.
{
    Node *node = NULL;

#ifdef DEBUG
    if ( !dofManagerList->includes(n) ) {
        _error2("giveNode: undefined dofManager (%d)", n);
    }

#endif
    node = dynamic_cast< Node * >( dofManagerList->at(n) );
    if ( node == NULL ) {
        _error2("giveNode: incompatible type of dofManager %d, can not convert", n);
    }

    return node;
}


ElementSide *
Domain :: giveSide(int n)
// Returns the n-th element side.
{
    DofManager *side = NULL;

    if ( dofManagerList->includes(n) ) {
        side = dofManagerList->at(n);
        if ( side->giveClassID() != ElementSideClass ) {
            _error2("giveSide: incompatible type of dofManager %d, can not convert", n);
        }
    } else {
        _error2("giveSide: undefined dofManager (%d)", n);
    }

    return ( ElementSide * ) side;
}


DofManager *
Domain :: giveDofManager(int n)
// Returns the n-th node. Creates this node if it does not exist yet.
{
    if ( dofManagerList->includes(n) ) {
        return dofManagerList->at(n);
    } else {
        _error2("giveDofManager: undefined dofManager (%d)", n);
        //      node = new Node(n,this) ;
        //      nodeList  -> put(n,node) ;}
    }

    return NULL;
}


CrossSection *
Domain :: giveCrossSection(int n)
// Returns the n-th cross section.
// yet.
{
    if ( crossSectionList->includes(n) ) {
        return crossSectionList->at(n);
    } else {
        _error2("giveCrossSection: undefined cross section (%d)", n);
    }

    return NULL;
}


NonlocalBarrier *
Domain :: giveNonlocalBarrier(int n)
// Returns the n-th NonlocalBarrier.
{
    if ( nonlocalBarierList->includes(n) ) {
        return nonlocalBarierList->at(n);
    } else {
        _error2("giveNonlocalBarrier: undefined barrier (%d)", n);
    }

    return NULL;
}


RandomFieldGenerator *
Domain :: giveRandomFieldGenerator(int n)
// Returns the n-th RandomFieldGenerator.
{
    if ( randomFieldGeneratorList->includes(n) ) {
        return randomFieldGeneratorList->at(n);
    } else {
        _error2("giveRandomFieldGenerator: undefined generator (%d)", n);
    }

    return NULL;
}


EngngModel *
Domain :: giveEngngModel()
// Returns the time integration algorithm. Creates it if it does not
// exist yet.
{
    if ( engineeringModel ) {
        return engineeringModel;
    } else {
        _error("giveEngngModel: Not defined");
    }

    return NULL;
}

void Domain :: resizeDofManagers(int _newSize) { dofManagerList->growTo(_newSize); }
void Domain :: resizeElements(int _newSize) { elementList->growTo(_newSize); }
void Domain :: resizeCrossSectionModels(int _newSize) { crossSectionList->growTo(_newSize); }
void Domain :: resizeMaterials(int _newSize) { materialList->growTo(_newSize); }
void Domain :: resizeNonlocalBarriers(int _newSize) { nonlocalBarierList->growTo(_newSize); }
void Domain :: resizeBoundaryConditions(int _newSize) { bcList->growTo(_newSize); }
void Domain :: resizeInitialConditions(int _newSize) { icList->growTo(_newSize); }
void Domain :: resizeLoadTimeFunctions(int _newSize) { loadTimeFunctionList->growTo(_newSize); }
void Domain :: resizeRandomFieldGenerators(int _newSize) { randomFieldGeneratorList->growTo(_newSize); }

void Domain :: setDofManager(int i, DofManager *obj) { dofManagerList->put(i, obj); }
void Domain :: setElement(int i, Element *obj) { elementList->put(i, obj); }
void Domain :: setCrossSection(int i, CrossSection *obj) { crossSectionList->put(i, obj); }
void Domain :: setMaterial(int i, Material *obj) { materialList->put(i, obj); }
void Domain :: setNonlocalBarrier(int i, NonlocalBarrier *obj) { nonlocalBarierList->put(i, obj); }
void Domain :: setBoundaryCondition(int i, GeneralBoundaryCondition *obj) { bcList->put(i, obj); }
void Domain :: setInitialCondition(int i, InitialCondition *obj) { icList->put(i, obj); }
void Domain :: setLoadTimeFunction(int i, LoadTimeFunction *obj) { loadTimeFunctionList->put(i, obj); }
void Domain :: setRandomFieldGenerator(int i, RandomFieldGenerator *obj) { randomFieldGeneratorList->put(i, obj); }

void Domain :: clearBoundaryConditions() { bcList->clear(true); };

int
Domain :: instanciateYourself(DataReader *dr)
// Creates all objects mentioned in the data file.
{
    const char *__keyword, *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                            // Required by IR_GIVE_FIELD macro

    int i, num;
    std :: string name, topologytype;
    int nnode, nelem, nmat, nload, nic, nloadtimefunc, ncrossSections, nbarrier, nrfg;
    DofManager *node;
    Element *elem;
    Material *mat;
    GeneralBoundaryCondition *load;
    InitialCondition *ic;
    LoadTimeFunction *ltf;
    CrossSection *crossSection;
    NonlocalBarrier *barrier;
    RandomFieldGenerator *rfg;

#ifdef __ENABLE_COMPONENT_LABELS
    // mapping from label to local numbers for dofmans and elements
    std :: map< int, int >dofManLabelMap, elemLabelMap;
#endif


    FILE *outputStream = this->giveEngngModel()->giveOutputStream();

    // read type of Domain to be solved
    InputRecord *ir = dr->giveInputRecord(DataReader :: IR_domainRec, 1);
    __keyword = "domain";
    result = ir->giveField(name, IFT_Domain_type, __keyword);
    if ( result != IRRT_OK ) {
        IR_IOERR(giveClassName(), __proc, IFT_Domain_type, __keyword, ir, result);
    }

    ir->finish();

#  ifdef VERBOSE
    VERBOSE_PRINT0("Instanciating domain ", this->number);
#  endif

    resolveDomainDofsDefaults( name.c_str() );
    fprintf( outputStream, "Domain type: %s, default ndofs per node is %d\n\n\n",
            name.c_str(), giveDefaultNodeDofIDArry().giveSize() );

    // read output manager record
    ir = dr->giveInputRecord(DataReader :: IR_outManRec, 1);
    outputManager->initializeFrom(ir);
    ir->finish();

    // read domain description
    ir = dr->giveInputRecord(DataReader :: IR_domainCompRec, 1);
    IR_GIVE_FIELD(ir, nnode, IFT_Domain_ndofman, "ndofman"); // Macro
    IR_GIVE_FIELD(ir, nelem, IFT_Domain_nelem, "nelem"); // Macro
    IR_GIVE_FIELD(ir, ncrossSections, IFT_Domain_ncrosssect, "ncrosssect"); // Macro
    IR_GIVE_FIELD(ir, nmat, IFT_Domain_nmat, "nmat"); // Macro
    IR_GIVE_FIELD(ir, nload, IFT_Domain_nbc, "nbc"); // Macro
    IR_GIVE_FIELD(ir, nic, IFT_Domain_nic, "nic"); // Macro
    IR_GIVE_FIELD(ir, nloadtimefunc, IFT_Domain_nloadtimefunct, "nltf"); // Macro
    IR_GIVE_OPTIONAL_FIELD(ir, topologytype, IFT_Domain_topology, "topology"); // Macro

    // read optional number of nonlocalBarriers
    nbarrier = 0;
    __keyword = "nbarrier";
    result = ir->giveOptionalField(nbarrier,  IFT_Domain_nbarrier, __keyword);
    // read optional number of RandomFieldGenerator
    nrfg = 0;
    result = ir->giveOptionalField(nrfg,  IFT_Domain_nrfg, "nrandgen");



    // read nodes
    dofManagerList->growTo(nnode);
    for ( i = 0; i < nnode; i++ ) {
        ir = dr->giveInputRecord(DataReader :: IR_dofmanRec, i + 1);
        // read type of dofManager
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

#ifdef __ENABLE_COMPONENT_LABELS
        // assign component number according to record order
        // component number (as given in input record) becomes label
        if ( ( node = CreateUsrDefDofManagerOfType(name.c_str(), i+1, this) ) == NULL ) {
            OOFEM_ERROR2("Domain :: instanciateYourself - Couldn't create node of type: %s\n", name.c_str());
        }
        node->initializeFrom(ir);
        if ( dofManLabelMap.find(num) == dofManLabelMap.end() ) {
            // label does not exist yet
            dofManLabelMap [ num ] = i + 1;
        } else {
            _error2("instanciateYourself: Dofmanager entry already exist (label=%d)", num);
        }

        node->setGlobalNumber(num);    // set label
        dofManagerList->put(i + 1, node);
#else
        // component numbers as given in input record
        if ( ( node = CreateUsrDefDofManagerOfType(name.c_str(), num, this) ) == NULL ) {
            OOFEM_ERROR2("Domain :: instanciateYourself - Couldn't create node of type: %s\n", name.c_str());
        }
        node->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > nnode ) ) {
            _error2("instanciateYourself: Invalid dofManager number (num=%d)", num);
        }

        if ( !dofManagerList->includes(num) ) {
            dofManagerList->put(num, node);
        } else {
            _error2("instanciateYourself: Dofmanager entry already exist (num=%d)", num);
        }

#endif

        //dofManagerList->put(i+1,node) ;
        ir->finish();
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Instanciated nodes & sides ", nnode)
#  endif

    // read elements
    elementList->growTo(nelem);
    for ( i = 0; i < nelem; i++ ) {
        ir = dr->giveInputRecord(DataReader :: IR_elemRec, i + 1);
        // read type of element
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

#ifdef __ENABLE_COMPONENT_LABELS
        if ( ( elem = CreateUsrDefElementOfType(name.c_str(), i + 1, this) ) == NULL ) {
            OOFEM_ERROR2("Domain :: instanciateYourself - Couldn't create element: %s", name.c_str());
        }
        elem->initializeFrom(ir);

        if ( elemLabelMap.find(num) == elemLabelMap.end() ) {
            // label does not exist yet
            elemLabelMap [ num ] = i + 1;
        } else {
            _error2("instanciateYourself: Element entry already exist (label=%d)", num);
        }

        elem->setGlobalNumber(num);
        elementList->put(i + 1, elem);
#else
        if ( ( elem = CreateUsrDefElementOfType(name.c_str(), num, this) ) == NULL ) {
            OOFEM_ERROR2("Domain :: instanciateYourself - Couldn't create element: %s", name.c_str());
        }
        elem->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > nelem ) ) {
            _error2("instanciateYourself: Invalid element number (num=%d)", num);
        }

        if ( !elementList->includes(num) ) {
            elementList->put(num, elem);
        } else {
            _error2("instanciateYourself: element entry already exist (num=%d)", num);
        }

#endif
        //elementList->put(i+1,elem) ;
        ir->finish();
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Instanciated elements ", nelem);
#  endif

    // read cross sections
    crossSectionList->growTo(ncrossSections);
    for ( i = 0; i < ncrossSections; i++ ) {
        ir = dr->giveInputRecord(DataReader :: IR_crosssectRec, i + 1);
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        if ( ( crossSection = CreateUsrDefCrossSectionOfType(name.c_str(), num, this) ) == NULL ) {
            OOFEM_ERROR2("Domain :: instanciateYourself - Couldn't create crosssection: %s", name.c_str());
        }
        crossSection->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > ncrossSections ) ) {
            _error2("instanciateYourself: Invalid crossSection number (num=%d)", num);
        }

        if ( !crossSectionList->includes(num) ) {
            crossSectionList->put(num, crossSection);
        } else {
            _error2("instanciateYourself: crossSection entry already exist (num=%d)", num);
        }

        //crossSectionList->put(i+1,crossSection) ;
        ir->finish();
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Instanciated cross sections ", ncrossSections)
#  endif

    // read materials
    materialList->growTo(nmat);
    for ( i = 0; i < nmat; i++ ) {
        ir = dr->giveInputRecord(DataReader :: IR_matRec, i + 1);
        // read type of material
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        if ( ( mat = CreateUsrDefMaterialOfType(name.c_str(), num, this) ) == NULL ) {
            OOFEM_ERROR2("Domain :: instanciateYourself - Couldn't create material: %s", name.c_str());
        }
        mat->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > nmat ) ) {
            _error2("instanciateYourself: Invalid material number (num=%d)", num);
        }

        if ( !materialList->includes(num) ) {
            materialList->put(num, mat);
        } else {
            _error2("instanciateYourself: material entry already exist (num=%d)", num);
        }

        //materialList->put(i+1,mat) ;
        ir->finish();
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Instanciated materials ", nmat)
#  endif

    if ( nbarrier ) {
        // read barriers
        nonlocalBarierList->growTo(nbarrier);
        for ( i = 0; i < nbarrier; i++ ) {
            ir = dr->giveInputRecord(DataReader :: IR_nlocBarRec, i + 1);
            // read type of load
            IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

            barrier = CreateUsrDefNonlocalBarrierOfType(name.c_str(), num, this);
            barrier->initializeFrom(ir);

            // check number
            if ( ( num < 1 ) || ( num > nbarrier ) ) {
                _error2("instanciateYourself: Invalid barrier number (num=%d)", num);
            }

            if ( !nonlocalBarierList->includes(num) ) {
                nonlocalBarierList->put(num, barrier);
            } else {
                _error2("instanciateYourself: barrier entry already exist (num=%d)", num);
            }

            //nonlocalBarierList->put(i+1,barrier) ;
            ir->finish();
        }

#  ifdef VERBOSE
        VERBOSE_PRINT0("Instanciated barriers ", nbarrier);
#  endif
    }

    if ( nrfg ) {
        // read random field generators
        randomFieldGeneratorList->growTo(nrfg);
        for ( i = 0; i < nrfg; i++ ) {
            ir = dr->giveInputRecord(DataReader :: IR_nRandomFieldGenRec, i + 1);
            // read type of load
            IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

            rfg = CreateUsrDefRandomFieldGenerator(name.c_str(), num, this);
            rfg->initializeFrom(ir);

            // check number
            if ( ( num < 1 ) || ( num > nrfg ) ) {
                _error2("instanciateYourself: Invalid generator number (num=%d)", num);
            }

            if ( !randomFieldGeneratorList->includes(num) ) {
                randomFieldGeneratorList->put(num, rfg);
            } else {
                _error2("instanciateYourself: generator entry already exist (num=%d)", num);
            }

            ir->finish();
        }

#  ifdef VERBOSE
        VERBOSE_PRINT0("Instanciated random generators ", nbarrier);
#  endif
    }



    // read boundary conditions
    bcList->growTo(nload);
    for ( i = 0; i < nload; i++ ) {
        ir = dr->giveInputRecord(DataReader :: IR_bcRec, i + 1);
        // read type of load
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        if ( ( load = CreateUsrDefBoundaryConditionOfType(name.c_str(), num, this) ) == NULL ) {
            OOFEM_ERROR2("Domain :: instanciateYourself - Couldn't create boundary condition: %s", name.c_str());
        }
        load->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > nload ) ) {
            _error2("instanciateYourself: Invalid boundary condition number (num=%d)", num);
        }

        if ( !bcList->includes(num) ) {
            bcList->put(num, load);
        } else {
            _error2("instanciateYourself: boundary condition entry already exist (num=%d)", num);
        }

        //loadList->put(i+1,load) ;
        ir->finish();
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Instanciated BCs ", nload)
#  endif

    // read initial conditions
    icList->growTo(nic);
    for ( i = 0; i < nic; i++ ) {
        ir = dr->giveInputRecord(DataReader :: IR_icRec, i + 1);
        // read type of load
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        ic = new InitialCondition(num, this);
        if ( ic ) {
            ic->initializeFrom(ir);
        } else {
            _error2("instanciateYourself: Creation of IC no. %d failed", num);
        }

        // check number
        if ( ( num < 1 ) || ( num > nic ) ) {
            _error2("instanciateYourself: Invalid initial condition number (num=%d)", num);
        }

        if ( !icList->includes(num) ) {
            icList->put(num, ic);
        } else {
            _error2("instanciateYourself: initial condition entry already exist (num=%d)", num);
        }

        ir->finish();
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Instanciated ICs ", nic)
#  endif


    // read load time functions
    loadTimeFunctionList->growTo(nloadtimefunc);
    for ( i = 0; i < nloadtimefunc; i++ ) {
        ir = dr->giveInputRecord(DataReader :: IR_ltfRec, i + 1);
        // read type of ltf
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);
        if ( ( ltf = CreateUsrDefLoadTimeFunctionOfType(name.c_str(), num, this) ) == NULL ) {
            OOFEM_ERROR2("Domain :: instanciateYourself - Couldn't create time function: %s", name.c_str());
        }
        ltf->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > nloadtimefunc ) ) {
            _error2("instanciateYourself: Invalid LoadTimeFunction number (num=%d)", num);
        }

        if ( !loadTimeFunctionList->includes(num) ) {
            loadTimeFunctionList->put(num, ltf);
        } else {
            _error2("instanciateYourself: LoadTimeFunction entry already exist (num=%d)", num);
        }

        //loadTimeFunctionList->put(i+1,ltf) ;
        ir->finish();
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Instanciated load-time fncts ", nloadtimefunc)
#  endif


#ifdef __ENABLE_COMPONENT_LABELS
    // change internal component references from labels to assigned local numbers
    MapBasedEntityRenumberingFunctor labelToLocNumFunctor(dofManLabelMap, elemLabelMap);
    for ( i = 1; i <= nnode; i++ ) {
        this->giveDofManager(i)->updateLocalNumbering(labelToLocNumFunctor);
    }

    for ( i = 1; i <= nelem; i++ ) {
        this->giveElement(i)->updateLocalNumbering(labelToLocNumFunctor);
    }

#endif
    this->topology = NULL;
    if ( topologytype.length() > 0 ) {
        this->topology = CreateUsrDefTopologyOfType(topologytype.c_str(), this);
        if ( !this->topology ) {
            OOFEM_ERROR2( "Domain :: instanciateYourself - Couldn't create topology of type '%s'", topologytype.c_str() );
        }

        return this->topology->instanciateYourself(dr);
    }

    return 1;
}


void
Domain :: error(const char *file, int line, const char *format, ...)
{
    char buffer [ MAX_ERROR_MSG_LENGTH ];
    va_list args;

    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);

    __OOFEM_ERROR3(file, line, "Class: Domain, number: %d\n%s", number, buffer);
}


void
Domain :: warning(const char *file, int line, const char *format, ...)
{
    char buffer [ MAX_ERROR_MSG_LENGTH ];
    va_list args;

    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);

    __OOFEM_WARNING3(file, line, "Class: Domain, number: %d\n%s", number, buffer);
}


const IntArray &
Domain :: giveDefaultNodeDofIDArry()
{
    // returns default DofID array, defining physical meaning of particular DOFs
    // in Node Dof collection
    if ( this->defaultNodeDofIDArry.giveSize() ) {
        return defaultNodeDofIDArry;
    }

    if ( dType == _2dPlaneStressRotMode ) {
        defaultNodeDofIDArry.setValues(3, D_u, D_v, R_w);
    } else if ( dType == _2dPlaneStressMode ) {
        defaultNodeDofIDArry.setValues(2, D_u, D_v);
    } else if ( dType == _PlaneStrainMode ) {
        defaultNodeDofIDArry.setValues(2, D_u, D_v);
    } else if  ( dType == _3dMode ) {
        defaultNodeDofIDArry.setValues(3, D_u, D_v, D_w);
    } else if ( dType == _3dAxisymmMode ) {
        defaultNodeDofIDArry.setValues(3, D_u, D_v, R_w);
    } else if  ( dType == _2dMindlinPlateMode ) {
        defaultNodeDofIDArry.setValues(3, D_w, R_u, R_v);
    } else if ( dType == _3dShellMode ) {
        defaultNodeDofIDArry.setValues(6, D_u, D_v, D_w, R_u, R_v, R_w);
    } else if  ( dType == _2dTrussMode ) {
        defaultNodeDofIDArry.setValues(2, D_u, D_w);
    } else if  ( dType == _1dTrussMode ) {
        defaultNodeDofIDArry.setValues(1, D_u);
    } else if  ( dType == _2dBeamMode ) {
        defaultNodeDofIDArry.setValues(3, D_u, D_w, R_v);
    } else if  ( dType == _2dLatticeMode ) {
        defaultNodeDofIDArry.setValues(3, D_u, D_v, R_w);
    } else if  ( dType == _HeatTransferMode ) {
        defaultNodeDofIDArry.setValues(1, T_f);
    } else if  ( dType == _HeatMass1Mode ) {
        defaultNodeDofIDArry.setValues(2, T_f, C_1);
    }  else if ( dType == _2dIncompressibleFlow ) {
        defaultNodeDofIDArry.setValues(3, V_u, V_v, P_f);
    }  else if ( dType == _3dIncompressibleFlow ) {
        defaultNodeDofIDArry.setValues(4, V_u, V_v, V_w, P_f);
    } else {
        _error2( "giveDefaultNodeDofIDArry : unknown domainType (%s)", __domainTypeToString(dType) );
    }

    return defaultNodeDofIDArry;
}


int
Domain :: giveNumberOfSpatialDimensions()
{
    //_HeatTransferMode _HeatMass1Mode // Are these deprecated?
    // Perhaps i shouldn't use the modes to determine this at all, but i couldn't see any other good way.
    if ( dType == _1dTrussMode ) {
        return 1;
    }

    if ( dType == _2dIncompressibleFlow || dType == _2dBeamMode || dType == _2dTrussMode || dType == _2dMindlinPlateMode || dType == _3dAxisymmMode || dType == _PlaneStrainMode || dType == _2dPlaneStressMode || dType == _2dPlaneStressRotMode ) {
        return 2;
    } else if ( dType == _3dIncompressibleFlow || dType == _3dShellMode || dType == _3dMode ) {
        return 3;
    } else {
        return 0;
    }
}


void
Domain :: resolveDomainDofsDefaults(const char *typeName)
//
// resolves default number of dofs per node according to domain type name.
// and also resolves default dof mask according to domain type.
//
{
    if ( !strncasecmp(typeName, "2dplanestressrot", 16) ) {
        dType = _2dPlaneStressRotMode;
    } else if ( !strncasecmp(typeName, "2dplanestress", 12) ) {
        dType = _2dPlaneStressMode;
    } else if ( !strncasecmp(typeName, "planestrain", 11) ) {
        dType = _PlaneStrainMode;
    } else if ( !strncasecmp(typeName, "3daxisymm", 9) ) {
        dType = _3dAxisymmMode;
    } else if  ( !strncasecmp(typeName, "2dmindlinplate", 14) ) {
        dType = _2dMindlinPlateMode;
    } else if ( !strncasecmp(typeName, "3dshell", 7) ) {
        dType = _3dShellMode;
    } else if  ( !strncasecmp(typeName, "2dtruss", 7) ) {
        dType = _2dTrussMode;
    } else if  ( !strncasecmp(typeName, "1dtruss", 7) ) {
        dType = _1dTrussMode;
    } else if  ( !strncasecmp(typeName, "2dbeam", 6) ) {
        dType = _2dBeamMode;
    } else if  ( !strncasecmp(typeName, "2dlattice", 9) ) {
        dType = _2dLatticeMode;
    } else if  ( !strncasecmp(typeName, "heattransfer", 11) )   {
        dType = _HeatTransferMode;
    } else if  ( !strncasecmp(typeName, "hema1", 5) ) {
        dType = _HeatMass1Mode;
    } else if ( !strncasecmp(typeName, "2dincompflow", 12) ) {
        dType = _2dIncompressibleFlow;
    } else if ( !strncasecmp(typeName, "3dincompflow", 12) ) {
        dType = _3dIncompressibleFlow;
    } else if  ( !strncasecmp(typeName, "3d", 2) ) {
        dType = _3dMode;
    } else {
        _error2("resolveDomainDofsDefaults : unknown domainType (%s)", typeName);
        return;
    }
}


#ifdef __OOFEG

void
Domain :: drawYourself(oofegGraphicContext &context)
//
// shows graphics representation of domain, respecting mode
//
{
    OGC_PlotModeType plotMode = context.giveIntVarPlotMode();
    if ( ( plotMode == OGC_nodeAnnotation ) || ( plotMode == OGC_nodeGeometry ) || ( plotMode == OGC_essentialBC ) ||
        ( plotMode == OGC_naturalBC ) || ( plotMode == OGC_nodeScalarPlot ) || ( plotMode == OGC_nodeVectorPlot ) ) {
        this->drawNodes(context);
    } else {
        this->drawElements(context);
    }
}


void
Domain :: drawElements(oofegGraphicContext &context)
{
    //
    // steps through element array and calls element(i)->show(mode,this);
    //
    for ( int i = 1; i <= this->giveNumberOfElements(); i++ ) {
        this->giveElement(i)->drawYourself(context);
    }
}


void
Domain :: drawNodes(oofegGraphicContext &context)
{
    //
    // steps through element array and calls element(i)->show(mode,this);
    //
    int nnodes = this->giveNumberOfDofManagers();
    for ( int i = 1; i <= nnodes; i++ ) {
        this->giveDofManager(i)->drawYourself(context);
    }
}

#endif


NodalRecoveryModel *
Domain :: giveSmoother()
{
    return this->smoother;
}


void
Domain :: setSmoother(NodalRecoveryModel *smoother, bool destroyOld)
{
    if ( destroyOld && this->smoother ) {
        delete this->smoother;
    }

    this->smoother = smoother;
}


void
Domain :: setTopology(TopologyDescription *topo, bool destroyOld)
{
    if ( destroyOld && this->topology ) {
        delete this->topology;
    }

    this->topology = topo;
}


ConnectivityTable *
Domain :: giveConnectivityTable()
//
// return connectivity Table - if no defined - creates new one
//
{
    if ( connectivityTable == NULL ) {
        connectivityTable = new ConnectivityTable(this);
    }

    return connectivityTable;
}


SpatialLocalizer *
Domain :: giveSpatialLocalizer()
//
// return connectivity Table - if no defined - creates new one
//
{
    //  if (spatialLocalizer == NULL) spatialLocalizer = new DummySpatialLocalizer(1, this);
    if ( spatialLocalizer == NULL ) {
        spatialLocalizer = new OctreeSpatialLocalizer(1, this);
    }

    return spatialLocalizer;
}


int
Domain ::  giveCorrespondingCoordinateIndex(int idof)
//
// find corresponding coordinate axis to idof
// if no - coordinate axis corespond to idof returns 0;
//
// if idof corresponds to displacement in direction of axis i then finction returns i
// otherwise 0;
//
{
    switch ( dType ) {
    case _2dBeamMode:
        if ( idof == 1 ) {
            return 1;
        } else if ( idof == 2 ) {
            return 3;
        }

        return 0;

    case _2dPlaneStressMode:
    case _PlaneStrainMode:
    case _2dLatticeMode:
        if ( idof == 1 ) {
            return 1;
        } else if ( idof == 2 ) {
            return 2;
        }

        return 0;

    case _2dPlaneStressRotMode:
        if ( idof == 1 ) {
            return 1;
        } else if ( idof == 2 ) {
            return 2;
        }

        return 0;

    case _2dTrussMode:
        if ( idof == 1 ) {
            return 1;
        } else if ( idof == 2 ) {
            return 3;
        }

        return 0;

    case _1dTrussMode:
        if ( idof == 1 ) {
            return 1;
        }

        return 0;

    case _2dMindlinPlateMode:
        if ( idof == 1 ) {
            return 3;
        }

        return 0;

    case _3dMode:
        if ( idof == 1 ) {
            return 1;
        } else if ( idof == 2 ) {
            return 2;
        } else if ( idof == 3 ) {
            return 3;
        }

        return 0;

    case _3dAxisymmMode:
        if ( idof == 1 ) {
            return 1;
        } else if ( idof == 2 ) {
            return 2;
        }

        return 0;

    case _3dShellMode:
        if ( idof == 1 ) {
            return 1;
        } else if ( idof == 2 ) {
            return 2;
        } else if ( idof == 3 ) {
            return 3;
        }

        return 0;

    case _3dIncompressibleFlow:
        if ( idof == 1 ) {
            return 1;
        } else if ( idof == 2 ) {
            return 2;
        } else if ( idof == 3 ) {
            return 3;
        }

        return 0;

    case _2dIncompressibleFlow:
        if ( idof == 1 ) {
            return 1;
        } else if ( idof == 2 ) {
            return 2;
        }

        return 0;

    default:
        _error("giveCorrespondingCoordinateIndex : unsupported domain type");
    }

    return 0;
}


int
Domain :: checkConsistency()
// this function transverse tree of all objects and invokes
// checkConsistency on this objects
// currently this function checks noly consistency
// of internal object structures, mainly whether referenced other objects
// are having required support
//
{
    int i, result = 1;
    int nnode, nelem, nmat;

    nnode = this->giveNumberOfDofManagers();
    nelem = this->giveNumberOfElements();
    nmat  = this->giveNumberOfMaterialModels();

    for ( i = 1; i <= nnode; i++ ) {
        result &= this->giveDofManager(i)->checkConsistency();
    }

    for ( i = 1; i <= nelem; i++ ) {
        result &= this->giveElement(i)->checkConsistency();
    }

    for ( i = 1; i <= nmat; i++ ) {
        result &= this->giveMaterial(i)->checkConsistency();
    }

    return result;
}

double Domain :: giveArea()
{
    double area = 0.0;
    for ( int i = 1; i <= this->giveNumberOfElements(); ++i ) {
        area += this->giveElement(i)->computeArea();
    }

    return area;
}

double Domain :: giveVolume()
{
    double volume = 0.0;
    for ( int i = 1; i <= this->giveNumberOfElements(); ++i ) {
        volume += this->giveElement(i)->computeVolume();
    }

    return volume;
}

ErrorEstimator *
Domain :: giveErrorEstimator() {
    return engineeringModel->giveDomainErrorEstimator(this->number);
}



#define SAVE_COMPONENTS(size,type,giveMethod)   \
  {                                             \
    type* obj;                                  \
    for ( i = 1; i <= size; i++ ) {             \
        obj = giveMethod(i);                    \
        if ( ( mode & CM_Definition ) ) {       \
            ct =  (int) obj->giveClassID();     \
            if ( !stream->write(& ct, 1) ) {    \
                THROW_CIOERR(CIO_IOERR);        \
            }                                   \
        }                                       \
        if ( ( iores = obj->saveContext(stream, mode) ) != CIO_OK ) { \
            THROW_CIOERR(iores);                \
        }                                       \
    }                                           \
  }

#define RESTORE_COMPONENTS(size,type,resizeMethod,creator,giveMethod,setMethod) \
  {                                         \
    type *obj;                              \
    if ( mode & CM_Definition ) {           \
        resizeMethod(size);                 \
    }                                       \
    for ( i = 1; i <= size; i++ ) {         \
        if ( mode & CM_Definition ) {       \
            if ( !stream->read(& ct, 1) ) { \
                THROW_CIOERR(CIO_IOERR);    \
            }                               \
            compId = ( classType ) ct;      \
            obj = creator(compId, 0, this); \
        } else {                            \
            obj = giveMethod(i);            \
        }                                   \
        if ( ( iores = obj->restoreContext(stream, mode) ) != CIO_OK ) { \
            THROW_CIOERR(iores);            \
        }                                   \
        if ( mode & CM_Definition ) {       \
            setMethod(i, obj);              \
        }                                   \
    }                                       \
  }

#define DOMAIN_NCOMP 9

contextIOResultType
Domain :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    int i, serNum;
    ErrorEstimator *ee;
    int ct;

    // save domain serial number
    serNum = this->giveSerialNumber();
    if ( !stream->write(& serNum, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( mode & CM_Definition ) ) {
      long ncomp[DOMAIN_NCOMP];
      ncomp[0]=this->giveNumberOfDofManagers();
      ncomp[1]=this->giveNumberOfElements();
      ncomp[2]=this->giveNumberOfMaterialModels();
      ncomp[3]=this->giveNumberOfCrossSectionModels();
      ncomp[4]=this->giveNumberOfBoundaryConditions();
      ncomp[5]=this->giveNumberOfInitialConditions();
      ncomp[6]=this->giveNumberOfLoadTimeFunctions();
      ncomp[7]=this->giveNumberOfNonlocalBarriers();
      ncomp[8]=this->giveNumberOfRandomFieldGenerators();

      // store number of components
      if ( !stream->write (ncomp, DOMAIN_NCOMP) ) {
        THROW_CIOERR(CIO_IOERR);
      }
    }

    // save dof managers
    SAVE_COMPONENTS(this->giveNumberOfDofManagers(),DofManager,this->giveDofManager);
    // elements and corresponding integration points
    SAVE_COMPONENTS(this->giveNumberOfElements(),Element,this->giveElement);
    // boundary conditions
    SAVE_COMPONENTS(this->giveNumberOfBoundaryConditions(),GeneralBoundaryCondition,this->giveBc);

    if ( ( mode & CM_Definition ) ) {
      // store material models
      SAVE_COMPONENTS(this->giveNumberOfMaterialModels(),Material,this->giveMaterial);
      // store cross sections
      SAVE_COMPONENTS(this->giveNumberOfCrossSectionModels(),CrossSection,this->giveCrossSection);
      // store initial conditions
      SAVE_COMPONENTS(this->giveNumberOfInitialConditions(),InitialCondition,this->giveIc);
      // store load time functions
      SAVE_COMPONENTS(this->giveNumberOfLoadTimeFunctions(),LoadTimeFunction,this->giveLoadTimeFunction);
      // store nonlocal barriers
      SAVE_COMPONENTS(this->giveNumberOfNonlocalBarriers(),NonlocalBarrier,this->giveNonlocalBarrier);
      // store random field generators
      SAVE_COMPONENTS(this->giveNumberOfRandomFieldGenerators(),RandomFieldGenerator,this->giveRandomFieldGenerator);
    } // end if ( ( mode & CM_Definition ) ) {


    // store error estimator data
    ee = this->giveErrorEstimator();
    if ( ee ) {
        if ( ( iores = ee->saveContext(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


contextIOResultType
Domain :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    int i, serNum;
    bool domainUpdated = false;
    ErrorEstimator *ee;
    int ct;
    classType compId;
    long ncomp[DOMAIN_NCOMP];

    int nnodes, nelem, nmat, ncs, nbc, nic, nltf, nnlb, nrfg;


    domainUpdated = false;
    serNum = this->giveSerialNumber();
    // restore domain serial number
    if ( !stream->read(& this->serialNumber, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( mode & CM_Definition ) ) {
        // read number of components
        if ( !stream->read(ncomp, DOMAIN_NCOMP) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        nnodes = ncomp[0];
        nelem = ncomp[1];
        nmat  = ncomp[2];
        ncs   = ncomp[3];
        nbc   = ncomp[4];
        nic   = ncomp[5];
        nltf  = ncomp[6];
        nnlb  = ncomp[7];
        nrfg  = ncomp[8];

        // clear receiver data
        dofManagerList->clear();
        elementList->clear();
        materialList->clear();
        bcList->clear();
        icList->clear();
        loadTimeFunctionList->clear();
        nonlocalBarierList->clear();
        randomFieldGeneratorList->clear();
        //this->clear();


        domainUpdated = true;
    } else {
        if ( serNum != this->giveSerialNumber() ) {
            // read corresponding domain
            OOFEM_LOG_INFO( "restoring domain %d.%d\n", this->number, this->giveSerialNumber() );
            DataReader *domainDr = this->engineeringModel->GiveDomainDataReader(1, this->giveSerialNumber(), contextMode_read);
            this->clear();

            if ( !this->instanciateYourself(domainDr) ) {
                _error("initializeAdaptive: domain Instanciation failed");
            }

            delete domainDr;
            domainUpdated = true;
        }

        nnodes=this->giveNumberOfDofManagers();
        nelem=this->giveNumberOfElements();
        nmat=this->giveNumberOfMaterialModels();
        ncs=this->giveNumberOfCrossSectionModels();
        nbc=this->giveNumberOfBoundaryConditions();
        nic=this->giveNumberOfInitialConditions();
        nltf=this->giveNumberOfLoadTimeFunctions();
        nnlb=this->giveNumberOfNonlocalBarriers();
        nrfg=this->giveNumberOfRandomFieldGenerators();
    }

    RESTORE_COMPONENTS(nnodes,DofManager,this->resizeDofManagers,CreateUsrDefDofManagerOfType,this->giveDofManager,this->setDofManager);
    RESTORE_COMPONENTS(nelem,Element,this->resizeElements,CreateUsrDefElementOfType,this->giveElement,this->setElement);
    RESTORE_COMPONENTS(nbc,GeneralBoundaryCondition,this->resizeBoundaryConditions,CreateUsrDefBoundaryConditionOfType,this->giveBc,this->setBoundaryCondition);

    if ( ( mode & CM_Definition ) ) {

        RESTORE_COMPONENTS(nmat,Material,this->resizeMaterials,CreateUsrDefMaterialOfType,this->giveMaterial,this->setMaterial);
        RESTORE_COMPONENTS(ncs,CrossSection,this->resizeCrossSectionModels,CreateUsrDefCrossSectionOfType,this->giveCrossSection,this->setCrossSection);
        RESTORE_COMPONENTS(nic,InitialCondition,this->resizeInitialConditions,CreateUsrDefInitialConditionOfType,this->giveIc,setInitialCondition);
        RESTORE_COMPONENTS(nltf,LoadTimeFunction,resizeLoadTimeFunctions,CreateUsrDefLoadTimeFunctionOfType,giveLoadTimeFunction,setLoadTimeFunction);
        RESTORE_COMPONENTS(nnlb,NonlocalBarrier,resizeNonlocalBarriers,CreateUsrDefNonlocalBarrierOfType, giveNonlocalBarrier,setNonlocalBarrier);
        RESTORE_COMPONENTS(nrfg,RandomFieldGenerator,resizeRandomFieldGenerators,CreateUsrDefRandomFieldGenerator,giveRandomFieldGenerator,setRandomFieldGenerator);
    }

    // restore error estimator data
    ee = this->giveErrorEstimator();
    if ( domainUpdated ) {
        ee->setDomain(this);
    }

    if ( ee ) {
        if ( ( iores = ee->restoreContext(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    if ( domainUpdated ) {
        if ( this->smoother ) {
            this->smoother->init();
        }
    }

    return CIO_OK;
}




#ifdef __PARALLEL_MODE

DomainTransactionManager *
Domain :: giveTransactionManager()
{
    if ( !transactionManager ) {
        if ( ( transactionManager = new DomainTransactionManager(this) ) == NULL ) {
            OOFEM_ERROR("Domain::giveTransactionManager: allocation failed");
        }
    }

    return transactionManager;
}


int Domain :: commitTransactions(DomainTransactionManager *tm)
{
    bool _exist;
    std :: map< int, FEMComponent * > :: const_iterator it;
    AList< DofManager > *dofManagerList_new = new AList< DofManager >(0);
    AList< Element > *elementList_new = new AList< Element >(0);


    if ( tm->dofmanTransactions.empty() && tm->elementTransactions.empty() ) {
        return 1;
    }

    this->initGlobalDofManMap();
    if ( !tm->dofmanTransactions.empty() ) {
        DofManager *dman;
        for ( it = tm->dofmanTransactions.begin(); it != tm->dofmanTransactions.end(); ++it ) {
            _exist = false;
            if ( dmanMap.find(it->first) != dmanMap.end() ) {
                _exist = true;
            }

            if ( _exist ) {
                int lnum = dmanMap [ it->first ]->giveNumber();
                dman = dofManagerList->unlink(lnum);
                dmanMap.erase(it->first);
                delete dman;
            }

            if ( it->second ) {
                dmanMap [ it->first ] = ( DofManager * ) it->second;
            }
        } // end loop over DofmanTransactions

    }

    this->initGlobalElementMap();
    if ( !tm->elementTransactions.empty() ) {
        int gen;
        Element *elem;

        for ( it = tm->elementTransactions.begin(); it != tm->elementTransactions.end(); ++it ) {
            gen = it->first;
            bool _exist = false;
            if ( elementMap.find(gen) != elementMap.end() ) {
                _exist = true;
            }

            if ( _exist ) {
                int lnum = elementMap [ gen ]->giveNumber();
                elem = elementList->unlink(lnum);
                elementMap.erase(gen);
                delete elem;
            }

            if ( it->second ) {
                elementMap [ gen ] = ( Element * ) it->second;
            }
        }
    }

    /*
     * if (tm->transactions.empty()) return 1;
     *
     * AList<DofManager> *dofManagerList_new = new AList<DofManager>(0) ;
     * AList<Element>    *elementList_new    = new AList<Element>(0) ;
     *
     * // put existing domain dofman and element records into domain maps
     * this->initGlobalDofManMap ();
     * this->initGlobalElemMap();
     *
     * // commit all transactions
     * while (!tm->transactions.empty()) {
     * DTM_Transaction& t = tm->transactions.front();
     * if (t._ttype == DTT_Remove) {
     *  if (t._ctype == DCT_DofManager) {
     *    dmanMap.erase (t._num);
     *    dman = dofManagerList->unlink (t._num);
     *    delete dman;
     *  } else if (t._ctype == DCT_Element) {
     *    elem = elementList->unlink (t._num);
     *    delete elem;
     *  } else {
     *    OOFEM_ERROR ("Domain::commitTransactions: unknown transaction component type");
     *  }
     * } else if (t._ttype == DTT_ADD) {
     *  if (t._ctype == DCT_DofManager) {
     *    dmanMap[t._num] = (*DofManager) t._obj;
     *  } else if (t._ctype == DCT_Element) {
     *    elemMap[t._num] = (*Element) t._obj;
     *  } else {
     *    OOFEM_ERROR ("Domain::commitTransactions: unknown transaction component type");
     *  }
     * } else {
     *  OOFEM_ERROR ("Domain::commitTransactions: unknown transaction type");
     * }
     *
     * // Pop new transaction
     * tm->transactions.pop_front();
     * } // while (!tm->transactions.empty()) {
     *
     */

    this->renumberDofManagers();
    this->renumberDofManData(tm);

    // initialize new dofman list
    int _i, _size = dmanMap.size();
    std :: map< int, DofManager * > :: iterator dmit;
    dofManagerList_new->clear();
    dofManagerList_new->growTo(_size);

    for ( _i = 0, dmit = dmanMap.begin(); dmit != dmanMap.end(); dmit++ ) {
        dofManagerList_new->put(++_i, dmit->second);
    }

    this->renumberElements();
    this->renumberElementData(tm);
    // initialize new element list
    _size = elementMap.size();
    std :: map< int, Element * > :: iterator elit;
    elementList_new->clear();
    elementList_new->growTo(_size);

    for ( _i = 0, elit = elementMap.begin(); elit != elementMap.end(); elit++ ) {
        elementList_new->put(++_i, elit->second);
    }


    tm->dofmanTransactions.clear();
    tm->elementTransactions.clear();


    this->dofManagerList->clear(false); // not the data
    delete dofManagerList;
    this->dofManagerList = dofManagerList_new;

    this->elementList->clear(false);
    delete elementList;
    this->elementList = elementList_new;

    this->giveConnectivityTable()->reset();
    this->giveSpatialLocalizer()->init(true);
    return 1;
}


void
Domain :: initGlobalDofManMap(bool forceinit)
{
    /*
     * Renumber here the master node number for rigid and hanging dofs, etc;
     * existing local nodes need mapping from old_local to new numbering,
     * but received nodes need mapping from global to new numbering
     *
     * -> we need to keep the list of received nodes! (now they are only introduced into globally indexed dmanMap!);
     */
    if ( forceinit || !dmanMapInitialized ) {
        int key, idofman, ndofman = this->giveNumberOfDofManagers();
        DofManager *dofman;
        dmanMap.clear();

        for ( idofman = 1; idofman <= ndofman; idofman++ ) {
            dofman = this->giveDofManager(idofman);
            key = dofman->giveGlobalNumber();
            dmanMap [ key ] = dofman;
        }
    }
}


void
Domain :: initGlobalElementMap(bool forceinit)
{
    if ( forceinit || !elementMapInitialized ) {
        int key, ielem, nelem = this->giveNumberOfElements();
        Element *elem;
        elementMap.clear();

        for ( ielem = 1; ielem <= nelem; ielem++ ) {
            elem = this->giveElement(ielem);
            key = elem->giveGlobalNumber();
            elementMap [ key ] = elem;
        }
    }
}


void
Domain :: renumberDofManData(DomainTransactionManager *tm)
{
    std :: map< int, DofManager * > :: iterator it;

    SpecificEntityRenumberingFunctor< Domain > domainGToLFunctor(this, &Domain :: LB_giveUpdatedGlobalNumber);
    SpecificEntityRenumberingFunctor< Domain > domainLToLFunctor(this, &Domain :: LB_giveUpdatedLocalNumber);


    for ( it = dmanMap.begin(); it != dmanMap.end(); it++ ) {
        if ( tm->dofmanTransactions.find(it->first) != tm->dofmanTransactions.end() ) {
            // received dof manager -> we map global numbers to new local number
            it->second->updateLocalNumbering(domainGToLFunctor); // g_to_l
        } else {
            // existing dof manager -> we map old local number to new local number
            it->second->updateLocalNumbering(domainLToLFunctor); // l_to_l
        }
    }
}


void
Domain :: renumberElementData(DomainTransactionManager *tm)
{
    std :: map< int, Element * > :: iterator it;

    SpecificEntityRenumberingFunctor< Domain > domainGToLFunctor(this, &Domain :: LB_giveUpdatedGlobalNumber);
    SpecificEntityRenumberingFunctor< Domain > domainLToLFunctor(this, &Domain :: LB_giveUpdatedLocalNumber);


    for ( it = elementMap.begin(); it != elementMap.end(); it++ ) {
        if ( tm->elementTransactions.find(it->first) != tm->elementTransactions.end() ) {
            // received dof manager -> we map global numbers to new local number
            it->second->updateLocalNumbering(domainGToLFunctor); // g_to_l
        } else {
            // existing dof manager -> we map old local number to new local number
            it->second->updateLocalNumbering(domainLToLFunctor); // l_to_l
        }
    }
}


void
Domain :: renumberDofManagers()
{
    int _locnum;
    std :: map< int, DofManager * > :: iterator it;

    for ( _locnum = 1, it = dmanMap.begin(); it != dmanMap.end(); it++ ) {
        it->second->setNumber(_locnum++);
    }
}


void
Domain :: renumberElements()
{
    int _locnum;
    std :: map< int, Element * > :: iterator it;

    for ( _locnum = 1, it = elementMap.begin(); it != elementMap.end(); it++ ) {
        it->second->setNumber(_locnum++);
    }
}


int
Domain :: LB_giveUpdatedLocalNumber(int num, EntityRenumberingScheme scheme)
{
    if ( scheme == ERS_DofManager ) {
        DofManager *dm = this->giveDofManager(num);
        if ( dm ) {
            return dm->giveNumber();
        } else {
            _error2("LB_giveUpdatedLocalNumber: dofman %d moved to remote partition, updated number not available", num);
        }
    } else {
        _error("LB_giveUpdatedLocalNumber: unsuported renumbering scheme");
    }

    return 0;
}


int
Domain :: LB_giveUpdatedGlobalNumber(int num, EntityRenumberingScheme scheme)
{
    if ( scheme == ERS_DofManager ) {
        DofManager *dm = dmanMap [ num ];
        if ( dm ) {
            return dm->giveNumber();
        } else {
            _error2("LB_giveUpdatedGlobalNumber: dofman [%d] not available on local partition, updated number not available", num);
            return 0;
        }
    } else {
        _error("LB_giveUpdatedGlobalNumber: unsuported renumbering scheme");
    }

    return 0;
}


int
Domain :: dofmanGlobal2Local(int _globnum)
{
    if ( dmanMap.find(_globnum) != dmanMap.end() ) {
        // dofman is already available -> update only
        return ( dmanMap [ _globnum ]->giveNumber() );
    } else {
        return 0;
    }
}


int
Domain :: elementGlobal2Local(int _globnum)
{
    if ( elementMap.find(_globnum) != elementMap.end() ) {
        // element is already available -> update only
        return ( elementMap [ _globnum ]->giveNumber() );
    } else {
        return 0;
    }
}


void
Domain :: setXfemManager(XfemManager *xfemManager)
{
    this->xfemManager = xfemManager;
}



#endif
} // end namespace oofem
