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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "domain.h"
#include "parametermanager.h"
#include "paramkey.h"

#include "element.h"
#include "timestep.h"
#include "node.h"
#include "elementside.h"
#include "material.h"
#include "crosssection.h"
#include "load.h"
#include "initialcondition.h"
#include "function.h"
#include "set.h"
#include "engngm.h"
#include "entityrenumberingscheme.h"
#include "datastream.h"
#include "contextioerr.h"
#include "verbose.h"
#include "connectivitytable.h"
#include "outputmanager.h"
#include "octreelocalizer.h"
#include "nodalrecoverymodel.h"
#include "nonlocalbarrier.h"
#include "classfactory.h"
#include "logger.h"
#include "xfem/xfemmanager.h"
#include "topologydescription.h"
#include "errorestimator.h"
#include "range.h"
#include "fracturemanager.h"
#include "datareader.h"
#include "oofemtxtdatareader.h"
#include "initmodulemanager.h"
#include "exportmodulemanager.h"
#include "xfem/enrichmentitem.h"
#include "xfem/nucleationcriterion.h"
#include "xfem/enrichmentfunction.h"
#include "xfem/propagationlaw.h"
#include "Contact/contactsurface.h"
#include "bctracker.h"

#include "boundarycondition.h"
#include "activebc.h"
#include "simpleslavedof.h"
#include "masterdof.h"


#ifdef __MPI_PARALLEL_MODE
 #include "parallel.h"
 #include "processcomm.h"
 #include "datastream.h"
 #include "communicator.h"
 #include "domaintransactionmanager.h"
#endif

#include <cstdarg>
#include <cstring>
#include <vector>
#include <set>

namespace oofem {
Domain :: Domain(int n, int serNum, EngngModel *e) : elementPPM(),
                                                    dofmanPPM(),
                                                    defaultNodeDofIDArry(),
                                                    bcTracker(this)
                                                     
    // Constructor. Creates a new domain.
{
    if ( !e->giveSuppressOutput() ) {
        outputManager = std::make_unique<OutputManager>(this);
    }

    this->engineeringModel = e;
    this->number = n;
    this->serialNumber = serNum;

    dType = _unknownMode;

    nonlocalUpdateStateCounter = 0;

    nsd = 0;
    axisymm = false;
    freeDofID = MaxDofID;

#ifdef __MPI_PARALLEL_MODE
    dmanMapInitialized = elementMapInitialized = false;
    transactionManager = NULL;
#endif
}

Domain :: ~Domain() { }

void
Domain :: clear()
// Clear receiver
{
    elementList.clear();
    elementGlobal2LocalMap.clear();
    dofmanGlobal2LocalMap.clear();
    dofManagerList.clear();
    materialList.clear();
    bcList.clear();
    icList.clear();
    functionList.clear();
    crossSectionList.clear();
    nonlocalBarrierList.clear();
    setList.clear();
    contactSurfaceList.clear();
    xfemManager = nullptr;
    if ( connectivityTable ) {
        connectivityTable->reset();
    }

    spatialLocalizer = nullptr;

    if ( smoother ) {
        smoother->clear();
    }

    ///@todo bp: how to clear/reset topology data?
    topology = nullptr;

#ifdef __MPI_PARALLEL_MODE
    transactionManager = nullptr;
#endif
}


Element *
Domain :: giveElement(int n)
{
#ifdef DEBUG
    if ( n < 1 || n > (int)elementList.size() ) {
        OOFEM_ERROR("undefined element (%d)", n);
    }
#endif
    return this->elementList[n-1].get();
}

Element *
Domain :: giveGlobalElement(int n)
{
    for ( auto &el: elementList ) {
        if ( el->giveGlobalNumber() == n ) {
            return el.get();
        }
    }

    return NULL;
}

int
Domain :: giveElementPlaceInArray(int iGlobalElNum) const
{
    auto res = elementGlobal2LocalMap.find(iGlobalElNum);

    if ( res != elementGlobal2LocalMap.end() ) {
        return res->second;
    } else {
        OOFEM_ERROR("returning -1 for iGlobalElNum: %d.", iGlobalElNum );
        return -1;
    }
}

int
Domain :: giveDofManPlaceInArray(int iGlobalDofManNum) const
{
    auto res = dofmanGlobal2LocalMap.find(iGlobalDofManNum);

    if ( res != dofmanGlobal2LocalMap.end() ) {
        return res->second;
    } else {
        OOFEM_ERROR("returning -1 for iGlobalDofManNum: %d.", iGlobalDofManNum );
        return -1;
    }
}

const IntArray &
Domain :: giveElementsWithMaterialNum(int iMaterialNum) const
{
    auto res = materialNum2ElMap.find(iMaterialNum);

    if ( res != materialNum2ElMap.end() ) {
        return res->second;
    } else {
        OOFEM_ERROR("Material not found.")
    }
}

Load *
Domain :: giveLoad(int n)
{
#ifdef DEBUG
    if ( n < 1 || n > (int)bcList.size() ) {
        OOFEM_ERROR("undefined load (%d)", n);
    }
    Load *answer = dynamic_cast< Load * >( bcList[n-1].get() );
    if ( answer ) {
        return answer;
    } else {
        OOFEM_ERROR("cannot cast boundary condition %d to Load class", n);
    }
#else
    return static_cast< Load * >( bcList[n-1].get() );

#endif
}


GeneralBoundaryCondition *
Domain :: giveBc(int n)
{
#ifdef DEBUG
    if ( n < 1 || n > (int)bcList.size() ) {
        OOFEM_ERROR("undefined bc (%d)", n);
    }
#endif
    return bcList[n-1].get();
}


InitialCondition *
Domain :: giveIc(int n)
{
#ifdef DEBUG
    if ( n < 1 || n > (int)icList.size() ) {
        OOFEM_ERROR("undefined ic (%d)", n);
    }
#endif

    return icList[n-1].get();
}


Function *
Domain :: giveFunction(int n)
{
#ifdef DEBUG
    if ( n < 1 || n > (int)functionList.size() ) {
        OOFEM_ERROR("undefined load-time function (%d)", n);
    }
#endif

    return functionList[n-1].get();
}


Material *
Domain :: giveMaterial(int n)
{
#ifdef DEBUG
    if ( n < 1 || n > (int)materialList.size() ) {
        OOFEM_ERROR("undefined material (%d)", n);
    }
#endif

    return materialList[n-1].get();
}

ElementSide *
Domain :: giveSide(int n)
{
#ifdef DEBUG
    if ( n < 1 || n > (int)dofManagerList.size() ) {
        OOFEM_ERROR("undefined dofManager (%d)", n);
    }

    ElementSide *side = dynamic_cast< ElementSide * >( dofManagerList[n-1].get() );
    if ( !side ) {
        OOFEM_ERROR("incompatible type of dofManager %d, can not convert", n);
    }
    return side;

#else
    return static_cast< ElementSide * >( dofManagerList[n-1].get() );

#endif
}


DofManager *
Domain :: giveDofManager(int n)
{
#ifdef DEBUG
    if ( n < 1 || n > (int)dofManagerList.size() ) {
        OOFEM_ERROR("undefined dofManager (%d)", n);
    }
#endif
    return this->dofManagerList[n-1].get();
}


DofManager *
Domain :: giveGlobalDofManager(int n)
{
    for ( auto &dman: dofManagerList ) {
        if ( dman->giveGlobalNumber() == n ) {
            return dman.get();
        }
    }

    return NULL;
}


CrossSection *
Domain :: giveCrossSection(int n)
{
#ifdef DEBUG
    if ( n < 1 || n > (int)crossSectionList.size() ) {
        OOFEM_ERROR("undefined cross section (%d)", n);
    }
#endif
    return crossSectionList[n-1].get();
}


NonlocalBarrier *
Domain :: giveNonlocalBarrier(int n)
{
#ifdef DEBUG
    if ( n < 1 || n > (int)nonlocalBarrierList.size() ) {
        OOFEM_ERROR("undefined barrier (%d)", n);
    }
#endif
    return nonlocalBarrierList[n-1].get();
}


Set *
Domain :: giveSet(int n)
{
#ifdef DEBUG
    if ( n < 1 || n > (int)setList.size() ) {
        OOFEM_ERROR("undefined set (%d)", n);
    }
#endif

    return setList[n-1].get();
}

XfemManager *
Domain :: giveXfemManager()
{
#ifdef DEBUG
    if ( !xfemManager ) {
        OOFEM_ERROR("undefined xfem manager");
    }
#endif
    return xfemManager.get();
}

bool
Domain :: hasXfemManager()
{
    return xfemManager.get() != NULL;
}


bool
Domain :: hasFractureManager()
{
    return fracManager.get() != NULL;
}

FractureManager *
Domain :: giveFractureManager()
{
#ifdef DEBUG
    if ( !fracManager ) {
        OOFEM_ERROR("undefined fracture manager");
    }
#endif
    return fracManager.get();
}

BCTracker*
Domain::giveBCTracker()
{
  return &bcTracker;
}
  
EngngModel *
Domain :: giveEngngModel()
{
#ifdef DEBUG
    if ( !engineeringModel ) {
        OOFEM_ERROR("Not defined");
    }
#endif

    return engineeringModel;
}


ContactSurface *
Domain::giveContactSurface(int n)
{
#ifdef DEBUG
    if ( n < 1 || n > ( int ) contactSurfaceList.size() ) {
        OOFEM_ERROR("undefined contact surface (%d)", n);
    }
#endif
    return this->contactSurfaceList [ n - 1 ].get();
}




void Domain :: resizeDofManagers(int _newSize) { dofManagerList.resize(_newSize); }
void Domain :: resizeElements(int _newSize) { elementList.resize(_newSize); }
void Domain :: resizeCrossSectionModels(int _newSize) { crossSectionList.resize(_newSize); }
void Domain :: resizeMaterials(int _newSize) { materialList.resize(_newSize); }
void Domain :: resizeNonlocalBarriers(int _newSize) { nonlocalBarrierList.resize(_newSize); }
void Domain :: resizeBoundaryConditions(int _newSize) { bcList.resize(_newSize); }
void Domain :: resizeInitialConditions(int _newSize) { icList.resize(_newSize); }
void Domain :: resizeFunctions(int _newSize) { functionList.resize(_newSize); }
void Domain :: resizeSets(int _newSize) { setList.resize(_newSize); }
void Domain :: resizeContactSurfaces(int _newSize) { contactSurfaceList.resize(_newSize); }

void Domain :: py_setDofManager(int i, DofManager *obj) { dofManagerList[i-1].reset(obj); dofmanGlobal2LocalMap[obj->giveGlobalNumber()] = i;}
void Domain :: py_setElement(int i, Element *obj) { elementList[i-1].reset(obj); elementGlobal2LocalMap[obj->giveGlobalNumber()] = i;}
void Domain :: py_setCrossSection(int i, CrossSection *obj) { crossSectionList[i-1].reset(obj); }
void Domain :: py_setMaterial(int i, Material *obj) { materialList[i-1].reset(obj); }
void Domain :: py_setNonlocalBarrier(int i, NonlocalBarrier *obj) { nonlocalBarrierList[i-1].reset(obj); }
void Domain :: py_setBoundaryCondition(int i, GeneralBoundaryCondition *obj) { bcList[i-1].reset(obj); }
void Domain :: py_setInitialCondition(int i, InitialCondition *obj) { icList[i-1].reset(obj); }
void Domain :: py_setFunction(int i, Function *obj) { functionList[i-1].reset(obj); }
void Domain :: py_setSet(int i, Set *obj) { setList[i-1].reset(obj); }

void Domain :: setDofManager(int i, std::unique_ptr<DofManager> obj) { dofmanGlobal2LocalMap[obj->giveGlobalNumber()] = i; dofManagerList[i-1] = std::move(obj); }
void Domain :: setElement(int i, std::unique_ptr<Element> obj) { elementGlobal2LocalMap[obj->giveGlobalNumber()] = i; elementList[i-1] = std::move(obj); }
void Domain :: setCrossSection(int i, std::unique_ptr<CrossSection> obj) { crossSectionList[i-1] = std::move(obj); }
void Domain :: setMaterial(int i, std::unique_ptr<Material> obj) { materialList[i-1] = std::move(obj); }
void Domain :: setNonlocalBarrier(int i, std::unique_ptr<NonlocalBarrier> obj) { nonlocalBarrierList[i-1] = std::move(obj); }
void Domain :: setBoundaryCondition(int i, std::unique_ptr<GeneralBoundaryCondition> obj) { bcList[i-1] = std::move(obj); }
void Domain :: setInitialCondition(int i, std::unique_ptr<InitialCondition> obj) { icList[i-1] = std::move(obj); }
void Domain :: setFunction(int i, std::unique_ptr<Function> obj) { functionList[i-1] = std::move(obj); }
void Domain :: setSet(int i, std::unique_ptr<Set> obj) { setList[i-1] = std::move(obj); }
void Domain :: setXfemManager(std::unique_ptr<XfemManager> obj) { xfemManager = std::move(obj); }

void Domain :: clearBoundaryConditions() { bcList.clear(); }
void Domain :: clearElements() { elementList.clear(); }

int Domain::instanciateYourself(DataReader& dr){
    auto &ir = dr.giveInputRecord(DataReader :: IR_domainRec, 1);
    return this->instanciateYourself(dr,ir);
}
int
Domain :: instanciateYourself(DataReader &dr, InputRecord& irDomain)
// Creates all objects mentioned in the data file.
{
    int num;
    std :: string name, topologytype;
    int componentRecPriority = 2;

    // read type of Domain to be solved
    {
        // ·auto &ir = dr.giveInputRecord(DataReader :: IR_domainRec, 1);
        IR_GIVE_FIELD(irDomain, name, _IFT_Domain_type); // This is inconsistent, "domain" isn't  exactly a field, but the actual record keyword.

        mDomainType = name;

        // irDomain.finish();
    }

    #  ifdef VERBOSE
        VERBOSE_PRINT0("Instanciating domain ", this->number);
    #  endif

    resolveDomainDofsDefaults( name.c_str() );

    // read output manager record
    {
        std :: string tmp;
        auto &ir = dr.giveInputRecord(DataReader :: IR_outManRec, 1);
        ir.giveRecordKeywordField(tmp);

        if ( !giveEngngModel()->giveSuppressOutput() ) {
            outputManager->initializeFrom(ir);
        }
        ir.finish();
    }

    // XML format (and perhaps others) does not contain DomainComp nested group, rather nests everything under domain directly
    std::shared_ptr<InputRecord> irdPtr(dr.hasFlattenedStructure()?irDomain.ptr():dr.giveInputRecord(DataReader :: IR_domainCompRec, 1).clone());
    IR_GIVE_OPTIONAL_FIELD(*irdPtr.get(), topologytype, _IFT_Domain_topology);
    this->nsd = -1; ///@todo Change this to default 0 when the domaintype record has been removed.
    IR_GIVE_OPTIONAL_FIELD(*irdPtr.get(), this->nsd, _IFT_Domain_numberOfSpatialDimensions);
    this->axisymm = irdPtr->hasField(_IFT_Domain_axisymmetric);


    // read nodes
    DataReader::GroupRecords dofManagerRecs=dr.giveGroupRecords(irdPtr,_IFT_Domain_ndofman,"Nodes",DataReader::IR_dofmanRec,/*optional*/false);
    DataReader::RecordGuard scope(dr,irdPtr.get()); // till the end of the scope
    dofManagerList.clear();
    dofManagerList.resize(dofManagerRecs.size());
    int dofManagerIndex0=0;
    for(auto& ir: dofManagerRecs){
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        // assign component number according to record order
        // component number (as given in input record) becomes label
        std :: unique_ptr< DofManager > dman( classFactory.createDofManager(name.c_str(), dofManagerIndex0+1, this) );
        if ( !dman ) {
            OOFEM_ERROR("Couldn't create node of type: %s\n", name.c_str());
        }

        dman->initializeFrom(ir, componentRecPriority);
        dman->setGlobalNumber(num);    // set label
        dofManagerList[dofManagerIndex0] = std :: move(dman);

        ir.finish();
        dofManagerIndex0++;
    }

    #  ifdef VERBOSE
        VERBOSE_PRINT0("Instantiated nodes & sides ", (int)dofManagerList.size())
    #  endif

    BuildDofManPlaceInArrayMap();

    // read elements
    DataReader::GroupRecords elemRecs=dr.giveGroupRecords(irdPtr,_IFT_Domain_nelem,"Elements",DataReader::IR_elemRec,/*optional*/false);
    elementList.clear();
    elementList.resize(elemRecs.size());
    int elemIndex0=0;
    for(auto& ir: elemRecs){
        // read type of element
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        std :: unique_ptr< Element >elem( classFactory.createElement(name.c_str(), elemIndex0+1, this) );
        if ( !elem ) {
            OOFEM_ERROR("Couldn't create element: %s", name.c_str());
        }
        elem->initializeFrom(ir, componentRecPriority);
        elem->setGlobalNumber(num);
        elementList[elemIndex0] = std :: move(elem);

        ir.finish();
        elemIndex0++;
    }

    BuildElementPlaceInArrayMap();

    // Support sets defined directly after the elements (special hack for backwards compatibility).
    setList.clear();
    if ( dr.peekNext("set") ) {
        DataReader::GroupRecords setRecs=dr.giveGroupRecords(irdPtr,_IFT_Domain_nset,"Sets",DataReader::IR_setRec,/*optional*/true);
        setList.resize(setRecs.size());
        for(InputRecord& ir: setRecs){
            // read type of set
            IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);
            // Only one set for now (i don't see any need to ever introduce any other version)
            std::unique_ptr< Set > set = std::make_unique<Set>(num, this); //classFactory.createSet(name.c_str(), num, this)
            if ( !set ) {
                OOFEM_ERROR("Couldn't create set: %s", name.c_str());
            }

            set->initializeFrom(ir);

            // check number
            if ( num < 1 || num > (int)setRecs.size() ) {
                OOFEM_ERROR("Invalid set number (num=%d)", num);
            }

            if ( !setList[num - 1] ) {
                setList[num - 1] = std :: move(set);
            } else {
                OOFEM_ERROR("Set entry already exist (num=%d)", num);
            }

            ir.finish();
        }
    }
    
    #  ifdef VERBOSE
        VERBOSE_PRINT0("Instantiated elements ", elementList.size());
    #  endif

	// read contact surfaces
	DataReader::GroupRecords contactSurfRecs=dr.giveGroupRecords(irdPtr,_IFT_Domain_ncontactsurf,"ContactSurfaces",DataReader::IR_contactSurfaceRec,/*optional*/true);
        contactSurfaceList.clear();
	contactSurfaceList.resize(contactSurfRecs.size());

	for(InputRecord& ir: contactSurfRecs){
	  // read type of contact surface
	  IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);
	  
	  std::unique_ptr< ContactSurface >surface(classFactory.createContactSurface(name.c_str(), num, this) );
	  if ( !surface ) {
            OOFEM_ERROR( "Couldn't create contact surface: %s", name.c_str() );
	  }
	  
	  surface->initializeFrom(ir);
	  
	  // check number
	  if ( num < 1 || num > (int)contactSurfRecs.size() ) {
            OOFEM_ERROR("Invalid contact surface number (num=%d)", num);
	  }
	  
	  if ( !contactSurfaceList [ num - 1 ] ) {
            contactSurfaceList [ num - 1 ] = std::move(surface);
	  } else {
            OOFEM_ERROR("Contact surface entry already exist (num=%d)", num);
	  }
	  
	  ir.finish();
	}
	


    // read cross sections
    DataReader::GroupRecords csRecs=dr.giveGroupRecords(irdPtr,_IFT_Domain_ncrosssect,"CrossSections",DataReader::IR_crosssectRec,/*optional*/false);
    crossSectionList.clear();
    crossSectionList.resize(csRecs.size());
    for(InputRecord& ir: csRecs){
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        std :: unique_ptr< CrossSection >crossSection( classFactory.createCrossSection(name.c_str(), num, this) );
        if ( !crossSection ) {
            OOFEM_ERROR("Couldn't create crosssection: %s", name.c_str());
        }

        crossSection->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > (int)csRecs.size() ) ) {
            OOFEM_ERROR("Invalid crossSection number (num=%d)", num);
        }

        if ( !crossSectionList[num - 1] ) {
            crossSectionList[num - 1] = std :: move(crossSection);
        } else {
            OOFEM_ERROR("crossSection entry already exist (num=%d)", num);
        }

        ir.finish();
    }

    #  ifdef VERBOSE
        VERBOSE_PRINT0("Instantiated cross sections ", crossSectionList.size())
    #  endif

    // read materials
    DataReader::GroupRecords matRecs=dr.giveGroupRecords(irdPtr,_IFT_Domain_nmat,"Materials",DataReader::IR_matRec,/*optional*/false);
    materialList.clear();
    materialList.resize(matRecs.size());
    for(InputRecord& ir: matRecs){
        // read type of material
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        std :: unique_ptr< Material >mat( classFactory.createMaterial(name.c_str(), num, this) );
        if ( !mat ) {
            OOFEM_ERROR("Couldn't create material: %s", name.c_str());
        }

        mat->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > (int)matRecs.size() ) ) {
            OOFEM_ERROR("Invalid material number (num=%d)", num);
        }

        if ( !materialList[num - 1] ) {
            materialList[num - 1] = std :: move(mat);
        } else {
            OOFEM_ERROR("material entry already exist (num=%d)", num);
        }

        ir.finish();
    }

    #  ifdef VERBOSE
        VERBOSE_PRINT0("Instantiated materials ", materialList.size())
    #  endif

    // read barriers
    DataReader::GroupRecords barrRecs=dr.giveGroupRecords(irdPtr,_IFT_Domain_nbarrier,"NonlocalBarriers",DataReader::IR_nlocBarRec,/*optional*/true);
    nonlocalBarrierList.clear();
    nonlocalBarrierList.resize(barrRecs.size());
    for(InputRecord& ir: barrRecs){
        // read type of load
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        std :: unique_ptr< NonlocalBarrier >barrier( classFactory.createNonlocalBarrier(name.c_str(), num, this) );
        if ( !barrier ) {
            OOFEM_ERROR("Couldn't create barrier: %s", name.c_str());
        }

        barrier->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > (int)barrRecs.size() ) ) {
            OOFEM_ERROR("Invalid barrier number (num=%d)", num);
        }

        if ( !nonlocalBarrierList[num - 1] ) {
            nonlocalBarrierList[num - 1] = std :: move(barrier);
        } else {
            OOFEM_ERROR("barrier entry already exist (num=%d)", num);
        }

        ir.finish();
    }

    #  ifdef VERBOSE
        if ( !nonlocalBarrierList.empty() ) {
            VERBOSE_PRINT0("Instantiated barriers ", nonlocalBarrierList.size());
        }
    #  endif

    // read boundary conditions
    DataReader::GroupRecords bcRecs=dr.giveGroupRecords(irdPtr,_IFT_Domain_nbc,"BoundaryConditions",DataReader::IR_bcRec,/*optional*/false);
    bcList.clear();
    bcList.resize(bcRecs.size());
    for(InputRecord& ir: bcRecs){
        // read type of bc
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        std :: unique_ptr< GeneralBoundaryCondition >bc( classFactory.createBoundaryCondition(name.c_str(), num, this) );
        if ( !bc ) {
            OOFEM_ERROR("Couldn't create boundary condition: %s", name.c_str());
        }

        bc->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > (int)bcRecs.size() ) ) {
            OOFEM_ERROR("Invalid boundary condition number (num=%d)", num);
        }

        if ( !bcList[num - 1] ) {
            bcList[num - 1] = std :: move(bc);
        } else {
            OOFEM_ERROR("boundary condition entry already exist (num=%d)", num);
        }

        ir.finish();
    }

    #  ifdef VERBOSE
        VERBOSE_PRINT0("Instantiated BCs ", bcList.size())
    #  endif

    // read initial conditions
    DataReader::GroupRecords icRecs=dr.giveGroupRecords(irdPtr,_IFT_Domain_nic,"InitialConditions",DataReader::IR_icRec,/*optional*/false);
    icList.clear();
    icList.resize(icRecs.size());
    for(InputRecord& ir: icRecs){
        // read type of load
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        std :: unique_ptr< InitialCondition >ic( new InitialCondition(num, this) );
        if ( !ic ) {
            OOFEM_ERROR("Creation of IC no. %d failed", num);
        }

        ic->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > (int)icRecs.size() ) ) {
            OOFEM_ERROR("Invalid initial condition number (num=%d)", num);
        }

        if ( !icList[num - 1] ) {
            icList[num - 1] = std :: move(ic);
        } else {
            OOFEM_ERROR("initial condition entry already exist (num=%d)", num);
        }

        ir.finish();
    }

    #  ifdef VERBOSE
        VERBOSE_PRINT0("Instantiated ICs ", icList.size())
    #  endif


    // read load time functions
    DataReader::GroupRecords ltfRecs=dr.giveGroupRecords(irdPtr,_IFT_Domain_nfunct,"LoadTimeFunctions",DataReader::IR_funcRec,/*optional*/false);
    functionList.clear();
    functionList.resize(ltfRecs.size());
    for(InputRecord& ir: ltfRecs){
        // read type of func
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        std :: unique_ptr< Function >func( classFactory.createFunction(name.c_str(), num, this) );
        if ( !func ) {
            OOFEM_ERROR("Couldn't create time function: %s", name.c_str());
        }

        func->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > (int)ltfRecs.size() ) ) {
            OOFEM_ERROR("Invalid Function number (num=%d)", num);
        }

        if ( !functionList[num - 1] ) {
            functionList[num - 1] = std :: move(func);
        } else {
            OOFEM_ERROR("Function entry already exist (num=%d)", num);
        }

        ir.finish();
    }

    #  ifdef VERBOSE
        VERBOSE_PRINT0("Instantiated load-time fncts ", functionList.size())
    #  endif

    // read sets
    if ( setList.size() == 0 ) {
        DataReader::GroupRecords setRecs=dr.giveGroupRecords(irdPtr,_IFT_Domain_nset,"Sets",DataReader::IR_setRec,/*optional*/true);
        setList.resize(setRecs.size());
        for(InputRecord& ir: setRecs){
            // read type of set
            IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);
            // Only one set for now (i don't see any need to ever introduce any other version)
            std :: unique_ptr< Set > set = std::make_unique<Set>(num, this); //classFactory.createSet(name.c_str(), num, this)
            if ( !set ) {
                OOFEM_ERROR("Couldn't create set: %s", name.c_str());
            }

            set->initializeFrom(ir);

            // check number
            if ( ( num < 1 ) || ( num > (int)setRecs.size() ) ) {
                OOFEM_ERROR("Invalid set number (num=%d)", num);
            }

            if ( !setList[num - 1] ) {
                setList[num - 1] = std :: move(set);
            } else {
                OOFEM_ERROR("Set entry already exist (num=%d)", num);
            }

            ir.finish();
        }
    }

    #  ifdef VERBOSE
        if ( setList.size() ) {
            VERBOSE_PRINT0("Instantiated sets ", setList.size());
        }
    #  endif

    InputRecord* xfmIr=dr.giveChildRecord(irdPtr,_IFT_Domain_nxfemman,DataReader::InputRecordTags[DataReader::IR_xfemManRec],DataReader::IR_xfemManRec,/*optional*/true);
    if(xfmIr){
        DataReader::RecordGuard scope(dr,xfmIr);
        IR_GIVE_RECORD_KEYWORD_FIELD(*xfmIr, name, num);
        xfemManager = classFactory.createXfemManager(name.c_str(), this);
        if ( !xfemManager ) {
            OOFEM_ERROR("Couldn't create xfemmanager: %s", name.c_str());
        }

        xfemManager->initializeFrom(*xfmIr);
        xfemManager->instanciateYourself(dr);
        #  ifdef VERBOSE
            VERBOSE_PRINT0("Instantiated xfem ", 1);
        #  endif
    }


    this->topology = NULL;
    if ( topologytype.length() > 0 ) {
        this->topology = classFactory.createTopology(topologytype.c_str(), this);
        if ( !this->topology ) {
            OOFEM_ERROR("Couldn't create topology of type '%s'", topologytype.c_str());
        }

        return this->topology->instanciateYourself(dr);
        #  ifdef VERBOSE
            VERBOSE_PRINT0("Instantiated topologies ", topologytype.length());
        #  endif
    }

    InputRecord* fmanIr=dr.giveChildRecord(irdPtr,_IFT_Domain_nfracman,"FractureManager",DataReader::IR_fracManRec,/*optional*/true);
    if (fmanIr) {
        fracManager = std::make_unique<FractureManager>(this);
        fracManager->initializeFrom(*fmanIr);
        fracManager->instanciateYourself(dr);
        #  ifdef VERBOSE
            VERBOSE_PRINT0("Instantiated fracture manager ", 1);
        #  endif
    }

    this->initializeFinish();
    return 1;
}

void
Domain::initializeFinish() {

    ///@todo Eventually remove this backwards compatibility:
    //_HeatTransferMode _HeatMass1Mode // Are these deprecated?
    // set the number of spatial dimensions
    if (this->nsd <=0) {
        if ( dType == _1dTrussMode ) {
            nsd = 1;
        } else if ( dType == _2dIncompressibleFlow || dType == _2dBeamMode || dType == _2dTrussMode || dType == _2dMindlinPlateMode || dType == _PlaneStrainMode || dType == _2dPlaneStressMode || dType == _2dPlaneStressRotMode || dType == _WarpingMode ) {
            nsd = 2;
        } else if ( dType == _3dIncompressibleFlow || dType == _3dShellMode || dType == _3dMode || dType == _3dDirShellMode ) {
            nsd = 3;
        } else if ( dType == _3dAxisymmMode ) {
            nsd = 2;
            axisymm = true;
        } 
    }
    

    spatialLocalizer = std::make_unique<OctreeSpatialLocalizer>(this);
    spatialLocalizer->init();
    connectivityTable = std::make_unique<ConnectivityTable>(this);
    OOFEM_LOG_INFO("Spatial localizer init done\n");

    // finish initialization of components (at present neded for dofManagers)
    for ( auto &dman: dofManagerList ) {
        dman->initializeFinish();
    }

    // change internal component references from labels to assigned local numbers
    MapBasedEntityRenumberingFunctor labelToLocNumFunctor(dofmanGlobal2LocalMap, elementGlobal2LocalMap);
    for ( auto &dman: this->dofManagerList ) {
        dman->updateLocalNumbering(labelToLocNumFunctor);
    }

    for ( auto &element: this->elementList ) {
        element->updateLocalNumbering(labelToLocNumFunctor);
    }

    for ( auto &set: setList ) {
        set->updateLocalNumbering(labelToLocNumFunctor);
    }

    BuildMaterialToElementMap();


    // New  - in development /JB
    // set element cross sections based on element set definition and set the corresponding
    // material based on the cs
    for ( int i = 1; i <= this->giveNumberOfCrossSectionModels(); i++ ) {
        if ( int setNum = this->giveCrossSection(i)->giveSetNumber() ) {
            Set *set = this->giveSet(setNum);
            for ( int ielem: set->giveElementList() ) {
                Element *element = this->giveElement( ielem );
                // hack to set the cross section from the cross section model to the element
                // this is obsolete feature, but supported here for backwards compatibility
                // the supported way is to set cross section using set element properties
                size_t indx = Element::IPK_Element_crosssect.getIndex();
                if (this->elementPPM.getPriority(ielem, indx) <= 1) {
                    this->elementPPM.setPriority(ielem, indx, 1);
                    element->setCrossSection(i);
                }  
            }
        }
    }

    // process sets and apply element and dofman properties
    for ( int i = 1; i <= this->giveNumberOfSets(); i++ ) {
        Set *set = this->giveSet(i);
        std::string elemProps = set->giveElementProperties();
        if (!elemProps.empty()) {
            OOFEMTXTInputRecord ir (-1, elemProps);
            for ( int ielem: set->giveElementList() ) {
                Element *element = this->giveElement( ielem );
                element->initializeFrom(ir, 1); // initialize with priority 1 (lower than component record priority)
            }
        }
        std::string dofmanProps = set->giveDofManProperties();
        if (!dofmanProps.empty()) {
            OOFEMTXTInputRecord ir (-1, dofmanProps);
            for ( int idofman: set->giveNodeList() ) {
                DofManager *dofman = this->giveDofManager( idofman );
                dofman->initializeFrom(ir, 1); // initialize with priority 1 (lower than component record priority)
            }
        }
    } 


    for ( auto &el: elementList ) {
        el->initializeFinish();
    }
}

void
Domain :: postInitialize()
{

    if (!spatialLocalizer) {
        spatialLocalizer= std::make_unique<OctreeSpatialLocalizer>(this);
        spatialLocalizer->init();
    }
    if (!connectivityTable) {
        connectivityTable = std::make_unique<ConnectivityTable>(this);
    }
    
    if ( this->hasXfemManager() ) {
        this->giveXfemManager()->postInitialize();
    }

    // Dofs must be created before dof managers due their post-initialization:
    this->createDofs();

    for ( auto &dman: dofManagerList ) {
        dman->postInitialize();
    }


    for ( auto &el: elementList ) {
        el->postInitialize();
    }
    //@TODO clear parameterPriorityManager

    for ( auto &cs: contactSurfaceList ) {
      cs->postInitialize();
    }

    
    for ( auto &bc: bcList ) {
        bc->postInitialize();
    }

    this->elementPPM.clear();
    this->dofmanPPM.clear();

}


std :: string
Domain :: errorInfo(const char *func) const
{
    return std::string("Domain::") + func + ", number: " + std::to_string(number);
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
        defaultNodeDofIDArry = {D_u, D_v, R_w};
    } else if ( dType == _2dPlaneStressMode ) {
        defaultNodeDofIDArry = {D_u, D_v};
    } else if ( dType == _PlaneStrainMode ) {
        defaultNodeDofIDArry = {D_u, D_v};
    } else if  ( dType == _3dMode ) {
        defaultNodeDofIDArry = {D_u, D_v, D_w};
    } else if ( dType == _3dAxisymmMode ) {
        defaultNodeDofIDArry = {D_u, D_v, R_w};
    } else if  ( dType == _2dMindlinPlateMode ) {
        defaultNodeDofIDArry = {D_w, R_u, R_v};
    } else if ( dType == _3dShellMode ) {
        defaultNodeDofIDArry = {D_u, D_v, D_w, R_u, R_v, R_w};
    } else if  ( dType == _2dTrussMode ) {
        defaultNodeDofIDArry = {D_u, D_w};
    } else if  ( dType == _1dTrussMode ) {
        defaultNodeDofIDArry = {D_u};
    } else if  ( dType == _2dBeamMode ) {
        defaultNodeDofIDArry = {D_u, D_w, R_v};
    } else if  ( dType == _2dLatticeMode ) {
        defaultNodeDofIDArry = {D_u, D_v, R_w};
    } else if  ( dType == _HeatTransferMode ) {
        defaultNodeDofIDArry = {T_f};
    } else if  ( dType == _Mass1TransferMode ) {
        defaultNodeDofIDArry = {C_1};
    } else if  ( dType == _HeatMass1Mode ) {
        defaultNodeDofIDArry = {T_f, C_1};
    }  else if ( dType == _2dIncompressibleFlow ) {
        defaultNodeDofIDArry = {V_u, V_v, P_f};
    }  else if ( dType == _3dIncompressibleFlow ) {
        defaultNodeDofIDArry = {V_u, V_v, V_w, P_f};
    }  else if ( dType == _3dDirShellMode ) {
        defaultNodeDofIDArry = {D_u, D_v, D_w, W_u, W_v, W_w, Gamma};
    }  else if ( dType == _2dLatticeMassTransportMode ) {
        defaultNodeDofIDArry = {P_f};
    }  else if ( dType == _3dLatticeMassTransportMode ) {
        defaultNodeDofIDArry = {P_f};
    }  else if ( dType == _3dLatticeMode ) {
        defaultNodeDofIDArry = {D_u, D_v, D_w, R_u, R_v, R_w};
    }  else if ( dType == _2dLatticeHeatTransferMode ) {
        defaultNodeDofIDArry = {T_f};
    }  else if ( dType == _3dLatticeHeatTransferMode ) {
        defaultNodeDofIDArry = {T_f};
    }  else if ( dType == _WarpingMode ) {
        defaultNodeDofIDArry = {D_w};
    } else {
        OOFEM_ERROR("unknown domainType (%s)", __domainTypeToString(dType));
    }

    return defaultNodeDofIDArry;
}


int
Domain :: giveNumberOfSpatialDimensions()
{
    return nsd;
}


bool
Domain :: isAxisymmetric()
{
    return axisymm;
}


void
Domain :: resolveDomainDofsDefaults(const char *typeName)
//
// resolves default number of dofs per node according to domain type name.
// and also resolves default dof mask according to domain type.
//
{
    if ( !strncmp(typeName, "2dplanestressrot", 16) ) {
        dType = _2dPlaneStressRotMode;
    } else if ( !strncmp(typeName, "2dplanestress", 12) ) {
        dType = _2dPlaneStressMode;
    } else if ( !strncmp(typeName, "planestrain", 11) ) {
        dType = _PlaneStrainMode;
    } else if ( !strncmp(typeName, "3daxisymm", 9) ) {
        dType = _3dAxisymmMode;
    } else if  ( !strncmp(typeName, "2dmindlinplate", 14) ) {
        dType = _2dMindlinPlateMode;
    } else if ( !strncmp(typeName, "3dshell", 7) ) {
        dType = _3dShellMode;
    } else if  ( !strncmp(typeName, "2dtruss", 7) ) {
        dType = _2dTrussMode;
    } else if  ( !strncmp(typeName, "1dtruss", 7) ) {
        dType = _1dTrussMode;
    } else if  ( !strncmp(typeName, "2dbeam", 6) ) {
        dType = _2dBeamMode;
    } else if  ( !strncmp(typeName, "2dlattice", 9) ) {
        dType = _2dLatticeMode;
    } else if  ( !strncmp(typeName, "heattransfer", 12) ) {
        dType = _HeatTransferMode;
    } else if  ( !strncmp(typeName, "mass1transfer", 13) ) {
        dType = _Mass1TransferMode;
    } else if  ( !strncmp(typeName, "hema1", 5) ) {
        dType = _HeatMass1Mode;
    } else if ( !strncmp(typeName, "2dincompflow", 12) ) {
        dType = _2dIncompressibleFlow;
    } else if ( !strncmp(typeName, "3dincompflow", 12) ) {
        dType = _3dIncompressibleFlow;
    } else if  ( !strncmp(typeName, "3ddirshell", 10) ) {
        dType = _3dDirShellMode;
    } else if  ( !strncmp(typeName, "2dmasslatticetransport", 22) ) {
        dType = _2dLatticeMassTransportMode;
    } else if  ( !strncmp(typeName, "3dlattice", 9) ) {
        dType = _3dLatticeMode;
    } else if  ( !strncmp(typeName, "3dmasslatticetransport", 22) ) {
        dType = _3dLatticeMassTransportMode;
    } else if  ( !strncmp(typeName, "2dheatlattice", 13) ) {
        dType = _3dLatticeMassTransportMode;
    } else if  ( !strncmp(typeName, "3dheatlattice", 13) ) {
        dType = _3dLatticeMassTransportMode;
    } else if  ( !strncmp(typeName, "3d", 2) ) {
        dType = _3dMode;
    } else if  ( !strncmp(typeName, "warping", 7) ) {
        dType = _WarpingMode;
    } else {
        OOFEM_ERROR("unknown domainType (%s)", typeName);
        return;
    }
}


NodalRecoveryModel *
Domain :: giveSmoother()
{
    return this->smoother.get();
}


void
Domain :: setSmoother(NodalRecoveryModel *newSmoother, bool destroyOld)
{
    if ( !destroyOld ) {
        this->smoother.release();
    }

    this->smoother.reset(newSmoother);
}


void
Domain :: setTopology(TopologyDescription *topo, bool destroyOld)
{
    if ( !destroyOld ) {
        this->topology.release();
    }

    this->topology.reset(topo);
}


ConnectivityTable *
Domain :: giveConnectivityTable()
//
// return connectivity Table - if no defined - creates new one
//
{
    if ( !connectivityTable ) {
        //connectivityTable = std::make_unique<ConnectivityTable>(this);
        OOFEM_LOG_ERROR("Connectivity table init error");
    }

    return connectivityTable.get();
}


SpatialLocalizer *
Domain :: giveSpatialLocalizer()
//
// return connectivity Table - if no defined - creates new one
//
{
//     if (spatialLocalizer == NULL) spatialLocalizer = new DummySpatialLocalizer(1, this);
    if ( spatialLocalizer ) {
        return spatialLocalizer.get();
    } else {
        OOFEM_LOG_ERROR("Spatial localizer init failure");
        return nullptr;      
    }
}

void Domain :: setSpatialLocalizer(std::unique_ptr<SpatialLocalizer> sl){
    spatialLocalizer = std::move(sl);
}


void
Domain :: createDofs()
{

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////// Step 1. Scan all required nodal dofs.
    std :: vector< std :: set< int > > node_dofs( this->giveNumberOfDofManagers() );
    for ( auto &element: this->elementList ) {
        IntArray dofids;
        // Scan for all dofs needed by element.
        for ( int j = 1; j <= element->giveNumberOfNodes(); ++j ) {
            element->giveDofManDofIDMask(j, dofids);
            for ( int k = 1; k <= dofids.giveSize(); k++ ) {
                node_dofs [ element->giveNode(j)->giveNumber() - 1 ].insert( dofids.at(k) );
            }
        }
    }
    for ( auto &dman: this->dofManagerList ) {
        // Nodes can also contain their own list of dofs (typical usecase: RigidArmNode )
        const IntArray *dofids = dman->giveForcedDofIDs();
        if ( dofids ) {
            for ( int k = 1; k <= dofids->giveSize(); ++k ) {
                node_dofs [ dman->giveNumber() - 1 ].insert( dofids->at(k) );
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Step 2. Scan all Dirichlet b.c.s (or active dofs). For every node we store a map from the dofid to it's b.c. number.
    // This loop won't check for slave dofs or so, and will give a bc id for every single relevant dof.
    // This must be a separate step since we store the inverse mapping (bc->dof instead of dof->bc) so we want to loop over all b.c.s to invert this.
    std :: vector< std :: map< int, int > > dof_bc( this->giveNumberOfDofManagers() );
    for ( int i = 1; i <= this->giveNumberOfBoundaryConditions(); ++i ) {
        GeneralBoundaryCondition *gbc = this->giveBc(i);
        if ( gbc->giveSetNumber() > 0 ) { ///@todo This will eventually not be optional.
            // Loop over nodes in set and store the bc number in each dof.
            Set *set = this->giveSet( gbc->giveSetNumber() );
            ActiveBoundaryCondition *active_bc = dynamic_cast< ActiveBoundaryCondition * >(gbc);
            BoundaryCondition *bc = dynamic_cast< BoundaryCondition * >(gbc);
            if ( bc || ( active_bc && active_bc->requiresActiveDofs() ) ) {
                const IntArray &appliedDofs = gbc->giveDofIDs();
                const IntArray &nodes = set->giveNodeList();
                for ( int inode = 1; inode <= nodes.giveSize(); ++inode ) {
                    for ( int idof = 1; idof <= appliedDofs.giveSize(); ++idof ) {
                        dof_bc [ nodes.at(inode) - 1 ] [ appliedDofs.at(idof) ] = i;
                    }
                }
            }
        }
    }
    // Step 2b. This step asks nodes for their bc-vector, which is the old approach to dirichlet b.c.s (i.e. this is for backwards compatibility)
    ///@todo Remove this input method whenever we decide on deprecating the old approach.
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); ++i ) {
        DofManager *dman = this->giveDofManager(i);
        const std :: map< int, int > *dmanBcs = dman->giveBcMap();
        if ( dmanBcs ) {
            dof_bc [ i - 1 ].insert( dmanBcs->begin(), dmanBcs->end() );     // This will ignore duplicated dofiditems.
        }
    }

    // Step 2c. THis will apply BC set to internal element nodes.
    for ( int i = 1; i <= this->giveNumberOfBoundaryConditions(); ++i ) {
        GeneralBoundaryCondition *gbc = this->giveBc(i);
        if ( gbc->giveSetNumber() > 0 ) { ///@todo This will eventually not be optional.
            // Loop over nodes in set and store the bc number in each dof.
            Set *set = this->giveSet( gbc->giveSetNumber() );
            BoundaryCondition *bc = dynamic_cast< BoundaryCondition * >(gbc);
            if ( bc ) {
                const IntArray &appliedDofs = gbc->giveDofIDs();
                const IntArray &internalelementnodes = set->giveInternalElementDofManagerList();
                for ( int indx = 1; indx <= internalelementnodes.giveSize()/2; ++indx ) {
                    int ielem = internalelementnodes.at(indx*2-1);
                    int idofman = internalelementnodes.at(indx*2);
                    for ( int idof = 1; idof <= appliedDofs.giveSize(); ++idof ) {
                        // test if element internal dofman has the dof
                        DofManager* idman = this->giveElement(ielem)->giveInternalDofManager(idofman);
                        if (idman && idman->hasDofID((DofIDItem) appliedDofs.at(idof))) {
                            // assign bc to existing DOF
                            idman->giveDofWithID((DofIDItem) appliedDofs.at(idof))->setBcId(i);
                        }
                    }
                }
            }
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Step 3. Same for initial conditions as for boundary conditions in step 2.
    std :: vector< std :: map< int, int > > dof_ic( this->giveNumberOfDofManagers() );
    for ( int i = 1; i <= this->giveNumberOfInitialConditions(); ++i ) {
        InitialCondition *ic = this->giveIc(i);
        if ( ic->giveSetNumber() > 0 ) { ///@todo This will eventually not be optional.
            // Loop over nodes in set and store the bc number in each dof.
            Set *set = this->giveSet( ic->giveSetNumber() );
            const IntArray &appliedDofs = ic->giveDofIDs();
            const IntArray &nodes = set->giveNodeList();
            for ( int inode = 1; inode <= nodes.giveSize(); ++inode ) {
                for ( int idof = 1; idof <= appliedDofs.giveSize(); ++idof ) {
                    dof_ic [ nodes.at(inode) - 1 ] [ appliedDofs.at(idof) ] = i;
                }
            }
        }
    }
    // Step 3b. This step asks nodes for their bc-vector, which is the old approach to dirichlet b.c.s (i.e. this is for backwards compatibility)
    ///@todo Remove this input method whenever we decide on deprecating the old approach.
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); ++i ) {
        DofManager *dman = this->giveDofManager(i);
        const std :: map< int, int > *dmanIcs = dman->giveIcMap();
        if ( dmanIcs ) {
            dof_ic [ i - 1 ].insert( dmanIcs->begin(), dmanIcs->end() );     // This will ignore duplicated dofiditems.
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Step 3. Create the dofs. This involves obtaining the correct
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); ++i ) {
        DofManager *dman = this->giveDofManager(i);
        //printf("Dofs in node %d (of %d) = %d\n", i, this->giveNumberOfDofManagers(), node_dofs[i-1].size());

        /* do not delete existing DOFs; that may be created during adaptive solution scheme (mesh generator applies DOFs) */
        if (0) dman->setNumberOfDofs(0);

        for ( int id: node_dofs [ i - 1 ] ) {
            // Find bc and ic if there are any, otherwise zero.
            int bcid = dof_bc [ i - 1 ].find(id) != dof_bc [ i - 1 ].end() ? dof_bc [ i - 1 ] [ id ] : 0;
            //int icid = dof_ic [ i - 1 ].find(id) != dof_ic [ i - 1 ].end() ? dof_ic [ i - 1 ] [ id ] : 0;

            // Determine the doftype:
            dofType dtype = DT_master;
            const std :: map< int, int > *dmanTypes = dman->giveDofTypeMap();
            if ( dmanTypes ) {
                std :: map< int, int > :: const_iterator it = dmanTypes->find(id);
                if ( it != dmanTypes->end() ) {
                    dtype = ( dofType ) it->second;
                }
            }
            // Check if active dofs are needed:
            if ( bcid > 0 ) {
                // What should take precedence here if there is a slave node?
                // Right now the active b.c. overrides anything set prior, if necessary.
                // This seems like the most suitable choice, but it could possibly be changed.
                ActiveBoundaryCondition *active_bc = dynamic_cast< ActiveBoundaryCondition * >( this->giveBc(bcid) );
                if ( active_bc && active_bc->requiresActiveDofs() ) {
                    dtype = DT_active;
                }
            }

            if ( !dman->isDofTypeCompatible(dtype) ) {
                OOFEM_ERROR("Incompatible dof type (%d) in node %d", dtype, i);
            }

            // Finally create the new DOF: 
            //printf("Creating: node %d, id = %d, dofType = %d, bc = %d, ic = %d\n", i, id, dtype, bcid, icid);
            if ( !dman->hasDofID((DofIDItem)id) ) {

                Dof *dof = classFactory.createDof(dtype, (DofIDItem)id, dman);
                // Slave dofs obtain their weights post-initialization, simple slave dofs must have their master node specified.
                if ( dtype == DT_simpleSlave ) {
                    static_cast< SimpleSlaveDof * >(dof)->setMasterDofManagerNum( ( * dman->giveMasterMap() ) [ id ] );
                }
                dman->appendDof(dof);
            }
        }
        // assign bc & ic data on all dofs (existing & newly created ones)
        for ( Dof *dof: *dman ) {
            // Find bc and ic if there are any, otherwise zero.
            DofIDItem id = dof->giveDofID();
            int bcid = dof_bc [ i - 1 ].find(id) != dof_bc [ i - 1 ].end() ? dof_bc [ i - 1 ] [ id ] : 0;
            int icid = dof_ic [ i - 1 ].find(id) != dof_ic [ i - 1 ].end() ? dof_ic [ i - 1 ] [ id ] : 0;
	    // preserve existing BC and IC if set (required for adaptivity, BC and IC set by mesher)
	    if (dof->isPrimaryDof()) {
	      if (dof->giveBcId()==0) dof->setBcId(bcid); // Note: slave dofs and such will simple ignore this.
	      if (dof->giveIcId()==0) dof->setIcId(icid);
	    }
        }
    }

    // XFEM manager create additional dofs themselves:
    if ( this->hasXfemManager() ) {
        xfemManager->createEnrichedDofs();
    }

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
    int result = 1;

    for ( auto &dman: this->dofManagerList ) {
        result &= dman->checkConsistency();
    }

    for ( auto &element: this->elementList ) {
        result &= element->checkConsistency();
    }

    for ( auto &material: this->materialList ) {
        result &= material->checkConsistency();
    }

    return result;
}

double
Domain :: giveArea()
{
    double area = 0.0;
    for ( auto &element: this->elementList ) {
        area += element->computeArea();
    }

    return area;
}

double
Domain :: giveVolume()
{
    double volume = 0.0;
    for ( auto &element: this->elementList ) {
        volume += element->computeVolume();
    }

    return volume;
}

double
Domain :: giveSize()
{
    double volume = 0.0;
    for ( auto &element: this->elementList ) {
        volume += element->computeVolumeAreaOrLength();
    }

    return volume;
}

int
Domain :: giveNextFreeDofID(int increment)
{
    if ( this->engineeringModel->isParallel() ) {
        OOFEM_ERROR("Additional dof id's not implemented/tested for parallel problems");
    }

    int freeID = this->freeDofID;
    this->freeDofID += increment;
    return freeID;
}

void
Domain :: resetFreeDofID()
{
    this->freeDofID = MaxDofID;
}

ErrorEstimator *
Domain :: giveErrorEstimator()
{
    return engineeringModel->giveDomainErrorEstimator(this->number);
}


OutputManager *
Domain :: giveOutputManager()
{
    return outputManager.get();
}


TopologyDescription *
Domain :: giveTopology()
{
    return topology.get();
}


template< typename T >
void save_components(T &list, DataStream &stream, ContextMode mode)
{
    if ( !stream.write((int)list.size()) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    for ( const auto &object: list ) {
        if ( ( mode & CM_Definition ) != 0 ) {
            if ( stream.write( std :: string( object->giveInputRecordName() ) ) == 0 ) {
                THROW_CIOERR(CIO_IOERR);
            }
        }
        object->saveContext(stream, mode);
    }
}


template< typename T, typename C >
void restore_components(T &list, DataStream &stream, ContextMode mode, const C &creator)
{
    int size = 0;
    if ( !stream.read(size) ) {
        THROW_CIOERR(CIO_IOERR);
    }
    if ( mode & CM_Definition ) {
        list.clear();
        list.resize(size);
    }
    for ( int i = 1; i <= size; i++ ) {
        if ( mode & CM_Definition ) {
            std :: string name;
            if ( !stream.read(name) ) {
                THROW_CIOERR(CIO_IOERR);
            }
            list[i-1] = creator(name, i);
        }
        list[i-1]->restoreContext(stream, mode);
    }
}


void
Domain :: saveContext(DataStream &stream, ContextMode mode)
{
    if ( !stream.write(this->giveSerialNumber()) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( mode & CM_Definition ) ) {
        save_components(this->setList, stream, mode);
        save_components(this->materialList, stream, mode);
        save_components(this->crossSectionList, stream, mode);
        save_components(this->icList, stream, mode);
        save_components(this->functionList, stream, mode);
        save_components(this->nonlocalBarrierList, stream, mode);
    }

    save_components(this->dofManagerList, stream, mode);
    save_components(this->elementList, stream, mode);
    save_components(this->bcList, stream, mode);

    auto ee = this->giveErrorEstimator();
    if ( ee ) {
        ee->saveContext(stream, mode);
    }
}


void
Domain :: restoreContext(DataStream &stream, ContextMode mode)
{
    bool domainUpdated = false;
    int serNum = this->giveSerialNumber();
    // restore domain serial number
    if ( !stream.read(this->serialNumber) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( mode & CM_Definition ) ) {
        // clear cached data:
        elementGlobal2LocalMap.clear();
        dofmanGlobal2LocalMap.clear();
 
        ///@todo Saving and restoring xfemmanagers.
        xfemManager = nullptr;

        restore_components(this->setList, stream, mode,
                           [this] (std::string &x, int i) { return std::make_unique<Set>(i, this); });
        restore_components(this->materialList, stream, mode,
                           [this] (std::string &x, int i) { return classFactory.createMaterial(x.c_str(), i, this); });
        restore_components(this->crossSectionList, stream, mode,
                           [this] (std::string &x, int i) { return classFactory.createCrossSection(x.c_str(), i, this); });
        restore_components(this->icList, stream, mode,
                           [this] (std::string &x, int i) { return classFactory.createInitialCondition(x.c_str(), i, this); });
        restore_components(this->functionList, stream, mode,
                           [this] (std::string &x, int i) { return classFactory.createFunction(x.c_str(), i, this); });
        restore_components(this->nonlocalBarrierList, stream, mode,
                           [this] (std::string &x, int i) { return classFactory.createNonlocalBarrier(x.c_str(), i, this); });

        domainUpdated = true;
    } else {
        if ( serNum != this->giveSerialNumber() ) {
            // read corresponding domain
            OOFEM_LOG_INFO("restoring domain %d.%d\n", this->number, this->giveSerialNumber());
            OOFEMTXTDataReader domainDr(this->engineeringModel->giveDomainFileName(1, this->giveSerialNumber()));
            this->clear();

            if ( !this->instanciateYourself(domainDr) ) {
                OOFEM_ERROR("domain Instanciation failed");
            }

            domainUpdated = true;
        }
    }

    restore_components(this->dofManagerList, stream, mode,
                       [this] (std::string &x, int i) { return classFactory.createDofManager(x.c_str(), i, this); });
    restore_components(this->elementList, stream, mode,
                       [this] (std::string &x, int i) { return classFactory.createElement(x.c_str(), i, this); });
    restore_components(this->bcList, stream, mode,
                       [this] (std::string &x, int i) { return classFactory.createBoundaryCondition(x.c_str(), i, this); });

    auto ee = this->giveErrorEstimator();
    if ( ee ) {
        if ( domainUpdated ) {
            ee->setDomain(this);
        }
        ee->restoreContext(stream, mode);
    }

    if ( domainUpdated ) {
        if ( this->smoother ) {
            this->smoother->clear();
        }
    }
}

#ifdef __MPI_PARALLEL_MODE

DomainTransactionManager *
Domain :: giveTransactionManager()
{
    if ( !transactionManager ) {
        transactionManager = std::make_unique<DomainTransactionManager>(this);
        if ( !transactionManager ) {
            OOFEM_ERROR("allocation failed");
        }
    }

    return transactionManager.get();
}


int Domain :: commitTransactions(DomainTransactionManager *tm)
{
    if ( tm->dofmanTransactions.empty() && tm->elementTransactions.empty() ) {
        return 1;
    }

    this->initGlobalDofManMap();
    for ( auto &dmanTrans: tm->dofmanTransactions ) {
        if ( dmanMap.find(dmanTrans.first) != dmanMap.end() ) {
            int lnum = dmanMap [ dmanTrans.first ]->giveNumber();
            DofManager *dman = dofManagerList[lnum-1].release();
            dmanMap.erase(dmanTrans.first);
            delete dman;
        }

        if ( dmanTrans.second ) {
            dmanMap [ dmanTrans.first ] = ( DofManager * ) dmanTrans.second;
        }
    }

    this->initGlobalElementMap();
    for ( auto elTrans: tm->elementTransactions ) {
        int gen = elTrans.first;
        if ( elementMap.find(gen) != elementMap.end() ) {
            int lnum = elementMap [ gen ]->giveNumber();
            Element *elem = elementList[lnum-1].release();
            elementMap.erase(gen);
            delete elem;
        }

        if ( elTrans.second ) {
            elementMap [ gen ] = ( Element * ) elTrans.second;
        }
    }

    /*
     * if (tm->transactions.empty()) return 1;
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
     *    dman = dofManagerList.release (t._num-1);
     *    delete dman;
     *  } else if (t._ctype == DCT_Element) {
     *    elem = elementList.release (t._num-1);
     *    delete elem;
     *  } else {
     *    OOFEM_ERROR("unknown transaction component type");
     *  }
     * } else if (t._ttype == DTT_ADD) {
     *  if (t._ctype == DCT_DofManager) {
     *    dmanMap[t._num] = (*DofManager) t._obj;
     *  } else if (t._ctype == DCT_Element) {
     *    elemMap[t._num] = (*Element) t._obj;
     *  } else {
     *    OOFEM_ERROR("unknown transaction component type");
     *  }
     * } else {
     *  OOFEM_ERROR("unknown transaction type");
     * }
     *
     * // Pop new transaction
     * tm->transactions.pop_front();
     * } // while (!tm->transactions.empty()) {
     *
     */

    this->renumberDofManagers();
    this->renumberDofManData(tm);

    this->renumberElements();
    this->renumberElementData(tm);

    ///@todo This is really quite bad. We shouldn't do unique ptrs like this (we have the same ptr in 2 different unique_ptrs!!!)
    // initialize new element list
    std :: vector< std :: unique_ptr< Element > > elementList_new;
    elementList_new.reserve(elementMap.size());
    for ( auto &map: elementMap ) {
        elementList_new.emplace_back(map.second);
    }
    for ( auto &el: this->elementList ) {
        el.release();
    }
    this->elementList = std :: move(elementList_new);

    // initialize new dofman list
    std :: vector< std :: unique_ptr< DofManager > > dofManagerList_new;
    dofManagerList_new.reserve(dmanMap.size());
    for ( auto &map: dmanMap ) {
        dofManagerList_new.emplace_back(map.second);
    }
    for ( auto &dman: this->dofManagerList ) {
        dman.release();
    }
    this->dofManagerList = std :: move(dofManagerList_new);

    BuildElementPlaceInArrayMap();
    BuildDofManPlaceInArrayMap();

    tm->dofmanTransactions.clear();
    tm->elementTransactions.clear();

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
        dmanMap.clear();

        for ( auto &dman: dofManagerList ) {
            int key = dman->giveGlobalNumber();
            dmanMap[key] = dman.get();
        }
    }
}


void
Domain :: initGlobalElementMap(bool forceinit)
{
    if ( forceinit || !elementMapInitialized ) {
        elementMap.clear();

        for ( auto &elem: elementList ) {
            int key = elem->giveGlobalNumber();
            elementMap[key] = elem.get();
        }
    }
}


void
Domain :: renumberDofManData(DomainTransactionManager *tm)
{
    SpecificEntityRenumberingFunctor< Domain > domainGToLFunctor(this, &Domain :: LB_giveUpdatedGlobalNumber);
    SpecificEntityRenumberingFunctor< Domain > domainLToLFunctor(this, &Domain :: LB_giveUpdatedLocalNumber);


    for ( auto &map : dmanMap ) {
        if ( tm->dofmanTransactions.find(map.first) != tm->dofmanTransactions.end() ) {
            // received dof manager -> we map global numbers to new local number
            map.second->updateLocalNumbering(domainGToLFunctor); // g_to_l
        } else {
            // existing dof manager -> we map old local number to new local number
            map.second->updateLocalNumbering(domainLToLFunctor); // l_to_l
        }
    }
}


void
Domain :: renumberElementData(DomainTransactionManager *tm)
{
    SpecificEntityRenumberingFunctor< Domain > domainGToLFunctor(this, &Domain :: LB_giveUpdatedGlobalNumber);
    SpecificEntityRenumberingFunctor< Domain > domainLToLFunctor(this, &Domain :: LB_giveUpdatedLocalNumber);

    for ( auto &map : elementMap ) {
        if ( tm->elementTransactions.find(map.first) != tm->elementTransactions.end() ) {
            // received dof manager -> we map global numbers to new local number
            map.second->updateLocalNumbering(domainGToLFunctor); // g_to_l
        } else {
            // existing dof manager -> we map old local number to new local number
            map.second->updateLocalNumbering(domainLToLFunctor); // l_to_l
        }
    }
}


void
Domain :: renumberDofManagers()
{
    int _locnum = 1;
    for ( auto &map: dmanMap ) {
        map.second->setNumber(_locnum++);
    }
}


void
Domain :: renumberElements()
{
    int _locnum = 1;
    for ( auto &map: elementMap ) {
        map.second->setNumber(_locnum++);
    }
}


int
Domain :: LB_giveUpdatedLocalNumber(int num, EntityRenumberingScheme scheme)
{
    if ( scheme == ERS_DofManager ) {
        auto dm = this->giveDofManager(num);
        if ( dm ) {
            return dm->giveNumber();
        } else {
            OOFEM_ERROR("dofman %d moved to remote partition, updated number not available", num);
        }
    } else {
        OOFEM_ERROR("unsuported renumbering scheme");
    }

    return 0;
}


int
Domain :: LB_giveUpdatedGlobalNumber(int num, EntityRenumberingScheme scheme)
{
    if ( scheme == ERS_DofManager ) {
        auto dm = dmanMap [ num ];
        if ( dm ) {
            return dm->giveNumber();
        } else {
            OOFEM_ERROR("dofman [%d] not available on local partition, updated number not available", num);
            return 0;
        }
    } else {
        OOFEM_ERROR("unsuported renumbering scheme");
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

#endif

void Domain :: BuildElementPlaceInArrayMap()
{
    elementGlobal2LocalMap.clear();

    int nelem = giveNumberOfElements();

    for ( int i = 1; i <= nelem; i++ ) {
        Element *elem = this->giveElement(i);
        int globnum = elem->giveGlobalNumber();
        if ( elementGlobal2LocalMap.find(globnum) == elementGlobal2LocalMap.end() ) {
            // label does not exist yet
            elementGlobal2LocalMap [ globnum ] = i;
        } else {
            OOFEM_ERROR("Element with same global number already exist (label=%d)", globnum);
        }
    }
}

void Domain :: BuildDofManPlaceInArrayMap()
{
    dofmanGlobal2LocalMap.clear();

    int ndman = giveNumberOfDofManagers();

    for ( int i = 1; i <= ndman; i++ ) {
        DofManager *dMan = this->giveDofManager(i);
        int globnum = dMan->giveGlobalNumber();
        if ( dofmanGlobal2LocalMap.find(globnum) == dofmanGlobal2LocalMap.end() ) {
            // label does not exist yet
            dofmanGlobal2LocalMap [ globnum ] = i;
        } else {
            OOFEM_ERROR("DofManagaer with same global number already exist (label=%d)", globnum);
        }
    }
}

void Domain :: BuildMaterialToElementMap()
{
    materialNum2ElMap.clear();

    int nelem = giveNumberOfElements();

    for ( int i = 1; i <= nelem; i++ ) {
        Element *elem = this->giveElement(i);
        int matNum = elem->giveMaterialNumber();
        materialNum2ElMap[ matNum ].followedBy(i);
    }
}

} // end namespace oofem
