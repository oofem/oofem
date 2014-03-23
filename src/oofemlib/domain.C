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

#include "domain.h"
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
#include "datareader.h"
#include "nodalrecoverymodel.h"
#include "nonlocalbarrier.h"
#include "classfactory.h"
#include "logger.h"
#include "xfem/xfemmanager.h"
#include "topologydescription.h"
#include "randomfieldgenerator.h"
#include "errorestimator.h"
#include "range.h"
#include "fracturemanager.h"
#include "dynamicinputrecord.h"
#include "dynamicdatareader.h"
#include "datareader.h"
#include "initmodulemanager.h"
#include "exportmodulemanager.h"
#include "xfem/enrichmentitem.h"
#include "xfem/enrichmentfunction.h"
#include "xfem/enrichmentdomain.h"
#include "xfem/propagationlaw.h"

#include "boundarycondition.h"
#include "activebc.h"
#include "simpleslavedof.h"
#include "masterdof.h"

#ifdef __PARALLEL_MODE
 #include "parallel.h"
 #include "processcomm.h"
 #include "datastream.h"
 #include "communicator.h"
 #include "domaintransactionmanager.h"
#endif

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

#include <cstdarg>
#include <cstring>
#include <vector>
#include <set>

namespace oofem {
Domain :: Domain(int n, int serNum, EngngModel *e) : defaultNodeDofIDArry()
    // Constructor. Creates a new domain.
{
    this->engineeringModel = e;
    this->number = n;
    this->serialNumber = serNum;

    elementList              = new AList< Element >(0);
    dofManagerList           = new AList< DofManager >(0);
    materialList             = new AList< Material >(0);
    bcList                   = new AList< GeneralBoundaryCondition >(0);
    icList                   = new AList< InitialCondition >(0);
    functionList     = new AList< Function >(0);
    crossSectionList         = new AList< CrossSection >(0);
    nonlocalBarierList       = new AList< NonlocalBarrier >(0);
    setList                  = new AList< Set >(0);
    randomFieldGeneratorList = new AList< RandomFieldGenerator >(0);

    dType                 = _unknownMode;

    xfemManager           = NULL;
    connectivityTable     = NULL;
    spatialLocalizer      = NULL;
    outputManager         = new OutputManager(this);
    smoother              = NULL;
    topology              = NULL;
    fracManager           = NULL;

    nonlocalUpdateStateCounter = 0;

    nsd = 0;
    axisymm = false;
    freeDofID = MaxDofID;

#ifdef __PARALLEL_MODE
    dmanMapInitialized = elementMapInitialized = false;
    transactionManager = NULL;
#endif
}


Domain *Domain :: Clone()
{
    /////////////////////////////////////////////////////////
    // Create a copy of the domain using
    // the dynamic data reader.

    EngngModel *eModel = this->giveEngngModel();

    int domNum = this->giveNumber();
    int serNum = this->giveSerialNumber();
    Domain *dNew = new Domain(domNum, serNum, eModel);


    DynamicDataReader dataReader;
    DynamicInputRecord *inputRec;


    //Domain
    inputRec = new DynamicInputRecord();
    inputRec->setField(mDomainType, _IFT_Domain_type);
    dataReader.insertInputRecord(DataReader :: IR_domainRec, inputRec);


    //Output
    inputRec = new DynamicInputRecord();
    inputRec->setRecordKeywordField(_IFT_OutputManager_Name, 1);
    inputRec->setField(_IFT_OutputManager_tstepall);
    inputRec->setField(_IFT_OutputManager_dofmanall);
    inputRec->setField(_IFT_OutputManager_elementall);
    dataReader.insertInputRecord(DataReader :: IR_outManRec, inputRec);

    //Components size record
    inputRec = new DynamicInputRecord();
    inputRec->setField(this->giveNumberOfDofManagers(),                 _IFT_Domain_ndofman);
    inputRec->setField(this->giveNumberOfElements(),                    _IFT_Domain_nelem);
    inputRec->setField(this->giveNumberOfCrossSectionModels(),  _IFT_Domain_ncrosssect);
    inputRec->setField(this->giveNumberOfMaterialModels(),              _IFT_Domain_nmat);
    inputRec->setField(this->giveNumberOfBoundaryConditions(),  _IFT_Domain_nbc);
    inputRec->setField(this->giveNumberOfInitialConditions(),   _IFT_Domain_nic);
    inputRec->setField(this->giveNumberOfFunctions(),   _IFT_Domain_nfunct);
    inputRec->setField(this->giveNumberOfSets(),                _IFT_Domain_nset);
    inputRec->setField(this->giveNumberOfSpatialDimensions(),   _IFT_Domain_numberOfSpatialDimensions);
    if ( this->isAxisymmetric() ) {
        inputRec->setField(_IFT_Domain_axisymmetric);
    }


    // fields to add:
    // inputRec->setField( , _IFT_Domain_nbarrier);
    // inputRec->setField( , _IFT_Domain_nrandgen);
    // inputRec->setField( , _IFT_Domain_topology);
    // inputRec->setField( , _IFT_Domain_nfracman);


    bool nxfemMan = 0;
    if ( this->hasXfemManager() ) {
        nxfemMan = 1;
    }
    inputRec->setField(nxfemMan, _IFT_Domain_nxfemman);


    dataReader.insertInputRecord(DataReader :: IR_domainCompRec, inputRec);


    //Nodes
    int nDofMan = this->giveNumberOfDofManagers();
    for ( int i = 1; i <= nDofMan; i++ ) {
        DynamicInputRecord *nodeRec = new DynamicInputRecord( *this->giveDofManager ( i ) );
        dataReader.insertInputRecord(DataReader :: IR_dofmanRec, nodeRec);
    }

    //Elements
    int nEl = this->giveNumberOfElements();
    for ( int i = 1; i <= nEl; i++ ) {
        DynamicInputRecord *elRec = new DynamicInputRecord( *this->giveElement ( i ) );
        dataReader.insertInputRecord(DataReader :: IR_elemRec, elRec);
    }


    //CrossSection
    int nCS = this->giveNumberOfCrossSectionModels();
    for ( int i = 1; i <= nCS; i++ ) {
        DynamicInputRecord *csRec = new DynamicInputRecord( *this->giveCrossSection ( i ) );
        dataReader.insertInputRecord(DataReader :: IR_crosssectRec, csRec);
    }


    //Material
    int nMat = this->giveNumberOfMaterialModels();
    for ( int i = 1; i <= nMat; i++ ) {
        DynamicInputRecord *matRec = new DynamicInputRecord( *this->giveMaterial ( i ) );
        dataReader.insertInputRecord(DataReader :: IR_matRec, matRec);
    }

    //Boundary Conditions
    int nBC = this->giveNumberOfBoundaryConditions();
    for ( int i = 1; i <= nBC; i++ ) {
        DynamicInputRecord *bcRec = new DynamicInputRecord( *this->giveBc ( i ) );
        dataReader.insertInputRecord(DataReader :: IR_bcRec, bcRec);
    }

    //Initial Conditions
    int nIC = this->giveNumberOfInitialConditions();
    for ( int i = 1; i <= nIC; i++ ) {
        DynamicInputRecord *icRec = new DynamicInputRecord( *this->giveIc ( i ) );
        dataReader.insertInputRecord(DataReader :: IR_icRec, icRec);
    }

    //Load-time functions
    int nLoads = this->giveNumberOfFunctions();
    for ( int i = 1; i <= nLoads; i++ ) {
        DynamicInputRecord *funcRec = new DynamicInputRecord( *this->giveFunction ( i ) );
        dataReader.insertInputRecord(DataReader :: IR_funcRec, funcRec);
    }


    //Sets
    int nSets = this->giveNumberOfSets();
    for ( int i = 1; i <= nSets; i++ ) {
        DynamicInputRecord *setRec = new DynamicInputRecord( *this->giveSet ( i ) );
        dataReader.insertInputRecord(DataReader :: IR_setRec, setRec);
    }

    //XFEM manager
    ///@todo Redesign this part (as well as this whole clone function); / Mikael
    if ( this->xfemManager != NULL ) {
        DynamicInputRecord *xmanRec = new DynamicInputRecord();
        xfemManager->giveInputRecord(* xmanRec);
        dataReader.insertInputRecord(DataReader :: IR_xfemManRec, xmanRec);


        // Enrichment items
        int nEI = xfemManager->giveNumberOfEnrichmentItems();
        for ( int i = 1; i <= nEI; i++ ) {
            EnrichmentItem *ei = xfemManager->giveEnrichmentItem(i);
            ei->appendInputRecords(dataReader);
        }
    }

    dNew->instanciateYourself(& dataReader);
    dNew->postInitialize();

    return dNew;
}

Domain :: ~Domain()
// Destructor.
{
    delete elementList;
    delete dofManagerList;
    delete materialList;
    delete bcList;
    delete icList;
    delete functionList;
    delete crossSectionList;
    delete nonlocalBarierList;
    delete setList;
    delete randomFieldGeneratorList;
    delete xfemManager;
    delete connectivityTable;
    delete spatialLocalizer;
    delete outputManager;
    delete smoother;
    delete topology;

#ifdef __PARALLEL_MODE
    delete transactionManager;
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
    functionList->clear();
    crossSectionList->clear();
    nonlocalBarierList->clear();
    setList->clear();
    randomFieldGeneratorList->clear();
    delete xfemManager;
    xfemManager = NULL;

    if ( connectivityTable ) {
        connectivityTable->reset();
    }

    delete spatialLocalizer;
    spatialLocalizer = NULL;

    if ( smoother ) {
        smoother->clear();
    }

    // bp: how to clear/reset topology data?
    delete topology;
    topology = NULL;

#ifdef __PARALLEL_MODE
    delete transactionManager;
    transactionManager = NULL;
#endif
}


Element *
Domain :: giveElement(int n)
// Returns the n-th element. Generates error if it is not defined yet.
{
#ifdef DEBUG
    if ( !elementList->includes(n) ) {
        OOFEM_ERROR("undefined element (%d)", n);
    }
#endif
    return elementList->at(n);
}

Element *
Domain :: giveGlobalElement(int n)
// Returns the global element with id n. Generates error if it is not defined yet.
{
    for ( int i = 1; i <= elementList->giveSize(); i++ ) {
        if ( elementList->at(i)->giveGlobalNumber() == n ) {
            return elementList->at(i);
        }
    }

    OOFEM_ERROR("undefined element id (%d)", n);
    return NULL;
}

Load *
Domain :: giveLoad(int n)
// Returns the n-th load. Generates the error if not defined.
{
#ifdef DEBUG
    if ( !bcList->includes(n) ) {
        OOFEM_ERROR("undefined load (%d)", n);
    }
    Load *answer = dynamic_cast< Load * >( bcList->at(n) );
    if ( answer ) {
        return answer;
    } else {
        OOFEM_ERROR("cannot cast boundary condition %d to Load class", n);
        return NULL;
    }
#else
    return static_cast< Load * >( bcList->at(n) );

#endif
}


GeneralBoundaryCondition *
Domain :: giveBc(int n)
// Returns the n-th bc. Generates the error if not defined.
{
#ifdef DEBUG
    if ( !bcList->includes(n) ) {
        OOFEM_ERROR("undefined bc (%d)", n);
    }
#endif
    return bcList->at(n);
}


InitialCondition *
Domain :: giveIc(int n)
// Returns the n-th ic. Generates the error if not defined.
{
#ifdef DEBUG
    if ( !icList->includes(n) ) {
        OOFEM_ERROR("undefined ic (%d)", n);
    }
#endif

    return icList->at(n);
}


Function *
Domain :: giveFunction(int n)
// Returns the n-th load-time function. Creates this fuction if it does
// not exist yet.
{
#ifdef DEBUG
    if ( !functionList->includes(n) ) {
        OOFEM_ERROR("undefined load-time function (%d)", n);
    }
#endif

    return functionList->at(n);
}


Material *
Domain :: giveMaterial(int n)
// Returns the n-th material. Creates this material if it does not exist
// yet.
{
#ifdef DEBUG
    if ( !materialList->includes(n) ) {
        OOFEM_ERROR("undefined material (%d)", n);
    }
#endif

    return materialList->at(n);
}


Node *
Domain :: giveNode(int n)
// Returns the n-th node if it exists.
{
#ifdef DEBUG
    if ( !dofManagerList->includes(n) ) {
        OOFEM_ERROR("undefined dofManager (%d)", n);
    }

    Node *node = dynamic_cast< Node * >( dofManagerList->at(n) );
    if ( node == NULL ) {
        OOFEM_ERROR("incompatible type of dofManager %d, can not convert", n);
    }

    return node;

#else
    return static_cast< Node * >( dofManagerList->at(n) );

#endif
}


ElementSide *
Domain :: giveSide(int n)
// Returns the n-th element side.
{
#ifdef DEBUG
    if ( !dofManagerList->includes(n) ) {
        OOFEM_ERROR("undefined dofManager (%d)", n);
    }

    ElementSide *side = dynamic_cast< ElementSide * >( dofManagerList->at(n) );
    if ( !side ) {
        OOFEM_ERROR("incompatible type of dofManager %d, can not convert", n);
    }
    return side;

#else
    return static_cast< ElementSide * >( dofManagerList->at(n) );

#endif
}


DofManager *
Domain :: giveDofManager(int n)
// Returns the n-th node. Creates this node if it does not exist yet.
{
#ifdef DEBUG
    if ( !dofManagerList->includes(n) ) {
        OOFEM_ERROR("undefined dofManager (%d)", n);
    }
#endif
    return dofManagerList->at(n);
}


CrossSection *
Domain :: giveCrossSection(int n)
// Returns the n-th cross section.
// yet.
{
#ifdef DEBUG
    if ( !crossSectionList->includes(n) ) {
        OOFEM_ERROR("undefined cross section (%d)", n);
    }
#endif
    return crossSectionList->at(n);
}


NonlocalBarrier *
Domain :: giveNonlocalBarrier(int n)
// Returns the n-th NonlocalBarrier.
{
#ifdef DEBUG
    if ( !nonlocalBarierList->includes(n) ) {
        OOFEM_ERROR("undefined barrier (%d)", n);
    }
#endif
    return nonlocalBarierList->at(n);
}


RandomFieldGenerator *
Domain :: giveRandomFieldGenerator(int n)
// Returns the n-th RandomFieldGenerator.
{
#ifdef DEBUG
    if ( !randomFieldGeneratorList->includes(n) ) {
        OOFEM_ERROR("undefined generator (%d)", n);
    }
#endif

    return randomFieldGeneratorList->at(n);
}


Set *
Domain :: giveSet(int n)
{
#ifdef DEBUG
    if ( !setList->includes(n) ) {
        OOFEM_ERROR("undefined set (%d)", n);
    }
#endif

    return setList->at(n);
}

XfemManager *
Domain :: giveXfemManager()
{
#ifdef DEBUG
    if ( !xfemManager ) {
        OOFEM_ERROR("undefined xfem manager");
    }
#endif
    return xfemManager;
}

bool
Domain :: hasXfemManager()
{
    return xfemManager != NULL;
}


bool
Domain :: hasFractureManager()
{
    return fracManager != NULL;
}

FractureManager *
Domain :: giveFractureManager()
{
#ifdef DEBUG
    if ( !fracManager ) {
        OOFEM_ERROR("undefined fracture manager");
    }
#endif
    return fracManager;
}

EngngModel *
Domain :: giveEngngModel()
// Returns the time integration algorithm. Creates it if it does not
// exist yet.
{
#ifdef DEBUG
    if ( !engineeringModel ) {
        OOFEM_ERROR("Not defined");
    }
#endif

    return engineeringModel;
}

void Domain :: resizeDofManagers(int _newSize) { dofManagerList->growTo(_newSize); }
void Domain :: resizeElements(int _newSize) { elementList->growTo(_newSize); }
void Domain :: resizeCrossSectionModels(int _newSize) { crossSectionList->growTo(_newSize); }
void Domain :: resizeMaterials(int _newSize) { materialList->growTo(_newSize); }
void Domain :: resizeNonlocalBarriers(int _newSize) { nonlocalBarierList->growTo(_newSize); }
void Domain :: resizeBoundaryConditions(int _newSize) { bcList->growTo(_newSize); }
void Domain :: resizeInitialConditions(int _newSize) { icList->growTo(_newSize); }
void Domain :: resizeFunctions(int _newSize) { functionList->growTo(_newSize); }
void Domain :: resizeRandomFieldGenerators(int _newSize) { randomFieldGeneratorList->growTo(_newSize); }
void Domain :: resizeSets(int _newSize) { setList->growTo(_newSize); }

void Domain :: setDofManager(int i, DofManager *obj) { dofManagerList->put(i, obj); }
void Domain :: setElement(int i, Element *obj) { elementList->put(i, obj); }
void Domain :: setCrossSection(int i, CrossSection *obj) { crossSectionList->put(i, obj); }
void Domain :: setMaterial(int i, Material *obj) { materialList->put(i, obj); }
void Domain :: setNonlocalBarrier(int i, NonlocalBarrier *obj) { nonlocalBarierList->put(i, obj); }
void Domain :: setBoundaryCondition(int i, GeneralBoundaryCondition *obj) { bcList->put(i, obj); }
void Domain :: setInitialCondition(int i, InitialCondition *obj) { icList->put(i, obj); }
void Domain :: setFunction(int i, Function *obj) { functionList->put(i, obj); }
void Domain :: setRandomFieldGenerator(int i, RandomFieldGenerator *obj) { randomFieldGeneratorList->put(i, obj); }
void Domain :: setSet(int i, Set *obj) { setList->put(i, obj); }

void Domain :: clearBoundaryConditions() { bcList->clear(true); }

int
Domain :: instanciateYourself(DataReader *dr)
// Creates all objects mentioned in the data file.
{
    IRResultType result;                            // Required by IR_GIVE_FIELD macro

    int num;
    std :: string name, topologytype;
    int nnode, nelem, nmat, nload, nic, nloadtimefunc, ncrossSections, nbarrier, nrfg, nset = 0;
    bool nxfemman = false;
    bool nfracman = false;
    RandomFieldGenerator *rfg;
    //XfemManager *xMan;
    // mapping from label to local numbers for dofmans and elements
    std :: map< int, int >dofManLabelMap, elemLabelMap;

    FILE *outputStream = this->giveEngngModel()->giveOutputStream();

    // read type of Domain to be solved
    InputRecord *ir = dr->giveInputRecord(DataReader :: IR_domainRec, 1);
    IR_GIVE_FIELD(ir, name, _IFT_Domain_type); // This is inconsistent, "domain" isn't  exactly a field, but the actual record keyword.

    mDomainType = name;

    ir->finish();

#  ifdef VERBOSE
    VERBOSE_PRINT0("Instanciating domain ", this->number);
#  endif

    resolveDomainDofsDefaults( name.c_str() );
    fprintf( outputStream, "Domain type: %s, default ndofs per node is %d\n\n\n",
            name.c_str(), giveDefaultNodeDofIDArry().giveSize() );

    // read output manager record
    std :: string tmp;
    ir = dr->giveInputRecord(DataReader :: IR_outManRec, 1);
    ir->giveRecordKeywordField(tmp);
    outputManager->initializeFrom(ir);
    ir->finish();

    // read domain description
    ir = dr->giveInputRecord(DataReader :: IR_domainCompRec, 1);
    IR_GIVE_FIELD(ir, nnode, _IFT_Domain_ndofman);
    IR_GIVE_FIELD(ir, nelem, _IFT_Domain_nelem);
    IR_GIVE_FIELD(ir, ncrossSections, _IFT_Domain_ncrosssect);
    IR_GIVE_FIELD(ir, nmat, _IFT_Domain_nmat);
    IR_GIVE_FIELD(ir, nload, _IFT_Domain_nbc);
    IR_GIVE_FIELD(ir, nic, _IFT_Domain_nic);
    IR_GIVE_FIELD(ir, nloadtimefunc, _IFT_Domain_nfunct);
    IR_GIVE_OPTIONAL_FIELD(ir, nset, _IFT_Domain_nset);
    IR_GIVE_OPTIONAL_FIELD(ir, nxfemman, _IFT_Domain_nxfemman);
    IR_GIVE_OPTIONAL_FIELD(ir, topologytype, _IFT_Domain_topology);
    this->nsd = -1; ///@todo Change this to default 0 when the domaintype record has been removed.
    IR_GIVE_OPTIONAL_FIELD(ir, this->nsd, _IFT_Domain_numberOfSpatialDimensions);
    this->axisymm = ir->hasField(_IFT_Domain_axisymmetric);
    IR_GIVE_OPTIONAL_FIELD(ir, nfracman, _IFT_Domain_nfracman);

    ///@todo Eventually remove this backwards compatibility:
    //_HeatTransferMode _HeatMass1Mode // Are these deprecated?
    if ( dType == _1dTrussMode ) {
        nsd = 1;
    } else if ( dType == _2dIncompressibleFlow || dType == _2dBeamMode || dType == _2dTrussMode || dType == _2dMindlinPlateMode || dType == _PlaneStrainMode || dType == _2dPlaneStressMode || dType == _2dPlaneStressRotMode ) {
        nsd = 2;
    } else if ( dType == _3dIncompressibleFlow || dType == _3dShellMode || dType == _3dMode || dType == _3dDirShellMode ) {
        nsd = 3;
    } else if ( dType == _3dAxisymmMode ) {
        nsd = 2;
        axisymm = true;
    }


    // read optional number of nonlocalBarriers
    nbarrier = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, nbarrier,  _IFT_Domain_nbarrier);
    // read optional number of RandomFieldGenerator
    nrfg = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, nrfg, _IFT_Domain_nrandgen);



    // read nodes
    dofManagerList->growTo(nnode);
    for ( int i = 1; i <= nnode; i++ ) {
        DofManager *node;
        ir = dr->giveInputRecord(DataReader :: IR_dofmanRec, i);
        // read type of dofManager
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        // assign component number according to record order
        // component number (as given in input record) becomes label
        if ( ( node = classFactory.createDofManager(name.c_str(), i, this) ) == NULL ) {
            OOFEM_ERROR("Couldn't create node of type: %s\n", name.c_str());
        }

        node->initializeFrom(ir);
        if ( dofManLabelMap.find(num) == dofManLabelMap.end() ) {
            // label does not exist yet
            dofManLabelMap [ num ] = i;
        } else {
            OOFEM_ERROR("iDofmanager entry already exist (label=%d)", num);
        }

        node->setGlobalNumber(num);    // set label
        dofManagerList->put(i, node);

        ir->finish();
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Instanciated nodes & sides ", nnode)
#  endif

    // read elements
    elementList->growTo(nelem);
    for ( int i = 1; i <= nelem; i++ ) {
        Element *elem;
        ir = dr->giveInputRecord(DataReader :: IR_elemRec, i);
        // read type of element
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        if ( ( elem = classFactory.createElement(name.c_str(), i, this) ) == NULL ) {
            OOFEM_ERROR("Couldn't create element: %s", name.c_str());
        }

        elem->initializeFrom(ir);

        if ( elemLabelMap.find(num) == elemLabelMap.end() ) {
            // label does not exist yet
            elemLabelMap [ num ] = i;
        } else {
            OOFEM_ERROR("Element entry already exist (label=%d)", num);
        }

        elem->setGlobalNumber(num);
        elementList->put(i, elem);

        ir->finish();
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Instanciated elements ", nelem);
#  endif

    // read cross sections
    crossSectionList->growTo(ncrossSections);
    for ( int i = 1; i <= ncrossSections; i++ ) {
        CrossSection *crossSection;
        ir = dr->giveInputRecord(DataReader :: IR_crosssectRec, i);
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        if ( ( crossSection = classFactory.createCrossSection(name.c_str(), num, this) ) == NULL ) {
            OOFEM_ERROR("Couldn't create crosssection: %s", name.c_str());
        }

        crossSection->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > ncrossSections ) ) {
            OOFEM_ERROR("Invalid crossSection number (num=%d)", num);
        }

        if ( !crossSectionList->includes(num) ) {
            crossSectionList->put(num, crossSection);
        } else {
            OOFEM_ERROR("crossSection entry already exist (num=%d)", num);
        }

        ir->finish();
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Instanciated cross sections ", ncrossSections)
#  endif

    // read materials
    materialList->growTo(nmat);
    for ( int i = 1; i <= nmat; i++ ) {
        Material *mat;
        ir = dr->giveInputRecord(DataReader :: IR_matRec, i);
        // read type of material
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        if ( ( mat = classFactory.createMaterial(name.c_str(), num, this) ) == NULL ) {
            OOFEM_ERROR("Couldn't create material: %s", name.c_str());
        }

        mat->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > nmat ) ) {
            OOFEM_ERROR("Invalid material number (num=%d)", num);
        }

        if ( !materialList->includes(num) ) {
            materialList->put(num, mat);
        } else {
            OOFEM_ERROR("material entry already exist (num=%d)", num);
        }

        ir->finish();
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Instanciated materials ", nmat)
#  endif

    // read barriers
    nonlocalBarierList->growTo(nbarrier);
    for ( int i = 1; i <= nbarrier; i++ ) {
        NonlocalBarrier *barrier;
        ir = dr->giveInputRecord(DataReader :: IR_nlocBarRec, i);
        // read type of load
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        if ( ( barrier = classFactory.createNonlocalBarrier(name.c_str(), num, this) ) == NULL ) {
            OOFEM_ERROR("Couldn't create barrier: %s", name.c_str());
        }

        barrier->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > nbarrier ) ) {
            OOFEM_ERROR("Invalid barrier number (num=%d)", num);
        }

        if ( !nonlocalBarierList->includes(num) ) {
            nonlocalBarierList->put(num, barrier);
        } else {
            OOFEM_ERROR("barrier entry already exist (num=%d)", num);
        }

        ir->finish();
    }

#  ifdef VERBOSE
    if ( nbarrier ) {
        VERBOSE_PRINT0("Instanciated barriers ", nbarrier);
    }
#  endif

    // read random field generators
    randomFieldGeneratorList->growTo(nrfg);
    for ( int i = 1; i <= nrfg; i++ ) {
        ir = dr->giveInputRecord(DataReader :: IR_nRandomFieldGenRec, i);
        // read type of load
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        rfg = classFactory.createRandomFieldGenerator(name.c_str(), num, this);
        rfg->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > nrfg ) ) {
            OOFEM_ERROR("Invalid generator number (num=%d)", num);
        }

        if ( !randomFieldGeneratorList->includes(num) ) {
            randomFieldGeneratorList->put(num, rfg);
        } else {
            OOFEM_ERROR("generator entry already exist (num=%d)", num);
        }

        ir->finish();
    }

#  ifdef VERBOSE
    if ( nrfg ) {
        VERBOSE_PRINT0("Instanciated random generators ", nbarrier);
    }
#  endif



    // read boundary conditions
    bcList->growTo(nload);
    for ( int i = 1; i <= nload; i++ ) {
        GeneralBoundaryCondition *load;
        ir = dr->giveInputRecord(DataReader :: IR_bcRec, i);
        // read type of load
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        if ( ( load = classFactory.createBoundaryCondition(name.c_str(), num, this) ) == NULL ) {
            OOFEM_ERROR("Couldn't create boundary condition: %s", name.c_str());
        }

        load->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > nload ) ) {
            OOFEM_ERROR("Invalid boundary condition number (num=%d)", num);
        }

        if ( !bcList->includes(num) ) {
            bcList->put(num, load);
        } else {
            OOFEM_ERROR("boundary condition entry already exist (num=%d)", num);
        }

        ir->finish();
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Instanciated BCs ", nload)
#  endif

    // read initial conditions
    icList->growTo(nic);
    for ( int i = 1; i <= nic; i++ ) {
        InitialCondition *ic;
        ir = dr->giveInputRecord(DataReader :: IR_icRec, i);
        // read type of load
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);

        if ( ( ic = new InitialCondition(num, this) ) == NULL ) {
            OOFEM_ERROR("Creation of IC no. %d failed", num);
        }

        ic->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > nic ) ) {
            OOFEM_ERROR("Invalid initial condition number (num=%d)", num);
        }

        if ( !icList->includes(num) ) {
            icList->put(num, ic);
        } else {
            OOFEM_ERROR("initial condition entry already exist (num=%d)", num);
        }

        ir->finish();
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Instanciated ICs ", nic)
#  endif


    // read load time functions
    functionList->growTo(nloadtimefunc);
    for ( int i = 1; i <= nloadtimefunc; i++ ) {
        Function *func;
        ir = dr->giveInputRecord(DataReader :: IR_funcRec, i);
        // read type of func
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);
        if ( ( func = classFactory.createFunction(name.c_str(), num, this) ) == NULL ) {
            OOFEM_ERROR("Couldn't create time function: %s", name.c_str());
        }

        func->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > nloadtimefunc ) ) {
            OOFEM_ERROR("Invalid Function number (num=%d)", num);
        }

        if ( !functionList->includes(num) ) {
            functionList->put(num, func);
        } else {
            OOFEM_ERROR("Function entry already exist (num=%d)", num);
        }

        ir->finish();
    }

#  ifdef VERBOSE
    VERBOSE_PRINT0("Instanciated load-time fncts ", nloadtimefunc)
#  endif

    // read load time functions
    setList->growTo(nset);
    for ( int i = 1; i <= nset; i++ ) {
        ir = dr->giveInputRecord(DataReader :: IR_setRec, i);
        // read type of set
        IR_GIVE_RECORD_KEYWORD_FIELD(ir, name, num);
        // Only one set for now (i don't see any need to ever introduce any other version)
        Set *set = new Set(num, this);
        /*if ( ( set = classFactory.createSet(name.c_str(), num, this) ) == NULL ) {
         *  OOFEM_ERROR("Couldn't create set: %s", name.c_str());
         * }*/

        set->initializeFrom(ir);

        // check number
        if ( ( num < 1 ) || ( num > nset ) ) {
            OOFEM_ERROR("Invalid set number (num=%d)", num);
        }

        if ( !setList->includes(num) ) {
            setList->put(num, set);
        } else {
            OOFEM_ERROR("Set entry already exist (num=%d)", num);
        }

        ir->finish();
    }

#  ifdef VERBOSE
    if ( nset ) {
        VERBOSE_PRINT0("Instanciated sets ", nset);
    }
#  endif

    if ( nxfemman ) {
        xfemManager = new XfemManager(this);
        ir = dr->giveInputRecord(DataReader :: IR_xfemManRec, 1);
        xfemManager->initializeFrom(ir);
        xfemManager->instanciateYourself(dr);
    }
#  ifdef VERBOSE
    if ( nxfemman ) {
        VERBOSE_PRINT0("Instanciated xfem ", nxfemman);
    }
#  endif

    this->topology = NULL;
    if ( topologytype.length() > 0 ) {
        this->topology = classFactory.createTopology(topologytype.c_str(), this);
        if ( !this->topology ) {
            OOFEM_ERROR("Couldn't create topology of type '%s'", topologytype.c_str());
        }

        return this->topology->instanciateYourself(dr);
    }
#  ifdef VERBOSE
    if ( topologytype.length() > 0 ) {
        VERBOSE_PRINT0("Instanciated topologies ", topologytype.length());
    }
#  endif


    if ( nfracman ) {
        fracManager = new FractureManager(this);
        ir = dr->giveInputRecord(DataReader :: IR_fracManRec, 1);
        fracManager->initializeFrom(ir);
        fracManager->instanciateYourself(dr);
    }
#  ifdef VERBOSE
    if ( nfracman ) {
        VERBOSE_PRINT0("Instanciated fracture manager ", nxfemman);
    }
#  endif

    // change internal component references from labels to assigned local numbers
    MapBasedEntityRenumberingFunctor labelToLocNumFunctor(dofManLabelMap, elemLabelMap);
    for ( int i = 1; i <= nnode; i++ ) {
        this->giveDofManager(i)->updateLocalNumbering(labelToLocNumFunctor);
    }

    for ( int i = 1; i <= nelem; i++ ) {
        this->giveElement(i)->updateLocalNumbering(labelToLocNumFunctor);
    }

    for ( int i = 1; i <= nset; i++ ) {
        this->giveSet(i)->updateLocalNumbering(labelToLocNumFunctor);
    }

    return 1;
}


void
Domain :: postInitialize()
{
    // Dofs must be created before dof managers due their post-initialization:
    this->createDofs();

    for ( int i = 1; i <= this->dofManagerList->giveSize(); i++ ) {
        this->dofManagerList->at(i)->postInitialize();
    }

    // New  - in development /JB
    // set element cross sections based on element set definition and set the corresponding
    // material based on the cs
    for ( int i = 1; i <= this->giveNumberOfCrossSectionModels(); i++ ) {
        if ( int setNum = this->giveCrossSection(i)->giveSetNumber() ) {
            Set *set = this->giveSet(setNum);
            const IntArray &elements = set->giveElementList();
            for ( int ielem = 1; ielem <= elements.giveSize(); ++ielem ) {
                Element *element = this->giveElement( elements.at(ielem) );
                element->setCrossSection(i);
            }
        }
    }


    //----

    for ( int i = 1; i <= this->elementList->giveSize(); i++ ) {
        this->elementList->at(i)->postInitialize();
    }
    for ( int i = 1; i <= this->bcList->giveSize(); i++ ) {
        this->bcList->at(i)->postInitialize();
    }
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
    } else if  ( !strncmp(typeName, "3d", 2) ) {
        dType = _3dMode;
    } else {
        OOFEM_ERROR("unknown domainType (%s)", typeName);
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
    if ( destroyOld ) {
        delete this->smoother;
    }

    this->smoother = smoother;
}


void
Domain :: setTopology(TopologyDescription *topo, bool destroyOld)
{
    if ( destroyOld ) {
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
        spatialLocalizer = new OctreeSpatialLocalizer(this);
    }

    return spatialLocalizer;
}


void
Domain :: createDofs()
{
    IntArray dofids;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////// Step 1. Scan all required nodal dofs.
    std :: vector< std :: set< int > > node_dofs( this->giveNumberOfDofManagers() );
    for ( int i = 1; i <= this->giveNumberOfElements(); ++i ) {
        // Scan for all dofs needed by element.
        Element *element = this->giveElement(i);
        for ( int j = 1; j <= element->giveNumberOfNodes(); ++j ) {
            element->giveDefaultDofManDofIDMask(j, dofids);
            for ( int k = 1; k <= dofids.giveSize(); k++ ) {
                node_dofs [ element->giveNode(j)->giveNumber() - 1 ].insert( dofids.at(k) );
            }
        }
    }
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); ++i ) {
        // Nodes can also contain their own list of dofs (typical usecase: RigidArmNode )
        DofManager *dman = this->giveDofManager(i);
        const IntArray *dofids = dman->giveForcedDofIDs();
        if ( dofids ) {
            for ( int k = 1; k <= dofids->giveSize(); ++k ) {
                node_dofs [ i - 1 ].insert( dofids->at(k) );
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
        int c = 0;
        //printf("Dofs in node %d (of %d) = %d\n", i, this->giveNumberOfDofManagers(), node_dofs[i-1].size());
        dman->setNumberOfDofs(0);
        for ( int id: node_dofs [ i - 1 ] ) {
            // Find bc and ic if there are any, otherwise zero.
            int bcid = dof_bc [ i - 1 ].find(id) != dof_bc [ i - 1 ].end() ? dof_bc [ i - 1 ] [ id ] : 0;
            int icid = dof_ic [ i - 1 ].find(id) != dof_ic [ i - 1 ].end() ? dof_ic [ i - 1 ] [ id ] : 0;

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
            Dof *dof = classFactory.createDof(dtype, ++c, dman);
            dof->setDofID((DofIDItem)id);
            dof->setBcId(bcid); // Note: slave dofs and such will simple ignore this.
            dof->setIcId(icid);
            // Slave dofs obtain their weights post-initialization, simple slave dofs must have their master node specified.
            if ( dtype == DT_simpleSlave ) {
                static_cast< SimpleSlaveDof * >(dof)->setMasterDofManagerNum( ( * dman->giveMasterMap() ) [ id ] );
            }
            dman->appendDof(dof);
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
    int nnode, nelem, nmat;

    nnode = this->giveNumberOfDofManagers();
    nelem = this->giveNumberOfElements();
    nmat  = this->giveNumberOfMaterialModels();

    for ( int i = 1; i <= nnode; i++ ) {
        result &= this->giveDofManager(i)->checkConsistency();
    }

    for ( int i = 1; i <= nelem; i++ ) {
        result &= this->giveElement(i)->checkConsistency();
    }

    for ( int i = 1; i <= nmat; i++ ) {
        result &= this->giveMaterial(i)->checkConsistency();
    }

    return result;
}

double
Domain :: giveArea()
{
    double area = 0.0;
    for ( int i = 1; i <= this->giveNumberOfElements(); ++i ) {
        area += this->giveElement(i)->computeArea();
    }

    return area;
}

double
Domain :: giveVolume()
{
    double volume = 0.0;
    for ( int i = 1; i <= this->giveNumberOfElements(); ++i ) {
        volume += this->giveElement(i)->computeVolume();
    }

    return volume;
}

double
Domain :: giveSize()
{
    double volume = 0.0;
    for ( int i = 1; i <= this->giveNumberOfElements(); ++i ) {
        volume += this->giveElement(i)->computeVolumeAreaOrLength();
    }

    return volume;
}

int
Domain :: giveNextFreeDofID(int increment)
{
#ifdef __PARALLEL_MODE
    if ( this->engineeringModel->isParallel() ) {
        OOFEM_ERROR("Additional dof id's not implemented/tested for parallel problems");
    }
#endif
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


#define SAVE_COMPONENTS(size, type, giveMethod)    \
    {                                               \
        for ( int i = 1; i <= size; i++ ) {         \
            type *obj = giveMethod(i);              \
            if ( ( mode & CM_Definition ) ) {       \
                if ( !stream->write( obj->giveInputRecordName() ) ) { \
                    THROW_CIOERR(CIO_IOERR);        \
                }                                   \
            }                                       \
            if ( ( iores = obj->saveContext(stream, mode) ) != CIO_OK ) { \
                THROW_CIOERR(iores);                \
            }                                       \
        }                                           \
    }

#define RESTORE_COMPONENTS(size, type, resizeMethod, creator, giveMethod, setMethod) \
    {                                           \
        if ( mode & CM_Definition ) {           \
            resizeMethod(size);                 \
        }                                       \
        for ( int i = 1; i <= size; i++ ) {     \
            type *obj;                          \
            if ( mode & CM_Definition ) {       \
                std :: string name;               \
                if ( !stream->read(name) ) {    \
                    THROW_CIOERR(CIO_IOERR);    \
                }                               \
                obj = creator(name.c_str(), 0, this); \
                if ( !obj ) {                   \
                    THROW_CIOERR(CIO_BADVERSION); \
                }                               \
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
    int serNum;
    ErrorEstimator *ee;

    // save domain serial number
    serNum = this->giveSerialNumber();
    if ( !stream->write(& serNum, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( mode & CM_Definition ) ) {
        long ncomp [ DOMAIN_NCOMP ];
        ncomp [ 0 ] = this->giveNumberOfDofManagers();
        ncomp [ 1 ] = this->giveNumberOfElements();
        ncomp [ 2 ] = this->giveNumberOfMaterialModels();
        ncomp [ 3 ] = this->giveNumberOfCrossSectionModels();
        ncomp [ 4 ] = this->giveNumberOfBoundaryConditions();
        ncomp [ 5 ] = this->giveNumberOfInitialConditions();
        ncomp [ 6 ] = this->giveNumberOfFunctions();
        ncomp [ 7 ] = this->giveNumberOfNonlocalBarriers();
        ncomp [ 8 ] = this->giveNumberOfRandomFieldGenerators();

        // store number of components
        if ( !stream->write(ncomp, DOMAIN_NCOMP) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        // Have to store materials (and possible other things first) before restoring integration points and such).
        SAVE_COMPONENTS(this->giveNumberOfMaterialModels(), Material, this->giveMaterial);
        SAVE_COMPONENTS(this->giveNumberOfCrossSectionModels(), CrossSection, this->giveCrossSection);
        SAVE_COMPONENTS(this->giveNumberOfInitialConditions(), InitialCondition, this->giveIc);
        SAVE_COMPONENTS(this->giveNumberOfFunctions(), Function, this->giveFunction);
        SAVE_COMPONENTS(this->giveNumberOfNonlocalBarriers(), NonlocalBarrier, this->giveNonlocalBarrier);
        SAVE_COMPONENTS(this->giveNumberOfRandomFieldGenerators(), RandomFieldGenerator, this->giveRandomFieldGenerator);
    }

    // save dof managers
    SAVE_COMPONENTS(this->giveNumberOfDofManagers(), DofManager, this->giveDofManager);
    // elements and corresponding integration points
    SAVE_COMPONENTS(this->giveNumberOfElements(), Element, this->giveElement);
    // boundary conditions
    SAVE_COMPONENTS(this->giveNumberOfBoundaryConditions(), GeneralBoundaryCondition, this->giveBc);

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
    int serNum;
    bool domainUpdated;
    ErrorEstimator *ee;
    long ncomp [ DOMAIN_NCOMP ];

    int nnodes, nelem, nmat, ncs, nbc, nic, nfunc, nnlb, nrfg;


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

        nnodes = ncomp [ 0 ];
        nelem = ncomp [ 1 ];
        nmat  = ncomp [ 2 ];
        ncs   = ncomp [ 3 ];
        nbc   = ncomp [ 4 ];
        nic   = ncomp [ 5 ];
        nfunc = ncomp [ 6 ];
        nnlb  = ncomp [ 7 ];
        nrfg  = ncomp [ 8 ];

        // clear receiver data
        dofManagerList->clear();
        elementList->clear();
        materialList->clear();
        bcList->clear();
        icList->clear();
        functionList->clear();
        nonlocalBarierList->clear();
        randomFieldGeneratorList->clear();
        setList->clear();
        ///@todo Saving and restoring xfemmanagers.
        delete xfemManager;
        xfemManager = NULL;
        //this->clear();

        RESTORE_COMPONENTS(nmat, Material, this->resizeMaterials, classFactory.createMaterial, this->giveMaterial, this->setMaterial);
        RESTORE_COMPONENTS(ncs, CrossSection, this->resizeCrossSectionModels, classFactory.createCrossSection, this->giveCrossSection, this->setCrossSection);
        RESTORE_COMPONENTS(nic, InitialCondition, this->resizeInitialConditions, classFactory.createInitialCondition, this->giveIc, setInitialCondition);
        RESTORE_COMPONENTS(nfunc, Function, resizeFunctions, classFactory.createFunction, giveFunction, setFunction);
        RESTORE_COMPONENTS(nnlb, NonlocalBarrier, resizeNonlocalBarriers, classFactory.createNonlocalBarrier, giveNonlocalBarrier, setNonlocalBarrier);
        RESTORE_COMPONENTS(nrfg, RandomFieldGenerator, resizeRandomFieldGenerators, classFactory.createRandomFieldGenerator, giveRandomFieldGenerator, setRandomFieldGenerator);

        domainUpdated = true;
    } else {
        if ( serNum != this->giveSerialNumber() ) {
            // read corresponding domain
            OOFEM_LOG_INFO("restoring domain %d.%d\n", this->number, this->giveSerialNumber());
            DataReader *domainDr = this->engineeringModel->GiveDomainDataReader(1, this->giveSerialNumber(), contextMode_read);
            this->clear();

            if ( !this->instanciateYourself(domainDr) ) {
                OOFEM_ERROR("domain Instanciation failed");
            }

            delete domainDr;
            domainUpdated = true;
        }

        nnodes = this->giveNumberOfDofManagers();
        nelem = this->giveNumberOfElements();
        nmat = this->giveNumberOfMaterialModels();
        ncs = this->giveNumberOfCrossSectionModels();
        nbc = this->giveNumberOfBoundaryConditions();
        nic = this->giveNumberOfInitialConditions();
        nfunc = this->giveNumberOfFunctions();
        nnlb = this->giveNumberOfNonlocalBarriers();
        nrfg = this->giveNumberOfRandomFieldGenerators();
    }

    RESTORE_COMPONENTS(nnodes, DofManager, this->resizeDofManagers, classFactory.createDofManager, this->giveDofManager, this->setDofManager);
    RESTORE_COMPONENTS(nelem, Element, this->resizeElements, classFactory.createElement, this->giveElement, this->setElement);
    RESTORE_COMPONENTS(nbc, GeneralBoundaryCondition, this->resizeBoundaryConditions, classFactory.createBoundaryCondition, this->giveBc, this->setBoundaryCondition);

    // restore error estimator data
    ee = this->giveErrorEstimator();
    if ( ee ) {
        if ( domainUpdated ) {
            ee->setDomain(this);
        }
        if ( ( iores = ee->restoreContext(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    if ( domainUpdated ) {
        if ( this->smoother ) {
            this->smoother->clear();
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
            OOFEM_ERROR("allocation failed");
        }
    }

    return transactionManager;
}


int Domain :: commitTransactions(DomainTransactionManager *tm)
{
    bool _exist;
    std :: map< int, DofManager * > :: const_iterator dit;
    std :: map< int, Element * > :: const_iterator eit;
    AList< DofManager > *dofManagerList_new = new AList< DofManager >(0);
    AList< Element > *elementList_new = new AList< Element >(0);


    if ( tm->dofmanTransactions.empty() && tm->elementTransactions.empty() ) {
        return 1;
    }

    this->initGlobalDofManMap();
    if ( !tm->dofmanTransactions.empty() ) {
        DofManager *dman;
        for ( dit = tm->dofmanTransactions.begin(); dit != tm->dofmanTransactions.end(); ++dit ) {
            _exist = false;
            if ( dmanMap.find(dit->first) != dmanMap.end() ) {
                _exist = true;
            }

            if ( _exist ) {
                int lnum = dmanMap [ dit->first ]->giveNumber();
                dman = dofManagerList->unlink(lnum);
                dmanMap.erase(dit->first);
                delete dman;
            }

            if ( dit->second ) {
                dmanMap [ dit->first ] = ( DofManager * ) dit->second;
            }
        } // end loop over DofmanTransactions
    }

    this->initGlobalElementMap();
    if ( !tm->elementTransactions.empty() ) {
        int gen;
        Element *elem;

        for ( eit = tm->elementTransactions.begin(); eit != tm->elementTransactions.end(); ++eit ) {
            gen = eit->first;
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

            if ( eit->second ) {
                elementMap [ gen ] = ( Element * ) eit->second;
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
        DofManager *dm = dmanMap [ num ];
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
} // end namespace oofem
