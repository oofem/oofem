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
 *               Copyright (C) 1993 - 2024   Borek Patzak
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

#include "additivemanufacturingproblem.h"
#include "GCodeParser.h"
#include "heavisidetimefunction.h"
#include "calculatorfunction.h"
#include "initialcondition.h"
#include "unknownnumberingscheme.h"
#include "masterdof.h"
#include "initialcondition.h"
#include "load.h"
#include "boundaryload.h"
#include "boundarycondition.h"
#include "node.h"
#include "dofmanager.h"
#include "freeconstantsurfaceload.h"
#include "dynamicinputrecord.h"
#include "tm/EngineeringModels/transienttransportproblem.h"
#include "connectivitytable.h"
#include "spatiallocalizer.h"
#include "primaryfield.h"
#include <vector>
#include <chrono>
#include <thread>

#include "engngm.h"
#include "timestep.h"
#include "function.h"
#include "metastep.h"
#include "exportmodulemanager.h"
#include "mathfem.h"
#include "oofemtxtdatareader.h"
#include "util.h"
#include "verbose.h"
#include "classfactory.h"
#include "domain.h"

#include <stdlib.h>

// #include "datarecord.h"


#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_EngngModel( AdditiveManufacturingProblem );

bool add_node_if_not_exists2( EngngModel *emodel, const VoxelNode &cn )
{
    Domain *d  = emodel->giveDomain( 1 );
    int nodeId = cn.id;

    std::unique_ptr<DofManager> dman( classFactory.createDofManager( "node", nodeId, d ) );
    if ( !dman ) {
        OOFEM_ERROR( "Couldn't create node %d\n", nodeId );
    }

    FloatArray coords = { cn.coords[0] / 1000, cn.coords[1] / 1000, cn.coords[2] / 1000 };

    dman->setCoordinates( coords );
    dman->setGlobalNumber( nodeId );


    // Set for node
    int nSetsBefore = d->giveNumberOfSets();
    d->resizeSets( nSetsBefore + 1 );
    int setId = nSetsBefore + 1;

    std::unique_ptr<Set> set = std::make_unique<Set>( setId, d );
    IntArray nodeIds         = { nodeId };
    set->setNumber( setId );
    set->setNodeList( nodeIds );
    d->setSet( setId, std::move( set ) );

    // create a dof
    Dof *dof = classFactory.createDof( DT_master, T_f, dman.get() );

    // Initial condition
    int nICBefore = d->giveNumberOfInitialConditions();
    int ICId      = nICBefore + 1;
    d->resizeInitialConditions( ICId );
    std::unique_ptr<InitialCondition> ic = std::make_unique<InitialCondition>( ICId, d );
    DynamicInputRecord icInput           = DynamicInputRecord( _IFT_InitialCondition_Name, ICId );

    Dictionary props = Dictionary();
    props.add( 'u', 235. );

    icInput.setField( props, _IFT_InitialCondition_conditions );
    IntArray dofIds = { T_f };
    icInput.setField( dofIds, _IFT_InitialCondition_dofs );
    icInput.setField( setId, _IFT_InitialCondition_set );
    ic->initializeFrom( icInput );

    d->setInitialCondition( ICId, std::move( ic ) );

    if ( coords[2] <= 1e-12 ) {
        // BC
        int nBCBefore = d->giveNumberOfBoundaryConditions();
        d->resizeBoundaryConditions( nBCBefore + 1 );
        int bcId = nBCBefore + 1;

        std::unique_ptr<GeneralBoundaryCondition> bc( classFactory.createBoundaryCondition( "BoundaryCondition", bcId, d ) );

        DynamicInputRecord myInput = DynamicInputRecord( _IFT_BoundaryCondition_Name, bcId );
        myInput.setField( 1, _IFT_GeneralBoundaryCondition_timeFunct );

        FloatArray vals = { 60.0 };
        IntArray dofs   = { T_f };
        myInput.setField( dofs, _IFT_GeneralBoundaryCondition_dofs );
        myInput.setField( vals, _IFT_BoundaryCondition_values );
        myInput.setField( setId, _IFT_GeneralBoundaryCondition_set );

        bc->initializeFrom( myInput );

        bc->setNumber( bcId );
        bc->postInitialize();

        d->setBoundaryCondition( bcId, std::move( bc ) );
        // set bcid to the created dof
        dof->setBcId( bcId ); // Note: slave dofs and such will simple ignore this.
    } else {

        // Fix node when its not activated
        int nTFBefore = d->giveNumberOfFunctions();
        d->resizeFunctions( nTFBefore + 1 );
        std::unique_ptr<CalculatorFunction> tf = std::make_unique<CalculatorFunction>( nTFBefore + 1, d );

        DynamicInputRecord myInput = DynamicInputRecord( _IFT_CalculatorFunction_Name, nTFBefore + 1 );
        myInput.setField( "1-h(" + std::to_string( cn.timeActivated ) + ")", _IFT_CalculatorFunction_f );
        tf->initializeFrom( myInput );

        d->setFunction( nTFBefore + 1, std::move( tf ) );

        int nBCBefore = d->giveNumberOfBoundaryConditions();
        d->resizeBoundaryConditions( nBCBefore + 1 );
        int bcId = nBCBefore + 1;

        std::unique_ptr<GeneralBoundaryCondition> bc( classFactory.createBoundaryCondition( "BoundaryCondition", bcId, d ) );

        DynamicInputRecord myInput2 = DynamicInputRecord( _IFT_BoundaryCondition_Name, bcId );
        myInput2.setField( 1, _IFT_GeneralBoundaryCondition_timeFunct );
        myInput2.setField( nTFBefore + 1, _IFT_GeneralBoundaryCondition_isImposedTimeFunct );

        FloatArray vals = { 235.0 };
        IntArray dofs   = { T_f };
        myInput2.setField( dofs, _IFT_GeneralBoundaryCondition_dofs );
        myInput2.setField( vals, _IFT_BoundaryCondition_values );
        myInput2.setField( setId, _IFT_GeneralBoundaryCondition_set );

        bc->initializeFrom( myInput2 );

        bc->setNumber( bcId );
        bc->postInitialize();

        d->setBoundaryCondition( bcId, std::move( bc ) );
        // set bcid to the created dof
        dof->setBcId( bcId ); // Note: slave dofs and such will simple ignore this.
    }


    dman->appendDof( dof );
    dman->postInitialize();
    d->setDofManager( nodeId, std::move( dman ) );

    return true;
}

void add_element_if_not_exists2( EngngModel *emodel, Voxel &cn )
{
    Domain *d = emodel->giveDomain( 1 );

    // Add new element
    int elId = cn.id;

    std::unique_ptr<Element> elem( classFactory.createElement( "brick1ht", elId, d ) );
    if ( !elem ) {
        OOFEM_ERROR( "Couldn't create element %d\n", elId );
    }

    IntArray enodes = {
        cn.nodes[0], cn.nodes[1], cn.nodes[2], cn.nodes[3],
        cn.nodes[4], cn.nodes[5], cn.nodes[6], cn.nodes[7]
    };

    elem->setDofManagers( enodes );
    elem->setGlobalNumber( elId );
    elem->setCrossSection( 1 );
    elem->setParallelMode( elementParallelMode::Element_local ); // TODO: WHY SOME ELEMENTS HAVE SET REMOTE????
    elem->postInitialize(); // this is needed to allocate the integration point array

    int nTFBefore = d->giveNumberOfFunctions();
    d->resizeFunctions( nTFBefore + 1 );
    std::unique_ptr<HeavisideTimeFunction> tf = std::make_unique<HeavisideTimeFunction>( nTFBefore + 1, d );

    DynamicInputRecord myInput = DynamicInputRecord( _IFT_HeavisideTimeFunction_Name, nTFBefore + 1 );
    myInput.setField( cn.time_activated(), _IFT_HeavisideTimeFunction_origin );
    myInput.setField( 1.0, _IFT_HeavisideTimeFunction_value );
    tf->initializeFrom( myInput );
    elem->setActivityTimeFunctionNumber( nTFBefore + 1 );

    d->setFunction( nTFBefore + 1, std::move( tf ) );
    d->setElement( elId, std::move( elem ) );

    // Add new node set
    // Used to apply initial boundary condition to the new DOFs
    int nSetsBefore = d->giveNumberOfSets();
    d->resizeSets( nSetsBefore + 6 );

    // Add 6 loads (1 per element face)
    int nBCBefore = d->giveNumberOfBoundaryConditions();
    d->resizeBoundaryConditions( nBCBefore + 6 );

    for ( int i = 0; i < 6; i++ ) {
        int setId = nSetsBefore + 1 + i;

        std::unique_ptr<Set> set = std::make_unique<Set>( setId, d );
        IntArray bounds          = { elId, 1 + i };
        set->setNumber( setId );
        set->setBoundaryList( bounds );
        // set->setNodeList(nNodeIds);
        d->setSet( setId, std::move( set ) );

        int bcid = nBCBefore + 1 + i;
        std::unique_ptr<GeneralBoundaryCondition> bc( classFactory.createBoundaryCondition( "freeconstantsurfaceload", bcid, d ) );

        DynamicInputRecord myInput = DynamicInputRecord( _IFT_FreeConstantSurfaceLoad_Name, bcid );
        myInput.setField( 1, _IFT_GeneralBoundaryCondition_timeFunct );
        myInput.setField( 3, _IFT_BoundaryLoad_loadtype );
        Dictionary props = Dictionary();
        props.add( 'a', 10. ); // 97 == 'a'
        myInput.setField( props, _IFT_BoundaryLoad_properties );

        // todo: UNCOMMENT
        FloatArray comps = { 30. }; // ambient temp
        myInput.setField( comps, _IFT_Load_components );
        myInput.setField( setId, _IFT_GeneralBoundaryCondition_set );

        bc->initializeFrom( myInput );

        bc->setNumber( bcid );
        bc->postInitialize();

        d->setBoundaryCondition( bcid, std::move( bc ) );
    }
}

bool add_sm_node_if_not_exists2( EngngModel *emodel, const VoxelNode &cn )
{
    Domain *d = emodel->giveDomain( 1 );

    int nodeId = cn.id;

    // was already added?
    if ( d->dofManagerList[cn.id - 1] != nullptr ) return false;

    std::unique_ptr<DofManager> dman( classFactory.createDofManager( "node", nodeId, d ) );
    if ( !dman ) {
        OOFEM_ERROR( "Couldn't create node %d\n", nodeId );
    }

    FloatArray coords = { cn.coords[0] / 1000, cn.coords[1] / 1000, cn.coords[2] / 1000 };

    /*DynamicInputRecord myInput = DynamicInputRecord(_IFT_Node_Name, nodeId);
    IntArray mBCs = {1,1,1};
    myInput.setField(mBCs, _IFT_DofManager_bc);
    myInput.setField(coords, _IFT_Node_coords);
    dman->initializeFrom(myInput);*/

    dman->setCoordinates( coords );
    dman->setGlobalNumber( nodeId );
    // dman->appendDof( new MasterDof(dman.get(), D_u) );
    // dman->appendDof( new MasterDof(dman.get(), D_v) );
    // dman->appendDof( new MasterDof(dman.get(), D_w) );

    // d->setDofManager(nodeId, std::move(dman));

    // create a dof
    Dof *dofU = classFactory.createDof( DT_master, D_u, dman.get() );
    Dof *dofV = classFactory.createDof( DT_master, D_v, dman.get() );
    Dof *dofW = classFactory.createDof( DT_master, D_w, dman.get() );

    if ( coords[2] <= 1e-12 ) {
        // Set for BC
        int nSetsBefore = d->giveNumberOfSets();
        int nBCBefore   = d->giveNumberOfBoundaryConditions();
        int bcId;

        /*if (nSetsBefore > 0 ){
            auto set = d->giveSet(nSetsBefore);
            auto nodeList = set->giveSpecifiedNodeList();
            nodeList.followedBy(nodeId,1);
            set->setNodeList(nodeList);
            bcId = nBCBefore;
        }
        else {*/
        d->resizeSets( nSetsBefore + 1 );
        int setId = nSetsBefore + 1;

        std::unique_ptr<Set> set = std::make_unique<Set>( setId, d );
        IntArray nodeIds         = { nodeId };
        set->setNumber( setId );
        set->setNodeList( nodeIds );
        d->setSet( setId, std::move( set ) );

        // BC
        // int nBCBefore = d->giveNumberOfBoundaryConditions();
        d->resizeBoundaryConditions( nBCBefore + 1 );
        bcId = nBCBefore + 1;

        std::unique_ptr<GeneralBoundaryCondition> bc( classFactory.createBoundaryCondition( "BoundaryCondition", bcId, d ) );

        DynamicInputRecord myInput = DynamicInputRecord( _IFT_BoundaryCondition_Name, bcId );
        myInput.setField( 1, _IFT_GeneralBoundaryCondition_timeFunct );

        FloatArray vals = { 0.0, 0.0, 0.0 };
        IntArray dofs   = { D_u, D_v, D_w };
        myInput.setField( dofs, _IFT_GeneralBoundaryCondition_dofs );
        myInput.setField( vals, _IFT_BoundaryCondition_values );
        myInput.setField( setId, _IFT_GeneralBoundaryCondition_set );

        bc->initializeFrom( myInput );

        bc->setNumber( bcId );
        bc->postInitialize();

        d->setBoundaryCondition( bcId, std::move( bc ) );

        //}

        // set bcid to the created dof
        dofU->setBcId( bcId ); // Note: slave dofs and such will simple ignore this.
        dofV->setBcId( bcId ); // Note: slave dofs and such will simple ignore this.
        dofW->setBcId( bcId ); // Note: slave dofs and such will simple ignore this.

    } else {
        // Set for node
        int nSetsBefore = d->giveNumberOfSets();
        d->resizeSets( nSetsBefore + 1 );
        int setId = nSetsBefore + 1;

        std::unique_ptr<Set> set = std::make_unique<Set>( setId, d );
        IntArray nodeIds         = { nodeId };
        set->setNumber( setId );
        set->setNodeList( nodeIds );
        d->setSet( setId, std::move( set ) );

        // Fix node when its not activated
        int nTFBefore = d->giveNumberOfFunctions();
        d->resizeFunctions( nTFBefore + 1 );
        std::unique_ptr<CalculatorFunction> tf = std::make_unique<CalculatorFunction>( nTFBefore + 1, d );

        DynamicInputRecord myInput = DynamicInputRecord( _IFT_CalculatorFunction_Name, nTFBefore + 1 );
        myInput.setField( "1-h(" + std::to_string( cn.timeActivated ) + ")", _IFT_CalculatorFunction_f );
        // myInput.setField("0", _IFT_CalculatorFunction_f);
        // std::cout << "Node " <<nodeId << " activated at " << cn.activatedAt << std::endl;
        // myInput.setField(1.0, _IFT_HeavisideTimeFunction_value);
        tf->initializeFrom( myInput );
        // elem->setActivityTimeFunctionNumber(nTFBefore+1);

        d->setFunction( nTFBefore + 1, std::move( tf ) );

        int nBCBefore = d->giveNumberOfBoundaryConditions();
        d->resizeBoundaryConditions( nBCBefore + 1 );
        int bcId = nBCBefore + 1;

        std::unique_ptr<GeneralBoundaryCondition> bc( classFactory.createBoundaryCondition( "BoundaryCondition", bcId, d ) );

        DynamicInputRecord myInput2 = DynamicInputRecord( _IFT_BoundaryCondition_Name, bcId );
        myInput2.setField( 1, _IFT_GeneralBoundaryCondition_timeFunct );
        myInput2.setField( nTFBefore + 1, _IFT_GeneralBoundaryCondition_isImposedTimeFunct );

        FloatArray vals = { 0.0, 0.0, 0.0 };
        IntArray dofs   = { D_u, D_v, D_w };
        myInput2.setField( dofs, _IFT_GeneralBoundaryCondition_dofs );
        myInput2.setField( vals, _IFT_BoundaryCondition_values );
        myInput2.setField( setId, _IFT_GeneralBoundaryCondition_set );

        bc->initializeFrom( myInput2 );

        bc->setNumber( bcId );
        // bc->setIsImposedTimeFunctionNumber(nTFBefore+1);
        bc->postInitialize();

        d->setBoundaryCondition( bcId, std::move( bc ) );
        // set bcid to the created dof
        dofU->setBcId( bcId ); // Note: slave dofs and such will simple ignore this.
        dofV->setBcId( bcId );
        dofW->setBcId( bcId );
    }

    dman->appendDof( dofU );
    dman->appendDof( dofV );
    dman->appendDof( dofW );

    d->setDofManager( nodeId, std::move( dman ) );
    return true;
}


void add_sm_element_if_not_exists2( EngngModel *emodel, Voxel &cn )
{
    Domain *d = emodel->giveDomain( 1 );

    // Add new element
    int elId = cn.id;

    // was already added?
    if ( d->elementList[cn.id - 1] != nullptr ) return;

    std::unique_ptr<Element> elem( classFactory.createElement( "LSpace", elId, d ) );
    if ( !elem ) {
        OOFEM_ERROR( "Couldn't create element %d\n", elId );
    }

    IntArray enodes = {
        cn.nodes[0], cn.nodes[1], cn.nodes[2], cn.nodes[3],
        cn.nodes[4], cn.nodes[5], cn.nodes[6], cn.nodes[7]
    };

    int nTFBefore = d->giveNumberOfFunctions();
    d->resizeFunctions( nTFBefore + 1 );
    std::unique_ptr<HeavisideTimeFunction> tf = std::make_unique<HeavisideTimeFunction>( nTFBefore + 1, d );

    DynamicInputRecord myInput = DynamicInputRecord( _IFT_HeavisideTimeFunction_Name, nTFBefore + 1 );
    myInput.setField( cn.time_activated(), _IFT_HeavisideTimeFunction_origin );
    myInput.setField( 1.0, _IFT_HeavisideTimeFunction_value );
    tf->initializeFrom( myInput );
    elem->setActivityTimeFunctionNumber( nTFBefore + 1 );

    d->setFunction( nTFBefore + 1, std::move( tf ) );

    elem->setDofManagers( enodes );
    elem->setGlobalNumber( elId );
    elem->setCrossSection( 1 );
    elem->setParallelMode( elementParallelMode::Element_local );
    elem->postInitialize(); // this is needed to allocate the integration point array

    d->setElement( elId, std::move( elem ) );
}

AdditiveManufacturingProblem ::AdditiveManufacturingProblem( int i, EngngModel *_master ) :
    EngngModel( i, _master ),
    adaptiveStepLength( false ),
    minStepLength( 0. ),
    maxStepLength( 0. ),
    reqIterations( 0. ),
    adaptiveStepSince( 0. ),
    endOfTimeOfInterest( 0. ),
    prevStepLength( 0. ),
    currentStepLength( 0. )
{
    ndomains = 1; // domain is needed to store the time step function

    dtFunction        = 0;
    stepMultiplier    = 1.;
    timeDefinedByProb = 0;

    // fclose(stdout);
    this->totalTimer.startTimer();
}

AdditiveManufacturingProblem ::~AdditiveManufacturingProblem()
{
    this->totalTimer.stopTimer();
    OOFEM_LOG_INFO( "TOTAL TIME SOLVING: %.2fs\n", this->totalTimer.getUtime() );
}

///////////
int AdditiveManufacturingProblem ::instanciateYourself( DataReader &dr, InputRecord &ir, const char *dataOutputFileName, const char *desc )
{
    int result;
    result = EngngModel ::instanciateYourself( dr, ir, dataOutputFileName, desc );
    ir.finish();
    // instanciate slave problems
    result &= this->instanciateSlaveProblems();
    return result;
}

int AdditiveManufacturingProblem ::instanciateDefaultMetaStep( InputRecord &ir )
{
    if ( timeDefinedByProb ) {
        /* just set a nonzero number of steps;
         * needed for instanciateDefaultMetaStep to pass; overall has no effect as time stepping is deteremined by slave
         */
        this->numberOfSteps = 1;
    }
    EngngModel ::instanciateDefaultMetaStep( ir );
    // there are no slave problems initiated so far, the overall metaStep will defined in a slave problem instantiation
    return 1;
}

int AdditiveManufacturingProblem ::instanciateSlaveProblems()
{
    // first instantiate master problem if defined
    // EngngModel *timeDefProb = NULL;
    emodelList.resize( inputStreamNames.size() );
    if ( timeDefinedByProb ) {
        OOFEMTXTDataReader dr( inputStreamNames[timeDefinedByProb - 1] );
        std ::unique_ptr<EngngModel> prob( InstanciateProblem( dr, this->pMode, this->contextOutputMode, this ) );
        // timeDefProb = prob.get();
        emodelList[timeDefinedByProb - 1] = std ::move( prob );
    }

    for ( int i = 1; i <= (int)inputStreamNames.size(); i++ ) {
        if ( i == timeDefinedByProb ) {
            continue;
        }

        OOFEMTXTDataReader dr( inputStreamNames[i - 1] );
        // the slave problem dictating time needs to have attribute master=NULL, other problems point to the dictating slave
        std ::unique_ptr<EngngModel> prob( InstanciateProblem( dr, this->pMode, this->contextOutputMode, this ) );
        emodelList[i - 1] = std ::move( prob );
    }

    return 1;
}


void AdditiveManufacturingProblem ::initializeFrom( InputRecord &ir )
{
    OOFEM_LOG_INFO( "Starting Additive Manufacturing solver\n" );

    IR_GIVE_FIELD( ir, this->gCodeFilePath, _IFT_AdditiveManufacturingProblem_gcode );
    IR_GIVE_FIELD( ir, this->stepX, _IFT_AdditiveManufacturingProblem_stepx );
    IR_GIVE_FIELD( ir, this->stepY, _IFT_AdditiveManufacturingProblem_stepy );
    IR_GIVE_FIELD( ir, this->stepZ, _IFT_AdditiveManufacturingProblem_stepz );

    IR_GIVE_OPTIONAL_FIELD( ir, this->minVOF, _IFT_AdditiveManufacturingProblem_minvof );

    IR_GIVE_OPTIONAL_FIELD( ir, this->skipSM, _IFT_AdditiveManufacturingProblem_skipsm );

    PrinterOptions po;
    po.steps            = { this->stepX, this->stepY, this->stepZ };
    po.sizes            = { (int)std::ceil( 250.0 / po.steps[0] ),
                   (int)std::ceil( 250.0 / po.steps[1] ),
                   (int)std::ceil( 250.0 / po.steps[2] ) };
    po.filamentDiameter = 1.75;
    po.layerHeightModel = LayerHeightModel::Constant;
    po.layerHeight      = 0.2;

    auto pr = Printer( po );

    OOFEM_LOG_INFO( "\n\nG-code file: %s\n", this->gCodeFilePath.c_str() );

    GCodeParser parser;
    std::vector<GCodeCommand> commands = parser.parseFile( this->gCodeFilePath );

    double queue_size = 3;

    for ( int i = 1; i < queue_size; i++ ) {
        pr.addCommandToQueue( commands[i] );
    }

    for ( size_t i = 0; i < commands.size(); i++ ) {
        size_t j = i + queue_size;

        pr.processCommand( commands[i] );

        if ( j < commands.size() )
            pr.addCommandToQueue( commands[j] );

        pr.popCommandFromQueue();
    }

    this->printer = pr;

    OOFEM_LOG_INFO( "\nFinished G-code parsing\n\n" );

    IR_GIVE_FIELD( ir, numberOfSteps, _IFT_EngngModel_nsteps );
    if ( numberOfSteps <= 0 ) {
        throw ValueInputException( ir, _IFT_EngngModel_nsteps, "nsteps must be > 0" );
    }
    if ( ir.hasField( _IFT_AdditiveManufacturingProblem_deltat ) ) {
        EngngModel ::initializeFrom( ir );
        IR_GIVE_FIELD( ir, deltaT, _IFT_AdditiveManufacturingProblem_deltat );
        dtFunction = 0;
    } else if ( ir.hasField( _IFT_AdditiveManufacturingProblem_prescribedtimes ) ) {
        EngngModel ::initializeFrom( ir );
        IR_GIVE_FIELD( ir, discreteTimes, _IFT_AdditiveManufacturingProblem_prescribedtimes );
        dtFunction = 0;
    } else if ( ir.hasField( _IFT_AdditiveManufacturingProblem_dtf ) ) {
        IR_GIVE_OPTIONAL_FIELD( ir, dtFunction, _IFT_AdditiveManufacturingProblem_dtf );
    } else {
        IR_GIVE_FIELD( ir, timeDefinedByProb, _IFT_AdditiveManufacturingProblem_timeDefinedByProb );
    }

    if ( ir.hasField( _IFT_AdditiveManufacturingProblem_adaptiveStepLength ) ) {
        adaptiveStepLength  = true;
        this->minStepLength = 0.;
        IR_GIVE_OPTIONAL_FIELD( ir, minStepLength, _IFT_AdditiveManufacturingProblem_minsteplength );
        this->maxStepLength = 1.e32;
        IR_GIVE_OPTIONAL_FIELD( ir, maxStepLength, _IFT_AdditiveManufacturingProblem_maxsteplength );
        this->reqIterations = 1;
        IR_GIVE_OPTIONAL_FIELD( ir, reqIterations, _IFT_AdditiveManufacturingProblem_reqiterations );
        this->endOfTimeOfInterest = 1.e32;
        IR_GIVE_OPTIONAL_FIELD( ir, endOfTimeOfInterest, _IFT_AdditiveManufacturingProblem_endoftimeofinterest );
        this->adaptiveStepSince = 0.;
        IR_GIVE_OPTIONAL_FIELD( ir, adaptiveStepSince, _IFT_AdditiveManufacturingProblem_adaptivestepsince );
    }


    IR_GIVE_OPTIONAL_FIELD( ir, stepMultiplier, _IFT_AdditiveManufacturingProblem_stepmultiplier );
    if ( stepMultiplier < 0 ) {
        throw ValueInputException( ir, _IFT_AdditiveManufacturingProblem_stepmultiplier, "stepMultiplier must be > 0" );
    }

    //    timeLag = 0.;
    //    IR_GIVE_OPTIONAL_FIELD(ir, timeLag, _IFT_AdditiveManufacturingProblem_timeLag);

    inputStreamNames.resize( 2 );
    if ( ir.hasField( _IFT_AdditiveManufacturingProblem_prob3 ) ) {
        inputStreamNames.resize( 3 );
    }

    IR_GIVE_FIELD( ir, inputStreamNames[0], _IFT_AdditiveManufacturingProblem_prob1 );
    IR_GIVE_FIELD( ir, inputStreamNames[1], _IFT_AdditiveManufacturingProblem_prob2 );
    if ( inputStreamNames.size() == 3 ) {
        IR_GIVE_OPTIONAL_FIELD( ir, inputStreamNames[2], _IFT_AdditiveManufacturingProblem_prob3 );
    }


    renumberFlag = true; // The staggered problem itself should always try to check if the sub-problems needs renumbering.

    coupledModels.resize( 3 );
    IR_GIVE_OPTIONAL_FIELD( ir, this->coupledModels, _IFT_AdditiveManufacturingProblem_coupling );


    if ( dtFunction < 1 ) {
        ndomains = 0;
        domainNeqs.clear();
        domainPrescribedNeqs.clear();
        domainList.clear();
    }

    suppressOutput = ir.hasField( _IFT_EngngModel_suppressOutput );

    if ( suppressOutput ) {
        printf( "Suppressing output.\n" );
    } else {

        if ( ( outputStream = fopen( this->dataOutputFileName.c_str(), "w" ) ) == NULL ) {
            throw ValueInputException( ir, "None", "can't open output file: " + this->dataOutputFileName );
        }

        fprintf( outputStream, "%s", PRG_HEADER );
        fprintf( outputStream, "\nStarting analysis on: %s\n", ctime( &this->startTime ) );
        fprintf( outputStream, "%s\n", simulationDescription.c_str() );

#ifdef __MPI_PARALLEL_MODE
        if ( this->isParallel() ) {
            fprintf( outputStream, "Problem rank is %d/%d on %s\n\n", this->rank, this->numProcs, this->processor_name );
        }
#endif
    }
}


void AdditiveManufacturingProblem ::updateAttributes( MetaStep *mStep )
{
    auto &ir = mStep->giveAttributesRecord();

    EngngModel ::updateAttributes( mStep );

    // update attributes of slaves
    for ( auto &emodel : emodelList ) {
        emodel->updateAttributes( mStep );
    }

    if ( !timeDefinedByProb ) {
        if ( ir.hasField( _IFT_AdditiveManufacturingProblem_deltat ) ) {
            IR_GIVE_FIELD( ir, deltaT, _IFT_AdditiveManufacturingProblem_deltat );
            IR_GIVE_OPTIONAL_FIELD( ir, dtFunction, _IFT_AdditiveManufacturingProblem_dtf );
            IR_GIVE_OPTIONAL_FIELD( ir, stepMultiplier, _IFT_AdditiveManufacturingProblem_stepmultiplier );
            if ( stepMultiplier < 0 ) {
                OOFEM_ERROR( "stepMultiplier must be > 0" )
            }
        } else if ( ir.hasField( _IFT_AdditiveManufacturingProblem_prescribedtimes ) ) {
            IR_GIVE_FIELD( ir, discreteTimes, _IFT_AdditiveManufacturingProblem_prescribedtimes );
        }
    }
}

Function *AdditiveManufacturingProblem ::giveDtFunction()
// Returns the load-time function of the receiver.
{
    if ( !dtFunction ) {
        return NULL;
    }

    return giveDomain( 1 )->giveFunction( dtFunction );
}

double
AdditiveManufacturingProblem ::giveDeltaT( int n )
{
    if ( giveDtFunction() ) {
        return giveDtFunction()->evaluateAtTime( n );
    }

    // in the first step the time increment is taken as the initial, user-specified value
    if ( stepMultiplier != 1 && currentStep ) {
        if ( currentStep->giveNumber() >= 2 ) {
            return ( currentStep->giveTargetTime() * ( stepMultiplier ) );
        }
    }

    if ( discreteTimes.giveSize() > 0 ) {
        return this->giveDiscreteTime( n ) - this->giveDiscreteTime( n - 1 );
    }

    if ( adaptiveStepLength ) {
        EngngModel *sp;
        int nite              = 1;
        double adjustedDeltaT = deltaT;

        if ( currentStep != NULL ) {
            if ( currentStep->giveNumber() != 0 ) {
                // return prescribed deltaT for times until time = adaptiveStepSince
                // can be used for consecutive force loading applied in a specified number of steps
                if ( !( currentStep->giveTargetTime() > this->adaptiveStepSince ) ) {
                    return adjustedDeltaT;
                }

                for ( int i = 1; i <= this->giveNumberOfSlaveProblems(); i++ ) {
                    sp   = this->giveSlaveProblem( i );
                    nite = max( sp->giveCurrentNumberOfIterations(), nite );
                }

                if ( nite > reqIterations ) {
                    adjustedDeltaT = this->prevStepLength * reqIterations / nite;
                } else {
                    adjustedDeltaT = this->prevStepLength * sqrt( sqrt( (double)reqIterations / (double)nite ) );
                }

                if ( adjustedDeltaT > maxStepLength ) {
                    adjustedDeltaT = maxStepLength;
                }

                if ( adjustedDeltaT < minStepLength ) {
                    adjustedDeltaT = minStepLength;
                }
            }
        }

        this->currentStepLength = adjustedDeltaT;

        return adjustedDeltaT;
    }

    return deltaT;
}

double
AdditiveManufacturingProblem ::giveDiscreteTime( int iStep )
{
    if ( ( iStep > 0 ) && ( iStep <= discreteTimes.giveSize() ) ) {
        return ( discreteTimes.at( iStep ) );
    }

    if ( ( iStep == 0 ) && ( iStep <= discreteTimes.giveSize() ) ) {
        return ( 0.0 );
    }

    OOFEM_ERROR( "invalid iStep" );
    return 0.0;
}

TimeStep *
AdditiveManufacturingProblem ::giveCurrentStep( bool force )
{
    if ( timeDefinedByProb ) {
        return emodelList[timeDefinedByProb - 1].get()->giveCurrentStep( true );
    } else {
        return EngngModel ::giveCurrentStep();
    }
}

TimeStep *
AdditiveManufacturingProblem ::givePreviousStep( bool force )
{
    if ( timeDefinedByProb ) {
        return emodelList[timeDefinedByProb - 1].get()->givePreviousStep( true );
    } else {
        return EngngModel ::givePreviousStep();
    }
}

TimeStep *
AdditiveManufacturingProblem ::giveSolutionStepWhenIcApply( bool force )
{
    if ( timeDefinedByProb ) {
        return emodelList[timeDefinedByProb - 1].get()->giveSolutionStepWhenIcApply( true );
    } else {
        if ( !stepWhenIcApply ) {
            int nFirst = giveNumberOfFirstStep();
            // stepWhenIcApply = std::make_unique<TimeStep>(giveNumberOfTimeStepWhenIcApply(), this, 0, -giveDeltaT(nFirst), giveDeltaT(nFirst), 0); //previous version for [-dt, 0]
            stepWhenIcApply = std::make_unique<TimeStep>( giveNumberOfTimeStepWhenIcApply(), this, 0, 0., giveDeltaT( nFirst ), 0 ); // now go from [0, dt]
        }

        return stepWhenIcApply.get();
    }
}


EngngModel *
AdditiveManufacturingProblem ::giveTimeControl()
{
    if ( !timeDefinedByProb ) {
        return this;
    } else { // time dictated by slave problem
        return this->giveSlaveProblem( timeDefinedByProb );
    }
}


int AdditiveManufacturingProblem ::giveNumberOfFirstStep( bool force )
{
    if ( timeDefinedByProb && !force ) {
        return emodelList[timeDefinedByProb - 1].get()->giveNumberOfFirstStep( true );
    } else {
        return EngngModel ::giveNumberOfFirstStep( force );
    }
}


TimeStep *
AdditiveManufacturingProblem ::giveNextStep()
{
    int istep                = this->giveNumberOfFirstStep();
    double totalTime         = 0;
    StateCounterType counter = 1;

    if ( !currentStep ) {
        // first step -> generate initial step
        currentStep = std::make_unique<TimeStep>( *giveSolutionStepWhenIcApply() );
    }

    double dt    = this->giveDeltaT( currentStep->giveNumber() + 1 );
    istep        = currentStep->giveNumber() + 1;
    totalTime    = currentStep->giveTargetTime() + this->giveDeltaT( istep );
    counter      = currentStep->giveSolutionStateCounter() + 1;
    previousStep = std ::move( currentStep );
    currentStep  = std::make_unique<TimeStep>( *previousStep, dt );

    if ( ( totalTime >= this->endOfTimeOfInterest ) && this->adaptiveStepLength ) {
        totalTime = this->endOfTimeOfInterest;
        OOFEM_LOG_INFO( "\n==================================================================\n" );
        OOFEM_LOG_INFO( "\nAdjusting time step length to: %lf \n\n", totalTime - previousStep->giveTargetTime() );
        currentStep         = std::make_unique<TimeStep>( istep, this, 1, totalTime, totalTime - previousStep->giveTargetTime(), counter );
        this->numberOfSteps = istep;
    } else {
        if ( this->adaptiveStepLength ) {
            OOFEM_LOG_INFO( "\n==================================================================\n" );
            OOFEM_LOG_INFO( "\nAdjusting time step length to: %lf \n\n", totalTime - previousStep->giveTargetTime() );
        }
        currentStep = std::make_unique<TimeStep>( istep, this, 1, totalTime, totalTime - previousStep->giveTargetTime(), counter );
    }

    // time and dt variables are set eq to 0 for statics - has no meaning
    return currentStep.get();
}


void AdditiveManufacturingProblem ::solveYourself()
{
    EngngModel *sp;
    sp = giveTimeControl();

    int smstep = 1, sjstep = 1;
    this->timer.startTimer( EngngModelTimer ::EMTT_AnalysisTimer );

    // Obtain eng models for TM and SM
    TransientTransportProblem *em_tm = dynamic_cast<TransientTransportProblem *>( this->emodelList.at( 0 ).get() );
    EngngModel *em_sm                = this->emodelList.at( 1 ).get();

    // Activate TM model
    {
        std::cout << "#3d printing TM elements " << printer.getGrid().active_elements() << std::endl;
        std::cout << "#3d printing TM nodes " << printer.getGrid().active_nodes() << std::endl;

        // Resize DOFs and Element lists
        em_tm->giveDomain( 1 )->resizeDofManagers( printer.getGrid().active_nodes() );
        em_tm->giveDomain( 1 )->resizeElements( printer.getGrid().active_elements() );

        for ( auto &n : this->printer.getGrid().giveNodes() ) {
            add_node_if_not_exists2( em_tm, n.second );
        }

        for ( auto &n : this->printer.getGrid().giveVoxels() ) {
            add_element_if_not_exists2( em_tm, n.second );
        }

        if ( em_tm->giveExportModuleManager()->giveNumberOfModules() > 0 ) {
            em_tm->giveExportModuleManager()->giveModule( 1 )->initialize();
        }
    }

    if ( !this->skipSM ) {
        std::cout << "#3d printing SM elements " << printer.getGrid().active_elements() << std::endl;
        std::cout << "#3d printing SM nodes " << printer.getGrid().active_nodes() << std::endl;

        // Resize DOFs and Element lists
        em_sm->giveDomain( 1 )->resizeDofManagers( printer.getGrid().active_nodes() );
        em_sm->giveDomain( 1 )->resizeElements( printer.getGrid().active_elements() );

        for ( auto &n : this->printer.getGrid().giveNodes() ) {
            add_sm_node_if_not_exists2( em_sm, n.second );
        }

        for ( auto &n : this->printer.getGrid().giveVoxels() ) {
            add_sm_element_if_not_exists2( em_sm, n.second );
        }

        if ( em_sm->giveExportModuleManager()->giveNumberOfModules() > 0 ) {
            em_sm->giveExportModuleManager()->giveModule( 1 )->initialize();
        }
    }

    em_tm->giveDomain( 1 )->giveConnectivityTable()->reset();

    if ( sp->giveCurrentStep() ) {
        smstep = sp->giveCurrentStep()->giveMetaStepNumber();
        sjstep = sp->giveMetaStep( smstep )->giveStepRelativeNumber( sp->giveCurrentStep()->giveNumber() ) + 1;
    }

    for ( int imstep = smstep; imstep <= sp->giveNumberOfMetaSteps(); imstep++ ) { // loop over meta steps
        MetaStep *activeMStep = sp->giveMetaStep( imstep );
        // update state according to new meta step in all slaves
        this->initMetaStepAttributes( activeMStep );

        int nTimeSteps = activeMStep->giveNumberOfSteps();
        for ( int jstep = sjstep; jstep <= nTimeSteps; jstep++ ) { // loop over time steps
            this->timer.startTimer( EngngModelTimer ::EMTT_SolutionStepTimer );
            this->timer.initTimer( EngngModelTimer ::EMTT_NetComputationalStepTimer );
            sp->preInitializeNextStep();
            sp->giveNextStep();

            // renumber equations if necessary. Ensure to call forceEquationNumbering() for staggered problems
            if ( this->requiresEquationRenumbering( sp->giveCurrentStep() ) ) {
                this->forceEquationNumbering();
            }

            Timer timer;
            timer.startTimer();
            this->initializeYourself( sp->giveCurrentStep() );
            timer.stopTimer();

            timer.startTimer();
            this->solveYourselfAt( sp->giveCurrentStep() );
            timer.stopTimer();

            timer.startTimer();
            this->updateYourself( sp->giveCurrentStep() );
            timer.stopTimer();

            timer.startTimer();
            this->terminate( sp->giveCurrentStep() );
            timer.stopTimer();

            this->timer.stopTimer( EngngModelTimer ::EMTT_SolutionStepTimer );
            double _steptime = this->timer.getUtime( EngngModelTimer ::EMTT_SolutionStepTimer );
            OOFEM_LOG_INFO( "EngngModel info: user time consumed by solution step %d: %.2fs\n",
                sp->giveCurrentStep()->giveNumber(), _steptime );

            if ( !suppressOutput ) {
                fprintf( this->giveOutputStream(), "\nUser time consumed by solution step %d: %.3f [s]\n\n",
                    sp->giveCurrentStep()->giveNumber(), _steptime );
            }

#ifdef __MPI_PARALLEL_MODE
            if ( loadBalancingFlag ) {
                this->balanceLoad( sp->giveCurrentStep() );
            }

#endif

            if ( ( sp->giveCurrentStep()->giveTargetTime() >= this->endOfTimeOfInterest ) && this->adaptiveStepLength ) {
                break;
            }
        }
    }
}

void AdditiveManufacturingProblem ::solveYourselfAt( TimeStep *tStep )
{
#ifdef VERBOSE
    OOFEM_LOG_RELEVANT( "Solving [step number %5d, time %e]\n", tStep->giveNumber(), tStep->giveTargetTime() );
#endif

    for ( auto &emodel : emodelList ) {
        emodel->solveYourselfAt( tStep );
    }

    tStep->incrementStateCounter();
}

int AdditiveManufacturingProblem ::forceEquationNumbering()
{
    int neqs = 0;
    for ( auto &emodel : emodelList ) {
        // renumber equations if necessary
        if ( emodel->requiresEquationRenumbering( emodel->giveCurrentStep() ) ) {
            neqs += emodel->forceEquationNumbering();
        }
    }

    return neqs;
}

void AdditiveManufacturingProblem ::updateYourself( TimeStep *tStep )
{
    if ( adaptiveStepLength ) {
        this->prevStepLength = this->currentStepLength;
    }

    for ( auto &emodel : emodelList ) {
        emodel->updateYourself( tStep );
    }

    EngngModel ::updateYourself( tStep );
}

void AdditiveManufacturingProblem ::terminate( TimeStep *tStep )
{
    for ( auto &emodel : emodelList ) {
        emodel->terminate( tStep );
    }
}

void AdditiveManufacturingProblem ::doStepOutput( TimeStep *tStep )
{
    for ( auto &emodel : emodelList ) {
        emodel->giveExportModuleManager()->doOutput( tStep );
    }
}


void AdditiveManufacturingProblem ::printOutputAt( FILE *file, TimeStep *tStep )
{
    // Subproblems handle the output themselves.
}


void AdditiveManufacturingProblem ::saveContext( DataStream &stream, ContextMode mode )
{
    EngngModel ::saveContext( stream, mode );
    for ( auto &emodel : emodelList ) {
        emodel->saveContext( stream, mode );
    }
}


void AdditiveManufacturingProblem ::restoreContext( DataStream &stream, ContextMode mode )
{
    EngngModel ::restoreContext( stream, mode );
    for ( auto &emodel : this->emodelList ) {
        emodel->restoreContext( stream, mode );
    }
}


EngngModel *
AdditiveManufacturingProblem ::giveSlaveProblem( int i )
{
    if ( ( i > 0 ) && ( i <= this->giveNumberOfSlaveProblems() ) ) {
        return this->emodelList[i - 1].get();
    } else {
        OOFEM_ERROR( "Undefined problem" );
    }

    return NULL;
}


int AdditiveManufacturingProblem ::checkProblemConsistency()
{
    // check internal consistency
    // if success returns nonzero
    int result = 1;
    for ( auto &emodel : emodelList ) {
        result &= emodel->checkProblemConsistency();
    }

#ifdef VERBOSE
    if ( result ) {
        OOFEM_LOG_DEBUG( "Consistency check:  OK\n" );
    } else {
        VERBOSE_PRINTS( "Consistency check", "failed" )
        exit( 1 );
    }

#endif

    return result;
}

void AdditiveManufacturingProblem ::updateDomainLinks()
{
    for ( auto &emodel : emodelList ) {
        emodel->updateDomainLinks();
    }
}

void AdditiveManufacturingProblem ::setRenumberFlag()
{
    for ( auto &emodel : emodelList ) {
        emodel->setRenumberFlag();
    }
}

#ifdef __OOFEG
void AdditiveManufacturingProblem ::drawYourself( oofegGraphicContext &gc )
{
    int ap = gc.getActiveProblemIndx();
    if ( ( ap > 0 ) && ( ap <= giveNumberOfSlaveProblems() ) ) {
        this->giveSlaveProblem( ap )->drawYourself( gc );
    }
}

void AdditiveManufacturingProblem ::drawElements( oofegGraphicContext &gc )
{
    int ap = gc.getActiveProblemIndx();
    if ( ( ap > 0 ) && ( ap <= giveNumberOfSlaveProblems() ) ) {
        this->giveSlaveProblem( ap )->drawElements( gc );
    }
}

void AdditiveManufacturingProblem ::drawNodes( oofegGraphicContext &gc )
{
    int ap = gc.getActiveProblemIndx();
    if ( ( ap > 0 ) && ( ap <= giveNumberOfSlaveProblems() ) ) {
        this->giveSlaveProblem( ap )->drawNodes( gc );
    }
}
#endif
} // end namespace oofem
