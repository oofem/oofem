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

// Milan ?????????????????
//#include "gpinitmodule.h"
// Milan ?????????????????

#include "nummet.h"
#include "sparsemtrx.h"
#include "engngm.h"
#include "timestep.h"
#include "metastep.h"
#include "element.h"
#include "set.h"
#include "load.h"
#include "bodyload.h"
#include "boundaryload.h"
#include "nodalload.h"
#include "oofemcfg.h"
#include "timer.h"
#include "dofmanager.h"
#include "node.h"
#include "activebc.h"
#include "timestep.h"
#include "verbose.h"
#include "datastream.h"
#include "oofemtxtdatareader.h"
#include "sloangraph.h"
#include "logger.h"
#include "errorestimator.h"
#include "contextioerr.h"
#include "outputmanager.h"
#include "exportmodulemanager.h"
#include "initmodulemanager.h"
#include "feinterpol3d.h"
#include "classfactory.h"
#include "xfem/xfemmanager.h"
#include "parallelcontext.h"
#include "unknownnumberingscheme.h"
#include "contact/contactmanager.h"


#ifdef __PARALLEL_MODE
 #include "problemcomm.h"
 #include "processcomm.h"
 #include "loadbalancer.h"
#endif

#include <cstdio>
#include <cstdarg>
#include <ctime>
#ifdef _OPENMP
    #include <omp.h>
#endif
#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
EngngModel :: EngngModel(int i, EngngModel *_master) : domainNeqs(), domainPrescribedNeqs(),
    exportModuleManager(this),
    initModuleManager(this),
    monitorManager(this)
{
    suppressOutput = false;

    number = i;
    numberOfSteps = 0;
    numberOfEquations = 0;
    numberOfPrescribedEquations = 0;
    renumberFlag = false;
    equationNumberingCompleted = 0;
    ndomains = 0;
    nMetaSteps = 0;
    profileOpt = false;
    nonLinFormulation = UNKNOWN;

    outputStream          = NULL;

    referenceFileName     = "";

    contextOutputMode     = COM_NoContext;
    contextOutputStep     = 0;
    pMode                 = _processor;  // for giveContextFile()
    pScale                = macroScale;

    master                = _master; // master mode by default
    // create context if in master mode; otherwise request context from master
    if ( master ) {
        context = master->giveContext();
    } else {
        context = new EngngModelContext();
    }

    parallelFlag = 0;
    numProcs = 1;
    rank = 0;
    nonlocalExt = 0;
#ifdef __PARALLEL_MODE
    loadBalancingFlag = false;
    force_load_rebalance_in_first_step = false;
    lb = NULL;
    lbm = NULL;
    communicator = NULL;
    nonlocCommunicator = NULL;
    commBuff = NULL;
 #ifdef __USE_MPI
    comm = MPI_COMM_SELF;
 #endif
#endif
}


EngngModel :: ~EngngModel()
{
    // master deletes the context
    if ( master == NULL ) {
        delete context;
    }

    //fclose (inputStream) ;
    if ( outputStream ) {
        fclose(outputStream);
    }

#ifdef __PARALLEL_MODE
    delete communicator;
    delete nonlocCommunicator;
    delete commBuff;
#endif
}


void EngngModel :: setParallelMode(bool newParallelFlag)
{
    parallelFlag = newParallelFlag;
    if ( parallelFlag ) {
        initParallel();
    }
}


void
EngngModel :: Instanciate_init()
{
    // create domains
    domainNeqs.clear();
    domainNeqs.resize(ndomains);
    domainPrescribedNeqs.clear();
    domainPrescribedNeqs.resize(ndomains);
    domainList.clear();
    domainList.reserve(ndomains);
    for ( int i = 1; i <= ndomains; i++ ) {
        domainList.push_back(std::make_unique<Domain>(i, 0, this));
    }

    this->initParallelContexts();
}


int EngngModel :: instanciateYourself(DataReader &dr, InputRecord &ir, const char *dataOutputFileName, const char *desc)
{
    referenceFileName = dr.giveReferenceName();

    bool inputReaderFinish = true;

    this->coreOutputFileName = std :: string(dataOutputFileName);
    this->dataOutputFileName = std :: string(dataOutputFileName);

    if ( this->giveProblemMode() == _postProcessor ) {
        // modify output file name to prevent output to be lost
        this->dataOutputFileName.append(".oofeg");
    }


    this->Instanciate_init(); // Must be done after initializeFrom


    this->startTime = time(NULL);

#  ifdef VERBOSE
    OOFEM_LOG_DEBUG( "Reading all data from \"%s\"\n", referenceFileName.c_str() );
#  endif

    simulationDescription = std::string(desc);

    try {
        // instanciate receiver
        this->initializeFrom(ir);
        exportModuleManager.initializeFrom(ir);
        initModuleManager.initializeFrom(ir);
        monitorManager.initializeFrom(ir);

        if ( this->nMetaSteps == 0 ) {
            inputReaderFinish = false;
            this->instanciateDefaultMetaStep(ir);
        } else {
            this->instanciateMetaSteps(dr);
        }

        // instanciate initialization module manager
        initModuleManager.instanciateYourself(dr, ir);
        // instanciate export module manager
        exportModuleManager.instanciateYourself(dr, ir);
        // instanciate monitor manager
        monitorManager.instanciateYourself(dr, ir);
        this->instanciateDomains(dr);

        exportModuleManager.initialize();

        // Milan ??????????????????
        //GPImportModule* gim = new GPImportModule(this);
        //gim -> getInput();
        // Milan ??????????????????

        // check emodel input record if no default metastep, since all has been read
        if ( inputReaderFinish ) {
            ir.finish();
        }
    } catch ( InputException &e ) {
        OOFEM_ERROR("Error initializing from user input: %s\n", e.what());
    }

    return 1;
}


void
EngngModel :: initializeFrom(InputRecord &ir)
{
    IR_GIVE_FIELD(ir, numberOfSteps, _IFT_EngngModel_nsteps);
    if ( numberOfSteps <= 0 ) {
        OOFEM_ERROR("nsteps not specified, bad format");
    }

    contextOutputStep =  0;
    IR_GIVE_OPTIONAL_FIELD(ir, contextOutputStep, _IFT_EngngModel_contextoutputstep);
    if ( contextOutputStep ) {
        this->setUDContextOutputMode(contextOutputStep);
    }

    renumberFlag = false;
    IR_GIVE_OPTIONAL_FIELD(ir, renumberFlag, _IFT_EngngModel_renumberFlag);
    profileOpt = false;
    IR_GIVE_OPTIONAL_FIELD(ir, profileOpt, _IFT_EngngModel_profileOpt);
    nMetaSteps   = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, nMetaSteps, _IFT_EngngModel_nmsteps);
    int _val = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, _IFT_EngngModel_nonLinFormulation);
    nonLinFormulation = ( fMode ) _val;

    int eeTypeId = -1;
    IR_GIVE_OPTIONAL_FIELD(ir, eeTypeId, _IFT_EngngModel_eetype);
    if ( eeTypeId >= 0 ) {
        this->defaultErrEstimator = classFactory.createErrorEstimator( ( ErrorEstimatorType ) eeTypeId, 1, this->giveDomain(1) );
        this->defaultErrEstimator->initializeFrom(ir);
    }

    IR_GIVE_OPTIONAL_FIELD(ir, parallelFlag, _IFT_EngngModel_parallelflag);
    // fprintf (stderr, "Parallel mode is %d\n", parallelFlag);

#ifdef __PARALLEL_MODE
    /* Load balancing support */
    _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, _IFT_EngngModel_loadBalancingFlag);
    loadBalancingFlag = _val;

    _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, _IFT_EngngModel_forceloadBalancingFlag);
    force_load_rebalance_in_first_step = _val;

#endif

    suppressOutput = ir.hasField(_IFT_EngngModel_suppressOutput);

    if ( suppressOutput ) {
        //printf("Suppressing output.\n");
    }
    else {

        if ( ( outputStream = fopen(this->dataOutputFileName.c_str(), "w") ) == NULL ) {
            OOFEM_ERROR("Can't open output file %s", this->dataOutputFileName.c_str());
        }

        fprintf(outputStream, "%s", PRG_HEADER);
        fprintf(outputStream, "\nStarting analysis on: %s\n", ctime(& this->startTime) );
        fprintf(outputStream, "%s\n", simulationDescription.c_str());

#ifdef __PARALLEL_MODE
        if ( this->isParallel() ) {
            fprintf(outputStream, "Problem rank is %d/%d on %s\n\n", this->rank, this->numProcs, this->processor_name);
        }
#endif
    }
}


int
EngngModel :: instanciateDomains(DataReader &dr)
{
    int result = 1;
    // read problem domains
    for ( auto &domain: domainList ) {
        result &= domain->instanciateYourself(dr);
    }
    this->postInitialize();

    return result;
}


int
EngngModel :: instanciateMetaSteps(DataReader &dr)
{
    // create meta steps
    metaStepList.clear();
    metaStepList.reserve(nMetaSteps);
    for ( int i = 1; i <= this->nMetaSteps; i++ ) {
        //MetaStep *mstep = new MetaStep(i, this);
        metaStepList.emplace_back(i, this);
    }

    // read problem domains
    for ( int i = 1; i <= this->nMetaSteps; i++ ) {
        auto &ir = dr.giveInputRecord(DataReader :: IR_mstepRec, i);
        metaStepList[i-1].initializeFrom(ir);
    }

    this->numberOfSteps = metaStepList.size();
    return 1;
}


int
EngngModel :: instanciateDefaultMetaStep(InputRecord &ir)
{
    if ( numberOfSteps == 0 ) {
        OOFEM_ERROR("nsteps cannot be zero");
    }

    // create default meta steps
    this->nMetaSteps = 1;
    metaStepList.clear();
    //MetaStep *mstep = new MetaStep(1, this, numberOfSteps, *ir);
    metaStepList.emplace_back(1, this, numberOfSteps, ir);

    return 1;
}


int
EngngModel :: giveNumberOfDomainEquations(int id, const UnknownNumberingScheme &num)
{
    //
    // returns number of equations of current problem
    // this method is implemented here, because some method may add some
    // conditions into the system and this may results into increased number of
    // equations.
    //
    if ( !equationNumberingCompleted ) {
        this->forceEquationNumbering();
    }

    return num.isDefault() ? domainNeqs.at(id) : domainPrescribedNeqs.at(id);
}


int
EngngModel :: forceEquationNumbering(int id)
{
    // forces equation renumbering for current time step
    // intended mainly for problems with changes of static system
    // during solution
    // OUTPUT:
    // sets this->numberOfEquations and this->numberOfPrescribedEquations and returns this value

    Domain *domain = this->giveDomain(id);
    TimeStep *currStep = this->giveCurrentStep();

    this->domainNeqs.at(id) = 0;
    this->domainPrescribedNeqs.at(id) = 0;

    if ( !this->profileOpt ) {
        for ( auto &node : domain->giveDofManagers() ) {
            node->askNewEquationNumbers(currStep);
        }

        for ( auto &elem : domain->giveElements() ) {
            int nnodes = elem->giveNumberOfInternalDofManagers();
            for ( int k = 1; k <= nnodes; k++ ) {
                elem->giveInternalDofManager(k)->askNewEquationNumbers(currStep);
            }
        }

        // For special boundary conditions;
        for ( auto &bc : domain->giveBcs() ) {
            int nnodes = bc->giveNumberOfInternalDofManagers();
            for ( int k = 1; k <= nnodes; k++ ) {
                bc->giveInternalDofManager(k)->askNewEquationNumbers(currStep);
            }
        }
    } else {
        // invoke profile reduction
        int initialProfile, optimalProfile;
        Timer timer;
        OOFEM_LOG_INFO("\nRenumbering DOFs with Sloan's algorithm...\n");
        timer.startTimer();

        SloanGraph graph(domain);
        graph.initialize();
        graph.tryParameters(0, 0);
        initialProfile = graph.giveOptimalProfileSize();
        graph.tryParameters(2, 1);
        graph.tryParameters(1, 0);
        graph.tryParameters(5, 1);
        graph.tryParameters(10, 1);
        optimalProfile = graph.giveOptimalProfileSize();

        timer.stopTimer();

        OOFEM_LOG_DEBUG( "Sloan's algorithm done in %.2fs\n", timer.getUtime() );
        OOFEM_LOG_DEBUG("Nominal profile %d (old) %d (new)\n", initialProfile, optimalProfile);

        //FILE* renTableFile = fopen ("rentab.dat","w");
        //graph.writeOptimalRenumberingTable (renTableFile);
        graph.askNewOptimalNumbering(currStep);
    }

    return domainNeqs.at(id);
}


int
EngngModel :: forceEquationNumbering()
{
    // set numberOfEquations counter to zero
    this->numberOfEquations = 0;
    this->numberOfPrescribedEquations = 0;

    OOFEM_LOG_DEBUG("Renumbering dofs in all domains\n");
    for ( int i = 1; i <= this->giveNumberOfDomains(); i++ ) {
        domainNeqs.at(i) = 0;
        this->numberOfEquations += this->forceEquationNumbering(i);
    }

    equationNumberingCompleted = 1;

    for ( int i = 1; i <= this->giveNumberOfDomains(); i++ ) {
        this->numberOfPrescribedEquations += domainPrescribedNeqs.at(i);
    }

    for ( std :: size_t i = 1; i <= parallelContextList.size(); i++ ) {
        this->parallelContextList[i-1].init((int)i);
    }


    return this->numberOfEquations;
}


void
EngngModel :: solveYourself()
{
    int smstep = 1, sjstep = 1;

    this->timer.startTimer(EngngModelTimer :: EMTT_AnalysisTimer);

    if ( this->currentStep ) {
        smstep = this->currentStep->giveMetaStepNumber();
        sjstep = this->giveMetaStep(smstep)->giveStepRelativeNumber( this->currentStep->giveNumber() ) + 1;
    }

    for ( int imstep = smstep; imstep <= nMetaSteps; imstep++, sjstep = 1 ) { //loop over meta steps
        auto activeMStep = this->giveMetaStep(imstep);
        // update state according to new meta step
        this->initMetaStepAttributes(activeMStep);

	for ( int jstep = sjstep; jstep <= activeMStep->giveNumberOfSteps(); jstep++ ) { //loop over time steps
            this->timer.startTimer(EngngModelTimer :: EMTT_SolutionStepTimer);
            this->timer.initTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);

            this->preInitializeNextStep();
            this->giveNextStep();

            // renumber equations if necessary. Ensure to call forceEquationNumbering() for staggered problems
            if ( this->requiresEquationRenumbering( this->giveCurrentStep() ) ) {
                this->forceEquationNumbering();
            }

            OOFEM_LOG_DEBUG("Number of equations %d\n", this->giveNumberOfDomainEquations( 1, EModelDefaultEquationNumbering()) );

            this->initializeYourself( this->giveCurrentStep() );
            this->solveYourselfAt( this->giveCurrentStep() );
            this->updateYourself( this->giveCurrentStep() );

            this->timer.stopTimer(EngngModelTimer :: EMTT_SolutionStepTimer);

            this->terminate( this->giveCurrentStep() );

            double _steptime = this->giveSolutionStepTime();
            OOFEM_LOG_INFO("EngngModel info: user time consumed by solution step %d: %.2fs\n",
                           this->giveCurrentStep()->giveNumber(), _steptime);

            if ( !suppressOutput ) {
                fprintf(this->giveOutputStream(), "\nUser time consumed by solution step %d: %.3f [s]\n\n",
                        this->giveCurrentStep()->giveNumber(), _steptime);
            }

#ifdef __PARALLEL_MODE
            if ( loadBalancingFlag ) {
                this->balanceLoad( this->giveCurrentStep() );
            }

#endif
        }
    }
}

TimeStep* EngngModel :: generateNextStep()
{
    int smstep = 1, sjstep = 1;
    if ( this->currentStep ) {
        smstep = this->currentStep->giveMetaStepNumber();
        sjstep = this->giveMetaStep(smstep)->giveStepRelativeNumber( this->currentStep->giveNumber() ) + 1;
    }

    // test if sjstep still valid for MetaStep
    if (sjstep > this->giveMetaStep(smstep)->giveNumberOfSteps())
        smstep++;
    if (smstep > nMetaSteps) return NULL; // no more metasteps

    this->initMetaStepAttributes(this->giveMetaStep(smstep));

    this->preInitializeNextStep();
    return this->giveNextStep();
}


void
EngngModel :: initMetaStepAttributes(MetaStep *mStep)
{
    // update attributes
    this->updateAttributes(mStep); // virtual function
    // finish data acquiring
    mStep->giveAttributesRecord().finish();
}

void
EngngModel :: updateAttributes(MetaStep *mStep)
{
    MetaStep *mStep1 = this->giveMetaStep( mStep->giveNumber() ); //this line ensures correct input file in staggered problem
    auto &ir = mStep1->giveAttributesRecord();

    if ( this->giveNumericalMethod(mStep1) ) {
        this->giveNumericalMethod(mStep1)->initializeFrom(ir);
    }

#ifdef __PARALLEL_MODE
    if ( this->giveLoadBalancer() ) {
        this->giveLoadBalancer()->initializeFrom(ir);
    }

    if ( this->giveLoadBalancerMonitor() ) {
        this->giveLoadBalancerMonitor()->initializeFrom(ir);
    }

#endif
}


void
EngngModel :: updateYourself(TimeStep *tStep)
{
    for ( auto &domain: domainList ) {
#  ifdef VERBOSE
        VERBOSE_PRINT0( "Updating domain ", domain->giveNumber() )
#  endif

        for ( auto &dman : domain->giveDofManagers() ) {
            dman->updateYourself(tStep);
        }

        // Update xfem manager if it is present
        if ( domain->hasXfemManager() ) {
            domain->giveXfemManager()->updateYourself(tStep);
        }

#  ifdef VERBOSE
        VERBOSE_PRINT0("Updated nodes ", domain->giveNumberOfDofManagers())
#  endif


        for ( auto &elem : domain->giveElements() ) {
            // skip remote elements (these are used as mirrors of remote elements on other domains
            // when nonlocal constitutive models are used. They introduction is necessary to
            // allow local averaging on domains without fine grain communication between domains).
            if ( elem->giveParallelMode() == Element_remote ) {
                continue;
            }

            elem->updateYourself(tStep);
        }

#  ifdef VERBOSE
        VERBOSE_PRINT0("Updated Elements ", domain->giveNumberOfElements())
#  endif
    }

    // if there is an error estimator, it should be updated so that values can be exported.
    if ( this->defaultErrEstimator ) {
        this->defaultErrEstimator->estimateError(equilibratedEM, tStep);
    }
}

void
EngngModel :: terminate(TimeStep *tStep)
{
    if ( !suppressOutput ) {
        this->doStepOutput(tStep);
        fflush( this->giveOutputStream() );
    } else {
        exportModuleManager.doOutput(tStep);
    }
    monitorManager.update(tStep, Monitor::MonitorEvent::TimeStepTermination);
    
    this->saveStepContext(tStep, CM_State | CM_Definition);
}


void
EngngModel :: doStepOutput(TimeStep *tStep)
{
    if ( !suppressOutput ) {
        this->printOutputAt(this->giveOutputStream(), tStep);
        fflush( this->giveOutputStream() );
    }

    // export using export manager
    exportModuleManager.doOutput(tStep);
}

void
EngngModel :: saveStepContext(TimeStep *tStep, ContextMode mode)
{
    if ( this->giveContextOutputMode() == COM_Always || this->giveContextOutputMode() == COM_Required || 
        ( this->giveContextOutputMode() == COM_UserDefined && tStep->giveNumber() % this->giveContextOutputStep() == 0 ) ) {

        auto fname = this->giveContextFileName(this->giveCurrentStep()->giveNumber(), this->giveCurrentStep()->giveVersion());
        FileDataStream stream(fname, true);
        this->saveContext(stream, mode);
    }
}


void
EngngModel :: printOutputAt(FILE *file, TimeStep *tStep)
{
    int domCount = 0;

    for ( auto &domain: domainList ) {
        domCount += domain->giveOutputManager()->testTimeStepOutput(tStep);
    }

    if ( domCount == 0 ) {
        return;              // do not print even Solution step header
    }

    fprintf(file, "\n==============================================================");
    fprintf(file, "\nOutput for time %.8e ", tStep->giveTargetTime() * this->giveVariableScale(VST_Time) );
    fprintf(file, "\n==============================================================\n");
    for ( auto &domain: domainList ) {
        fprintf( file, "Output for domain %3d\n", domain->giveNumber() );

        domain->giveOutputManager()->doDofManOutput(file, tStep);
        domain->giveOutputManager()->doElementOutput(file, tStep);
    }
}


void
EngngModel :: printOutputAt(FILE *file, TimeStep *tStep, const IntArray &nodeSets, const IntArray &elementSets)
{
    for ( auto &domain: domainList ) {
        int dnum = domain->giveNumber();
        fprintf( file, "Output for domain %3d\n", dnum );
        int nset = nodeSets.giveSize() < dnum ? 0 : nodeSets.at(dnum);
        int eset = elementSets.giveSize() < dnum ? 0 : elementSets.at(dnum);

        this->outputNodes(file, *domain, tStep, nset);
        this->outputElements(file, *domain, tStep, eset);
        ///@todo Add general support for reaction forces
#if 0
        this->outputReactionForces(file, *domain, tStep, nset);
#endif
    }
}


void
EngngModel :: outputNodes(FILE *file, Domain &domain, TimeStep *tStep, int setNum)
{
    fprintf(file, "\n\nNode output:\n------------------\n");

    if ( setNum == 0 ) { // No set specified, export all
        for ( auto &dman : domain.giveDofManagers() ) {
            if ( dman->giveParallelMode() == DofManager_null ) {
                continue;
            }
            dman->printOutputAt(file, tStep);
        }
    } else {
        auto &nodes = domain.giveSet(setNum)->giveNodeList();

        for ( int inode : nodes ) {
            auto dman = domain.giveDofManager(inode);
            if ( dman->giveParallelMode() == DofManager_null ) {
                continue;
            }
            dman->printOutputAt(file, tStep);
        }
    }
    fprintf(file, "\n\n");
}


void
EngngModel :: outputElements(FILE *file, Domain &domain, TimeStep *tStep, int setNum)
{
    fprintf(file, "\n\nElement output:\n---------------\n");

    if ( setNum == 0 ) {
        for ( auto &elem : domain.giveElements() ) {
            if ( elem->giveParallelMode() == Element_remote ) {
                continue;
            }
            elem->printOutputAt(file, tStep);
        }
    } else {
        auto &elements = domain.giveSet(setNum)->giveElementList();
        for ( int ielem : elements ) {
            auto element = domain.giveElement(ielem);
            if ( element->giveParallelMode() == Element_remote ) {
                continue;
            }
            element->printOutputAt(file, tStep);
        }
    }
    fprintf(file, "\n\n");
}


void EngngModel :: printYourself()
{
    printf("\nEngineeringModel: instance %s\n", this->giveClassName() );
    printf("number of steps: %d\n", this->giveNumberOfSteps() );
    printf("number of eq's : %d\n", numberOfEquations);
}

void EngngModel :: printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep)
{
    iDof->printSingleOutputAt(stream, tStep, 'd', VM_Total);
}

void EngngModel :: assemble(SparseMtrx &answer, TimeStep *tStep, const MatrixAssembler &ma,
                            const UnknownNumberingScheme &s, Domain *domain)
{
    IntArray loc;
    FloatMatrix mat, R;
#ifdef _OPENMP
    omp_lock_t writelock;
    omp_init_lock(&writelock);
#endif

    this->timer.resumeTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
    int nelem = domain->giveNumberOfElements();
#ifdef _OPENMP
#pragma omp parallel for shared(answer) private(mat, R, loc)
#endif
    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
        auto element = domain->giveElement(ielem);
        // skip remote elements (these are used as mirrors of remote elements on other domains
        // when nonlocal constitutive models are used. They introduction is necessary to
        // allow local averaging on domains without fine grain communication between domains).
        if ( element->giveParallelMode() == Element_remote || !element->isActivated(tStep) || !this->isElementActivated(element) ) {
            continue;
        }

        ma.matrixFromElement(mat, *element, tStep);

        if ( mat.isNotEmpty() ) {
            ma.locationFromElement(loc, *element, s);
            ///@todo This rotation matrix is not flexible enough.. it can only work with full size matrices and doesn't allow for flexibility in the matrixassembler.
            if ( element->giveRotationMatrix(R) ) {
                mat.rotatedWith(R);
            }

#ifdef _OPENMP
 #pragma omp critical
#endif
            if ( answer.assemble(loc, mat) == 0 ) {
                OOFEM_ERROR("sparse matrix assemble error");
            }
        }
    }

#ifdef _OPENMP
#pragma omp parallel for shared(answer) private(mat, R, loc)
#endif
    for ( auto &bc : domain->giveBcs() ) {
        auto abc = dynamic_cast< ActiveBoundaryCondition * >(bc.get());

        if ( abc ) {
            /// @note: Some active bcs still make changes even when they are not applied
            /// We should probably reconsider this approach, so that they e.g. just prescribe their lagrange mult. instead.
#ifdef _OPENMP
            ma.assembleFromActiveBC(answer, *abc, tStep, s, s, &writelock);
#else
            ma.assembleFromActiveBC(answer, *abc, tStep, s, s);
#endif
        } else if ( bc->giveSetNumber() ) {
            if ( !bc->isImposed(tStep) ) continue;
            auto load = dynamic_cast< Load * >(bc.get());
            if ( !load ) continue;
            // Now we assemble the corresponding load type for the respective components in the set:
            IntArray loc, bNodes;
            FloatMatrix mat, R;
            BodyLoad *bodyLoad;
            SurfaceLoad* sLoad;
            EdgeLoad* eLoad;
            Set *set = domain->giveSet( bc->giveSetNumber() );

            if ( ( bodyLoad = dynamic_cast< BodyLoad * >(load) ) ) { // Body load:
                const IntArray &elements = set->giveElementList();
                for ( auto ielem : elements ) {
                    auto element = domain->giveElement( ielem );
                    mat.clear();
                    ma.matrixFromLoad(mat, *element, bodyLoad, tStep);

                    if ( mat.isNotEmpty() ) {
                        if ( element->giveRotationMatrix(R) ) {
                            mat.rotatedWith(R);
                        }

                        ma.locationFromElement(loc, *element, s);
#ifdef _OPENMP
            			omp_set_lock(&writelock);
#endif
                        answer.assemble(loc, mat);
#ifdef _OPENMP
			            omp_unset_lock(&writelock);
#endif
                    }
                }
            } else if ( ( sLoad = dynamic_cast< SurfaceLoad * >(load) ) ) {
                const auto &surfaces = set->giveBoundaryList();
                for ( int ibnd = 1; ibnd <= surfaces.giveSize() / 2; ++ibnd ) {
                    auto element = domain->giveElement( surfaces.at(ibnd * 2 - 1) );
                    int boundary = surfaces.at(ibnd * 2);
                    mat.clear();
                    ma.matrixFromSurfaceLoad(mat, *element, sLoad, boundary, tStep);

                    if ( mat.isNotEmpty() ) {
                        bNodes = element->giveInterpolation()->boundaryGiveNodes(boundary);
                        if ( element->computeDofTransformationMatrix(R, bNodes, false) ) {
                            mat.rotatedWith(R);
                        }

                        ma.locationFromElementNodes(loc, *element, bNodes, s);
 
 #ifdef _OPENMP
            			omp_set_lock(&writelock);
#endif			
                        answer.assemble(loc, mat);
#ifdef _OPENMP
			            omp_unset_lock(&writelock);
#endif
                    }
                }
            } else if ( ( eLoad = dynamic_cast< EdgeLoad * >(load) ) ) {
                const auto &edges = set->giveEdgeList();
                for ( int ibnd = 1; ibnd <= edges.giveSize() / 2; ++ibnd ) {
                    auto element = domain->giveElement( edges.at(ibnd * 2 - 1) );
                    int boundary = edges.at(ibnd * 2);
                    mat.clear();
                    ma.matrixFromEdgeLoad(mat, *element, eLoad, boundary, tStep);

                    if ( mat.isNotEmpty() ) {
                        bNodes = element->giveInterpolation()->boundaryEdgeGiveNodes(boundary);
                        if ( element->computeDofTransformationMatrix(R, bNodes, false) ) {
                            mat.rotatedWith(R);
                        }

                        ma.locationFromElementNodes(loc, *element, bNodes, s);
#ifdef _OPENMP
            			omp_set_lock(&writelock);
#endif			
                        answer.assemble(loc, mat);
#ifdef _OPENMP
			            omp_unset_lock(&writelock);
#endif
                    }
                }
            }
        }
    }

    if ( domain->hasContactManager() ) {
        OOFEM_ERROR("Contact problems temporarily deactivated");
        //domain->giveContactManager()->assembleTangentFromContacts(answer, tStep, type, s, s);
    }

    this->timer.pauseTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);

    answer.assembleBegin();
    answer.assembleEnd();
}


void EngngModel :: assemble(SparseMtrx &answer, TimeStep *tStep, const MatrixAssembler &ma,
                            const UnknownNumberingScheme &rs, const UnknownNumberingScheme &cs,
                            Domain *domain)
// Same as assemble, but with different numbering for rows and columns
{
    IntArray r_loc, c_loc, dofids(0);
    FloatMatrix mat, R;
#ifdef _OPENMP
    omp_lock_t writelock;
    omp_init_lock(&writelock);
#endif

    this->timer.resumeTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
    int nelem = domain->giveNumberOfElements();
#ifdef _OPENMP
#pragma omp parallel for shared(answer) private(mat, R, r_loc, c_loc)
#endif
    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
        Element *element = domain->giveElement(ielem);

        if ( element->giveParallelMode() == Element_remote || !element->isActivated(tStep) || !this->isElementActivated(element) ) {
            continue;
        }

        ma.matrixFromElement(mat, *element, tStep);
        if ( mat.isNotEmpty() ) {
            ma.locationFromElement(r_loc, *element, rs);
            ma.locationFromElement(c_loc, *element, cs);
            // Rotate it
            ///@todo This rotation matrix is not flexible enough.. it can only work with full size matrices and doesn't allow for flexibility in the matrixassembler.
            if ( element->giveRotationMatrix(R) ) {
                mat.rotatedWith(R);
            }

#ifdef _OPENMP
 #pragma omp critical
#endif
            if ( answer.assemble(r_loc, c_loc, mat) == 0 ) {
                OOFEM_ERROR("sparse matrix assemble error");
            }
        }
    }

#ifdef _OPENMP
#pragma omp parallel for shared(answer) private(mat, R, r_loc, c_loc)
#endif
    for ( auto &gbc : domain->giveBcs() ) {
        ActiveBoundaryCondition *bc = dynamic_cast< ActiveBoundaryCondition * >( gbc.get() );
        if ( bc != NULL ) {
#ifdef _OPENMP
            ma.assembleFromActiveBC(answer, *bc, tStep, rs, cs, &writelock);
#else
            ma.assembleFromActiveBC(answer, *bc, tStep, rs, cs);
#endif
        }
    }

    if ( domain->hasContactManager() ) {
        OOFEM_ERROR("Contant problems temporarily deactivated");
        //domain->giveContactManager()->assembleTangentFromContacts(answer, tStep, type, rs, cs);
    }

    this->timer.pauseTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);

    answer.assembleBegin();
    answer.assembleEnd();
}


void EngngModel :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                  const VectorAssembler &va, ValueModeType mode,
                                  const UnknownNumberingScheme &s, Domain *domain, FloatArray *eNorms)
{
    if ( eNorms ) {
        int maxdofids = domain->giveMaxDofID();
#ifdef __PARALLEL_MODE
        if ( this->isParallel() ) {
            int val;
            MPI_Allreduce(& maxdofids, & val, 1, MPI_INT, MPI_MAX, this->comm);
            maxdofids = val;
        }
#endif
        eNorms->resize(maxdofids);
        eNorms->zero();
    }

    this->assembleVectorFromDofManagers(answer, tStep, va, mode, s, domain, eNorms);
    this->assembleVectorFromElements(answer, tStep, va, mode, s, domain, eNorms);
    this->assembleVectorFromBC(answer, tStep, va, mode, s, domain, eNorms);

    if ( this->isParallel() ) {
        if ( eNorms ) {
            FloatArray localENorms = * eNorms;
            this->giveParallelContext(domain->giveNumber())->accumulate(localENorms, *eNorms);
        }
    }
}


void EngngModel :: assembleVectorFromDofManagers(FloatArray &answer, TimeStep *tStep, const VectorAssembler &va, ValueModeType mode,
                                                 const UnknownNumberingScheme &s, Domain *domain, FloatArray *eNorms)
{
    ///@todo This should be removed when loads are given through sets.
    IntArray loc, dofids;
    FloatArray charVec;
    FloatMatrix R;
    IntArray dofIDarry;
    int nnode = domain->giveNumberOfDofManagers();

    this->timer.resumeTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
    // Note! For normal master dofs, loc is unique to each node, but there can be slave dofs, so we must keep it shared, unfortunately.
#ifdef _OPENMP
#pragma omp parallel for shared(answer, eNorms) private(R, charVec, loc, dofids)
#endif
    for ( int i = 1; i <= nnode; i++ ) {
        DofManager *node = domain->giveDofManager(i);

        charVec.clear();
        for ( int iload : *node->giveLoadArray() ) {   // to more than one load
            Load *load = domain->giveLoad(iload);
            va.vectorFromNodeLoad(charVec, *node, static_cast< NodalLoad* >(load), tStep, mode);

            if ( node->giveParallelMode() == DofManager_shared ) {
                charVec.times( 1. / ( node->givePartitionsConnectivitySize() ) );
            }

            if ( charVec.isNotEmpty() ) {
                if ( node->computeM2LTransformation(R, dofIDarry) ) {
                    charVec.rotatedWith(R, 't');
                }

                if ( load->giveDofIDs().giveSize() ) {
                    node->giveLocationArray(load->giveDofIDs(), loc, s);
                } else {
                    node->giveCompleteLocationArray(loc, s);
                }
#ifdef _OPENMP
 #pragma omp critical
#endif
                {
                    answer.assemble(charVec, loc);
                    if ( eNorms ) {
                        node->giveCompleteMasterDofIDArray(dofids);
                        eNorms->assembleSquared(charVec, dofids);
                    }
                }
            }
        }
    }

    this->timer.pauseTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
}


void EngngModel :: assembleVectorFromBC(FloatArray &answer, TimeStep *tStep,
                                        const VectorAssembler &va, ValueModeType mode,
                                        const UnknownNumberingScheme &s, Domain *domain, FloatArray *eNorms)
{
    int nbc = domain->giveNumberOfBoundaryConditions();
#ifdef _OPENMP
    omp_lock_t writelock;
    omp_init_lock(&writelock);
#endif

    this->timer.resumeTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
#ifdef _OPENMP
#pragma omp parallel for shared(answer, eNorms)
#endif
    for ( int i = 1; i <= nbc; ++i ) {
        GeneralBoundaryCondition *bc = domain->giveBc(i);
        ActiveBoundaryCondition *abc;
        Load *load;

        if ( ( abc = dynamic_cast< ActiveBoundaryCondition * >(bc) ) ) {
#ifdef _OPENMP
            va.assembleFromActiveBC(answer, *abc, tStep, mode, s, eNorms, &writelock);
#else
            va.assembleFromActiveBC(answer, *abc, tStep, mode, s, eNorms);
#endif
        } else if ( bc->giveSetNumber() && ( load = dynamic_cast< Load * >(bc) ) && bc->isImposed(tStep) ) {
            // Now we assemble the corresponding load type fo the respective components in the set:
            IntArray dofids, loc;
            FloatArray charVec;
            FloatMatrix R;
            BodyLoad *bodyLoad;
            SurfaceLoad *sLoad;
            EdgeLoad *eLoad;
            //BoundaryEdgeLoad *eLoad;
            NodalLoad *nLoad;
            Set *set = domain->giveSet( bc->giveSetNumber() );

            if ( ( bodyLoad = dynamic_cast< BodyLoad * >(load) ) ) { // Body load:
                const IntArray &elements = set->giveElementList();
                for ( int ielem = 1; ielem <= elements.giveSize(); ++ielem ) {
                    Element *element = domain->giveElement( elements.at(ielem) );
                    if ( element->isActivated(tStep) && this->isElementActivated(element) ) {
                        charVec.clear();
                        va.vectorFromLoad(charVec, *element, bodyLoad, tStep, mode);

                        if ( charVec.isNotEmpty() ) {
                            if ( element->giveRotationMatrix(R) ) {
                                charVec.rotatedWith(R, 't');
                            }

                            va.locationFromElement(loc, *element, s, & dofids);
#ifdef _OPENMP
                			omp_set_lock(&writelock);
#endif
                            answer.assemble(charVec, loc);
                            if ( eNorms ) {
                                eNorms->assembleSquared(charVec, dofids);
                            }
#ifdef _OPENMP
            			    omp_unset_lock(&writelock);
#endif
                        }
                    }
                }
            } else if ( ( sLoad = dynamic_cast< SurfaceLoad * >(load) ) ) { // Surface load:
                const IntArray &boundaries = set->giveBoundaryList();
                for ( int ibnd = 1; ibnd <= boundaries.giveSize() / 2; ++ibnd ) {
                    Element *element = domain->giveElement( boundaries.at(ibnd * 2 - 1) );
                    if ( element->isActivated(tStep) && this->isElementActivated(element) ) {

                        int boundary = boundaries.at(ibnd * 2);
                        charVec.clear();
                        va.vectorFromSurfaceLoad(charVec, *element, sLoad, boundary, tStep, mode);

                        if ( charVec.isNotEmpty() ) {
                            //element->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);
                            auto bNodes = element->giveBoundarySurfaceNodes(boundary);
                            if ( element->computeDofTransformationMatrix(R, bNodes, false) ) {
                                charVec.rotatedWith(R, 't');
                            }

                            va.locationFromElementNodes(loc, *element, bNodes, s, & dofids);
#ifdef _OPENMP
                  			omp_set_lock(&writelock);
#endif                            
                            answer.assemble(charVec, loc);

                            if ( eNorms ) {
                                eNorms->assembleSquared(charVec, dofids);
                            }
#ifdef _OPENMP
            			    omp_unset_lock(&writelock);
#endif                            
                        }
                    }
                }
            } else if ( ( eLoad = dynamic_cast< EdgeLoad * >(load) ) ) { // Edge load:
                const IntArray &edgeBoundaries = set->giveEdgeList();
                for ( int ibnd = 1; ibnd <= edgeBoundaries.giveSize() / 2; ++ibnd ) {
                    Element *element = domain->giveElement( edgeBoundaries.at(ibnd * 2 - 1) );
                    if ( element->isActivated(tStep) && this->isElementActivated(element) ) {
                        int boundary = edgeBoundaries.at(ibnd * 2);
                        charVec.clear();
                        va.vectorFromEdgeLoad(charVec, *element, eLoad, boundary, tStep, mode);

                        if ( charVec.isNotEmpty() ) {
                            //element->giveInterpolation()->boundaryEdgeGiveNodes(bNodes, boundary);
                            auto bNodes = element->giveBoundaryEdgeNodes(boundary);
                            if ( element->computeDofTransformationMatrix(R, bNodes, false) ) {
                                charVec.rotatedWith(R, 't');
                            }

                            va.locationFromElementNodes(loc, *element, bNodes, s, & dofids);
#ifdef _OPENMP
            			    omp_set_lock(&writelock);
#endif                            
                            answer.assemble(charVec, loc);

                            if ( eNorms ) {
                                eNorms->assembleSquared(charVec, dofids);
                            }
#ifdef _OPENMP
            			    omp_unset_lock(&writelock);
#endif                            
                        }
                    }
                }
            } else if ( ( nLoad = dynamic_cast< NodalLoad * >(load) ) ) { // Nodal load:
                const IntArray &nodes = set->giveNodeList();
                for ( int idman = 1; idman <= nodes.giveSize(); ++idman ) {
                    DofManager *node = domain->giveDofManager( nodes.at(idman) );
                    charVec.clear();
                    va.vectorFromNodeLoad(charVec, *node, nLoad, tStep, mode);

                    if ( charVec.isNotEmpty() ) {
                        if ( node->computeM2LTransformation(R, nLoad->giveDofIDs()) ) {
                            charVec.rotatedWith(R, 't');
                        }

                        node->giveLocationArray(nLoad->giveDofIDs(), loc, s);
#ifdef _OPENMP
            			omp_set_lock(&writelock);
#endif                        
                        answer.assemble(charVec, loc);

                        if ( eNorms ) {
                            eNorms->assembleSquared(charVec, dofids);
                        }
#ifdef _OPENMP
            			omp_unset_lock(&writelock);
#endif                        
                    }
                }
            }
        }
    }

    this->timer.pauseTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
}


void EngngModel :: assembleVectorFromElements(FloatArray &answer, TimeStep *tStep,
                                              const VectorAssembler &va, ValueModeType mode,
                                              const UnknownNumberingScheme &s, Domain *domain, FloatArray *eNorms)
//
// for each element in domain
// and assembling every contribution to answer
//
{
    IntArray loc, dofids;
    FloatMatrix R;
    FloatArray charVec;
    int nelem = domain->giveNumberOfElements();
    bool assembleFlag = false;

    ///@todo Checking the chartype is not since there could be some other chartype in the future. We need to try and deal with chartype in a better way.
    /// For now, this is the best we can do.
    if ( this->isParallel() ) {
        // Copies internal (e.g. Gauss-Point) data from remote elements to make sure they have all information necessary for nonlocal averaging.
        this->exchangeRemoteElementData(RemoteElementExchangeTag);
    }

    this->timer.resumeTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
    ///@todo Consider using private answer variables and sum them up at the end, but it just might be slower then a shared variable.
#ifdef _OPENMP
#pragma omp parallel for shared(answer, eNorms) private(R, charVec, loc, dofids)
#endif
    for ( int i = 1; i <= nelem; i++ ) {

      Element *element = domain->giveElement(i);

        // skip remote elements (these are used as mirrors of remote elements on other domains
        // when nonlocal constitutive models are used. They introduction is necessary to
        // allow local averaging on domains without fine grain communication between domains).
        if ( element->giveParallelMode() == Element_remote ) {
            continue;
        }

        if ( !element->isActivated(tStep) || !this->isElementActivated(element) ) {
            continue;
        }

        va.vectorFromElement(charVec, *element, tStep, mode);

        if ( charVec.isNotEmpty() ) {
            if ( element->giveRotationMatrix(R) ) {
                charVec.rotatedWith(R, 't');
            }
            va.locationFromElement(loc, *element, s, & dofids);
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                answer.assemble(charVec, loc);
                if ( eNorms ) {
                    eNorms->assembleSquared(charVec, dofids);
                }
            }
        }
    }

#ifdef _OPENMP
#pragma omp parallel for shared(answer, eNorms) private(R, charVec, loc, dofids)
#endif
    for ( int i = 1; i <= nelem; i++ ) {
        Element *element = domain->giveElement(i);

        // skip remote elements (these are used as mirrors of remote elements on other domains
        // when nonlocal constitutive models are used. They introduction is necessary to
        // allow local averaging on domains without fine grain communication between domains).
        if ( element->giveParallelMode() == Element_remote ) {
            continue;
        }

        if ( !element->isActivated(tStep) || !this->isElementActivated(element) ) {
            continue;
        }

        // obtain form element its body, surface, edge, and point loads
        const IntArray& list = element->giveBodyLoadList();
        if (!list.isEmpty()) {
          for (int iload=1; iload<=list.giveSize(); iload++) { // loop over body loads
            BodyLoad *bodyLoad;
            if ((bodyLoad = dynamic_cast< BodyLoad * >(domain->giveLoad(list.at(iload))))) {
              charVec.clear();
              va.vectorFromLoad(charVec, *element, bodyLoad, tStep, mode);

              if ( charVec.isNotEmpty() ) {
                if ( element->giveRotationMatrix(R) ) {
                  charVec.rotatedWith(R, 't');
                }

                va.locationFromElement(loc, *element, s, & dofids);
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    answer.assemble(charVec, loc);            
                    if ( eNorms ) {
                        eNorms->assembleSquared(charVec, dofids);
                    }
                }
              }
            }
            
          } // loop over body load list
        } // if (!(list = element->giveBodyLoadList()).isEmpty())
    }

#ifdef _OPENMP
#pragma omp parallel for shared(answer, eNorms) private(R, charVec, loc, dofids, assembleFlag)
#endif
    for ( int i = 1; i <= nelem; i++ ) {
        Element *element = domain->giveElement(i);

        // skip remote elements (these are used as mirrors of remote elements on other domains
        // when nonlocal constitutive models are used. They introduction is necessary to
        // allow local averaging on domains without fine grain communication between domains).
        if ( element->giveParallelMode() == Element_remote ) {
            continue;
        }

        if ( !element->isActivated(tStep) || !this->isElementActivated(element) ) {
            continue;
        }

        // obtain from element its boundaryloads (surface+edge)
        const IntArray& list2 = element->giveBoundaryLoadList();

        for (int j=1; j<=list2.giveSize()/2; j++) { // loop over boundary loads
            int iload = list2.at(j * 2 - 1) ;
            int boundary = list2.at(j * 2);
            SurfaceLoad *sLoad;
            EdgeLoad *eLoad;
            assembleFlag = false;
            IntArray bNodes;

            if ((eLoad = dynamic_cast< EdgeLoad * >(domain->giveLoad(iload)))) {
                charVec.clear();
                va.vectorFromEdgeLoad(charVec, *element, eLoad, boundary, tStep, mode);
            
                if ( charVec.isNotEmpty() ) {
                    //element->giveInterpolation()->boundaryEdgeGiveNodes(bNodes, boundary);
                    bNodes = element->giveBoundaryEdgeNodes(boundary);
                    if ( element->computeDofTransformationMatrix(R, bNodes, false) ) {
                        charVec.rotatedWith(R, 't');
                    }
                    assembleFlag = true;
                }
            } else if ((sLoad = dynamic_cast< SurfaceLoad * >(domain->giveLoad(iload)))) {
                charVec.clear();
                va.vectorFromSurfaceLoad(charVec, *element, sLoad, boundary, tStep, mode);
            
                if ( charVec.isNotEmpty() ) {
                    //element->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);
                    bNodes = element->giveBoundarySurfaceNodes(boundary);
                    if ( element->computeDofTransformationMatrix(R, bNodes, false) ) {
                        charVec.rotatedWith(R, 't');
                    }
                    assembleFlag = true;
                }
            } else {
                OOFEM_ERROR ("Unsupported element boundary load type");
            }

            if ( assembleFlag ) {
                // assemble the contribution
                va.locationFromElementNodes(loc, *element, bNodes, s, & dofids);
#ifdef _OPENMP
#pragma omp critical
#endif
                {
                    answer.assemble(charVec, loc);
                    if ( eNorms ) {
                        eNorms->assembleSquared(charVec, dofids);
                    }
                }
            } // end loop over lement boundary loads
        }

    } // end loop over elements

    this->timer.pauseTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
}

void
EngngModel :: assembleExtrapolatedForces(FloatArray &answer, TimeStep *tStep, CharType type, Domain *domain)
{
    // Simply assembles contributions from each element in domain
    IntArray loc;
    FloatArray charVec, delta_u;
    FloatMatrix charMatrix, R;
    int nelems = domain->giveNumberOfElements();
    EModelDefaultEquationNumbering dn;

    answer.resize( this->giveNumberOfDomainEquations( domain->giveNumber(), EModelDefaultEquationNumbering() ) );
    answer.zero();

    this->timer.resumeTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);

#ifdef _OPENMP
//#pragma omp parallel for shared(answer) private(R, charMatrix, charVec, loc, delta_u)
#endif
    for ( int i = 1; i <= nelems; i++ ) {
        Element *element = domain->giveElement(i);

        // Skip remote elements (these are used as mirrors of remote elements on other domains
        // when nonlocal constitutive models are used. Their introduction is necessary to
        // allow local averaging on domains without fine grain communication between domains).
        if ( element->giveParallelMode() == Element_remote ) {
            continue;
        }

        if ( !element->isActivated(tStep) || !this->isElementActivated(element) ) {
            continue;
        }

            element->giveLocationArray(loc, dn);

        // Take the tangent from the previous step
        ///@todo This is not perfect. It is probably no good for viscoelastic materials, and possibly other scenarios that are rate dependent
        ///(tangent will be computed for the previous step, with whatever deltaT it had)

        element->giveCharacteristicMatrix(charMatrix, type, tStep);
        if ( charMatrix.isNotEmpty() ) {
            ///@note Temporary work-around for active b.c. used in multiscale (it can't support VM_Incremental easily).
            
#if 0
            element->computeVectorOf(VM_Incremental, tStep, delta_u);
#else
            element->computeVectorOf(VM_Total, tStep, delta_u);
            FloatArray tmp;

            if ( tStep->isTheFirstStep() ) {
                tmp = delta_u;
                tmp.zero();
            } else {
                element->computeVectorOf(VM_Total, tStep->givePreviousStep(), tmp);
            }

            delta_u.subtract(tmp);
#endif

            charVec.beProductOf(charMatrix, delta_u);
            if ( element->giveRotationMatrix(R) ) {
                charVec.rotatedWith(R, 't');
            }

            ///@todo Deal with element deactivation and reactivation properly.
#ifdef _OPENMP
#pragma omp critical
#endif
            {
                answer.assemble(charVec, loc);
            }
        }
    }


    this->timer.pauseTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
}


void
EngngModel :: assemblePrescribedExtrapolatedForces(FloatArray &answer, TimeStep *tStep, CharType type, Domain *domain)
{
    // Simply assembles contributions from each element in domain
    IntArray loc;
    FloatArray charVec, delta_u;
    FloatMatrix charMatrix, R;
    int nelems = domain->giveNumberOfElements();
    EModelDefaultEquationNumbering dn;

    answer.resize( this->giveNumberOfDomainEquations( domain->giveNumber(), EModelDefaultEquationNumbering() ) );
    answer.zero();

    this->timer.resumeTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);

#ifdef _OPENMP
//#pragma omp parallel for shared(answer) private(R, charMatrix, charVec, loc, delta_u)
#endif
    for ( int i = 1; i <= nelems; i++ ) {
        Element *element = domain->giveElement(i);

        // Skip remote elements (these are used as mirrors of remote elements on other domains
        // when nonlocal constitutive models are used. Their introduction is necessary to
        // allow local averaging on domains without fine grain communication between domains).
        if ( element->giveParallelMode() == Element_remote ) {
            continue;
        }

        if ( !element->isActivated(tStep) ) {
            continue;
        }

        element->giveLocationArray(loc, dn);

        // Take the tangent from the previous step
        ///@todo This is not perfect. It is probably no good for viscoelastic materials, and possibly other scenarios that are rate dependent
        ///(tangent will be computed for the previous step, with whatever deltaT it had)
        element->giveCharacteristicMatrix(charMatrix, type, tStep);
        element->computeVectorOfPrescribed(VM_Incremental, tStep, delta_u);
        if ( charMatrix.isNotEmpty() ) {
            charVec.beProductOf(charMatrix, delta_u);
            if ( element->giveRotationMatrix(R) ) {
                charVec.rotatedWith(R, 't');
            }

            ///@todo Deal with element deactivation and reactivation properly.
    #ifdef _OPENMP
    #pragma omp critical
    #endif
            {
                answer.assemble(charVec, loc);
            }
        }
    }

    this->timer.pauseTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
}

void
EngngModel :: assembleVectorFromContacts(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode,
                                    const UnknownNumberingScheme &s, Domain *domain, FloatArray *eNorms)
{
    if( domain->hasContactManager()) {
        domain->giveContactManager()->assembleVectorFromContacts(answer, tStep, type, mode, s, domain, eNorms);
    }
}


void
EngngModel :: updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
{
    OOFEM_ERROR("Unknown Type of component.");
}


void
EngngModel :: updateSolution(FloatArray &solutionVector, TimeStep *tStep, Domain *d)
{
    OOFEM_ERROR("updateSolution is not implemented.");
}


void
EngngModel :: updateInternalRHS(FloatArray &answer, TimeStep *tStep, Domain *d, FloatArray *eNorms)
{
    OOFEM_ERROR("updateInternalRHS is not implemented.");
}


void
EngngModel :: updateMatrix(SparseMtrx &mat, TimeStep *tStep, Domain *d)
{
    OOFEM_ERROR("updateMatrix is not implemented.");
}


void
EngngModel :: initStepIncrements()
//
// resets all temp data in elements and their gps so only
// non temp variables remain.
// this function can be used if called before this->terminate()
// for current step restart, because it causes complete reset
// of temp variables.
//
{
    for ( auto &domain: domainList ) {
        for ( auto &elem : domain->giveElements() ) {
            // skip remote elements (these are used as mirrors of remote elements on other domains
            // when nonlocal constitutive models are used. They introduction is necessary to
            // allow local averaging on domains without fine grain communication between domains).
            if ( elem->giveParallelMode() == Element_remote ) {
                continue;
            }

            elem->initForNewStep();
        }
    }
}


void
EngngModel :: updateDomainLinks()
{
    this->giveExportModuleManager()->initialize();
}


void EngngModel :: saveContext(DataStream &stream, ContextMode mode)
//
// this procedure is used mainly for two reasons:
//
// - to save context of current model in order to be able to backtrace
//   computational history (useful in nonlinear analysis, when you want
//   explore another, previously detected load path
//
// - to save context of current model in order to enable
//   possibility of post-processing.
//
// saving context means: (if needed may be enhanced)
//
// - save EngngModel state varialbles as displacement, velocity, .. vectors
// - save Elements stress, strain and material history.
//
// This version saves only Element and Material properties.
//
{
    contextIOResultType iores;
 
    if ( !stream.write(giveCurrentStep()->giveNumber()) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    giveCurrentStep()->saveContext(stream);

    if ( !stream.write(numberOfEquations) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = domainNeqs.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.write(numberOfPrescribedEquations) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = domainPrescribedNeqs.storeYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( !stream.write(renumberFlag) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    for ( auto &domain: domainList ) {
        domain->saveContext(stream, mode);
    }

    // store nMethod
    NumericalMethod *nmethod = this->giveNumericalMethod( this->giveMetaStep( giveCurrentStep()->giveMetaStepNumber() ) );
    if ( nmethod ) {
        nmethod->saveContext(stream, mode);
    }
}


void EngngModel :: restoreContext(DataStream &stream, ContextMode mode)
//
// this procedure is used mainly for two reasons:
//
// - to restore context of current model in order to be able to backtrace
//   computational history (useful in nonlinear analysis, when you want
//   explore another, previously detected load path
//
// - to restore context of current model in order to enable
//   possibility of post-processing.
//
// restoring context means: (if needed may be enhanced)
//
// - rest. EngngModel state varialbles as displacement, velocity, .. vectors
// - rest. Elements stress, strain and material history.
//
// This version loads only Element and Material properties.
//
// This function is inverse to the saveContext() member function
{
    contextIOResultType iores;

    // restore solution step
    int istep;
    if ( !stream.read(istep) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( !currentStep ) {
        currentStep = std::make_unique<TimeStep>(istep, this, 0, 0., 0., 0);
    }

    currentStep->restoreContext(stream);

    // this->updateAttributes (currentStep);

    int pmstep = currentStep->giveMetaStepNumber();
    if ( nMetaSteps ) {
        if ( !this->giveMetaStep(pmstep)->isStepValid(istep - 1) ) {
            pmstep--;
        }
    }

    previousStep = std::make_unique<TimeStep>(istep - 1, this, pmstep, currentStep->giveTargetTime ( ) - currentStep->giveTimeIncrement(),
                                currentStep->giveTimeIncrement(), currentStep->giveSolutionStateCounter() - 1);

    // restore numberOfEquations and domainNeqs array
    if ( !stream.read(numberOfEquations) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = domainNeqs.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // restore numberOfPrescribedEquations and domainNeqs array
    if ( !stream.read(numberOfPrescribedEquations) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = domainPrescribedNeqs.restoreYourself(stream) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // restore renumber flag
    if ( !stream.read(renumberFlag) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    for ( auto &domain: domainList ) {
        domain->restoreContext(stream, mode);
    }

    // restore nMethod
    NumericalMethod *nmethod = this->giveNumericalMethod( this->giveCurrentMetaStep() );
    if ( nmethod ) {
        nmethod->restoreContext(stream, mode);
    }

    this->updateDomainLinks();
    this->updateAttributes( this->giveCurrentMetaStep() );
    this->initStepIncrements();
}


MetaStep *
EngngModel :: giveCurrentMetaStep()
{
    return this->giveMetaStep( this->giveCurrentStep()->giveMetaStepNumber() );
}


std ::string
EngngModel :: giveContextFileName(int tStepNumber, int stepVersion) const
{
    std :: string fname = this->coreOutputFileName;
    char fext [ 100 ];
    sprintf(fext, ".%d.%d.osf", tStepNumber, stepVersion);
    return fname + fext;
}


std :: string
EngngModel :: giveDomainFileName(int domainNum, int domainSerNum) const
{
    std :: string fname = this->coreOutputFileName;
    char fext [ 100 ];
    sprintf(fext, ".domain.%d.%d.din", domainNum, domainSerNum);
    return fname + fext;
}

std :: string
EngngModel :: errorInfo(const char *func) const
{
    if ( this->isParallel() ) {
        return std::string(this->giveClassName()) + "::" + func + ", Rank: " + std::to_string(rank);
    } else {
        return std::string(this->giveClassName()) + "::" + func;
    }
}

Domain *
EngngModel :: giveDomain(int i)
{
    if ( ( i > 0 ) && ( i <= (int)this->domainList.size() ) ) {
        return this->domainList[i-1].get();
    } else {
        OOFEM_ERROR("Undefined domain");
    }

    return NULL;
}

void
EngngModel :: setDomain(int i, Domain *ptr, bool iDeallocateOld)
{
    if ( i < 1 || i > (int)this->domainList.size() ) {
        OOFEM_ERROR("Domain index %d out of range [1,%d]", i, (int)this->domainList.size());
    }
    if ( !iDeallocateOld ) {
        this->domainList[i-1].release();
    }
    this->domainList[i-1].reset(ptr);
}


ParallelContext *
EngngModel :: giveParallelContext(int i)
{
    if  ( i > (int)parallelContextList.size() ) {
        OOFEM_ERROR("context not initialized for this problem");
    }

    return &this->parallelContextList[i-1];
}

void
EngngModel :: initParallelContexts()
{
    parallelContextList.clear();
    for ( int i = 0; i < this->giveNumberOfDomains(); ++i ) {
        parallelContextList.emplace_back(this);
    }
}


MetaStep *
EngngModel :: giveMetaStep(int i)
{
    if ( ( i > 0 ) && ( i <= this->nMetaSteps ) ) {
        return &this->metaStepList[i-1];
    } else {
        OOFEM_ERROR("undefined metaStep (%d)", i);
    }

    return NULL;
}

void
EngngModel :: letOutputBaseFileNameBe(const std :: string &src)
{
    this->dataOutputFileName = src;

    if ( outputStream ) fclose(outputStream);

    if ( !suppressOutput ) {
        if ( ( outputStream = fopen(this->dataOutputFileName.c_str(), "w") ) == NULL ) {
            OOFEM_ERROR("Can't open output file %s", this->dataOutputFileName.c_str());
        }
    }
}

FILE *
EngngModel :: giveOutputStream()
// Returns an output stream on the data file of the receiver.
{
    if ( !outputStream ) {
        OOFEM_ERROR("No output stream opened!");
    }

    return outputStream;
}

double
EngngModel :: giveSolutionStepTime()
{
    return this->timer.getUtime(EngngModelTimer :: EMTT_SolutionStepTimer);
}

void
EngngModel :: giveAnalysisTime(int &rhrs, int &rmin, int &rsec, int &uhrs, int &umin, int &usec)
{
    double rtsec = this->timer.getWtime(EngngModelTimer :: EMTT_AnalysisTimer);
    double utsec = this->timer.getUtime(EngngModelTimer :: EMTT_AnalysisTimer);
    rsec = rmin = rhrs = 0;
    usec = umin = uhrs = 0;
    this->timer.convert2HMS(rhrs, rmin, rsec, rtsec);
    this->timer.convert2HMS(uhrs, umin, usec, utsec);
}

void
EngngModel :: terminateAnalysis()
{
    int rsec = 0, rmin = 0, rhrs = 0;
    int usec = 0, umin = 0, uhrs = 0;
    time_t endTime = time(NULL);
    this->timer.stopTimer(EngngModelTimer :: EMTT_AnalysisTimer);


    // compute real time consumed
    this->giveAnalysisTime(rhrs, rmin, rsec, uhrs, umin, usec);

    if(!suppressOutput) {
        FILE *out = this->giveOutputStream();
        fprintf(out, "\nFinishing analysis on: %s\n", ctime(& endTime) );
        fprintf(out, "Real time consumed: %03dh:%02dm:%02ds\n", rhrs, rmin, rsec);
        fprintf(out, "User time consumed: %03dh:%02dm:%02ds\n\n\n", uhrs, umin, usec);
    }

    OOFEM_LOG_FORCED("\n\nANALYSIS FINISHED\n\n\n");
    OOFEM_LOG_FORCED("Real time consumed: %03dh:%02dm:%02ds\n", rhrs, rmin, rsec);
    OOFEM_LOG_FORCED("User time consumed: %03dh:%02dm:%02ds\n", uhrs, umin, usec);
    exportModuleManager.terminate();
}

int
EngngModel :: checkProblemConsistency()
{
    int result = 1;

    result &= this->checkConsistency();

    for ( auto &domain: domainList ) {
        result &= domain->checkConsistency();
    }

#  ifdef VERBOSE
    if ( result ) {
        OOFEM_LOG_DEBUG("Consistency check:  OK\n");
    } else {
        VERBOSE_PRINTS("Consistency check", "failed")
        exit(1);
    }

#  endif

    return result;
}


void
EngngModel :: postInitialize()
{
    // set meta step bounds
    int istep = this->giveNumberOfFirstStep(true);
    for ( auto &metaStep: metaStepList ) {
        istep = metaStep.setStepBounds(istep);
    }

    for ( auto &domain: domainList ) {
        domain->postInitialize();
    }
}

void
EngngModel :: init()
{
    initModuleManager.doInit();
}


void
EngngModel :: initParallel()
{
 #ifdef __USE_MPI
    int len;
    MPI_Get_processor_name(processor_name, & len);
    this->comm = MPI_COMM_WORLD;
    MPI_Comm_rank(this->comm, & this->rank);
    MPI_Comm_size(this->comm, & numProcs);
 #else
    OOFEM_ERROR("Can't do it, only compiled for sequential runs");
 #endif
 #ifdef __VERBOSE_PARALLEL
    OOFEM_LOG_RELEVANT("[%d/%d] Running on %s\n", rank, numProcs, processor_name);
 #endif
}


#ifdef __OOFEG
void EngngModel :: drawYourself(oofegGraphicContext &gc)
{
    OGC_PlotModeType plotMode = gc.giveIntVarPlotMode();

    if ( ( plotMode == OGC_nodeAnnotation ) || ( plotMode == OGC_nodeGeometry ) || ( plotMode == OGC_essentialBC ) ||
        ( plotMode == OGC_naturalBC ) || ( plotMode == OGC_nodeScalarPlot ) || ( plotMode == OGC_nodeVectorPlot ) ) {
        this->drawNodes(gc);
    } else {
        this->drawElements(gc);
    }
}

void EngngModel :: drawElements(oofegGraphicContext &gc)
{
    Domain *d = this->giveDomain( gc.getActiveDomain() );
    TimeStep *tStep = this->giveCurrentStep();
    for ( auto &elem : d->giveElements() ) {
        elem->drawYourself(gc, tStep);
    }
}

void EngngModel :: drawNodes(oofegGraphicContext &gc)
{
    Domain *d = this->giveDomain( gc.getActiveDomain() );
    TimeStep *tStep = this->giveCurrentStep();
    for ( auto &dman : d->giveDofManagers() ) {
        dman->drawYourself(gc, tStep);
    }
}

#endif


void
EngngModel :: initializeCommMaps(bool forceInit)
{
#ifdef __PARALLEL_MODE
    // Set up communication patterns.
    communicator->setUpCommunicationMaps(this, true, forceInit);
    if ( nonlocalExt ) {
        nonlocCommunicator->setUpCommunicationMaps(this, true, forceInit);
    }
#else
    OOFEM_ERROR("Can't set up comm maps, parallel support not compiled");
#endif
}


int
EngngModel :: updateSharedDofManagers(FloatArray &answer, const UnknownNumberingScheme &s, int ExchangeTag)
{
    if ( isParallel() ) {
#ifdef __PARALLEL_MODE
        int result = 1;
 #ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT( "EngngModel :: updateSharedDofManagers", "Packing data", this->giveRank() );
 #endif

        ArrayWithNumbering tmp;
        tmp.array = & answer;
        tmp.numbering = & s;
        result &= communicator->packAllData(this, & tmp, & EngngModel :: packDofManagers);

 #ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT( "EngngModel :: updateSharedDofManagers", "Exchange started", this->giveRank() );
 #endif

        result &= communicator->initExchange(ExchangeTag);

 #ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT( "EngngModel :: updateSharedDofManagers", "Receiving and unpacking", this->giveRank() );
 #endif

        result &= communicator->unpackAllData(this, & tmp, & EngngModel :: unpackDofManagers);
        result &= communicator->finishExchange();
        return result;
#else
        OOFEM_ERROR("Support for parallel mode not compiled in.");
        return 0;
#endif
    } else {
        return 1;
    }

}


int
EngngModel :: exchangeRemoteElementData(int ExchangeTag)
{

    if ( isParallel() && nonlocalExt ) {
#ifdef __PARALLEL_MODE
        int result = 1;
 #ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT( "EngngModel :: exchangeRemoteElementData", "Packing remote element data", this->giveRank() );
 #endif

        result &= nonlocCommunicator->packAllData(this, & EngngModel :: packRemoteElementData);

 #ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT( "EngngModel :: exchangeRemoteElementData", "Remote element data exchange started", this->giveRank() );
 #endif

        result &= nonlocCommunicator->initExchange(ExchangeTag);

 #ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT( "EngngModel :: exchangeRemoteElementData", "Receiveng and Unpacking remote element data", this->giveRank() );
 #endif

        if ( !( result &= nonlocCommunicator->unpackAllData(this, & EngngModel :: unpackRemoteElementData) ) ) {
            OOFEM_ERROR("Receiveng and Unpacking remote element data");
        }

        result &= nonlocCommunicator->finishExchange();
        return result;
#else
        OOFEM_ERROR("Support for parallel mode not compiled in.");
        return 0;
#endif
    } else {
        return 1;
    }
}

#ifdef __PARALLEL_MODE
void
EngngModel :: balanceLoad(TimeStep *tStep)
{
    this->giveLoadBalancerMonitor();
    this->giveLoadBalancer();
    if ( !lb ) {
        OOFEM_WARNING("No load balancer found, skipping load balancing step");
        return;
    }

    //print statistics for current step
    lb->printStatistics();

    if ( tStep->isNotTheLastStep() ) {
        LoadBalancerMonitor :: LoadBalancerDecisionType _d = lbm->decide(tStep);
        if ( ( _d == LoadBalancerMonitor :: LBD_RECOVER ) ||
            ( ( tStep->isTheFirstStep() ) && force_load_rebalance_in_first_step ) ) {
            this->timer.startTimer(EngngModelTimer :: EMTT_LoadBalancingTimer);

            // determine nwe partitioning
            lb->calculateLoadTransfer();
            // pack e-model solution data into dof dictionaries
            this->packMigratingData(tStep);
            // migrate data
            lb->migrateLoad( this->giveDomain(1) );
            // renumber itself
            this->forceEquationNumbering();
            // re-init export modules
            this->giveExportModuleManager()->initialize();
 #ifdef __VERBOSE_PARALLEL
            // debug print
            EModelDefaultEquationNumbering dn;
            int nnodes = giveDomain(1)->giveNumberOfDofManagers();
            int myrank = this->giveRank();
            fprintf(stderr, "\n[%d] Nodal Table\n", myrank);
            for ( int i = 1; i <= nnodes; i++ ) {
                if ( giveDomain(1)->giveDofManager(i)->giveParallelMode() == DofManager_local ) {
                    fprintf( stderr, "[%d]: %5d[%d] local ", myrank, i, giveDomain(1)->giveDofManager(i)->giveGlobalNumber() );
                } else if ( giveDomain(1)->giveDofManager(i)->giveParallelMode() == DofManager_shared ) {
                    fprintf( stderr, "[%d]: %5d[%d] shared ", myrank, i, giveDomain(1)->giveDofManager(i)->giveGlobalNumber() );
                }

                for ( Dof *dof: *giveDomain(1)->giveDofManager(i) ) {
                    fprintf( stderr, "(%d)", dof->giveEquationNumber(dn) );
                }

                fprintf(stderr, "\n");
            }

 #endif
            // unpack (restore) e-model solution data from dof dictionaries
            this->unpackMigratingData(tStep);

            this->timer.stopTimer(EngngModelTimer :: EMTT_LoadBalancingTimer);
            double _steptime = this->timer.getUtime(EngngModelTimer :: EMTT_LoadBalancingTimer);
            OOFEM_LOG_INFO("[%d] EngngModel info: user time consumed by load rebalancing %.2fs\n",
                           this->giveRank(), _steptime);
        }
    }
}


int
EngngModel :: packRemoteElementData(ProcessCommunicator &processComm)
{
    int result = 1;
    const IntArray &toSendMap = processComm.giveToSendMap();
    CommunicationBuffer &send_buff = processComm.giveProcessCommunicatorBuff()->giveSendBuff();
    Domain *domain = this->giveDomain(1);


    for ( int ielem : toSendMap ) {
        result &= domain->giveElement( ielem )->packUnknowns( send_buff, this->giveCurrentStep() );
    }

    return result;
}


int
EngngModel :: unpackRemoteElementData(ProcessCommunicator &processComm)
{
    int result = 1;
    const IntArray &toRecvMap = processComm.giveToRecvMap();
    CommunicationBuffer &recv_buff = processComm.giveProcessCommunicatorBuff()->giveRecvBuff();
    Domain *domain = this->giveDomain(1);


    for ( int ielem : toRecvMap ) {
        Element *element = domain->giveElement( ielem );
        if ( element->giveParallelMode() == Element_remote ) {
            result &= element->unpackAndUpdateUnknowns( recv_buff, this->giveCurrentStep() );
        } else {
            OOFEM_ERROR("element is not remote");
        }
    }

    return result;
}


int
EngngModel :: packDofManagers(ArrayWithNumbering *srcData, ProcessCommunicator &processComm)
{
    FloatArray *src = srcData->array;
    const UnknownNumberingScheme &s = * srcData->numbering;
    int result = 1;
    Domain *domain = this->giveDomain(1);
    const IntArray &toSendMap = processComm.giveToSendMap();
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();

    ///@todo Shouldn't hardcode domain number 1
    ///@todo Must fix: Internal dofmanagers in xfem and bc
    for ( int inode : toSendMap ) {
        DofManager *dman = domain->giveDofManager( inode );
        for ( auto &jdof: *dman ) {
            if ( jdof->isPrimaryDof() ) {
                int eqNum = jdof->giveEquationNumber(s);
                if ( eqNum ) {
                    result &= pcbuff->write( src->at(eqNum) );
                }
            }
        }
    }

    return result;
}


int
EngngModel :: unpackDofManagers(ArrayWithNumbering *destData, ProcessCommunicator &processComm)
{
    FloatArray *dest = destData->array;
    const UnknownNumberingScheme &s = * destData->numbering;
    int result = 1;
    Domain *domain = this->giveDomain(1);
    const IntArray &toRecvMap = processComm.giveToRecvMap();
    ProcessCommunicatorBuff *pcbuff = processComm.giveProcessCommunicatorBuff();
    double value;

    ///@todo Shouldn't hardcode domain number 1
    ///@todo Must fix: Internal dofmanagers in bc
    for ( int inode : toRecvMap ) {
        DofManager *dman = domain->giveDofManager( inode );
        dofManagerParallelMode dofmanmode = dman->giveParallelMode();
        for ( auto &jdof: *dman ) {
            int eqNum = jdof->giveEquationNumber(s);
            if ( jdof->isPrimaryDof() && eqNum ) {
                result &= pcbuff->read(value);
                if ( dofmanmode == DofManager_shared ) {
                    dest->at(eqNum) += value;
                } else if ( dofmanmode == DofManager_remote ) {
                    dest->at(eqNum)  = value;
                } else {
                    OOFEM_ERROR("unknown dof namager parallel mode");
                }
            }
        }
    }

    return result;
}

#endif
} // end namespace oofem
