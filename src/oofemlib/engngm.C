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

// Milan ?????????????????
//#include "gpinitmodule.h"
// Milan ?????????????????

#include "nummet.h"
#include "sparsemtrx.h"
#include "engngm.h"
#include "timestep.h"
#include "metastep.h"
#include "element.h"
#include "oofemdef.h"
#include "clock.h"
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
#include "xfemmanager.h"
#include "outputmanager.h"
#include "exportmodulemanager.h"
#include "initmodulemanager.h"
#include "usrdefsub.h"

#ifdef __CEMHYD_MODULE
 #include "cemhydmat.h"
#endif

#ifndef __MAKEDEPEND
 #include <cstdio>
 #include <cstdarg>
// include unistd.h; needed for access
 #ifdef HAVE_UNISTD_H
  #include <unistd.h>
 #elif _MSC_VER
  #include <io.h>
 #endif
#endif

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

#ifdef TIME_REPORT
 #ifndef __MAKEDEPEND
  #include <time.h>
 #endif
#endif

namespace oofem {
EngngModel :: EngngModel(int i, EngngModel *_master) : domainNeqs(), domainPrescribedNeqs()
// constructor
{

    number = i;
    currentStep = NULL;
    previousStep = NULL;
    stepWhenIcApply = NULL;
    defaultErrEstimator = NULL;
    numberOfSteps = 0;
    numberOfEquations = 0;
    numberOfPrescribedEquations = 0;
    renumberFlag = false;
    equationNumberingCompleted = 0;
    ndomains = 0;
    nMetaSteps = 0;
    nxfemman = 0;
    nonLinFormulation = UNKNOWN;

    outputStream          = NULL;

    domainList            = new AList< Domain >(0);
    metaStepList          = new AList< MetaStep >(0);
    xfemManagerList       = new AList< XfemManager >(0);

    contextOutputMode     = COM_NoContext;
    contextOutputStep     = 0;
    pMode                 = _processor;  // for giveContextFile()
    pScale                = macroScale;

    exportModuleManager   = new ExportModuleManager(this);
    initModuleManager     = new InitModuleManager(this);
    master                = _master; // master mode by default
    // create context if in master mode; otherwise request context from master
    if ( master ) {
        context = master->giveContext();
    } else {
        context = new EngngModelContext();
    }

    parallelFlag = 0;
#ifdef __PARALLEL_MODE
    loadBalancingFlag = false;
    force_load_rebalance_in_first_step = false;
    lb = NULL;
    lbm = NULL;
#endif

#ifdef __PETSC_MODULE
    petscContextList = new AList< PetscContext >(0);
#endif
}

EngngModel :: EngngModel(int i, char *s, EngngModel *_master) : domainNeqs(), domainPrescribedNeqs()
// constructor
{
    number = i;
    currentStep = NULL;
    previousStep = NULL;
    stepWhenIcApply = NULL;
    defaultErrEstimator = NULL;
    numberOfSteps = 0;
    numberOfEquations = 0;
    numberOfPrescribedEquations = 0;
    renumberFlag = false;
    equationNumberingCompleted = 0;
    ndomains = 0;
    nMetaSteps = 0;
    nxfemman = 0;

    outputStream          = NULL;

    domainList            = new AList< Domain >(0);
    metaStepList          = new AList< MetaStep >(0);
    xfemManagerList       = new AList< XfemManager >(0);

    exportModuleManager   = new ExportModuleManager(this);
    initModuleManager     = new InitModuleManager(this);
    master                = _master; // master mode by default
    // create context if in master mode; otherwise request context from master
    if ( master ) {
        context = master->giveContext();
    } else {
        context = new EngngModelContext();
    }

    parallelFlag = 0;
#ifdef __PARALLEL_MODE
    loadBalancingFlag = false;
    force_load_rebalance_in_first_step = false;
    lb = NULL;
    lbm = NULL;
#endif

#ifdef __PETSC_MODULE
    petscContextList = new AList< PetscContext >(0);
#endif
}


EngngModel ::  ~EngngModel()
// destructor
{
    if (previousStep == currentStep) {
        if(previousStep != NULL) {
            delete this->currentStep;
        }
    } else {
        if(currentStep != NULL) {
            delete currentStep;
        }
        if(previousStep != NULL) {
            delete previousStep;
        }
    }

    if(stepWhenIcApply != NULL){
      delete stepWhenIcApply;
    }

    delete domainList;
    delete metaStepList;
    delete xfemManagerList;

#ifdef __PETSC_MODULE
    delete petscContextList;
#endif

    if ( exportModuleManager ) {
        delete exportModuleManager;
    }

    if ( initModuleManager ) {
        delete initModuleManager;
    }

    // master deletes the context
    if ( master == NULL ) {
        delete context;
    }

    //fclose (inputStream) ;
    if ( outputStream ) {
        fclose(outputStream);
    }

    if ( defaultErrEstimator ) {
        delete defaultErrEstimator;
    }

#ifdef __PARALLEL_MODE
    if ( loadBalancingFlag ) {
        if ( lb ) {
            delete lb;
        }

        if ( lbm ) {
            delete lbm;
        }
    }
#endif
}


void EngngModel :: setParallelMode(bool parallelFlag)
{
    this->parallelFlag = parallelFlag;
if ( this->parallelFlag ) {
#ifndef __PARALLEL_MODE
        OOFEM_ERROR("EngngModel :: setParallelMode - Can't do it, only compiled for sequential runs");
#else
        initParallel();
#endif
    }
}


int EngngModel :: instanciateYourself(DataReader *dr, InputRecord *ir, const char *dataOutputFileName, const char *desc)
// simple input - only number of steps variable is read
{
    Domain *domain;
    int i;
    bool inputReaderFinish = true;

    this->coreOutputFileName = std::string(dataOutputFileName);
    this->dataOutputFileName = std::string(dataOutputFileName);

    if ( this->giveProblemMode() == _postProcessor ) {
        // modify output file name to prevent output to be lost
        this->dataOutputFileName.append(".oofeg");
    }

    if ( ( outputStream = fopen(this->dataOutputFileName.c_str(), "w") ) == NULL ) {
        _error2("instanciateYourself: Can't open output file %s", this->dataOutputFileName.c_str());
    }

    fprintf(outputStream, "%s", PRG_HEADER);
    this->startTime = getTime();
    //this->startClock= this-> getClock();
    fprintf( outputStream, "\nStarting analysis on: %s\n", ctime(& this->startTime) );

    fprintf(outputStream, "%s\n", desc);

#  ifdef VERBOSE
    OOFEM_LOG_DEBUG("Reading all data from input file \n");
#  endif
#ifdef __PARALLEL_MODE
    if ( this->isParallel() )
        fprintf(outputStream, "Problem rank is %d/%d on %s\n\n", this->rank, this->numProcs, this->processor_name);
#endif

    // create domains
    domainNeqs.resize(this->ndomains);
    domainPrescribedNeqs.resize(this->ndomains);
    domainList->growTo(ndomains);
    for ( i = 0; i < this->ndomains; i++ ) {
        domain =  new Domain(i + 1, 0, this);
        domainList->put(i + 1, domain);
    }

#ifdef __PETSC_MODULE
    this->initPetscContexts();
#endif

    // instanciate receiver
    this->initializeFrom(ir);
    exportModuleManager->initializeFrom(ir);
    initModuleManager->initializeFrom(ir);

    if ( this->nMetaSteps == 0 ) {
        inputReaderFinish = false;
        this->instanciateDefaultMetaStep(ir);
    } else {
        this->instanciateMetaSteps(dr);
    }

    // instanciate initialization module manager
    initModuleManager->instanciateYourself(dr, ir);
    // instanciate export module manager
    exportModuleManager->instanciateYourself(dr, ir);
    this->instanciateDomains(dr);

    exportModuleManager->initialize();

    // Milan ??????????????????
    //GPImportModule* gim = new GPImportModule(this);
    //gim -> getInput();
    // Milan ??????????????????

    // instantiate xfemmanager
    XfemManager *xm;
    xfemManagerList->growTo(nxfemman);
    for ( i = 0; i < this->nxfemman; i++ ) {
        xm =  new XfemManager(this, i + 1);
        ir = dr->giveInputRecord(DataReader :: IR_xfemManRec, 1);
        // XfemManager has to be put into xfemManagerList before xm->initializeFrom, otherwise Enrichmentitem cannot access XfemManager
        // or we have to make a reference from EnrichmentItem also
        xfemManagerList->put(i + 1, xm);
        xm->initializeFrom(ir);
        xm->instanciateYourself(dr);
        int last = xm->computeFictPosition();
        this->setNumberOfEquations(1, last);
        xm->updateIntegrationRule();
    }

    // check emodel input record if no default metastep, since all has been read
    if ( inputReaderFinish ) {
        ir->finish();
    }

    return 1;
}

IRResultType
EngngModel :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, numberOfSteps, IFT_EngngModel_nsteps, "nsteps"); // Macro
    if ( numberOfSteps <= 0 ) {
        _error("instanciateFrom: nsteps not specified, bad format");
    }

    contextOutputStep =  0;
    IR_GIVE_OPTIONAL_FIELD(ir, contextOutputStep, IFT_EngngModel_contextoutputstep, "contextoutputstep"); // Macro
    if ( contextOutputStep ) {
        this->setUDContextOutputMode(contextOutputStep);
    }

    renumberFlag = false;
    IR_GIVE_OPTIONAL_FIELD(ir, renumberFlag, IFT_EngngModel_renumberFlag, "renumber");                // Macro
    profileOpt = false;
    IR_GIVE_OPTIONAL_FIELD(ir, profileOpt, IFT_EngngModel_profileOpt, "profileopt");                // Macro
    nMetaSteps   = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, nMetaSteps, IFT_EngngModel_nmsteps, "nmsteps");                // Macro
    nxfemman   = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, nxfemman, IFT_EngngModel_nxfemman, "nxfemman");          // Macro
    int _val = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, IFT_EngngModel_nonLinFormulation, "nonlinform");
    nonLinFormulation = ( fMode ) _val;

    int eeTypeId = -1;
    IR_GIVE_OPTIONAL_FIELD(ir, eeTypeId, IFT_EngngModel_eetype, "eetype");
    if (eeTypeId >= 0) {
        this->defaultErrEstimator = CreateUsrDefErrorEstimator( ( ErrorEstimatorType ) eeTypeId, 1, this->giveDomain(1) );
        this->defaultErrEstimator->initializeFrom(ir);
    }

#ifdef __PARALLEL_MODE
    IR_GIVE_OPTIONAL_FIELD(ir, parallelFlag, IFT_EngngModel_parallelflag, "parallelflag"); // Macro
    // fprintf (stderr, "Parallel mode is %d\n", parallelFlag);

    /* Load balancing support */
    _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, IFT_NonLinearStatic_loadBalancingFlag, "lbflag"); // Macro
    loadBalancingFlag = _val;

    _val = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _val, IFT_NonLinearStatic_forceloadBalancingFlag, "forcelb1"); // Macro
    force_load_rebalance_in_first_step = _val;

#endif
    return IRRT_OK;
}

int
EngngModel :: instanciateDomains(DataReader *dr)
{
    int i, result = 1;
    // read problem domains
    for ( i = 0; i < this->ndomains; i++ ) {
        result &= domainList->at(i + 1)->instanciateYourself(dr);
    }

    return result;
}


int
EngngModel :: instanciateMetaSteps(DataReader *dr)
{
    int i, result = 1;
    MetaStep *mstep;

    // create meta steps
    metaStepList->growTo(nMetaSteps);
    for ( i = 0; i < this->nMetaSteps; i++ ) {
        mstep =  new MetaStep(i + 1, this);
        metaStepList->put(i + 1, mstep);
    }

    // read problem domains
    for ( i = 0; i < this->nMetaSteps; i++ ) {
        InputRecord *ir = dr->giveInputRecord(DataReader :: IR_mstepRec, i + 1);
        result &= metaStepList->at(i + 1)->initializeFrom(ir);
    }

    // set meta step bounds
    int istep = this->giveNumberOfFirstStep();
    for ( i = 0; i < this->nMetaSteps; i++ ) {
        istep = metaStepList->at(i + 1)->setStepBounds(istep);
    }

    this->numberOfSteps = istep - 1;
    OOFEM_LOG_RELEVANT("Total number of solution steps     %d\n", numberOfSteps);
    return result;
}


int
EngngModel :: instanciateDefaultMetaStep(InputRecord *ir)
{
    MetaStep *mstep;

    if ( numberOfSteps == 0 ) {
        _error("instanciateDefaultMetaStep: nsteps cannot be zero");
    }

    // create default meta steps
    this->nMetaSteps = 1;
    metaStepList->growTo(nMetaSteps);
    mstep =  new MetaStep(1, this, numberOfSteps, *ir);
    metaStepList->put(1, mstep);

    // set meta step bounds
    int istep = this->giveNumberOfFirstStep() - 1;
    metaStepList->at(1)->setStepBounds(istep + 1);

    OOFEM_LOG_RELEVANT("Total number of solution steps     %d\n",  numberOfSteps);
    return 1;
}


int
EngngModel :: giveNumberOfEquations(EquationID)
{
    //
    // returns number of equations of current problem
    // this method is implemented here, because some method may add some
    // conditions in to system and this may results into increased number of
    // equations.
    //
    if ( equationNumberingCompleted ) {
        return numberOfEquations;
    }

    return this->forceEquationNumbering();
}

int
EngngModel :: giveNumberOfPrescribedEquations(EquationID)
{
    //
    // returns number of equations of current problem
    // this method is implemented here, because some method may add some
    // conditions in to system and this may results into increased number of
    // equations.
    //
    if ( !equationNumberingCompleted ) {
        this->forceEquationNumbering();
    }

    return numberOfPrescribedEquations;
}

int
EngngModel :: giveNumberOfDomainEquations(int id, EquationID)
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

    return domainNeqs.at(id);
}


int
EngngModel :: giveNumberOfPrescribedDomainEquations(int id, EquationID)
{
    //
    // returns number of equations of current problem
    // this method is implemented here, because some method may add some
    // conditions in to system and this may results into increased number of
    // equations.
    //
    if ( equationNumberingCompleted ) {
        return domainPrescribedNeqs.at(id);
    }

    this->forceEquationNumbering();
    return domainPrescribedNeqs.at(id);
}

int
EngngModel :: forceEquationNumbering(int id)
{
    // forces equation renumbering for current time step
    // intended mainly for problems with changes of static system
    // during solution
    // OUTPUT:
    // sets this->numberOfEquations and this->numberOfPrescribedEquations and returns this value

    int i, k, nnodes, nelem, nbc;
    Element *elem;
    Domain *domain = this->giveDomain(id);
    TimeStep *currStep = this->giveCurrentStep();
    IntArray loc;

    this->domainNeqs.at(id) = 0;
    this->domainPrescribedNeqs.at(id) = 0;

    nnodes = domain->giveNumberOfDofManagers();
    nelem  = domain->giveNumberOfElements();
    nbc    = domain->giveNumberOfBoundaryConditions();

    if ( !this->profileOpt ) {
        for ( i = 1; i <= nnodes; i++ ) {
            domain->giveDofManager(i)->askNewEquationNumbers(currStep);
        }
        for ( i = 1; i <= nelem; ++i ) {
            elem = domain->giveElement(i);
            nnodes = elem->giveNumberOfInternalDofManagers();
            for ( k = 1; k <= nnodes; k++ ) {
                elem->giveInternalDofManager(k)->askNewEquationNumbers(currStep);
            }
        }
        // For special boundary conditions;
        for ( i = 1; i <= nbc; ++i ) {
            GeneralBoundaryCondition *bc = domain->giveBc(i);
            nnodes = bc->giveNumberOfInternalDofManagers();
            for ( k = 1; k <= nnodes; k++ ) {
                bc->giveInternalDofManager(k)->askNewEquationNumbers(currStep);
            }
        }
    } else {
        // invoke profile reduction
        int initialProfile, optimalProfile;
        //clock_t time_0 = this->getClock(), time_1;
        oofem_timeval tstart;
        getUtime(tstart);
        OOFEM_LOG_INFO("\nRenumbering DOFs with Sloan's algorithm...\n");

        SloanGraph graph(domain);
        graph.initialize();
        graph.tryParameters(0, 0);
        initialProfile = graph.giveOptimalProfileSize();
        graph.tryParameters(2, 1);
        graph.tryParameters(1, 0);
        graph.tryParameters(5, 1);
        graph.tryParameters(10, 1);
        optimalProfile = graph.giveOptimalProfileSize();

        oofem_timeval ut;
        getRelativeUtime(ut, tstart);

        OOFEM_LOG_DEBUG("done in %.2fs\n", ( double ) ( ut.tv_sec + ut.tv_usec / ( double ) OOFEM_USEC_LIM ) );
        OOFEM_LOG_DEBUG("Nominal profile %d (old) %d (new)\n", initialProfile, optimalProfile);

        //FILE* renTableFile = fopen ("rentab.dat","w");
        //graph.writeOptimalRenumberingTable (renTableFile);
        graph.askNewOptimalNumbering(currStep);
    }

    // invalidate element local copies of location arrays
    nelem = domain->giveNumberOfElements();
    for ( i = 1; i <= nelem; i++ ) {
        domain->giveElement(i)->invalidateLocationArray();
    }

    return domainNeqs.at(id);
}


int
EngngModel :: forceEquationNumbering()
{
    int i;
    // set numberOfEquations counter to zero
    this->numberOfEquations = 0;
    this->numberOfPrescribedEquations = 0;

    OOFEM_LOG_DEBUG("Renumbering dofs in all domains\n");
    for ( i = 1; i <= this->ndomains; i++ ) {
        domainNeqs.at(i) = 0;
        this->numberOfEquations += this->forceEquationNumbering(i);
    }

    equationNumberingCompleted = 1;

    for ( i = 1; i <= this->ndomains; i++ ) {
        //this->numberOfPrescribedEquations+=giveNumberOfPrescribedDomainEquations(i);
        this->numberOfPrescribedEquations += domainPrescribedNeqs.at(i);
    }

#ifdef __PETSC_MODULE
    for ( i = 1; i <= petscContextList->giveSize(); i++ ) {
        this->petscContextList->at(i)->init(i);
    }

#endif


    return this->numberOfEquations;
}


void
EngngModel :: solveYourself()
{
    int imstep, jstep, nTimeSteps;
    int smstep = 1, sjstep = 1;
    MetaStep *activeMStep;
    FILE *out = this->giveOutputStream();
    //#ifdef TIME_REPORT
    this->timer.startTimer(EngngModelTimer :: EMTT_AnalysisTimer);
    //#endif

    if ( this->currentStep ) {
        smstep = this->currentStep->giveMetaStepNumber();
        sjstep = this->giveMetaStep(smstep)->giveStepRelativeNumber( this->currentStep->giveNumber() ) + 1;
    }

    for ( imstep = smstep; imstep <= nMetaSteps; imstep++, sjstep = 1 ) { //loop over meta steps
        activeMStep = this->giveMetaStep(imstep);
        // update state according to new meta step
        this->initMetaStepAttributes( activeMStep );

        nTimeSteps = activeMStep->giveNumberOfSteps();
        for ( jstep = sjstep; jstep <= nTimeSteps; jstep++ ) { //loop over time steps
            //#ifdef TIME_REPORT
            this->timer.startTimer(EngngModelTimer :: EMTT_SolutionStepTimer);
            this->timer.initTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
            //#endif
            this->giveNextStep();

            // renumber equations if necessary. Ensure to call forceEquationNumbering() for staggered problems
            if ( this->requiresEquationRenumbering( this->giveCurrentStep() ) || this->giveClassID()== classType(StaggeredProblemClass) ) {
                this->forceEquationNumbering();
            }

            this->solveYourselfAt( this->giveCurrentStep() );
            this->updateYourself( this->giveCurrentStep() );
            this->terminate( this->giveCurrentStep() );

            //#ifdef TIME_REPORT
            this->timer.stopTimer(EngngModelTimer :: EMTT_SolutionStepTimer);
            //this->timer.stopTimer(EngngModelTimer::EMTT_NetComputationalStepTimer);
            double _steptime = this->timer.getUtime(EngngModelTimer :: EMTT_SolutionStepTimer);
            OOFEM_LOG_INFO("EngngModel info: user time consumed by solution step %d: %.2fs\n",
                           this->giveCurrentStep()->giveNumber(), _steptime);

            fprintf(out, "\nUser time consumed by solution step %d: %.3f [s]\n\n",
                    this->giveCurrentStep()->giveNumber(), _steptime);
            //#endif

#ifdef __PARALLEL_MODE
            if ( loadBalancingFlag ) {
                this->balanceLoad( this->giveCurrentStep() );
            }

#endif
        }
    }
}

void
EngngModel :: initMetaStepAttributes(MetaStep *mStep)
{
    // update attributes
    this->updateAttributes(mStep); // virtual function
    // finish data acquiring
    mStep->giveAttributesRecord()->finish();
}

void
EngngModel :: updateAttributes(MetaStep *mStep)
{
    InputRecord *ir = mStep->giveAttributesRecord();

    if ( this->giveNumericalMethod(mStep) ) {
        this->giveNumericalMethod(mStep)->initializeFrom(ir);
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
EngngModel :: updateYourself(TimeStep *stepN)
{
    int idomain, ndomains = this->giveNumberOfDomains();
    int j, nnodes;
    Domain *domain;

    for ( idomain = 1; idomain <= ndomains; idomain++ ) {
        domain = this->giveDomain(idomain);

#  ifdef VERBOSE
        VERBOSE_PRINT0( "Updating domain ", domain->giveNumber() )
#  endif

        nnodes = domain->giveNumberOfDofManagers();
        for ( j = 1; j <= nnodes; j++ ) {
            domain->giveDofManager(j)->updateYourself(stepN);
            domain->giveNode(j)->updateYourself(stepN);
            //domain->giveDofManager(j)->printOutputAt(File, stepN);
        }

#  ifdef VERBOSE
        VERBOSE_PRINT0("Updated nodes & sides ", nnodes)
#  endif


        Element *elem;

        int nelem = domain->giveNumberOfElements();
        for ( j = 1; j <= nelem; j++ ) {
            elem = domain->giveElement(j);
#ifdef __PARALLEL_MODE
            // skip remote elements (these are used as mirrors of remote elements on other domains
            // when nonlocal constitutive models are used. They introduction is necessary to
            // allow local averaging on domains without fine grain communication between domains).
            if ( elem->giveParallelMode() == Element_remote ) {
                continue;
            }

#endif
            elem->updateYourself(stepN);
            //elem -> printOutputAt(File, stepN) ;

#ifdef __CEMHYD_MODULE
            //store temperature and associated volume on each GP before performing averaging
            if ( elem->giveMaterial()->giveClassID() == CemhydMatClass ) {
                TransportElement *element = ( TransportElement * ) elem;
                CemhydMat *cem = ( CemhydMat * ) element->giveMaterial();
                cem->clearWeightTemperatureProductVolume(element);
                cem->storeWeightTemperatureProductVolume(element, stepN);
            }
#endif
        }

#  ifdef VERBOSE
        VERBOSE_PRINT0("Updated Elements ", nelem)
#  endif

#ifdef __CEMHYD_MODULE
        //perform averaging on each material instance
        for ( j = 1; j <= domain->giveNumberOfMaterialModels(); j++ ) {
            if ( domain->giveMaterial(j)->giveClassID() == CemhydMatClass ) {
                CemhydMat *cem = ( CemhydMat * ) domain->giveMaterial(j);
                cem->averageTemperature();
            }
        }
#endif

#  ifdef VERBOSE
        VERBOSE_PRINT0("Updated Materials ", nelem)
#  endif
    }

    // if there is an error estimator, it should be updated so that values can be exported.
    if ( this->defaultErrEstimator ) {
        this->defaultErrEstimator->estimateError( equilibratedEM, stepN );
    }
}

void
EngngModel :: terminate(TimeStep *stepN)
{
    this->doStepOutput(stepN);
    this->saveStepContext(stepN);
}


void
EngngModel :: doStepOutput(TimeStep *stepN)
{
    FILE *File = this->giveOutputStream();

    // print output
    this->printOutputAt(File, stepN);
    // export using export manager
    exportModuleManager->doOutput(stepN);
}

void
EngngModel :: saveStepContext(TimeStep *stepN)
{
    // save context if required
    // default - save only if ALWAYS is set ( see cltypes.h )

    if ( ( this->giveContextOutputMode() == COM_Always ) ||
        ( this->giveContextOutputMode() == COM_Required ) ) {
        this->saveContext(NULL, CM_State);
    } else if ( this->giveContextOutputMode() == COM_UserDefined ) {
        if ( stepN->giveNumber() % this->giveContextOutputStep() == 0 ) {
            this->saveContext(NULL, CM_State);
        }
    }
}


void
EngngModel :: printOutputAt(FILE *File, TimeStep *stepN)
{
    //FILE* File = this -> giveDomain() -> giveOutputStream() ;
    int domCount = 0, idomain;
    Domain *domain;

    // fprintf (File,"\nOutput for time step number %d \n\n",stepN->giveNumber());
    for ( idomain = 1; idomain <= this->ndomains; idomain++ ) {
        domain = this->giveDomain(idomain);
        domCount += domain->giveOutputManager()->testTimeStepOutput(stepN);
    }

    if ( domCount == 0 ) {
        return;              // do not print even Solution step header
    }

    fprintf(File, "\n==============================================================");
    fprintf( File, "\nOutput for time % .8e ", stepN->giveTargetTime() * this->giveVariableScale(VST_Time) );
    fprintf(File, "\n==============================================================\n");
    for ( idomain = 1; idomain <= this->ndomains; idomain++ ) {
        domain = this->giveDomain(idomain);
        fprintf( File, "Output for domain %3d\n", domain->giveNumber() );

        domain->giveOutputManager()->doDofManOutput(File, stepN);
        domain->giveOutputManager()->doElementOutput(File, stepN);
    }
}

void EngngModel :: printYourself()
{
    printf( "\nEngineeringModel: instance %s\n", this->giveClassName() );
    printf( "number of steps: %d\n", this->giveNumberOfSteps() );
    printf("number of eq's : %d\n", numberOfEquations);
}

void EngngModel :: assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
                            CharType type, const UnknownNumberingScheme &s, Domain *domain)
//
// assembles matrix answer by calling
// element(i) -> giveCharacteristicMatrix ( type, tStep );
// for each element in domain
// and assembling every contribution to answer
//
//
{
    int ielem;
    IntArray loc;
    FloatMatrix mat, R;
    Element *element;

    if ( answer == NULL ) {
        _error("assemble: NULL pointer encountered.");
    }

    this->timer.resumeTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
    int nelem = domain->giveNumberOfElements();
    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        element = domain->giveElement(ielem);
#ifdef __PARALLEL_MODE
        // skip remote elements (these are used as mirrors of remote elements on other domains
        // when nonlocal constitutive models are used. They introduction is necessary to
        // allow local averaging on domains without fine grain communication between domains).
        if ( element->giveParallelMode() == Element_remote ) {
            continue;
        }

#endif
        this->giveElementCharacteristicMatrix(mat, ielem, type, tStep, domain);

        if ( mat.isNotEmpty() ) {
            element->giveLocationArray(loc, eid, s);
            if (element->giveRotationMatrix(R, eid)) {
                mat.rotatedWith(R);
            }
            if ( answer->assemble(loc, mat) == 0 ) {
                _error("assemble: sparse matrix assemble error");
            }
        }
    }

    this->timer.pauseTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);

    answer->assembleBegin();
    answer->assembleEnd();
}

void EngngModel :: assemble(SparseMtrx *answer, TimeStep *tStep, EquationID r_id, EquationID c_id,
                            CharType type, const UnknownNumberingScheme &ns,
                            Domain *domain)
//
// assembles matrix answer by  calling
// element(i) -> giveCharacteristicMatrix ( type, tStep );
// for each element in domain
// and assembling every contribution to answer
//
//
{
    int ielem;
    IntArray r_loc, c_loc;
    FloatMatrix mat, Rr, Rc;
    Element *element;

    if ( answer == NULL ) {
        _error("assemble: NULL pointer encountered.");
    }

    this->timer.resumeTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
    int nelem = domain->giveNumberOfElements();
    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        element = domain->giveElement(ielem);
#ifdef __PARALLEL_MODE
        // skip remote elements (these are used as mirrors of remote elements on other domains
        // when nonlocal constitutive models are used. They introduction is necessary to
        // allow local averaging on domains without fine grain communication between domains).
        if ( element->giveParallelMode() == Element_remote ) {
            continue;
        }

#endif
        this->giveElementCharacteristicMatrix(mat, ielem, type, tStep, domain);

        if ( mat.isNotEmpty() ) {
            element->giveLocationArray(r_loc, r_id, ns);
            element->giveLocationArray(c_loc, c_id, ns);
            // Rotate it
            if (element->giveRotationMatrix(Rr, r_id)) {
                FloatMatrix tmpMat;
                tmpMat.beTProductOf(Rr, mat); ///@todo Check transpose here
                mat = tmpMat;
            }
            if (element->giveRotationMatrix(Rc, c_id)) {
                FloatMatrix tmpMat;
                tmpMat.beProductOf(mat, Rc); ///@todo Check transpose here
                mat = tmpMat;
            }
            if ( answer->assemble(r_loc, c_loc, mat) == 0 ) {
                _error("assemble: sparse matrix assemble error");
            }
        }
    }

    this->timer.pauseTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);

    answer->assembleBegin();
    answer->assembleEnd();
}

void EngngModel :: assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
                            CharType type, const UnknownNumberingScheme &rs, const UnknownNumberingScheme &cs,
                            Domain *domain)
// Same as assemble, but with different numbering for rows and columns
{
    int ielem;
    IntArray r_loc, c_loc;
    FloatMatrix mat, R;
    Element *element;
    if ( answer == NULL ) OOFEM_ERROR("EngngModel :: assemble: NULL pointer encountered.");

    this->timer.resumeTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
    int nelem = domain->giveNumberOfElements();
    for ( ielem = 1; ielem <= nelem; ielem++ ) {
        element = domain->giveElement(ielem);
#ifdef __PARALLEL_MODE
        if ( element->giveParallelMode() == Element_remote ) continue;
#endif
        this->giveElementCharacteristicMatrix(mat, ielem, type, tStep, domain);
        if ( mat.isNotEmpty() ) {
            element->giveLocationArray(r_loc, eid, rs);
            element->giveLocationArray(c_loc, eid, cs);
            // Rotate it
            if (element->giveRotationMatrix(R, eid)) {
                mat.rotatedWith(R);
            }
            if ( answer->assemble(r_loc, c_loc, mat) == 0 ) {
                OOFEM_ERROR("EngngModel :: assemble: sparse matrix assemble error");
            }
        }
    }

    this->timer.pauseTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);

    answer->assembleBegin();
    answer->assembleEnd();
}


double EngngModel :: assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                  CharType type, ValueModeType mode,
                                  const UnknownNumberingScheme &s, Domain *domain)
{
    double localNorm = 0.0;
    localNorm += this->assembleVectorFromDofManagers(answer, tStep, eid, type, mode, s, domain);
    localNorm += this->assembleVectorFromElements(answer, tStep, eid, type, mode, s, domain);
    localNorm += this->assembleVectorFromActiveBC(answer, tStep, eid, type, mode, s, domain);
#ifdef __PARALLEL_MODE
 #ifdef __PETSC_MODULE
    return this->givePetscContext(domain->giveNumber(), eid)->accumulate(localNorm);
 #else
    if (this->isParallel()) {
        double globalNorm;
        MPI_Allreduce(&localNorm, &globalNorm, 1, MPI_DOUBLE, MPI_SUM, comm);
        return globalNorm;
    }
 #endif
#else
    return localNorm;
#endif
}


double EngngModel :: assembleVectorFromDofManagers(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                                 CharType type, ValueModeType mode,
                                                 const UnknownNumberingScheme &s, Domain *domain)
//
// assembles matrix answer by  calling
// node(i) -> computeLoadVectorAt (tStep);
// for each element in domain
// and assembling every contribution to answer
//
//
{
    if ( type != ExternalForcesVector ) { // Dof managers can only have external loads.
        return 0.0;
    }
    IntArray loc;
    FloatArray charVec;
    FloatMatrix R;
    IntArray dofIDarry(0);
    DofManager *node;
    double norm = 0.0;
    int nnode = domain->giveNumberOfDofManagers();

    this->timer.resumeTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
    for ( int i = 1; i <= nnode; i++ ) {
        node = domain->giveDofManager(i);
        node->computeLoadVectorAt(charVec, tStep, mode);
#ifdef __PARALLEL_MODE
        if ( node->giveParallelMode() == DofManager_shared ) {
            charVec.times( 1./( node->givePartitionsConnectivitySize() ) );
        }

#endif
        if ( charVec.isNotEmpty() ) {
            node->giveCompleteLocationArray(loc, s);
            if (node->computeM2LTransformation(R, dofIDarry)) {
                charVec.rotatedWith(R, 't');
            }
            answer.assemble(charVec, loc);
            norm += charVec.computeSquaredNorm();
        }
    }

    this->timer.pauseTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
    return norm;
}


double EngngModel :: assembleVectorFromActiveBC(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, Domain *domain)
{
    double norm = 0.0;
    int nbc = domain->giveNumberOfBoundaryConditions();

    this->timer.resumeTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
    for ( int i = 1; i <= nbc; ++i ) {
        ActiveBoundaryCondition *bc = dynamic_cast<ActiveBoundaryCondition*>(domain->giveBc(i));
        if (bc != NULL)
            norm += bc->assembleVector(answer,tStep, eid, type, mode, s, domain);
    }
    this->timer.pauseTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
    return norm;
}


double EngngModel :: assembleVectorFromElements(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                              CharType type, ValueModeType mode,
                                              const UnknownNumberingScheme &s, Domain *domain)
//
// for each element in domain
// and assembling every contribution to answer
//
{

    IntArray loc;
    FloatMatrix R;
    FloatArray charVec;
    Element *element;
    double norm = 0.0;
    int nelem = domain->giveNumberOfElements();

    this->timer.resumeTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
    for ( int i = 1; i <= nelem; i++ ) {
        element = domain->giveElement(i);
#ifdef __PARALLEL_MODE
        // skip remote elements (these are used as mirrors of remote elements on other domains
        // when nonlocal constitutive models are used. They introduction is necessary to
        // allow local averaging on domains without fine grain communication between domains).
        if ( element->giveParallelMode() == Element_remote ) {
            continue;
        }

#endif
        element->giveLocationArray(loc, eid, s);
        this->giveElementCharacteristicVector(charVec, i, type, mode, tStep, domain);
        if ( charVec.isNotEmpty() ) {
            if ( element->giveRotationMatrix(R, eid) ) {
                charVec.rotatedWith(R, 't');
            }
            answer.assemble(charVec, loc);
            norm += charVec.computeSquaredNorm();
        }
    }

    this->timer.pauseTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);
    return norm;
}


void
EngngModel :: assembleExtrapolatedForces(FloatArray &answer, TimeStep *tStep, EquationID eid, CharType type, Domain *domain)
{
    // Simply assembles contributions from each element in domain
    Element *element;
    IntArray loc;
    FloatArray charVec, delta_u;
    FloatMatrix charMatrix;
    int nelems;
    EModelDefaultEquationNumbering dn;

    answer.resize( this->giveNumberOfEquations( eid ) );
    answer.zero();

#ifdef __PARALLEL_MODE
    if ( isParallel() ) {
        exchangeRemoteElementData( RemoteElementExchangeTag  );
    }
#endif

    nelems = domain->giveNumberOfElements();
    this->timer.resumeTimer(EngngModelTimer :: EMTT_NetComputationalStepTimer);

    for ( int i = 1; i <= nelems; i++ ) {
        element = domain->giveElement(i);
#ifdef __PARALLEL_MODE
        // Skip remote elements (these are used as mirrors of remote elements on other domains
        // when nonlocal constitutive models are used. Their introduction is necessary to
        // allow local averaging on domains without fine grain communication between domains).
        if ( element->giveParallelMode() == Element_remote ) continue;
#endif
        element->giveLocationArray( loc, eid, dn );

        // Take the tangent from the previous step
        ///@todo This is not perfect. It is probably no good for viscoelastic materials, and possibly other scenarios that are rate dependent
        ///(tangent will be computed for the previous step, with whatever deltaT it had)
        element->giveCharacteristicMatrix( charMatrix, type, tStep );
        element->computeVectorOf(eid, VM_Incremental, tStep, delta_u);
        charVec.beProductOf(charMatrix, delta_u);

        ///@todo Deal with element deactivation and reactivation properly.

        answer.assemble(charVec, loc);
    }

    this->timer.pauseTimer( EngngModelTimer :: EMTT_NetComputationalStepTimer );
}


void
EngngModel :: giveElementCharacteristicMatrix(FloatMatrix &answer, int num, CharType type, TimeStep *tStep, Domain *domain)
{
    domain->giveElement(num)->giveCharacteristicMatrix(answer, type, tStep);
}


void
EngngModel :: giveElementCharacteristicVector(FloatArray &answer, int num, CharType type, ValueModeType mode, TimeStep *tStep, Domain *domain)
{
    domain->giveElement(num)->giveCharacteristicVector(answer, type, mode, tStep);
}


void
EngngModel ::  updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d)
//
// updates some component, which is used by numerical method
// to newly reached state
//
{
    _error("updateComponent: Unknown Type of component.");
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
    int j;
    Element *elem;
    Domain *domain;

    for ( int idomain = 1; idomain <= this->ndomains; idomain++ ) {
        domain = this->giveDomain(idomain);

        int nelem = domain->giveNumberOfElements();
        for ( j = 1; j <= nelem; j++ ) {
            elem = domain->giveElement(j);
#ifdef __PARALLEL_MODE
            // skip remote elements (these are used as mirrors of remote elements on other domains
            // when nonlocal constitutive models are used. They introduction is necessary to
            // allow local averaging on domains without fine grain communication between domains).
            if ( elem->giveParallelMode() == Element_remote ) {
                continue;
            }

#endif
            elem->initForNewStep();
        }
    }
}


void
EngngModel :: updateDomainLinks()
{
    this->giveExportModuleManager()->initialize();
};




contextIOResultType EngngModel :: saveContext(DataStream *stream, ContextMode mode, void *obj)
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
    int closeFlag = 0;
    FILE *file;


    OOFEM_LOG_INFO("Storing context\n");
    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, this->giveCurrentStep()->giveNumber(),
                                    this->giveCurrentStep()->giveVersion(), contextMode_write) ) {
            THROW_CIOERR(CIO_IOERR); // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    // store solution step
    if ( ( iores = giveCurrentStep()->saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // store numberOfEquations and domainNeqs array
    if ( !stream->write(& numberOfEquations, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = domainNeqs.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // store numberOfPrescribedEquations and domainNeqs array
    if ( !stream->write(& numberOfPrescribedEquations, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = domainPrescribedNeqs.storeYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // store renumber flag
    if ( !stream->write(& renumberFlag, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    for ( int idomain = 1; idomain <= this->ndomains; idomain++ ) {
        this->giveDomain(idomain)->saveContext(stream, mode, obj);
    }


    // store nMethod
    NumericalMethod *nmethod = this->giveNumericalMethod( this->giveMetaStep( giveCurrentStep()->giveMetaStepNumber() ) );
    if ( nmethod ) {
        if ( ( iores = nmethod->saveContext(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }


    if ( closeFlag ) {
        fclose(file);
        delete(stream);
        stream = NULL;
    }                                                         // ensure consistent records

    return CIO_OK;
}


contextIOResultType EngngModel :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
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
//
// WARNING obj is cast into int pointer  to time step for which to seek Context.
//
{
    contextIOResultType iores;
    int serNum, closeFlag = 0, istep, iversion;
    bool domainUpdated = false;
    Domain *domain;
    FILE *file;

    this->resolveCorrespondingStepNumber(istep, iversion, obj);
    OOFEM_LOG_RELEVANT("Restoring context for time step %d.%d\n", istep, iversion);

    if ( stream == NULL ) {
        if ( !this->giveContextFile(& file, istep, iversion, contextMode_read) ) {
            THROW_CIOERR(CIO_IOERR);                                                              // override
        }

        stream = new FileDataStream(file);
        closeFlag = 1;
    }

    // restore solution step
    if ( currentStep == NULL ) {
        currentStep = new TimeStep(istep, this, 0, 0., 0., 0);
    }

    if ( ( iores = currentStep->restoreContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // this->updateAttributes (currentStep);

    int pmstep = currentStep->giveMetaStepNumber();
    if ( nMetaSteps ) {
        if ( !this->giveMetaStep(pmstep)->isStepValid(istep - 1) ) {
            pmstep--;
        }
    }

    if ( previousStep ) {
        delete previousStep;
    }

    previousStep = new TimeStep(istep - 1, this, pmstep, currentStep->giveTargetTime ( ) - currentStep->giveTimeIncrement(),
                                currentStep->giveTimeIncrement(), currentStep->giveSolutionStateCounter() - 1);

    // restore numberOfEquations and domainNeqs array
    if ( !stream->read(& numberOfEquations, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = domainNeqs.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // restore numberOfPrescribedEquations and domainNeqs array
    if ( !stream->read(& numberOfPrescribedEquations, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    if ( ( iores = domainPrescribedNeqs.restoreYourself(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    // restore renumber flag
    if ( !stream->read(& renumberFlag, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    for ( int idomain = 1; idomain <= this->ndomains; idomain++ ) {
        domain = this->giveDomain(idomain);
        serNum = domain->giveSerialNumber();
        domain->restoreContext(stream, mode, obj);
        if ( serNum != domain->giveSerialNumber() ) {
            domainUpdated = true;
        }
    }

    // restore nMethod
    NumericalMethod *nmethod = this->giveNumericalMethod( this->giveCurrentMetaStep() );
    if ( nmethod ) {
        if ( ( iores = nmethod->restoreContext(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    this->updateDomainLinks();
    this->updateAttributes( this->giveCurrentMetaStep() );
    this->initStepIncrements();

    if ( closeFlag ) {
        fclose(file);
        delete stream;
        stream = NULL;
    }                                                           // ensure consistent records

    return CIO_OK;
}


void
EngngModel :: resolveCorrespondingStepNumber(int &istep, int &iversion, void *obj)
{
    //
    // returns corresponding step number
    //
    if ( obj == NULL ) {
        istep = 1;
        iversion = 0;
        return;
    }

    istep = * ( int * ) obj;
    iversion = * ( ( ( int * ) obj ) + 1 );

    if ( istep > this->giveNumberOfSteps() ) {
        istep = this->giveNumberOfSteps();
    }

    if ( istep <= 0 ) {
        istep = 1;
    }
}


MetaStep *
EngngModel :: giveCurrentMetaStep()
{
    return this->giveMetaStep( this->giveCurrentStep()->giveMetaStepNumber() );
}


int
EngngModel :: giveContextFile(FILE **contextFile, int stepNumber, int stepVersion, ContextFileMode cmode, int errLevel)
//
//
// assigns context file of given step number to stream
// returns nonzero on success
//
{
    std::string fname = this->coreOutputFileName;
    char fext[100];
    sprintf(fext, ".%d.%d.osf", stepNumber, stepVersion);
    fname += fext;

    if ( cmode ==  contextMode_read ) {
        * contextFile = fopen(fname.c_str(), "rb"); // open for reading
    } else {
        * contextFile = fopen(fname.c_str(), "wb"); // open for writing,
    }

    //  rewind (*contextFile); // seek the beginning
    // // overwrite if exist
    // else *contextFile = fopen(fname,"r+"); // open for reading and writing

    if ( ( * contextFile == NULL ) && errLevel > 0 ) {
        _error2("giveContextFile : can't open %s", fname.c_str());
    }

    return 1;
}

bool
EngngModel :: testContextFile(int stepNumber, int stepVersion)
{
    std::string fname = this->coreOutputFileName;
    char fext[100];
    sprintf(fext, ".%d.%d.osf", stepNumber, stepVersion);
    fname.append(fext);

#ifdef HAVE_ACCESS
    if ( access(fname.c_str(), R_OK) == 0 ) {
        return true;
    } else {
        return false;
    }

#elif _MSC_VER
    if ( _access(fname.c_str(), 4) == 0 ) {
        return true;
    } else {
        return false;
    }

#else
    return true;

#endif
}

DataReader *
EngngModel :: GiveDomainDataReader(int domainNum, int domainSerNum, ContextFileMode cmode)
//
//
// returns domain i/o file
// returns nonzero on success
//
{
    std::string fname = this->coreOutputFileName;
    char fext [ 100 ];
    sprintf(fext, ".domain.%d.%d.din", domainNum, domainSerNum);
    fname += fext;

    DataReader *dr;

    if ( ( dr = new OOFEMTXTDataReader(fname.c_str()) ) == NULL ) {
        _error("Creation of DataReader failed");
    }

    return dr;
}

void
EngngModel :: error(const char *file, int line, const char *format, ...) const
{
    char buffer [ MAX_ERROR_MSG_LENGTH ];
    va_list args;

    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);
#ifdef __PARALLEL_MODE
    __OOFEM_ERROR4(file, line, "Class: %s, Rank: %d\n%s", giveClassName(), rank, buffer);
#else
    __OOFEM_ERROR3(file, line, "Class: %s\n%s", giveClassName(), buffer);
#endif
}


void EngngModel :: warning(const char *file, int line, const char *format, ...) const
//
// this function handles error reporting
// prints errorMsg enriched by ClasName and Number
// to the standard error stream
{
    char buffer [ MAX_ERROR_MSG_LENGTH ];
    va_list args;

    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);
#ifdef __PARALLEL_MODE
    __OOFEM_WARNING4(file, line, "Class: %s, Rank: %d\n%s", giveClassName(), rank, buffer);
#else
    __OOFEM_WARNING3(file, line, "Class: %s\n%s", giveClassName(), buffer);
#endif
}

Domain *
EngngModel :: giveDomain(int i)
{
    if ( ( i > 0 ) && ( i <= this->ndomains ) ) {
        return this->domainList->at(i);
    } else {
        _error("giveDomain: Undefined domain");
    }

    return NULL;
}

XfemManager *
EngngModel :: giveXfemManager(int i)
{
    if ( xfemManagerList->includes(i) ) {
        return this->xfemManagerList->at(i);
    } else {
        _error2("giveXfemManager: undefined xfem manager (%d)", i);
        return NULL; // return NULL to prevent compiler warnings
    }
}

bool
EngngModel :: hasXfemManager(int i)
{
    return xfemManagerList->includes(i);
}


#ifdef __PETSC_MODULE
PetscContext *
EngngModel :: givePetscContext(int i, EquationID eid)
{
    if ( ( i > 0 ) && ( i <= this->ndomains ) ) {
        if  ( i > petscContextList->giveSize() ) {
            _error("givePetscContext: petsc context not initialized for this problem");
        }

        return this->petscContextList->at(i);
    } else {
        _error2("givePetscContext: Undefined domain index %d ", i);
    }

    return NULL;
}

void
EngngModel :: initPetscContexts() { }
#endif



MetaStep *
EngngModel :: giveMetaStep(int i)
{
    if ( ( i > 0 ) && ( i <= this->nMetaSteps ) ) {
        return this->metaStepList->at(i);
    } else {
        _error2("giveMetaStep: undefined metaStep (%d)", i);
    }

    return NULL;
}

FILE *
EngngModel :: giveOutputStream()
// Returns an output stream on the data file of the receiver.
{
    if ( !outputStream ) {
#ifdef _MSC_VER

        char *tmp = tmpnam(NULL);
        _warning2("giveOutputStream: using default output stream %s", tmp);
        outputStream = fopen(tmp, "w");
#else

        char sfn[] = "oofem.out.XXXXXX";
        int fd = -1;
        FILE *sfp;

        if ( ( fd = mkstemp(sfn) ) == -1 ||
            ( sfp = fdopen(fd, "w+") ) == NULL ) {
            if ( fd != -1 ) {
                unlink(sfn);
                close(fd);
            }

            OOFEM_ERROR2("Failed to create temporary file %s\n", sfn);
            return ( NULL );
        }

        outputStream = sfp;
#endif
    }

    return outputStream;
}

void
EngngModel :: terminateAnalysis()
{
    double tsec;
    int nsec = 0, nmin = 0, nhrs = 0;
    FILE *out = this->giveOutputStream();
    time_t endTime = getTime();
    this->timer.stopTimer(EngngModelTimer :: EMTT_AnalysisTimer);


    fprintf( out, "\nFinishing analysis on: %s\n", ctime(& endTime) );
    // compute real time consumed
    tsec = this->timer.getWtime(EngngModelTimer :: EMTT_AnalysisTimer);
    this->timer.convert2HMS(nhrs, nmin, nsec, tsec);
    fprintf(out, "Real time consumed: %03dh:%02dm:%02ds\n", nhrs, nmin, nsec);
    LOG_FORCED_MSG(oofem_logger, "\n\nANALYSIS FINISHED\n\n\n");
    LOG_FORCED_MSG(oofem_logger, "Real time consumed: %03dh:%02dm:%02ds\n", nhrs, nmin, nsec);
    // compute processor time used by the program
    // nsec = (endClock - startClock) / CLOCKS_PER_SEC;
    tsec = this->timer.getUtime(EngngModelTimer :: EMTT_AnalysisTimer);
    this->timer.convert2HMS(nhrs, nmin, nsec, tsec);
    fprintf(out, "User time consumed: %03dh:%02dm:%02ds\n\n\n", nhrs, nmin, nsec);
    LOG_FORCED_MSG(oofem_logger, "User time consumed: %03dh:%02dm:%02ds\n", nhrs, nmin, nsec);

    exportModuleManager->terminate();
}

int
EngngModel :: checkProblemConsistency()
{
    int result = 1;

    result &= this->checkConsistency();

    int ndomains = this->giveNumberOfDomains();
    for ( int i = 1; i <= ndomains; i++ ) {
        result &= this->giveDomain(i)->checkConsistency();
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
EngngModel :: init()
{
    initModuleManager->doInit();
}


#ifdef __PARALLEL_MODE
void
EngngModel :: initParallel()
{
 #ifdef __USE_MPI
    int len;
    MPI_Get_processor_name(processor_name, & len);
    MPI_Comm_rank(MPI_COMM_WORLD, & this->rank);
    MPI_Comm_size(MPI_COMM_WORLD, & numProcs);
 #endif
 #ifdef __VERBOSE_PARALLEL
    OOFEM_LOG_RELEVANT("[%d/%d] Running on %s\n", rank, numProcs, processor_name);
 #endif
}

#endif

#ifdef __OOFEG
void EngngModel :: drawYourself(oofegGraphicContext &context) {
    this->giveDomain( context.getActiveDomain() )->drawYourself(context);
}

void EngngModel :: drawElements(oofegGraphicContext &context) {
    this->giveDomain( context.getActiveDomain() )->drawElements(context);
}

void EngngModel :: drawNodes(oofegGraphicContext &context) {
    this->giveDomain( context.getActiveDomain() )->drawNodes(context);
}

#endif


#ifdef __PARALLEL_MODE
void
EngngModel :: balanceLoad(TimeStep *atTime)
{
    LoadBalancerMonitor :: LoadBalancerDecisionType _d;
    this->giveLoadBalancerMonitor();
    this->giveLoadBalancer();
    EModelDefaultEquationNumbering dn;

    //print statistics for current step
    lb->printStatistics();

    if ( atTime->isNotTheLastStep() ) {
        _d = lbm->decide(atTime);
        if ( ( _d == LoadBalancerMonitor :: LBD_RECOVER ) ||
            ( ( atTime->isTheFirstStep() ) && force_load_rebalance_in_first_step ) ) {
            this->timer.startTimer(EngngModelTimer :: EMTT_LoadBalancingTimer);

            // determine nwe partitioning
            lb->calculateLoadTransfer();
            // pack e-model solution data into dof dictionaries
            this->packMigratingData(atTime);
            // migrate data
            lb->migrateLoad( this->giveDomain(1) );
            // renumber itself
            this->forceEquationNumbering();
 #ifdef __VERBOSE_PARALLEL
            // debug print
            int i, j, nnodes = giveDomain(1)->giveNumberOfDofManagers();
            int myrank = this->giveRank();
            fprintf(stderr, "\n[%d] Nodal Table\n", myrank);
            for ( i = 1; i <= nnodes; i++ ) {
                if ( giveDomain(1)->giveDofManager(i)->giveParallelMode() == DofManager_local ) {
                    fprintf( stderr, "[%d]: %5d[%d] local ", myrank, i, giveDomain(1)->giveDofManager(i)->giveGlobalNumber() );
                } else if ( giveDomain(1)->giveDofManager(i)->giveParallelMode() == DofManager_shared ) {
                    fprintf( stderr, "[%d]: %5d[%d] shared ", myrank, i, giveDomain(1)->giveDofManager(i)->giveGlobalNumber() );
                }

                for ( j = 1; j <= giveDomain(1)->giveDofManager(i)->giveNumberOfDofs(); j++ ) {
                    fprintf( stderr, "(%d)", giveDomain(1)->giveDofManager(i)->giveDof(j)->giveEquationNumber(dn) );
                }

                fprintf(stderr, "\n");
            }

 #endif
            // unpack (restore) e-model solution data from dof dictionaries
            this->unpackMigratingData(atTime);

            this->timer.stopTimer(EngngModelTimer :: EMTT_LoadBalancingTimer);
            double _steptime = this->timer.getUtime(EngngModelTimer :: EMTT_LoadBalancingTimer);
            OOFEM_LOG_INFO("[%d] EngngModel info: user time consumed by load rebalancing %.2fs\n",
                           this->giveRank(), _steptime);
        }
    }
}
#endif
} // end namespace oofem
