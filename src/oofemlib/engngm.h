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

#ifndef engngm_h
#define engngm_h

#include "inputrecord.h"

#include "alist.h"
#include "intarray.h"
#include "fieldmanager.h"
#include "timer.h"
#include "chartype.h"
#include "classtype.h"
#include "unknowntype.h"
#include "varscaletype.h"
#include "equationid.h"
#include "numericalcmpn.h"
#include "valuemodetype.h"
#include "problemmode.h"
#include "fmode.h"
#include "dofiditem.h"
#include "contextoutputmode.h"
#include "contextfilemode.h"
#include "contextioresulttype.h"
#include "xfemmanager.h"
#include "unknownnumberingscheme.h"

#ifdef __PARALLEL_MODE
 #include "parallel.h"
 #include "problemcommunicatormode.h"
 #include "loadbalancer.h"
#endif

#ifdef __PETSC_MODULE
 #include "petsccontext.h"
 #include "petscordering.h"
#endif

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

#ifndef __MAKEDEPEND
 #include <string>
#endif


namespace oofem {

class Domain;
class TimeStep;
class Dof;
class DofManager;
class DataReader;
class DataStream;
class ErrorEstimator;
class MetaStep;
class MaterialInterface;
class SparseMtrx;
class NumericalMethod;
class XfemManager;
class InitModuleManager;
class ExportModuleManager;
class FloatMatrix;
class FloatArray;

#ifdef __PARALLEL_MODE
class ProblemCommunicator;
#endif


/**
 * Class EngngModelContext represents a context, which is shared by all problem engng sub-models.
 * In principle every problem (represented by the EngngModel class) can be made part of more
 * complex problem, providing only part of the solution, possibly depending on the results of
 * other sub-problems. Typical example is staggered heat and structural analysis.
 * The context provides the common (shared) resources within the problem. The context is created by the
 * master problem (the class representing staggered problem, for example). The subproblems are
 * then created in so-called maintained (or slave) mode and they request the context from master.
 */
class EngngModelContext
{
protected:
    /// Common fieldManager providing shared field register for the problem.
    FieldManager fieldManager;
public:

    EngngModelContext() { }
    FieldManager *giveFieldManager() { return & ( this->fieldManager ); }
};


/**
 * Abstract base class representing the "problem" under consideration.
 * The engineering model describes the problem and type of analysis to be done.
 * It "knows" the type and form of governing equation, knows how to assemble this
 * problem from local element contributions and uses suitable instance of
 * numerical method class to solve problem. It possesses and manages one or more problem domains.
 * Concept of time steps is introduced. For problems discretized in time the introduction of
 * time step is natural. In other cases, time steps can represent for example
 * load increments or different load cases.
 *
 * The solution steps are grouped together into so called meta steps. The meta step can be thought as an sequence of
 * solution steps, with same set of attributes used to drive behavior of engng model.
 * For each metastep, the engng model typically updates its control attributes according to
 * metastep engng attributes (see initMetaStepAttributes and updateAttributes services)
 * and creates the solution steps accordingly. This allows to switch to
 * different time increment, different solution control, etc. If no metastep is specified, the engng model
 * creates default one for all required solution steps. There are two services, where attributes are updated, the first one,
 * used for those attributes, which do not vary during solution of problem are set in initializeForm service.
 * The second service is updateAttributes, where the attributes allowed to change (with metastep validity) are updated.
 * If no metastep is introduced, default one is created (with attributes set to engng model init record).
 * Then there is no difference, whether attributes are read in initializeFrom or updateAttributes, but
 * preferred scheme is to read all attributes in initializeFrom and left updateAttributes service empty.
 *
 * The basic EngngModel tasks are following
 * - Assembling governing equations by summing contributions from problem domains (typically from nodes and elements),
 * - Solving the problem described by governing equation(s) using suitable instance of
 *   numerical method. This requires interfacing numericalMethod characteristic elements
 *   with components in governing equation
 *   EngngModel must map each component of governing
 *   equation(s) (which has physical meaning) to corresponding numerical component of Numerical
 *   method. This mapping between physical components to independent numerical components
 *   (understood by numerical method) is important, because it allows numerical method to be used by
 *   many EngngModel with different meaning of particular components.
 * - Returning unknown values (according to requested type and mode). Used by Dofs to
 *   access their corresponding unknowns.
 * - Terminating time step by updating nodal and element values (including integration points update).
 * - Updating dofs unknowns dictionaries if model supports changes of static system (see Dof class
 *   documentation for detailed explanation). In general if static system changes are not supported,
 *   when dof are requested for unknowns, they use their associate equation number to ask EngngModel
 *   for this unknown. Unknowns are therefore stored in EngngModel and are requested by dofs.
 *   On the other hand, when static system changes are supported, the equation numbers of dofs
 *   can vary during solution. Therefore, so called unknowns dictionary at dof level are introduced.
 *   All unknowns are stored on dof level and dofs will use in such case their own dictionaries
 *   instead of requesting EngngModel. The EngngModel is fully responsible to update this
 *   dictionary for each dof with all necessary unknowns (see updateDofUnknownsDictionary function).
 */
class EngngModel
{
public:
#ifdef __PARALLEL_MODE
    enum EngngModel_UpdateMode { EngngModel_Unknown_Mode, EngngModel_SUMM_Mode, EngngModel_SET_Mode };
    enum EngngModelCommType { PC_default, PC_nonlocal };
#endif

    /**
     * Means to choose methods for finding a good initial guess.
     * This is ad-hoc methods, often problem-specific for obtaining a point from which Newton iterations work.
     * Only nonlinear analysis needs to worry about this.
     */
    enum InitialGuess {
        IG_None = 0, ///< No special treatment for new iterations. Probably means ending up using @f$ {}^{n+1}x = {}^{n}x @f$ for all free dofs.
        IG_Tangent = 1, ///< Solves an approximated tangent problem from the last iteration. Useful for changing Dirichlet boundary conditions.
        //IG_Extrapolated = 2, ///< Assumes constant increment extrapolating @f$ {}^{n+1}x = {}^{n}x + \Delta t\delta{x}'@f$, where @f$ \delta x' = ({}^{n}x - {}^{n-1}x)/{}^{n}Delta t@f$.
    };

protected:
    /// Number of receiver domains.
    int ndomains;
    /// List of problem domains.
    AList< Domain > *domainList;
    /// Total number of time steps.
    int numberOfSteps;
    /// Total number of equation in current time step.
    int numberOfEquations;
    /// Total number or prescribed equations in current time step.
    int numberOfPrescribedEquations;
    /// Number of equations per domain.
    IntArray domainNeqs;
    /// Number of prescribed equations per domain.
    IntArray domainPrescribedNeqs;
    /// Renumbering flag (renumbers equations after each step, necessary if Dirichlet BCs change).
    int renumberFlag;
    /// Profile optimized numbering flag (using Sloan's algorithm).
    int profileOpt;
    /// Equation numbering completed flag.
    int equationNumberingCompleted;
    /// Number of meta steps.
    int nMetaSteps;
    /// List of problem metasteps.
    AList< MetaStep > *metaStepList;
    /// Solution step when IC (initial conditions) apply.
    TimeStep *stepWhenIcApply;
    /// Current time step.
    TimeStep *currentStep;
    /// Previous time step.
    TimeStep *previousStep;
    /// Receivers id.
    int number;

    //DefaultEquationNumbering
    EModelDefaultEquationNumbering defaultNumberingScheme;
    // DefaultPrescribedEquationNumbering
    EModelDefaultPrescribedEquationNumbering defaultPrescribedNumberingScheme;

    /// Path to output stream.
    std::string dataOutputFileName;
    /// String with core output file name
    std::string coreOutputFileName;
    /// Output stream.
    FILE *outputStream;
    /// Domain context output mode.
    ContextOutputMode contextOutputMode;
    int contextOutputStep;

    /// Export module manager.
    ExportModuleManager *exportModuleManager;
    /// Initialization module manager.
    InitModuleManager *initModuleManager;

    /// Domain mode.
    problemMode pMode;
    /// Multiscale mode.
    problemScale pScale;
    /// Solution start time.
    time_t startTime;
    // initial value of processor time used by program
    // clock_t startClock;

    /// Master e-model; if defined receiver is in maintained (slave) mode.
    EngngModel *master;
    /// Context.
    EngngModelContext *context;
    /// E-model timer.
    EngngModelTimer timer;
    /// Flag indicating that the receiver runs in parallel.
    int parallelFlag;
    /// List of Xfemmanagers.
    AList< XfemManager > *xfemManagerList;
    /// Number of Xfemmanagers.
    int nxfemman;
    /// Type of non linear formulation (total or updated formulation).
    enum fMode nonLinFormulation;
    /// Error estimator. Useful for adaptivity, or simply printing errors output.
    ErrorEstimator *defaultErrEstimator;

#ifdef __PARALLEL_MODE
    /// Domain rank in a group of collaborating processes (0..groupSize-1).
    int rank;
    /// Total number of collaborating processes.
    int numProcs;
    /// Processor name.
    char processor_name [ PROCESSOR_NAME_LENGTH ];
    /// Communicator mode. Determines current strategy used.
    ProblemCommunicatorMode commMode;


    /**@name Load balancing attributes */
    //@{
    /// Load Balancer.
    LoadBalancer *lb;
    LoadBalancerMonitor *lbm;
    /// If set to true, load balancing is active.
    bool loadBalancingFlag;
    /// Debug flag forcing load balancing after first step.
    bool force_load_rebalance_in_first_step;
    //@}

#endif // __PARALLEL_MODE

#ifdef __PETSC_MODULE
    /// List where PETSc contexts are stored.
    AList< PetscContext > *petscContextList;
#endif

public:
    /**
     * Constructor. Creates Engng model with number i.
     */
    EngngModel(int i, EngngModel *_master = NULL);
    /// Destructor.
    virtual ~EngngModel();

    /**
     * Service for accessing particular problem domain.
     * Generates error if no such domain is defined.
     * @param n Pointer to n-th domain is returned.
     * @return Domain number n.
     */
    Domain *giveDomain(int n);
    /**
     * Sets i-th domain of receiver. Given domain is assumed to be owned (and deleted) by receiver.
     * The old domain, if defined, will be deleted.
     * @param i Domain index.
     * @param ptr Pointer to valid domain instance.
     */
    void setDomain (int i, Domain *ptr);
    /// Returns number of domains in problem.
    int giveNumberOfDomains() { return ndomains; }
    /** Service for accessing ErrorEstimator corresponding to particular domain */
    virtual ErrorEstimator *giveDomainErrorEstimator(int n) { return defaultErrEstimator; }
    /** Returns material interface representation for given domain */
    virtual MaterialInterface *giveMaterialInterface(int n) { return NULL; }
    /** Returns XfemManager at a particular position */
    XfemManager *giveXfemManager(int n);
    /** Return true if XfemManager at a particular position is available */
    bool hasXfemManager(int n);
    void setNumberOfEquations(int id, int neq) {
        numberOfEquations = neq;
        domainNeqs.at(id) = neq;
    }
    // input / output
    /// Returns file descriptor of output file.
    FILE *giveOutputStream();
    /**
     * Returns base output file name
     * to which extensions, like .out .vtu .osf should be added.
     * In current implementation, output file name is simply returned.
     * @param dest Path and base file name will be copied into the array pointed to by dest.
     * @param n Not more than n bytes of src are copied
     */
    std::string giveOutputBaseFileName() { return dataOutputFileName; }
    /**
     * Sets the base output file name.
     * @see giveOutputBaseFileName
     * @param src New output file name.
     */
    void letOutputBaseFileNameBe(const std::string &src) { dataOutputFileName = src; }
    /**
     * Returns domain context output mode.
     */
    ContextOutputMode giveContextOutputMode() { return contextOutputMode; }
    /**
     * Returns domain context output step.
     */
    int giveContextOutputStep() { return contextOutputStep; }
    /**
     * Sets context output mode of receiver.
     * @param contextMode domain context mode.
     */
    void setContextOutputMode(ContextOutputMode contextMode)
    { contextOutputMode = contextMode; }
    /**
     * Sets user defined context output mode (it sets contextOutputMode to contextOutputMode),
     * setting contextOutputStep to given value.
     * @param cStep new context output step
     */
    void setUDContextOutputMode(int cStep)
    {
        contextOutputMode = COM_UserDefined;
        contextOutputStep = cStep;
    }
    /**
     * Sets domain mode to given mode.
     * @param pmode Problem mode.
     */
    void setProblemMode(problemMode pmode) { pMode = pmode; }
    /**
     * Sets the problem to run in parallel (or not).
     * @param parallelFlag Determines parallel mode.
     */
    void setParallelMode(bool parallelFlag);
    /// Returns domain mode.
    problemMode giveProblemMode() { return pMode; }
    /**
     * Sets scale in multiscale simulation.
     * @param pscale Problem scale.
     */
    void setProblemScale(problemScale pscale) { pScale = pscale; }
    /// Returns scale in multiscale simulation
    problemScale giveProblemScale() { return pScale; }
    /// Sets the renumber flag to true.
    virtual void setRenumberFlag() { this->renumberFlag = true; }
    /// Sets the renumber flag to false.
    virtual void resetRenumberFlag() { this->renumberFlag = false; }

    /**
     * Performs analysis termination after finishing analysis.
     */
    void terminateAnalysis();

    // solving
    /**
     * Starts solution process. Implementation should invoke for each time step
     * solveYourselfAt function with time step as parameter. Time steps are created
     * using giveNextStep function (this will set current time step to newly created,
     * and updates previous step).
     */
    virtual void solveYourself();
    /**
     * Solves problem for given time step. Should assemble characteristic matrices and vectors
     * if necessary and solve problem using appropriate numerical method. After finishing solution,
     * this->updateYourself function for updating solution state and then this->terminate
     * function (for updating nodal and element values) should be called.
     */
    virtual void solveYourselfAt(TimeStep *tStep) { }
    /**
     * Terminates the solution of time step. Default implementation calls prinOutput() service and if specified,
     * context of whole domain is stored and output for given time step is printed.
     */
    virtual void terminate(TimeStep *tStep);
    /**
     * Prints the ouput of the solution step (using virtual this->printOutputAtservice)
     * to the stream detemined using this->giveOutputStream() method
     * and calls exportModuleManager to do output.
     */
    virtual void doStepOutput(TimeStep *tStep);
    /**
     * Saves context of given solution step, if required (determined using this->giveContextOutputMode() method).
     */
    void saveStepContext(TimeStep *tStep);
    /**
     * Updates internal state after finishing time step. (for example total values may be
     * updated according to previously solved increments). Then element values are also updated
     * (together with related integration points and material statuses).
     */
    virtual void updateYourself(TimeStep *tStep);
    /**
     * Provides the opportunity to initialize state variables stored in element
     * integration points according to
     * initial conditions using function initializeYourself() on element level.
     * Should be called when current time step is time step when IC will apply
     * (see EngngModel::giveNumberOfTimeStepWhenIcApply)
     * somewhere from solveYourselfAt function). Implementation must be provided.
     * Default implementation is empty.
     */
    virtual void initializeYourself(TimeStep *tStep) { }
    /**
     * Initializes the newly generated discretization state according to previous solution.
     * This process should typically include restoring old solution, instanciating newly
     * generated domain(s) and by mapping procedure.
     */
    virtual int initializeAdaptive(int stepNumber) { return 0; }

    /**
     * Returns total number of equations in active (current time step) time step.
     * The UnknownType parameter allows to distinguish between several possible governing equations, that
     * can be numbered separately.
     */
    virtual int giveNumberOfEquations(EquationID eid);
    /**
     * Returns total number of prescribed equations in active (current time step) time step.
     * The UnknownType parameter allows to distinguish between several possible governing equations, that
     * can be numbered separately.
     */
    virtual int giveNumberOfPrescribedEquations(EquationID eid);
    /**
     * Returns number of equations for given domain in active (current time step) time step.
     * The EquationID parameter allows to distinguish between several possible governing equations, that
     * can be numbered separately.
     */
    virtual int giveNumberOfDomainEquations(int di, EquationID eid);
    /**
     * Returns number of prescribed equations for given domain in active (current time step) time step.
     * The UnknownType parameter allows to distinguish between several possible governing equations, that
     * can be numbered separately.
     */
    virtual int giveNumberOfPrescribedDomainEquations(int di, EquationID eid);
    //virtual IntArray* GiveBanWidthVector ();


    // management components
    /**
     * Provides backward mapping between numerical component and characteristic
     * component on EngngModel level.
     */
    virtual CharType giveTypeOfComponent(NumericalCmpn) { return UnknownCharType; }
    /**
     * Returns requested unknown. Unknown at give time step is characterized by its type and mode
     * and by its equation number. This function is used by Dofs, when they are requested for
     * their associated unknowns.
     * @see Dof::giveUnknown
     */
    virtual double giveUnknownComponent(EquationID, ValueModeType, TimeStep *, Domain *, Dof *) { return 0.0; }
    virtual double giveUnknownComponent(UnknownType, ValueModeType, TimeStep *, Domain *, Dof *) { return 0.0; }

#ifdef __PARALLEL_MODE
    /**
     * Updates unknown. Unknown at give time step is characterized by its type and mode
     * and by its equation number. This function is used by Dofs, when they are requested for
     * their update of associated unknowns. Declared only in __PARALLEL_MODE
     * @see Dof::giveUnknown method
     */
    virtual void updateUnknownComponent(EquationID, ValueModeType, TimeStep *, int,
                                        double, EngngModel_UpdateMode) { return; }

    virtual ProblemCommunicator *giveProblemCommunicator(EngngModelCommType t) { return NULL; }
#endif
    /**
     * Initializes whole problem according to its description stored in inputStream.
     * Prints header, opens the outFileName, instanciate itself the receiver using
     * using virtual initializeFrom service and instanciates all problem domains.
     */
    virtual int instanciateYourself(DataReader *dr, InputRecord *ir, const char *outFileName, const char *desc);
    /**
     * Initializes receiver according to object description in input reader.
     * InitString can be imagined as data record in component database
     * belonging to receiver. Receiver may use value-name extracting functions
     * to extract particular field from record.*/
    virtual IRResultType initializeFrom(InputRecord *ir);
    /// Instanciate problem domains by calling their instanciateYourself() service
    int instanciateDomains(DataReader *dr);
    /// Instanciate problem meta steps by calling their instanciateYourself() service
    int instanciateMetaSteps(DataReader *dr);
    /// Instanciate default metastep, if nmsteps is zero
    int instanciateDefaultMetaStep(InputRecord *ir);

    /**
     * Update receiver attributes according to step metaStep attributes.
     * Allows the certain parameters or attributes to be updated for particular metastep.
     * The metastep provides the attributes record, from which the corresponding attributes can
     * be read. The service takes a MetaStep parameter. It is recommended, to implement this service in such way, that multiple calls
     * for steps belonging to same MetaStep does not change response.
     * The default implementation updates the numerical method attributes.
     * @param mStep Meta step.
     */
    virtual void updateAttributes(MetaStep *mStep);
    /**
     * Update e-model attributes attributes according to step metaStep attributes.
     * Calls updateAttributes. At the end the meta step input reader finish() service
     * is called in order to allow for unread attribute check.
     */
    void initMetaStepAttributes(MetaStep *mStep);
    /**
     * Stores the state of model to output stream. Stores not only the receiver state,
     * but also same function is invoked for all DofManagers and Elements in associated
     * domain. Note that by storing element context also contexts of all associated
     * integration points (and material statuses) are stored.
     * Stored context is associated with current time step. One time step can have only
     * one associated context. Multiple call to saveContext within same time step
     * override previously saved context for this step.
     * By default the stream parameter is used to store data and is not closed.
     * If stream is NULL, new file descriptor is created and this must be also closed at the end.
     * @param stream Context stream. If NULL then new file descriptor will be opened and closed
     * at the end else the stream given as parameter will be used and not closed at the end.
     * @param mode Determines amount of info in stream.
     * @param obj Void pointer to an int array containing two values:time step number and
     * version of a context file to be restored.
     * @return contextIOResultType.
     * @exception ContextIOERR If error encountered.
     */
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores the state of model from output stream. Restores not only the receiver state,
     * but also same function is invoked for all DofManagers and Elements in associated
     * domain. Note that by restoring element context also contexts of all associated
     * integration points (and material statuses) are restored.
     * Each context is associated with unique time step. Only one context per time step is
     * allowed. Restore context function will restore such context, which is related
     * (through its step number) to time step number and version given in obj parameter.
     * Restoring context will change current time step in order to correspond to newly restored
     * context.
     * @param stream Context file.
     * @param mode Determines amount of info in stream.
     * @param obj Void pointer to an int array containing two values:time step number and
     * version of a context file to be restored.
     * @return contextIOResultType.
     * @exception ContextIOERR exception if error encountered.
     */
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Updates domain links after the domains of receiver have changed. Used mainly after
     * restoring context - the domains may change and this service is then used
     * to update domain variables in all components belonging to receiver
     * like error estimators, solvers, etc, having domains as attributes.
     */
    virtual void updateDomainLinks();
    void resolveCorrespondingStepNumber(int &, int &, void *obj);
    /// Returns current meta step.
    MetaStep *giveCurrentMetaStep();
    /// Returns current time step.
    TimeStep *giveCurrentStep() { if ( master ) { return master->giveCurrentStep(); } else { return currentStep; } }
    /// Returns previous time step.
    TimeStep *givePreviousStep() { if ( master ) { return master->givePreviousStep(); } else { return previousStep; } }
    /// Returns next time step (next to current step) of receiver.
    virtual TimeStep *giveNextStep() { return NULL; }
    /// Returns the solution step when Initial Conditions (IC) apply.
    virtual TimeStep *giveSolutionStepWhenIcApply() { if ( master ) { return master->giveCurrentStep(); } else { return stepWhenIcApply; } }
    /// Returns number of first time step used by receiver.
    virtual int giveNumberOfFirstStep() { if ( master ) { return master->giveNumberOfFirstStep(); } else { return 1; } }
    /// Return number of meta steps.
    int giveNumberOfMetaSteps() { return nMetaSteps; }
    /// Returns the i-th meta step.
    MetaStep *giveMetaStep(int i);
    /// Returns total number of steps.
    int giveNumberOfSteps() { if ( master ) { return master->giveNumberOfSteps(); } else { return numberOfSteps; } }
    /// Returns end of time interest (time corresponding to end of time integration).
    virtual double giveEndOfTimeOfInterest() { return 0.; }
    /// Returns the time step number, when initial conditions should apply.
    virtual int giveNumberOfTimeStepWhenIcApply() {
        if ( master ) { return master->giveNumberOfTimeStepWhenIcApply(); } else { return 0; } }
    /// Returns reference to receiver's numerical method.
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep) { return NULL; }
    /// Returns receiver's export module manager.
    ExportModuleManager *giveExportModuleManager() { return exportModuleManager; }
    /// Returns reference to receiver timer (EngngModelTimer).
    EngngModelTimer *giveTimer() { return & timer; }

    /**
     * Increases number of equations of receiver's domain and returns newly created equation number.
     * Used mainly by DofManagers to allocate their corresponding equation number if it
     * is not currently allocated.
     * The DofIDItem parameter allows to distinguish between several possible governing equations, that
     * can be numbered separately.
     */
    virtual int giveNewEquationNumber(int domain, DofIDItem) { return ++domainNeqs.at(domain); }
    /**
     * Increases number of prescribed equations of receiver's domain and returns newly created equation number.
     * Used mainly by DofManagers to allocate their corresponding equation number if it
     * is not currently allocated.
     * The DofIDItem parameter allows to distinguish between several possible governing equations, that
     * can be numbered separately.
     */
    virtual int giveNewPrescribedEquationNumber(int domain, DofIDItem) { return ++domainPrescribedNeqs.at(domain); }
    /**
     * Assigns context file-descriptor for given step number to stream.
     * Returns nonzero on success.
     * @param contextFile Assigned file descriptor.
     * @param stepNumber Solution step number to store/restore.
     * @param stepVersion Version of step.
     * @param cmode Determines the i/o mode of context file.
     * @param errLevel Determines the amount of warning messages if errors are encountered, level 0 no warnings reported.
     */
    int giveContextFile(FILE **contextFile, int stepNumber, int stepVersion,
                                     ContextFileMode cmode, int errLevel = 1);
    /** Returns true if context file for given step and version is available */
    bool testContextFile(int stepNumber, int stepVersion);
    /**
     * Creates new DataReader for given domain.
     * Returns nonzero on success.
     * @param domainNum Domain number.
     * @param domainSerNum Domain serial number.
     * @param cmode Determines the i/o mode of context file.
     */
    DataReader *GiveDomainDataReader(int domainNum, int domainSerNum, ContextFileMode cmode);
    /**
     * Updates components mapped to numerical method if necessary during solution process.
     * Some numerical methods may require updating
     * mapped components during solution process (e.g., updating of tangent stiffness
     * when using updated Newton-Raphson method).
     * @param tStep Time when component is updated.
     * @param cmpn Numerical component to update.
     * @param d Domain.
     */
    virtual void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d);
    /**
     * Initializes solution of new time step. Default implementation
     * resets all internal history variables (in integration points of elements)
     * to previously reached equilibrium values.
     * Can be used for time step restart.
     */
    virtual void initStepIncrements();
    /**
     * Forces equation renumbering on given domain. All equation numbers in all dofManagers are invalidated,
     * and new equation numbers are generated starting from domainNeqs entry corresponding to given domain.
     * It will update numberOfEquations variable accordingly.
     * Should be used at startup to force equation numbering and therefore sets numberOfEquations.
     * Must be used if model supports changes of static system to assign new valid equation numbers
     * to dofManagers.
     */
    virtual int forceEquationNumbering(int i);
    /**
     * Forces equation renumbering on all domains associated to engng model.
     * All equation numbers in all domains for all dofManagers are invalidated,
     * and new equation numbers are generated starting from 1 on each domain.
     * It will update numberOfEquations variable accordingly.
     * Should be used at startup to force equation numbering and therefore sets numberOfEquations.
     * Must be used if model supports changes of static system to assign new valid equation numbers
     * to dofManagers.
     */
    virtual int forceEquationNumbering();
    /**
     * Indicates if EngngModel requires Dofs dictionaries to be updated.
     * If EngngModel does not support changes
     * of static system, the dof
     * forwards the requests for its unknowns to EngngModel, where unknowns are naturally kept.
     * This is possible, because dof equation number is same during whole solution.
     * But when changes of static system are allowed, several problem arise. For example
     * by solving simple incremental static with allowed static changes, the incremental displacement
     * vector of structure can not be added to total displacement vector of structure, because
     * equation numbers may have changed, and one can not simply add these vector to obtain new
     * total displacement vector, because incompatible displacement will be added.
     * To solve this problem, unknown dictionary at dof level has been assumed. Dof then keeps
     * its unknowns in its own private dictionary.
     * After computing increment of solution, engngModel updates for each dof its unknowns in its
     * dictionary (using updateUnknownsDictionary function). For aforementioned example
     * engngModel updates incremental values but also total value by asking dof for previous total
     * value (dof will use its dictionary, does not asks back EngngModel) adds corresponding increment
     * and updates total value in dictionary.
     */
    virtual int requiresUnknownsDictionaryUpdate() { return false; }
    /**
     * Returns true if equation renumbering is required for given solution step.
     * This may of course change the number of equation and in general there is no guarantee
     * that for a certain dof the same equation will be assigned. So the use of
     * DOF unknowns dictionaries is generally recommended.
     */
    virtual bool requiresEquationRenumbering(TimeStep *tStep) { return renumberFlag; }
    //virtual int supportsBoundaryConditionChange () {return 0;}
    /**
     * Updates necessary values in Dofs unknown dictionaries.
     * @see EngngModel::requiresUnknownsDictionaryUpdate
     * @see Dof::updateUnknownsDictionary
     */
    virtual void updateDofUnknownsDictionary(DofManager *, TimeStep *) { }
    /**
     * This method is responsible for computing unique dictionary id (ie hash value) from
     * given equationId, valueModeType and time step. This function is used by particular dofs
     * to access unknown identified by given parameters from its dictionary using computed index.
     * Usually the hash algorithm should produce index that depend on time step relatively to
     * actual one to avoid storage of complete history.
     */
    virtual int giveUnknownDictHashIndx(EquationID type, ValueModeType mode, TimeStep *stepN) { return 0; }
    /**
     * Returns UnknownNUmberingScheme related to given EquationID
     */
    virtual UnknownNumberingScheme& giveUnknownNumberingScheme(EquationID type) {
        return this->defaultNumberingScheme;
    }

    /**
     * Returns characteristic matrix of element. The Element::giveCharacteristicMatrix function
     * should not be called directly, because EngngModel may require some special modification
     * of characteristic matrices supported on element level. But default implementation does
     * the direct call to element level.
     * @param answer Characteristic matrix.
     * @param num Element number.
     * @param type Type of matrix requested.
     * @param tStep Time step when response is computed.
     * @param domain Source domain.
     */
    virtual void giveElementCharacteristicMatrix(FloatMatrix &answer, int num, CharType type, TimeStep *tStep, Domain *domain);
    /**
     * Returns characteristic vector of element. The Element::giveCharacteristicVector function
     * should not be called directly, because EngngModel may require some special modification
     * of characteristic vectors supported on element level. But default implementation does
     * the direct call to element level.
     * @param answer Characteristic vector.
     * @param num Element number.
     * @param type Type of vector requested.
     * @param mode Mode of unknown (total, incremental, rate of change).
     * @param tStep Time step when response is computed.
     * @param domain Source domain.
     */
    virtual void giveElementCharacteristicVector(FloatArray &answer, int num, CharType type, ValueModeType mode, TimeStep *tStep, Domain *domain);

#ifdef __PETSC_MODULE
    /**
     * Returns the PETSc context corresponding to given domain (n) and unknown type
     * Default implementation returns i-th context from petscContextList.
     */
    virtual PetscContext *givePetscContext(int n, EquationID eid);
    /**
     * Creates PETSc contexts. Must be implemented by derived classes since the governing equation type is required
     * for context creation.
     */
    virtual void initPetscContexts();
#endif

    /**
     * Assembles characteristic matrix of required type into given sparse matrix.
     * @param answer Assembled matrix.
     * @param tStep Time step, when answer is assembled.
     * @param eid Determines type of equation and corresponding element code numbers.
     * @param s Determines the equation numbering scheme.
     * @param type Characteristic components of type type are requested from elements and assembled.
     * @param domain Source domain.
     */
    virtual void assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
            CharType type, const UnknownNumberingScheme &s, Domain *domain);
    /**
     * Assembles characteristic matrix of required type into given sparse matrix.
     * @param answer assembled matrix
     * @param tStep Time step, when answer is assembled.
     * @param eid Determines type of equation and corresponding element code numbers.
     * @param r_s Determines the equation numbering scheme for the rows.
     * @param c_s Determines the equation numbering scheme for the columns.
     * @param type Characteristic components of type type are requested from elements and assembled.
     * @param domain Source domain.
     */
    virtual void assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
            CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, Domain *domain);
    /**
     * Assembles characteristic matrix of required type into given sparse matrix.
     * @param answer Assembled matrix.
     * @param tStep Time step, when answer is assembled.
     * @param r_id Determines type of equation and corresponding element code numbers for matrix rows.
     * @param c_id Determines type of equation and corresponding element code numbers for matrix columns.
     * @param s Determines the equation numbering scheme.
     * @param type Characteristic components of type type are requested from elements and assembled.
     * @param domain Source domain.
     */
    virtual void assemble(SparseMtrx *answer, TimeStep *tStep, EquationID r_id, EquationID c_id,
            CharType type, const UnknownNumberingScheme &s, Domain *domain);

    /**
     * Assembles characteristic vector of required type from dofManagers, element, and active boundary conditions, into given vector.
     * This routine is simple a convenient call to all three subroutines, since this is most likely what any engineering model will want to do.
     * The return value is used to normalize the residual when checking for convergence in nonlinear problems.
     * For parallel problems, the returned norm is also summed over all processes.
     * @param answer Assembled vector.
     * @param eid Determines type of equation and corresponding element code numbers.
     * @param mode Mode of unknown (total, incremental, rate of change).
     * @param tStep Time step, when answer is assembled.
     * @param type Characteristic components of type type are requested.
     * @param s Determines the equation numbering scheme.
     * @param domain Domain to assemble from.
     * @return Sum of element/node norm (squared) of assembled vector.
     */
    double assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                          CharType type, ValueModeType mode,
                          const UnknownNumberingScheme &s, Domain *domain);
    /**
     * Assembles characteristic vector of required type from dofManagers into given vector.
     * @param answer Assembled vector.
     * @param eid Determines type of equation and corresponding element code numbers.
     * @param mode Mode of unknown (total, incremental, rate of change).
     * @param tStep Time step, when answer is assembled.
     * @param type Characteristic components of type type are requested.
     * @param s Determines the equation numbering scheme.
     * @param domain Domain to assemble from.
     * @return Sum of element norm (squared) of assembled vector.
     */
    double assembleVectorFromDofManagers(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                         CharType type, ValueModeType mode,
                                         const UnknownNumberingScheme &s, Domain *domain);
    /**
     * Assembles characteristic vector of required type from elements into given vector.
     * @param answer Assembled vector.
     * @param tStep Time step, when answer is assembled.
     * @param eid Determines type of equation and corresponding element code numbers.
     * @param mode Mode of unknown (total, incremental, rate of change).
     * @param type Characteristic components of type type are requested
     * from elements and assembled using prescribed eqn numbers.
     * @param s Determines the equation numbering scheme.
     * @param domain Domain to assemble from.
     * @return Sum of element norm (squared) of assembled vector.
     */
    double assembleVectorFromElements(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                      CharType type, ValueModeType mode,
                                      const UnknownNumberingScheme &s, Domain *domain);

    /**
     * Assembles characteristic vector of required type from active boundary conditions.
     * @param answer Assembled vector.
     * @param tStep Time step, when answer is assembled.
     * @param eid Determines type of equation and corresponding element code numbers.
     * @param mode Mode of unknown (total, incremental, rate of change).
     * @param type Characteristic components of type type are requested
     * from elements and assembled using prescribed eqn numbers.
     * @param s Determines the equation numbering scheme.
     * @param domain Domain to assemble from.
     * @return Sum of element norm (squared) of assembled vector.
     */
    double assembleVectorFromActiveBC(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                      CharType type, ValueModeType mode,
                                      const UnknownNumberingScheme &s, Domain *domain);

    /**
     * Assembles the extrapolated internal forces vector,
     * useful for obtaining a good initial guess in nonlinear analysis with Dirichlet boundary conditions.
     * @param answer Assembled vector.
     * @param tStep Time step, when answer is assembled.
     * @param eid Determines type of equation and corresponding element code numbers.
     * @param type Determines the type of matrix to use, typically the tangent matrix or possibly the elastic tangent.
     * @param domain Domain to assemble from.
     * @return Sum of element norm (squared) of assembled vector.
     */
    void assembleExtrapolatedForces(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                    CharType type, Domain *domain);

protected:
#ifdef __PARALLEL_MODE
    /**
     * Packs receiver data when rebalancing load. When rebalancing happens, the local numbering will be lost on majority of processors.
     * Instead of identifying values of solution vectors that have to be send/received and then performing renumbering, all solution vectors
     * are assumed to be stored in dof dictionaries before data migration. Then dofs will take care themselves for packing and unpacking. After
     * data migration and local renumbering, the solution vectors will be restored from dof dictionary data back.
     */
    virtual void packMigratingData(TimeStep *) { }
    /**
     * Unpacks receiver data when rebalancing load. When rebalancing happens, the local numbering will be lost on majority of processors.
     * Instead of identifying values of solution vectors that have to be send/received and then performing renumbering, all solution vectors
     * are assumed to be stored in dof dictionaries before data migration. Then dofs will take care themselves for packing and unpacking. After
     * data migration and local renumbering, the solution vectors will be restored from dof dictionary data back.
     */
    virtual void unpackMigratingData(TimeStep *) { }
#endif

public:
    /**
     * Allows programmer to test some receiver's internal data, before computation begins.
     * @return Nonzero if receiver check is o.k.
     */
    virtual int checkConsistency() { return 1; }
    /**
     * Allows programmer to test problem its internal data, before computation begins.
     * @return Nonzero if receiver check is o.k.
     */
    int checkProblemConsistency();
    /**
     * Initializes the receiver state. Default implementation calls initModuleManager::doInit service to
     * invoke initialization by individual init modules.
     */
    virtual void init();

    /**
     * Prints output of receiver to output domain stream, for given time step.
     * Corresponding function for element gauss points is invoked
     * (gaussPoint::printOutputAt).
     */
    virtual void printOutputAt(FILE *file, TimeStep *tStep);


    // input / output
    /// Prints state of receiver. Useful for debugging.
    void printYourself();

    /**
     * DOF printing routine. Called by DofManagers to print Dof specific part.
     * Dof class provides component printing routines, but emodel is responsible
     * for what will be printed at DOF level.
     * @param stream output stream
     * @param iDof dof to be processed
     * @param atTime solution step
     */
    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime) = 0;


    // identification
    /// Returns class name of the receiver.
    virtual const char *giveClassName() const { return "EngngModel"; }
    /// Returns classType id of receiver.
    virtual classType giveClassID() const { return EngngModelClass; }
    /// Returns nonzero if receiver does incremental analysis.
    virtual int isIncremental() { return 0; }
    /// Returns nonzero if nonlocal stiffness option activated.
    virtual int useNonlocalStiffnessOption() { return 0; }
    /// Returns true if receiver in parallel mode
    bool isParallel() { return ( parallelFlag != 0 ); }

    /**
     * Indicates type of non linear computation (total or updated formulation).
     * This is used for example on Nodal level to update coordinates
     * if updated formulation
     * is done, or on element level, when non linear contributions are computed.
     */
    virtual fMode giveFormulation() { return nonLinFormulation; }
    /*
     * Returns Load Response Mode of receiver.
     * This value indicates, whether nodes and elements should assemble
     * total or incremental load vectors.
     *
     * virtual LoadResponseMode giveLoadResponseMode () {return TotalLoad;}
     */
    /// Context requesting service
    EngngModelContext *giveContext() { return this->context; }
    /// Returns number of slave problems.
    virtual int giveNumberOfSlaveProblems() { return 0; }
    /// Returns i-th slave problem.
    virtual EngngModel *giveSlaveProblem(int i) { return NULL; }

    /// Returns the Equation scaling flag, which is used to indicate that governing equation(s) are scaled, or non-dimensionalized.
    virtual bool giveEquationScalingFlag() { return false; }
    /// Returns the scale factor for given variable type.
    virtual double giveVariableScale(VarScaleType varId) { return 1.0; }


#ifdef __PARALLEL_MODE
    /**
     * Determines the space necessary for send/receive buffer.
     * It uses related communication map pattern to determine the maximum size needed.
     * @param commMap Communication map used to send/receive messages.
     * @param buff Communication buffer.
     * @param packUnpackType Determines the type of packed quantity, used by receiver
     * to estimate the size of pack/unpack buffer accordingly.
     * @return Upper bound of space needed.
     */
    virtual int estimateMaxPackSize(IntArray &commMap, CommunicationBuffer &buff, int packUnpackType) { return 0; }
    /**
     * Recovers the load balance between processors, if needed. Uses load balancer monitor and load balancer
     * instances to decide if rebalancing is needed (monitor) and to repartition the domain (load balancer).
     * Method is responsible for packing all relevant data (the use of dof dictionaries is assumed to store e-model
     * dof related staff, which can later help in renumbering after rebalancing) and to send/receive all data.
     * Then the local update and renumbering is necessary to get consistent data structure.
     */
    virtual void balanceLoad(TimeStep *tStep);
    /** Returns reference to receiver's load balancer. */
    virtual LoadBalancer *giveLoadBalancer() { return NULL; }
    /** Returns reference to receiver's load balancer monitor. */
    virtual LoadBalancerMonitor *giveLoadBalancerMonitor() { return NULL; }

    /// Returns domain rank in a group of collaborating processes (0..groupSize-1)
    int giveRank() { return rank; }
    /// Returns the number of collaborating processes.
    int giveNumberOfProcesses() { return numProcs; }
    /// Request domain rank and problem size
    void initParallel();
    /// Returns reference to itself -> required by communicator.h
    EngngModel *giveEngngModel() { return this; }
    /// Returns Communicator mode. Determines current domain-decomposition strategy used.
    ProblemCommunicatorMode giveProblemCommMode() { return this->commMode; }
#endif

#ifdef __OOFEG
    virtual void drawYourself(oofegGraphicContext &context);
    virtual void drawElements(oofegGraphicContext &context);
    virtual void drawNodes(oofegGraphicContext &context);
    /**
     * Shows the sparse structure of required matrix, type == 1 stiffness.
     */
    virtual void showSparseMtrxStructure(int type, oofegGraphicContext &context, TimeStep *atTime) { }
#endif

    /**@name Error and warning reporting methods.
     * These methods will print error (or warning) message using oofem default loggers.
     * Do not use these methods directly, to avoid specify file and line parameters.
     * More preferably, use these methods via corresponding OOFEM_CLASS_ERROR and OOFEM_CLASS_WARNING macros,
     * that will include file and line parameters automatically.
     *
     * Uses variable number of arguments, so a format string followed by optional arguments is expected
     * (according to printf conventions).
     * @param file Source file name, where error encountered (where error* function called).
     * @param line Source file line number, where error encountered.
     */
    //@{
    /// prints error message and exits
    void error(const char *file, int line, const char *format, ...) const;
    /// prints warning message
    void warning(const char *file, int line, const char *format, ...) const;
    //@}
};

} // end namespace oofem
#endif // engngm_h

