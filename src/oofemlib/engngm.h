/* $Header: /home/cvs/bp/oofem/oofemlib/src/engngm.h,v 1.32.4.1 2004/04/05 15:19:43 bp Exp $ */
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


//   ************************
//   *** CLASS ENGNGMODEL ***
//   ************************


#ifndef engngm_h
#define engngm_h

#ifdef __PARALLEL_MODE
#include "parallel.h"
#include "loadbalancer.h"
#endif

#include "alist.h"

#include "inputrecord.h"
#include "datareader.h"

#include "sparsemtrx.h"
#include "flotarry.h"
#include "element.h"
#include "dofmanager.h"
#include "exportmodulemanager.h"
#include "field.h"
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

#ifndef __MAKEDEPEND
#include <stdio.h>

#ifdef __PETSC_MODULE
#include "petsccontext.h"
#include "petscordering.h"
#endif
#endif


#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#endif

class Domain;
class NumericalMethod;
class TimeStep;
class ErrorEstimator;
class MetaStep;
class MaterialInterface;

/**
 * Class EngngModelContext represents a context, which is shared by all problem engng sub-models.
 * In principle every problem (represented by the EngngModel class) can be made part of more
 * complex problem, providing only part of the solution, possibly depending on the results of
 * other sub-problems. Typical example is staggered heat and structural analysis.
 * The context provides the common (shared) resources within the problem. The context is created by the
 * master problem (the class represendting staggered problem, for example). The subproblems are
 * then created in so-called maintained (or slave) mode and they request the context from master.
 */
class EngngModelContext
{
protected:
    /// Common fieldManager providing shared field regisry for the problem
    FieldManager fieldManager;
public:

    EngngModelContext() { }
    FieldManager *giveFieldManager() { return & ( this->fieldManager ); }
};





/**
 * Abstract base class representing the "problem" under consideration.
 * The enginerring model describes the problem and type of analysis to be done.
 * It "knows" the type and form of governing equation, knows how to assemble this
 * problem from local element contributions and uses suitable instance of
 * numerical method class to solve problem. It posses and manages one or more problem domains.
 * Concept of time steps is introduced. For problems discretized in time the introduction of
 * time step is natural. In other cases, time steps can represent for example
 * load increments or different load cases.
 *
 * The solution steps are grouped together into so called meta steps. The meta step can be thought as an sequence of
 * solution steps, with same set of attributes used to drive behaviour of engng model.
 * For each metastep, the enegng model typically updates its controll attributes according to
 * metastep engng attributes (see initMetaStepAttributes and updateAttributes services)
 * and creates the solution steps acoordingly. This allows to switch to
 * different time increment, different solution controll, etc. If no metastep is specified, the engng model
 * creates default one for all required solution steps. There are two services, where attributes are updated, the first one,
 * used for those attributes, which do not vary during solution of problem are set in initializeForm service.
 * The second service is updateAttributes, where the attributes allowed to change (with metastep validity) are updated.
 * If no metastep is introduced, default one is created (with attributes set to engng model init record).
 * Then there is no difference, whether attributes are read in initializeFrom or updateAttributes, but
 * preffered scheme is to read all attributes in initializeFrom and left updateAttributes service empty.
 *
 * The basic EngngModel tasks are following
 * <UL>
 * <LI>
 * assembling governing equations by summing contributions from problem domains (typically from nodes and elements),</LI>
 * <LI>
 * Solving the problem descibed by governing equation(s) using suitable instance of
 * numerical method. This requires interfacing numericalMethod characteristic elements
 * with components in governing equation
 * EngngModel must map each component of governing
 * equation(s) (which has physical meaning) to corresponding numerical component of Numerical
 * method. This mapping between physical components to independent numerical components
 * (uderstand by numerical method) is important, because it allows numerical method to be used by
 * many EngngModel with diferent meaning of particular components.</LI>
 * <LI>
 * Returning unknown values (acording to requsted type and mode). Used by Dofs to
 * access their corresponding unknowns.</LI>
 * <LI>
 * Terminating time step by updating nodal and element values (including integration points update).</LI>
 * <LI>
 * Updating dofs unknowns dictionaries if  model supports changes of static system (see Dof class
 * documentation for detailed explanation). In general if static system changes are not supported,
 * when dof are requested for unknowns, they use their associate equation number to ask EngngModel
 * for this unknown. Unknowns are therefore stored in EngngModel and are requested by dofs.
 * On the other hand, when static system changes are supported, the equation numbers of dofs
 * can vary during solution. Therefore, so called unknowns dictionary at dof level are introduced.
 * All unknowns are stored on dof level and dofs will use in such case their own dictionaries
 * instead of requsting EngngModel. The EngngModel is fully responsible to update this
 * dictionary for each dof with all necessary unknowns (see updateDofUnknownsDictionary function).</LI>
 * </UL>
 *
 */
class EngngModel
{
    /*
     * This class implements the abstract class of problem to be solved by using FEM.
     * DESCRIPTION :
     * The Engineering Model ( like Static, EigenValue , Dynamic, Thermal, ... problems)
     * handles the problem to be solved and  governing equation of this problem
     * ( Kr=f, Mr"+Cr'+Kr=f, ... ).
     * The model contains key-component numericalMethod, that will be used to
     * solve governing equations. Model is responsible to recognize allowed numerical
     * method, to give type of characteristics matrices that it need.
     * The model is also responsible for selecting best numerical model according to
     * number of unknowns, prescribed accuracy, type of load, previously computed
     * results .... .
     *
     * if updating of DOFs uknowns dictionary is necessary the funcion
     * requiresNodeUnknowsDictionaryUpdate() returns nonzero. This may be used
     * for models, which supports dyynamic changes of static system, where renumbering
     * of equatins is necessary. Then is necessary to overload
     * updateDofUnknownsDictionary() function to update DOF unknowns dictionary, where
     * unknowns are hold instead of keeping them in global unknowns vectors in engng instancies.
     * Then after finishing loading step, unknowns are updated in DOFs dictionaries
     * using updateDofUnknownsDictionary() function called from this->terminate().
     * Then global vectors can be deleted, new equation numbers are set (because we change
     * static system) and new solution for next step can begin.
     * Keeping unknowns in DOFs dictionary results to independence on equation numbering
     * changes during overall solution process.
     *
     * By default requiresNodeUnknowsDictionaryUpdate() returns 0 and unknowns are
     * stored in global vectors owned by engng model, and dofs are asking for them
     * using this->giveUnknownComponent() member function.
     *
     * This function is not necesarry - we use unknowns DOF dictionary instead
     * so particular dofs are not asking emodel, but use their own dictionary
     * - this dictionary is set maintained by emodel->updateDofUnknownsDictionary()
     *
     *
     * TASKS :
     * - Interfacing numericalMethod to elements. for example if Numerical Method solve
     *   following type of problem Au"+Bu=f, this class is responsible to say:
     *   yes for my type of problem A is mass matrix and B is stiffness matrix.
     * - returning end of time of interest (usefull for some time depeemndent models).
     * - assembling governing equations by summing contributions from nodes and elements
     * - if updaating of DOFs uknowns dictionary is necessary the funcion
     * requiresNodeUnknowsDictionaryUpdate() returns nonzero, and is necessary to overload
     * updateDofUnknownsDictionary() function to update this map;
     * - initializing state variables stored in elements gp's acording to
     * initial conditions using function initializeYourself()
     */

public:

#ifdef __PARALLEL_MODE
    enum EngngModel_UpdateMode { EngngModel_Unknown_Mode, EngngModel_SUMM_Mode, EngngModel_SET_Mode };
#endif

protected:
    /// number of receiver domains
    int ndomains;
    /// List of problem domains
    AList< Domain > *domainList;
    /// Total number of time steps
    int numberOfSteps;
    /// total number of equation in cuurent time step
    int numberOfEquations;
    /// total number or prescribed equations in current time step
    int numberOfPrescribedEquations;
    /// number of equations per domain
    IntArray domainNeqs;
    /// number of prescribed equations per domain
    IntArray domainPrescribedNeqs;
    /// renumbering flag
    int renumberFlag;
    /// equation numbering completed flag
    int equationNumberingCompleted;
    /// number of meta steps
    int nMetaSteps;
    /// List of problem metasteps
    AList< MetaStep > *metaStepList;
    /// Solution step when IC (initial conditions) apply
    TimeStep *stepWhenIcApply;
    /// Currnet time step
    TimeStep *currentStep;
    /// Previous time step
    TimeStep *previousStep;
    /// receivers id
    int number;

    /// Path to input stream
    //char*                 dataInputFileName ;
    /// Path to output stream
    char *dataOutputFileName;
    /// Output stream
    FILE *outputStream;
    /// Input stream
    //FILE* inputStream;
    /// Domain context output mode
    ContextOutputMode contextOutputMode;
    int contextOutputStep;

    ///Export module manager
    ExportModuleManager *exportModuleManager;


    /// Domain mode
    problemMode pMode;
    /// solution start time
    time_t startTime;
    // initial value of processor time used by program
    // clock_t startClock;

    /// master e-model; if defined receiver is in maintained (slave) mode
    EngngModel *master;
    /// context
    EngngModelContext *context;
    /// e-model timer
    EngngModelTimer timer;
    /**
     * flag indicating that the receiver runs in parallel.
     */
    int parallelFlag;


#ifdef __PARALLEL_MODE
    /// domain rank in a group of colaborating processes (0..groupSize-1)
    int rank;
    /// total number of colaborating processes
    int numProcs;
    /// processor name
    char processor_name [ PROCESSOR_NAME_LENGTH ];

    /**@name Load balancing attributes */
    //@{
    /// Load Balancer
    LoadBalancer *lb;
    LoadBalancerMonitor *lbm;
    /// if set to true, load balancing is active
    bool loadBalancingFlag;
    /// debug flag forcing load balancing after first step
    bool force_load_rebalance_in_first_step;
    //@}

#endif // __PARALLEL_MODE

#ifdef __PETSC_MODULE
    /// list where petsc contexts are stored
    AList< PetscContext > *petscContextList;
#endif


public:
    /**
     * Constructor. Creates Engng model with number i belonging to domain d.
     */
    EngngModel(int i, EngngModel *_master = NULL);   // constructor
    /**
     * Constructor. Creates Engng model with number i and input file given by path.
     */
    EngngModel(int i, char *s, EngngModel *_master = NULL);
    /// Destructor.
    virtual ~EngngModel();      // destructor

    /**
     * Service for accessing particular problem domain.
     * Generates error if no such domain is defined.
     * @param: n pointer to n-th domain is returned
     */
    Domain *giveDomain(int n);
    /// Returns number of domains in problem.
    int         giveNumberOfDomains() { return ndomains; }
    /** Service for accessing ErrorEstimator corresponding to particular domain */
    virtual ErrorEstimator *giveDomainErrorEstimator(int n) { return NULL; }
    /** Returns material interface representation for given domain */
    virtual MaterialInterface *giveMaterialInterface(int n) { return NULL; }

    // input / output
    /// Returns input file path.
    //char*              giveInputDataFileName () ;
    /// Returns file descriptor of output file
    FILE *giveOutputStream();
    /** Returns base output file name
     *  to which extensions, like .out .vtk .osf should be added.
     *  In current implementation, output file name is simply returned.
     *  @param path and base file name will be copied into the array pointed to by  dest
     *  @param not more than n bytes of src  are copied
     */
    char *giveOutputBaseFileName(char *dest, size_t n) { return strncpy(dest, dataOutputFileName, n); }

    //FILE*              giveInputStream () ;

    /*
     * Returns current time in seconds as returned by time call.
     * @return current time in time_t structure.
     */
    //time_t             getTime ();
    /*
     * Returns an approximation of processor time used by the program.
     * The value returned is the  CPU  time  used  so  far  as  a
     * clock_t;  to  get  the  number  of seconds used, divide by
     * CLOCKS_PER_SEC. Calls clock ANSI C function.
     * The C standard allows for arbitrary values at the start of
     * the   program;  take  the  difference  between  the  value
     * returned from a call to this method at the start of  the  pro-
     * gram and the end to get maximum portability.
     */
    //clock_t            getClock ();
    /**
     * Returns domain context output mode.
     */
    ContextOutputMode  giveContextOutputMode() { return contextOutputMode; }
    /**
     * Returns domain context output step.
     */
    int                giveContextOutputStep() { return contextOutputStep; }
    /**
     * Sets context output mode of receiver.
     * @param contextMode domain context mode.
     */
    void               setContextOutputMode(ContextOutputMode contextMode)
    { contextOutputMode = contextMode; }
    /**
     * Sets user defined context output mode (it sets contextOutputMode to contextOutputMode),
     * setting contextOutputStep to given value.
     * @param cStep new context output step
     */
    void               setUDContextOutputMode(int cStep)
    { contextOutputMode = USERDEFINED;
      contextOutputStep = cStep; }
    /**
     * Sets domain mode to given mode.
     * @param mode domain mode.
     */
    void               setProblemMode(problemMode mode) { pMode = mode; }
    /// Returns domain mode.
    problemMode         giveProblemMode()        { return pMode; }
    /// Sets the renumber flag to TRUE
    virtual void                setRenumberFlag() { this->renumberFlag = 1; }
    /// Sets the renumber flag to FALSE
    virtual void                resetRenumberFlag() { this->renumberFlag = 0; }

    /**
     * Performs analysis termination after finishing analysis.
     */
    void               terminateAnalysis();

    // solving
    /**
     * Starts solution process. Implementation should invoke for each time step
     * solveYourselfAt function with time step as parameter. Time steps are created
     * using giveNextStep function (this will set current time step to newly created,
     * and updates previous step).
     */
    virtual void               solveYourself();
    /**
     * Solves problem for given time step. Should assemble characteristic matrices and vectors
     * if necessary and solve problem using appropriate numerical method. After finishing solution,
     * this->updateYourself function for updating solution state and then this->terminate
     * function (for updating nodal and element values) should be called.
     */
    virtual void               solveYourselfAt(TimeStep *) { }
    //virtual int                requiresNewLhs () {return 1;}
    /**
     * Terminates the solution of time step. Default implementation calls prinOutput() service and if specified,
     * context of whole domain is stored and output for given time step is printed.
     */
    virtual void               terminate(TimeStep *);
    /**
     * Prints the ouput of the solution step (using virtual this->printOutputAtservice)
     * to the stream detemined using this->giveOutputStream() method
     * and calls exportModuleManager to do output.
     */
    virtual void              doStepOutput(TimeStep *);
    /**
     * Saves context of given solution step, if required (determined using this->giveContextOutputMode() method).
     */
    void                       saveStepContext(TimeStep *);
    /**
     * Updates internal state after finishing time step. (for example total values may be
     * updated according to previously solved increments).  Then element values are also updated
     * (together with related integration points and material statuses).
     */
    virtual void               updateYourself(TimeStep *stepN);
    /**
     * Provides the oportunity to initialize state variables stored in element
     * integration points acording to
     * initial conditions using function initializeYourself() on element level.
     * Should be called when curent time step is time step when IC will aply
     * (see EngngModel::giveNumberOfTimeStepWhenIcApply)
     * somewhere from solveYourselfAt function). Implementation must be provided.
     * Default implementation is empty.
     */
    virtual void               initializeYourself(TimeStep *) { }
    /**
     * Initializes the newly generated discretization state acording to previous solution.
     * This process should typically include restoring old solution, instanciating newly
     * generated domain(s) and by mapping procedure.
     */
    virtual int                initializeAdaptive(int stepNumber) { return 0; }

    /**
     * Returns total number of equations in active (current time steep) time step.
     * The UnknownType parameter allows to distinguis between several possible governing equations, that
     * can be numbered separately.
     */
    virtual int giveNumberOfEquations(EquationID);
    /**
     * Returns total number of prescribed equations in active (current time steep) time step.
     * The UnknownType parameter allows to distinguis between several possible governing equations, that
     * can be numbered separately.
     */
    virtual int giveNumberOfPrescribedEquations(EquationID);
    /**
     * Returns number of equations for given domain in active (current time steep) time step.
     * The UnknownType parameter allows to distinguis between several possible governing equations, that
     * can be numbered separately.
     */
    virtual int giveNumberOfDomainEquations(int, EquationID);
    /**
     * Returns number of prescribed equations for given domain in active (current time steep) time step.
     * The UnknownType parameter allows to distinguis between several possible governing equations, that
     * can be numbered separately.
     */
    virtual int giveNumberOfPrescribedDomainEquations(int, EquationID);
    //virtual IntArray*          GiveBanWidthVector ();


    // management  components
    /**
     * Provides backward mapping between numerical component and characteristic
     * component on EngngModel level.
     */
    virtual CharType  giveTypeOfComponent(NumericalCmpn) { return UnknownCharType; }
    /**
     * Returns requested unknown. Unknown at give time step is characterized by its type and mode
     * and by its equation number. This function is used by Dofs, when they are requsted for
     * their associated unknowns.
     * @see Dof::giveUnknown method
     */
    virtual double    giveUnknownComponent(EquationID, ValueModeType, TimeStep *, Domain *, Dof *) { return 0.0; }
    virtual double    giveUnknownComponent(UnknownType, ValueModeType, TimeStep *, Domain *, Dof *) { return 0.0; }

#ifdef __PARALLEL_MODE
    /**
     * Updates unknown. Unknown at give time step is characterized by its type and mode
     * and by its equation number. This function is used by Dofs, when they are requsted for
     * their update of associated unknowns. Declared only in __PARALLEL_MODE
     * @see Dof::giveUnknown method
     */
    virtual void updateUnknownComponent(EquationID, ValueModeType, TimeStep *, int,
                                        double, EngngModel_UpdateMode) { return; }

#endif
    /**
     * Initializes whole problem acording to its description stored in inputStream.
     * Prints header, opens the outFileName, instanciate itself the receicer using
     * using virtual initializeFrom service and instancites all problem domains.
     */
    virtual int instanciateYourself(DataReader *dr, InputRecord *ir, char *outFileName, char *desc);
    /**
     * Initializes receiver acording to object description in input reader.
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
     * be read. The service takes TimeStep as parameter, from which corresponding MetaStep is
     * requested. It is recomended, to implement this service in such way, that multiple calls
     * for steps belonging to same MetaStep does not change response.
     * The default implementation updates the numerical method attributes.
     * @param TimeStep time step.
     */
    virtual void updateAttributes(TimeStep *);
    /**
     * Update e-model attributes attributes according to step metaStep attributes.
     * Calls updateAttributes. At the end the meta step input reader finish() service
     * is called in order to allow for unread attribute check.
     */
    void initMetaStepAttributes(TimeStep *tStep);
    /**
     * Stores the  state of model to output stream. Stores not only the receiver state,
     * but also same function is invoked for all DofManagers and Elements in associated
     * domain. Note that by storing element  context also contexts of all associated
     * integration points (and material statuses) are stored.
     * Stored context is associated with current time step. One time step can have only
     * one associated context. Multiple call to saveContext within same time step
     * owerride previously saved context for this step.
     * By default the stream paprameter is used to store data and is not closed.
     * If stream is NULL, new file descriptor is created and this must be also closed at the end.
     * @param stream - context stream. If NULL then new file descriptor will be openned and closed
     * at the end else the stream given as parameter will be used and not closed at the end.
     * @param mode determines ammount of info in stream
     * @return contextIOResultType.
     * @exception throws an ContextIOERR exception if error encountered
     */
    virtual contextIOResultType                saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores the  state of model from output stream. Restores not only the receiver state,
     * but also same function is invoked for all DofManagers and Elements in associated
     * domain. Note that by restoring element  context also contexts of all associated
     * integration points (and material statuses) are restored.
     * Each context is associated with unique time step. Only one context per time step is
     * allowed. Restore context function will restore such contex, which is related
     * (through its step number) to time step number and version given in obj parameter.
     * Restoring context will change current time step in order to correspond to newly restored
     * context.
     * @param stream context file
     * @param mode determines ammount of info in stream
     * @param obj is a void pointer to an int array containing two values:time step number and
     * version of a context file to be restored.
     * @return contextIOResultType.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    virtual contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Updates domain links after the domains of receiver have changed. Used mainly after
     * restoring context - the domains may change and this service is then used
     * to update domain variables in all components belonging to receiver
     * like errorestimators, solvers, etc, having domains as attributes.
     */
    virtual void updateDomainLinks() { };
    void               resolveCorrespondingStepNumber(int &, int &, void *obj);
    /// Returns current time step.
    TimeStep *giveCurrentStep() { if ( master ) { return master->giveCurrentStep(); } else { return currentStep; } }
    /// Returns previous time step.
    TimeStep *givePreviousStep() { if ( master ) { return master->givePreviousStep(); } else { return previousStep; } }
    /// Returns next time step (next to current step) of receiver.
    virtual TimeStep *giveNextStep() { return NULL; }
    /// Returns the solution step when Initial Conditions (IC) apply
    virtual TimeStep *giveSolutionStepWhenIcApply() { if ( master ) { return master->giveCurrentStep(); } else { return stepWhenIcApply; } }
    /// Returns number of first time step used by receiver.
    virtual int        giveNumberOfFirstStep() { if ( master ) { return master->giveNumberOfFirstStep(); } else { return 1; } }
    /// Return number of meta steps
    int                 giveNumberOfMetaSteps() { return nMetaSteps; }
    /// Returns the i-th meta step
    MetaStep *giveMetaStep(int i);
    /// Returns total number of steps.
    int                giveNumberOfSteps() { if ( master ) { return master->giveNumberOfSteps(); } else { return numberOfSteps; } }
    /// Returns end of time interest (time corresponding to end of time integration).
    virtual double     giveEndOfTimeOfInterest() { return 0.; }
    /// Returns the time step number, when initial conditions should apply.
    virtual int       giveNumberOfTimeStepWhenIcApply() {
        if ( master ) { return master->giveNumberOfTimeStepWhenIcApply(); } else { return 0; } }
    /// Returns reference to receiver's numerical method
    virtual NumericalMethod *giveNumericalMethod(TimeStep *) { return NULL; }
    /// Returns receiver's export mudule manager
    ExportModuleManager *giveExportModuleManager() { return exportModuleManager; }
    /// Returns reference to recever timer (EngngModelTimer)
    EngngModelTimer *giveTimer() { return & timer; }

    /**
     * Increases number of equations of receiver's domain and returns newly created equation number.
     * Used mainly by DofManagers to allocate their corresponding equation number if it
     * is not currently allocated.
     * The DofIDItem parameter allows to distinguis between several possible governing equations, that
     * can be numbered separately.
     */
    virtual int      giveNewEquationNumber(int domain, DofIDItem) { return ++domainNeqs.at(domain); }
    /**
     * Increases number of prescribed equations of receiver's domain and returns newly created equation number.
     * Used mainly by DofManagers to allocate their corresponding equation number if it
     * is not currently allocated.
     * The DofIDItem parameter allows to distinguis between several possible governing equations, that
     * can be numbered separately.
     */
    virtual int      giveNewPrescribedEquationNumber(int domain, DofIDItem) { return ++domainPrescribedNeqs.at(domain); }
    /**
     * Assigns context file-descriptor for given step number to stream.
     * Returns nonzero on success.
     * @param stepNumber solution step number to store/restore
     * @param stepVersion version of step
     * @param cmode determines the i/o mode of context file
     * @param errLevel determines the amout of warning messages if errors are encountered, level 0 no warnings reported.
     */
    int              giveContextFile(FILE **contextFile, int stepNumber, int stepVersion,
                                     ContextFileMode cmode, int errLevel = 1);
    /** Returns true if context file for given step and version is available */
    bool             testContextFile(int stepNumber, int stepVersion);
    /**
     * Creates new DataReader for given domain.
     * Returns nonzero on success.
     * @param domainNum domain number
     * @param domainSerNum domain seerial number
     * @param cmode determines the i/o mode of context file
     */
    DataReader *GiveDomainDataReader(int domainNum, int domainSerNum, ContextFileMode cmode);
    /**
     * Updates components mapped to numerical method if necessary during solution process.
     * Some numerical methods may require updating
     * mapped components during solution process (e.g., updating of tanget stiffness
     * when using updated Newton-Raphson method).
     * @param tStep time when component is updated.
     * @param cmpn Numerical component to update.
     * @param d domain
     */
    virtual void      updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d);
    /**
     * Initializes solution of new time step. Default implementation
     * resets all internal history variables (in integration points of elements)
     * to previously reached equilibrium values.
     * Can be used for time step restart.
     */
    virtual void      initStepIncrements();
    /**
     * Forces equation renumbering on given domain. All equation numbers in all dofManagers are invalidated,
     * and new equation numbers are generated starting from domainNeqs entry corresponding to given domain.
     * It will update numberOfEquations variable accordingly.
     * Should be used at startup to force equation numbering and therefore sets numberOfEquations.
     * Must be used if model supports changes of static system to assign  new valid equation numbers
     * to dofManagers.
     */
    virtual int       forceEquationNumbering(int i);
    /**
     * Forces equation renumbering on all domains associated to engng model.
     * All equation numbers in all domains for all dofManagers are invalidated,
     * and new equation numbers are generated starting from 1 on each domain.
     * It will update numberOfEquations variable accordingly.
     * Should be used at startup to force equation numbering and therefore sets numberOfEquations.
     * Must be used if model supports changes of static system to assign  new valid equation numbers
     * to dofManagers.
     */
    virtual int       forceEquationNumbering();
    /**
     * Indicates if Engngmodel requires Dofs dictionaries to be updated.
     * If EngngModel does not support changes
     * of static system, the dof
     * frowards the requests for its unknowns to EngngModel, where unknowns are naturaly kept.
     * This is posible, because dof equation number is same during whole solution.
     * But when changes of static system are allowed, several problem arise. For example
     * by solving simple  incremental static with allowed static changes, the incremetal displacement
     * vector of structure can not be added to total displacement vector of structure, because
     * equation numbers may have changed, and one can not simply add these vector to obtain new
     * total displacement vector, because uncompatible displacement will be added.
     * To solve this problem, uknown dictionary at dof level has been assumed. Dof then keeps
     * its unknowns in its onw private dictionary.
     * After computing increment of solution, engngModel updates for each dof its unknowns  in its
     * dictionary (using updateUnknownsDictionary function). For aforementioned example
     * engngModel updates incremental values but also total value by asking dof for previous total
     * value (dof will use its dictionary, does not asks back EngngModel) adds corresponding increment
     * and updates total value in dictionary.
     */
    virtual int       requiresUnknownsDictionaryUpdate() { return 0; }
    /**
     * Returns true if equation renumbering is required for given solution step.
     * This may of course change the number of equation and in general there is no gauarantee
     * that for a certain dof the same eautiaon will be assigned. So the use of
     * DOF unknowns dictionaries is generally recomended.
     */
    virtual bool requiresEquationRenumbering(TimeStep *) { return false; }
    //virtual int       supportsBoundaryConditionChange () {return 0;}
    /**
     * Updates necessary values in Dofs unknown dictionaries.
     * @see EngngModel::requiresUnknownsDictionaryUpdate
     * @see Dof::updateUnknownsDictionary
     */
    virtual void      updateDofUnknownsDictionary(DofManager *, TimeStep *) { }
    /**
     *  This method is responsible for computing unique dictionary id (ie hash value) from
     *  given equationId, valueModeType and timestep. This function is used by particular dofs
     *  to access unknown identified by given params from its dictionary using computed index.
     *  Usually the hash algorithm shoud produce index that depend on timestep relativelly to
     *  actual one to avoid storage of complete history.
     */
    virtual int       giveUnknownDictHashIndx(EquationID type, ValueModeType mode, TimeStep *stepN) { return 0; }

    // we don't directlt call element ->GiveCharacteristicMatrix() function, because some
    // engngm classes may require special modification of base types supported on
    // element class level
    /**
     * Returns characteristic matrix of element. The Element::GiveCharacteristicMatrix function
     * should not be called directly, because EngngModel may require some special modification
     * of characteristic matrices supported on element level. But default implementation does
     * the direct call to element level.
     * @param answer characteristic matrix
     * @param num element number
     * @param type type of CharMatrix requsted
     * @param tStep time step when response is computed
     * @param domain source domain
     */
    virtual void giveElementCharacteristicMatrix(FloatMatrix &answer, int num, CharType type, TimeStep *tStep, Domain *domain)
    { domain->giveElement(num)->giveCharacteristicMatrix(answer, type, tStep); }
    /**
     * Returns characteristic vector of element. The Element::GiveCharacteristicVector function
     * should not be called directly, because EngngModel may require some special modification
     * of characteristic vectors supported on element level. But default implementation does
     * the direct call to element level.
     * @param answer characteristic vector
     * @param num element number
     * @param type type of vector requsted
     * @param tStep time step when response is computed
     * @param domain source domain
     */
    virtual void giveElementCharacteristicVector(FloatArray &answer, int num, CharType type, ValueModeType mode, TimeStep *tStep, Domain *domain)
    { domain->giveElement(num)->giveCharacteristicVector(answer, type, mode, tStep); }

#ifdef __PETSC_MODULE
    /**
     * Returns the petsc context corresponding to given domain (n) and unknown type
     * Default implementation returns i-th context from petscContextList.
     */
    virtual PetscContext *givePetscContext(int n, EquationID ut);
    /**
     * Creates Petsc contexts. Must be implemented by derived classes since the governing equation type is reqired
     * for context creation.
     */
    virtual void initPetscContexts();
#endif


protected:
    /**
     * Assembles characteristic matrix of required type into given sparse matrix.
     * @param answer assembled matrix
     * @param tStep time step, when answer is assembled.
     * @param ut determines type of equation and corresponding element code numbers
     * @param type characterisctic components of type type are requsted from elements and assembled.
     * @param domain source domain
     */
    virtual void       assemble(SparseMtrx *answer, TimeStep *tStep, EquationID ut, CharType type, Domain *domain);
    /**
     * Assembles characteristic matrix of required type into given sparse matrix.
     * @param answer assembled matrix
     * @param tStep time step, when answer is assembled.
     * @param r_id determines type of equation and corresponding element code numbers for matrix rows
     * @param c_id determines type of equation and corresponding element code numbers for matrix columns
     * @param type characterisctic components of type type are requsted from elements and assembled.
     * @param domain source domain
     */
    virtual void       assemble(SparseMtrx *answer, TimeStep *tStep, EquationID r_id, EquationID c_id, CharType type, Domain *domain);
    /**
     * Assembles characteristic vector of required type into given vector.
     * @param answer assembled vector
     * @param tStep time step, when answer is assembled.
     * @param type characterisctic components of type type are requsted
     * from dofManagers/elements and assembled.
     */
    //virtual void       assemble (FloatArray&, TimeStep*, CharType type, Domain* domain) ;
    /**
     * Assembles characteristic vector of required type from dofManagers into given vector.
     * @param answer assembled vector
     * @param tStep time step, when answer is assembled.
     * @param type characterisctic components of type type are requsted
     * from dofManagers and assembled using code numbers.
     */
    virtual void assembleVectorFromDofManagers(FloatArray &, TimeStep *, EquationID ut, CharType type, ValueModeType mode, Domain *domain);
    /**
     * Assembles prescribed characteristic vector of required type from dofManagers into given vector.
     * @param answer assembled vector
     * @param tStep time step, when answer is assembled.
     * @param type characterisctic components of type type are requsted
     * from dofManagers and assembled using prescribed eqn numbers.
     */
    void assemblePrescribedVectorFromDofManagers(FloatArray &, TimeStep *, EquationID, CharType type, ValueModeType mode, Domain * domain);
    /**
     * Assembles characteristic vector of required type from elements into given vector.
     * @param answer assembled vector
     * @param tStep time step, when answer is assembled.
     * @param type characterisctic components of type type are requsted
     * from elements and assembled using  using code numbers.
     */
    void assembleVectorFromElements(FloatArray &, TimeStep *, EquationID, CharType type, ValueModeType mode, Domain * domain);
    /**
     * Assembles prescribed characteristic vector of required type from elements into given vector.
     * @param answer assembled vector
     * @param tStep time step, when answer is assembled.
     * @param type characterisctic components of type type are requsted
     * from elements and assembled using prescribed eqn numbers.
     */
    void assemblePrescribedVectorFromElements(FloatArray &, TimeStep *, EquationID, CharType type, ValueModeType mode, Domain * domain);

#ifdef __PETSC_MODULE
    void petsc_assembleVectorFromDofManagers(Vec, TimeStep *, EquationID ut, CharType type, ValueModeType mode, Domain * domain);
    void petsc_assemblePrescribedVectorFromDofManagers(Vec, TimeStep *, EquationID ut, CharType type, ValueModeType mode, Domain * domain);
    void petsc_assembleVectorFromElements(Vec, TimeStep *, EquationID ut, CharType type, ValueModeType mode, Domain * domain);
    void petsc_assemblePrescribedVectorFromElements(Vec, TimeStep *, EquationID ut, CharType type, ValueModeType mode, Domain * domain);
#endif
#ifdef __PARALLEL_MODE
    /** Packs receiver data when rebalancing load. When rebalancing happens, the local numbering will be lost on majority of processors.
     *  Instead of identifying values of solution vectors that have to be send/received and then performing renumbering, all solution vectors
     *  are assumed to be stored in dof dictionaries before data migration. Then dofs will take care themselves for packing and unpacking. After
     *  data migration and local renubering, the solution vectors will be restored from dof dictionary data back.
     */
    virtual void packMigratingData(TimeStep *) { }
    /** Unpacks receiver data when rebalancing load. When rebalancing happens, the local numbering will be lost on majority of processors.
     *  Instead of identifying values of solution vectors that have to be send/received and then performing renumbering, all solution vectors
     *  are assumed to be stored in dof dictionaries before data migration. Then dofs will take care themselves for packing and unpacking. After
     *  data migration and local renubering, the solution vectors will be restored from dof dictionary data back.
     */
    virtual void unpackMigratingData(TimeStep *) { }
#endif
public:

    // consistency check
    /** Allows programmer to test some receiver's internal data, before computation begins.
     * @return nonzero if receiver check is o.k. */
    virtual int checkConsistency() { return 1; }    // returns nonzero if o.k.
    /** Allows programmer to test problem its internal data, before computation begins.
     * @return nonzero if receiver check is o.k. */
    int checkProblemConsistency();      // returns nonzero if o.k.

    /**
     * Prints output of receiver to ouput domain stream, for given time step.
     * Corresponding function for element gauss points is invoked
     * (gaussPoint::printOutputAt).
     */
    virtual void                  printOutputAt(FILE *, TimeStep *);


    // input / output
    /// Prints stete of receiver. Usefull for debugging.
    void printYourself();

    /** DOF printing routine. Called by DofManagers to print Dof specific part.
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
    /// retun true if receiver in parallel mode
    bool isParallel() { return ( parallelFlag != 0 ); }

    /**
     * Indicates type of non linear computation (total or updated formulation).
     * This is used for example on Nodal level to update coordinates
     * if updated formulation
     * is done, or on element level, when non linear contributions are computed.
     */
    virtual fMode giveFormulation() { return UNKNOWN; } // for non-linear computation
    /*
     * Returns Load Response Mode of receiver.
     * This value indicates, whether nodes and elements should assemble
     * total or incremental load vectors.
     *
     * virtual  LoadResponseMode giveLoadResponseMode () {return TotalLoad;}
     */
    /// Context requesting service
    EngngModelContext *giveContext() { return this->context; }
    /**
     * Returns number of slave problems */
    virtual int giveNumberOfSlaveProblems() { return 0; }
    /**Returns i-th slave problem */
    virtual EngngModel *giveSlaveProblem(int i) { return NULL; }

    /// Returns the Equation scaling flag, which is used to indicate that governing equation(s) are scaled, or non-dimensionalized
    virtual bool giveEquationScalingFlag() { return false; }
    /// Returns the scale factor for given variable type
    virtual double giveVariableScale(VarScaleType varId) { return 1.0; }


#ifdef __PARALLEL_MODE
    /**
     * Determines the space necessary for send/receive buffer.
     * It uses related communication map pattern to determine the maximum size needed.
     * @param commMap communication map used to send/receive messages
     * @param buff communication buffer
     * @param packUnpackType determines the type of packed quantity, used by receiver
     * to estimate the size of pack/unpack buffer accordingly.
     * @return upper bound of space needed
     */
    virtual int estimateMaxPackSize(IntArray &commMap, CommunicationBuffer &buff, int packUnpackType) { return 0; }
    /**
     * Recovers the load balance between processors, if needed. Uses load balancer monitor and load balancer
     * instances to decide if rebalancing is needed (monitor) and to repartition the domain (load balancer).
     * Method is responsible for packing all relevant data (the use of dof dictionaries is assumed to store e-model
     * dof related staff, which can later help in renumbering after rebalancing) and to send/receive all data.
     * Then the local update and renumbering is necessary to get consistent data structure.
     */
    virtual void balanceLoad(TimeStep *);
    /** returns reference to receiver's load balancer*/
    virtual LoadBalancer *giveLoadBalancer() { return NULL; }
    /** returns reference to receiver's load balancer monitor*/
    virtual LoadBalancerMonitor *giveLoadBalancerMonitor() { return NULL; }


#endif

#ifdef __PARALLEL_MODE

    /// Returns domain rank in a group of colaborating processes (0..groupSize-1)
    int                giveRank() { return rank; }
    /// Returns the number of colaborating processis
    int                giveNumberOfProcesses() { return numProcs; }
    /// Request domain rank and problem size
    void               initParallel();
    /// Returns reference to itself -> required by comunicator.h
    EngngModel *giveEngngModel() { return this; }
#endif

#ifdef __OOFEG
    virtual void               drawYourself(oofegGraphicContext &context);
    virtual void               drawElements(oofegGraphicContext &context);
    virtual void               drawNodes(oofegGraphicContext &context);
    /**
     * Shows the sparse structure of required matrix, type == 1 stiffness.
     */
    virtual void       showSparseMtrxStructure(int type, oofegGraphicContext &context, TimeStep *atTime) { }
#endif

    /**@name error and warning reporting methods
     * These methods will print error (or warning) message using oofem default loggers.
     * Do not use these methods directly, to avoid specify file and line parameters.
     * More preferably, use these methods via corresponding OOFEM_CLASS_ERROR and OOFEM_CLASS_WARNING macros,
     * that will include file and line parameters automatically.
     *
     * Uses variable number of arguments, so a format string followed by optional argumens is expected
     * (according to printf conventions).
     * @param file  source file name, where error encountered (where error* function called)
     * @param line  source file line number, where error encountered
     */
    //@{
    /// prints error message and exits
    void error(const char *file, int line, const char *format, ...) const;
    /// prints warning message
    void warning(const char *file, int line, const char *format, ...) const;
    //@}
};

typedef EngngModel Problem;

#endif // engngm_h









