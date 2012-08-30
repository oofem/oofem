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

#ifndef staggeredproblem_h
#define staggeredproblem_h

#include "engngm.h"
#include "alist.h"

#include "inputrecord.h"

namespace oofem {
/**
 * Implementation of general sequence (staggered) problem. The problem consists in sequence of
 * low level problems (slaves) which are executed sequentially and where the results
 * of particular slave depends on the results of previous slaves in sequence.
 * Typical example is heat&mass transfer analysis followed by mechanical one, which
 * takes into account the temperature field from the first analysis.
 *
 * The sequence problem is represented by this class. It maintains list
 * of subsequent (slave) problems and it is executes the slave problems. It is responsible
 * for solution step generation and synchronization between slave problems.
 * The transfer of required state variables is done by mapping of corresponding variables
 * between problem domains. This allows to to transfer primary (nodal) values of one problem to
 * integration points of subsequent problem or to use completely different discretizations for
 * slave problems.
 *
 * Since the master problem is responsible for synchronization, it is responsible for
 * generation the solution steps. Therefore, the solution step specification, as well as
 * relevant meta step attributes are specified at master level.
 *
 * @note To avoid confusion,
 * the slaves are treated in so-called maintained mode. In this mode, the attributes and
 * meta step attributes are taken from the master. The local attributes, even if specified,
 * are ignored.
 *
 * @todo Move to oofemlib
 */
class StaggeredProblem : public EngngModel
{
protected:
    /// Number of engineering models to run.
    int nModels;
    /// List of engineering models to solve sequentially.
    AList< EngngModel > *emodelList;
    double deltaT;
    std :: string *inputStreamNames;
    /// Associated time function for time step increment
    int dtTimeFunction;
    /**
     * Constant multiplier, optional input parameter. This parameter determines the ratio of
     * two consecutive time steps. Efficient for creep and relaxation analyses.
     */
    double stepMultiplier;

    /// Optional parameter which specify problems to define load time functions
    int timeDefinedByProb;

public:
    /**
     * Constructor. Creates an engineering model with number i belonging to domain d.
     */
    StaggeredProblem(int i, EngngModel *_master = NULL);
    /// Destructor.
    virtual ~StaggeredProblem();

    void setContextOutputMode(ContextOutputMode contextMode);
    void setUDContextOutputMode(int cStep);
    void setProblemMode(problemMode pmode);
    virtual void setRenumberFlag();

    virtual void solveYourself();
    virtual void solveYourselfAt(TimeStep *tStep);
    virtual int forceEquationNumbering();
    virtual void updateYourself(TimeStep *stepN);
    virtual void initializeYourself(TimeStep *tStep) { }
    virtual int initializeAdaptive(int stepNumber) { return 0; }
    virtual void terminate(TimeStep *tStep);
    virtual void doStepOutput(TimeStep *tStep);

    virtual int instanciateYourself(DataReader *dr, InputRecord *ir, const char *outFileName, const char *desc);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void updateAttributes(MetaStep *mStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual void updateDomainLinks();

    void printYourself();
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime) { }
    virtual TimeStep *giveNextStep();
    virtual TimeStep *giveSolutionStepWhenIcApply();

    // identification
    virtual const char *giveClassName() const { return "StaggeredProblem"; }
    virtual classType giveClassID() const { return StaggeredProblemClass; }
    virtual int isIncremental() { return 0; }
    virtual int useNonlocalStiffnessOption() { return 0; }

    virtual fMode giveFormulation() { return UNKNOWN; }
    /**
     * Returns time function for time step increment.
     * Used time function should provide step lengths as function of step number.
     * Initial step with number 0 is considered as [ -dt(0), 0 ], first step is [ 0, dt(1) ], ...
     */
    LoadTimeFunction *giveDtTimeFunction();

    /**
     * Returns the timestep length for given step number n, initial step is number 0
     */
    double giveDeltaT(int n);

#ifdef __OOFEG
    void drawYourself(oofegGraphicContext &context);
    void drawElements(oofegGraphicContext &context);
    void drawNodes(oofegGraphicContext &context);
    /**
     * Shows the sparse structure of required matrix, type == 1 stiffness.
     */
    virtual void showSparseMtrxStructure(int type, oofegGraphicContext &context, TimeStep *atTime) { }
#endif

    virtual int checkConsistency();

    virtual EngngModel *giveSlaveProblem(int i);
    virtual int giveNumberOfSlaveProblems() { return nModels; }

    virtual int giveNumberOfFirstStep() { if ( master ) { return master->giveNumberOfFirstStep(); } else { return 1; } }
    virtual int giveNumberOfTimeStepWhenIcApply() {
        if ( master ) { return master->giveNumberOfTimeStepWhenIcApply(); } else { return 0; }
    }
    virtual int instanciateDefaultMetaStep(InputRecord *ir);

protected:
    int instanciateSlaveProblems();
};
} // end namespace oofem
#endif // staggeredproblem_h
