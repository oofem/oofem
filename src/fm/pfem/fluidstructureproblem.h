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

#ifndef fluidstructureproblem_h
#define fluidstructureproblem_h

#include "staggeredproblem.h"
#include "inputrecord.h"

///@name Input fields for StaggeredProblem
//@{
#define _IFT_FluidStructureProblem_Name "fluidstuctureproblem"
#define _IFT_StaggeredProblem_deltat "deltat"
#define _IFT_StaggeredProblem_dtf "dtf"
#define _IFT_StaggeredProblem_timeDefinedByProb "timedefinedbyprob"
#define _IFT_StaggeredProblem_stepmultiplier "stepmultiplier"
#define _IFT_StaggeredProblem_prescribedtimes "prescribedtimes"
#define _IFT_StaggeredProblem_prob1 "prob1"
#define _IFT_StaggeredProblem_prob2 "prob2"
#define _IFT_StaggeredProblem_coupling "coupling"
//@}

namespace oofem {
/**
 * Implementation of fluid-structure interaction ) problem based on Dirichlet-Neumann
 * approach. The problem consists in pair of structural (DIIDynamic) and fluid (PFEM)
 * problem.
 *
 * The FluidStructureProblem provides an iterative synchronization of sub-problems.
 * The convergence criterion is based on the difference of the pressure and velocity
 * values on the interface from the subsequent iterative steps.
 */
class OOFEM_EXPORT FluidStructureProblem : public StaggeredProblem
{
protected:
    IntArray interactionParticles;
    double tol;
    int iterationNumber;

public:
    /**
     * Constructor. Creates an engineering model with number i belonging to domain d.
     */
    FluidStructureProblem(int i, EngngModel *_master = NULL);
    /// Destructor.
    virtual ~FluidStructureProblem();

    void setContextOutputMode(ContextOutputMode contextMode);
    void setUDContextOutputMode(int cStep);
    void setProblemMode(problemMode pmode);

    virtual void solveYourselfAt(TimeStep *tStep);
    virtual void initializeYourself(TimeStep *tStep);
    virtual int initializeAdaptive(int stepNumber) { return 0; }

    void printYourself();
    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime) { }

    virtual void preInitializeNextStep();

    // identification
    virtual const char *giveClassName() const { return "FluidStructureProblem"; }
    virtual const char *giveInputRecordName() const { return _IFT_FluidStructureProblem_Name; }
    virtual int isIncremental() { return 0; }
    virtual int useNonlocalStiffnessOption() { return 0; }

    virtual fMode giveFormulation() { return UNKNOWN; }


    /// Returns list of model number that this model is coupled with. Used for staggered approach.
    void giveCoupledModels(IntArray &answer) { answer = coupledModels; }

#ifdef __OOFEG
    /**
     * Shows the sparse structure of required matrix, type == 1 stiffness.
     */
    virtual void showSparseMtrxStructure(int type, oofegGraphicContext &context, TimeStep *atTime) { }
#endif

    virtual int giveNumberOfSlaveProblems() { return ( int ) inputStreamNames.size(); }

    virtual int giveNumberOfFirstStep() { if ( master ) { return master->giveNumberOfFirstStep(); } else { return 1; } }
    virtual int giveNumberOfTimeStepWhenIcApply() {
        if ( master ) { return master->giveNumberOfTimeStepWhenIcApply(); } else { return 0; }
    }

    int giveIterationNumber() { return iterationNumber; }
};
} // end namespace oofem
#endif // fluidstructureproblem_h
