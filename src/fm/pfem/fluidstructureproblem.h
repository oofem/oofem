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

///@name Input fields for FluidStructureProblem
//@{
#define _IFT_FluidStructureProblem_Name "fluidstuctureproblem"
#define _IFT_FluidStructureProblem_rtolv "rtolv"
#define _IFT_FluidStructureProblem_rtolp "rtolp"
#define _IFT_FluidStructureProblem_maxiter "maxiter"

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
    /// Iteration counter
    int iterationNumber;

    /// Convergence tolerance.
    double rtolv, rtolp;
    /// Max number of iterations.
    int maxiter;

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

    void solveYourselfAt(TimeStep *tStep) override;
    void initializeYourself(TimeStep *tStep) override;
    int initializeAdaptive(int stepNumber) override { return 0; }

    void initializeFrom(InputRecord &ir) override;

    void printYourself();
    void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime) override { }

    void preInitializeNextStep() override;

    // identification
    const char *giveClassName() const override { return "FluidStructureProblem"; }
    const char *giveInputRecordName() const override { return _IFT_FluidStructureProblem_Name; }
    int isIncremental() { return 0; }
    int useNonlocalStiffnessOption() override { return 0; }

    fMode giveFormulation() override { return UNKNOWN; }

    /// Returns list of model number that this model is coupled with. Used for staggered approach.
    void giveCoupledModels(IntArray &answer) { answer = coupledModels; }

#ifdef __OOFEG
    /**
     * Shows the sparse structure of required matrix, type == 1 stiffness.
     */
    void showSparseMtrxStructure(int type, oofegGraphicContext &context, TimeStep *atTime) override { }
#endif

    int giveNumberOfSlaveProblems() override { return ( int ) inputStreamNames.size(); }

    int giveIterationNumber() { return iterationNumber; }
};
} // end namespace oofem
#endif // fluidstructureproblem_h
