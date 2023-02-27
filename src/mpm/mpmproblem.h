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

#ifndef mpmproblem_h
#define mpmproblem_h

#include "engngm.h"
#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "sparsemtrx.h"
#include "primaryfield.h"
#include "function.h"

///@name Input fields for mpmproblem
//@{
#define _IFT_MPMProblem_Name "mpmproblem"
#define _IFT_MPMProblem_nsmax "nsmax"
#define _IFT_MPMProblem_rtol "rtol"
#define _IFT_MPMProblem_manrmsteps "manrmsteps"
#define _IFT_MPMProblem_initt "initt"
#define _IFT_MPMProblem_deltat "deltat"
#define _IFT_MPMProblem_deltatfunction "deltatfunction"
#define _IFT_MPMProblem_prescribedtimes "prescribedtimes"
#define _IFT_MPMProblem_alpha "alpha"
#define _IFT_MPMProblem_changingproblemsize "changingproblemsize"


//@}

namespace oofem {

/**
 * Callback class for assembling mid point effective tangents. 
 * @todo Need to parametrize individual contributing terms, ther locations and multilication factors.
 */
class MPMLhsAssembler : public MatrixAssembler
{
protected:
    double alpha;
    double deltaT;

public:
    MPMLhsAssembler(double alpha, double deltaT);
    void matrixFromElement(FloatMatrix &mat, Element &element, TimeStep *tStep) const override;
};

/**
 * Callback class for assembling rhs external forces
 */
class MPMRhsAssembler : public VectorAssembler
{
protected:
    double alpha;
    double deltaT;
public:
    MPMRhsAssembler(double alpha, double deltaT) : VectorAssembler(), alpha(alpha), deltaT(deltaT) {}
    void vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const override;
};

/**
 * Callback class for assembling residuals
 */
class MPMResidualAssembler : public VectorAssembler
{
    protected:
    double alpha;
    double deltaT;
public:
    MPMResidualAssembler(double alpha, double deltaT) : VectorAssembler(), alpha(alpha), deltaT(deltaT) {}
    void vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const override;
};


/**
 * This class represents generic nonlinear multi-physics problem. The problem can be growing/decreasing, signalized by flag "changingproblemsize"
 * in the problem description. The solution is stored in UnknownsField, which can obtain/ project solution from/to DOFs (nodes). If the problem
 * keeps the same equation numbers, solution is taken from UnknownsField without any projection, which is more efficient. See the matlibmanual
 * for solution strategy of balance equations and the solution algorithm.
 *
 * @todo Documentation errors (there is no "UnknownsField" used here).
 */
class MPMProblem : public EngngModel
{
protected:
    enum nlttp_ModeType { nrsolverModifiedNRM, nrsolverFullNRM, nrsolverAccelNRM };

    SparseMtrxType sparseMtrxType = SMT_Skyline;
    /// This field stores solution vector. For fixed size of problem, the PrimaryField is used, for growing/decreasing size, DofDistributedPrimaryField applies.
    std :: unique_ptr< PrimaryField > UnknownsField;

    std :: unique_ptr< SparseMtrx > jacobianMatrix;
    FloatArray internalForces;
    FloatArray eNorm;

    /// Numerical method used to solve the problem
    std::unique_ptr<SparseNonLinearSystemNM> nMethod;
    LinSystSolverType solverType = ST_Direct; ///@todo Remove this and use nonlinear methods.
    std :: unique_ptr< SparseLinearSystemNM > linSolver; ///@todo Remove this and use nonlinear methods.

    bool keepTangent = false;

    double rtol = 0.;
    int nsmax = 0;
    nlttp_ModeType NR_Mode = nrsolverModifiedNRM;
    int MANRMSteps = 0;
    int currentIterations = 0;

    /// Initial time from which the computation runs. Default is zero.
    double initT = 0.;
    /// Length of time step.
    double deltaT = 0.;
    double alpha = 0.;
    /// Associated time function for time step increment.
    int dtFunction = 0;
    /// Specified times where the problem is solved
    FloatArray discreteTimes;
    /// Determines if there are change in the problem size (no application/removal of Dirichlet boundary conditions).
    bool changingProblemSize = false;
    /**
     * Contains last time stamp of internal variable update.
     * This update is made via various services
     * (like those for computing real internal forces or updating the internal state).
     */
    StateCounterType internalVarUpdateStamp;

public:
    MPMProblem(int i, EngngModel * _master);

    TimeStep* giveNextStep() override;
    void solveYourselfAt(TimeStep *tStep) override;
    void updateYourself(TimeStep *tStep) override;
    double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof) override;

    void initializeFrom(InputRecord &ir) override;

    // identification
    const char *giveInputRecordName() const { return _IFT_MPMProblem_Name; }
    const char *giveClassName() const override { return "MPMProblem"; }
    fMode giveFormulation() override { return nonLinFormulation; }
    int giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep) override;
    void updateDofUnknownsDictionary(DofManager *dman, TimeStep *tStep) override;

    int giveCurrentNumberOfIterations() override { return currentIterations; }
    TimeStep *giveSolutionStepWhenIcApply(bool force = false) override;
    NumericalMethod *giveNumericalMethod(MetaStep *mStep) override;


protected:
    void updateInternalState(TimeStep *tStep) ;
    void applyIC(TimeStep *tStep) ;
    void createPreviousSolutionInDofUnknownsDictionary(TimeStep *tStep);
    /**
     * Returns the time step length for given step number n, initial step is number 0.
     */
    double giveDeltaT(int n);
    /**
     * Copy unknowns in DOF's from previous to current position.
     * @param mode What the unknown describes (increment, total value etc.).
     * @param fromTime From which time step to obtain value.
     * @param toTime To which time to copy.
     */
    virtual void copyUnknownsInDictionary(ValueModeType mode, TimeStep *fromTime, TimeStep *toTime);
    Function * giveDtFunction();
    /**
     * Returns time for time step number n (array discreteTimes must be specified)
     */
    double giveDiscreteTime(int n);
};
} // end namespace oofem
#endif // mpmproblem_h
