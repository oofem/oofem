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

#ifndef nlinearstatic_h
#define nlinearstatic_h

#include "sm/EngineeringModels/linearstatic.h"
#include "sparsenonlinsystemnm.h"
#include "function.h"

///@name Input fields for NonLinearStatic
//@{
#define _IFT_NonLinearStatic_Name "nonlinearstatic"
#define _IFT_NonLinearStatic_controlmode "controlmode"
#define _IFT_NonLinearStatic_deltat "deltat"
#define _IFT_NonLinearStatic_deltatfunction "deltatfunction"
#define _IFT_NonLinearStatic_stiffmode "stiffmode"
#define _IFT_NonLinearStatic_refloadmode "refloadmode"
#define _IFT_NonLinearStatic_keepll "keepll"
#define _IFT_NonLinearStatic_donotfixload "donotfixload"
#define _IFT_NonLinearStatic_nonlocstiff "nonlocstiff"
#define _IFT_NonLinearStatic_nonlocalext "nonlocalext"
#define _IFT_NonLinearStatic_loadBalancingFlag "loadbalancingflag"
#define _IFT_NonLinearStatic_forceloadBalancingFlag "forceloadbalancingflag"
#define _IFT_NonLinearStatic_updateElasticStiffnessFlag "updateelasticstiffnessflag"
//@}

namespace oofem {
/// Type determining the stiffness mode.
enum NonLinearStatic_stiffnessMode {
    nls_tangentStiffness = 0, ///< The tangent stiffness is used and updated whenever requested.
    nls_secantStiffness = 1, ///< The secant stiffness is used and updated whenever requested.
    nls_elasticStiffness = 2, ///< The initial elastic stiffness is used in the whole solution.
    nls_secantInitialStiffness = 3, ///< The secant stiffness is used and updated only at the beginning of new load step.
};
/// Type determining type of loading control. This type determines the solver to be used.
enum NonLinearStatic_controlType {
    nls_indirectControl = 0, ///< A generalized norm of displacement and loading vectors is controlled. In current implementation, the CALM solver is used, the reference load vector is FIXED.
    nls_directControl = 1,   ///< Describes the direct control where load or displacement (or both) are controlled.
};

/**
 * This class implements nonlinear static engineering problem.
 * Solution of this problem is performed  as a series of increments (loading or displacement).
 * At start of Each increment we assemble new tangent stiffness, and iteratively trying
 * to fulfill balance of external and real internal forces
 * at end of load step (see numerical method ).
 * The loading applied can bo of two types:
 * - proportional incremental  loading
 * - non-proportional fixed loading, reflecting the previous history,
 *   but could not be scaled (like dead weight).
 *
 * Tasks:
 * - Creating Numerical method for solving nonlinear problem.
 * - Assembling tangent stiffness matrix.
 * - Interfacing Numerical method to Elements.
 * - Managing time steps.
 */
class NonLinearStatic : public LinearStatic
{
protected:
    double prevStepLength, currentStepLength;
    FloatArray totalDisplacement,  incrementOfDisplacement;
    FloatArray internalForces;

    /// A load vector already applied, which does not scales.
    FloatArray initialLoadVector;
    /**
     * Incremental load vector which is applied in loading step (as a whole for direct control
     * or proportionally for indirect control).
     */
    FloatArray incrementalLoadVector;
    /// A load vector which does not scale for prescribed DOFs.
    FloatArray initialLoadVectorOfPrescribed;
    /// Incremental Load Vector for prescribed DOFs.
    FloatArray incrementalLoadVectorOfPrescribed;

    double loadLevel, cumulatedLoadLevel;
    bool mstepCumulateLoadLevelFlag;
    int currentIterations;
    NonLinearStatic_stiffnessMode stiffMode;
    int loadInitFlag;
    int nonlocalStiffnessFlag;
    bool updateElasticStiffnessFlag;
    NM_Status numMetStatus;
    /// Numerical method used to solve the problem.
    SparseNonLinearSystemNM *nMethod;
    /// Characterizes the type of control used.
    NonLinearStatic_controlType controlMode;
    /// Intrinsic time increment.
    double deltaT;
    /// Associated time function for time step increment.
    int dtFunction;

   /**
     * The following parameter allows to specify how the reference load vector
     * is obtained from given totalLoadVector and initialLoadVector.
     * The initialLoadVector desribes the part of loading which does not scale.
     * If refLoadInputMode is rlm_total (default) then the reference incremental load vector is defined as
     * totalLoadVector assembled at given time.
     * If refLoadInputMode is rlm_inceremental then the reference load vector is
     * obtained as incremental load vector at given time.
     */
    SparseNonLinearSystemNM :: referenceLoadInputModeType refLoadInputMode;

    /// The initial guess type to use before starting the nonlinear solver.
    InitialGuess initialGuessType;

public:
    NonLinearStatic(int i, EngngModel *master = nullptr);
    virtual ~NonLinearStatic();

    Function *giveDtFunction();
    double giveDeltaT(int n);
    void solveYourself() override;
    void solveYourselfAt(TimeStep *tStep) override;
    void terminate(TimeStep *tStep) override;

    void printOutputAt(FILE *file, TimeStep *tStep) override;

    void updateComponent(TimeStep *tStep, NumericalCmpn, Domain *d) override;
    void updateSolution(FloatArray &solutionVector, TimeStep *tStep, Domain *d) override;
    void updateMatrix(SparseMtrx &mat, TimeStep *tStep, Domain *d) override;
    void updateAttributes(MetaStep *mStep) override;

    double giveUnknownComponent(ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof) override;
    void initializeFrom(InputRecord &ir) override;
    TimeStep *giveSolutionStepWhenIcApply(bool force = false) override;
    TimeStep *giveNextStep() override;
    NumericalMethod *giveNumericalMethod(MetaStep *mStep) override;

    double giveLoadLevel() override { return cumulatedLoadLevel + loadLevel; }

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    void updateDomainLinks() override;

    // identification
    const char *giveClassName() const override { return "NonLinearStatic"; }
    fMode giveFormulation() override { return nonLinFormulation; }
    int useNonlocalStiffnessOption() override { return this->nonlocalStiffnessFlag; }
    /// For load balancing purposes we store all values so hash is computed from mode value only
    int giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep) override { return ( int ) mode; }
    int giveCurrentNumberOfIterations() override { return currentIterations; }

#ifdef __OOFEG
    void showSparseMtrxStructure(int type, oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

    int estimateMaxPackSize(IntArray &commMap, DataStream &buff, int packUnpackType) override;
#ifdef __PARALLEL_MODE
    LoadBalancer *giveLoadBalancer() override;
    LoadBalancerMonitor *giveLoadBalancerMonitor() override;

#endif

protected:
    void assemble(SparseMtrx &answer, TimeStep *tStep, const MatrixAssembler &ma,
                  const UnknownNumberingScheme &, Domain *domain) override;
    void proceedStep(int di, TimeStep *tStep);
    virtual void updateLoadVectors(TimeStep *tStep);
    void computeExternalLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di) override;
    void assembleIncrementalReferenceLoadVectors(FloatArray &_incrementalLoadVector,
                                                 FloatArray &_incrementalLoadVectorOfPrescribed,
                                                 SparseNonLinearSystemNM :: referenceLoadInputModeType _refMode,
                                                 Domain *sourceDomain, TimeStep *tStep);

    void packMigratingData(TimeStep *tStep) override;
    void unpackMigratingData(TimeStep *tStep) override;
};
} // end namespace oofem
#endif // nlinearstatic_h
