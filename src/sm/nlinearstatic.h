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

#ifndef nlinearstatic_h
#define nlinearstatic_h

#include "linearstatic.h"
#include "sparsenonlinsystemnm.h"

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
    nls_directControl2 = 2,  ///<
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
    NM_Status numMetStatus;
    /// Numerical method used to solve the problem.
    SparseNonLinearSystemNM *nMethod;
    /// Characterizes the type of control used.
    NonLinearStatic_controlType controlMode;
    /// Intrinsic time increment.
    double deltaT;
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
    NonLinearStatic(int i, EngngModel *_master = NULL);
    virtual ~NonLinearStatic();

    virtual void solveYourself();
    virtual void solveYourselfAt(TimeStep *tStep);
    virtual void terminate(TimeStep *tStep);

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void updateYourself(TimeStep *tStep);
    virtual void updateComponent(TimeStep *tStep, NumericalCmpn, Domain *d);
    virtual void updateAttributes(MetaStep *mStep);

    virtual double giveUnknownComponent(EquationID eid, ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);
    virtual double giveUnknownComponent(UnknownType ut, ValueModeType type, TimeStep *tStep, Domain *d, Dof *dof);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual TimeStep *giveNextStep();
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual void updateDomainLinks();

    // identification
    virtual const char *giveClassName() const { return "NonLinearStatic"; }
    virtual classType giveClassID() const { return NonLinearStaticClass; }
    virtual int isIncremental() { return 1; }
    virtual fMode giveFormulation() { return nonLinFormulation; }
    virtual int useNonlocalStiffnessOption() { return this->nonlocalStiffnessFlag; }
    /// For load balancing purposes we store all values with same EquationID; so hash is computed from mode value only
    virtual int giveUnknownDictHashIndx(EquationID type, ValueModeType mode, TimeStep *stepN)
    { return ( int ) mode; }

#ifdef __OOFEG
    void showSparseMtrxStructure(int type, oofegGraphicContext &context, TimeStep *atTime);
#endif

#ifdef __PARALLEL_MODE
    virtual int estimateMaxPackSize(IntArray &commMap, CommunicationBuffer &buff, int packUnpackType);

    virtual LoadBalancer *giveLoadBalancer();
    virtual LoadBalancerMonitor *giveLoadBalancerMonitor();

#endif

#ifdef __PETSC_MODULE
    virtual void initPetscContexts();
#endif

protected:
    virtual void assemble(SparseMtrx *answer, TimeStep *tStep, EquationID ut, CharType type,
                  const UnknownNumberingScheme &, Domain *domain);
    void proceedStep(int di, TimeStep *tStep);
    void updateLoadVectors(TimeStep *tStep);
    virtual void computeExternalLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di);
    void assembleIncrementalReferenceLoadVectors(FloatArray &_incrementalLoadVector,
                                                 FloatArray &_incrementalLoadVectorOfPrescribed,
                                                 SparseNonLinearSystemNM :: referenceLoadInputModeType _refMode,
                                                 Domain *sourceDomain, EquationID ut, TimeStep *tStep);
#ifdef __PARALLEL_MODE
    virtual void packMigratingData(TimeStep *tStep);
    virtual void unpackMigratingData(TimeStep *tStep);
#endif
};
} // end namespace oofem
#endif // nlinearstatic_h
