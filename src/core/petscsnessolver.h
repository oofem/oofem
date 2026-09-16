/*
 * PETSc SNES adapter for OOFEM nonlinear problems.
 */

#ifndef petscsnessolver_h
#define petscsnessolver_h

#include "sparsenonlinsystemnm.h"

#include <petscsnes.h>

#include <map>
#include <memory>
#include <string>

#define _IFT_PetscSNESSolver_Name "petscsnes"
#define _IFT_PetscSNESSolver_maxIterations "snesmaxiter"
#define _IFT_PetscSNESSolver_maxFunctionEvaluations "snesmaxfunc"
#define _IFT_PetscSNESSolver_absoluteTolerance "snesatol"
#define _IFT_PetscSNESSolver_relativeTolerance "snesrtol"
#define _IFT_PetscSNESSolver_stepTolerance "snesstol"
#define _IFT_PetscSNESSolver_type "snestype"
#define _IFT_PetscSNESSolver_optionsPrefix "snesoptionsprefix"
#define _IFT_PetscSNESSolver_monitor "snesmonitor"
#define _IFT_PetscSNESSolver_solutionDependentExternalForces "soldepextforces"
#define _IFT_PetscSNESSolver_fieldScaling "snesfieldscaling"
#define _IFT_PetscSNESSolver_forceScale "forcescale"
#define _IFT_PetscSNESSolver_forceScaleDofs "forcescaledofs"
#define _IFT_PetscSNESSolver_boundedTrustRegion "snestrustregion"
#define _IFT_PetscSNESSolver_trustInitialRadius "snestrdelta0"
#define _IFT_PetscSNESSolver_trustMinimumRadius "snestrdeltamin"
#define _IFT_PetscSNESSolver_trustMaximumRadius "snestrdeltamax"
#define _IFT_PetscSNESSolver_trustMaximumTrials "snestrmaxtrials"
#define _IFT_PetscSNESSolver_minimumFieldScaleSquared "snesminscalesquared"
#define _IFT_PetscSNESSolver_boundDofIds "snesbounddofs"
#define _IFT_PetscSNESSolver_lowerBounds "sneslowerbounds"
#define _IFT_PetscSNESSolver_upperBounds "snesupperbounds"
#define _IFT_PetscSNESSolver_irreversibleDofIds "snesirreversibledofs"

namespace oofem {

class PetscSparseMtrx;

/**
 * Delegates nonlinear globalization to PETSc SNES while retaining OOFEM's
 * residual and tangent assembly.  The solved PETSc residual is
 *
 *     f(X) = F_internal(X) - R_external(X).
 *
 * Its Jacobian therefore has the same sign as the tangent assembled by the
 * existing OOFEM engineering models.
 *
 * The initial implementation supports load control for serial OOFEM analyses
 * and requires a PetscSparseMtrx. Runtime PETSc options (for example
 * -snes_type newtontr) override the defaults read from the OOFEM input record.
 */
class OOFEM_EXPORT PetscSNESSolver : public SparseNonLinearSystemNM
{
protected:
    int maxIterations;
    int maxFunctionEvaluations;
    double absoluteTolerance;
    double relativeTolerance;
    double stepTolerance;
    std :: string snesType;
    std :: string optionsPrefix;
    bool monitor;
    bool solutionDependentExternalForces;
    bool fieldScaling;
    /// OOFEM linear solver used only by an engineering-model predictor.
    std :: unique_ptr< SparseLinearSystemNM > predictorLinearSolver;
    std :: map< int, double > fieldForceScale;
    double minimumFieldScaleSquared = 1.e-6;
    // Optional trust-region truncation of the reduced-space Newton direction.
    // PETSc retains VI active-set selection and projects every trial to bounds.
    bool boundedTrustRegion = false;
    double trustInitialRadius = 1.0, trustMinimumRadius = 1.e-8;
    double trustMaximumRadius = 1.e8, trustRadius = 1.0;
    int trustMaximumTrials = 20;
    FloatArray trustStepWeights;
    static PetscErrorCode applyBoundedTrustRegion(SNESLineSearch search, void *context);
    PetscErrorCode boundedTrustRegionStep(SNESLineSearch search);
    PetscReal boundTolerance = 1.e-8;
    IntArray boundDofIds, irreversibleDofIds;
    FloatArray lowerBounds, upperBounds;
    // Freeze the start-of-step values across cutbacks and staggered sweeps.
    // Captured before the first solve: do not predict irreversible fields.
    int boundReferenceStep = -1;
    FloatArray boundReferenceSolution;
    FloatArray equationLowerBounds, equationUpperBounds;
    /// Numbering of the vectors/matrix passed to the current compact solve.
    const UnknownNumberingScheme *customEquationNumbering = nullptr;

    PetscSparseMtrx *activeMatrix;
    FloatArray *activeReferenceLoad;
    FloatArray *activeInitialLoad;
    FloatArray *activeSolution;
    FloatArray *activeInternalForces;
    const FloatArray *activeInternalForcesEBENorm;
    TimeStep *activeTimeStep;
    FloatArray residualRowScale;
    bool residualRowScaleInitialized;
    int functionEvaluations;
    int jacobianEvaluations;
    int domainErrors;
    std :: string callbackFailure;
    std :: string lastDomainError;

    SNESConvergedReason lastReason;
    int lastIterations;
    double lastFunctionNorm;

    static PetscErrorCode formFunction(SNES snes, Vec x, Vec f, void *context);
    static PetscErrorCode formJacobian(SNES snes, Vec x, Mat jacobian,
                                       Mat preconditioner, void *context);
    static PetscErrorCode monitorIteration(SNES snes, PetscInt iteration,
                                           PetscReal functionNorm, void *context);

    PetscErrorCode evaluateFunction(SNES snes, Vec x, Vec f);
    PetscErrorCode evaluateJacobian(SNES snes, Vec x, Mat jacobian);
    void applyTrialSolution(SNES snes, Vec x, bool updateInternalForces);
    void updateOofemState(int iteration, bool updateInternalForces);
    void formTotalExternalLoad(FloatArray &answer);
    void buildResidualRowScale(const FloatArray &externalForces, const FloatArray &elementNorm);
    void buildEquationDofIdMap(IntArray &equationDofIds) const;
    void buildEquationBounds(const FloatArray &referenceSolution);
    PetscErrorCode applyResidualRowScale(FloatArray &residual) const;
    PetscErrorCode applyJacobianRowScale(Mat jacobian) const;
    void checkPetscError(PetscErrorCode error, const char *operation) const;
    void invalidateLayoutCaches();

public:
    PetscSNESSolver(Domain *domain, EngngModel *model);
    ~PetscSNESSolver() override;

    ConvergedReason solve(SparseMtrx &matrix, FloatArray &referenceLoad,
                          FloatArray *initialLoad, FloatArray &solution,
                          FloatArray &increment, FloatArray &internalForces,
                          const FloatArray &internalForcesEBENorm,
                          double &loadLevel, referenceLoadInputModeType loadMode,
                          int &iterations, TimeStep *timeStep) override;

    void initializeFrom(const std :: shared_ptr< InputRecord > &record) override;
    SparseLinearSystemNM *giveLinearSolver() override;
    void setDomain(Domain *domain) override;
    void reinitialize() override;
    void printState(FILE *outputStream) override;
    void setStepReferenceSolution(const FloatArray &solution, TimeStep *step) override;
    bool prepareConvergenceResidual(FloatArray &residual, const FloatArray &solution,
                                    const FloatArray &external, const FloatArray &elementNorm) override;
    bool supportsCustomEquationNumbering() const override {
        // StationaryMPM does not currently provide a compact ExternalRhs
        // callback.  Fixed external vectors assembled before solve() are fine.
        return !solutionDependentExternalForces;
    }
    void setCustomEquationNumbering(const UnknownNumberingScheme *numbering) override;

    const char *giveClassName() const override { return "PetscSNESSolver"; }
    const char *giveInputRecordName() const { return _IFT_PetscSNESSolver_Name; }
};

} // namespace oofem

#endif // petscsnessolver_h
