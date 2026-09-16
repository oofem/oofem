/*
 * PETSc SNES adapter for OOFEM nonlinear problems.
 */

#include "petscsnessolver.h"

#include "classfactory.h"
#include "convergenceexception.h"
#include "dof.h"
#include "dofiditem.h"
#include "element.h"
#include "engngm.h"
#include "error.h"
#include "generalboundarycondition.h"
#include "parallelcontext.h"
#include "petscsparsemtrx.h"
#include "sparselinsystemnm.h"
#include "timestep.h"
#include "unknownnumberingscheme.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <sstream>
#include <utility>

namespace oofem {

namespace {

template< typename Function >
class PetscScopeExit
{
    Function action;

public:
    explicit PetscScopeExit(Function &&action) : action(std :: move(action)) { }
    PetscScopeExit(const PetscScopeExit &) = delete;
    PetscScopeExit &operator=(const PetscScopeExit &) = delete;
    ~PetscScopeExit() noexcept { action(); }
};

} // namespace

REGISTER_SparseNonLinearSystemNM(PetscSNESSolver)

PetscSNESSolver :: PetscSNESSolver(Domain *domain, EngngModel *model) :
    SparseNonLinearSystemNM(domain, model),
    maxIterations(50),
    maxFunctionEvaluations(10000),
    absoluteTolerance(1.e-50),
    relativeTolerance(1.e-8),
    stepTolerance(1.e-12),
    snesType(SNESNEWTONLS),
    monitor(true),
    solutionDependentExternalForces(false),
    fieldScaling(true),
    activeMatrix(nullptr),
    activeReferenceLoad(nullptr),
    activeInitialLoad(nullptr),
    activeSolution(nullptr),
    activeInternalForces(nullptr),
    activeInternalForcesEBENorm(nullptr),
    activeTimeStep(nullptr),
    residualRowScaleInitialized(false),
    functionEvaluations(0),
    jacobianEvaluations(0),
    domainErrors(0),
    lastReason(SNES_CONVERGED_ITERATING),
    lastIterations(0),
    lastFunctionNorm(0.0)
{
    deltaL = 1.0;
}


PetscSNESSolver :: ~PetscSNESSolver() = default;


void
PetscSNESSolver :: initializeFrom(const std :: shared_ptr< InputRecord > &record)
{
    SparseNonLinearSystemNM :: initializeFrom(record);

    // A NumericalMethod instance is reused across metasteps.  Each metastep
    // record is self-contained, so an omitted option must select the solver
    // default rather than silently inherit the preceding metastep's value.
    maxIterations = 50;
    maxFunctionEvaluations = 10000;
    absoluteTolerance = 1.e-50;
    relativeTolerance = 1.e-8;
    stepTolerance = 1.e-12;
    snesType = SNESNEWTONLS;
    optionsPrefix.clear();
    monitor = true;
    solutionDependentExternalForces = false;
    fieldScaling = true;
    minimumFieldScaleSquared = 1.e-6;
    boundedTrustRegion = false;
    trustInitialRadius = 1.0;
    trustMinimumRadius = 1.e-8;
    trustMaximumRadius = 1.e8;
    trustMaximumTrials = 20;
    boundDofIds.clear();
    lowerBounds.clear();
    upperBounds.clear();
    irreversibleDofIds.clear();
    fieldForceScale.clear();

    // Bounds and metrics are tied to both the option set and equation layout.
    // Rebuild them after every metastep reconfiguration.
    invalidateLayoutCaches();

    IR_GIVE_OPTIONAL_FIELD(record, maxIterations, _IFT_PetscSNESSolver_maxIterations);
    IR_GIVE_OPTIONAL_FIELD(record, maxFunctionEvaluations,
                           _IFT_PetscSNESSolver_maxFunctionEvaluations);
    IR_GIVE_OPTIONAL_FIELD(record, absoluteTolerance,
                           _IFT_PetscSNESSolver_absoluteTolerance);
    IR_GIVE_OPTIONAL_FIELD(record, relativeTolerance,
                           _IFT_PetscSNESSolver_relativeTolerance);
    IR_GIVE_OPTIONAL_FIELD(record, stepTolerance,
                           _IFT_PetscSNESSolver_stepTolerance);
    IR_GIVE_OPTIONAL_FIELD(record, snesType, _IFT_PetscSNESSolver_type);
    IR_GIVE_OPTIONAL_FIELD(record, optionsPrefix, _IFT_PetscSNESSolver_optionsPrefix);

    int monitorValue = monitor ? 1 : 0;
    IR_GIVE_OPTIONAL_FIELD(record, monitorValue, _IFT_PetscSNESSolver_monitor);
    monitor = monitorValue != 0;
    solutionDependentExternalForces =
        record->hasField(_IFT_PetscSNESSolver_solutionDependentExternalForces);

    int fieldScalingValue = fieldScaling ? 1 : 0;
    IR_GIVE_OPTIONAL_FIELD(record, fieldScalingValue,
                           _IFT_PetscSNESSolver_fieldScaling);
    fieldScaling = fieldScalingValue != 0;

    IR_GIVE_OPTIONAL_FIELD(record, boundedTrustRegion,
                           _IFT_PetscSNESSolver_boundedTrustRegion);
    IR_GIVE_OPTIONAL_FIELD(record, trustInitialRadius,
                           _IFT_PetscSNESSolver_trustInitialRadius);
    IR_GIVE_OPTIONAL_FIELD(record, trustMinimumRadius,
                           _IFT_PetscSNESSolver_trustMinimumRadius);
    IR_GIVE_OPTIONAL_FIELD(record, trustMaximumRadius,
                           _IFT_PetscSNESSolver_trustMaximumRadius);
    IR_GIVE_OPTIONAL_FIELD(record, trustMaximumTrials,
                           _IFT_PetscSNESSolver_trustMaximumTrials);
    if ( !std::isfinite(trustInitialRadius) || !std::isfinite(trustMinimumRadius) ||
         !std::isfinite(trustMaximumRadius) || trustMinimumRadius <= 0.0 ||
         trustInitialRadius < trustMinimumRadius || trustMaximumRadius < trustInitialRadius ||
         trustMaximumTrials < 1 ) {
        OOFEM_ERROR("invalid SNES trust-region radii or trial limit");
    }

    IR_GIVE_OPTIONAL_FIELD(record, minimumFieldScaleSquared,
                           _IFT_PetscSNESSolver_minimumFieldScaleSquared);
    IR_GIVE_OPTIONAL_FIELD(record, boundDofIds,
                           _IFT_PetscSNESSolver_boundDofIds);
    IR_GIVE_OPTIONAL_FIELD(record, lowerBounds,
                           _IFT_PetscSNESSolver_lowerBounds);
    IR_GIVE_OPTIONAL_FIELD(record, upperBounds,
                           _IFT_PetscSNESSolver_upperBounds);
    IR_GIVE_OPTIONAL_FIELD(record, irreversibleDofIds,
                           _IFT_PetscSNESSolver_irreversibleDofIds);
    if ( boundDofIds.giveSize() != lowerBounds.giveSize() ||
         boundDofIds.giveSize() != upperBounds.giveSize() ) {
        OOFEM_ERROR("SNES bound DOFs, lower bounds and upper bounds must have equal lengths");
    }
    for ( int i = 0; i < boundDofIds.giveSize(); ++i ) {
        if ( boundDofIds[i] <= 0 || !std::isfinite(lowerBounds[i]) ||
             !std::isfinite(upperBounds[i]) || lowerBounds[i] > upperBounds[i] ) {
            OOFEM_ERROR("invalid SNES variable bounds");
        }
        for ( int j = 0; j < i; ++j ) {
            if ( boundDofIds[i] == boundDofIds[j] ) {
                OOFEM_ERROR("duplicate SNES bound DOF");
            }
        }
    }
    for ( int id : irreversibleDofIds ) {
        if ( !boundDofIds.contains(id) ) {
            OOFEM_ERROR("every irreversible SNES DOF must have variable bounds");
        }
    }
    if ( !std::isfinite(minimumFieldScaleSquared) || minimumFieldScaleSquared <= 0.0 ) {
        OOFEM_ERROR(_IFT_PetscSNESSolver_minimumFieldScaleSquared
                    " must be finite and positive");
    }
    const bool hasForceScale = record->hasField(_IFT_PetscSNESSolver_forceScale);
    const bool hasForceScaleDofs =
        record->hasField(_IFT_PetscSNESSolver_forceScaleDofs);
    if ( hasForceScale != hasForceScaleDofs ) {
        OOFEM_ERROR(_IFT_PetscSNESSolver_forceScale " and "
                    _IFT_PetscSNESSolver_forceScaleDofs
                    " must be specified together");
    }
    if ( hasForceScale ) {
        IntArray dofs;
        FloatArray scales;
        IR_GIVE_FIELD(record, scales, _IFT_PetscSNESSolver_forceScale);
        IR_GIVE_FIELD(record, dofs, _IFT_PetscSNESSolver_forceScaleDofs);
        if ( dofs.giveSize() != scales.giveSize() ) {
            OOFEM_ERROR(_IFT_PetscSNESSolver_forceScale " and "
                        _IFT_PetscSNESSolver_forceScaleDofs
                        " must have equal lengths");
        }
        for ( int i = 0; i < dofs.giveSize(); ++i ) {
            if ( scales [ i ] <= 0.0 || !std :: isfinite(scales [ i ]) ) {
                OOFEM_ERROR(_IFT_PetscSNESSolver_forceScale
                            " values must be finite and positive");
            }
            fieldForceScale [ dofs [ i ] ] = scales [ i ];
        }
    }

    if ( maxIterations < 1 || maxFunctionEvaluations < 1 ) {
        OOFEM_ERROR("SNES iteration and function-evaluation limits must be positive");
    }
    if ( absoluteTolerance < 0.0 || relativeTolerance < 0.0 ||
         stepTolerance < 0.0 ) {
        OOFEM_ERROR("SNES tolerances must be non-negative");
    }
    if ( snesType.empty() ) {
        OOFEM_ERROR(_IFT_PetscSNESSolver_type " must not be empty");
    }

    // StaticStructural asks the nonlinear method for an OOFEM linear solver
    // before solve() when its tangent predictor is enabled.  This is separate
    // from the KSP owned internally by SNES and uses the same PetscSparseMtrx.
    this->giveLinearSolver()->initializeFrom(record);
}


SparseLinearSystemNM *
PetscSNESSolver :: giveLinearSolver()
{
    if ( !predictorLinearSolver ) {
        predictorLinearSolver = classFactory.createSparseLinSolver(
            ST_Petsc, domain, engngModel
        );
        if ( !predictorLinearSolver ) {
            OOFEM_ERROR("failed to create PETSc linear solver for tangent predictor");
        }
    }
    return predictorLinearSolver.get();
}


void
PetscSNESSolver :: setDomain(Domain *newDomain)
{
    SparseNonLinearSystemNM :: setDomain(newDomain);
    if ( predictorLinearSolver ) {
        predictorLinearSolver->setDomain(newDomain);
    }
    invalidateLayoutCaches();
}


void
PetscSNESSolver :: reinitialize()
{
    SparseNonLinearSystemNM :: reinitialize();
    if ( predictorLinearSolver ) {
        predictorLinearSolver->reinitialize();
    }
    invalidateLayoutCaches();
}


void
PetscSNESSolver :: invalidateLayoutCaches()
{
    boundReferenceStep = -1;
    boundReferenceSolution.clear();
    equationLowerBounds.clear();
    equationUpperBounds.clear();
    residualRowScale.clear();
    residualRowScaleInitialized = false;
    trustStepWeights.clear();
    customEquationNumbering = nullptr;
}


void
PetscSNESSolver :: buildEquationDofIdMap(IntArray &equationDofIds) const
{
    EModelDefaultEquationNumbering defaultNumbering;
    const UnknownNumberingScheme &numbering = customEquationNumbering ?
        *customEquationNumbering : defaultNumbering;
    ParallelContext *parallelContext =
        engngModel->giveParallelContext(domain->giveNumber());
    const int numberOfEquations = customEquationNumbering ?
        customEquationNumbering->giveRequiredNumberOfDomainEquation() :
        engngModel->giveNumberOfDomainEquations(
            domain->giveNumber(), defaultNumbering
        );
    equationDofIds.resize(numberOfEquations);
    equationDofIds.zero();

    auto recordDofManager = [&equationDofIds, &numbering](DofManager *manager) {
        for ( Dof *dof : *manager ) {
            if ( !dof->isPrimaryDof() ) {
                continue;
            }
            const int equation = dof->giveEquationNumber(numbering);
            if ( equation > 0 ) {
                equationDofIds.at(equation) = dof->giveDofID();
            }
        }
    };

    for ( auto &manager : domain->giveDofManagers() ) {
        if ( parallelContext->isLocal(manager.get()) ) {
            recordDofManager(manager.get());
        }
    }
    for ( auto &element : domain->giveElements() ) {
        if ( element->giveParallelMode() != Element_local ) {
            continue;
        }
        for ( int i = 1; i <= element->giveNumberOfInternalDofManagers(); ++i ) {
            recordDofManager(element->giveInternalDofManager(i));
        }
    }
    for ( auto &bc : domain->giveBcs() ) {
        for ( int i = 1; i <= bc->giveNumberOfInternalDofManagers(); ++i ) {
            recordDofManager(bc->giveInternalDofManager(i));
        }
    }
}


void
PetscSNESSolver :: buildEquationBounds(const FloatArray &referenceSolution)
{
    if ( !boundDofIds.giveSize() ) {
        equationLowerBounds.clear();
        equationUpperBounds.clear();
        return;
    }

    IntArray dofIds;
    buildEquationDofIdMap(dofIds);
    if ( dofIds.giveSize() != referenceSolution.giveSize() ) {
        OOFEM_ERROR("SNES bound reference layout is incompatible with the active numbering");
    }

    equationLowerBounds.resize(referenceSolution.giveSize());
    equationUpperBounds.resize(referenceSolution.giveSize());
    for ( int i = 0; i < referenceSolution.giveSize(); ++i ) {
        equationLowerBounds[i] = PETSC_NINFINITY;
        equationUpperBounds[i] = PETSC_INFINITY;
        for ( int j = 0; j < boundDofIds.giveSize(); ++j ) {
            if ( dofIds[i] != boundDofIds[j] ) {
                continue;
            }
            equationLowerBounds[i] = lowerBounds[j];
            equationUpperBounds[i] = upperBounds[j];
            if ( irreversibleDofIds.contains(dofIds[i]) ) {
                equationLowerBounds[i] = std::max(
                    equationLowerBounds[i], referenceSolution[i]
                );
            }
            if ( equationLowerBounds[i] > equationUpperBounds[i] ) {
                OOFEM_ERROR("irreversible SNES reference exceeds its upper bound");
            }
            break;
        }
    }
}


void
PetscSNESSolver :: buildResidualRowScale(const FloatArray &externalForces, const FloatArray &elementNorm)
{
    IntArray equationDofIds;
    buildEquationDofIdMap(equationDofIds);
    if ( equationDofIds.giveSize() != externalForces.giveSize() ) {
        OOFEM_ERROR("SNES equation/DOF map has incompatible size");
    }

    const int numberOfDofIds = domain->giveMaxDofID();
    FloatArray externalNorm2(numberOfDofIds);
    FloatArray equationCount(numberOfDofIds);
    externalNorm2.zero();
    equationCount.zero();
    for ( int equation = 1; equation <= equationDofIds.giveSize(); ++equation ) {
        const int dofId = equationDofIds.at(equation);
        if ( dofId <= 0 ) {
            continue;
        }
        const double value = externalForces.at(equation);
        externalNorm2.at(dofId) += value * value;
        equationCount.at(dofId) += 1.0;
    }

    ParallelContext *parallelContext =
        engngModel->giveParallelContext(domain->giveNumber());
    FloatArray collective(numberOfDofIds);
    parallelContext->accumulate(externalNorm2, collective);
    externalNorm2 = collective;
    parallelContext->accumulate(equationCount, collective);
    equationCount = collective;

    FloatArray dofRowScale(numberOfDofIds);
    dofRowScale.zero();
    for ( int dofId = 1; dofId <= numberOfDofIds; ++dofId ) {
        if ( equationCount.at(dofId) == 0.0 ) {
            continue;
        }
        double scale2 = externalNorm2.at(dofId);
        if ( dofId <= elementNorm.giveSize() ) {
            scale2 += elementNorm.at(dofId);
        }
        auto explicitScale = fieldForceScale.find(dofId);
        if ( explicitScale != fieldForceScale.end() ) {
            scale2 += equationCount.at(dofId) *
                      explicitScale->second * explicitScale->second;
        }
        if ( !std :: isfinite(scale2) || scale2 < 0.0 ) {
            OOFEM_ERROR("invalid SNES residual scale for DOF group %d", dofId);
        }
        // This is sqrt(buildResidualWeights()) from NRSolver. Applying it
        // to both f and the rows of J preserves the exact Newton correction.
        dofRowScale.at(dofId) =
            scale2 >= minimumFieldScaleSquared ? 1.0 / std :: sqrt(scale2) : 1.0;
        if ( this->giveIterationLogEnabled() &&
             engngModel->giveProblemScale() == macroScale ) {
            OOFEM_LOG_INFO(
                "PetscSNESSolver: field %-8s equations=%d scale2=%.6e "
                "row-scale=%.6e\n",
                DofIDItemToString(static_cast< DofIDItem >(dofId)).c_str(),
                static_cast< int >(equationCount.at(dofId)), scale2,
                dofRowScale.at(dofId)
            );
        }
    }

    residualRowScale.resize(equationDofIds.giveSize());
    for ( int equation = 1; equation <= equationDofIds.giveSize(); ++equation ) {
        const int dofId = equationDofIds.at(equation);
        residualRowScale.at(equation) = dofId > 0 ? dofRowScale.at(dofId) : 1.0;
    }
}


PetscErrorCode
PetscSNESSolver :: applyResidualRowScale(FloatArray &residual) const
{
    if ( !fieldScaling ) {
        return PETSC_SUCCESS;
    }
    if ( residual.giveSize() != residualRowScale.giveSize() ) {
        return PETSC_ERR_ARG_SIZ;
    }
    for ( int equation = 1; equation <= residual.giveSize(); ++equation ) {
        residual.at(equation) *= residualRowScale.at(equation);
    }
    return PETSC_SUCCESS;
}


PetscErrorCode
PetscSNESSolver :: applyJacobianRowScale(Mat jacobian) const
{
    if ( !fieldScaling ) {
        return PETSC_SUCCESS;
    }
    Vec scale = nullptr;
    activeMatrix->createVecGlobal(&scale);
    // PetscSparseMtrx follows the historical OOFEM convention and returns
    // one on successful scatter, not PETSC_SUCCESS (zero).
    activeMatrix->scatterL2G(residualRowScale, scale);
    PetscErrorCode error = MatDiagonalScale(jacobian, scale, nullptr);
    VecDestroy(&scale);
    return error;
}


void
PetscSNESSolver :: checkPetscError(PetscErrorCode error,
                                   const char *operation) const
{
    if ( error != PETSC_SUCCESS ) {
        OOFEM_ERROR("PETSc operation '%s' failed with error code %d",
                    operation, static_cast< int >(error));
    }
}


void
PetscSNESSolver :: formTotalExternalLoad(FloatArray &answer)
{
    answer = *activeReferenceLoad;
    if ( activeInitialLoad ) {
        answer.add(*activeInitialLoad);
    }
}


void
PetscSNESSolver :: applyTrialSolution(SNES snes, Vec x,
                                      bool updateInternalForces)
{
    activeMatrix->scatterG2L(x, *activeSolution);

    PetscInt iteration = 0;
    SNESGetIterationNumber(snes, &iteration);
    updateOofemState(static_cast< int >(iteration), updateInternalForces);
}


void
PetscSNESSolver :: updateOofemState(int iteration,
                                    bool updateInternalForces)
{

    // SNES may evaluate several line-search/trust-region trial points in one
    // outer iteration. Change the state counter before updating fields and
    // contact so no cached state can be reused at another trial point.
    activeTimeStep->incrementStateCounter();

    engngModel->initForNewIteration(
        domain, activeTimeStep, iteration, *activeSolution
    );

    if ( updateInternalForces ) {
        engngModel->updateComponent(activeTimeStep, InternalRhs, domain);
    }
}


PetscErrorCode
PetscSNESSolver :: evaluateFunction(SNES snes, Vec x, Vec f)
{
    applyTrialSolution(snes, x, true);

    if ( solutionDependentExternalForces ) {
        engngModel->updateComponent(activeTimeStep, ExternalRhs, domain);
    }

    FloatArray externalForces;
    formTotalExternalLoad(externalForces);
    FloatArray residual = *activeInternalForces;
    residual.subtract(externalForces);
    // MPM updates the element-by-element force norm together with InternalRhs.
    // At entry to solve() that norm can still be empty, so freeze the field
    // metric only after the first admissible residual evaluation.
    if ( fieldScaling && !residualRowScaleInitialized ) {
        buildResidualRowScale(externalForces, *activeInternalForcesEBENorm);
        residualRowScaleInitialized = true;
    }
    PetscErrorCode error = applyResidualRowScale(residual);
    if ( error != PETSC_SUCCESS ) {
        return error;
    }
    activeMatrix->scatterL2G(residual, f);
    ++functionEvaluations;
    return PETSC_SUCCESS;
}


PetscErrorCode
PetscSNESSolver :: evaluateJacobian(SNES snes, Vec x, Mat jacobian)
{
    // Re-evaluate the internal state at x before asking materials for their
    // tangent. This is deliberately conservative: SNES is allowed to request
    // a Jacobian without a preceding function evaluation at the same x.
    applyTrialSolution(snes, x, true);
    if ( fieldScaling && !residualRowScaleInitialized ) {
        FloatArray externalForces;
        formTotalExternalLoad(externalForces);
        buildResidualRowScale(externalForces, *activeInternalForcesEBENorm);
        residualRowScaleInitialized = true;
    }
    engngModel->updateComponent(activeTimeStep, NonLinearLhs, domain);
    // Symbolic MPM integrals insert entries directly and do not pass through
    // EngngModel::assemble, which normally finalizes a PETSc matrix. KSP may
    // only inspect or factor a matrix after this collective final assembly.
    PetscErrorCode error = activeMatrix->assembleBegin();
    if ( error != PETSC_SUCCESS ) {
        return error;
    }
    error = activeMatrix->assembleEnd();
    if ( error != PETSC_SUCCESS ) {
        return error;
    }
    error = MatCopy(*activeMatrix->giveMtrx(), jacobian, SAME_NONZERO_PATTERN);
    if ( error != PETSC_SUCCESS ) {
        return error;
    }
    error = applyJacobianRowScale(jacobian);
    if ( error != PETSC_SUCCESS ) {
        return error;
    }
    ++jacobianEvaluations;
    return PETSC_SUCCESS;
}


PetscErrorCode
PetscSNESSolver :: formFunction(SNES snes, Vec x, Vec f, void *context)
{
    auto *solver = static_cast< PetscSNESSolver * >(context);
    try {
        return solver->evaluateFunction(snes, x, f);
    } catch ( const ConvergenceException &exception ) {
        solver->lastDomainError = exception.what();
        ++solver->domainErrors;
        VecSet(f, PETSC_INFINITY);
        return SNESSetFunctionDomainError(snes);
    } catch ( const std :: exception &exception ) {
        solver->callbackFailure = exception.what();
        return PETSC_ERR_USER;
    } catch ( ... ) {
        solver->callbackFailure = "unknown exception in OOFEM residual callback";
        return PETSC_ERR_USER;
    }
}


PetscErrorCode
PetscSNESSolver :: formJacobian(SNES snes, Vec x, Mat jacobian,
                                Mat preconditioner, void *context)
{
    ( void ) jacobian;
    ( void ) preconditioner;
    auto *solver = static_cast< PetscSNESSolver * >(context);
    try {
        return solver->evaluateJacobian(snes, x, jacobian);
    } catch ( const ConvergenceException &exception ) {
        solver->lastDomainError = exception.what();
        ++solver->domainErrors;
        return SNESSetJacobianDomainError(snes);
    } catch ( const std :: exception &exception ) {
        solver->callbackFailure = exception.what();
        return PETSC_ERR_USER;
    } catch ( ... ) {
        solver->callbackFailure = "unknown exception in OOFEM Jacobian callback";
        return PETSC_ERR_USER;
    }
}


PetscErrorCode
PetscSNESSolver :: monitorIteration(SNES, PetscInt iteration,
                                    PetscReal functionNorm, void *context)
{
    auto *solver = static_cast< PetscSNESSolver * >(context);
    if ( solver->giveIterationLogEnabled() &&
         solver->engngModel->giveProblemScale() == macroScale ) {
        OOFEM_LOG_INFO("PetscSNESSolver: iter %d, ||f|| = %.6e\n",
                       static_cast< int >(iteration),
                       static_cast< double >(functionNorm));
    }
    return PETSC_SUCCESS;
}


PetscErrorCode
PetscSNESSolver :: applyBoundedTrustRegion(SNESLineSearch search, void *context)
{
    auto *solver = static_cast<PetscSNESSolver *>(context);
    try {
        return solver->boundedTrustRegionStep(search);
    } catch (const std::exception &exception) {
        solver->callbackFailure = exception.what();
        return PETSC_ERR_USER;
    } catch (...) {
        solver->callbackFailure = "unknown exception in bounded trust-region callback";
        return PETSC_ERR_USER;
    }
}


PetscErrorCode
PetscSNESSolver :: boundedTrustRegionStep(SNESLineSearch search)
{
    SNES snes;
    Vec x, f, newton, trial, trialResidual;
    Mat jacobian;
#if PETSC_VERSION_LT(3, 21, 0)
    SNESLineSearchVIProjectFunc project = nullptr;
    SNESLineSearchVINormFunc viNorm = nullptr;
#else
    SNESLineSearchVIProjectFn *project = nullptr;
    SNESLineSearchVINormFn *viNorm = nullptr;
#endif
    PetscCall(SNESLineSearchGetSNES(search, &snes));
    PetscCall(SNESLineSearchGetVecs(search, &x, &f, &newton, &trial, &trialResidual));
#if PETSC_VERSION_LT(3, 23, 0)
    PetscCall(SNESLineSearchGetVIFunctions(search, &project, &viNorm));
#else
    PetscCall(SNESLineSearchGetVIFunctions(search, &project, &viNorm, nullptr));
#endif
    PetscCheck(project && viNorm, PetscObjectComm((PetscObject)snes), PETSC_ERR_ARG_WRONGSTATE,
               "bounded trust region requires the PETSc VI projection and norm");
    PetscCall(SNESGetJacobian(snes, &jacobian, nullptr, nullptr, nullptr));
    FloatArray direction;
    activeMatrix->scatterG2L(newton, direction);
    auto *parallel = engngModel->giveParallelContext(domain->giveNumber());
    if (trustStepWeights.giveSize() == 0) {
        // Freeze a balanced, dimensionless metric for this nonlinear solve.
        // As in NRSolver, the initial full Newton step has metric norm one.
        IntArray ids;
        buildEquationDofIdMap(ids);
        FloatArray sums(domain->giveMaxDofID()), global;
        sums.zero();
        for (int i=0; i<ids.giveSize(); ++i) {
            if (ids[i] > 0) sums.at(ids[i]) += direction[i]*direction[i];
        }
        parallel->accumulate(sums, global);
        int fields=0;
        for (double value : global) {
            if (std::isfinite(value) && value > std::numeric_limits<double>::min()) ++fields;
        }
        trustStepWeights.resize(direction.giveSize());
        trustStepWeights.zero();
        for (int i=0; i<ids.giveSize(); ++i) {
            if (fields && ids[i]>0 && global.at(ids[i]) > std::numeric_limits<double>::min()) {
                trustStepWeights[i]=1.0/(fields*global.at(ids[i]));
            }
        }
    }
    auto metricNorm = [&](const FloatArray &step) {
        FloatArray weighted=step;
        for (int i=0; i<step.giveSize(); ++i) weighted[i]*=trustStepWeights[i];
        return std::sqrt(std::max(0.0, parallel->localDotProduct(step, weighted)));
    };
    const double newtonNorm=metricNorm(direction);
    PetscReal baseNorm;
    PetscCall(viNorm(snes, f, x, &baseNorm));
    const double baseMerit=0.5*baseNorm*baseNorm;
    // RAII also releases work storage if a PETSc operation or callback fails.
    struct WorkVector { Vec v=nullptr; ~WorkVector() { VecDestroy(&v); } } actualStep;
    PetscCall(VecDuplicate(x, &actualStep.v));
    for (int attempt=0; attempt<trustMaximumTrials; ++attempt) {
        const double attemptedRadius=trustRadius;
        const double lambda=newtonNorm > attemptedRadius ? attemptedRadius/newtonNorm : 1.0;
        PetscCall(VecWAXPY(trial, -lambda, newton, x));
        PetscCall(project(snes, trial)); // bounds before any material evaluation
        PetscCall(VecWAXPY(actualStep.v, -1.0, trial, x));
        FloatArray actual;
        activeMatrix->scatterG2L(actualStep.v, actual);
        const double stepNorm=metricNorm(actual);
        // Predict using the actual projected step and the same VI norm at
        // the trial point. The physical Jacobian retained by OOFEM is untouched.
        PetscCall(MatMult(jacobian, actualStep.v, trialResidual));
        // VecWAXPY above gives actualStep = x - trial, i.e. PETSc's
        // (positive) correction. Thus f(trial) ~= f(x) - J*actualStep.
        PetscCall(VecAYPX(trialResidual, -1.0, f));
        PetscReal modelNorm;
        PetscCall(viNorm(snes, trialResidual, trial, &modelNorm));
        const double predicted=baseMerit-0.5*modelNorm*modelNorm;
        PetscCall(SNESComputeFunction(snes, trial, trialResidual));
        // Our residual callback fills f with infinity on a domain error.
        // Check the full vector before VI masking can hide an invalid entry.
        PetscReal fullNorm, trialNorm=PETSC_INFINITY;
        PetscCall(VecNorm(trialResidual, NORM_2, &fullNorm));
        if (std::isfinite(fullNorm)) PetscCall(viNorm(snes, trialResidual, trial, &trialNorm));
        const double reduction=baseMerit-0.5*trialNorm*trialNorm;
        const double ratio=predicted>0.0 && std::isfinite(trialNorm) ? reduction/predicted : -1.0;
        const bool accepted=std::isfinite(ratio) && ratio>1.e-4 && reduction>0.0;
        if (!std::isfinite(ratio) || ratio<0.25) {
            trustRadius=std::max(trustMinimumRadius, 0.25*attemptedRadius);
        } else if (ratio>0.75 && stepNorm>=0.8*attemptedRadius) {
            trustRadius=std::min(trustMaximumRadius, 2.0*attemptedRadius);
        }
        if ( monitor && this->giveIterationLogEnabled() &&
             engngModel->giveProblemScale() == macroScale ) {
            OOFEM_LOG_INFO("SNES VI trust region: %s, lambda=%.6e, radius=%.6e, "
                           "merit=%.6e -> %.6e, rho=%.6e\n",
                           accepted ? "accepted" : "rejected", lambda, attemptedRadius,
                           baseMerit, 0.5*trialNorm*trialNorm, ratio);
        }
        if (accepted) {
            PetscCall(VecCopy(trial, x));
            PetscCall(VecCopy(trialResidual, f));
            PetscCall(VecCopy(actualStep.v, newton));
            PetscCall(SNESLineSearchSetLambda(search, lambda));
            PetscCall(SNESLineSearchSetReason(search, SNES_LINESEARCH_SUCCEEDED));
            PetscCall(SNESLineSearchComputeNorms(search));
            return PETSC_SUCCESS;
        }
        if (attemptedRadius<=trustMinimumRadius) break;
    }
    // Restore the baseline material state as well as x/f before reporting
    // failure; a failed solve will subsequently restore the physical step base.
    PetscCall(SNESComputeFunction(snes, x, f));
    PetscCall(SNESLineSearchSetReason(search, SNES_LINESEARCH_FAILED_REDUCT));
    PetscCall(SNESLineSearchComputeNorms(search));
    return PETSC_SUCCESS;
}


ConvergedReason
PetscSNESSolver :: solve(SparseMtrx &matrix, FloatArray &referenceLoad,
                         FloatArray *initialLoad, FloatArray &solution,
                         FloatArray &increment, FloatArray &internalForces,
                         const FloatArray &internalForcesEBENorm,
                         double &loadLevel, referenceLoadInputModeType loadMode,
                         int &iterations, TimeStep *timeStep)
{
    ( void ) loadMode;

    auto *petscMatrix = dynamic_cast< PetscSparseMtrx * >(&matrix);
    if ( !petscMatrix ) {
        OOFEM_ERROR("PetscSNESSolver requires smtype %d (PetscSparseMtrx)",
                    static_cast< int >(SMT_PetscMtrx));
    }
    if ( matrix.giveNumberOfRows() != solution.giveSize() ||
         referenceLoad.giveSize() != solution.giveSize() ||
         internalForces.giveSize() != solution.giveSize() ) {
        OOFEM_ERROR("incompatible matrix/vector dimensions passed to PetscSNESSolver");
    }
    if ( increment.giveSize() != 0 &&
         increment.giveSize() != solution.giveSize() ) {
        OOFEM_ERROR("incompatible solution-increment dimension passed to PetscSNESSolver");
    }
    if ( engngModel->isParallel() ) {
        // A rank-local material exception can occur before OOFEM's residual or
        // matrix assembly reaches its own collectives. PETSc only synchronizes
        // domain-error flags after every callback has returned, which is too
        // late to prevent mismatched collectives inside that callback.
        OOFEM_ERROR("PetscSNESSolver currently supports serial OOFEM analyses only");
    }

    auto clearActiveState = [this]() noexcept {
        activeMatrix = nullptr;
        activeReferenceLoad = nullptr;
        activeInitialLoad = nullptr;
        activeSolution = nullptr;
        activeInternalForces = nullptr;
        activeInternalForcesEBENorm = nullptr;
        activeTimeStep = nullptr;
    };
    PetscScopeExit< decltype(clearActiveState) > activeStateCleanup(
        std :: move(clearActiveState)
    );

    activeMatrix = petscMatrix;
    activeReferenceLoad = &referenceLoad;
    activeInitialLoad = initialLoad;
    activeSolution = &solution;
    activeInternalForces = &internalForces;
    activeInternalForcesEBENorm = &internalForcesEBENorm;
    activeTimeStep = timeStep;
    functionEvaluations = 0;
    jacobianEvaluations = 0;
    domainErrors = 0;
    callbackFailure.clear();
    lastDomainError.clear();
    lastReason = SNES_CONVERGED_ITERATING;
    lastIterations = 0;
    lastFunctionNorm = 0.0;
    residualRowScaleInitialized = false;

    Vec x = nullptr;
    Vec residual = nullptr;
    Vec lower = nullptr, upper = nullptr;
    Mat jacobian = nullptr;
    SNES snes = nullptr;

    struct PetscObjectCleanup {
        SNES &snes;
        Vec &lower;
        Vec &upper;
        Mat &jacobian;
        Vec &residual;
        Vec &x;
        bool cleaned = false;

        void cleanup() noexcept
        {
            if ( cleaned ) {
                return;
            }
            SNESDestroy(&snes);
            VecDestroy(&lower);
            VecDestroy(&upper);
            MatDestroy(&jacobian);
            VecDestroy(&residual);
            VecDestroy(&x);
            cleaned = true;
        }

        ~PetscObjectCleanup() noexcept { cleanup(); }
    } petscObjectCleanup { snes, lower, upper, jacobian, residual, x };

    // The engineering model may already have added a tangent predictor to X
    // and increment.  SNES starts from that predicted X, but successful output
    // and failed-solve rollback are both relative to the accepted step base.
    FloatArray stepBaseSolution = solution;
    if ( increment.giveSize() ) {
        stepBaseSolution.subtract(increment);
    }
    auto restoreStepBase = [&]() {
        // Do not write x here: a PETSc callback that returned an error may
        // leave that vector read-locked until SNES is destroyed.
        solution = stepBaseSolution;
        increment.resize(stepBaseSolution.giveSize());
        increment.zero();
        engngModel->initStepIncrements();
        updateOofemState(0, true);
        if ( solutionDependentExternalForces ) {
            engngModel->updateComponent(activeTimeStep, ExternalRhs, domain);
        }
    };

    try {
        if ( boundDofIds.giveSize() && boundReferenceStep != timeStep->giveNumber() ) {
            boundReferenceSolution = stepBaseSolution;
            boundReferenceStep = timeStep->giveNumber();
            equationLowerBounds.clear();
            equationUpperBounds.clear();
        }
        loadLevel = 1.0;
        deltaL = 1.0;

        if ( !fieldScaling ) {
            residualRowScale.resize(solution.giveSize());
            residualRowScale.zero();
            residualRowScale.add(1.0);
            residualRowScaleInitialized = true;
        } else {
            residualRowScale.resize(0);
        }

        activeMatrix->createVecGlobal(&x);
        activeMatrix->scatterL2G(solution, x);
        checkPetscError(VecDuplicate(x, &residual), "VecDuplicate");

        MPI_Comm communicator = PetscObjectComm(
            reinterpret_cast< PetscObject >(*activeMatrix->giveMtrx())
        );

        // Populate the OOFEM matrix once before duplicating its PETSc sparsity
        // pattern. SNES uses the duplicate so row scaling never contaminates the
        // physical tangent retained by the engineering model.
        activeTimeStep->incrementStateCounter();
        engngModel->initForNewIteration(domain, activeTimeStep, 0, solution);
        engngModel->updateComponent(activeTimeStep, NonLinearLhs, domain);
        checkPetscError(activeMatrix->assembleBegin(), "initial MatAssemblyBegin");
        checkPetscError(activeMatrix->assembleEnd(), "initial MatAssemblyEnd");
        checkPetscError(MatDuplicate(*activeMatrix->giveMtrx(),
                                     MAT_DO_NOT_COPY_VALUES, &jacobian),
                        "MatDuplicate(SNES Jacobian)");

        checkPetscError(SNESCreate(communicator, &snes), "SNESCreate");
        // Optimized PETSc builds do not check SNESSetJacobianDomainError() by
        // default. The OOFEM callback uses it for inadmissible material states,
        // so enable the check explicitly (SNESSetFromOptions may still override).
        checkPetscError(SNESSetCheckJacobianDomainError(snes, PETSC_TRUE),
                        "SNESSetCheckJacobianDomainError");
        checkPetscError(SNESSetFunction(snes, residual, formFunction, this),
                        "SNESSetFunction");
        checkPetscError(SNESSetJacobian(snes, jacobian, jacobian,
                                        formJacobian, this),
                        "SNESSetJacobian");
        checkPetscError(SNESSetType(snes, snesType.c_str()), "SNESSetType");
        checkPetscError(SNESSetTolerances(
            snes, absoluteTolerance, relativeTolerance, stepTolerance,
            maxIterations, maxFunctionEvaluations
        ), "SNESSetTolerances");
        checkPetscError(SNESSetErrorIfNotConverged(snes, PETSC_FALSE),
                        "SNESSetErrorIfNotConverged");
        if ( !optionsPrefix.empty() ) {
            checkPetscError(SNESSetOptionsPrefix(snes, optionsPrefix.c_str()),
                            "SNESSetOptionsPrefix");
        }
        if ( monitor ) {
            checkPetscError(SNESMonitorSet(snes, monitorIteration, this, nullptr),
                            "SNESMonitorSet");
        }
        // Command-line PETSc options deliberately come last and override all
        // defaults above, including SNES/KSP/PC type and tolerances.
        checkPetscError(SNESSetFromOptions(snes), "SNESSetFromOptions");
        boundTolerance = 1.e-8; // PETSc VI default; honor the same options prefix.
        checkPetscError(PetscOptionsGetReal(nullptr, optionsPrefix.c_str(),
                                          "-snes_vi_zero_tolerance", &boundTolerance, nullptr),
                        "PetscOptionsGetReal(VI bound tolerance)");

        if ( boundDofIds.giveSize() ) {
            SNESType selectedType;
            checkPetscError(SNESGetType(snes, &selectedType), "SNESGetType");
            if ( std::string(selectedType) != SNESVINEWTONRSLS &&
                 std::string(selectedType) != SNESVINEWTONSSLS ) {
                OOFEM_ERROR("SNES variable bounds require vinewtonrsls or vinewtonssls");
            }
            if ( boundReferenceSolution.giveSize() != solution.giveSize() ) {
                OOFEM_ERROR("SNES bound reference layout changed within a step");
            }
            buildEquationBounds(boundReferenceSolution);
            checkPetscError(VecDuplicate(x, &lower), "VecDuplicate(lower bounds)");
            checkPetscError(VecDuplicate(x, &upper), "VecDuplicate(upper bounds)");
            activeMatrix->scatterL2G(equationLowerBounds, lower);
            activeMatrix->scatterL2G(equationUpperBounds, upper);
            checkPetscError(SNESVISetVariableBounds(snes, lower, upper), "SNESVISetVariableBounds");
        }

        if (boundedTrustRegion) {
            SNESType selectedType;
            checkPetscError(SNESGetType(snes, &selectedType), "SNESGetType(trust region)");
            if (!boundDofIds.giveSize() || std::string(selectedType)!=SNESVINEWTONRSLS) {
                OOFEM_ERROR(_IFT_PetscSNESSolver_boundedTrustRegion
                            " requires bounded vinewtonrsls");
            }
            SNESLineSearch search;
            checkPetscError(SNESGetLineSearch(snes, &search), "SNESGetLineSearch(trust region)");
            checkPetscError(SNESLineSearchSetType(search, SNESLINESEARCHSHELL), "SNESLineSearchSetType(shell)");
#if PETSC_VERSION_LT(3, 21, 0)
            checkPetscError(SNESLineSearchShellSetUserFunc(search, applyBoundedTrustRegion, this),
                            "SNESLineSearchShellSetUserFunc(trust region)");
#else
            checkPetscError(SNESLineSearchShellSetApply(search, applyBoundedTrustRegion, this),
                            "SNESLineSearchShellSetApply(trust region)");
#endif
            trustRadius=trustInitialRadius;
            trustStepWeights.clear();
        }

        const PetscErrorCode solveError = SNESSolve(snes, nullptr, x);
        if ( solveError == PETSC_SUCCESS ) {
            checkPetscError(SNESGetConvergedReason(snes, &lastReason),
                            "SNESGetConvergedReason");
            PetscInt petscIterations = 0;
            checkPetscError(SNESGetIterationNumber(snes, &petscIterations),
                            "SNESGetIterationNumber");
            lastIterations = static_cast< int >(petscIterations);
            PetscReal functionNorm = 0.0;
            checkPetscError(SNESGetFunctionNorm(snes, &functionNorm),
                            "SNESGetFunctionNorm");
            lastFunctionNorm = static_cast< double >(functionNorm);
        }

        // SNES owns x during the solve. Leave OOFEM at the accepted state only
        // after convergence. Every rejected solve and every exception below is
        // restored to the exact step base by the common catch block.
        if ( solveError == PETSC_SUCCESS && lastReason > 0 ) {
            activeMatrix->scatterG2L(x, solution);
            updateOofemState(lastIterations, true);
            if ( solutionDependentExternalForces ) {
                engngModel->updateComponent(activeTimeStep, ExternalRhs, domain);
            }
            increment = solution;
            increment.subtract(stepBaseSolution);
        }
        iterations = lastIterations;

        if ( solveError != PETSC_SUCCESS ) {
            std :: ostringstream message;
            message << "PETSc SNESSolve failed with error code " << solveError;
            if ( !callbackFailure.empty() ) {
                message << ": " << callbackFailure;
            }
            throw ConvergenceException(message.str());
        }

        const char *reasonName = SNESConvergedReasons[lastReason];
        // In compact staggered output the outer sweep table reports this solve's
        // iteration count. Keep a divergent reason visible, but avoid repeating a
        // successful subproblem summary before every outer table.
        if ( this->giveIterationLogEnabled() || lastReason <= 0 ) {
            OOFEM_LOG_INFO(
                "PetscSNESSolver: %s after %d iterations, %d residual and %d Jacobian "
                "evaluations, final ||f||=%.6e, domain errors=%d\n",
                reasonName, lastIterations, functionEvaluations, jacobianEvaluations,
                lastFunctionNorm, domainErrors
            );
        }

        if ( lastReason <= 0 ) {
            std :: ostringstream message;
            message << "PETSc SNES diverged: " << reasonName;
            if ( !lastDomainError.empty() ) {
                message << "; last inadmissible state: " << lastDomainError;
            }
            if ( !callbackFailure.empty() ) {
                message << "; callback failure: " << callbackFailure;
            }
            throw ConvergenceException(message.str());
        }

        return CR_CONVERGED;
    } catch ( ... ) {
        std :: exception_ptr originalFailure = std :: current_exception();
        try {
            restoreStepBase();
        } catch ( const std :: exception &exception ) {
            OOFEM_WARNING("PetscSNESSolver failed to restore the step-base state: %s",
                          exception.what());
        } catch ( ... ) {
            OOFEM_WARNING("PetscSNESSolver failed to restore the step-base state");
        }
        std :: rethrow_exception(originalFailure);
    }
}


void
PetscSNESSolver :: setStepReferenceSolution(const FloatArray &solution, TimeStep *step)
{
    boundReferenceSolution = solution;
    boundReferenceStep = step->giveNumber();
    equationLowerBounds.clear();
    equationUpperBounds.clear();
}


void
PetscSNESSolver :: setCustomEquationNumbering(const UnknownNumberingScheme *numbering)
{
    customEquationNumbering = numbering;
    // These arrays are indexed by the selected equation layout and must not
    // survive a switch between a compact group and the default full system.
    equationLowerBounds.clear();
    equationUpperBounds.clear();
    residualRowScale.clear();
    residualRowScaleInitialized = false;
    trustStepWeights.clear();
}

bool
PetscSNESSolver :: prepareConvergenceResidual(FloatArray &residual, const FloatArray &solution,
                                             const FloatArray &external, const FloatArray &elementNorm)
{
    if ( boundDofIds.giveSize() &&
         equationLowerBounds.giveSize() != solution.giveSize() ) {
        if ( boundReferenceSolution.giveSize() != solution.giveSize() ) {
            OOFEM_ERROR("SNES convergence bound reference has incompatible size");
        }
        buildEquationBounds(boundReferenceSolution);
    }
    for ( int i = 0; boundDofIds.giveSize() && i < solution.giveSize(); ++i ) {
        if ( (solution[i] <= equationLowerBounds[i] + boundTolerance && residual[i] > 0.0) ||
             (solution[i] >= equationUpperBounds[i] - boundTolerance && residual[i] < 0.0) ) {
            residual[i] = 0.0;
        }
    }
    if ( fieldScaling ) {
        buildResidualRowScale(external, elementNorm);
        checkPetscError(applyResidualRowScale(residual), "outer residual scaling");
    }
    return fieldScaling;
}

void
PetscSNESSolver :: printState(FILE *outputStream)
{
    if ( lastReason == SNES_CONVERGED_ITERATING ) {
        return;
    }
    std :: fprintf(outputStream,
                   "PETSc SNES reason: %s; iterations: %d; final norm: %.6e\n",
                   SNESConvergedReasons[lastReason], lastIterations,
                   lastFunctionNorm);
}

} // namespace oofem
