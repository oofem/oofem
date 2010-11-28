/**
 * @date 2010-09-06
 * @author Mikael Öhman
 */

#ifndef stokesflowstressrve_h
#define stokesflowstressrve_h

#include "inputrecord.h"
#include "stokesflow.h"

#include "sparsemtrxtype.h"
#include "sparselinsystemnm.h"
#include "linsystsolvertype.h"

namespace oofem {
enum HomogenizationType {
    HT_None = 0,
    HT_StressDirichlet = 1, /// Dirichlet boundary condition
    HT_StressNeumann = 2, /// Neumann boundary condition
    HT_StressPeriodic = 3, /// Periodic boundary condition
};


/** Stokes flow model with homogenization for stress.
 * @author Mikael Öhman
 */
class StokesFlowStressHomogenization : public StokesFlow
{
protected:
    HomogenizationType ht;

    /// Non-free submatrices of the total stiffness matrix.
    SparseMtrx *K_cf, *K_cc;
    /// Coefficient matrix for values of prescribed dofs.
    FloatMatrix *C;
    /// Solver type for homogenization
    LinSystSolverType solverType;
    /// Numerical method for homogenization
    SparseLinearSystemNM *linNumericalMethod;

public:
    StokesFlowStressHomogenization(int i, EngngModel *_master = NULL);
    virtual ~StokesFlowStressHomogenization();

    /** Initialization from given input record.
     * Reads
     * - lstype Linear solver type (enum, optional, default ST_Petsc). Should conform to the non-linear solver.
     * - homogenizationtype Sets which type of homogenization to use for multiscale analysis (enum, optional, default HT_None)
     * @see StokesFlow::initializeFrom
     */
    IRResultType initializeFrom(InputRecord *ir);

    /** Updates the coefficient matrix C.
     * @param tStep Current time step
     */
    virtual void updateYourself(TimeStep *tStep) { StokesFlow :: updateYourself(tStep);
                                                   this->updateC(); };

    /**
     * Computes the macroscopic stress.
     * @param answer Macroscopic stress
     * @param input Strain rate
     * @param tStep Time step to evaluate for
     * @return True if successfull.
     */
    virtual bool computeMacroStress(FloatArray &answer, const FloatArray &input, TimeStep *tStep);

    /** Computes the macroscopic tangent.
     * @param answer Macroscopic tangent
     * @param tStep Time step to evaluate for
     */
    virtual void computeMacroTangent(FloatMatrix &answer, TimeStep *tStep);

    /** Computes the macroscopic stress for the given input.
     * @param answer Macroscopic stress
     * @param tStep Time step to evaluate for
     */
    virtual void computeMacroStressFromDirichlet(FloatArray &answer, TimeStep *tStep);

    /** Computes the macroscopic tangent for the previously last stress computated stress.
     * May not be called before computeMacroStressFromDirichlet.
     * @param answer The macroscopic stress-strainrate tangent
     * @param tStep Time step to evaluate for
     */
    virtual void computeMacroStressTangentFromDirichlet(FloatMatrix &answer, TimeStep *tStep);

    /** Activates type of homogenization.
     * @return True if mode is supported
     */
    virtual bool activateHomogenizationMode(HomogenizationType ht);

    /// Updates coefficient vectors for new nodal positions.
    void updateC();

    const char *giveClassName() const { return "StokesFlowStressRVE"; }
    classType giveClassID() const { return StokesFlowStressHomogenizationClass; }
};
} // end namespace oofem

#endif // stokesflowstressrve_h


