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

#ifndef calmls_h
#define calmls_h

#include <set>
#include <vector>
#include <memory>

#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "floatarray.h"
#include "intarray.h"
#include "dofiditem.h"

///@name Input fields for CylindricalALM
//@{
#define _IFT_CylindricalALM_Name "calm"
#define _IFT_CylindricalALM_psi "psi"
#define _IFT_CylindricalALM_maxiter "maxiter"
#define _IFT_CylindricalALM_maxrestarts "maxrestarts"
#define _IFT_CylindricalALM_minsteplength "minsteplength"
#define _IFT_CylindricalALM_steplength "steplength"
#define _IFT_CylindricalALM_initialsteplength "initialsteplength"
#define _IFT_CylindricalALM_forcedinitialsteplength "forcedinitialsteplength"
#define _IFT_CylindricalALM_reqiterations "reqiterations"
#define _IFT_CylindricalALM_miniterations "miniter"
#define _IFT_CylindricalALM_manrmsteps "manrmsteps"
#define _IFT_CylindricalALM_hpcmode "hpcmode"
#define _IFT_CylindricalALM_hpc "hpc"
#define _IFT_CylindricalALM_hpcw "hpcw"
#define _IFT_CylindricalALM_lstype "lstype"
#define _IFT_CylindricalALM_linesearch "linesearch"
#define _IFT_CylindricalALM_lsearchtol "lsearchtol"
#define _IFT_CylindricalALM_lsearchamp "lsearchamp"
#define _IFT_CylindricalALM_lsearchmaxeta "lsearchmaxeta"
#define _IFT_CylindricalALM_nccdg "nccdg"
#define _IFT_CylindricalALM_ccdg "ccdg"
#define _IFT_CylindricalALM_rtolv "rtolv"
#define _IFT_CylindricalALM_rtolf "rtolf"
#define _IFT_CylindricalALM_rtold "rtold"
#define _IFT_CylindricalALM_rootselectiontype "rootselectiontype"
//@}


namespace oofem {
class Domain;
class EngngModel;
class ParallelContext;

#define calm_SMALL_NUM 1.e-20
#define calm_SMALL_ERROR_NUM 1.e-6

/**
 * Implementation of sparse nonlinear solver with indirect control.
 * It uses Cylindrical Arc-Length Method algorithm.
 *
 * DESCRIPTION :
 * Perform solution of non-linear problem in the form
 * Kt * deltaR = g
 * where g is defined as g = g(Lambda,X) = Lambda* R - F(X).
 * During iteration process, we want g became zero vector.
 *
 * This method uses Modified Newton Raphson iteration scheme
 *
 * If we solve non-linear static we can interpret symbols as follows:
 *
 * Kt     - tangential stiffness
 * dX     - increment of displacements
 * g      - vector of unbalanced forces (at the end should be zero one)
 * Lambda - Load level
 * Rt     - Quasi Total Load vector Rt = R0 + Lambda*R
 * X      - total displacement vector
 * F(X)   - Nodal representation of (real) internal forces.
 * Psi    - control parameter (=0. means displacement control)
 * Bergan_k0 - variable used to compute bergan's parameters of current stiffness.
 * calm_NR_Mode - variable controlling the mode of NRM (ModifiedNR, Full NRM (stiffness update after each iteration),
 *              Modified Accelerated NRM (we perform iteration with stiffness matrix updated only after calm_MANRMSteps)
 * calm_NR_OldMode - variable containing the old mode of NRM, which will be restored after
 *              calm_NR_ModeTick iterations.
 * calm_NR_ModeTick - see calm_NR_OldMode.
 * calm_MANRMSteps - if calm_NR_Mode == calm_accelNRM, it specifies, that new updated
 *                 stiffness matrix is assembled after calm_MANRMSteps.
 * calm_Control - variable indicating the ALM control.
 * calm_HPCIndirectDofMask - Mask, telling which dofs are used for HPC.
 * calm_HPCWeights - dofs weights in constrain.
 *
 * Tasks:
 * - solving problem: solveYourselfAt.
 * - returning results (increment of displacement, reached level of loading and so on)
 *
 * Variable description  :
 * - K(N,N) = Stiffness matrix (assumed positive definite)
 * - dX(N) = Iterative increment of displacement
 * - Rt = "quasi" total load vector
 * - R0 = Fixed (Initial) load vector
 * - deltaL = Step length
 * - deltaLambda = Increment of load level
 * - Psi = Determines loading or displacement control
 * - dX = Current total increment
 * - F = Nodal representation of (real) internal forces
 * - DeltaLambda = Current total increment of load level
 * - Lambda = Total load level
 * - RTOL= Convergence tolerance
 * - NSMAX = Maximum number of sweeps allowed
 *
 * OUTPUT : (after call solveYourselfAt)
 * - dX = Reached displacement increment
 * - Lambda = Total load level
 * - nite = Number of iterations required to fulfill balance
 * - status = NM_status with flags set to reached state
 */
class OOFEM_EXPORT CylindricalALM : public SparseNonLinearSystemNM
{
protected:
    /// CALM mode type; determines the calm step length control.
    enum calm_ControlType {
        calm_hpc_off = 0, ///< Full ALM with quadratic constrain and all dofs.
        calm_hpc_on, ///< Full ALM with quadratic constrain, taking into account only selected dofs.
        calml_hpc, ///<  Linearized ALM (only displacements), taking into account only selected dofs with given weight.
    };

    /// Controlling mode of Newton-Raphson Method.
    enum calm_NR_ModeType {
        calm_modifiedNRM, ///< Keeps the old tangent.
        calm_fullNRM, ///< Updates the tangent every iteration.
        calm_accelNRM, ///< Updates the tangent after a few steps.
    };

    typedef std :: set< DofIDItem > __DofIDSet;

    int nsmax;
    int maxRestarts;
    double Psi;
    double deltaL, minStepLength, maxStepLength;
    int solved, numberOfRequiredIterations;
    calm_NR_ModeType calm_NR_Mode, calm_NR_OldMode;
    int calm_NR_ModeTick;
    int calm_MANRMSteps;

    /// Minimum hard number of iteration.s
    int minIterations;

    /// Variables for HyperPlaneControl.
    int calm_hpc_init;
    calm_ControlType calm_Control;
    FloatArray calm_HPCWeights;
    /// Array containing equation numbers of dofs under indirect control.
    IntArray calm_HPCIndirectDofMask;
    /// Input array containing dofmanagers and corresponding dof numbers under indirect control.
    IntArray calm_HPCDmanDofSrcArray;
    /// Input array of dofman weights (for hpcmode 2).
    FloatArray calm_HPCDmanWeightSrcArray;

    /// Linear system solver.
    std :: unique_ptr< SparseLinearSystemNM > linSolver;
    /// linear system solver ID.
    LinSystSolverType solverType;

    /// Line search flag.
    int lsFlag;
    /// Line search tolerance.
    double ls_tolerance;
    /// Line search amplification factor.
    double amplifFactor;
    /// Line search parameters (limits).
    double maxEta, minEta;

    // Support for evaluation of error norms for user defined dof-groups.
    /// Number of convergence criteria dof groups.
    int nccdg;
    /// Convergence criteria dof groups.
    std :: vector< __DofIDSet >ccDofGroups;
    /// Relative unbalanced force tolerance for each group.
    FloatArray rtolf;
    /// Relative iterative displacement change tolerance for each group.
    FloatArray rtold;
    /// Parallel context for computing norms, dot products and such.
    ParallelContext *parallel_context;

    /** Type characterizing how the roots of the quadratic equation 
     *  for the arc-length parameter are selected
     *  0 -- cosinus based criterion
     *  1 -- dot product based criterion
     */
    enum RootSelectionType {
      RST_Cos=0,
      RST_Dot=1,
    };

    /// Root selection type
    RootSelectionType rootselectiontype;
    /// previous increment of dX, needed by root selection type 1
    FloatArray old_dX;

public:
    CylindricalALM(Domain * d, EngngModel * m);
    virtual ~CylindricalALM();

    // Overloaded methods:
    NM_Status solve(SparseMtrx &K, FloatArray &R, FloatArray *R0,
                    FloatArray &X, FloatArray &dX, FloatArray &F,
                    const FloatArray &internalForcesEBENorm, double &ReachedLambda, referenceLoadInputModeType rlm,
                    int &nite, TimeStep *) override;
    double giveCurrentStepLength() override { return deltaL; }
    void setStepLength(double s) override { deltaL = s; }
    void initializeFrom(InputRecord &ir) override;
    bool referenceLoad() const override { return true; }
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
    void setDomain(Domain *d) override {
        this->domain = d;
        if ( linSolver ) {
            linSolver->setDomain(d);
        }
    }
    void reinitialize() override {
        calm_hpc_init = 1;
        if ( linSolver ) {
            linSolver->reinitialize();
        }
    }
    const char *giveClassName() const override { return "CylindricalALM"; }
    const char *giveInputRecordName() const { return _IFT_CylindricalALM_Name; }

    SparseLinearSystemNM *giveLinearSolver() override;

protected:
    void convertHPCMap();

    /**
     * @return If 0 then ok, 1 then failure (restart).
     */
    int computeDeltaLambda(double &deltaLambda, const FloatArray &dX, const FloatArray &deltaXt,
                           const FloatArray &deltaX_, const FloatArray &R, double RR, double eta,
                           double deltaL, double DeltaLambda0, int neq);

    void search(int istep, FloatArray &prod, FloatArray &eta, double amp,
                double maxeta, double mineta, int &status);

    /// Evaluates the convergence criteria.
    bool checkConvergence(const FloatArray &R, const FloatArray *R0, const FloatArray &F,
                          const FloatArray &X, const FloatArray &ddX,
                          double Lambda, double RR0, double RR, double drProduct,
                          const FloatArray &internalForcesEBENorm, int nite, bool &errorOutOfRange);

    /// Perform line search optimization of step length
    void do_lineSearch(FloatArray &X, const FloatArray &XInitial, const FloatArray &deltaX_, const FloatArray &deltaXt,
                       const FloatArray &dXm1, FloatArray &dX, FloatArray &ddX,
                       const FloatArray &R, const FloatArray *R0, const FloatArray &F,
                       double &DeltaLambda, double &DeltaLambdam1, double &deltaLambda,
                       double &Lambda, double &ReachedLambda, double RR, double &drProduct, TimeStep *tStep);
};
} // end namespace oofem
#endif // calmls_h
