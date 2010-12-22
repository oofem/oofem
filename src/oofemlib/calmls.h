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
 *               Copyright (C) 1993 - 2010   Borek Patzak
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

//   ********************************************
//   *** CLASS CYLINDRICAL ARC LENGTH METHOD  ***
//   ********************************************


#ifndef calmls_h
#define calmls_h

#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <set>
 #include <vector>
#endif
#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "sparsemtrx.h"
#include "flotarry.h"

namespace oofem {
class Domain;
class EngngModel;

#define calm_SMALL_NUM 1.e-20
#define calm_SMALL_ERROR_NUM 1.e-6

/**
 * Implementation of sparse nonlinear solver with indirect control.
 * It uses Cylindrical Arc-Length Method algorithm.
 */
class CylindricalALM : public SparseNonLinearSystemNM
{
    /*
     * This class implements the class NumericalMethod instance Cylindrical Arc-Length Method
     * for solving non-linear problems.
     *
     * DESCRIPTION :
     * Perform solution of non-linear problem in the form
     * Kt * deltaR = g
     * where g is defined as g = g(Lambda,r) = Lambda* R - F(r).
     * During iteration process, we want g became zero vector.
     *
     * =======>   This method uses Modified Newton Raphson iteration scheme  <======
     *
     * If we solve non-linear static we can interpret symbols as follows:
     *
     * Kt     - tangential stiffness
     * deltaR - increment of displacements
     * g      - vector of unbalanced forces (at the end should be zero one)
     * Lambda - Load level
     * Rt     - Quasi Total Load vector Rt = R0 + Lambda*R
     * r      - total displacement vector
     * F(r)   - Nodal representation of (real) internal forces.
     * Psi    - control parameter (=0. means displacement control)
     * Bergan_k0 - variable used to compute bergan's parameters of current stiffness.
     * calm_NR_Mode - variable controlling the mode of NRM (ModifiedNR, Full NRM (stiffness update after each iteration),
     *              Modified Accelerated NRM (we perform iteration with stiffness matrix updated only after calm_MANRMSteps)
     * calm_NR_OldMode - variable containing the old mode of NRM, which will be restored after
     *              calm_NR_ModeTick iterations.
     * calm_NR_ModeTick - see calm_NR_OldMode.
     * calm_MANRMSteps - if calm_NR_Mode == calm_accelNRM, it specifies, that new updated
     *                 stiffness matrix is assembled after calm_MANRMSteps.
     * calm_Controll - variable indicating the ALM control.
     * calm_HPCIndirectDofMask - Mask, telling which dofs are used for HPC.
     * calm_HPCWeights - dofs weights in constrain.
     * TASKS :
     *
     * - solving problem
     *   solveYourselfAt.
     * - returning results (increment of displacement,
     *                      reached level of loading and so on)
     *
     * Variable description  :
     *
     *       K(N,N)    = STIFFNESS MATRIX (ASSUMED POZITIVE DEFINITE)        *
     *       deltaR(N) = ITERATIVE INCREMENT OF DISPLACEMENT                 *
     *       Rt        = "quasi" total LOAD VECTOR                           *
     *   R0        = Fixed (Initial) load vector                         *
     * deltaL    = STEP LENGTH                                         *
     * deltaLambda= INCREMENT OF LOAD LEVEL                            *
     * Psi       = DETERMINES LOADING OR DISPLACEMENT CONTROL          *
     * DeltaR    = CURRENT TOTAL INCREMENT                             *
     *       F         = NODAL REPRESENTATION OF (REAL) INTERNAL FORCES      *
     * DeltaLambda= CURRENT TOTAL INCEREMENT OF LOAD LEVEL             *
     * Lambda    = TOTAL LOAD LEVEL
     *
     *       RTOL      = CONVERGENCE TOLERANCE                               *
     *       NSMAX     = MAXIMUM NUMBER OF SWEEPS ALLOVED                    *
     *
     * OUTPUT : (after call solveYourselfAt)
     *       K(N,N)    = DIAGONALIZED STIFFNESS MATRIX                       *
     *       DeltaR    = REACHED DISPLACEMENT INCREMENT                      *
     * Lambda    = TOTAL LOAD LEVEL                                    *
     * nite      = NUMBER OF ITERATIONS REQUIRED TO FULLFIL BALANCE    *
     * status    = NM_status with flags set to reached state (see cltypes.h) *
     *
     */
protected:
    /// CALM mode type; determines the calm step length control.
    enum calm_ControllType {
        calm_hpc_off = 0, /// Full ALM with quadratic constrain and all dofs.
        calm_hpc_on, /// Full ALM with quadratic constrain, taking into account only selected dofs.
        calml_hpc, ///  Linearized ALM (only displacements), taking into account only selected dofs with given weight.
    };

    // TODO: Document me
    enum calm_NR_ModeType {
        calm_modifiedNRM,
        calm_fullNRM,
        calm_accelNRM,
    };

    typedef std :: set< DofID >__DofIDSet;

    int nsmax;
    double Psi;
    double deltaL, minStepLength, maxStepLength;
    int solved, numberOfRequiredIterations;
    calm_NR_ModeType calm_NR_Mode, calm_NR_OldMode;
    int calm_NR_ModeTick;
    int calm_MANRMSteps;

    /// Minimum hard number of iteration.s
    int minIterations;

    /// Variables for HyperPlaneControll.
    int calm_hpc_init;
    calm_ControllType calm_Controll;
    FloatArray calm_HPCWeights;
    /// Array containing equation numbers of dofs under indirect control.
    IntArray calm_HPCIndirectDofMask;
    /// Input array containing dofmanagers and corresponding dof numbers under indirect control.
    IntArray calm_HPCDmanDofSrcArray;
    /// Input array of dofman weights (for hpcmode 2).
    FloatArray calm_HPCDmanWeightSrcArray;

    /// Linear system solver.
    SparseLinearSystemNM *linSolver;
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


public:
    CylindricalALM(int i, Domain *d, EngngModel *m, EquationID ut);
    ~CylindricalALM();

    // Overloaded methods:
    virtual NM_Status solve(SparseMtrx *k, FloatArray *Ri, FloatArray *R0,
                            FloatArray *Rr, FloatArray *r, FloatArray *DeltaR, FloatArray *F,
                            double &internalForcesEBENorm, double &ReachedLambda, referenceLoadInputModeType rlm,
                            int &nite, TimeStep *);
    virtual double giveCurrentStepLength() { return deltaL; }
    virtual void setStepLength(double l) { deltaL = l; }
    IRResultType initializeFrom(InputRecord *ir);
    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual void setDomain(Domain *d) {
        this->domain = d;
        if ( linSolver ) {
            linSolver->setDomain(d);
        }
    }
    virtual void reinitialize() {
        calm_hpc_init = 1;
        if ( linSolver ) {
            linSolver->reinitialize();
        }
    }
    const char *giveClassName() const { return "CylindricalALM"; }
    classType giveClassID() const { return CylindricalALMSolverClass; }

protected:
    void convertHPCMap();
    SparseLinearSystemNM *giveLinearSolver();
    int  computeDeltaLambda(double &deltaLambda, FloatArray &DeltaR, FloatArray &deltaRt,
                            FloatArray &deltaR_, FloatArray &R, double RR, double eta,
                            double deltaL, double DeltaLambda0, int neq);

    void search(int istep, FloatArray &prod, FloatArray &eta, double amp,
                double maxeta, double mineta, int &status);

    /// Evaluates the convergence criteria.
    bool checkConvergence(FloatArray &R, FloatArray *R0, FloatArray &F,
                          FloatArray &r, FloatArray &rIterIncr,
                          double Lambda, double RR0, double RR, double drProduct,
                          double internalForcesEBENorm, int nite, bool &errorOutOfRange);

    /// Perform line search optimization of step length
    void do_lineSearch(FloatArray &r, FloatArray &rInitial, FloatArray &deltaR_, FloatArray &deltaRt,
                       FloatArray &DeltaRm1, FloatArray &DeltaR, FloatArray &deltaR,
                       FloatArray &R, FloatArray *R0, FloatArray &F,
                       double &DeltaLambda, double &DeltaLambdam1, double &deltaLambda,
                       double &Lambda, double &ReachedLambda, double RR, double &drProduct, TimeStep *tNow);
};
} // end namespace oofem
#endif // calmls_h









