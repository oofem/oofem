/* $Header: /home/cvs/bp/oofem/oofemlib/src/calmls.h,v 1.1.4.1 2004/04/13 11:28:15 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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
#endif
#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "sparsemtrx.h"
#include "flotarry.h"


class Domain;
class EngngModel;

#define calm_SMALL_NUM 1.e-20
#define calm_SMALL_ERROR_NUM 1.e-6

/**
 * Implementation of sparse nonlinear solver with indirect controll.
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
     * During iteration proces, we want g became zero vector.
     *
     * =======>   This method uses Modified Newton Raphson iteration scheme  <======
     *
     * If we solve non-linear static we can interprete symbols as follows:
     *
     * Kt     - tangential stiffness
     * deltaR - increment of displacements
     * g      - vector of unbalanced forces (at the end should be zero one)
     * Lambda - Load level
     * Rt     - Quasi Total Load vector Rt = R0 + Lambda*R
     * r      - total displacement vector
     * F(r)   - Nodal representation of (real) internal forces.
     * Psi    - control parametr (=0. means displacement control)
     * Bergan_k0 - variable used to compute bergan's parameters of current stiffness.
     * calm_NR_Mode - variable controlling the mode of NRM (ModifiedNR, Full NRM (stifnees update after each iteration),
     *              Modified Accelerated NRM (we perform iteration with stiffness matrix updated only after calm_MANRMSteps)
     * calm_NR_OldMode - variable containing the old mode of NRM, which will be restored after
     *              calm_NR_ModeTick iterations.
     * calm_NR_ModeTick - see calm_NR_OldMode.
     * calm_MANRMSteps - if calm_NR_Mode == calm_accelNRM, it specifies, that new updated
     *                 stiffness matrix is assembled after calm_MANRMSteps.
     * calm_Controll - variable indicating the ALM controll.
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
    /**
     * CALM mode type; determines the calm step length control.
     * calm_hpc_off - full ALM with quadratic constrain and all dofs
     * calm_hpc_on  - full ALM with quadratic constrain, taking into account only selected dofs
     * calml_hpc - linearized ALM (only displacements), taking into account only selected dofs with given weight
     */
    enum    calm_ControllType { calm_hpc_off = 0, calm_hpc_on, calml_hpc };
    enum    calm_NR_ModeType { calm_modifiedNRM, calm_fullNRM, calm_accelNRM };

    //FloatArray     *F;
    int nsmax;
    double Psi;
    double deltaL, minStepLength, maxStepLength;
    int solved, numberOfRequiredIterations;
    calm_NR_ModeType calm_NR_Mode, calm_NR_OldMode;
    int calm_NR_ModeTick;
    int calm_MANRMSteps;

    // variables for HyperPlaneControll
    int calm_hpc_init;
    calm_ControllType calm_Controll;
    FloatArray calm_HPCWeights;
    // array containing equation numbers of dofs under indirect controll
    IntArray calm_HPCIndirectDofMask;
    // input array containing dofmanagers and corresponding dof numbers under indirect controll
    IntArray calm_HPCDmanDofSrcArray;
    // input arry of dofman weights (for hpcmode 2)
    FloatArray calm_HPCDmanWeightSrcArray;

    // linear system solver
    SparseLinearSystemNM *linSolver;
    // linear system solver ID
    LinSystSolverType solverType;

    /// lineserach flag
    int lsFlag;
    /// line search tolerance
    double ls_tolerance;
    /// line serch aplification factor
    double amplifFactor;
    /// line search parameters (limits)
    double maxEta, minEta;

public:
    CylindricalALM(int i, Domain *d, EngngModel *m, EquationID ut);
    // constructor
    ~CylindricalALM();              // destructor


    /**
     * Solves the given sparse linear system of equations g(x,l)=l-F(x); dx=K^{-1}g+ dl K^{-1}R.
     * Total load vector not passed, it is defined as l*R+R0, where l is scale factor
     * @param K coefficient matrix (K = dF/dx; stiffness matrix)
     * @param R  incremental Rhs (incremental load)
     * @param R0 initial Rhs (initial load)
     * @param Rr linearization of K*rri, where rri is increment of prescribed displacements
     * @param r  total solution (total displacement)
     * @param dr increment of solution (incremental displacaments)
     * @param l  Rhs scale factor (load level)
     * @param rtol prescribed tolerance (g residual and iterative r change;)
     * @param F  InternalRhs (real internal forces)
     * @param rlm - reference load mode
     * @return NM_Status value
     */
    virtual NM_Status solve(SparseMtrx *k, FloatArray *Ri, FloatArray *R0,
                            FloatArray *Rr, FloatArray *r, FloatArray *DeltaR, FloatArray *F,
                            double &ReachedLambda, double rtol, referenceLoadInputModeType rlm,
                            int &nite, TimeStep *);

    virtual double giveCurrentStepLength() { return deltaL; }
    virtual void   setStepLength(double l) { deltaL = l; }

    // management  components
    IRResultType initializeFrom(InputRecord *ir);

    /** Stores receiver state to output stream.
     *  Receiver should write class-id first in order to allow test
     *  whether correct data are then restored.
     *  @param stream output stream
     *  @param mode determines ammount of info in stream (state, definition,...)
     *  @param obj special parameter, used only to send particular integration
     *  point to material class version of this method. Except this
     *  case, obj parameter is always NULL pointer.*/
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /** Restores the receiver state previously written in stream.
     *  @see saveContext member function.*/
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    // identification
    const char *giveClassName() const { return "CylindricalALM"; }
    classType giveClassID() const { return CylindricalALMSolverClass; }
    /// sets associated Domain
    virtual void         setDomain(Domain *d) { this->domain = d;
                                                if ( linSolver ) { linSolver->setDomain(d); } }
    /// This method clears receiver cached data dependent on topology, when it changes.
    virtual void reinitialize() { calm_hpc_init = 1;
                                  if ( linSolver ) { linSolver->reinitialize(); } }
protected:
    void convertHPCMap();
    SparseLinearSystemNM *giveLinearSolver();
    int  computeDeltaLambda(double &deltaLambda, FloatArray &DeltaR, FloatArray &deltaRt,
                            FloatArray &deltaR_, FloatArray &R, double RR, double eta,
                            double deltaL, double DeltaLambda0, int neq);

    void search(int istep, FloatArray &prod, FloatArray &eta, double amp,
                double maxeta, double mineta, int &status);
};

#endif // calmls_h









