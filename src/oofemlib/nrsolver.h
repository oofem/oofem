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


#ifndef nrsolver_h
#define nrsolver_h

#include <set>
#include <vector>

#include "sparselinsystemnm.h"
#include "sparsenonlinsystemnm.h"
#include "sparsemtrx.h"
#include "floatarray.h"
#include "linesearch.h"

#ifdef __PETSC_MODULE
 #include <petscksp.h>
#endif

///@name Input fields for NRSolver
//@{
#define _IFT_NRSolver_maxiter "maxiter"
#define _IFT_NRSolver_miniterations "miniter"
#define _IFT_NRSolver_minsteplength "minsteplength"
#define _IFT_NRSolver_manrmsteps "manrmsteps"
#define _IFT_NRSolver_lstype "lstype"
#define _IFT_NRSolver_ddm "ddm"
#define _IFT_NRSolver_ddv "ddv"
#define _IFT_NRSolver_ddltf "ddltf"
#define _IFT_NRSolver_linesearch "linesearch"
#define _IFT_NRSolver_rtolv "rtolv"
#define _IFT_NRSolver_rtolf "rtolf"
#define _IFT_NRSolver_rtold "rtold"
//@}

namespace oofem {
class Domain;
class EngngModel;

/**
 * This class implements Newton-Raphson Method, derived from abstract NumericalMethod class
 * for solving non-linear problems.
 * Except the traditional load control, it also provides direct displacement control without
 * requiring BC applied,
 * so that the equation renumbering is not required when combining arc-length and Newton-Raphson
 * solvers within a simulation.
 *
 * The direct displacement control is achieved by adding a large number alpha to the corresponding
 * diagonal member of K and replacing the right-hand side by alpha*prescribed_value.
 * If alpha is very much larger than other stiffness coefficients then this alteration
 * effectively replaces the corresponding equation by
 * alpha*unknown_value = alpha*prescribed_value
 * that is, the required condition, but the whole system remains symmetric and minimal
 * changes are necessary in the computational sequence.
 * The above artifice has been introduced by Payne and Irons.
 */
class NRSolver : public SparseNonLinearSystemNM
{
private:
    enum nrsolver_ModeType { nrsolverModifiedNRM, nrsolverFullNRM, nrsolverAccelNRM };

    int nite, nsmax, minIterations;
    double minStepLength;
    int solved;
    nrsolver_ModeType NR_Mode, NR_OldMode;
    int NR_ModeTick;
    int MANRMSteps;

    /// linear system solver
    SparseLinearSystemNM *linSolver;
    /// linear system solver ID
    LinSystSolverType solverType;
    /// sparse matrix version, used to control constrains application to stiffness
    SparseMtrx :: SparseMtrxVersionType smConstraintVersion;
    /// number of prescribed displacement
    int numberOfPrescribedDofs;
    /**
     * Flag indicating that some dofs are controlled under displacement control.
     * In parallel mode, numberOfPrescribedDofs is local (related to specific partition)
     * so its nonzero value does not mean that there are no prescribed dofs on
     * other partitions.
     */
    bool prescribedDofsFlag;

    /// Array of pairs identifying prescribed dofs (node, dof)
    IntArray prescribedDofs;
    /// Array of prescribed values
    FloatArray prescribedDofsValues;
    /// Load Time Function of prescribed values
    int prescribedDisplacementLTF;
    /// Array of prescribed equations
    IntArray prescribedEqs;
    /// Flag indicating that prescribedEqs were initialized
    bool prescribedEqsInitFlag;
    /// Computed reactions. They are stored in order to print them in printState method.
    FloatArray lastReactions;
    /// Flag indicating whether to use line-search
    int lsFlag;
    /// Line search solver
    LineSearchNM *linesearchSolver;

#ifdef __PETSC_MODULE
    IS prescribedEgsIS;
    bool prescribedEgsIS_defined;
#endif

    /// Relative unbalanced force tolerance for each group
    FloatArray rtolf;
    /// Relative iterative displacement change tolerance for each group
    FloatArray rtold;

public:
    NRSolver(Domain *d, EngngModel *m);
    virtual ~NRSolver();

    // Overloaded methods:
    virtual NM_Status solve(SparseMtrx *k, FloatArray *R, FloatArray *R0,
                            FloatArray *X, FloatArray *dX, FloatArray *F,
                            const FloatArray &internalForcesEBENorm, double &l, referenceLoadInputModeType rlm,
                            int &nite, TimeStep *);
    virtual void printState(FILE *outputStream);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveClassName() const { return "NRSolver"; }

    virtual void setDomain(Domain *d) {
        this->domain = d;
        if ( linSolver ) {
            linSolver->setDomain(d);
        }
        if ( linesearchSolver ) {
            linesearchSolver->setDomain(d);
        }
    }
    virtual void reinitialize() {
        if ( linSolver ) {
            linSolver->reinitialize();
        }
    }

    virtual SparseLinearSystemNM *giveLinearSolver();

protected:
    /// Constructs and returns a line search solver.
    LineSearchNM *giveLineSearchSolver();
    /// Initiates prescribed equations
    void initPrescribedEqs();
    void applyConstraintsToStiffness(SparseMtrx *k);
    void applyConstraintsToLoadIncrement(int nite, const SparseMtrx *k, FloatArray &R,
                                         referenceLoadInputModeType rlm, TimeStep *atTime);

    /**
     * Determines whether or not the solution has reached convergence.
     * @return True if solution has converged, otherwise false.
     */
    bool checkConvergence(FloatArray &RT, FloatArray &F, FloatArray &rhs, FloatArray &ddX, FloatArray &X,
                          double RRT, const FloatArray &internalForcesEBENorm, int nite, bool &errorOutOfRange, TimeStep *tNow);
};
} // end namespace oofem
#endif // nrsolver_h
