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


//   *************************************
//   *** CLASS SparseNonLinearSystemNM ***
//   *************************************


#ifndef sparsenonlinsystemnm_h
#define sparsenonlinsystemnm_h

#include "nummet.h"
#include "nmstatus.h"
#include "floatarray.h"
#include "intarray.h"

///@name Input fields for SparseNonLinearSystemNM
//@{
#define _IFT_NonLinearStatic_randPertAmplitude "rpa"
#define _IFT_NonLinearStatic_randSeed "rseed"
#define _IFT_NonLinearStatic_pert "pert"
#define _IFT_NonLinearStatic_pertw "pertw"
//@}

namespace oofem {
class EngngModel;
class SparseMtrx;
class FloatArray;
class TimeStep;
class SparseLinearSystemNM;

/**
 * This base class is an abstraction for all numerical methods solving sparse
 * nonlinear system of equations. The purpose of this class is to declare
 * the general interface to all numerical methods solving this kind of
 * problem. This interface allows to use any suitable
 * instance of the Numerical method class to the solve problem,
 * and leave the  whole engineering model code,
 * including mapping, unchanged, because all instances of this class
 * provide the common interface.
 */
class OOFEM_EXPORT SparseNonLinearSystemNM : public NumericalMethod
{
protected:
    /// Load level
    double deltaL;

    /// Amplitude of a random perturbation applied on the solution before the iteration process
    double randPertAmplitude;
    int randSeed;
    bool pert_init_needed;
    IntArray igp_PertDmanDofSrcArray;
    FloatArray igp_PertWeightArray;
    IntArray igp_Map;
    FloatArray igp_Weight;
public:
    /**
     * The following parameter allows to specify how the reference load vector
     * is obtained from given totalLoadVector and initialLoadVector.
     * The initialLoadVector describes the part of loading which does not scale.
     */
    enum referenceLoadInputModeType {
        rlm_total=0, ///< The reference incremental load vector is defined as totalLoadVector assembled at given time.
        rlm_incremental=1, ///< The reference load vector is obtained as incremental load vector at given time.
    };    

    /// Constructor
    SparseNonLinearSystemNM(Domain * d, EngngModel * m) : NumericalMethod(d, m), igp_PertDmanDofSrcArray(), igp_PertWeightArray(), igp_Map(), igp_Weight() { pert_init_needed = false; }
    /// Destructor
    virtual ~SparseNonLinearSystemNM() { }

    /**
     * Solves the given sparse linear system of equations @f$ s  R + R_0 - F(X) = 0 @f$.
     * Total load vector is not passed, it is defined as @f$ s R + R_0 @f$, where s is scale factor.
     * @see EngngModel::updateComponent Used to update the stiffness matrix and load vector.
     * @param K  Coefficient matrix (@f$\displaystyle K = \frac{\partial F}{\partial X} @f$; stiffness matrix).
     * @param R  Reference incremental RHS (incremental load).
     * @param R0 Initial RHS (initial load).
     * @param X  Total solution (total displacement).
     * @param dX Increment of solution (incremental displacements).
     * @param F InternalRhs (real internal forces).
     * @param internalForcesEBENorm Norm of internal nodal forces (evaluated on element by element basis) (split into each DOF id).
     * @param s RHS scale factor (load level).
     * @param rlm Reference load mode.
     * @param nite Number of iterations needed.
     * @param tStep Time step to solve for.
     * @return NM_Status value.
     */
    virtual NM_Status solve(SparseMtrx &K, FloatArray &R, FloatArray *R0, FloatArray *iR,
                            FloatArray &X, FloatArray &dX, FloatArray &F,
                            const FloatArray &internalForcesEBENorm, double &s, referenceLoadInputModeType rlm,
                            int &nite, TimeStep *tStep) = 0;

    /**
     * Returns step length.
     * @see solve For more details on the step length s.
     * @return Current step length
     */
    virtual double giveCurrentStepLength() { return this->deltaL; }
    /**
     * Sets the step length.
     * @see solve For more details on the step length s.
     * @param s New step length.
     */
    virtual void setStepLength(double s) { this->deltaL = s; }
    /**
     * Prints status message of receiver to output stream.
     * Prints the message corresponding to last solve.
     * @param outputStream Stream to print state to.
     */
    virtual void printState(FILE *outputStream) { }

    /**
     * Constructs (if necessary) and returns a linear solver.
     * Public method because some problems require it for sensitivity analysis, etc. even for nonlinear problems (e.g. tangent relations in multiscale simulations).
     */
    virtual SparseLinearSystemNM *giveLinearSolver() { return NULL; }

    IRResultType initializeFrom(InputRecord *ir);
    virtual void convertPertMap();
    virtual void applyPerturbation(FloatArray* displacement);

    virtual const char *giveClassName() const { return "SparseNonLinearSystemNM"; }
    /// Error printing helper.
    std :: string errorInfo(const char *func) const { return std :: string(giveClassName()) + "::" + func; }
};
} // end namespace oofem
#endif // sparsenonlinsystemnm_h
