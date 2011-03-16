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


//   *************************************
//   *** CLASS SparseNonLinearSystemNM ***
//   *************************************


#ifndef sparsenonlinsystemnm_h
#define sparsenonlinsystemnm_h


#include "nummet.h"

#include "equationid.h"
#include "nmstatus.h"

#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

namespace oofem {
class EngngModel;
class SparseMtrx;
class FloatArray;


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
class SparseNonLinearSystemNM : public NumericalMethod
{
public:
    /**
     * The following parameter allows to specify how the reference load vector
     * is obtained from given totalLoadVector and initialLoadVector.
     * The initialLoadVector describes the part of loading which does not scale.
     */
    enum referenceLoadInputModeType {
        rlm_total=0, /// The reference incremental load vector is defined as totalLoadVector assembled at given time.
        rlm_inceremental=1, /// The reference load vector is obtained as incremental load vector at given time.
    };

protected:
    /// Equation to solve for.
    EquationID ut;

    /// Load level
    double deltaL;

public:
    /// Constructor
    SparseNonLinearSystemNM(int i, Domain *d, EngngModel *m, EquationID ut) : NumericalMethod(i, d, m) { this->ut = ut; }
    /// Destructor
    ~SparseNonLinearSystemNM() { }

    /**
     * Solves the given sparse linear system of equations @f$s  R + R_0 - F(r) = 0@f$; @f$ dr=K^{-1}(s R + R_0 - F({}^n r)) @f$.
     * Total load vector not passed, it is defined as @f$ s R + R_0 @f$, where s is scale factor.
     * @see EngngModel::updateComponent Which is used to update the stiffness matrix and load vector.
     * @param K  Coefficient matrix (@f$\displaystyle K = \frac{\partial F}{\partial r} @f$; stiffness matrix).
     * @param R  Reference incremental RHS (incremental load).
     * @param R0 Initial RHS (initial load).
     * @param Rr Linearization of K*rri, where rri is increment of prescribed displacements.
     * @param r  Total solution (total displacement).
     * @param dr Increment of solution (incremental displacements).
     * @param s  RHS scale factor (load level).
     * @param F  InternalRhs (real internal forces).
     * @param nite Number of iterations needed.
     * @param rlm Reference load mode.
     * @param internalForcesEBENorm Norm of internal nodal forces (evaluated on element by element basis).
     * @param tStep Time step to solve for.
     * @return NM_Status value.
     */
    virtual NM_Status solve(SparseMtrx *K, FloatArray *R, FloatArray *R0,
                            FloatArray *Rr, FloatArray *r, FloatArray *dr, FloatArray *F,
                            double &internalForcesEBENorm, double &s, referenceLoadInputModeType rlm,
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

    // Overloaded methods:
    const char *giveClassName() const { return "SparseNonLinearSystemNM"; }
    classType giveClassID() const { return SparseNonLinearSystemNMClass; }

public:
};
} // end namespace oofem
#endif // sparsenonlinsystemnm_h
