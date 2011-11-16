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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef stokesflowstressrve_h
#define stokesflowstressrve_h

#include "inputrecord.h"
#include "stokesflow.h"

#include "sparsemtrxtype.h"
#include "sparselinsystemnm.h"
#include "linsystsolvertype.h"

namespace oofem {
class PrescribedGradient;

/**
 * Stokes flow model with homogenization for stress-strainrate.
 * @author Mikael Öhman
 */
class StokesFlowStressHomogenization : public StokesFlow
{
protected:
    /// Non-free sub-matrices of the total stiffness matrix.
    SparseMtrx *K_cf, *K_cc;
    /// Coefficient matrix for values of prescribed dofs.
    FloatMatrix C;
    /// Solver type for homogenization.
    LinSystSolverType solverType;
    /// Numerical method for homogenization.
    SparseLinearSystemNM *linNumericalMethod;

public:
    StokesFlowStressHomogenization(int i, EngngModel *_master = NULL);
    virtual ~StokesFlowStressHomogenization();

    virtual void updateYourself(TimeStep *tStep);

    PrescribedGradient *giveDirichletBC();

    /**
     * Computes the volume of the whole microscopic domain (pores included).
     * @param di Domain index.
     * @return RVE size (volume or area).
     */
    virtual double computeSize(int di);
    /**
     * Computes the macroscopic stress.
     * @param answer Macroscopic stress.
     * @param input Strain rate.
     * @param tStep Time step to evaluate at.
     * @return True if successful.
     */
    virtual bool computeMacroStress(FloatArray &answer, const FloatArray &input, TimeStep *tStep);

    /** Computes the macroscopic tangent.
     * @param answer Macroscopic tangent.
     * @param tStep Time step to evaluate at.
     */
    virtual void computeMacroTangent(FloatMatrix &answer, TimeStep *tStep);

    /**
     * Computes the macroscopic stress for the given input.
     * @param answer Macroscopic stress.
     * @param tStep Time step to evaluate at.
     */
    virtual void computeMacroStressFromDirichlet(FloatArray &answer, TimeStep *tStep);

    /**
     * Computes the macroscopic tangent for the previously last stress computed stress.
     * May not be called before computeMacroStressFromDirichlet.
     * @param answer The macroscopic stress-strain rate tangent.
     * @param tStep Time step to evaluate at.
     */
    virtual void computeMacroStressTangentFromDirichlet(FloatMatrix &answer, TimeStep *tStep);

    virtual int forceEquationNumbering(int di);

    virtual SparseLinearSystemNM *giveLinearNumericalMethod();

    const char *giveClassName() const { return "StokesFlowStressStressHomogenization"; }
    classType giveClassID() const { return StokesFlowStressHomogenizationClass; }
};
} // end namespace oofem

#endif // stokesflowstressrve_h


