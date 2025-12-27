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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef stokesflow_h
#define stokesflow_h

#include "fluidmodel.h"
#include "sparsemtrxtype.h"
#include "topologydescription.h"
#include "linsystsolvertype.h"
#include "floatmatrix.h"
#include "primaryfield.h"
#include "floatarray.h"

///@name Input fields for Stokes' Flow
//@{
#define _IFT_StokesFlow_Name "stokesflow"
#define _IFT_StokesFlow_deltat "deltat"
//@}

namespace oofem {
class SparseNonLinearSystemNM;
class MeshQualityErrorEstimator;

/**
 * Implements the engineering model to solve incompressible Stokes flow.
 * Stokes flow means acceleration is ignored.
 *
 * @author Carl Sandström
 * @author Mikael Öhman
 */
class StokesFlow : public FluidModel
{
protected:
    /// Time increment read from input record.
    double deltaT = 1.;
    /// Primary unknowns.
    std :: unique_ptr< PrimaryField > velocityPressureField;
    /// Sparse matrix type.
    SparseMtrxType sparseMtrxType;
    /// Numerical method.
    std :: unique_ptr< SparseNonLinearSystemNM > nMethod;
    /// Linear solver type.
    LinSystSolverType solverType = ST_Direct;
    /// Element norm for nonlinear analysis (squared)
    FloatArray eNorm;

    /// Used for determining if a new mesh must be created.
    std :: unique_ptr< MeshQualityErrorEstimator > meshqualityee;
    /// Maximum deformation allowed
    double maxdef = 1.;

    std :: unique_ptr< SparseMtrx > stiffnessMatrix;
    FloatArray solutionVector;
    FloatArray internalForces;

    /// Topology state, most notably used for determining if there is a need to remesh.
    TopologyState ts;

public:
    StokesFlow(int i, EngngModel * _master = nullptr);
    ~StokesFlow();

    void solveYourselfAt(TimeStep *tStep) override;

    /**
     * Updates everything for the problem.
     * Updates the internal state for the elements.
     * Also calls updateNodalPositions.
     */
    void updateYourself(TimeStep *tStep) override;

    double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *domain, Dof *dof) override;
    bool newDofHandling() override { return true; }

    double giveReynoldsNumber() override;

    int forceEquationNumbering(int id) override;

    void initializeFrom(InputRecord &ir) override;

    int checkConsistency() override;
    void doStepOutput(TimeStep *tStep) override;
    void updateInternalState(TimeStep *tStep);
    void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d) override;
    void updateSolution(FloatArray &solutionVector, TimeStep *tStep, Domain *d) override;
    void updateInternalRHS(FloatArray &answer, TimeStep *tStep, Domain *d, FloatArray *eNorm) override;
    void updateMatrix(SparseMtrx &mat, TimeStep *tStep, Domain *d) override;
    NumericalMethod *giveNumericalMethod(MetaStep *mStep) override;
    TimeStep *giveNextStep() override;

    const char *giveClassName() const override { return "StokesFlow"; }
    const char *giveInputRecordName() const { return _IFT_StokesFlow_Name; }
};
} // end namespace oofem

#endif // stokesflow_h
