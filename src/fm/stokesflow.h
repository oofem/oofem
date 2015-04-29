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
    double deltaT;
    /// Primary unknowns.
    std :: unique_ptr< PrimaryField > velocityPressureField;
    /// Sparse matrix type.
    SparseMtrxType sparseMtrxType;
    /// Numerical method.
    std :: unique_ptr< SparseNonLinearSystemNM > nMethod;
    /// Linear solver type.
    LinSystSolverType solverType;
    /// Element norm for nonlinear analysis (squared)
    FloatArray eNorm;

    /// Used for determining if a new mesh must be created.
    std :: unique_ptr< MeshQualityErrorEstimator > meshqualityee;
    /// Maximum deformation allowed
    double maxdef;

    std :: unique_ptr< SparseMtrx > stiffnessMatrix;
    FloatArray solutionVector;
    FloatArray internalForces;

    /// Topology state, most notably used for determining if there is a need to remesh.
    TopologyState ts;

public:
    StokesFlow(int i, EngngModel * _master = NULL);
    virtual ~StokesFlow();

    virtual void solveYourselfAt(TimeStep *tStep);

    /**
     * Updates everything for the problem.
     * Updates the internal state for the elements.
     * Also calls updateNodalPositions.
     */
    virtual void updateYourself(TimeStep *tStep);

    virtual double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *domain, Dof *dof);

    virtual double giveReynoldsNumber();

    virtual int forceEquationNumbering(int id);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int checkConsistency();
    virtual void doStepOutput(TimeStep *tStep);
    void updateInternalState(TimeStep *tStep);
    virtual void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d);
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);
    virtual TimeStep *giveNextStep();

    virtual const char *giveClassName() const { return "StokesFlow"; }
    virtual const char *giveInputRecordName() const { return _IFT_StokesFlow_Name; }
};
} // end namespace oofem

#endif // stokesflow_h
