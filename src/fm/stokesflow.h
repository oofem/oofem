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

#ifndef stokesflow_h
#define stokesflow_h

#include "fluidmodel.h"
#include "sparsemtrxtype.h"
#include "topologydescription.h"
#include "linsystsolvertype.h"
#include "floatmatrix.h"

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
    PrimaryField *velocityPressureField;
    /// Sparse matrix type.
    SparseMtrxType sparseMtrxType;
    /// Numerical method.
    SparseNonLinearSystemNM *nMethod;
    /// Linear solver type.
    LinSystSolverType solverType;
    /// Element norm for nonlinear analysis (squared)
    FloatArray eNorm;

    /// Used for determining if a new mesh must be created.
    MeshQualityErrorEstimator *meshqualityee;
    /// Maximum deformation allowed
    double maxdef;

    SparseMtrx *stiffnessMatrix;
    FloatArray internalForces;
    FloatArray externalForces;
    FloatArray incrementOfSolution;

    /// Topology state, most notably used for determining if there is a need to remesh.
    TopologyState ts;

    /**
     * Boolean value to keep track of advancing the primary field.
     * This is necessary when solveYourselfAt might be called several times per time step. Shouldn't need to do this.
     */
    bool hasAdvanced;

public:
    StokesFlow(int i, EngngModel *_master = NULL);
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

    /**
     * Initialization from given input record.
     * Reads
     * - mstype Sparse matrix type (enum, optional, default SMT_PetscMtrx).
     * - deltat Time increment (real, optional, default 1.0).
     * - nodalupdatescheme How to update nodal positions (enum, optional, default US_None).
     */
    virtual IRResultType initializeFrom(InputRecord *ir);

#ifdef __PETSC_MODULE
    virtual void initPetscContexts();
#endif

    virtual int checkConsistency();
    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep);
    virtual void doStepOutput(TimeStep *tStep);
    void updateInternalState(TimeStep *tStep);
    virtual void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d);
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);
    virtual TimeStep *giveNextStep();

    virtual const char *giveClassName() const { return "StokesFlow"; }
    virtual classType giveClassID() const { return StokesFlowClass; }
};
} // end namespace oofem

#endif // stokesflow_h
