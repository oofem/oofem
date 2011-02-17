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

#ifndef stokesflow_h
#define stokesflow_h

#include "engngm.h"
#include "sparsemtrxtype.h"
#include "sparsenonlinsystemnm.h"
//#include "linsystsolvertype.h"

namespace oofem {
/**
 * Implements the engineering model to solve incompressible Stokes flow.
 * Stokes flow means acceleration is ignored.
 * @author Carl Sandström
 * @author Mikael Öhman
 */
class StokesFlow : public EngngModel
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

    SparseMtrx *stiffnessMatrix;
    FloatArray internalForces;
    FloatArray externalForces;
    FloatArray incrementOfSolution;

    // Related to homogenization
    /** Boolean value to keep track of advancing the primary field.
     * This is necessary when solveYourselfAt might be called several times per timestep. Shouldn't need to do this.
     */
    bool hasAdvanced;

public:
    StokesFlow(int i, EngngModel *_master = NULL);
    virtual ~StokesFlow();

    /// @see EngngModel::solveYourselfAt
    void solveYourselfAt(TimeStep *tStep);

    /** Updates everything for the problem.
     * Updates the internal state for the elements.
     * Also calls updateNodalPositions.
     */
    virtual void updateYourself(TimeStep *tStep);

    /// @see EngngModel::giveUnknownComponent
    virtual double giveUnknownComponent(EquationID eid, ValueModeType mode, TimeStep *tStep, Domain *domain, Dof *dof);
    virtual double    giveUnknownComponent(UnknownType, ValueModeType, TimeStep *, Domain *, Dof *);

    /** Numbers all equations.
     * Numbers velocities first for each node in order, then pressures.
     * @param id Domain id
     * @return Number of equations for domain id
     */
    virtual int forceEquationNumbering(int id);

    /** Initialization from given input record.
     * Reads
     * - mstype Sparse matrix type (enum, optional, default SMT_PetscMtrx)
     * - deltat Time increment (real, optional, default 1.0)
     * - nodalupdatescheme How to update nodal positions (enum, optional, default US_None)
     */
    virtual IRResultType initializeFrom(InputRecord *ir);

    /// @see EngngModel::checkConsistency
    virtual int checkConsistency();

    /// @see EngngModel::printDofOutput
    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *tStep);

    /** Updates internal state of all gauss points.
     * @param tStep Current time step
     */
    virtual void updateInternalState(TimeStep *tStep);

    /// @see EngngModel::updateComponent
    virtual void updateComponent(TimeStep *tStep, NumericalCmpn cmpn, Domain *d);

    /// @see EngngModel::giveNumericalMethod
    NumericalMethod *giveNumericalMethod(TimeStep *tStep);

    /// @see EngngModel::giveNextStep  (i think there should be a default)
    TimeStep *giveNextStep();

    // Shouldn't need to overload. Should be moved to EngngModel::terminate TODO
    void terminate(TimeStep *tStep) { this->updateYourself(tStep);
                                      this->doStepOutput(tStep);
                                      this->saveStepContext(tStep); }

    const char *giveClassName() const { return "StokesFlow"; }
    classType giveClassID() const { return StokesFlowClass; }
};
} // end namespace oofem

#endif // stokesflow_h


