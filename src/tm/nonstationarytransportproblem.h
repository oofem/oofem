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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef nonstationarytransportproblem_h
#define nonstationarytransportproblem_h

#include "structengngmodel.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "primaryfield.h"
#include "dofdistributedprimaryfield.h"
#include "transportmaterial.h"

namespace oofem {
/**
 * This class represents linear nonstationary transport problem.
 */
class NonStationaryTransportProblem : public EngngModel
{
protected:
    /**
     * Contains last time stamp of internal variable update.
     * This update is made via various services
     * (like those for computing real internal forces or updating the internal state).
     */
    StateCounterType internalVarUpdateStamp;

    SparseMtrx *lhs;
    /// Right hand side vector from boundary conditions.
    FloatArray bcRhs;
    /// This field stores solution vector. For fixed size of problem, the PrimaryField is used, for growing/decreasing size, DofDistributedPrimaryField applies.
    PrimaryField *UnknownsField;

    LinSystSolverType solverType;
    SparseMtrxType sparseMtrxType;
    /// Numerical method used to solve the problem.
    SparseLinearSystemNM *nMethod;

    /// Initial time from which the computation runs. Default is zero.
    double initT;
    double deltaT;
    double alpha;

    /// If set then stabilization using lumped capacity will be used.
    int lumpedCapacityStab;

    // if set, the receiver flux field will be exported using FieldManager
    //int exportFieldFlag;

    /// Associated time function for time step increment.
    int dtTimeFunction;

    /// Determines if there are change in the problem size (no application/removal of Dirichlet boundary conditions).
    bool changingProblemSize;

public:
    /// Constructor.
    NonStationaryTransportProblem(int i, EngngModel *_master);
    /// Destructor.
    virtual ~NonStationaryTransportProblem();

    virtual void solveYourselfAt(TimeStep *tStep);
    virtual void updateYourself(TimeStep *tStep);
    virtual double giveUnknownComponent(EquationID, ValueModeType, TimeStep *, Domain *, Dof *);
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual void updateDomainLinks();

    virtual TimeStep *giveNextStep();
    virtual TimeStep *giveSolutionStepWhenIcApply();
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int checkConsistency();

    // identification
    virtual const char *giveClassName() const { return "NonStationaryTransportProblem"; }
    virtual classType giveClassID() const { return NonStationaryTransportProblemClass; }
    virtual fMode giveFormulation() { return TL; }

    /// Allows to change number of equations during solution.
    virtual int requiresUnknownsDictionaryUpdate() { return changingProblemSize; }
    virtual bool requiresEquationRenumbering(TimeStep *) { return changingProblemSize; }
    //Store solution vector to involved DoFs
    //virtual void updateDofUnknownsDictionary(DofManager *dman, TimeStep *tStep);

    virtual int giveUnknownDictHashIndx(EquationID type, ValueModeType mode, TimeStep *stepN);

    virtual void giveElementCharacteristicMatrix(FloatMatrix &answer, int num,
                                                 CharType type, TimeStep *tStep, Domain *domain);

    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime);

    /**
     * Returns time function for time step increment.
     * Used time function should provide step lengths as function of step number.
     * Initial step with number 0 is considered as [ -dt(0), 0 ], first step is [ 0, dt(1) ], ...
     */
    LoadTimeFunction *giveDtTimeFunction();

    /**
     * Returns the time step length for given step number n, initial step is number 0.
     */
    double giveDeltaT(int n);

    void averageOverElements(TimeStep *tStep);

#ifdef __PETSC_MODULE
    virtual void initPetscContexts();
#endif

protected:
    virtual void assembleAlgorithmicPartOfRhs(FloatArray &rhs, EquationID ut,
                                              const UnknownNumberingScheme &s, TimeStep *tStep);

    /**
     * This function is normally called at the first time to project initial conditions to previous (0^th) solution vector.
     * @param tStep Previous solution step.
     */
    virtual void applyIC(TimeStep *tStep);

    /**
     * Assembles part of RHS due to Dirichlet boundary conditions.
     * @param answer Global vector where the contribution will be added.
     * @param tStep Solution step.
     * @param eid Equation ID.
     * @param mode Mode of result.
     * @param lhsType Type of element matrix to be multiplied by vector of prescribed.
     * The giveElementCharacteristicMatrix service is used to get/compute element matrix.
     * @param s A map of non-default equation numbering if required.
     * @param d Domain.
     */
    virtual void assembleDirichletBcRhsVector(FloatArray &answer, TimeStep *tStep, EquationID eid, ValueModeType mode,
                                              CharType lhsType, const UnknownNumberingScheme &s, Domain *d);

    /**
     * Updates IP values on elements.
     * @param stepN Solution step.
     */
    virtual void updateInternalState(TimeStep *stepN);
};
} // end namespace oofem
#endif // nonstationarytransportproblem_h
