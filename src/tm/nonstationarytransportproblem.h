/* $Header: /home/cvs/bp/oofem/tm/src/nonstationarytransportproblem.h,v 1.2 2003/05/19 13:04:10 bp Exp $ */
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

//
// Class NonStartionaryTransportProblem
//

#ifndef nonstationarytransportproblem_h
#define nonstationarytransportproblem_h

#include "structengngmodel.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "primaryfield.h"
#include "dofdistributedprimaryfield.h"
#include "transportmaterial.h"
#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

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
    ///Right hand side vector from boundary conditions
    FloatArray bcRhs;
    ///This field stores solution vector. For fixed size of problem, the PrimaryField is used, for growing/decreasing size, DofDistributedPrimaryField applies
    PrimaryField *UnknownsField;

    LinSystSolverType solverType;
    SparseMtrxType sparseMtrxType;
    /// Numerical method used to solve the problem
    SparseLinearSystemNM *nMethod;

    double deltaT;
    double alpha;

    /// if set then stabilization using lumped capacity will be used
    int lumpedCapacityStab;

    //// if set, the receiver flux field will be exported using FieldManager
    //int exportFieldFlag;

    /// Associated time function for time step increment
    int dtTimeFunction;

    /// changingProblemSize=0 means no change in the problem size (no application/removal of Dirichet boundary conditions)
    bool changingProblemSize;

public:
    ///constructor
    NonStationaryTransportProblem(int i, EngngModel *_master);
    ///destructor
    ~NonStationaryTransportProblem();

    void solveYourselfAt(TimeStep *);
    /**
     * Updates nodal values. The method calls also this->updateDofUnknownsDictionary for updating DOFs unknowns dictionaries,
     * because the model supports assignment of Dirichlet b.c. at various time steps. The number of equations may be different
     * in each timeStep. The element internal state update is also forced using updateInternalState service.
     */
    virtual void updateYourself(TimeStep *);
    double giveUnknownComponent(EquationID, ValueModeType, TimeStep *, Domain *, Dof *);
    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    void updateDomainLinks();

    TimeStep *giveNextStep();
    TimeStep *giveSolutionStepWhenIcApply();
    NumericalMethod *giveNumericalMethod(TimeStep *);

    /// Initialization from given input record
    IRResultType initializeFrom(InputRecord *ir);

    /// Consistency check
    virtual int checkConsistency(); // returns nonzero if o.k.

    // identification
    const char *giveClassName() const { return "NonStationaryTransportProblem"; }
    classType giveClassID() const { return NonStationaryTransportProblemClass; }
    fMode giveFormulation() { return TL; }

    ///Allows to change number of equations during solution
    virtual int       requiresUnknownsDictionaryUpdate() { return changingProblemSize; }
    virtual bool      requiresEquationRenumbering(TimeStep *) { return changingProblemSize; }
    //Store solution vector to involved DoFs
    //virtual void      updateDofUnknownsDictionary(DofManager *, TimeStep *);

    /**
     * The array of MasterDof::unknowns contains, in this case, one previous solution and previous RHS.
     * The hash index tells position in the array MasterDof::unknowns, depending on @param mode
     */
    virtual int giveUnknownDictHashIndx(EquationID type, ValueModeType mode, TimeStep *stepN);

    virtual void giveElementCharacteristicMatrix(FloatMatrix &answer, int num,
                                                 CharType type, TimeStep *tStep, Domain *domain);

    /** DOF printing routine. Called by DofManagers to print Dof specific part.
     * Dof class provides component printing routines, but emodel is responsible
     * for what will be printed at DOF level.
     * @param stream output stream
     * @param iDof dof to be processed
     * @param atTime solution step
     */
    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime);

    /**
     * Returns time function for time step increment.
     * Used time function should provide step lengths as function of step number.
     * Initial step with number 0 is considered as [ -dt(0), 0 ], first step is [ 0, dt(1) ], ...
     */
    LoadTimeFunction *giveDtTimeFunction();

    /**
     * Returns the timestep length for given step number n, initial step is number 0
     */
    double giveDeltaT(int n);


    void averageOverElements(TimeStep *tStep);


#ifdef __PETSC_MODULE
    /**
     * Creates Petsc contexts. Must be implemented by derived classes since the governing equation type is reqired
     * for context creation.
     */
    virtual void initPetscContexts();
#endif


protected:
    virtual void assembleAlgorithmicPartOfRhs(FloatArray &rhs, EquationID ut,
                                              const UnknownNumberingScheme &s, TimeStep *tStep);

    /**
     * This function is normally called at the first time to project initial conditions to previous (0^th) solution vector.
     * @param tStep previous solution step
     */
    virtual void applyIC(TimeStep *);

    /**
     * Assembles part of rhs due to Dirichlet boundary conditions.
     * @param answer global vector where the contribution will be added
     * @param tStep solution step
     * @param mode CharTypeMode of result
     * @param lhsType type of element matrix to be multiplied by vector of prescribed.
     * The giveElementCharacteristicMatrix service is used to get/compute element matrix.
     * @param s a map of non-default equation numbering if required
     * @param d domain
     */
    virtual void assembleDirichletBcRhsVector(FloatArray &answer, TimeStep *tStep, EquationID ut, ValueModeType mode,
                                              CharType lhsType, const UnknownNumberingScheme &s, Domain *d);

    /**
     * Updates IP values on elements
     * @param TimeStep solution step
     */
    virtual void updateInternalState(TimeStep *stepN);
};
} // end namespace oofem
#endif // nonstationarytransportproblem_h
