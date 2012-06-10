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

#ifndef stationarytransportproblem_h
#define stationarytransportproblem_h

#include "structengngmodel.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "primaryfield.h"

namespace oofem {
/**
 * This class represents stationary transport problem.
 */
class StationaryTransportProblem : public EngngModel
{
protected:
    SparseMtrx *conductivityMatrix;
    FloatArray rhsVector;

    //FloatArray solutionVector;
    PrimaryField FluxField;

    LinSystSolverType solverType;
    SparseMtrxType sparseMtrxType;
    /// Numerical method used to solve the problem
    SparseLinearSystemNM *nMethod;
    //// if set, the receiver flux field will be exported using FieldManager
    //int exportFieldFlag;

public:
    StationaryTransportProblem(int i, EngngModel *_master = NULL) : EngngModel(i, _master), rhsVector(),
        FluxField(this, 1, FT_TransportProblemUnknowns, EID_ConservationEquation, 0)
    {
        conductivityMatrix = NULL;
        ndomains = 1;
        nMethod = NULL;
    }
    virtual ~StationaryTransportProblem()
    {
        delete  conductivityMatrix;
        if ( nMethod ) { delete nMethod; } }

    virtual void solveYourselfAt(TimeStep *tStep);
    virtual void updateYourself(TimeStep *tStep);
    virtual double giveUnknownComponent(EquationID eid, ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof);
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual void updateDomainLinks();

    virtual TimeStep *giveNextStep();
    virtual NumericalMethod *giveNumericalMethod(MetaStep *mStep);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int checkConsistency();

    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime);

    // identification
    virtual const char *giveClassName() const { return "StationaryTransportProblem"; }
    virtual classType giveClassID() const { return StationaryTransportProblemClass; }
    virtual fMode giveFormulation() { return TL; }

#ifdef __PETSC_MODULE
    virtual void initPetscContexts();
#endif

protected:
    /**
     * Assembles part of rhs due to Dirichlet boundary conditions.
     * @param answer Global vector where the contribution will be added.
     * @param tStep Solution step.
     * @param eid Equation ID.
     * @param mode CharTypeMode of result.
     * @param lhsType Type of element matrix to be multiplied by vector of prescribed.
     * The giveElementCharacteristicMatrix service is used to get/compute element matrix.
     * @param s Numbering scheme for unknowns.
     * @param d Domain.
     */
    void assembleDirichletBcRhsVector(FloatArray &answer, TimeStep *tStep, EquationID eid, ValueModeType mode,
                                      CharType lhsType, const UnknownNumberingScheme &s, Domain *d);

    /**
     * Updates IP values on elements
     * @param stepN Solution step.
     */
    virtual void updateInternalState(TimeStep *stepN);
};
} // end namespace oofem
#endif // stationarytransportproblem_h
