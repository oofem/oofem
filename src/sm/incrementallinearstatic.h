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

#ifndef incrementallinearstatic_h
#define incrementallinearstatic_h

#include "structengngmodel.h"
#include "sparselinsystemnm.h"
#include "sparsemtrxtype.h"

namespace oofem {
/**
 * This class implements Incremental LinearStatic Engineering problem.
 * Problem is solved as series of explicit linear solutions.
 * It will solve any nonlinear problem by writing the expressions explicitly
 * @f$ K({}^{n-1} r)\cdot \delta r = R({}^{n-1}r) - F @f$.
 * Explicit solver has requirements on the time step size to obtain a stable solution.
 *
 * This class can be used for solving linear creep problems.
 *
 * Supports the changes of static scheme (applying, removing and changing  boundary conditions)
 * during the analysis.
 *
 */
class IncrementalLinearStatic : public StructuralEngngModel
{
protected:
    SparseMtrx *stiffnessMatrix;
    FloatArray loadVector;
    FloatArray internalLoadVector;
    FloatArray incrementOfDisplacementVector;
    FloatArray discreteTimes;
    bool fixedSteps;
    double deltaT;
    double endOfTimeOfInterest;
    /// Numerical method used to solve the problem
    SparseLinearSystemNM *nMethod;

    LinSystSolverType solverType;
    SparseMtrxType sparseMtrxType;

public:
    IncrementalLinearStatic(int i, EngngModel *_master = NULL);
    virtual ~IncrementalLinearStatic();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void solveYourself();
    virtual void solveYourselfAt(TimeStep *tStep);

    virtual void updateYourself(TimeStep *tStep);
    virtual double giveUnknownComponent(EquationID type, ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof);
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual TimeStep *giveNextStep();

    /**
     * This function returns time valid for iStep time step, used in integration
     * of structure response.
     * This functions, when invoked for the first time, generates table of times.
     * Times in this table are generated according to:
     * - if there exists some load time function with abrupt change,
     *   then time just before and just at abrupt is included in table.
     * - between these steps under constant loads (or when load changes continuously)
     *   we use progressively increasing time step. They are best chosen so that time
     *   step be kept constant in the log (t-t') scale
     * @param iStep Time step number.
     */
    double giveDiscreteTime(int iStep);
    virtual double giveEndOfTimeOfInterest() { return endOfTimeOfInterest; }

    virtual NumericalMethod *giveNumericalMethod(TimeStep *tStep);

    virtual void printDofOutputAt(FILE *stream, Dof *iDof, TimeStep *atTime);

    virtual void terminate(TimeStep *tStep);

    virtual fMode giveFormulation() { return TL; }

    virtual const char *giveClassName() const { return "IncrementalLinearStatic"; }
    virtual classType giveClassID() const { return IncrementalLinearStaticClass; }

    virtual int requiresUnknownsDictionaryUpdate() { return true; }
    virtual bool requiresEquationRenumbering(TimeStep *) { return true; }
    virtual void updateDofUnknownsDictionary(DofManager *, TimeStep *);
    // Here we store only total and incremental value; so hash is computed from mode value only
    virtual int giveUnknownDictHashIndx(EquationID type, ValueModeType mode, TimeStep *stepN) { return (int) mode; }
};
} // end namespace oofem
#endif // incrementallinearstatic_h
