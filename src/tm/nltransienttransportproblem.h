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

#ifndef nltransienttransportproblem_h
#define nltransienttransportproblem_h

#include "nonstationarytransportproblem.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"

namespace oofem {
/**
 * This class represents nonlinear transient transport problem. The problem can be growing/decreasing, signalized by flag "changingproblemsize"
 * in the problem description. The solution is stored in UnknownsField, which can obtain/ project solution from/to DOFs (nodes). If the problem
 * keeps the same equation numbers, solution is taken from UnknownsField without any projection, which is more efficient. See the matlibmanual
 * for solution strategy of balance equations and the solution algorithm.
 *
 * @todo Documentation errors (there is no "UnknownsField" used here).
 */
class NLTransientTransportProblem : public NonStationaryTransportProblem
{
protected:
    enum nlttp_ModeType { nrsolverModifiedNRM, nrsolverFullNRM, nrsolverAccelNRM };

    double rtol;
    int nsmax;
    nlttp_ModeType NR_Mode;
    int MANRMSteps;

public:
    /// Constructor.
    NLTransientTransportProblem(int i, EngngModel *_master);
    /// Destructor.
    virtual ~NLTransientTransportProblem();

    virtual void solveYourselfAt(TimeStep *tStep);
    virtual void updateYourself(TimeStep *tStep);
    virtual double giveUnknownComponent(EquationID eid, ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof);

    virtual IRResultType initializeFrom(InputRecord *ir);

    // identification
    virtual const char *giveClassName() const { return "NLTransientTransportProblem"; }
    virtual classType giveClassID() const { return NLTransientTransportProblemClass; }
    virtual fMode giveFormulation() { return nonLinFormulation; }
    virtual int giveUnknownDictHashIndx(EquationID type, ValueModeType mode, TimeStep *stepN);
    virtual void updateDofUnknownsDictionary(DofManager *dman, TimeStep *tStep);
    /**
     * Copy unknowns in DOF's from previous to current position.
     * @param type Determines type of equation
     * @param mode What the unknown describes (increment, total value etc.).
     * @param fromTime From which time step to obtain value.
     * @param toTime To which time to copy.
     */
    virtual void copyUnknownsInDictionary(EquationID type, ValueModeType mode, TimeStep *fromTime, TimeStep *toTime);

protected:
    virtual void updateInternalState(TimeStep *tStep);
    virtual void applyIC(TimeStep *tStep);
    void createPreviousSolutionInDofUnknownsDictionary(TimeStep *tStep);
    void assembleAlgorithmicPartOfRhs(FloatArray &rhs, EquationID ut, const UnknownNumberingScheme &s, TimeStep *tStep, int nite);
};
} // end namespace oofem
#endif // nltransienttransportproblem_h
