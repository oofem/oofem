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

#ifndef nltransienttransportproblem_h
#define nltransienttransportproblem_h

#include "tm/EngineeringModels/nonstationarytransportproblem.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"

///@name Input fields for NLTransientTransportProblem
//@{
#define _IFT_NLTransientTransportProblem_Name "nltransienttransportproblem"
#define _IFT_NLTransientTransportProblem_nsmax "nsmax"
#define _IFT_NLTransientTransportProblem_rtol "rtol"
#define _IFT_NLTransientTransportProblem_manrmsteps "manrmsteps"
//@}

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

    double rtol = 0.;
    int nsmax = 0;
    nlttp_ModeType NR_Mode = nrsolverModifiedNRM;
    int MANRMSteps = 0;
    int currentIterations = 0;

public:
    NLTransientTransportProblem(int i, EngngModel * _master);

    TimeStep* giveNextStep() override;
    void solveYourselfAt(TimeStep *tStep) override;
    void updateYourself(TimeStep *tStep) override;
    double giveUnknownComponent(ValueModeType mode, TimeStep *tStep, Domain *d, Dof *dof) override;

    void initializeFrom(InputRecord &ir) override;

    // identification
    const char *giveInputRecordName() const { return _IFT_NLTransientTransportProblem_Name; }
    const char *giveClassName() const override { return "NLTransientTransportProblem"; }
    fMode giveFormulation() override { return nonLinFormulation; }
    int giveUnknownDictHashIndx(ValueModeType mode, TimeStep *tStep) override;
    void updateDofUnknownsDictionary(DofManager *dman, TimeStep *tStep) override;

    int giveCurrentNumberOfIterations() override { return currentIterations; }

protected:
    void updateInternalState(TimeStep *tStep) override;
    void applyIC(TimeStep *tStep) override;
    void createPreviousSolutionInDofUnknownsDictionary(TimeStep *tStep);
    void assembleAlgorithmicPartOfRhs(FloatArray &rhs, const UnknownNumberingScheme &s, TimeStep *tStep) override;
};
} // end namespace oofem
#endif // nltransienttransportproblem_h
