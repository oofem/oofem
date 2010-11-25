/* $Header: /home/cvs/bp/oofem/tm/src/nltransienttransportproblem.h,v 1.1 2003/04/14 16:01:39 bp Exp $ */
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
// Class NonLinearTransientTransportProblem
//

#ifndef nltransienttransportproblem_h
#define nltransienttransportproblem_h

#include "nonstationarytransportproblem.h"
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#ifndef __MAKEDEPEND
 #include <stdio.h>
#endif

namespace oofem {
/**
 * This class represents nonlinear transient transport problem.
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
    NLTransientTransportProblem(int i, EngngModel *_master = NULL) : NonStationaryTransportProblem(i, _master)  { }
    ~NLTransientTransportProblem() { }

    void solveYourselfAt(TimeStep *);
    /**
     * Updates nodal values
     * (calls also this->updateDofUnknownsDictionary for updating dofs unknowns dictionaries
     * if model supports changes of static system). The element internal state update is also forced using
     * updateInternalState service.
     */
    virtual void               updateYourself(TimeStep *);
    double giveUnknownComponent(EquationID, ValueModeType, TimeStep *, Domain *, Dof *);

    /// Initialization from given input record
    IRResultType initializeFrom(InputRecord *ir);

    // identification
    const char *giveClassName() const { return "NLTransientTransportProblem"; }
    classType giveClassID()      const { return NLTransientTransportProblemClass; }
    fMode giveFormulation() { return nonLinFormulation; }

protected:
    /**
     * Updates nodal values and IP values on elements
     * (calls also this->updateDofUnknownsDictionary for updating dofs unknowns dictionaries
     * if model supports changes of static system). The element internal state update is also forced using
     * updateInternalState service.
     */
    void updateInternalState(TimeStep *);
    void applyIC(TimeStep *);
    void assembleAlgorithmicPartOfRhs(FloatArray &rhs, EquationID ut,
                                      const UnknownNumberingScheme &s, TimeStep *tStep, int nite);
};
} // end namespace oofem
#endif // nltransienttransportproblem_h
