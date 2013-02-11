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

#ifndef scalarerrorindicator_h
#define scalarerrorindicator_h

#include "compiler.h"
#include "errorestimator.h"
#include "internalstatetype.h"

namespace oofem {
class RemeshingCriteria;

/**
 * The class representing scalar error indicator.
 * It indicates element error based on the value of some suitable scalar value obtained from the
 * element integration points and corresponding material model.
 */
class ScalarErrorIndicator : public ErrorEstimator
{
protected:
    /// Type of internal variable to be indicator (type for temp and nontemp version).
    int indicatorType;
    /// Corresponding internal state type.
    InternalStateType varType;

public:
    /// Constructor
    ScalarErrorIndicator(int n, Domain *d) : ErrorEstimator(n, d) { eeType = EET_SEI; }
    /// Destructor
    virtual ~ScalarErrorIndicator() { }

    virtual double giveElementError(EE_ErrorType type, Element *elem, TimeStep *tStep);
    virtual double giveValue(EE_ValueType type, TimeStep *tStep) { return 0.0; }
    virtual int estimateError(EE_ErrorMode mode, TimeStep *tStep);
    virtual RemeshingCriteria *giveRemeshingCrit();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual const char *giveClassName() const { return "ScalarErrorIndicator"; }
    virtual classType giveClassID() const { return ScalarErrorIndicatorClass; }
};
} // end namespace oofem
#endif // scalarerrorindicator_h
