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

#ifndef PRESCRIBEDGRADIENTBCNEUMANN_H_
#define PRESCRIBEDGRADIENTBCNEUMANN_H_

#include "prescribedgradientbcweakperiodic.h"

#define _IFT_PrescribedGradientBCNeumann_Name   "prescribedgradientbcneumann"

namespace oofem {
/**
 * Imposes a prescribed gradient weakly on the boundary
 * with a Neumann boundary condition.
 *
 * To reduce duplication of code, we exploit the fact that
 * Neumann boundary conditions can be obtained as a special
 * case of weakly periodic boundary conditions.
 *
 * @author Erik Svenning
 * @date Mar 5, 2014
 */

class PrescribedGradientBCNeumann : public PrescribedGradientBCWeakPeriodic {
public:
	PrescribedGradientBCNeumann(int n, Domain * d);
	virtual ~PrescribedGradientBCNeumann();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual const char *giveClassName() const { return "PrescribedGradientBCNeumann"; }
    virtual const char *giveInputRecordName() const { return _IFT_PrescribedGradientBCNeumann_Name; }
};

} /* namespace oofem */

#endif /* PRESCRIBEDGRADIENTBCNEUMANN_H_ */
