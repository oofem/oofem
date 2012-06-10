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

#include "doublepowerlaw.h"
#include "mathfem.h"

namespace oofem {
IRResultType
DoublePowerLawMaterial :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    MaxwellChainMaterial :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, E28, IFT_DoublePowerLawMaterial_e28, "e28"); // Macro
    IR_GIVE_FIELD(ir, fi1, IFT_DoublePowerLawMaterial_fi1, "fi1"); // Macro
    IR_GIVE_FIELD(ir, m, IFT_DoublePowerLawMaterial_m, "m"); // Macro
    IR_GIVE_FIELD(ir, n, IFT_DoublePowerLawMaterial_n, "n"); // Macro
    IR_GIVE_FIELD(ir, alpha, IFT_DoublePowerLawMaterial_alpha, "alpha"); // Macro

    return IRRT_OK;
}


double
DoublePowerLawMaterial :: computeCreepFunction(GaussPoint *gp, double atTime, double ofAge)
{
    // computes the value of creep function at time ofAge
    // when load is acting from atTime
    // WARNING: Area returned by crossSection is assumed to be in [m^2].

    double e0;
    double h1, h2, h3;

    e0 = 1.50 * E28;

    h1 = __OOFEM_POW(atTime - ofAge, n);
    h2 = __OOFEM_POW(ofAge, -m) + alpha;
    h3 = fi1 / e0;

    return 1. / e0 + h1 * h2 * h3;
}
} // end namespace oofem
