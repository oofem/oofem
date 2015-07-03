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

#include "doublepowerlaw.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(DoublePowerLawMaterial);

IRResultType
DoublePowerLawMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, E28, _IFT_DoublePowerLawMaterial_e28);
    IR_GIVE_FIELD(ir, fi1, _IFT_DoublePowerLawMaterial_fi1);
    IR_GIVE_FIELD(ir, m, _IFT_DoublePowerLawMaterial_m);
    IR_GIVE_FIELD(ir, n, _IFT_DoublePowerLawMaterial_n);
    IR_GIVE_FIELD(ir, alpha, _IFT_DoublePowerLawMaterial_alpha);

    return MaxwellChainMaterial :: initializeFrom(ir);
}


double
DoublePowerLawMaterial :: computeCreepFunction(double t, double t_prime)
{

    // computes the value of creep function at time t
    // when load is acting from time t_prime
    // t-t_prime = duration of loading

    double e0;
    double h1, h2, h3;

    e0 = 1.50 * E28;

    h1 = pow(t - t_prime, n);
    h2 = pow(t_prime, -m) + alpha;
    h3 = fi1 / e0;

    return 1. / e0 + h1 * h2 * h3;
}
} // end namespace oofem
