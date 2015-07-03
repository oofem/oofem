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

#include "cebfip78.h"
#include "mathfem.h"
#include "gausspoint.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Material(CebFip78Material);

IRResultType
CebFip78Material :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, E28, _IFT_CebFip78Material_e28);
    IR_GIVE_FIELD(ir, fibf, _IFT_CebFip78Material_fibf);
    IR_GIVE_FIELD(ir, kap_a_per_area, _IFT_CebFip78Material_kap_a_per_area);
    IR_GIVE_FIELD(ir, kap_c, _IFT_CebFip78Material_kap_c);
    IR_GIVE_FIELD(ir, kap_tt, _IFT_CebFip78Material_kap_tt);
    IR_GIVE_FIELD(ir, u, _IFT_CebFip78Material_u);

    return MaxwellChainMaterial :: initializeFrom(ir);
}


double
CebFip78Material :: computeCreepFunction(double t, double t_prime)
{
    // computes the value of creep function at time t
    // when load is acting from time t_prime
    // t-t_prime = duration of loading

    double e0;
    double fi, fi0, firv, fiir, hd, alpha, beta;
    double tt, t0;

    t0 = this->kap_tt * this->kap_c * t_prime;
    tt  = this->kap_tt * this->kap_c * t;

    e0 = E28 * sqrt( 1.36 * t0 / ( t0 + 10. ) );

    fi0 = 0.95 * pow(t0, -0.3) - 0.1;
    if ( fi0 < 0. ) {
        fi0 = 0.;
    }

    firv = 0.11 + 0.2 * atan( 0.05 * pow( ( tt - t0 ), 2. / 3. ) );
    if ( firv < 0.4 ) {
        firv = 0.4;
    }

    hd = this->kap_a_per_area * 1000. * 1000. / this->u;
    if ( hd < 50. ) {
        hd = 50.;
    }

    if ( hd > 1500. ) {
        hd = 1500.;
    }

    alpha = 0.078 * exp( -1.20 * log10(hd) );
    beta  = 0.530 * exp( -0.13 * log10(hd) );

    fiir = 0.25 * ( 5.4 - log10(hd) ) * ( exp( -pow(alpha * t0, beta) ) - exp( -pow(alpha * tt, beta) ) );
    fiir = fibf * fiir;

    fi = fi0 + firv + fiir;

    return ( 1. / e0 ) + fi / E28;
}
} // end namespace oofem
