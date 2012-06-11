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

#include "cebfip78.h"
#include "mathfem.h"
#include "gausspnt.h"
#include "crosssection.h"

namespace oofem {
IRResultType
CebFip78Material :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    MaxwellChainMaterial :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, E28, IFT_CebFip78Material_e28, "e28"); // Macro
    IR_GIVE_FIELD(ir, fibf, IFT_CebFip78Material_fibf, "fibf"); // Macro
    IR_GIVE_FIELD(ir, kap_a, IFT_CebFip78Material_kap_a, "kap_a"); // Macro
    IR_GIVE_FIELD(ir, kap_c, IFT_CebFip78Material_kap_c, "kap_c"); // Macro
    IR_GIVE_FIELD(ir, kap_tt, IFT_CebFip78Material_kap_tt, "kap_tt"); // Macro
    IR_GIVE_FIELD(ir, u, IFT_CebFip78Material_u, "u"); // Macro

    // t0   = readDouble (initString,"curringendtime");
    return IRRT_OK;
}


double
CebFip78Material :: computeCreepFunction(GaussPoint *gp, double atTime, double ofAge)
{
    // computes the value of creep function at time ofAge
    // when load is acting from atTime
    // WARNING: Area returned by crossSection is assumed to be in [m^2].

    double e0;
    double fi, fi0, firv, fiir, hd, alpha, beta;
    double t, t0;
    CrossSection *cs = gp->giveCrossSection();

    t0 = this->kap_tt * this->kap_c * ofAge;
    t  = this->kap_tt * this->kap_c * atTime;

    e0 = E28 * sqrt( 1.36 * t0 / ( t0 + 10. ) );

    fi0 = 0.95 * pow(t0, -0.3) - 0.1;
    if ( fi0 < 0. ) {
        fi0 = 0.;
    }

    firv = 0.11 + 0.2 * atan( 0.05 * pow( ( t - t0 ), 2. / 3. ) );
    if ( firv < 0.4 ) {
        firv = 0.4;
    }

    hd = this->kap_a * cs->give(CS_Area) * 1000. * 1000. / this->u;
    if ( hd < 50. ) {
        hd = 50.;
    }

    if ( hd > 1500. ) {
        hd = 1500.;
    }

    alpha = 0.078 * exp( -1.20 * log10(hd) );
    beta  = 0.530 * exp( -0.13 * log10(hd) );

    fiir = 0.25 * ( 5.4 - log10(hd) ) * ( exp( -pow(alpha * t0, beta) ) - exp( -pow(alpha * t, beta) ) );
    fiir = fibf * fiir;

    fi = fi0 + firv + fiir;

    return ( 1. / e0 ) + fi / E28;
}
} // end namespace oofem
