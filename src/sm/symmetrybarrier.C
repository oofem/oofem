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

#include "symmetrybarrier.h"
#include "intarray.h"
#include "flotarry.h"
#include "mathfem.h"

namespace oofem {
SymmetryBarrier :: SymmetryBarrier(int n, Domain *aDomain) :
    NonlocalBarrier(n, aDomain), origin(), normals(), mask(), lcs(3, 3)
    // Constructor. Creates an element with number n, belonging to aDomain.
{ }


SymmetryBarrier :: ~SymmetryBarrier()
// Destructor.
{ }

void
SymmetryBarrier :: applyConstraint(const FloatArray &c1, const FloatArray &c2, double &weight,
                                   bool &shieldFlag, NonlocalMaterialExtensionInterface *nei)
{
    // compute node coordinates in barrrier lcs
    FloatArray mc2(3), help(3);
    int i, ii, jj, kk, mi, dim = c1.giveSize();
    double d;

    shieldFlag = false;


    for ( i = 1; i <= dim; i++ ) {
        help.at(i) = c2.at(i) - origin.at(i);
    }

    // first compute mirrors to active planes
    // loop over active planes, mirror source point and compute weights
    for ( mi = 1; mi <= 3; mi++ ) {
        if ( mask.at(mi) ) {
            // plane active
            // mirror source
            for ( d = 0.0, i = 1; i <= dim; i++ ) {
                d += help.at(i) * lcs.at(mi, i);
            }

            for ( i = 1; i <= dim; i++ ) {
                mc2.at(i) = c2.at(i) - 2.0 *d *lcs.at(mi, i);
            }

            // compute weight of mirrored source
            weight += nei->computeWeightFunction(c1, mc2);
        }
    }

    // compute mirrors to lines common to two active planes
    for ( mi = 0; mi < 3; mi++ ) {
        ii = mi + 1;
        jj = ( mi + 1 ) % 3 + 1;
        kk = ( mi + 2 ) % 3 + 1;

        if ( mask.at(ii) && mask.at(jj) ) {
            // compute mirror point
            for ( d = 0.0, i = 1; i <= dim; i++ ) {
                d += help.at(i) * lcs.at(kk, i);
            }

            d = sqrt(d);
            for ( i = 1; i <= dim; i++ ) {
                mc2.at(i) = c2.at(i) - 2.0 * ( c2.at(i) - ( origin.at(i) + d * lcs.at(kk, i) ) );
            }

            // compute weight of mirrored source
            weight += nei->computeWeightFunction(c1, mc2);
        }
    }

    // finally compute mirror to origin if all three planes are active
    if ( mask.at(1) && mask.at(2) && mask.at(3) ) {
        // mirror source
        for ( i = 1; i <= dim; i++ ) {
            mc2.at(i) = c2.at(i) - 2.0 * ( c2.at(i) - origin.at(i) );
        }

        // compute weight of mirrored source
        weight += nei->computeWeightFunction(c1, mc2);
    }
}

IRResultType
SymmetryBarrier :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    FloatArray normals;
    int j;

    IR_GIVE_FIELD(ir, origin, IFT_SymmetryBarrier_origin, "origin"); // Macro
    IR_GIVE_FIELD(ir, normals, IFT_SymmetryBarrier_normals, "normals"); // Macro


    lcs.resize(3, 3);
    int size = normals.giveSize();
    if ( !( ( size == 0 ) || ( size == 6 ) ) ) {
        _warning2( "initializeFrom: lcs in node %d is not properly defined, will be ignored", this->giveNumber() );
    }

    if ( size == 6 ) {
        double n1 = 0.0, n2 = 0.0;
        // compute transformation matrix
        for ( j = 1; j <= 3; j++ ) {
            lcs.at(1, j) = normals.at(j);
            n1 += normals.at(j) * normals.at(j);
            lcs.at(2, j) = normals.at(j + 3);
            n2 += normals.at(j + 3) * normals.at(j + 3);
        }

        n1 = sqrt(n1);
        n2 = sqrt(n2);
        if ( ( n1 <= 1.e-6 ) || ( n2 <= 1.e-6 ) ) {
            _error("instanciateFrom : lcs input error");
        }

        for ( j = 1; j <= 3; j++ ) { // normalize e1' e2'
            lcs.at(1, j) /= n1;
            lcs.at(2, j) /= n2;
        }

        // vector e3' computed from vector product of e1', e2'
        lcs.at(3, 1) = lcs.at(1, 2) * lcs.at(2, 3) - lcs.at(1, 3) * lcs.at(2, 2);
        lcs.at(3, 2) = lcs.at(1, 3) * lcs.at(2, 1) - lcs.at(1, 1) * lcs.at(2, 3);
        lcs.at(3, 3) = lcs.at(1, 1) * lcs.at(2, 2) - lcs.at(1, 2) * lcs.at(2, 1);
    }

    IR_GIVE_FIELD(ir, mask, IFT_SymmetryBarrier_activemask, "activemask"); // Macro
    if ( mask.giveSize() != 3 ) {
        _error("instanciateFrom: activemask size should be 3");
    }

    return IRRT_OK;
}
} // end namespace oofem
