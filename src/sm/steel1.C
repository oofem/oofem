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

#include "steel1.h"
#include "isolinearelasticmaterial.h"
#include "mathfem.h"

namespace oofem {
Steel1 :: Steel1(int n, Domain *d) : PerfectlyPlasticMaterial(n, d)
    // constructor
{
    linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
}



double
Steel1 :: computeYCValueAt(GaussPoint *gp,
                           FloatArray *currentStress,
                           FloatArray *currentPlasticStrain)

//
// computes the value of receiver.
// testing if yield condition is fulfilled
// function double      computeValueAt (GaussPoint*)
// should return zero if yield surface is reached,
//               less than zero if surface still not reached
//               > 0. if yield criteria previously violated.
//
// ig gp->status are posibly stored hardenining variables
{
    // double answer;
    return this->computeJ2InvariantAt(currentStress) - this->give('k', gp);
}


double
Steel1 :: computeJ2InvariantAt(FloatArray *currentStress)
//
// computes the J2 value of receiver.
//
{
    double answer;
    double v1, v2, v3;

    if ( currentStress == NULL ) {
        return 0.0;
    }

    v1 = ( ( currentStress->at(1) - currentStress->at(2) ) * ( currentStress->at(1) - currentStress->at(2) ) );
    v2 = ( ( currentStress->at(2) - currentStress->at(3) ) * ( currentStress->at(2) - currentStress->at(3) ) );
    v3 = ( ( currentStress->at(3) - currentStress->at(1) ) * ( currentStress->at(3) - currentStress->at(1) ) );

    answer = ( 1. / 6. ) * ( v1 + v2 + v3 ) + currentStress->at(4) * currentStress->at(4) +
             currentStress->at(5) * currentStress->at(5) + currentStress->at(6) * currentStress->at(6);

    return sqrt(answer);
}




FloatArray *
Steel1 :: GiveYCStressGradient(GaussPoint *gp,
                               FloatArray *currentStress,
                               FloatArray *currentPlasticStrain)

//
// - returning vector of derivatives of yield surface with respect to stresses.
//
// ig gp->status are posibly stored hardening variables
{
    double f, sigx, sigy, sigz, sx, sy, sz;
    FloatArray *answer = new FloatArray(6);

    if ( currentStress == NULL ) {
        return answer;
    }

    f = this->computeJ2InvariantAt(currentStress);
    sigx = currentStress->at(1);
    sigy = currentStress->at(2);
    sigz = currentStress->at(3);

    sx = ( 2. / 3. ) * sigx - ( 1. / 3. ) * sigy - ( 1. / 3. ) * sigz;
    sy = ( 2. / 3. ) * sigy - ( 1. / 3. ) * sigx - ( 1. / 3. ) * sigz;
    sz = ( 2. / 3. ) * sigz - ( 1. / 3. ) * sigy - ( 1. / 3. ) * sigx;

    answer->at(1) = 0.5 * sx / f;
    answer->at(2) = 0.5 * sy / f;
    answer->at(3) = 0.5 * sz / f;
    answer->at(4) = currentStress->at(4) / f;
    answer->at(5) = currentStress->at(5) / f;
    answer->at(6) = currentStress->at(6) / f;

    return answer;
}


FloatArray *
Steel1 :: GiveLCStressGradient(GaussPoint *gp,
                               FloatArray *currentStress,
                               FloatArray *currentPlasticStrain)

//
// - returning vector of derivatives of loading surface with respect to stresses.
//
// ig gp->status are posibly stored hardenining variables
{
    return this->GiveYCStressGradient(gp, currentStress, currentPlasticStrain);
}


FloatArray *
Steel1 :: GiveYCPlasticStrainGradient(GaussPoint *gp,
                                      FloatArray *currentStress,
                                      FloatArray *currentPlasticStrain)

//
//returning vector of derivatives of yield surface with respects to plastic strains
//
// ig gp->status are posibly stored hardenining variables
{
    return new FloatArray(6);
}


FloatArray *
Steel1 :: GiveLCPlasticStrainGradient(GaussPoint *gp,
                                      FloatArray *currentStress,
                                      FloatArray *currentPlasticStrain)

//
//returning vector of derivatives of loading  surface with respects to plastic strains
//
// ig gp->status are posibly stored hardenining variables
{
    return this->GiveYCPlasticStrainGradient(gp, currentStress, currentPlasticStrain);
}


IRResultType
Steel1 :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    double value;

    PerfectlyPlasticMaterial :: initializeFrom(ir);
    IR_GIVE_FIELD(ir, value, IFT_Steel1_ry, "ry"); // Macro
    propertyDictionary->add( 'k', value / sqrt(3.) );

    return IRRT_OK;
}
} // end namespace oofem
