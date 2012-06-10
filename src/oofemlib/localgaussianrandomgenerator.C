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

#include "localgaussianrandomgenerator.h"
#include "mathfem.h"

#include <time.h>

namespace oofem {
LocalGaussianRandomGenerator :: LocalGaussianRandomGenerator(int num, Domain *d) : RandomFieldGenerator(num, d)
{}

LocalGaussianRandomGenerator :: ~LocalGaussianRandomGenerator()
{}

void
LocalGaussianRandomGenerator :: generateRandomValue(double &value, FloatArray *position)
{
    value = normalCdfInverse(ran1(& randomInteger), mean, variance * mean);
}

IRResultType
LocalGaussianRandomGenerator :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, mean, IFT_LocalGaussianRandomGenerator_mean, "mean"); // Macro
    IR_GIVE_FIELD(ir, variance, IFT_LocalGaussianRandomGenerator_variance, "variance"); // Macro
    randomInteger = -time(NULL);
    int seed = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, seed, IFT_LocalGaussianRandomGenerator_seed, "seed"); // Macro
    if ( seed ) {
        randomInteger = seed;
    }

    return IRRT_OK;
}


#define IA 16807
#define IM 2147483647
#define AM ( 1.0 / IM )
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV ( 1 + ( IM - 1 ) / NTAB )
#define EPS 1.2e-7
#define RNMX ( 1.0 - EPS )

double LocalGaussianRandomGenerator :: ran1(long *idum)
{
    int j;
    long k;
    static long iy = 0;
    static long iv [ NTAB ];
    float temp;

    if ( * idum <= 0 || !iy ) {
        if ( -( * idum ) < 1 ) {
            * idum = 1;
        } else {
            * idum = -( * idum );
        }

        for ( j = NTAB + 7; j >= 0; j-- ) {
            k = ( * idum ) / IQ;
            * idum = IA * ( * idum - k * IQ ) - IR * k;
            if ( * idum < 0 ) {
                * idum += IM;
            }

            if ( j < NTAB ) {
                iv [ j ] = * idum;
            }
        }

        iy = iv [ 0 ];
    }

    k = ( * idum ) / IQ;
    * idum = IA * ( * idum - k * IQ ) - IR * k;
    if ( * idum < 0 ) {
        * idum += IM;
    }

    j = iy / NDIV;
    iy = iv [ j ];
    iv [ j ] = * idum;
    if ( ( temp = AM * iy ) > RNMX ) {
        return RNMX;
    } else {
        return temp;
    }
}

double LocalGaussianRandomGenerator :: normalCdfInverse(double cdf, double a, double b)
{
    double x;
    double x2;
    if ( cdf < 0.0 || 1.0 < cdf ) {
        OOFEM_ERROR("LocalGaussianRandomGenerator :: normalCdfInverse - NORMAL_CDF_INV - Fatal error!\nCDF < 0 or 1 < CDF.");
    }

    x2 = normal01CdfInverse(cdf);
    x = a + b * x2;
    return x;
}

double LocalGaussianRandomGenerator :: normal01CdfInverse(double p)
{
    double a [ 8 ] = {
        3.3871328727963666080,     1.3314166789178437745e+2,
        1.9715909503065514427e+3,  1.3731693765509461125e+4,
        4.5921953931549871457e+4,  6.7265770927008700853e+4,
        3.3430575583588128105e+4,  2.5090809287301226727e+3
    };
    double b [ 8 ] = {
        1.0,                       4.2313330701600911252e+1,
        6.8718700749205790830e+2,  5.3941960214247511077e+3,
        2.1213794301586595867e+4,  3.9307895800092710610e+4,
        2.8729085735721942674e+4,  5.2264952788528545610e+3
    };
    double c [ 8 ] = {
        1.42343711074968357734,     4.63033784615654529590,
        5.76949722146069140550,     3.64784832476320460504,
        1.27045825245236838258,     2.41780725177450611770e-1,
        2.27238449892691845833e-2,  7.74545014278341407640e-4
    };
    double const1 = 0.180625;
    double const2 = 1.6;
    double d [ 8 ] = {
        1.0,                        2.05319162663775882187,
        1.67638483018380384940,     6.89767334985100004550e-1,
        1.48103976427480074590e-1,  1.51986665636164571966e-2,
        5.47593808499534494600e-4,  1.05075007164441684324e-9
    };
    double e [ 8 ] = {
        6.65790464350110377720,     5.46378491116411436990,
        1.78482653991729133580,     2.96560571828504891230e-1,
        2.65321895265761230930e-2,  1.24266094738807843860e-3,
        2.71155556874348757815e-5,  2.01033439929228813265e-7
    };
    double f [ 8 ] = {
        1.0,                        5.99832206555887937690e-1,
        1.36929880922735805310e-1,  1.48753612908506148525e-2,
        7.86869131145613259100e-4,  1.84631831751005468180e-5,
        1.42151175831644588870e-7,  2.04426310338993978564e-15
    };
    double q;
    double r;
    double split1 = 0.425;
    double split2 = 5.0;
    double value;

    if ( p <= 0.0 ) {
        value = -HUGE_VAL;
        return value;
    }

    if ( 1.0 <= p ) {
        value = HUGE_VAL;
        return value;
    }

    q = p - 0.5;

    if ( fabs(q) <= split1 ) {
        r = const1 - q * q;
        value = q * dpolyValue(8, a, r) / dpolyValue(8, b, r);
    } else   {
        if ( q < 0.0 ) {
            r = p;
        } else   {
            r = 1.0 - p;
        }

        if ( r <= 0.0 ) {
            value = -1.0;
            OOFEM_ERROR("LocalGaussianRandomGenerator :: normal01CdfInverse - r < 0.0!");
        }

        r = sqrt( -log(r) );
        if ( r <= split2 ) {
            r = r - const2;
            value = dpolyValue(8, c, r) / dpolyValue(8, d, r);
        } else   {
            r = r - split2;
            value = dpolyValue(8, e, r) / dpolyValue(8, f, r);
        }

        if ( q < 0.0 ) {
            value = -value;
        }
    }

    return value;
}


double LocalGaussianRandomGenerator :: dpolyValue(int n, double a[], double x)
{
    int i;
    double value;
    value = 0.0;
    for ( i = n - 1; 0 <= i; i-- ) {
        value = value * x + a [ i ];
    }

    return value;
}
} // end namespace oofem
