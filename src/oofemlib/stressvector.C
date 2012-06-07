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

#include "stressvector.h"
#include "strainvector.h"
#include "mathfem.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "error.h"

namespace oofem {
StressVector :: StressVector(MaterialMode m) : StressStrainBaseVector(m)
{ }

StressVector :: StressVector(const FloatArray &src, MaterialMode m) : StressStrainBaseVector(src, m)
{ }

void
StressVector :: computeDeviatoricVolumetricSplit(StressVector &dev, double &vol) const
{
    MaterialMode myMode = this->giveStressStrainMode();

    if ( myMode == _1dMat ) {
        // 1D model
        OOFEM_ERROR("StressVector::computeDeviatoricVolumetricSplit: No Split for 1D!");
        //  dev.resize(1); dev.at(1) = 0.0;
        // vol = this->at (1);
    } else if ( myMode == _PlaneStress ) {
        // plane stress problem
        OOFEM_ERROR("StressVector::computeDeviatoricVolumetricSplit: No Split for plane stress!");
        //  dev = *this;
        // vol = (this->at(1)+this->at(2))/2.0;
        //    dev.at(1) -= vol;
        // dev.at(2) -= vol;
    } else {
        // 3d, plane strain or axisymmetric problem
        dev = * this;
        vol = ( this->at(1) + this->at(2) + this->at(3) ) / 3.0;
        dev.at(1) -= vol;
        dev.at(2) -= vol;
        dev.at(3) -= vol;
    }
}

void
StressVector :: computeDeviatoricVolumetricSum(StressVector &answer, double vol) const
{
    MaterialMode myMode = this->giveStressStrainMode();

    if ( myMode == _1dMat ) {
        // 1D model
        OOFEM_ERROR("StressVector::computeDeviatoricVolumetricSum : No sum for 1D!");
        //  dev.resize(1); dev.at(1) = 0.0;
        // vol = this->at (1);
    } else if ( myMode == _PlaneStress ) {
        // plane stress problem
        OOFEM_ERROR("StressVector::computeDeviatoricVolumetricSum : No sum for plane stress!");
        //  dev = *this;
        // vol = (this->at(1)+this->at(2))/2.0;
        //    dev.at(1) -= vol;
        // dev.at(2) -= vol;
    } else {
        // 3d, plane strain or axisymmetric problem
        answer = * this;
        for ( int i = 0; i < 3; i++ ) {
            answer(i) += vol;
        }
    }
}

void
StressVector :: computePrincipalValues(FloatArray &answer) const
{
    //
    // This function cumputes Principal values of strain vector.
    // Engineering notation is used.
    //
    MaterialMode myMode = this->giveStressStrainMode();
    double swap;
    int size = this->giveSize();
    int nonzeroFlag = 0;

    if ( myMode == _1dMat ) {
        answer = * this;
    } else if ( myMode == _PlaneStress ) {
        // 2D problem
        double ast, dst, D = 0.0;
        answer.resize(2);

        for ( int i = 1; i <= size; i++ ) {
            if ( fabs( this->at(i) ) > 1.e-20 ) {
                nonzeroFlag = 1;
            }
        }

        if ( nonzeroFlag == 0 ) {
            answer.zero();
            return;
        }

        ast = this->at(1) + this->at(2);
        dst = this->at(1) - this->at(2);
        D = dst * dst + 4.0 * this->at(3) * this->at(3);

        if ( D < 0. ) {
            OOFEM_ERROR("StressVector::computePrincipalValues: Imaginar roots ");
        }

        D = sqrt(D);
        answer.at(1) = 0.5 * ( ast - D );
        answer.at(2) = 0.5 * ( ast + D );

        // sort result
        if ( answer.at(1) > answer.at(2) ) {
            return;
        } else {
            swap = answer.at(1);
            answer.at(1) = answer.at(2);
            answer.at(2) = swap;
            return;
        }
    } else {
        // 3D problem
        double I1 = 0.0, I2 = 0.0, I3 = 0.0, s1, s2, s3;
        int i, j;
        FloatArray s;

        this->convertToFullForm(s);
        for ( i = 1; i <= size; i++ ) {
            if ( fabs( s.at(i) ) > 1.e-20 ) {
                nonzeroFlag = 1;
            }
        }

        answer.resize(3);
        answer.zero();
        if ( nonzeroFlag == 0 ) {
            return;
        }

        I1 = s.at(1) + s.at(2) + s.at(3);
        I2 = s.at(1) * s.at(2) + s.at(2) * s.at(3) + s.at(3) * s.at(1) -
             ( s.at(4) * s.at(4) + s.at(5) * s.at(5) + s.at(6) * s.at(6) );
        I3 = s.at(1) * s.at(2) * s.at(3) + 2. * s.at(4) * s.at(5) * s.at(6) -
             ( s.at(1) * s.at(4) * s.at(4) + s.at(2) * s.at(5) * s.at(5) +
              s.at(3) * s.at(6) * s.at(6) );

        /*
         * Call cubic3r to ensure, that all three real eigenvalues will be found, because we have symmetric tensor.
         * This aloows to overcome various rounding errors when solving general cubic equation.
         */
        cubic3r( ( double ) -1., I1, -I2, I3, & s1, & s2, & s3, & i );

        if ( i > 0 ) {
            answer.at(1) = s1;
        }

        if ( i > 1 ) {
            answer.at(2) = s2;
        }

        if ( i > 2 ) {
            answer.at(3) = s3;
        }

        // sort results
        for ( i = 1; i < 3; i++ ) {
            for ( j = 1; j < 3; j++ ) {
                if ( answer.at(j + 1) > answer.at(j) ) {
                    swap = answer.at(j + 1);
                    answer.at(j + 1) = answer.at(j);
                    answer.at(j) = swap;
                }
            }
        }

        return;
    }
}

void
StressVector :: computePrincipalValDir(FloatArray &answer, FloatMatrix &dir) const
{
    //
    // This function cumputes Principal values & directions corresponding to receiver.
    //
    // Return Values:
    //
    // matrix dir -> principal directions of strains or stresses
    // array sp -> principal strains or stresses
    //

    FloatMatrix ss;
    FloatArray sp;
    double swap;
    int i, ii, jj, kk, nval;
    int nonzeroFlag = 0;
    int size = this->giveSize();
    MaterialMode myMode = this->giveStressStrainMode();

    if ( myMode == _1dMat ) {
        answer = * this;
        dir.resize(1, 1);
        dir.at(1, 1) = 1.0;
        return;
    } else if ( myMode == _PlaneStress ) {
        // 2D problem
        ss.resize(2, 2);
        answer.resize(2);

        for ( i = 1; i <= size; i++ ) {
            if ( fabs( this->at(i) ) > 1.e-20 ) {
                nonzeroFlag = 1;
            }
        }

        if ( nonzeroFlag == 0 ) {
            answer.zero();
            ss.zero();
            return;
        }

        ss.at(1, 1) = this->at(1);
        ss.at(2, 2) = this->at(2);
        ss.at(1, 2) = ss.at(2, 1) = this->at(3);
    } else {
        // 3D problem
        ss.resize(3, 3);
        FloatArray s;

        answer.resize(3);

        this->convertToFullForm(s);
        for ( i = 1; i <= size; i++ ) {
            if ( fabs( s.at(i) ) > 1.e-20 ) {
                nonzeroFlag = 1;
            }
        }

        if ( nonzeroFlag == 0 ) {
            answer.zero();
            ss.zero();
            return;
        }

        ss.at(1, 1) = s.at(1);
        ss.at(2, 2) = s.at(2);
        ss.at(3, 3) = s.at(3);
        ss.at(1, 2) = ss.at(2, 1) = s.at(6);
        ss.at(1, 3) = ss.at(3, 1) = s.at(5);
        ss.at(2, 3) = ss.at(3, 2) = s.at(4);
    }

#if 0
    ss.Jacobi(& answer, & dir, & i);
#else
    ss.jaco_(answer, dir, 10);
#endif
    // sort results
    nval = 3;
    if ( myMode == _PlaneStress ) {
        nval = 2;
    }

    for ( ii = 1; ii < nval; ii++ ) {
        for ( jj = 1; jj < nval; jj++ ) {
            if ( answer.at(jj + 1) > answer.at(jj) ) {
                // swap eigen values and eigen vectors
                swap = answer.at(jj + 1);
                answer.at(jj + 1) = answer.at(jj);
                answer.at(jj) = swap;
                for ( kk = 1; kk <= nval; kk++ ) {
                    swap = dir.at(kk, jj + 1);
                    dir.at(kk, jj + 1) = dir.at(kk, jj);
                    dir.at(kk, jj) = swap;
                }
            }
        }
    }
}

void
StressVector :: printYourself() const
{
    printf("StressVector (MaterialMode %d)\n", mode);
    for ( int i = 1; i <= size; ++i ) {
        printf( "%10.3e  ", this->at(i) );
    }

    printf("\n");
}



double
StressVector :: computeFirstInvariant() const
{
    //
    // This function computes the first invariant
    // of the stress vector
    //
    MaterialMode myMode = giveStressStrainMode();
    if ( myMode == _1dMat ) {
        // 1d problem
        return values [ 0 ];
    } else if ( myMode == _PlaneStress ) {
        // 2d problem: plane stress
        return values [ 0 ] + values [ 1 ];
    } else {
        // plane strain, axisymmetry or full 3d
        return values [ 0 ] + values [ 1 ] + values [ 2 ];
    }
}

double
StressVector :: computeSecondInvariant() const
{
    //
    // This function computes the second invariant of
    // of the deviatoric stress vector
    //
    MaterialMode myMode = giveStressStrainMode();
    if ( myMode == _1dMat ) {
        // 1d problem
        return .5 * values [ 0 ] * values [ 0 ];
    } else if ( myMode == _PlaneStress ) {
        // 2d problem: plane stress
        return .5 * ( values [ 0 ] * values [ 0 ] + values [ 1 ] * values [ 1 ] ) + values [ 2 ] * values [ 2 ];
    } else if ( myMode == _PlaneStrain || myMode == _3dRotContinuum ) {
        //  plane strain or axisymmetry
        return .5 * ( values [ 0 ] * values [ 0 ] + values [ 1 ] * values [ 1 ] + values [ 2 ] * values [ 2 ] ) +
               values [ 3 ] * values [ 3 ];
    } else {
        // 3d problem
        return .5 * ( values [ 0 ] * values [ 0 ] + values [ 1 ] * values [ 1 ] + values [ 2 ] * values [ 2 ] ) +
               values [ 3 ] * values [ 3 ] + values [ 4 ] * values [ 4 ] + values [ 5 ] * values [ 5 ];
    }
}

double
StressVector :: computeThirdInvariant() const
{
    //
    // This function computes the third invariant
    // of the deviatoric stress vector.
    //
    MaterialMode myMode = giveStressStrainMode();
    if ( myMode == _1dMat ) {
        // 1d problem
        return ( 1. / 3. ) * values [ 0 ] * values [ 0 ] * values [ 0 ];
    } else if ( myMode == _PlaneStress ) {
        // 2d problem: plane stress
        return ( 1. / 3. ) * ( values [ 0 ] * values [ 0 ] * values [ 0 ] + 3. * values [ 1 ] * values [ 2 ] * values [ 2 ]
                               + 3. * values [ 0 ] * values [ 2 ] * values [ 2 ] + values [ 1 ] * values [ 1 ] * values [ 1 ] );
    } else if ( myMode == _PlaneStrain || myMode == _3dRotContinuum ) {
        // plane strain or axisymmetry
        return ( 1. / 3. ) * ( values [ 0 ] * values [ 0 ] * values [ 0 ] + 3. * values [ 0 ] * values [ 3 ] * values [ 3 ] +
                               3. * values [ 1 ] * values [ 3 ] * values [ 3 ] + values [ 1 ] * values [ 1 ] * values [ 1 ] +
                               values [ 2 ] * values [ 2 ] * values [ 2 ] );
    } else {
        // 3d problem
        return ( 1. / 3. ) * ( values [ 0 ] * values [ 0 ] * values [ 0 ] + 3. * values [ 0 ] * values [ 5 ] * values [ 5 ] +
                               3. * values [ 0 ] * values [ 4 ] * values [ 4 ] + 6. * values [ 3 ] * values [ 5 ] * values [ 4 ] +
                               3. * values [ 1 ] * values [ 5 ] * values [ 5 ] + 3 * values [ 2 ] * values [ 4 ] * values [ 4 ] +
                               values [ 1 ] * values [ 1 ] * values [ 1 ] + 3. * values [ 1 ] * values [ 3 ] * values [ 3 ] +
                               3. * values [ 2 ] * values [ 3 ] * values [ 3 ] + values [ 2 ] * values [ 2 ] * values [ 2 ] );
    }
}


void
StressVector :: computeAllThreeHWCoordinates(double &xsi,
                                             double &rho,
                                             double &theta) const
{
    xsi = this->computeFirstCoordinate();

    StressVector deviatoricStress( this->giveStressStrainMode() );
    double volumetricStress;
    this->computeDeviatoricVolumetricSplit(deviatoricStress, volumetricStress);
    rho = deviatoricStress.computeSecondCoordinate();
    theta = deviatoricStress.computeThirdCoordinate();
}


double
StressVector :: computeFirstCoordinate() const
{
    //
    // This function computes the first Haigh-Westergaard coordinate
    // from the stress state
    //
    return computeFirstInvariant() / sqrt(3.);
}

double
StressVector :: computeSecondCoordinate() const
{
    //
    // This function computes the second Haigh-Westergaard coordinate
    // from the deviatoric stress state
    //
    ;
    return sqrt( 2. * computeSecondInvariant() );
}

double
StressVector :: computeThirdCoordinate() const
{
    //
    // This function computes the third Haigh-Westergaard coordinate
    // from the deviatoric stress state
    //
    double c1 = 0.0;
    if ( computeSecondInvariant() == 0. ) {
        c1 = 0.0;
    } else {
        c1 = ( 3. * sqrt(3.) / 2. ) * computeThirdInvariant() / ( pow( computeSecondInvariant(), ( 3. / 2. ) ) );
    }

    if ( c1 > 1.0 ) {
        c1 = 1.0;
    }

    if ( c1 < -1.0 ) {
        c1 = -1.0;
    }

    return 1. / 3. * acos(c1);
}

void
StressVector :: applyElasticCompliance(StrainVector &strain, const double EModulus, const double nu) const
{
    //
    // This function multiplies the receiver by the elastic compliance matrix
    // and stores the result in strain
    //
    MaterialMode myMode = giveStressStrainMode();
    if ( myMode == _1dMat ) {
        strain(0) = values [ 0 ] / EModulus;
    } else if ( myMode == _PlaneStress ) {
        strain(0) = ( values [ 0 ] - nu * values [ 1 ] ) / EModulus;
        strain(1) = ( -nu * values [ 0 ] + values [ 1 ] ) / EModulus;
        strain(2) = ( ( 2. + 2. * nu ) * values [ 2 ] ) / EModulus;
    } else if ( myMode == _PlaneStrain ) {
        strain(0) = ( values [ 0 ] - nu * values [ 1 ] - nu * values [ 2 ] ) / EModulus;
        strain(1) = ( -nu * values [ 0 ] + values [ 1 ] - nu * values [ 2 ] ) / EModulus;
        strain(2) = ( -nu * values [ 0 ] - nu * values [ 1 ] + values [ 2 ] ) / EModulus;
        strain(3) = 2. * ( 1. + nu ) * values [ 3 ] / EModulus;
    } else if ( myMode == _3dRotContinuum ) {
        // Order: r,  theta,  z,  zr
        strain(0) = ( values [ 0 ] - nu * values [ 1 ] - nu * values [ 2 ] ) / EModulus;
        strain(1) = ( -nu * values [ 0 ] + values [ 1 ] - nu * values [ 2 ] ) / EModulus;
        strain(2) = ( -nu * values [ 0 ] - nu * values [ 1 ] + values [ 2 ] ) / EModulus;
        strain(3) = 2. * ( 1. + nu ) * values [ 3 ] / EModulus;
    } else {
        strain(0) = ( values [ 0 ] - nu * values [ 1 ] - nu * values [ 2 ] ) / EModulus;
        strain(1) = ( -nu * values [ 0 ] + values [ 1 ] - nu * values [ 2 ] ) / EModulus;
        strain(2) = ( -nu * values [ 0 ] - nu * values [ 1 ] + values [ 2 ] ) / EModulus;
        strain(3) = ( 2. * ( 1. + nu ) * values [ 3 ] ) / EModulus;
        strain(4) = ( 2. * ( 1. + nu ) * values [ 4 ] ) / EModulus;
        strain(5) = ( 2. * ( 1. + nu ) * values [ 5 ] ) / EModulus;
    }
}

void
StressVector :: applyDeviatoricElasticCompliance(StrainVector &strain,
                                                 const double EModulus,
                                                 const double nu) const
{
    //
    // This function applies the elastic compliance to the deviatoric strain vector
    //
    applyDeviatoricElasticCompliance( strain, EModulus / 2. / ( 1. + nu ) );
}

void
StressVector :: applyDeviatoricElasticCompliance(StrainVector &strain,
                                                 const double GModulus) const
{
    //
    // This function applies the elastic compliance to the deviatoric strain vector
    //
    MaterialMode myMode = giveStressStrainMode();
    if ( myMode == _1dMat ) {
        OOFEM_ERROR("StressVector::applyDeviatoricElasticCompliance: No Split for 1D");
    } else if ( myMode == _PlaneStress ) {
        OOFEM_ERROR("StressVector::applyDeviatoricElasticCompliance: No Split for Plane Stress");
    } else if ( myMode == _PlaneStrain ) {
        strain(0) = 1. / ( 2. * GModulus ) * values [ 0 ];
        strain(1) = 1. / ( 2. * GModulus ) * values [ 1 ];
        strain(2) = 1. / ( 2. * GModulus ) * values [ 2 ];
        strain(3) = 1. / GModulus * values [ 3 ];
    } else if ( myMode == _PlaneStrainGrad ) {
        strain(0) = 1. / ( 2. * GModulus ) * values [ 0 ];
        strain(1) = 1. / ( 2. * GModulus ) * values [ 1 ];
        strain(2) = 1. / ( 2. * GModulus ) * values [ 2 ];
        strain(3) = 1. / GModulus * values [ 3 ];
        strain(4) = values[4];
    } else if ( myMode == _3dRotContinuum ) {
        // Order: r,  theta,  z,  zr
        strain(0) = 1. / ( 2. * GModulus ) * values [ 0 ];
        strain(1) = 1. / ( 2. * GModulus ) * values [ 1 ];
        strain(2) = 1. / ( 2. * GModulus ) * values [ 2 ];
        strain(3) = 1. / GModulus * values [ 3 ];
    } else {
        strain(0) = 1. / ( 2. * GModulus ) * values [ 0 ];
        strain(1) = 1. / ( 2. * GModulus ) * values [ 1 ];
        strain(2) = 1. / ( 2. * GModulus ) * values [ 2 ];
        strain(3) = 1. / GModulus * values [ 3 ];
        strain(4) = 1. / GModulus * values [ 4 ];
        strain(5) = 1. / GModulus * values [ 5 ];
    }
}

double
StressVector :: computeStressNorm() const
{
    //
    // This function computes the tensorial Norm of the stress in engineering notation
    //
    MaterialMode myMode = giveStressStrainMode();
    if ( myMode == _1dMat ) {
        // 1d problem
        return abs(values [ 0 ]);
    } else if ( myMode == _PlaneStress ) {
        // 2d problem: plane stress
        return sqrt(values [ 0 ] * values [ 0 ] + values [ 1 ] * values [ 1 ] + 2. * values [ 2 ] * values [ 2 ]);
    } else if ( myMode == _PlaneStrain || myMode == _3dRotContinuum ) {
        //  plane strain or axisymmetry
        return sqrt(values [ 0 ] * values [ 0 ] + values [ 1 ] * values [ 1 ] + values [ 2 ] * values [ 2 ] +
                    2. * values [ 3 ] * values [ 3 ]);
    } else if ( myMode == _PlaneStrainGrad) {
        //  plane strain or axisymmetry
        return sqrt(values [ 0 ] * values [ 0 ] + values [ 1 ] * values [ 1 ] + values [ 2 ] * values [ 2 ] +
                    2. * values [ 3 ] * values [ 3 ]);
    } else {
        // 3d problem
        return sqrt(values [ 0 ] * values [ 0 ] + values [ 1 ] * values [ 1 ] + values [ 2 ] * values [ 2 ] +
                    2. * values [ 3 ] * values [ 3 ] + 2. * values [ 4 ] * values [ 4 ] + 2. * values [ 5 ] * values [ 5 ]);
    }
}


void
StressVector :: giveTranformationMtrx(FloatMatrix &answer,
                                      const FloatMatrix &base,
                                      int transpose) const
//
// returns transformation matrix for 3d - stress to another system of axes,
// given by base.
// In base (FloatMatrix[3,3]) there are on each column stored vectors of
// coordinate system to which we do transformation.
//
// If transpose == 1 we transpose base matrix before transforming
//
{
    FloatMatrix t;
    answer.resize(6, 6);
    answer.zero();

    if ( transpose == 1 ) {
        t.beTranspositionOf(base);
    } else {
        t = base;
    }

    answer.at(1, 1) = t.at(1, 1) * t.at(1, 1);
    answer.at(1, 2) = t.at(2, 1) * t.at(2, 1);
    answer.at(1, 3) = t.at(3, 1) * t.at(3, 1);
    answer.at(1, 4) = 2.0 * t.at(2, 1) * t.at(3, 1);
    answer.at(1, 5) = 2.0 * t.at(1, 1) * t.at(3, 1);
    answer.at(1, 6) = 2.0 * t.at(1, 1) * t.at(2, 1);

    answer.at(2, 1) = t.at(1, 2) * t.at(1, 2);
    answer.at(2, 2) = t.at(2, 2) * t.at(2, 2);
    answer.at(2, 3) = t.at(3, 2) * t.at(3, 2);
    answer.at(2, 4) = 2.0 * t.at(2, 2) * t.at(3, 2);
    answer.at(2, 5) = 2.0 * t.at(1, 2) * t.at(3, 2);
    answer.at(2, 6) = 2.0 * t.at(1, 2) * t.at(2, 2);

    answer.at(3, 1) = t.at(1, 3) * t.at(1, 3);
    answer.at(3, 2) = t.at(2, 3) * t.at(2, 3);
    answer.at(3, 3) = t.at(3, 3) * t.at(3, 3);
    answer.at(3, 4) = 2.0 * t.at(2, 3) * t.at(3, 3);
    answer.at(3, 5) = 2.0 * t.at(1, 3) * t.at(3, 3);
    answer.at(3, 6) = 2.0 * t.at(1, 3) * t.at(2, 3);

    answer.at(4, 1) = t.at(1, 2) * t.at(1, 3);
    answer.at(4, 2) = t.at(2, 2) * t.at(2, 3);
    answer.at(4, 3) = t.at(3, 2) * t.at(3, 3);
    answer.at(4, 4) = ( t.at(2, 2) * t.at(3, 3) + t.at(3, 2) * t.at(2, 3) );
    answer.at(4, 5) = ( t.at(1, 2) * t.at(3, 3) + t.at(3, 2) * t.at(1, 3) );
    answer.at(4, 6) = ( t.at(1, 2) * t.at(2, 3) + t.at(2, 2) * t.at(1, 3) );

    answer.at(5, 1) = t.at(1, 1) * t.at(1, 3);
    answer.at(5, 2) = t.at(2, 1) * t.at(2, 3);
    answer.at(5, 3) = t.at(3, 1) * t.at(3, 3);
    answer.at(5, 4) = ( t.at(2, 1) * t.at(3, 3) + t.at(3, 1) * t.at(2, 3) );
    answer.at(5, 5) = ( t.at(1, 1) * t.at(3, 3) + t.at(3, 1) * t.at(1, 3) );
    answer.at(5, 6) = ( t.at(1, 1) * t.at(2, 3) + t.at(2, 1) * t.at(1, 3) );

    answer.at(6, 1) = t.at(1, 1) * t.at(1, 2);
    answer.at(6, 2) = t.at(2, 1) * t.at(2, 2);
    answer.at(6, 3) = t.at(3, 1) * t.at(3, 2);
    answer.at(6, 4) = ( t.at(2, 1) * t.at(3, 2) + t.at(3, 1) * t.at(2, 2) );
    answer.at(6, 5) = ( t.at(1, 1) * t.at(3, 2) + t.at(3, 1) * t.at(1, 2) );
    answer.at(6, 6) = ( t.at(1, 1) * t.at(2, 2) + t.at(2, 1) * t.at(1, 2) );
}
} // end namespace oofem
