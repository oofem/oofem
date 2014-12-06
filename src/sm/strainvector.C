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

#include "strainvector.h"
#include "stressvector.h"
#include "mathfem.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "error.h"

namespace oofem {
StrainVector :: StrainVector(MaterialMode m) : StressStrainBaseVector(m)
{ }

StrainVector :: StrainVector(const FloatArray &src, MaterialMode m) : StressStrainBaseVector(src, m)
{ }

void
StrainVector :: computeDeviatoricVolumetricSplit(StrainVector &dev, double &vol) const
{
    MaterialMode myMode = this->giveStressStrainMode();

    if ( myMode == _1dMat ) {
        // 1D model
        OOFEM_ERROR("No Split for 1D!");
        //  dev.resize(1); dev.at(1) = 0.0;
        // vol = this->at (1);
    } else if ( myMode == _PlaneStress ) {
        // plane stress problem
        OOFEM_ERROR("No Split for plane stress!");
        //  dev = *this;
        // vol = (this->at(1)+this->at(2))/2.0;
        //    dev.at(1) -= vol;
        // dev.at(2) -= vol;
    } else {
        // 3d problem or plane strain problem
        dev = * this;
        vol = ( this->at(1) + this->at(2) + this->at(3) ) / 3.0;
        dev.at(1) -= vol;
        dev.at(2) -= vol;
        dev.at(3) -= vol;
    }
}

void
StrainVector :: computeDeviatoricVolumetricSum(StrainVector &answer, const double vol) const
{
    MaterialMode myMode = this->giveStressStrainMode();

    if ( myMode == _1dMat ) {
        // 1D model
        OOFEM_ERROR("No sum for 1D!");
        //  dev.resize(1); dev.at(1) = 0.0;
        // vol = this->at (1);
    } else if ( myMode == _PlaneStress ) {
        // plane stress problem
        OOFEM_ERROR("No sum for plane stress!");
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
StrainVector :: computePrincipalValues(FloatArray &answer) const
{
    //
    // This function computes sorted principal values of a strain vector.
    // Engineering notation is used.
    //

    MaterialMode myMode = this->giveStressStrainMode();
    int size = this->giveSize();

    if ( myMode == _1dMat ) {
        answer = * this;
    } else if ( myMode == _PlaneStress ) {
        // 2D problem
        double ast, dst, D;
        answer.resize(2);

        ast = this->at(1) + this->at(2);
        dst = this->at(1) - this->at(2);
        D = sqrt( dst * dst + this->at(3) * this->at(3) );
        answer.at(1) = 0.5 * ( ast + D );
        answer.at(2) = 0.5 * ( ast - D );
    } else {
        // 3D problem
        double I1 = 0.0, I2 = 0.0, I3 = 0.0, s1, s2, s3;
        FloatArray s;
        int nonzeroFlag = 0;

        this->convertToFullForm(s);
        for ( int i = 1; i <= size; i++ ) {
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
             0.25 * ( s.at(4) * s.at(4) + s.at(5) * s.at(5) + s.at(6) * s.at(6) );
        I3 = s.at(1) * s.at(2) * s.at(3) +
             0.25 * ( s.at(4) * s.at(5) * s.at(6) - s.at(1) * s.at(4) * s.at(4) -
                      s.at(2) * s.at(5) * s.at(5) - s.at(3) * s.at(6) * s.at(6) );

        /*
         * Call cubic3r to ensure, that all three real eigenvalues will be found, because we have symmetric tensor.
         * This allows to overcome various rounding errors when solving general cubic equation.
         */
        int num;
        cubic3r( ( double ) -1., I1, -I2, I3, & s1, & s2, & s3, &num);

        if (num > 0 ) {
            answer.at(1) = s1;
        }

        if (num > 1 ) {
            answer.at(2) = s2;
        }

        if (num > 2 ) {
            answer.at(3) = s3;
        }

        // sort results
        for ( int i = 1; i < 3; i++ ) {
            for ( int j = 1; j < 3; j++ ) {
                if ( answer.at(j + 1) > answer.at(j) ) {
                    std :: swap( answer.at(j + 1), answer.at(j) );
                }
            }
        }
    }
}

void
StrainVector :: computeMaxPrincipalDir(FloatArray &answer) const
//
// This function computes the principal direction of the receiver
// associated with the maximum principal value.
//
{
    FloatArray princval;
    FloatMatrix princdir;
    this->computePrincipalValDir(princval, princdir);
    int n = princval.giveSize();
    answer.resize(n);
    for ( int i = 1; i <= n; i++ ) {
        answer.at(i) = princdir.at(i, 1);
    }
}

void
StrainVector :: computePrincipalValDir(FloatArray &answer, FloatMatrix &dir) const
{
    //
    // This function computes principal values & directions of the receiver.
    //
    // Return Values:
    //
    // answer -> principal strains (ordered from largest to smallest)
    // dir -> principal strain directions
    //

    FloatMatrix ss;
    FloatArray sp;
    int nval;
    int nonzeroFlag = 0;
    int size = this->giveSize();
    MaterialMode myMode = this->giveStressStrainMode();

    if ( myMode == _1dMat ) {
        // 1D problem
        answer = * this;
        dir.resize(1, 1);
        dir.at(1, 1) = 1.0;
        return;
    } else if ( myMode == _PlaneStress ) {
        // 2D problem
        ss.resize(2, 2);
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

        ss.at(1, 1) = this->at(1);
        ss.at(2, 2) = this->at(2);

        ss.at(1, 2) = ss.at(2, 1) = 0.5 * this->at(3);
    } else {
        // 3D problem
        ss.resize(3, 3);
        FloatArray s;

        answer.resize(3);

        this->convertToFullForm(s);
        for ( int i = 1; i <= size; i++ ) {
            if ( fabs( s.at(i) ) > 1.e-20 ) {
                nonzeroFlag = 1;
            }
        }

        if ( nonzeroFlag == 0 ) {
            answer.zero();
            return;
        }

        ss.at(1, 1) = s.at(1);
        ss.at(2, 2) = s.at(2);
        ss.at(3, 3) = s.at(3);
        ss.at(1, 2) = ss.at(2, 1) = 0.5 * s.at(6);
        ss.at(1, 3) = ss.at(3, 1) = 0.5 * s.at(5);
        ss.at(2, 3) = ss.at(3, 2) = 0.5 * s.at(4);
    }

#if 0
    ss.Jacobi(& answer, & dir, & i);
#else
    ss.jaco_(answer, dir, 10);
#endif
    // sort results (from the largest to the smallest eigenvalue)
    nval = 3;
    if ( myMode == _PlaneStress ) {
        nval = 2;
    }

    for ( int ii = 1; ii < nval; ii++ ) {
        for ( int jj = 1; jj < nval; jj++ ) {
            if ( answer.at(jj + 1) > answer.at(jj) ) {
                // swap eigenvalues and eigenvectors
                std :: swap( answer.at(jj + 1), answer.at(jj) );
                for ( int kk = 1; kk <= nval; kk++ ) {
                    std :: swap( dir.at(kk, jj + 1), dir.at(kk, jj) );
                }
            }
        }
    }
}

void
StrainVector :: printYourself() const
{
    printf("StrainVector (MaterialMode %d)\n", mode);
    for ( double x : this->values ) {
        printf("%10.3e  ", x);
    }

    printf("\n");
}


double
StrainVector :: computeVolumeChange() const
{
    //
    // This function computes the change of volume. This is different from the volumetric strain.
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
StrainVector :: computeStrainNorm() const
{
    //
    // This function computes the tensorial Norm of the strain in engineering notation
    //
    MaterialMode myMode = giveStressStrainMode();
    if ( myMode == _1dMat ) {
        // 1d problem
        return sqrt(values [ 0 ] * values [ 0 ]);
    } else if ( myMode == _PlaneStress ) {
        // 2d problem: plane stress
        return sqrt( .5 * ( 2. * values [ 0 ] * values [ 0 ] + 2 * values [ 1 ] * values [ 1 ] + values [ 2 ] * values [ 2 ] ) );
    } else if ( myMode == _PlaneStrain ) {
        //  plane strain or axisymmetry
        return sqrt( .5 * ( 2. * values [ 0 ] * values [ 0 ] + 2. * values [ 1 ] * values [ 1 ] + 2. * values [ 2 ] * values [ 2 ] +
                            values [ 3 ] * values [ 3 ] ) );
    } else {
        // 3d problem
        return sqrt( .5 * ( 2. * values [ 0 ] * values [ 0 ] + 2. * values [ 1 ] * values [ 1 ] + 2. * values [ 2 ] * values [ 2 ] +
                            values [ 3 ] * values [ 3 ] + values [ 4 ] * values [ 4 ] + values [ 5 ] * values [ 5 ] ) );
    }
}

void
StrainVector :: applyElasticStiffness(StressVector &stress, const double EModulus, const double nu) const
{
    //
    // This function applies the elastic stiffness to the total strain vector
    //

    MaterialMode myMode = giveStressStrainMode();
    if ( myMode == _1dMat ) {
        stress(0) = EModulus * values [ 0 ];
    } else if ( myMode == _PlaneStress ) {
        double factor = EModulus / ( 1. - nu * nu );
        stress(0) = factor * ( values [ 0 ] + nu * values [ 1 ] );
        stress(1) = factor * ( nu * values [ 0 ] + values [ 1 ] );
        stress(2) = factor * ( ( ( 1 - nu ) / 2. ) * values [ 2 ] );
    } else if ( myMode == _PlaneStrain ) {
        if ( nu >= 0.5 ) {
            OOFEM_ERROR("nu must be less than 0.5");
        }

        double factor = EModulus / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
        stress(0) = factor * ( ( 1. - nu ) * values [ 0 ] + nu * values [ 1 ] );
        stress(1) = factor * ( nu * values [ 0 ] + ( 1. - nu ) * values [ 1 ] );
        stress(2) = nu * ( stress(0) + stress(1) );
        stress(3) = factor * ( ( ( 1. - 2. * nu ) / 2. ) * values [ 3 ] );
    } else {
        if ( nu >= 0.5 ) {
            OOFEM_ERROR("nu must be less than 0.5");
        }

        double factor = EModulus / ( ( 1. + nu ) * ( 1. - 2. * nu ) );
        stress(0) = factor * ( ( 1. - nu ) * values [ 0 ] + nu * values [ 1 ] + nu * values [ 2 ] );
        stress(1) = factor * ( nu * values [ 0 ] + ( 1. - nu ) * values [ 1 ] + nu * values [ 2 ] );
        stress(2) = factor * ( nu * values [ 0 ] + nu * values [ 1 ] + ( 1. - nu ) * values [ 2 ] );
        stress(3) = factor * ( ( ( 1. - 2. * nu ) / 2. ) * values [ 3 ] );
        stress(4) = factor * ( ( ( 1. - 2. * nu ) / 2. ) * values [ 4 ] );
        stress(5) = factor * ( ( ( 1. - 2. * nu ) / 2. ) * values [ 5 ] );
    }
}

void
StrainVector :: applyDeviatoricElasticStiffness(StressVector &stress,
                                                const double EModulus,
                                                const double nu) const
{
    //
    // This function applies the elastic stiffness to the deviatoric strain vector
    //

    applyDeviatoricElasticStiffness( stress, EModulus / ( 2. * ( 1. + nu ) ) );
}

void
StrainVector :: applyDeviatoricElasticStiffness(StressVector &stress,
                                                const double GModulus) const
{
    //
    // This function applies the elastic stiffness to the deviatoric strain vector
    //

    MaterialMode myMode = giveStressStrainMode();
    if ( myMode == _1dMat ) {
        OOFEM_ERROR("No Split for 1D");
    } else if ( myMode == _PlaneStress ) {
        OOFEM_ERROR("No Split for Plane Stress");
    } else if ( myMode == _PlaneStrain ) {
        if ( stress.giveStressStrainMode() != _PlaneStrain ) {
            stress.letStressStrainModeBe(_PlaneStrain);
        }

        stress(0) = 2. * GModulus * values [ 0 ];
        stress(1) = 2. * GModulus * values [ 1 ];
        stress(2) = 2. * GModulus * values [ 2 ];
        stress(3) = GModulus * values [ 3 ];
    } else if ( myMode == _PlaneStrainGrad ) {
        if ( stress.giveStressStrainMode() != _PlaneStrainGrad ) {
            stress.letStressStrainModeBe(_PlaneStrainGrad);
        }

        stress(0) = 2. * GModulus * values [ 0 ];
        stress(1) = 2. * GModulus * values [ 1 ];
        stress(2) = 2. * GModulus * values [ 2 ];
        stress(3) = GModulus * values [ 3 ];
        stress(4) = values [ 4 ];
    } else {
        if ( stress.giveStressStrainMode() != _3dMat ) {
            stress.letStressStrainModeBe(_3dMat);
        }

        stress(0) = 2. * GModulus * values [ 0 ];
        stress(1) = 2. * GModulus * values [ 1 ];
        stress(2) = 2. * GModulus * values [ 2 ];
        stress(3) = GModulus * values [ 3 ];
        stress(4) = GModulus * values [ 4 ];
        stress(5) = GModulus * values [ 5 ];
    }
}



void
StrainVector :: giveTranformationMtrx(FloatMatrix &answer,
                                      const FloatMatrix &base,
                                      int transpose) const
//
// returns transformation matrix for 3d - strains to another system of axes,
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
    answer.at(1, 4) = t.at(2, 1) * t.at(3, 1);
    answer.at(1, 5) = t.at(1, 1) * t.at(3, 1);
    answer.at(1, 6) = t.at(1, 1) * t.at(2, 1);

    answer.at(2, 1) = t.at(1, 2) * t.at(1, 2);
    answer.at(2, 2) = t.at(2, 2) * t.at(2, 2);
    answer.at(2, 3) = t.at(3, 2) * t.at(3, 2);
    answer.at(2, 4) = t.at(2, 2) * t.at(3, 2);
    answer.at(2, 5) = t.at(1, 2) * t.at(3, 2);
    answer.at(2, 6) = t.at(1, 2) * t.at(2, 2);

    answer.at(3, 1) = t.at(1, 3) * t.at(1, 3);
    answer.at(3, 2) = t.at(2, 3) * t.at(2, 3);
    answer.at(3, 3) = t.at(3, 3) * t.at(3, 3);
    answer.at(3, 4) = t.at(2, 3) * t.at(3, 3);
    answer.at(3, 5) = t.at(1, 3) * t.at(3, 3);
    answer.at(3, 6) = t.at(1, 3) * t.at(2, 3);

    answer.at(4, 1) = 2.0 * t.at(1, 2) * t.at(1, 3);
    answer.at(4, 2) = 2.0 * t.at(2, 2) * t.at(2, 3);
    answer.at(4, 3) = 2.0 * t.at(3, 2) * t.at(3, 3);
    answer.at(4, 4) = ( t.at(2, 2) * t.at(3, 3) + t.at(3, 2) * t.at(2, 3) );
    answer.at(4, 5) = ( t.at(1, 2) * t.at(3, 3) + t.at(3, 2) * t.at(1, 3) );
    answer.at(4, 6) = ( t.at(1, 2) * t.at(2, 3) + t.at(2, 2) * t.at(1, 3) );

    answer.at(5, 1) = 2.0 * t.at(1, 1) * t.at(1, 3);
    answer.at(5, 2) = 2.0 * t.at(2, 1) * t.at(2, 3);
    answer.at(5, 3) = 2.0 * t.at(3, 1) * t.at(3, 3);
    answer.at(5, 4) = ( t.at(2, 1) * t.at(3, 3) + t.at(3, 1) * t.at(2, 3) );
    answer.at(5, 5) = ( t.at(1, 1) * t.at(3, 3) + t.at(3, 1) * t.at(1, 3) );
    answer.at(5, 6) = ( t.at(1, 1) * t.at(2, 3) + t.at(2, 1) * t.at(1, 3) );

    answer.at(6, 1) = 2.0 * t.at(1, 1) * t.at(1, 2);
    answer.at(6, 2) = 2.0 * t.at(2, 1) * t.at(2, 2);
    answer.at(6, 3) = 2.0 * t.at(3, 1) * t.at(3, 2);
    answer.at(6, 4) = ( t.at(2, 1) * t.at(3, 2) + t.at(3, 1) * t.at(2, 2) );
    answer.at(6, 5) = ( t.at(1, 1) * t.at(3, 2) + t.at(3, 1) * t.at(1, 2) );
    answer.at(6, 6) = ( t.at(1, 1) * t.at(2, 2) + t.at(2, 1) * t.at(1, 2) );
}
} // end namespace oofem
