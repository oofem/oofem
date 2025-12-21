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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "index.h"
#include "tensor1.h"
#include "tensor4.h"
#pragma once

using namespace FTensor;



/**
 * Functions providing basic operations with tensors
 * based on FTensor library
 * It provides transformation from the second and fourth order tensors in 3d to FloatArrayF<9> and FloatArray<9,9>
 * It provides calculation of determinant, cofactor, inverse, and tensor cross product
 * More to be added gradually
 * @author Martin Horak
 **/

namespace oofem {
class Tensor2_3d : public Tensor2< double, 3, 3 >
{
public:
    using Tensor2< double, 3, 3 >::Tensor2;

    /**
     * Creates second-order tensor in 3d
     * Initialized with zeros
     */

    Tensor2_3d() {
        this->data [ 0 ] [ 0 ] = 0;
        this->data [ 0 ] [ 1 ] = 0.;
        this->data [ 0 ] [ 2 ] = 0.;
        this->data [ 1 ] [ 0 ] = 0.;
        this->data [ 1 ] [ 1 ] = 0.;
        this->data [ 1 ] [ 2 ] = 0.;
        this->data [ 2 ] [ 0 ] = 0.;
        this->data [ 2 ] [ 1 ] = 0.;
        this->data [ 2 ] [ 2 ] = 0.;
    }

    /**
     * Creates a second-order order tensor in 3d from floatarrayf<9>
     */

    Tensor2_3d(const oofem::FloatArrayF< 9 > &array) {
        this->data [ 0 ] [ 0 ] = array.at(1);
        this->data [ 0 ] [ 1 ] = array.at(6);
        this->data [ 0 ] [ 2 ] = array.at(5);
        this->data [ 1 ] [ 0 ] = array.at(9);
        this->data [ 1 ] [ 1 ] = array.at(2);
        this->data [ 1 ] [ 2 ] = array.at(4);
        this->data [ 2 ] [ 0 ] = array.at(8);
        this->data [ 2 ] [ 1 ] = array.at(7);
        this->data [ 2 ] [ 2 ] = array.at(3);
    }

    /**
     * Transforms a second-order tensor into a floatarrayf<9>,  using the Voigt notation
     */

    const inline FloatArrayF< 9 >to_voigt_form()
    {
        return {
                   this->operator()(0, 0),
                   this->operator()(1, 1),
                   this->operator()(2, 2),
                   this->operator()(1, 2),
                   this->operator()(0, 2),
                   this->operator()(0, 1),
                   this->operator()(2, 1),
                   this->operator()(2, 0),
                   this->operator()(1, 0)
        };
    }


    /**
     * Computes power of a second order tensor
     * @param pow power
     * @param nf tolerance for eigenvalue calculation
     * @return Power of a Second-order tensor
     */
    Tensor2_3d computeTensorPower(double pow, int nf = 10)
    {
        auto [ eVals, eVecs ] = this->eigs(nf);
        return computeTensorPowerFromEigs(eVals, eVecs, pow);
    }

    /**
     * Computes power of a second order tensor from eigenvalues and eigenvectors
     * @param eVals eigen values
     * @param eVecs eigenvectors
     * @param m power
     * @return Power of a Second-order tensor
     */

    static Tensor2_3d computeTensorPowerFromEigs(const FloatArray &eVals, const FloatMatrix &eVecs, double m)
    {
        Tensor2_3d Cpow;
        for ( int i = 0; i < 3; i++ ) {
            for ( int j = 0; j < 3; j++ ) {
                Cpow(i, j) = pow(eVals.at(1), m) * eVecs(i, 0) * eVecs(j, 0) + pow(eVals.at(2), m) * eVecs(i, 1) * eVecs(j, 1) + pow(eVals.at(3), m) * eVecs(i, 2) * eVecs(j, 2);
            }
        }
        return Cpow;
    }

    /**
     * Computes power of a second order tensor from eigenvalues and eigenvectors and given coeffcinets for each eigenvalue
     * @param eVals eigen values
     * @param eVecs eigenvectors
     * @param m power
     * @param coeff coefficients
     * @return Power of a Second-order tensor
     */
    static Tensor2_3d computeTensorPowerFromEigs(const FloatArray &eVals, const FloatMatrix &eVecs, double m, const FloatArray &coeff)
    {
        Tensor2_3d Cpow;
        for ( int i = 0; i < 3; i++ ) {
            for ( int j = 0; j < 3; j++ ) {
                Cpow(i, j) -= coeff.at(1) * pow(eVals.at(1), m) * eVecs.at(i, 0) * eVecs.at(j, 0) + coeff.at(2) * pow(eVals.at(2), m) * eVecs.at(i, 1) * eVecs.at(j, 1) + coeff.at(3) * pow(eVals.at(3), m) * eVecs.at(i, 2) * eVecs.at(j, 2);
            }
        }
        return Cpow;
    }



    /**
     * Solves the eigenvalues and eigenvectors of symmetric second-order tensor
     *  adapted from FloatMatrix function jacobi_
     *
     * @param nf - number of significant figures
     * @return pair of eigenvalues and eigenvectors
     */
    std::pair< FloatArrayF< 3 >, FloatMatrixF< 3, 3 > >eigs(int nf = 15) const
    {
        Tensor2_3d copy(* this);
        // Local variables
        int neq = 3;
        double c_b2 = .10;

        FloatArrayF< 3 >eval;
        FloatMatrixF< 3, 3 >v;

        for ( int i = 1; i <= neq; i++ ) {
            eval.at(i) = copy(i - 1, i - 1);
        }

        double tol = pow(c_b2, nf);
        double sum = 0.0;
        for ( int i = 1; i <= neq; ++i ) {
            for ( int j = 1; j <= neq; ++j ) {
                sum += fabs(copy(i - 1, j - 1) );
                v.at(i, j) = 0.0;
            }

            v.at(i, i) = 1.0;
        }

        if ( sum > 0.0 ) {
            // ---- REDUCE MATRIX TO DIAGONAL ----------------
            int ite = 0;
            double ssum;
            do {
                ssum = 0.0;
                for ( int j = 2; j <= neq; ++j ) {
                    int ih = j - 1;
                    for ( int i = 1; i <= ih; ++i ) {
                        if ( ( fabs(copy(i - 1, j - 1) ) / sum ) > tol ) {
                            ssum += fabs(copy(i - 1, j - 1) );
                            // ---- CALCULATE ROTATION ANGLE -----------------
                            double aa = atan2(copy(i - 1, j - 1) * 2.0, eval.at(i) - eval.at(j) ) /  2.0;
                            double si = sin(aa);
                            double co = cos(aa);
                            // ---- MODIFY "I" AND "J" COLUMNS OF "A" AND "V"
                            for ( int k = 1; k < i; ++k ) {
                                double tt = copy(k - 1, i - 1);
                                copy(k - 1, i - 1) = co * tt + si * copy(k - 1, j - 1);
                                copy(k - 1, j - 1) = -si * tt + co * copy(k - 1, j - 1);
                                tt = v.at(k, i);
                                v.at(k, i) = co * tt + si * v.at(k, j);
                                v.at(k, j) = -si * tt + co * v.at(k, j);
                            }

                            // diagonal term (i,i)
                            double tt = eval.at(i);
                            eval.at(i) = co * tt + si * copy(i - 1, j - 1);
                            double aij = -si * tt + co * copy(i - 1, j - 1);
                            tt = v.at(i, i);
                            v.at(i, i) = co * tt + si * v.at(i, j);
                            v.at(i, j) = -si * tt + co * v.at(i, j);

                            for ( int k = i + 1; k < j; ++k ) {
                                double tt = copy(i - 1, k - 1);
                                copy(i - 1, k - 1) = co * tt + si * copy(k - 1, j - 1);
                                copy(k - 1, j - 1) = -si * tt + co * copy(k - 1, j - 1);
                                tt = v.at(k, i);
                                v.at(k, i) = co * tt + si * v.at(k, j);
                                v.at(k, j) = -si * tt + co * v.at(k, j);
                            }

                            // diagonal term (j,j)
                            tt = copy(i - 1, j - 1);
                            double aji = co * tt + si * eval.at(j);
                            eval.at(j) = -si * tt + co * eval.at(j);

                            tt = v.at(j, i);
                            v.at(j, i) = co * tt + si * v.at(j, j);
                            v.at(j, j) = -si * tt + co * v.at(j, j);
                            //
                            for ( int k = j + 1; k <= neq; ++k ) {
                                double tt = copy(i - 1, k - 1);
                                copy(i - 1, k - 1) = co * tt + si * copy(j - 1, k - 1);
                                copy(j - 1, k - 1) = -si * tt + co * copy(j - 1, k - 1);
                                tt = v.at(k, i);
                                v.at(k, i) = co * tt + si * v.at(k, j);
                                v.at(k, j) = -si * tt + co * v.at(k, j);
                            }

                            // ---- MODIFY DIAGONAL TERMS --------------------
                            eval.at(i) = co * eval.at(i) + si * aji;
                            eval.at(j) = -si * aij + co * eval.at(j);
                            copy(i - 1, j - 1) = 0.0;
                        } else {
                            // ---- A(I,J) MADE ZERO BY ROTATION -------------
                            ;
                        }
                    }
                }

                // ---- CHECK FOR CONVERGENCE --------------------
                if ( ++ite > 50 ) {
                    OOFEM_ERROR("too many iterations");
                }
            } while ( fabs(ssum) / sum > tol );
        }

        return { eval, v };
    }


    /**
     * Calculated derivative of C^m wrt C
     * @param m - power
     * @param nf - number of significant figures
     * @return Second-order tensor dC^m_dC
     */
    Tensor4_3d  compute_dCm_dC(double m, int nf = 10)
    {
        auto [ eVals, eVecs ] = this->eigs(nf);
        return compute_dCm_dC_fromEigs(m, eVals, eVecs);
    }


    /**
     * Calculated derivative of C^m wrt C from given eigenvalues, eigenvectors and power m
     * @param m - power
     * @param nf - number of significant figures
     * @return Second-order tensor dC^m_dC
     */
    static Tensor4_3d compute_dCm_dC_fromEigs(double m, const FloatArray &lam, const FloatMatrix &N)
    {
        FloatArray d(3), c(3);

        c.at(1) = pow(lam.at(1),  m);
        c.at(2) = pow(lam.at(2),  m);
        c.at(3) = pow(lam.at(3),  m);

        d.at(1) = m * pow(lam.at(1), m - 1.);
        d.at(2) = m * pow(lam.at(2), m - 1.);
        d.at(3) = m * pow(lam.at(3), m - 1.);

        FloatMatrix theta(3, 3);


        // compute auxiliary variables
        // the computation differes depends on if the eigenvalues of C are equal or not
        if ( lam.at(1) != lam.at(2) ) {
            if ( lam.at(2) != lam.at(3) ) {
                if ( lam.at(1) != lam.at(3) ) {
                    // all eigenvalues are different
                    for ( int i = 1; i <= 3; i++ ) {
                        for ( int j = 1; j <= 3; j++ ) {
                            if ( i == j ) {
                                continue;
                            } else {
                                theta.at(i, j) = ( c.at(i) - c.at(j) ) / ( lam.at(i) - lam.at(j) );
                            }
                        }
                    }
                } else { //l1 == l3 && l1 != l2
                    for ( int i = 1; i <= 3; i++ ) {
                        for ( int j = 1; j <= 3; j++ ) {
                            if ( i == j ) {
                                continue;
                            } else {
                                if ( ( i == 1 && j == 3 ) || ( i == 3 && j == 1 ) ) {
                                    theta.at(i, j) = 1. / 2. * d.at(i);
                                } else {
                                    theta.at(i, j) = ( c.at(i) - c.at(j) ) / ( lam.at(i) - lam.at(j) );
                                }
                            }
                        }
                    }
                }
            } else { //l2 == l3 && l1 != l2
                for ( int i = 1; i <= 3; i++ ) {
                    for ( int j = 1; j <= 3; j++ ) {
                        if ( i == j ) {
                            continue;
                        } else {
                            if ( ( i == 2 && j == 3 ) || ( i == 3 && j == 2 ) ) {
                                theta.at(i, j) = 1. / 2. * d.at(i);
                            } else {
                                theta.at(i, j) = ( c.at(i) - c.at(j) ) / ( lam.at(i) - lam.at(j) );
                            }
                        }
                    }
                }
            }
        } else if ( lam.at(1) != lam.at(3) ) { // l1 == l2  && l1 != l3
            for ( int i = 1; i <= 3; i++ ) {
                for ( int j = 1; j <= 3; j++ ) {
                    if ( i == j ) {
                        continue;
                    } else {
                        if ( ( i == 1 && j == 2 ) || ( i == 2 && j == 1 ) ) {
                            theta.at(i, j) = 1. / 2. * d.at(i);
                        } else {
                            theta.at(i, j) = ( c.at(i) - c.at(j) ) / ( lam.at(i) - lam.at(j) );
                        }
                    }
                }
            }
        } else { // l1 == l2 == l3
            for ( int i = 1; i <= 3; i++ ) {
                for ( int j = 1; j <= 3; j++ ) {
                    theta.at(i, j) = 1. / 2. * d.at(i);
                }
            }
        }


        FloatMatrix M(9, 9);
        std::vector< std::vector< int > >vIindex = {{ 1, 6, 5 }, { 9, 2, 4 }, { 8, 7, 3 } };
        for ( int i = 1; i <= 3; i++ ) {
            for ( int j = 1; j <= 3; j++ ) {
                for ( int k = 1; k <= 3; k++ ) {
                    for ( int l = 1; l <= 3; l++ ) {
                        M.at(vIindex [ i - 1 ] [ j - 1 ], vIindex [ k - 1 ] [ l - 1 ]) = N.at(k, i) * N.at(l, j) + N.at(k, j) * N.at(l, i);
                    }
                }
            }
        }

        Tensor4_3d dCm_dC;
        for ( int k = 1; k <= 3; k++ ) {
            for ( int l = 1; l <= 3; l++ ) {
                for ( int m = 1; m <= 3; m++ ) {
                    for ( int n = 1; n <= 3; n++ ) {
                        for ( int i = 1; i <= 3; i++ ) {
                            dCm_dC(k - 1, l - 1, m - 1, n - 1) += 0.5 * d.at(i) * N.at(k, i) * N.at(l, i) * M.at(vIindex [ i - 1 ] [ i - 1 ], vIindex [ m - 1 ] [ n - 1 ]);
                            for ( int j  = 1; j <= 3; j++ ) {
                                if ( j != i ) {
                                    dCm_dC(k - 1, l - 1, m - 1, n - 1) += 0.5 * theta.at(i, j) * N.at(k, i) * N.at(l, j) * M.at(vIindex [ i - 1 ] [ j - 1 ], vIindex [ m - 1 ] [ n - 1 ]);
                                }
                            }
                        }
                    }
                }
            }
        }
        return dCm_dC;
    }





    /**
     * Creates second-order unit tensor
     * @return Second-order unit tensor
     */
    static const Tensor2_3d UnitTensor()
    {
        return { 1.0, 0., 0., 0., 1.0, 0., 0., 0., 1.0 };
    }


    /**
     * Computes cofactor
     * @return Second-order cofactor tensor
     */
    inline Tensor2_3d compute_cofactor() const
    {
        Tensor2_3d cofF;
        cofF(i_3, j_3) = 0.5 * this->compute_tensor_cross_product(* this)(i_3, j_3);
        return cofF;
    }

    /**
     * Computes determinant
     * @return determinant
     */
    inline double compute_determinant() const
    {
        return ( 1. / 6. *  this->compute_tensor_cross_product(* this)(m_3, n_3) * this->operator()(m_3, n_3) );
    }

    /**
     * Computes determinant and cofactor
     * @return determinant and cofactor
     */
    inline std::pair< double, Tensor2_3d >compute_determinant_and_cofactor() const
    {
        Tensor2_3d cofF;
        cofF(i_3, j_3) =  0.5 * this->compute_tensor_cross_product(* this)(i_3, j_3);
        auto detF =  1. / 3. * cofF(i_3, j_3) * this->operator()(i_3, j_3);
        return { detF, cofF };
    }

    /**
     * Computes inverse
     * @return inverse second-order tensor
     */
    inline Tensor2_3d compute_inverse() const
    {
        Tensor2_3d iF;
        auto [ J, cofF ] = this->compute_determinant_and_cofactor();
        iF(i_3, j_3) = 1. / J * cofF(j_3, i_3);
        return iF;
    }

    /**
     * Computes second-order tensor cross product
     * @param tensor to calculate tensor cross product with
     * @return tensor cross product, eps_{ikm}eps_{jln}this_{kl}B_{mn}
     */
    inline Tensor2_3d compute_tensor_cross_product(const Tensor2_3d &B) const
    {
        Tensor2_3d C(this->operator()(1, 1) * B(2, 2) - this->operator()(1, 2) * B(2, 1) - this->operator()(2, 1) * B(1, 2) + this->operator()(2, 2) * B(1, 1), this->operator()(1, 2) * B(2, 0) - this->operator()(1, 0) * B(2, 2) + this->operator()(2, 0) * B(1, 2) - this->operator()(2, 2) * B(1, 0), this->operator()(1, 0) * B(2, 1) - this->operator()(1, 1) * B(2, 0) - this->operator()(2, 0) * B(1, 1) + this->operator()(2, 1) * B(1, 0), this->operator()(0, 2) * B(2, 1) - this->operator()(0, 1) * B(2, 2) + this->operator()(2, 1) * B(0, 2) - this->operator()(2, 2) * B(0, 1), this->operator()(0, 0) * B(2, 2) - this->operator()(0, 2) * B(2, 0) - this->operator()(2, 0) * B(0, 2) + this->operator()(2, 2) * B(0, 0), this->operator()(0, 1) * B(2, 0) - this->operator()(0, 0) * B(2, 1) + this->operator()(2, 0) * B(0, 1) - this->operator()(2, 1) * B(0, 0), this->operator()(0, 1) * B(1, 2) - this->operator()(0, 2) * B(1, 1) - this->operator()(1, 1) * B(0, 2) + this->operator()(1, 2) * B(0, 1), this->operator()(0, 2) * B(1, 0) - this->operator()(0, 0) * B(1, 2) + this->operator()(1, 0) * B(0, 2) - this->operator()(1, 2) * B(0, 0), this->operator()(0, 0) * B(1, 1) - this->operator()(0, 1) * B(1, 0) - this->operator()(1, 0) * B(0, 1) + this->operator()(1, 1) * B(0, 0) );
        return C;
    }


    /**
     * Computes fourth-order tensor cross product
     * @return tensor cross product, eps_{ikm}eps_{jln}this_{mn}
     */
    inline Tensor4_3d compute_tensor_cross_product() const
    {
        Tensor4_3d Ax(0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);
        Ax(1, 1, 0, 0) =  this->operator()(2, 2);
        Ax(1, 2, 0, 0) = -this->operator()(2, 1);
        Ax(2, 1, 0, 0) = -this->operator()(1, 2);
        Ax(2, 2, 0, 0) =  this->operator()(1, 1);

        Ax(0, 1, 1, 0) = -this->operator()(2, 2);
        Ax(0, 2, 1, 0) =  this->operator()(2, 1);
        Ax(2, 1, 1, 0) =  this->operator()(0, 2);
        Ax(2, 2, 1, 0) = -this->operator()(0, 1);

        Ax(0, 1, 2, 0) =  this->operator()(1, 2);
        Ax(0, 2, 2, 0) = -this->operator()(1, 1);
        Ax(1, 1, 2, 0) = -this->operator()(0, 2);
        Ax(1, 2, 2, 0) =  this->operator()(0, 1);

        Ax(1, 0, 0, 1) = -this->operator()(2, 2);
        Ax(1, 2, 0, 1) =  this->operator()(2, 0);
        Ax(2, 0, 0, 1) =  this->operator()(1, 2);
        Ax(2, 2, 0, 1) = -this->operator()(1, 0);

        Ax(0, 0, 1, 1) =  this->operator()(2, 2);
        Ax(0, 2, 1, 1) = -this->operator()(2, 0);
        Ax(2, 0, 1, 1) = -this->operator()(0, 2);
        Ax(2, 2, 1, 1) =  this->operator()(0, 0);

        Ax(0, 0, 2, 1) = -this->operator()(1, 2);
        Ax(0, 2, 2, 1) =  this->operator()(1, 0);
        Ax(1, 0, 2, 1) =  this->operator()(0, 2);
        Ax(1, 2, 2, 1) = -this->operator()(0, 0);

        Ax(1, 0, 0, 2) =  this->operator()(2, 1);
        Ax(1, 1, 0, 2) = -this->operator()(2, 0);
        Ax(2, 0, 0, 2) = -this->operator()(1, 1);
        Ax(2, 1, 0, 2) =  this->operator()(1, 0);

        Ax(0, 0, 1, 2) = -this->operator()(2, 1);
        Ax(0, 1, 1, 2) =  this->operator()(2, 0);
        Ax(2, 0, 1, 2) =  this->operator()(0, 1);
        Ax(2, 1, 1, 2) = -this->operator()(0, 0);

        Ax(0, 0, 2, 2) =  this->operator()(1, 1);
        Ax(0, 1, 2, 2) = -this->operator()(1, 0);
        Ax(1, 0, 2, 2) = -this->operator()(0, 1);
        Ax(1, 1, 2, 2) =  this->operator()(0, 0);

        return Ax;
    }
};

// not used for now...
class Tensor2sym_3d : public Tensor2_symmetric< double, 3 >
{
public:
    using Tensor2_symmetric< double, 3 >::Tensor2_symmetric;
    /*    Tensor2sym_3d(const oofem::FloatArrayF<6> &array){
     * this->data[0] = array.at(1);
     * this->data[1] = array.at(6);
     * this->data[2] = array.at(5);
     * this->data[3] = array.at(6);
     * this->data[4] = array.at(2);
     * this->data[5] = array.at(4);
     * this->data[6] = array.at(5);
     * this->data[7] = array.at(4);
     * this->data[8] = array.at(3);
     * }
     */


    const inline FloatArrayF< 6 >to_voigt_form()
    {
        return {
                   this->operator()(0, 0),
                   this->operator()(1, 1),
                   this->operator()(2, 2),
                   this->operator()(1, 2),
                   this->operator()(0, 2),
                   this->operator()(0, 1)
        };
    }
};
} // end namespace oofem
