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

#ifndef floatmatrixf_h
#define floatmatrixf_h

#include "oofemcfg.h"
#include "contextioresulttype.h"
#include "contextmode.h"
#include "floatmatrix.h"
#include "floatarrayf.h"

#include <array>
#include <initializer_list>
#include <algorithm>
#include <utility>
#include <iosfwd>

namespace oofem {

/**
 * Implementation of matrix containing floating point numbers.
 * @author Mikael Öhman
 */
template<int N, int M>
class OOFEM_EXPORT FloatMatrixF
{
protected:
    /// Values of matrix stored column wise.
    std::array< double, N*M > values;

public:
    /// @name Iterator for for-each loops:
    //@{
    auto begin() { return this->values.begin(); }
    auto end() { return this->values.end(); }
    auto begin() const { return this->values.begin(); }
    auto end() const { return this->values.end(); }
    //@}

    /**
     * Constructor (values are specified column-wise)
     * @note The syntax {{x,y,z},{...}} can be achieved by nested initializer_list, but 
     */
    template<typename... V> FloatMatrixF(V... x) : values{x...} { }
    /// Copy constructor.
    FloatMatrixF(const FloatMatrixF<N,M> &mat) : values(mat.values) {}

    /// Assignment operator
    FloatMatrixF &operator=(const FloatMatrixF<N,M> &mat)
    {
        values = mat.values;
        return * this;
    }

    /**
     * Checks size of receiver towards requested bounds.
     * @param i Required number of rows.
     * @param j Required number of columns.
     */
    void checkBounds(int i, int j) const
    {
        if ( i <= 0 ) {
            OOFEM_ERROR("matrix error on rows : %d < 0", i);
        }
        if ( j <= 0 ) {
            OOFEM_ERROR("matrix error on columns : %d < 0", j);
        }
        if ( i > N ) {
            OOFEM_ERROR("matrix error on rows : %d > %d", i, M);
        }
        if ( j > M ) {
            OOFEM_ERROR("matrix error on columns : %d > %d", j, M);
        }
    }

    /// Returns number of rows of receiver.
    int rows() const { return N; }
    /// Returns number of columns of receiver.
    int cols() const { return M; }

    /**
     * Coefficient access function. Returns value of coefficient at given
     * position of the receiver. Implements 1-based indexing.
     * @param i Row position of coefficient.
     * @param j Column position of coefficient.
     */
    double at(int i, int j) const
    {
#ifndef NDEBUG
        this->checkBounds(i, j);
#endif
        return values [ ( j - 1 ) * N + i - 1 ];
    }
    /**
     * Coefficient access function. Returns value of coefficient at given
     * position of the receiver. Implements 1-based indexing.
     * @param i Row position of coefficient.
     * @param j Column position of coefficient.
     */
    inline double &at(int i, int j)
    {
#ifndef NDEBUG
        this->checkBounds(i, j);
#endif
        return values [ ( j - 1 ) * N + i - 1 ];
    }

    /**
     * Direct value access (column major). Implements 0-based indexing.
     * @param i Position in data.
     */
    double &operator[](int i)
    {
        return values[ i ];
    }
    /**
     * Direct value access (column major). Implements 0-based indexing.
     * @param i Position in data.
     */
    double operator[](int i) const
    {
        return values[ i ];
    }

    /**
     * Coefficient access function. Returns l-value of coefficient at given
     * position of the receiver. Implements 0-based indexing.
     * @param i Row position of coefficient.
     * @param j Column position of coefficient.
     */
    double &operator()(int i, int j)
    {
#ifndef NDEBUG
        this->checkBounds(i + 1, j + 1);
#endif
        return values [ j * N + i ];
    }
    /**
     * Coefficient access function. Implements 0-based indexing.
     * @param i Row position of coefficient.
     * @param j Column position of coefficient.
     */
    double operator()(int i, int j) const
    {
#ifndef NDEBUG
        this->checkBounds(i + 1, j + 1);
#endif 
        return values [ j * N + i ];
    }

    /**
     * Sets the values of the matrix in specified column. If matrix size is zero, the size is adjusted.
     * @param src Array to set at column c.
     * @param c Column to set.
     */
    void setColumn(const FloatArrayF<N> &src, int c)
    {
        for ( int i = 0; i < N; i++ ) {
            (*this)(i, c) = src[i];
        }
    }
    /**
     * Adds to the receiver the product @f$ a^{\mathrm{T}}\cdot b \mathrm{d}V @f$. If the receiver has zero size, it is expanded.
     * Assumes that receiver and product @f$ a^{\mathrm{T}}\cdot b \mathrm{d}V @f$ are symmetric matrices. Computes only the
     * upper half of receiver.
     * @param a Matrix a in equation.
     * @param b Matrix b in equation.
     * @param dV Scaling factor.
     */
    template<int P>
    void plusProductSymmUpper(const FloatMatrixF<P,N> &a, const FloatMatrixF<P,M> &b, double dV)
    {
        ///@todo Larger matrices needs optimized BLAS.
        for ( int i = 0; i < N; ++i ) {
            for ( int j = i; j < M; ++j ) {
                double summ = 0.;
                for ( int k = 0; k < P; ++k ) {
                    summ += a(k, i) * b(k, j);
                }
                (*this)(i, j) += summ * dV;
            }
        }
    }
    /**
     * Adds to the receiver the dyadic product @f$ a \otimes a \mathrm{d}V @f$.
     * Computes only the upper half of receiver.
     * @param a Array a in equation.
     * @param dV Scaling factor.
     */
    void plusDyadSymmUpper(const FloatArrayF<N> &a, double dV)
    {
        ///@todo This method should only exist if we have N == M, how do I enforce this?
        for ( int i = 0; i < N; ++i ) {
            for ( int j = i; j < N; ++j ) {
                (*this)(i, j) = a[i] * a[j] * dV;
            }
        }
    }
    /**
     * Adds to the receiver the product @f$a^{\mathrm{T}} \cdot b \mathrm{d}V@f$.
     * @param a Matrix a in equation.
     * @param b Matrix b in equation.
     * @param dV Scaling factor.
     */
    template<int P>
    void plusProductUnsym(const FloatMatrixF<P,N> &a, const FloatMatrixF<P,M> &b, double dV)
    {
        ///@todo Larger matrices needs optimized BLAS.
        for ( int i = 0; i < N; ++i ) {
            for ( int j = 0; j < M; ++j ) {
                double summ = 0.;
                for ( int k = 0; k < P; ++k ) {
                    summ += a(k, i) * b(k, j);
                }
                (*this)(i, j) += summ * dV;
            }
        }
    }
    /**
     * Adds to the receiver the product @f$a \otimes b \mathrm{d}V@f$.
     * @param a Array a in equation.
     * @param b Array b in equation.
     * @param dV Scaling factor.
     */
    void plusDyadUnsym(const FloatArray &a, const FloatArray &b, double dV)
    {
        ///@todo This method should only exist if we have N == M, how do I enforce this?
        for ( int i = 0; i < N; ++i ) {
            for ( int j = 0; j < M; ++j ) {
                (*this)(i, j) = a[i] * b[j] * dV;
            }
        }
    }
    /**
     * Initializes the lower half of the receiver according to the upper half.
     */
    void symmetrized()
    {
        for ( int i = 2; i <= this->rows(); i++ ) {
            for ( int j = 1; j < i; j++ ) {
                this->at(i, j) = this->at(j, i);
            }
        }
    }

    /// Prints matrix to stdout. Useful for debugging.
    void printYourself(const std::string &name="FloatMatrixF") const
    {
        printf("%s (%d x %d): \n", name.c_str(), N, M);
        if ( this->rows() <= 250 && this->cols() <= 250 ) {
            for ( int i = 0; i < this->rows(); ++i ) {
                for ( int j = 0; j < this->cols() && j < 100; ++j ) {
                    printf( "%10.3e  ", (*this)(i, j) );
                }
                printf("\n");
            }
        } else {
            for ( int i = 0; i < this->rows() && i < 20; ++i ) {
                for ( int j = 0; j < this->cols() && j < 10; ++j ) {
                    printf( "%10.3e  ", (*this)(i, j) );
                }
                if ( this->cols() > 10 ) printf(" ...");
                printf("\n");
            }
            if ( this->cols() > 20 )  printf(" ...\n");
        }
    }

    /**
     * Exposes the internal values of the matrix. Should typically not be used outside of matrix classes.
     * @return Pointer to the values of the matrix.
     */
    const double *givePointer() const { return values.data(); }
    double *givePointer() { return values.data(); }

    contextIOResultType storeYourself(DataStream &stream) const
    {
        if ( !stream.write(values.data(), N*M) ) {
            return CIO_IOERR;
        }
        return CIO_OK;
    }

    contextIOResultType restoreYourself(DataStream &stream)
    {
        if ( !stream.read(values.data(), N*M) ) {
            return CIO_IOERR;
        }
        return CIO_OK;
    }

    int givePackSize(DataStream &buff) const
    {
        return buff.givePackSizeOfDouble(N*M);
    }
};

template<int N, int M>
std::ostream & operator << ( std::ostream & out, const FloatMatrixF<N,M> & x )
{
    out << N << " " << M << " {";
    for ( int i = 0; i < N; ++i ) {
        for ( int j = 0; j < M; ++j ) {
            out << " " << x(i, j);
        }
        out << ";";
    }
    out << "}";
    return out;
}

template<int N, int M>
FloatMatrixF<N,M> operator * ( double a, const FloatMatrixF<N,M> & x )
{
    FloatMatrixF<N,M> out;
    for ( int i = 0; i < N*M; ++i ) {
        out[i] = x[i] * a;
    }
    return out;
}

template<int N, int M>
FloatMatrixF<N,M> operator * ( const FloatMatrixF<N,M> & x, double a )
{
    return a*x;
}

template<int N, int M>
FloatMatrixF<N,M> operator / ( const FloatMatrixF<N,M> & x, double a )
{
    FloatMatrixF<N,M> out;
    for ( int i = 0; i < N*M; ++i ) {
        out[i] = x[i] / a;
    }
    return out;
}

template<int N, int M>
FloatMatrixF<N,M> operator + ( const FloatMatrixF<N,M> & x, const FloatMatrixF<N,M> & y )
{
    FloatMatrixF<N,M> out;
    for ( int i = 0; i < N*M; ++i ) {
        out[i] = x[i] + y[i];
    }
    return out;
}

template<int N, int M>
FloatMatrixF<N,M> operator - ( const FloatMatrixF<N,M> & x, const FloatMatrixF<N,M> & y )
{
    FloatMatrixF<N,M> out;
    for ( int i = 0; i < N*M; ++i ) {
        out[i] = x[i] - y[i];
    }
    return out;
}

template<int N, int M>
FloatMatrixF<N,M> &operator += ( FloatMatrixF<N,M> & x, const FloatMatrixF<N,M> & y )
{
    for ( int i = 0; i < N*M; ++i ) {
        x[i] += y[i];
    }
    return x;
}

template<int N, int M>
FloatMatrixF<N,M> &operator -= ( FloatMatrixF<N,M> & x, const FloatMatrixF<N,M> & y )
{
    for ( int i = 0; i < N*M; ++i ) {
        x[i] -= y[i];
    }
    return x;
}

template<int N, int M>
FloatMatrixF<N,M> &operator *= ( FloatMatrixF<N,M> & x, double a )
{
    for ( int i = 0; i < N*M; ++i ) {
        x[i] *= a;
    }
    return x;
}

/// Returns true if no element is not finite (NAN or infinite)
template<int N, int M>
bool isfinite(const FloatMatrixF<N,M> &mat)
{
    for ( auto &val : mat ) {
        if ( !std::isfinite(val) ) {
            return false;
        }
    }
    return true;
}

/**
 * Constructs the N matrix
 * @param n Vector with components which will appear in respective diagonal.
 */
template<int N, int NSD>
FloatMatrixF<N*3,NSD> Nmatrix(const FloatArrayF<N> &n)
{
    FloatMatrixF<N*3,NSD> out;
    for ( int i = 0; i < N; ++i ) {
        for ( int j = 0; j < NSD; ++j ) {
            out( j, i * NSD + j ) = n[j];
        }
    }
    return out;
}

/**
 * Constructs a local coordinary system for the given normal.
 * @param normal Normal (normalized).
 */
inline FloatMatrixF<2,2> local_cs(const FloatArrayF<2> &normal)
{
    return {
        normal[0], normal[1],
        normal[1], -normal[0]
    };
}

inline FloatMatrixF<3,3> local_cs(const FloatArrayF<3> &normal)
{
    // Create a permutated vector of n, *always* length 1 and significantly different from n.
    FloatArrayF<3> tangent = {normal[1], -normal[2], normal[0]};

    // Construct orthogonal vector
    double npn = dot(tangent, normal);
    tangent += -npn * normal;
    tangent /= norm(tangent);
    FloatArrayF<3> bigent;
    cross(tangent, normal);

    FloatMatrixF<3,3> out;
    out(0, 0) = normal[0];
    out(0, 1) = normal[1];
    out(0, 2) = normal[2];

    out(1, 0) = bigent[0];
    out(1, 1) = bigent[1];
    out(1, 2) = bigent[2];

    out(2, 0) = tangent[0];
    out(2, 1) = tangent[1];
    out(2, 2) = tangent[2];
    return out;
}

/// Constructs transposed matrix
template<int N, int M>
FloatMatrixF<M,N> transpose(const FloatMatrixF<N,M> &mat)
{
    FloatMatrixF<M,N> out;
    for ( int i = 0; i < N; ++i ) {
        for ( int j = 0; j < M; ++i ) {
            out(j, i) = mat(i, j);
        }
    }
    return out;
}

/// Computes @f$ a \cdot b @f$.
template<int N, int M, int P>
FloatMatrixF<N,P> dot(const FloatMatrixF<N,M> &a, const FloatMatrixF<M,P> &b)
{
    FloatMatrixF<N,P> out;
    ///@todo BLAS for larger matrix sizes (maybe)
    for ( int i = 0; i < N; i++ ) {
        for ( int j = 0; j < P; j++ ) {
            double coeff = 0.;
            for ( int k = 0; k < M; k++ ) {
                coeff += a(i, k) * b(k, j);
            }
            out(i, j) = coeff;
        }
    }
    return out;
}

/// Computes @f$ a \cdot b^{\mathrm{T}} @f$.
template<int N, int M, int P>
FloatMatrixF<N,P> dotT(const FloatMatrixF<N,M> &a, const FloatMatrixF<P,M> &b)
{
    FloatMatrixF<N,P> out;
    ///@todo BLAS for larger matrix sizes (maybe)
    for ( int i = 0; i < N; i++ ) {
        for ( int j = 0; j < P; j++ ) {
            double coeff = 0.;
            for ( int k = 0; k < M; k++ ) {
                coeff += a(i, k) * b(j, k);
            }
            out(i, j) = coeff;
        }
    }
    return out;
}

/// Computes @f$ a^{\mathrm{T}} \cdot b @f$.
template<int N, int M, int P>
FloatMatrixF<N,P> Tdot(const FloatMatrixF<M,N> &a, const FloatMatrixF<M,P> &b)
{
    FloatMatrixF<N,P> out;
    ///@todo BLAS for larger matrix sizes (maybe)
    for ( int i = 0; i < N; i++ ) {
        for ( int j = 0; j < P; j++ ) {
            double coeff = 0.;
            for ( int k = 0; k < M; k++ ) {
                coeff += a(k, i) * b(k, j);
            }
            out(i, j) = coeff;
        }
    }
    return out;
}

/// Computes @$ m_ij x_j = m \cdot x @$
template<int N, int M>
FloatArrayF<N> dot(const FloatMatrixF<N,M> &m, const FloatArrayF<M> &x)
{
    FloatArrayF<N> out;
    for ( int i = 0; i < N; i++ ) {
        double sum = 0.;
        for ( int j = 0; j < M; j++ ) {
            sum += m(i, j) * x[j];
        }
        out[i] = sum;
    }
    return out;
}

/// Computes @$ m_ji x_j = x \cdot m @$
template<int N, int M>
FloatArrayF<N> dot(const FloatMatrixF<M,N> &m, const FloatArrayF<M> &x)
{
    FloatArrayF<N> out;
    for ( int i = 0; i < M; i++ ) {
        double sum = 0.;
        for ( int j = 0; j < N; j++ ) {
            sum += x[j] * m(j, i);
        }
        out[i] = sum;
    }
    return out;
}

/// Computes the dyadic product @f$ m_{ij} = a_i b_j @f$.
template<int N, int M>
FloatMatrixF<N,M> dyad(const FloatArrayF<N> &a, const FloatArrayF<M> &b)
{
    FloatMatrixF<N,M> out;
    for ( int i = 0; i < N; ++i ) {
        for ( int j = 0; j < M; ++j ) {
            out(i, j) = a[i] * b[j];
        }
    }
    return out;
}

/// Computes @f$ a = r^{\mathrm{T}} \cdot a \cdot r @f$
template<int N, int M>
FloatMatrixF<M,M> rotate(FloatMatrixF<N,N> &a, const FloatMatrixF<N,M> &r)
{
    return dot(dotT(r, a), r);
}

/// Computes @f$ a = r \cdot a \cdot r^{\mathrm{T}} @f$
template<int N, int M>
FloatMatrixF<M,M> unrotate(FloatMatrixF<N,N> &a, const FloatMatrixF<M,N> &r)
{
    return dotT(dot(r, a), r);
}

/**
 * Returns a copy of the given matrix column,
 */
template<int N, int M>
FloatArrayF<N> column(const FloatMatrixF<N,M> &mat, int col)
{
    FloatArrayF<N> ans;
    for ( int row = 0; row < N; ++row ) {
        ans[row] = mat(row, col);
    }
    return ans;
}

/**
 * Symmetrizes and stores a matrix in Voigt form:
 * x_11, x_22, x_33, x_23, x_13, x_12, x_32, x_31, x_21
 */
inline FloatMatrixF<3,3> from_voigt_form(const FloatArrayF<9> &v)
{
    return {
        v[0], v[8], v[7], 
        v[5], v[1], v[6], 
        v[4], v[3], v[2]
    }; // note: column major order
}

/**
 * Symmetrizes and stores a matrix in Voigt form:
 * x_11, x_22, x_33, x_23, x_13, x_12, x_32, x_31, x_21
 */
inline FloatArrayF<9> to_voigt_form(const FloatMatrixF<3,3> &t)
{
    return {
        t(0, 0),
        t(1, 1),
        t(2, 2),
        t(1, 2),
        t(0, 2),
        t(0, 1),
        t(2, 1),
        t(2, 0),
        t(1, 0)
    };
}

/**
 * Symmetrizes and stores a matrix in stress Voigt form:
 * s_11, s_22, s_33, s_23, s_13, s_12
 */
inline FloatMatrixF<3,3> from_voigt_stress(const FloatArrayF<6> &v)
{
    return {
        v[0], v[5], v[4], 
        v[5], v[1], v[3], 
        v[4], v[3], v[2]
    }; // note: column-major order
}

/**
 * Symmetrizes and stores a matrix in stress Voigt form:
 * s_11, s_22, s_33, s_23, s_13, s_12
 */
inline FloatArrayF<6> to_voigt_stress(const FloatMatrixF<3,3> &t)
{
    return {
        t(0, 0),
        t(1, 1),
        t(2, 2),
        0.5 * ( t(1, 2) + t(2, 1) ),
        0.5 * ( t(0, 2) + t(2, 0) ),
        0.5 * ( t(0, 1) + t(1, 0) )
    };
}

/**
 * Converts Voigt vector to matrix form:
 * e_11, e_22, e_33, v_23, v_13, v_12
 */
inline FloatMatrixF<3,3> from_voigt_strain(const FloatArrayF<6> &v)
{
    return {
        v[0], 0.5*v[5], 0.5*v[4],
        0.5*v[5], v[1], 0.5*v[3],
        0.5*v[4], 0.5*v[3], v[2]
    }; // note: column-major order
}

/**
 * Symmetrizes and stores a matrix in strain Voigt form:
 * e_11, e_22, e_33, v_23, v_13, v_12
 */
inline FloatArrayF<6> to_voigt_strain(const FloatMatrixF<3,3> &t)
{
    return {
        t(0, 0),
        t(1, 1),
        t(2, 2),
        ( t(1, 2) + t(2, 1) ),
        ( t(0, 2) + t(2, 0) ),
        ( t(0, 1) + t(1, 0) )
    };
}

/**
 * Constructs diagonal matrix from vector.
 */
template<int N>
FloatMatrixF<N,N> diag(const FloatArrayF<N> &v)
{
    FloatMatrixF<N,N> out;
    for ( int i = 0; i < N; ++i ) {
        out(i, i) = v[i];
    }
    return out;
}

/**
 * Constructs diagonal matrix from vector.
 */
template<int N, int M>
FloatArrayF<N*M> flatten(const FloatMatrixF<N,M> &m)
{
    FloatArrayF<N*M> out;
    for (int i = 0; i < N*M; ++i) {
        out[i] = m[i];
    }
    return out;
}

/**
 * Swaps the fourth and sixth index in the array. This is to reorder the indices
 * from OOFEM's order to Abaqus' and vice versa.
 */
#if 0
void swap_46(FloatMatrixF<6,6> &t)
{
    std::swap( t(3, 0), t(5, 0) );
    std::swap( t(3, 1), t(5, 1) );
    std::swap( t(3, 2), t(5, 2) );
    std::swap( t(0, 3), t(0, 5) );
    std::swap( t(1, 3), t(1, 5) );
    std::swap( t(2, 3), t(2, 5) );
    std::swap( t(3, 3), t(5, 5) );
    std::swap( t(4, 3), t(4, 5) );
    std::swap( t(5, 3), t(3, 5) );
    std::swap( t(3, 4), t(5, 4) );
}

void swap_46(FloatMatrixF<9,9> &t)
{
    // OOFEM: 11, 22, 33, 23, 13, 12, 32, 31, 21
    // UMAT:  11, 22, 33, 12, 13, 23, 32, 21, 31

    std::array<int, 9> abq2oo = {0, 1, 2, 5, 4, 3, 6, 8, 7};

    FloatMatrixF<9,9> tmp;
    for ( int i = 0; i < 9; i++ ) {
        for ( int j = 0; j < 9; j++ ) {
            tmp(i, j) = t(abq2oo[ i ], abq2oo[ j ]);
        }
    }

    t = tmp;
}
#endif

/// Constructs an identity matrix
template<int N>
FloatMatrixF<N,N> eye()
{
    FloatMatrixF<N,N> out;
    for ( int i = 0; i < N; ++i ) {
        out(i, i) = 1.;
    }
    return out;
}

/// I_dev matrix in Voigt (stress) form
const FloatMatrixF<6,6> I_dev6 = {
    2./3., -1./3., -1/3., 0., 0., 0.,
    -1./3., 2./3., -1/3., 0., 0., 0.,
    -1./3., -1./3., 2/3., 0., 0., 0.,
    0., 0., 0., 0.5, 0., 0.,
    0., 0., 0., 0., 0.5, 0.,
    0., 0., 0., 0., 0., 0.5,
};

/// I(x)I expressed in Voigt form
const FloatMatrixF<6,6> I6_I6 = {
    1., 1., 1., 0., 0., 0.,
    1., 1., 1., 0., 0., 0.,
    1., 1., 1., 0., 0., 0.,
    0., 0., 0., 0., 0., 0.,
    0., 0., 0., 0., 0., 0.,
    0., 0., 0., 0., 0., 0.,
};


/**
 * Constructs "N" matrix
 * N_1 0   0   N_2 ...
 * 0   N_1 0   0   ...
 * 0   0   N_1 0   ... 
 */
template<int N, int NSD>
FloatMatrixF<NSD,N> Nmatrix(const FloatArrayF<N> &n)
{
    FloatMatrixF<NSD,N*NSD> out;
    for ( int i = 0; i < N; ++i ) {
        for ( int j = 0; j < NSD; ++j ) {
            out( j, i * NSD + j ) = n[j];
        }
    }
    return out;
}

/// Constructs the B matrix for 3D momentum balance problems.
template<int N>
FloatMatrixF<N,6> Bmatrix_3d(const FloatMatrixF<3,N> &dN)
{
    FloatMatrixF<N,6> B;
    for ( int j = 0, k = 0; j < dN.rows(); j++, k += 3 ) {
        B(0, k + 0) = B(3, k + 1) = B(4, k + 2) = dN(0, j);
        B(1, k + 1) = B(3, k + 0) = B(5, k + 2) = dN(1, j);
        B(2, k + 2) = B(4, k + 0) = B(5, k + 1) = dN(2, j);
    }
    return B;
}

/// Constructs the B matrix for plane stress momentum balance problems.
template<int N>
FloatMatrixF<3,N> Bmatrix_2d(const FloatMatrixF<2,N> &dN)
{
    FloatMatrixF<3,N*2> B;
    for ( int j = 0, k = 0; j < N; j++, k += 2 ) {
        B(0, k + 0) = B(2, k + 1) = dN(j, 0);
        B(1, k + 1) = B(2, k + 0) = dN(j, 1);
    }
    return B;
}

///@todo Move these more costly operations to a seperate header?
///@{

/**
 * Computes the  trace of the matrix.
 */
template<int N>
double trace(const FloatMatrixF<N,N> &mat)
{
    double s = 0.;
    for ( int i = 0; i < N; ++i ) {
        s += mat(i, i);
    }
    return s;
}

/**
 * Computes the Frobenius norm of the receiver.
 * The Frobenius norm is defined as the square root of the sum of the absolute squares of its elements.
 * @return Frobenius norm.
 */
template<int N>
double frobeniusNorm(const FloatMatrixF<N,N> &mat)
{
    double n = 0.;
    for ( int i = 0; i < N*N; ++i ) {
        n += mat[i] * mat[i];
    }
    return std::sqrt( n );
    //return std::sqrt( std::inner_product(mat.values.begin(), mat.values.end(), mat.values.begin(), 0.) );
}

/**
 * Computes the operator norm of the receiver.
 * @param p Norm type, 1 norm, else 2 norm.
 * @return Norm.
 */
template<int N>
double norm(const FloatMatrixF<N,N> &mat, int p=1)
{
    if ( p == 1 ) {
        double max_col = 0.;
        for ( int j = 0; j < N; j++ ) {
            double col_sum  = 0.;
            for ( int i = 0; i < N; i++ ) {
                col_sum += fabs( mat(i, j) );
            }
            if ( col_sum > max_col ) {
                max_col = col_sum;
            }
        }
        return max_col;
    } else if ( p == 2 ) {
        auto e = eig(transpose(mat) * mat, 1);
        return sqrt(e.first(0));
    }
}

/**
 * Computes the reciprocal conditioning of the receiver. From 0 to 1, where 0 is singular and 1 best.
 * The receiver must be square.
 * Works identically as MATLAB/Octaves rcond().
 * @param p Norm type, 1 norm, else 2 norm.
 * @return Conditioning of receiver.
 */
template<int N>
double rcond(const FloatMatrixF<N,N> &mat, int p=1)
{
    ///@todo Do we need lapack here?
    double anorm = norm(mat, p);
    if ( det(mat) <= 1e-6 * anorm ) {
        return 0.;
    }
    return 1. / ( norm(inv(mat)) * anorm, p );
}

/// Computes the determinant
template<int N>
double det(const FloatMatrixF<N,N> &mat)
{
    OOFEM_ERROR("TODO");
}

/// Computes the determinant
template<>
inline double det(const FloatMatrixF<2,2> &mat)
{
    return mat(0, 0) * mat(1, 1) - mat(0, 1) * mat(1, 0);
}

/// Computes the determinant
template<>
inline double det(const FloatMatrixF<3,3> &mat)
{
    return mat(0, 0) * mat(1, 1) * mat(2, 2) + mat(0, 1) * mat(1, 2) * mat(2, 0) +
           mat(0, 2) * mat(1, 0) * mat(2, 1) - mat(0, 2) * mat(1, 1) * mat(2, 0) -
           mat(1, 2) * mat(2, 1) * mat(0, 0) - mat(2, 2) * mat(0, 1) * mat(1, 0);
}

/// Computes the inverse
template<int N>
FloatMatrixF<N,N> inv(const FloatMatrixF<N,N> &mat)
{
    // gaussian elimination - slow but safe
    auto tmp = mat;
    // initialize answer to be unity matrix;
    auto out = eye<N>();
    // lower triangle elimination by columns
    for ( int i = 1; i < N; i++ ) {
        double piv = tmp.at(i, i);
        if ( std::abs(piv) < 1.e-24 ) {
            OOFEM_ERROR("pivot (%d,%d) to close to small (< 1.e-24)", i, i);
        }
        for ( int j = i + 1; j <= N; j++ ) {
            double linkomb = tmp.at(j, i) / tmp.at(i, i);
            for ( int k = i; k <= N; k++ ) {
                tmp.at(j, k) -= tmp.at(i, k) * linkomb;
            }
            for ( int k = 1; k <= N; k++ ) {
                out.at(j, k) -= out.at(i, k) * linkomb;
            }
        }
    }
    // upper triangle elimination by columns
    for ( int i = N; i > 1; i-- ) {
        double piv = tmp.at(i, i);
        for ( int j = i - 1; j > 0; j-- ) {
            double linkomb = tmp.at(j, i) / piv;
            for ( int k = i; k > 0; k-- ) {
                tmp.at(j, k) -= tmp.at(i, k) * linkomb;
            }
            for ( int k = N; k > 0; k-- ) {
                out.at(j, k) -= out.at(i, k) * linkomb;
            }
        }
    }
    // diagonal scaling
    for ( int i = 1; i <= N; i++ ) {
        for ( int j = 1; j <= N; j++ ) {
            out.at(i, j) /= tmp.at(i, i);
        }
    }
    return out;
}

/// Computes the inverse
template<>
inline FloatMatrixF<2,2> inv(const FloatMatrixF<2,2> &mat)
{
    FloatMatrixF<2,2> out;
    double d = det(mat);
    out(0, 0) = mat(1, 1) / d;
    out(1, 0) = -mat(1, 0) / d;
    out(0, 1) = -mat(0, 1) / d;
    out(1, 1) = mat(0, 0) / d;
    return out;
}

/// Computes the inverse
template<>
inline FloatMatrixF<3,3> inv(const FloatMatrixF<3,3> &mat)
{
    FloatMatrixF<3,3> out;
    double d = det(mat);
    out(0, 0) = ( mat(1, 1) * mat(2, 2) - mat(1, 2) * mat(2, 1) ) / d;
    out(1, 0) = ( mat(1, 2) * mat(2, 0) - mat(1, 0) * mat(2, 2) ) / d;
    out(2, 0) = ( mat(1, 0) * mat(2, 1) - mat(1, 1) * mat(2, 0) ) / d;
    out(0, 1) = ( mat(0, 2) * mat(2, 1) - mat(0, 1) * mat(2, 2) ) / d;
    out(1, 1) = ( mat(0, 0) * mat(2, 2) - mat(0, 2) * mat(2, 0) ) / d;
    out(2, 1) = ( mat(0, 1) * mat(2, 0) - mat(0, 0) * mat(2, 1) ) / d;
    out(0, 2) = ( mat(0, 1) * mat(1, 2) - mat(0, 2) * mat(1, 1) ) / d;
    out(1, 2) = ( mat(0, 2) * mat(1, 0) - mat(0, 0) * mat(1, 2) ) / d;
    out(2, 2) = ( mat(0, 0) * mat(1, 1) - mat(0, 1) * mat(1, 0) ) / d;
    return out;
}

/**
 * Computes (real) eigenvalues and eigenvectors of receiver (must be symmetric)
 * @param mat Matrix.
 * @param nf Number of significant figures.
 * @return Pair of eigenvalues and vectors..
 */
template<int N>
std::pair<FloatArrayF<N>, FloatMatrixF<N,N>> eig(FloatMatrixF<N,N> &mat, int nf)
{
    OOFEM_ERROR("TODO");
}

template<>
inline std::pair<FloatArrayF<2>, FloatMatrixF<2,2>> eig(FloatMatrixF<2,2> &mat, int nf)
{
    double b = - (mat(0,0) + mat(1,1))/2.;
    double c = mat(0,0) * mat(1,1) - mat(0,1) * mat(1,0);
    FloatArrayF<2> vals = {
        -b + std::sqrt(b*b - c),
        -b - std::sqrt(b*b - c)
    };
    FloatMatrixF<2,2> vecs;
    OOFEM_ERROR("TODO"); // is this even worth specializing? I don't think we ever use it for 2x2. Typically just 3x3, 6x6, 9x9
    return {vals, vecs};
}

/**
 * Solves the  system of linear equations @f$ K\cdot a = b @f$ .
 * Uses Gaussian elimination with pivoting directly on receiver.
 * @param b RHS of linear system.
 * @param answer Solution of linear equations.
 * @param transpose Solves for the transpose of K.
 * @return False if K is singular, otherwise true.
 */
template<int N>
bool solve(FloatMatrixF<N,N> &k, const FloatArrayF<N> &b, FloatArrayF<N> &answer, bool transpose=false)
{
    OOFEM_ERROR("TODO");
}

/**
 * Solves the  system of linear equations @f$ K\cdot A = B @f$ .
 * Uses Gaussian elimination with pivoting directly on receiver.
 * @param B RHS of linear system.
 * @param answer Solution of linear equations, each column corresponding to columns in B.
 * @param transpose Solves for the transpose of K.
 * @return False if K is singular, otherwise true.
 */
template<int N, int M>
bool solve(FloatMatrixF<N,N> &k, FloatMatrixF<N,M> &B, FloatMatrixF<N,M> &answer, bool transpose=false)
{
    OOFEM_ERROR("TODO");
}
///@}

} // end namespace oofem
#endif // floatmatrixf_h
