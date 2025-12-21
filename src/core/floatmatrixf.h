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

#pragma once

#ifndef floatmatrixf_h
#define floatmatrixf_h

#include "oofemenv.h"
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


#define _LOOP_FLOATMATRIX(M,r,c) for(std::size_t c=0; c<M.cols(); c++) for(std::size_t r=0; r<M.rows(); r++)

/**
 * Implementation of matrix containing floating point numbers.
 * @author Mikael Öhman
 */
template<std::size_t N, std::size_t M>
class OOFEM_EXPORT FloatMatrixF
#ifdef _USE_EIGEN
: public MatrixRCd<N,M>
#endif
{
#ifdef _USE_EIGEN
public:
    typedef MatrixRCd<N,M> MatrixRCd_NM;
    OOFEM_EIGEN_DERIVED(FloatMatrixF,MatrixRCd_NM);
    const double *data() const { return MatrixRCd_NM::data(); }
    double *data() { return MatrixRCd_NM::data(); }
    template<
        typename... Args,
        class = typename std::enable_if_t<sizeof...(Args) == M*N>,
        class = std::enable_if_t<(std::conjunction_v<std::is_same<double, Args>...>)>
    >
    FloatMatrixF(Args... args) {
        std::initializer_list<double> ini{ args ... };
        assert(ini.size()==rows()*cols());
        int i=0;
        for(const double& x: ini){
            (*this)(i%rows(),i/rows())=x;
            i++;
        }
    }
    // redefine because of signedness (move to Index in the future pehaps)
    std::size_t rows() const{ return MatrixRCd_NM::rows(); }
    std::size_t cols() const{ return MatrixRCd_NM::cols(); }
#else
protected:
    /// Values of matrix stored column wise.
    std::array< double, N*M > values;
public:
    /**
     * Constructor (values are specified column-wise)
     * @note The syntax {{x,y,z},{...}} can be achieved by nested initializer_list, but 
     */
    template< typename... V, class = typename std::enable_if_t<sizeof...(V) == N*M> >
    FloatMatrixF(V... x) noexcept : values{x...} { }
    /**
     * Empty ctor, initializes to zero.
     */
    FloatMatrixF() noexcept : values{} { }
    /// Copy constructor.
    FloatMatrixF(const FloatMatrixF<N,M> &mat) noexcept : values(mat.values) {}
    /// FloatMatrix conversion constructor.
    FloatMatrixF(const FloatMatrix &mat)
    {
        #ifndef NDEBUG
            if ( mat.rows() != N || mat.cols() != M ) {
                throw std::out_of_range("Can't convert dynamic float matrix of size " +
                    std::to_string(mat.rows()) + "x" + std::to_string(mat.cols()) +
                    " to fixed size " + std::to_string(N) + "x" + std::to_string(M));
            }
        #endif
        #if 0
             std::copy_n(mat.begin(), N*M, values.begin());
        #else
             for(Index c=0; c<mat.cols(); c++) for(Index r=0; r<mat.rows(); r++) (*this)(r,c)=mat(r,c);
        #endif
    }

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
    void checkBounds(std::size_t i, std::size_t j) const
    {
        if ( i <= 0 ) {
            throw std::out_of_range("matrix error on rows : " + std::to_string(i) + " <= 0");
        } else if ( j <= 0 ) {
            throw std::out_of_range("matrix error on rows : " + std::to_string(j) + " <= 0");
        } else if ( i > N ) {
            throw std::out_of_range("matrix error on rows : " + std::to_string(i) + " < " + std::to_string(N));
        } else if ( j > M ) {
            throw std::out_of_range("matrix error on rows : " + std::to_string(j) + " < " + std::to_string(M));
        }
    }

    /// Returns number of rows of receiver.
    std::size_t rows() const { return N; }
    /// Returns number of columns of receiver.
    std::size_t cols() const { return M; }
    const double *data() const { return values.data(); }
    double *data() { return values.data(); }



    /**
     * Direct value access (column major). Implements 0-based indexing.
     * @param i Position in data.
     */
    double &operator[](std::size_t i)
    {
        return values[ i ];
    }
    /**
     * Direct value access (column major). Implements 0-based indexing.
     * @param i Position in data.
     */
    double operator[](std::size_t i) const
    {
        return values[ i ];
    }

    /**
     * Coefficient access function. Returns l-value of coefficient at given
     * position of the receiver. Implements 0-based indexing.
     * @param i Row position of coefficient.
     * @param j Column position of coefficient.
     */
    double &operator()(std::size_t i, std::size_t j)
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
    double operator()(std::size_t i, std::size_t j) const
    {
        #ifndef NDEBUG
            this->checkBounds(i + 1, j + 1);
        #endif
        return values [ j * N + i ];
    }

    /**
     * Extract sub matrix (can also reorder the matrix). Implements 0-based indexing.
     * @param r Rows to extract.
     * @param c Columns to extract.
     */
    template<std::size_t R, std::size_t C>
    FloatMatrixF<R,C> operator()(std::size_t const (&r)[R], std::size_t const (&c)[C]) const
    {
        FloatMatrixF<R,C> x;
        for ( std::size_t i = 0; i < R; ++i ) {
            for ( std::size_t j = 0; j < C; ++j ) {
                x(i, j) = (*this)(r[i], c[j]);
            }
        }
        return x;
    }
#endif
    /**
    * Coefficient access function. Returns value of coefficient at given
    * position of the receiver. Implements 1-based indexing.
    * @param i Row position of coefficient.
    * @param j Column position of coefficient.
    */
    double at(std::size_t i, std::size_t j) const
    {
        #ifndef NDEBUG
        this->checkBounds(i, j);
        #endif
        return (*this)(i-1, j-1);
    }
    /**
    * Coefficient access function. Returns value of coefficient at given
    * position of the receiver. Implements 1-based indexing.
    * @param i Row position of coefficient.
    * @param j Column position of coefficient.
    */
    inline double &at(std::size_t i, std::size_t j)
    {
        #ifndef NDEBUG
        this->checkBounds(i, j);
        #endif
        return (*this)(i-1, j-1);
    }

    static FloatMatrixF<N,M> fromColumns(FloatArrayF<N> const (&x)[M]) noexcept
    {
        FloatMatrixF<N,M> ret;
        for ( std::size_t r = 0; r < N; ++r) {
            for ( std::size_t c = 0; c < M; ++c) {
                ret(r,c) = x[c][r];
            }
        }
        return ret;
    }


    /// Assemble x into self.
    template<std::size_t R, std::size_t C>
    inline void assemble(const FloatMatrixF<R,C> &x, int const (&r)[R], int const (&c)[C] )
    {
        for ( std::size_t i = 0; i < R; ++i ) {
            for ( std::size_t j = 0; j < C; ++j ) {
                (*this)(r[i], c[j]) += x(i, j);
            }
        }
    }

    
    /**
     * Sets the values of the matrix in specified column. If matrix size is zero, the size is adjusted.
     * @param src Array to set at column c.
     * @param c Column to set.
     */
    void setColumn(const FloatArrayF<N> &src, std::size_t c)
    {
        for ( std::size_t i = 0; i < N; i++ ) {
            (*this)(i, c) = src[i];
        }
    }

    /**
     * Sets the values of the matrix in specified row.
     * @param src Array to set at row c.
     * @param r Column to set.
     */
    void setRow(const FloatArrayF<M> &src, int r)
    {
        for ( std::size_t i = 0; i < M; i++ ) {
            (*this)(r, i) = src[i];
        }
    }
    
    /**
     * Sets the values of the matrix in specified column. If matrix size is zero, the size is adjusted.
     * @param src Array to set at column c.
     * @param c Column to set.
     */
    FloatArrayF<N> column(std::size_t j) const
    {
        FloatArrayF<N> c;
        for ( std::size_t i = 0; i < N; i++ ) {
            c[i] = (*this)(i, j);
        }
        return c;
    }

    /**
     * Extract column from matrix
     */
    template<std::size_t C, class = typename std::enable_if_t<C < M>>
    FloatArrayF<N> column() const
    {
        FloatArrayF<N> c;
        for ( std::size_t i = 0; i < N; i++ ) {
            c[i] = (*this)(i, C);
        }
        return c;
    }

    /**
     * Extract row from matrix
     */
    template<std::size_t R, class = typename std::enable_if_t<R < N>>
    FloatArrayF<M> row() const
    {
        FloatArrayF<M> r;
        for ( std::size_t j = 0; j < M; j++ ) {
            r[j] = (*this)(R, j);
        }
        return r;
    }

    /**
     * Adds to the receiver the product @f$ a^{\mathrm{T}}\cdot b \mathrm{d}V @f$. If the receiver has zero size, it is expanded.
     * Assumes that receiver and product @f$ a^{\mathrm{T}}\cdot b \mathrm{d}V @f$ are symmetric matrices. Computes only the
     * upper half of receiver.
     * @param a Matrix a in equation.
     * @param b Matrix b in equation.
     * @param dV Scaling factor.
     */
    template<std::size_t P>
    void plusProductSymmUpper(const FloatMatrixF<P,N> &a, const FloatMatrixF<P,M> &b, double dV)
    {
        ///@todo Larger matrices needs optimized BLAS.
        for ( std::size_t i = 0; i < N; ++i ) {
            for ( std::size_t j = i; j < M; ++j ) {
                double summ = 0.;
                for ( std::size_t k = 0; k < P; ++k ) {
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
        for ( std::size_t i = 0; i < N; ++i ) {
            for ( std::size_t j = i; j < N; ++j ) {
                (*this)(i, j) += a[i] * a[j] * dV;
            }
        }
    }
    /**
     * Adds to the receiver the product @f$a^{\mathrm{T}} \cdot b \mathrm{d}V@f$.
     * @param a Matrix a in equation.
     * @param b Matrix b in equation.
     * @param dV Scaling factor.
     */
    template<std::size_t P>
    void plusProductUnsym(const FloatMatrixF<P,N> &a, const FloatMatrixF<P,M> &b, double dV)
    {
        ///@todo Larger matrices needs optimized BLAS.
        for ( std::size_t i = 0; i < N; ++i ) {
            for ( std::size_t j = 0; j < M; ++j ) {
                double summ = 0.;
                for ( std::size_t k = 0; k < P; ++k ) {
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
        for ( std::size_t i = 0; i < N; ++i ) {
            for ( std::size_t j = 0; j < M; ++j ) {
                (*this)(i, j) += a[i] * b[j] * dV;
            }
        }
    }
    /**
     * Initializes the lower half of the receiver according to the upper half.
     */
    void symmetrized()
    {
        for ( std::size_t i = 2; i <= this->rows(); i++ ) {
            for ( std::size_t j = 1; j < i; j++ ) {
                this->at(i, j) = this->at(j, i);
            }
        }
    }

    /// Prints matrix to stdout. Useful for debugging.
    void printYourself(const std::string &name="FloatMatrixF") const
    {
        printf("%s (%zu x %zu): \n", name.c_str(), N, M);
        if ( this->rows() <= 250 && this->cols() <= 250 ) {
            for ( std::size_t i = 0; i < this->rows(); ++i ) {
                for ( std::size_t j = 0; j < this->cols() && j < 100; ++j ) {
                    printf( "%10.3e  ", (*this)(i, j) );
                }
                printf("\n");
            }
        } else {
            for ( std::size_t i = 0; i < this->rows() && i < 20; ++i ) {
                for ( std::size_t j = 0; j < this->cols() && j < 10; ++j ) {
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
    const double *givePointer() const { return data(); }
    double *givePointer() { return data(); }

    contextIOResultType storeYourself(DataStream &stream) const
    {
        if ( !stream.write(data(), N*M) ) {
            return CIO_IOERR;
        }
        return CIO_OK;
    }

    contextIOResultType restoreYourself(DataStream &stream)
    {
        if ( !stream.read(data(), N*M) ) {
            return CIO_IOERR;
        }
        return CIO_OK;
    }

    int givePackSize(DataStream &buff) const
    {
        return buff.givePackSizeOfDouble(N*M);
    }
};


/// Assemble components into zero matrix.
template<std::size_t N, std::size_t M, std::size_t R, std::size_t C>
inline FloatMatrixF<N,M> assemble(const FloatMatrixF<R,C> &x, int const (&r)[R], int const (&c)[C] )
{
    FloatMatrixF<N,M> out;
    for ( std::size_t i = 0; i < R; ++i ) {
        for ( std::size_t j = 0; j < C; ++j ) {
            out(r[i], c[j]) = x(i, j);
        }
    }
    return out;
}


template<std::size_t N, std::size_t M>
std::ostream & operator << ( std::ostream & out, const FloatMatrixF<N,M> & x )
{
    out << N << " " << M << " {";
    for ( std::size_t i = 0; i < N; ++i ) {
        for ( std::size_t j = 0; j < M; ++j ) {
            out << " " << x(i, j);
        }
        out << ";";
    }
    out << "}";
    return out;
}

template<std::size_t N, std::size_t M>
FloatMatrixF<N,M> operator * ( double a, const FloatMatrixF<N,M> & x )
{
    FloatMatrixF<N,M> out;
    _LOOP_FLOATMATRIX(x,r,c) out(r,c)=x(r,c)*a;
    return out;
}

template<std::size_t N, std::size_t M>
FloatMatrixF<N,M> operator * ( const FloatMatrixF<N,M> & x, double a )
{
    return a*x;
}

template<std::size_t N, std::size_t M>
FloatMatrixF<N,M> operator / ( const FloatMatrixF<N,M> & x, double a )
{
    FloatMatrixF<N,M> out;
    _LOOP_FLOATMATRIX(x,r,c) out(r,c)=x(r,c)/a;
    return out;
}

template<std::size_t N, std::size_t M>
FloatMatrixF<N,M> operator + ( const FloatMatrixF<N,M> & x, const FloatMatrixF<N,M> & y )
{
    FloatMatrixF<N,M> out;
    _LOOP_FLOATMATRIX(x,r,c) out(r,c)=x(r,c)+y(r,c);
    return out;
}

template<std::size_t N, std::size_t M>
FloatMatrixF<N,M> operator - ( const FloatMatrixF<N,M> & x, const FloatMatrixF<N,M> & y )
{
    FloatMatrixF<N,M> out;
    _LOOP_FLOATMATRIX(x,r,c) out(r,c)=x(r,c)-y(r,c);
    return out;
}

template<std::size_t N, std::size_t M>
FloatMatrixF<N,M> &operator += ( FloatMatrixF<N,M> & x, const FloatMatrixF<N,M> & y )
{
    _LOOP_FLOATMATRIX(x,r,c) x(r,c)+=y(r,c);
    return x;
}

template<std::size_t N, std::size_t M>
FloatMatrixF<N,M> &operator -= ( FloatMatrixF<N,M> & x, const FloatMatrixF<N,M> & y )
{
    _LOOP_FLOATMATRIX(x,r,c) x(r,c)-=y(r,c);
    return x;
}

template<std::size_t N, std::size_t M>
FloatMatrixF<N,M> &operator *= ( FloatMatrixF<N,M> & x, double a )
{
    _LOOP_FLOATMATRIX(x,r,c) x(r,c)*=a;
    return x;
}

/// Returns true if no element is not finite (NAN or infinite)
template<std::size_t N, std::size_t M>
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
template<std::size_t N, std::size_t NSD>
FloatMatrixF<N*3,NSD> Nmatrix(const FloatArrayF<N> &n)
{
    FloatMatrixF<N*3,NSD> out;
    for ( std::size_t i = 0; i < N; ++i ) {
        for ( std::size_t j = 0; j < NSD; ++j ) {
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
    auto bigent = cross(tangent, normal);

    FloatMatrixF<3,3> out;
    //out.setColumn(normal, 0);
    //out.setColumn(bigent, 1);
    //out.setColumn(tangent, 2);
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
template<std::size_t N, std::size_t M>
FloatMatrixF<M,N> transpose(const FloatMatrixF<N,M> &mat)
{
    FloatMatrixF<M,N> out;
    for ( std::size_t i = 0; i < N; ++i ) {
        for ( std::size_t j = 0; j < M; ++j ) {
            out(j, i) = mat(i, j);
        }
    }
    return out;
}

/// Computes @f$ a \cdot b @f$.
template<std::size_t N, std::size_t M, std::size_t P>
FloatMatrixF<N,P> dot(const FloatMatrixF<N,M> &a, const FloatMatrixF<M,P> &b)
{
    FloatMatrixF<N,P> out;
    ///@todo BLAS for larger matrix sizes (maybe)
    for ( std::size_t i = 0; i < N; i++ ) {
        for ( std::size_t j = 0; j < P; j++ ) {
            double coeff = 0.;
            for ( std::size_t k = 0; k < M; k++ ) {
                coeff += a(i, k) * b(k, j);
            }
            out(i, j) = coeff;
        }
    }
    return out;
}

/// Computes @f$ a \cdot b^{\mathrm{T}} @f$.
template<std::size_t N, std::size_t M, std::size_t P>
FloatMatrixF<N,P> dotT(const FloatMatrixF<N,M> &a, const FloatMatrixF<P,M> &b)
{
    FloatMatrixF<N,P> out;
    ///@todo BLAS for larger matrix sizes (maybe)
    for ( std::size_t i = 0; i < N; i++ ) {
        for ( std::size_t j = 0; j < P; j++ ) {
            double coeff = 0.;
            for ( std::size_t k = 0; k < M; k++ ) {
                coeff += a(i, k) * b(j, k);
            }
            out(i, j) = coeff;
        }
    }
    return out;
}

/// Computes @f$ a^{\mathrm{T}} \cdot b @f$.
template<std::size_t N, std::size_t M, std::size_t P>
FloatMatrixF<N,P> Tdot(const FloatMatrixF<M,N> &a, const FloatMatrixF<M,P> &b)
{
    FloatMatrixF<N,P> out;
    ///@todo BLAS for larger matrix sizes (maybe)
    for ( std::size_t i = 0; i < N; i++ ) {
        for ( std::size_t j = 0; j < P; j++ ) {
            double coeff = 0.;
            for ( std::size_t k = 0; k < M; k++ ) {
                coeff += a(k, i) * b(k, j);
            }
            out(i, j) = coeff;
        }
    }
    return out;
}

/// Computes @$ m_ij x_j = m \cdot x @$
template<std::size_t N, std::size_t M>
FloatArrayF<N> dot(const FloatMatrixF<N,M> &m, const FloatArrayF<M> &x)
{
    FloatArrayF<N> out;
    for ( std::size_t i = 0; i < N; i++ ) {
        double sum = 0.;
        for ( std::size_t j = 0; j < M; j++ ) {
            sum += m(i, j) * x[j];
        }
        out[i] = sum;
    }
    return out;
}

/// Computes @$ x_j m_ji = x \cdot m = m^{\mathrm{T}} \cdot x @$
template<std::size_t N, std::size_t M>
FloatArrayF<N> dot(const FloatArrayF<M> &x, const FloatMatrixF<M,N> &m)
{
    FloatArrayF<N> out;
    for ( std::size_t i = 0; i < N; i++ ) {
        double sum = 0.;
        for ( std::size_t j = 0; j < M; j++ ) {
            sum += x[j] * m(j, i);
        }
        out[i] = sum;
    }
    return out;
}

/// Computes @$ x_j m_ji = x \cdot m = m^{\mathrm{T}} \cdot x @$
template<std::size_t N, std::size_t M>
FloatArrayF<N> Tdot(const FloatMatrixF<M,N> &m, const FloatArrayF<M> &x)
{
    return dot(x, m);
}

/// Computes the dyadic product @f$ m_{ij} = a_i b_j @f$.
template<std::size_t N, std::size_t M>
FloatMatrixF<N,M> dyad(const FloatArrayF<N> &a, const FloatArrayF<M> &b)
{
    FloatMatrixF<N,M> out;
    for ( std::size_t i = 0; i < N; ++i ) {
        for ( std::size_t j = 0; j < M; ++j ) {
            out(i, j) = a[i] * b[j];
        }
    }
    return out;
}

/// Computes @f$ a = r^{\mathrm{T}} \cdot a \cdot r @f$
template<std::size_t N, std::size_t M>
FloatMatrixF<M,M> rotate(FloatMatrixF<N,N> &a, const FloatMatrixF<N,M> &r)
{
    return dot(Tdot(r, a), r);
}

/// Computes @f$ a = r \cdot a \cdot r^{\mathrm{T}} @f$
template<std::size_t N, std::size_t M>
FloatMatrixF<M,M> unrotate(FloatMatrixF<N,N> &a, const FloatMatrixF<M,N> &r)
{
    return dotT(dot(r, a), r);
}

/**
 * Returns a copy of the given matrix column,
 */
template<std::size_t N, std::size_t M>
FloatArrayF<N> column(const FloatMatrixF<N,M> &mat, std::size_t col)
{
    FloatArrayF<N> ans;
    for ( std::size_t row = 0; row < N; ++row ) {
        ans[row] = mat(row, col);
    }
    return ans;
}

/**
 * Symmetrizes and stores a matrix in Voigt form:
 * x_11, x_22, x_33, x_23, x_13, x_12, x_32, x_31, x_21
 */
inline FloatMatrixF<3,3> from_voigt_form_9(const FloatArrayF<9> &v)
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
inline FloatArrayF<9> to_voigt_form_33(const FloatMatrixF<3,3> &t)
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
inline FloatMatrixF<3,3> from_voigt_stress_6(const FloatArrayF<6> &v)
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
inline FloatArrayF<6> to_voigt_stress_33(const FloatMatrixF<3,3> &t)
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
inline FloatMatrixF<3,3> from_voigt_strain_6(const FloatArrayF<6> &v)
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
inline FloatArrayF<6> to_voigt_strain_33(const FloatMatrixF<3,3> &t)
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
template<std::size_t N>
FloatMatrixF<N,N> diag(const FloatArrayF<N> &v)
{
    FloatMatrixF<N,N> out;
    for ( std::size_t i = 0; i < N; ++i ) {
        out(i, i) = v[i];
    }
    return out;
}


/**
 * Extract diagonal from matrix.
 */
template<std::size_t N>
FloatArrayF<N> diag(const FloatMatrixF<N,N> &m)
{
    FloatArrayF<N> out;
    for ( std::size_t i = 0; i < N; ++i ) {
        out[i] = m(i, i);
    }
    return out;
}

/**
 * Constructs diagonal matrix from vector.
 */
template<std::size_t N, std::size_t M>
FloatArrayF<N*M> flatten(const FloatMatrixF<N,M> &m)
{
    FloatArrayF<N*M> out;
    for(std::size_t c=0; c<m.cols(); c++) for(std::size_t r=0; r<m.rows(); r++) out[c*m.rows()+r]=m(r,c);
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

    std::array<std::size_t, 9> abq2oo = {0, 1, 2, 5, 4, 3, 6, 8, 7};

    FloatMatrixF<9,9> tmp;
    for ( std::size_t i = 0; i < 9; i++ ) {
        for ( std::size_t j = 0; j < 9; j++ ) {
            tmp(i, j) = t(abq2oo[ i ], abq2oo[ j ]);
        }
    }

    t = tmp;
}
#endif

/// Constructs a zero matrix (this is the default behavior when constructing a matrix, this is just for nicer syntax)
template<std::size_t N,std::size_t M>
FloatMatrixF<N,M> zero()
{
    return FloatMatrixF<N,M>();
}

/// Constructs an identity matrix
template<std::size_t N>
FloatMatrixF<N,N> eye()
{
    FloatMatrixF<N,N> out;
    for ( std::size_t i = 0; i < N; ++i ) {
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
template<std::size_t N, std::size_t NSD>
FloatMatrixF<NSD,N> Nmatrix(const FloatArrayF<N> &n)
{
    FloatMatrixF<NSD,N*NSD> out;
    for ( std::size_t i = 0; i < N; ++i ) {
        for ( std::size_t j = 0; j < NSD; ++j ) {
            out( j, i * NSD + j ) = n[j];
        }
    }
    return out;
}

/// Constructs the B matrix for 3D momentum balance problems.
template<std::size_t N>
FloatMatrixF<6, N*3> Bmatrix_3d(const FloatMatrixF<3,N> &dN)
{
    FloatMatrixF<6, N*3> B;
    for ( std::size_t j = 0, k = 0; j < dN.rows(); j++, k += 3 ) {
        B(0, k + 0) = B(3, k + 1) = B(4, k + 2) = dN(0, j);
        B(1, k + 1) = B(3, k + 0) = B(5, k + 2) = dN(1, j);
        B(2, k + 2) = B(4, k + 0) = B(5, k + 1) = dN(2, j);
    }
    return B;
}

/// Constructs the B matrix for plane stress momentum balance problems.
template<std::size_t N>
FloatMatrixF<3,N*2> Bmatrix_2d(const FloatMatrixF<2,N> &dN)
{
    FloatMatrixF<3,N*2> B;
    for ( std::size_t j = 0, k = 0; j < N; j++, k += 2 ) {
        B(0, k + 0) = B(2, k + 1) = dN(0, j);
        B(1, k + 1) = B(2, k + 0) = dN(1, j);
    }
    return B;
}

///@todo Move these more costly operations to a seperate header?
///@{

/**
 * Computes the  trace of the matrix.
 */
template<std::size_t N>
double trace(const FloatMatrixF<N,N> &mat)
{
    double s = 0.;
    for ( std::size_t i = 0; i < N; ++i ) {
        s += mat(i, i);
    }
    return s;
}

/**
 * Computes the Frobenius norm of the receiver.
 * The Frobenius norm is defined as the square root of the sum of the absolute squares of its elements.
 * @return Frobenius norm.
 */
template<std::size_t N>
double frobeniusNorm(const FloatMatrixF<N,N> &mat)
{
    double n = 0.;
    _LOOP_FLOATMATRIX(mat,r,c) n+=mat(r,c)*mat(r,c);
    return std::sqrt( n );
    //return std::sqrt( std::inner_product(mat.values.begin(), mat.values.end(), mat.values.begin(), 0.) );
}

/**
 * Computes the operator norm of the receiver.
 * @param p Norm type, 1 norm, else 2 norm.
 * @return Norm.
 */
template<std::size_t N>
double norm(const FloatMatrixF<N,N> &mat, int p=1)
{
    if ( p == 1 ) {
        double max_col = 0.;
        for ( std::size_t j = 0; j < N; j++ ) {
            double col_sum  = 0.;
            for ( std::size_t i = 0; i < N; i++ ) {
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
template<std::size_t N>
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
inline double det(const FloatMatrixF<2,2> &mat)
{
    return mat(0, 0) * mat(1, 1) - mat(0, 1) * mat(1, 0);
}

/// Computes the determinant
inline double det(const FloatMatrixF<3,3> &mat)
{
    return mat(0, 0) * mat(1, 1) * mat(2, 2) + mat(0, 1) * mat(1, 2) * mat(2, 0) +
           mat(0, 2) * mat(1, 0) * mat(2, 1) - mat(0, 2) * mat(1, 1) * mat(2, 0) -
           mat(1, 2) * mat(2, 1) * mat(0, 0) - mat(2, 2) * mat(0, 1) * mat(1, 0);
}


/// Computes the inverse
inline FloatMatrixF<2,2> inv_22(const FloatMatrixF<2,2> &mat, double /*zeropiv*/)
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
inline FloatMatrixF<3,3> inv_33(const FloatMatrixF<3,3> &mat, double /*zeropiv*/)
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


/// Computes the inverse
template<std::size_t N>
FloatMatrixF<N,N> inv(const FloatMatrixF<N,N> &mat, double zeropiv=1e-24)
{
    if constexpr (N==2) return inv_22(mat,zeropiv);
    if constexpr (N==3) return inv_33(mat,zeropiv);
    // gaussian elimination - slow but safe
    auto tmp = mat;
    // initialize answer to be unity matrix;
    auto out = eye<N>();
    // lower triangle elimination by columns
    for ( std::size_t  i = 1; i < N; i++ ) {
        double piv = tmp.at(i, i);
        if ( std::abs(piv) <= zeropiv ) {
            OOFEM_ERROR("pivot (%d,%d) to close to small", (int)i, (int)i);
        }
        for ( std::size_t j = i + 1; j <= N; j++ ) {
            double linkomb = tmp.at(j, i) / tmp.at(i, i);
            for ( std::size_t  k = i; k <= N; k++ ) {
                tmp.at(j, k) -= tmp.at(i, k) * linkomb;
            }
            for ( std::size_t  k = 1; k <= N; k++ ) {
                out.at(j, k) -= out.at(i, k) * linkomb;
            }
        }
    }
    // upper triangle elimination by columns
    for ( std::size_t  i = N; i > 1; i-- ) {
        double piv = tmp.at(i, i);
        for ( std::size_t  j = i - 1; j > 0; j-- ) {
            double linkomb = tmp.at(j, i) / piv;
            for ( std::size_t  k = i; k > 0; k-- ) {
                tmp.at(j, k) -= tmp.at(i, k) * linkomb;
            }
            for ( std::size_t  k = N; k > 0; k-- ) {
                out.at(j, k) -= out.at(i, k) * linkomb;
            }
        }
    }
    // diagonal scaling
    for ( std::size_t  i = 1; i <= N; i++ ) {
        for ( std::size_t  j = 1; j <= N; j++ ) {
            out.at(i, j) /= tmp.at(i, i);
        }
    }
    return out;
}


/**
 * Computes (real) eigenvalues and eigenvectors of receiver (must be symmetric)
 * @param mat Matrix.
 * @param nf Number of significant figures.
 * @return Pair of eigenvalues and vectors..
 */
template<std::size_t N>
std::pair<FloatArrayF<N>, FloatMatrixF<N,N>> eig(const FloatMatrixF<N,N> &mat, int nf=9)
{
    FloatArrayF<N> eval = diag(mat);
    FloatMatrixF<N,N> v = eye<N>();

    double sum = 0.0;
    for ( std::size_t i = 0; i < N; ++i ) {
        for ( std::size_t j = 0; j < N; ++j ) {
            sum += fabs( mat(i, j) );
        }
    }

    if ( sum <= 0.0 ) {
        return {zeros<N>(), eye<N>()};
    }

    auto m = mat;
    /* ---- REDUCE MATRIX TO DIAGONAL ---------------- */
    double c_b2 = .10;
    double tol = pow(c_b2, nf);
    int ite = 0;
    double ssum;
    do {
        ssum = 0.0;
        for ( std::size_t j = 1; j < N; ++j ) {
            for ( std::size_t i = 0; i < j; ++i ) {
                if ( ( fabs( m(i, j) ) / sum ) > tol ) {
                    ssum += fabs( m(i, j) );
                    /* ---- CALCULATE ROTATION ANGLE ----------------- */
                    double aa = atan2( m(i, j) * 2.0, eval[i] - eval[j] ) /  2.0;
                    double si = sin(aa);
                    double co = cos(aa);
                    /*
                     *   // ---- MODIFY "I" AND "J" COLUMNS OF "A" AND "V"
                     *   for (k = 0; k < neq; ++k) {
                     *    tt = m(k, i);
                     *    m(k, i) = co * tt + si * m(k, j);
                     *    m(k, j) = -si * tt + co * m(k, j);
                     *    tt = v(k, i);
                     *    v(k, i) = co * tt + si * v(k, j);
                     *    // L500:
                     *    v(k, j) = -si * tt + co * v(k, j);
                     *   }
                     *   // ---- MODIFY DIAGONAL TERMS --------------------
                     *   m(i, i) = co * m(i, i) + si * m(j, i);
                     *   m(j, j) = -si * m(i, j) + co * m(j, j);
                     *   m(i, j) = 0.0;
                     *   // ---- MAKE "A" MATRIX SYMMETRICAL --------------
                     *   for (k = 1; k <= neq; ++k) {
                     *    m(i, k) = m(k, i);
                     *    m(j, k) = m(k, j);
                     *    // L600:
                     *   }
                     */
                    // ---- MODIFY "I" AND "J" COLUMNS OF "A" AND "V"
                    for ( std::size_t k = 0; k < i; ++k ) {
                        double tt2 = m(k, i);
                        m(k, i) = co * tt2 + si * m(k, j);
                        m(k, j) = -si * tt2 + co * m(k, j);
                        tt2 = v(k, i);
                        v(k, i) = co * tt2 + si * v(k, j);
                        v(k, j) = -si * tt2 + co * v(k, j);
                    }

                    // diagonal term (i,i)
                    double tt = eval[i];
                    eval[i] = co * tt + si * m(i, j);
                    double aij = -si * tt + co * m(i, j);
                    tt = v(i, i);
                    v(i, i) = co * tt + si * v(i, j);
                    v(i, j) = -si * tt + co * v(i, j);

                    for ( std::size_t k = i + 1; k < j; ++k ) {
                        double tt2 = m(i, k);
                        m(i, k) = co * tt2 + si * m(k, j);
                        m(k, j) = -si * tt2 + co * m(k, j);
                        tt2 = v(k, i);
                        v(k, i) = co * tt2 + si * v(k, j);
                        v(k, j) = -si * tt2 + co * v(k, j);
                    }

                    // diagonal term (j,j)
                    tt = m(i, j);
                    double aji = co * tt + si * eval[j];
                    eval[j] = -si * tt + co * eval[j];

                    tt = v(j, i);
                    v(j, i) = co * tt + si * v(j, j);
                    v(j, j) = -si * tt + co * v(j, j);
                    //
                    for ( std::size_t k = j + 1; k < N; ++k ) {
                        double tt2 = m(i, k);
                        m(i, k) = co * tt2 + si * m(j, k);
                        m(j, k) = -si * tt2 + co * m(j, k);
                        tt2 = v(k, i);
                        v(k, i) = co * tt2 + si * v(k, j);
                        v(k, j) = -si * tt2 + co * v(k, j);
                    }

                    // ---- MODIFY DIAGONAL TERMS --------------------
                    eval[i] = co * eval[i] + si * aji;
                    eval[j] = -si * aij + co * eval[j];
                    m(i, j) = 0.0;
                } else {
                    /* ---- A(I,J) MADE ZERO BY ROTATION ------------- */
                }
            }
        }
        /* ---- CHECK FOR CONVERGENCE -------------------- */
        if ( ++ite > 50 ) {
            throw std::runtime_error("eig computation failed");
        }
    } while ( fabs(ssum) / sum > tol );
    return {eval, v};
}


/**
 * Computes (real) eigenvalues and eigenvectors of receiver (must be symmetric) using inverse iterations.
 * @param mat Matrix.
 * @return Pair of eigenvalues and vectors..
 */
template<std::size_t N>
std::pair<FloatArrayF<N>, FloatMatrixF<N,N>> eig_inverse(const FloatMatrixF<N,N> &mat)
{
    int nitem = 100;
    double rtol = 1e-12;

    FloatArrayF<N> w;
    std::array<FloatArrayF<N>, N> x;
    for ( std::size_t i = 0; i < N; ++i ) {
        w[i] += 1.;
        x[i][i] = 1.;
    }
    auto invmat = inv(mat);

    for ( int it = 0; it < nitem; it++ ) {
        auto z = x;
        //  solve matrix equation K.X = M.X
        for ( std::size_t i = 0; i < N; ++i ) {
            x[i] = dot(invmat, z[i]);
        }

        //  evaluation of Rayleigh quotients
        auto old_w = w;
        for ( std::size_t i = 0; i < N; i++ ) {
            w[i] = dot(z[i], x[i]) / dot(x[i], x[i]);
        }

        orthogonalize(x);

        //  check convergence
        if ( norm2(old_w - w) <= norm2(w) * rtol ) {
            break;
        }
    }

    return {w, x};
}


template<std::size_t N>
void orthogonalize(std::array<FloatArrayF<N>, N> &x)
{
    for ( std::size_t j = 0; j < N; j++ ) {
        auto t = x[j];
        for ( std::size_t ii = 0; ii < j; ii++ ) {
            x[j] -= dot(x[ii], t) * x[ii];
        }
        x[j] *= 1.0 / norm(x[j]);
    }
}


#if 0
template<>
inline std::pair<FloatArrayF<2>, FloatMatrixF<2,2>> eig(const FloatMatrixF<2,2> &mat, int nf)
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
#endif

/**
 * Solves the  system of linear equations @f$ K\cdot a = b @f$ .
 * Uses Gaussian elimination with pivoting directly on receiver.
 * @param b RHS of linear system.
 * @return Solution of linear equations.
 */
template<std::size_t N>
std::pair<bool, FloatArrayF<N>> solve_check(FloatMatrixF<N,N> mtrx, const FloatArrayF<N> &b, double zeropiv = 1e-20)
{
    auto answer = b;
    // initialize answer to be unity matrix;
    // lower triangle elimination by columns
    for ( std::size_t i = 0; i < N - 1; i++ ) {
        // find the suitable row and pivot
        double piv = fabs( mtrx(i, i) );
        std::size_t pivRow = i;
        for ( std::size_t j = i + 1; j < N; j++ ) {
            if ( fabs( mtrx(j, i) ) > piv ) {
                pivRow = j;
                piv = fabs( mtrx(j, i) );
            }
        }

        if ( piv <= zeropiv ) {
            return {false, zeros<N>()};
        }

        // exchange rows
        if ( pivRow != i ) {
            for ( std::size_t j = i; j < N; j++ ) {
                double help = mtrx(i, j);
                mtrx(i, j) = mtrx(pivRow, j);
                mtrx(pivRow, j) = help;
            }
            double help = answer[i];
            answer[i] = answer[pivRow];
            answer[pivRow] = help;
        }

        for ( std::size_t j = i + 1; j < N; j++ ) {
            double linkomb = mtrx(j, i) / mtrx(i, i);
            for ( std::size_t k = i; k < N; k++ ) {
                mtrx(j, k) -= mtrx(i, k) * linkomb;
            }

            answer[j] -= answer[i] * linkomb;
        }
    }

    // back substitution
    for ( int i = N - 1; i >= 0; i-- ) {
        double help = 0.;
        for ( std::size_t j = i + 1; j < N; j++ ) {
            help += mtrx(i, j) * answer[j];
        }

        answer[i] = ( answer[i] - help ) / mtrx(i, i);
    }
    return {true, answer};
}


/**
 * Solves the  system of linear equations @f$ K\cdot a = b @f$ .
 * Uses Gaussian elimination with pivoting directly on receiver.
 * @param b RHS of linear system.
 * @return Solution of linear equations.
 */
template<std::size_t N>
FloatArrayF<N> solve(FloatMatrixF<N,N> mtrx, const FloatArrayF<N> &b, double zeropiv=1e-20)
{
    auto tmp = solve_check(mtrx, b, zeropiv);
    if ( tmp.first ) {
        return tmp.second;
    } else {
        throw std::runtime_error("Singular pivot encountered");
    }
}

/**
 * Solves the  system of linear equations @f$ K\cdot A = B @f$ .
 * Uses Gaussian elimination with pivoting directly on receiver.
 * @param B RHS of linear system.
 * @return Solution of linear equations, each column corresponding to columns in B.
 */
template<std::size_t N, std::size_t M>
FloatMatrixF<N,M> solve(FloatMatrixF<N,N> mtrx, const FloatMatrixF<N,M> &B, double zeropiv=1e-20)
{
    auto answer = B;
    // initialize answer to be unity matrix;
    // lower triangle elimination by columns
    for ( std::size_t i = 0; i < N - 1; i++ ) {
        // find the suitable row and pivot
        double piv = fabs( mtrx(i, i) );
        std::size_t pivRow = i;
        for ( std::size_t j = i + 1; j < N; j++ ) {
            if ( fabs( mtrx(j, i) ) > piv ) {
                pivRow = j;
                piv = fabs( mtrx(j, i) );
            }
        }

        if ( fabs(piv) < zeropiv ) {
            throw std::runtime_error("pivot too small, matrix problem could not be solved");
        }

        // exchange rows
        if ( pivRow != i ) {
            for ( std::size_t j = i; j < N; j++ ) {
                double help = mtrx(i, j);
                mtrx(i, j) = mtrx(pivRow, j);
                mtrx(pivRow, j) = help;
            }

            for ( std::size_t j = 0; j < M; j++ ) {
                double help = answer(i, j);
                answer(i, j) = answer(pivRow, j);
                answer(pivRow, j) = help;
            }
        }

        for ( std::size_t j = i + 1; j < N; j++ ) {
            double linkomb = mtrx(j, i) / mtrx(i, i);
            for ( std::size_t k = i; k < N; k++ ) {
                mtrx(j, k) -= mtrx(i, k) * linkomb;
            }

            for ( std::size_t k = 0; k < M; k++ ) {
                answer(j, k) -= answer(i, k) * linkomb;
            }
        }
    }

    // back substitution
    for ( std::size_t i = N - 1; i >= 0; i-- ) {
        for ( std::size_t k = 0; k < M; k++ ) {
            double help = 0.;
            for ( std::size_t j = i + 1; j < N; j++ ) {
                help += mtrx(i, j) * answer(j, k);
            }

            answer(i, k) = ( answer(i, k) - help ) / mtrx(i, i);
        }
    }
    return answer;
}
///@}

} // end namespace oofem
#endif // floatmatrixf_h
