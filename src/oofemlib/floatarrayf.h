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

#ifndef floatarrayf_h
#define floatarrayf_h

#include "oofemcfg.h"
#include "contextioresulttype.h"
#include "contextmode.h"
#include "error.h"

#include <array>
#include <initializer_list>
#include <cmath>
#include <algorithm>
#include <numeric>

namespace oofem {

/**
 * Class representing vector of real numbers with compile time fixed size.
 * @author Mikael Öhman
 */
template<int N>
class OOFEM_EXPORT FloatArrayF
{
protected:
    /// Stored values.
    std::array< double, N > values;

public:
    /// @name Iterator for for-each loops:
    //@{
    auto begin() { return this->values.begin(); }
    auto end() { return this->values.end(); }
    auto begin() const { return this->values.begin(); }
    auto end() const { return this->values.end(); }
    //@}

    /// Ctor
    template<typename... V> FloatArrayF(V... x) : values{x...} { }
    /// Copy ctor
    FloatArrayF(FloatArrayF<N> &x) : values{x.values} { }
    /// Ctor from dynamic array (size must match)
    FloatArrayF(FloatArray &x)
    {
#ifndef NDEBUG
        if ( x.giveSize() != N ) {
            OOFEM_ERROR("Can't convert dynamic float array of size %d to fixed size %d\n", x.giveSize(), N);
        }
#endif
        std::copy(x.begin(), x.end(), values.begin());
    }

    /// Assignment operator
    void operator = (const FloatArrayF<N> &src) { values = src.values; }

    /**
     * Coefficient access function. Returns value of coefficient at given
     * position of the receiver. Provides 1-based indexing access.
     * @param i Position of coefficient in array.
     */
    inline double &at(int i)
    {
#ifndef NDEBUG
        this->checkBounds( i );
#endif
        return values [ i - 1 ];
    }
    /**
     * Coefficient access function. Returns l-value of coefficient at given
     * position of the receiver. Provides 1-based indexing access.
     * @param i Position of coefficient in array.
     */
    inline double at(int i) const
    {
#ifndef NDEBUG
        this->checkBounds( i );
#endif
        return values [ i - 1 ];
    }

    /**
     * Coefficient access function. Returns value of coefficient at given
     * position of the receiver. Provides 0-based indexing access.
     * @param i Position of coefficient in array.
     */
    inline double &operator[] (int i)
    {
#ifndef NDEBUG
        if ( i >= this->size() ) {
            OOFEM_ERROR( "array error on index : %d >= %d", i, this->size() );
        }
#endif
        return values [ i ];
    }
    /**
     * Coefficient access function. Returns value of coefficient at given
     * position of the receiver. Provides 0-based indexing access.
     * @param i Position of coefficient in array.
     */
    inline const double &operator() (int i) const { return this->operator[](i); } 
    inline const double &operator[] (int i) const { 
#ifndef NDEBUG
        if ( i >= this->size() ) {
            OOFEM_ERROR( "array error on index : %d >= %d", i, this->size() );
        }
#endif
        return values [ i ];
    }
    /**
     * Checks size of receiver towards requested bounds.
     * Current implementation will call exit(1), if dimension
     * mismatch found.
     * @param i Required size of receiver.
     */
    void checkBounds(int i) const
    {
        if ( i <= 0 ) {
            OOFEM_ERROR("array error on index : %d <= 0", i);
        } else if ( i > this->size() ) {
            OOFEM_ERROR("array error on index : %d > %d", i, this->size());
        }
    }
    /// Returns the size of receiver.
    int size() const { return N; }
    /**
     * Print receiver on stdout with custom name.
     * @param name Display name of reciever.
     */
    void printYourself(const std::string &name="FloatArrayF") const
    {
        printf("%s (%d): \n", name.c_str(), this->size());
        for ( double x: *this ) {
            printf( "%10.3e  ", x );
        }
        printf("\n");
    }
    /**
     * Gives the pointer to the raw data, breaking encapsulation.
     * @return Pointer to values of array
     */
    const double *givePointer() const { return values.data(); }
    double *givePointer() { return values.data(); }

    contextIOResultType storeYourself(DataStream &stream) const
    {
        if ( !stream.write(this->givePointer(), size) ) {
            return CIO_IOERR;
        }
        return CIO_OK;
    }

    contextIOResultType restoreYourself(DataStream &stream)
    {
        if ( !stream.read(this->givePointer(), size) ) {
            return CIO_IOERR;
        }
        return CIO_OK;
    }

    int givePackSize(DataStream &buff) const
    {
        return buff.givePackSizeOfDouble(this->size());
    }

    //friend class FloatMatrixF;
};

/// Convert to dynamic size matrix
template<int N>
FloatArray to_dynamic( const FloatArrayF<N> & x )
{
    FloatArray out(N);
    std::copy(x.begin(), x.end(), out.begin());
    return x;
}

/// Print to stream
template<int N>
std::ostream & operator << ( std::ostream & out, const FloatArrayF<N> & x )
{
    out << x.size();
    for ( double xi : x ) {
        out << " " << xi;
    }
    return out;
}

/// Simple math operations
template<int N>
FloatArrayF<N> operator * ( double a, const FloatArrayF<N> & x )
{
    FloatArrayF<N> out;
    for ( int i = 0; i < N; ++i ) {
        out[i] = x[i] * a;
    }
    return out;
}

template<int N>
FloatArrayF<N> operator * ( const FloatArrayF<N> & x, double a )
{
    return a*x;
}

template<int N>
FloatArrayF<N> operator / ( const FloatArrayF<N> & x, double a )
{
    FloatArrayF<N> out;
    for ( int i = 0; i < N; ++i ) {
        out[i] = x[i] * a;
    }
    return out;
}

template<int N>
FloatArrayF<N> operator ^ ( const FloatArrayF<N> & x, double a)
{
    FloatArrayF<N> out;
    for ( int i = 0; i < N; ++i ) {
        out[i] = std::pow(x[i], a);
    }
    return out;
}

template<int N>
FloatArrayF<N> operator + ( const FloatArrayF<N> & x, const FloatArrayF<N> & y )
{
    FloatArrayF<N> out;
    for ( int i = 0; i < N; ++i ) {
        out[i] = x[i] + y[i];
    }
    return out;
}

template<int N>
FloatArrayF<N> operator - ( const FloatArrayF<N> & x, const FloatArrayF<N> & y )
{
    FloatArrayF<N> out;
    for ( int i = 0; i < N; ++i ) {
        out[i] = x[i] - y[i];
    }
    return out;
}

template<int N>
FloatArrayF<N> &operator += ( FloatArrayF<N> & x, const FloatArrayF<N> & y )
{
    for ( int i = 0; i < N; ++i ) {
        x[i] += y[i];
    }
    return x;
}

template<int N>
FloatArrayF<N> &operator -= ( FloatArrayF<N> & x, const FloatArrayF<N> & y )
{
    for ( int i = 0; i < N; ++i ) {
        x[i] -= y[i];
    }
    return x;
}

template<int N>
FloatArrayF<N> &operator *= ( FloatArrayF<N> & x, double a )
{
    for ( int i = 0; i < N; ++i ) {
        x[i] *= a;
    }
    return x;
}

template<int N>
FloatArrayF<N> &operator /= ( FloatArrayF<N> & x, double a )
{
    for ( int i = 0; i < N; ++i ) {
        x[i] /= a;
    }
    return x;
}

template<int N>
FloatArrayF<N> operator ^= ( FloatArrayF<N> & x, double a)
{
    for ( int i = 0; i < N; ++i ) {
        x[i] = std::pow(x[i], a);
    }
}


/// Returns true if all coefficients of the receiver are 0, else false.
template<int N>
bool isZero(const FloatArrayF<N> &x)
{
    for ( auto &x : x ) {
        if ( x != 0. ) {
            return false;
        }
    }
    return true;
}

/// Returns true if all coefficients of the receiver are finite, else false.
template<int N>
bool isfinite( const FloatArrayF<N> &x )
{
    for ( auto &val : x ) {
        if ( !std::isfinite(val) ) {
            return false;
        }
    }
    return true;
}

/// Computes the L2 norm of x
template<int N>
double norm( const FloatArrayF<N> & x )
{
    double ans = 0.;
    for ( auto &val : x ) {
        ans += val * val;
    }
    return std::sqrt(ans);
}

/// Computes the sum of x
template<int N>
double sum( const FloatArrayF<N> & x )
{
    return std::accumulate(x.begin(), x.end(), 0.);
}

/// Computes the product of x
template<int N>
double product( const FloatArrayF<N> & x )
{
    return std::accumulate(x.begin(), x.end(), 1.0, [](double a, double b) { return a*b; });
}

/// Computes @$ x \cross y @$
FloatArrayF<3> cross( const FloatArrayF<3> & x, const FloatArrayF<3> & y )
{
    return {
        x[1] * y[2] - x[2] * y[1],
        x[2] * y[0] - x[0] * y[2],
        x[0] * y[1] - x[1] * y[0],
    };
}

/// Computes @$ x_i y_i @$
template<int N>
double dot( const FloatArrayF<N> & x, const FloatArrayF<N> & y )
{
    double ans = 0.;
    for ( int i = 0; i < N; ++i ) {
        ans += x[i] * y[i];
    }
    return ans;
}

/**
 * Swaps the fourth and sixth index in the array. This is to reorder the indices
 * from OOFEM's order to Abaqus' and vice versa.
 */
void swap_46(FloatArrayF<6> &t)
{
    std::swap( t[3], t[5] );
}

void swap_46(FloatArrayF<9> &t)
{
    // OOFEM: 11, 22, 33, 23, 13, 12, 32, 31, 21
    // UMAT:  11, 22, 33, 12, 13, 23, 32, 21, 31
    std::swap( t[3], t[5] );
    std::swap( t[7], t[8] );
}

} // end namespace oofem
#endif // floatarrayf_h
