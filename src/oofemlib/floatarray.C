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

#include "floatarray.h"
#include "intarray.h"
#include "floatmatrix.h"
#include "mathfem.h"
#include "error.h"
#include "datastream.h"
#include "mathfem.h"

#include <cstdarg>
#include <cstdlib>
#include <cstring>
#include <ostream>
#include <memory>
#include <numeric>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

#define FAST_RESIZE(newsize) \
    if ( (newsize) < this->size() ) { \
        this->values.resize((newsize)); \
    } else if ( (newsize) > this->size() ) { \
        this->values.assign((newsize), 0.); \
    }

#ifdef __LAPACK_MODULE
extern "C" {
    extern void dgemv_(const char *trans, const int *m, const int *n, const double *alpha, const double *a, const int *lda, const double *x,
                       const int *incx, const double *beta, double *y, const int *incy, int aColumns, int xSize, int ySize);
    // Y = Y + alpha * X
    extern void daxpy_(const int *n, const double *alpha, const double *x, const int *incx, double *y, const int *incy, int xsize, int ysize);
}
#endif

namespace oofem {

bool FloatArray :: isFinite() const
{
    for(double val : values) {
        if(!std::isfinite(val)) {
            return false;
        }
    }

    return true;
}


void
FloatArray :: beScaled(double s, const FloatArray &b)
{
    FAST_RESIZE(b.size());

    for ( std::size_t i = 0; i < this->size(); ++i ) {
        (*this) [ i ] = s * b [ i ];
    }
}


void FloatArray :: add(const FloatArray &b)
// Performs the operation a=a+b, where a stands for the receiver. If the
// receiver's size is 0, adjusts its size to that of b. Returns the
// receiver.
{
    if ( b.isEmpty() ) {
        return;
    }

    if ( !this->size() ) {
        * this = b;
        return;
    }

#  ifndef NDEBUG
    if ( this->size() != b.size() ) {
        OOFEM_ERROR("dimension mismatch in a[%d]->add(b[%d])", this->giveSize(), b.giveSize());
    }

#  endif

#ifdef __LAPACK_MODULE
    int inc = 1;
    double s = 1.;
    int size = this->giveSize();
    daxpy_(& size, & s, b.givePointer(), & inc, this->givePointer(), & inc, b.giveSize(), this->giveSize());
#else
    for ( std::size_t i = 0; i < this->size(); i++ ) {
        (*this) [ i ] += b [ i ];
    }
#endif
}


void FloatArray :: add(double offset)
{
    for ( double &x: *this ) {
        x += offset;
    }
}


void FloatArray :: add(double factor, const FloatArray &b)
// Performs the operation a=a+factor*b, where a stands for the receiver. If the
// receiver's size is 0, adjusts its size to that of b.
{
    if ( this->isEmpty() ) {
        this->beScaled(factor, b);
        return;
    }

#  ifndef NDEBUG
    if ( this->size() != b.size() ) {
        OOFEM_ERROR("dimension mismatch in a[%d]->add(b[%d])", this->giveSize(), b.giveSize());
    }

#  endif

#ifdef __LAPACK_MODULE
    int inc = 1;
    int size = this->giveSize();
    daxpy_(& size, & factor, b.givePointer(), & inc, this->givePointer(), & inc, b.giveSize(), this->giveSize());
#else
    for ( std::size_t i = 0; i < this->size(); ++i ) {
        (*this) [ i ] += factor * b [ i ];
    }
#endif
}


void FloatArray :: plusProduct(const FloatMatrix &b, const FloatArray &s, double dV)
// Performs the operation a += b^T . s * dV
{
    std::size_t nRows = b.giveNumberOfRows();
    std::size_t nColumns = b.giveNumberOfColumns();

    if ( this->isEmpty() ) {
        this->values.assign( nColumns, 0. );
    }

#  ifndef NDEBUG
    if ( this->giveSize() != b.giveNumberOfColumns() ) {
        OOFEM_ERROR( "dimension mismatch in a[%d] and b[%d, *]", this->giveSize(), b.giveNumberOfColumns() );
    }
#  endif

#ifdef __LAPACK_MODULE
    double beta = 1.;
    int inc = 1;
    dgemv_("t", & nRows, & nColumns, & dV, b.givePointer(), & nRows, s.givePointer(), & inc, & beta, this->givePointer(), & inc, nColumns, nColumns, nRows);
#else
    for ( std::size_t i = 1; i <= nColumns; i++ ) {
        double sum = 0.;
        for ( std::size_t j = 1; j <= nRows; j++ ) {
            sum += b.at(j, i) * s.at(j);
        }
        this->at(i) += sum * dV;
    }
#endif
}


void FloatArray :: subtract(const FloatArray &src)
// Performs the operation a=a-src, where a stands for the receiver. If the
// receiver's size is 0, adjusts its size to that of src.
{
    if ( src.isEmpty() ) {
        return;
    }

    if ( this->isEmpty() ) {
        FAST_RESIZE(src.size());
        for ( std::size_t i = 0; i < this->size(); ++i ) {
            (*this) [ i ] = -src [ i ];
        }

        return;
    }

#  ifndef NDEBUG
    if ( this->size() != src.size() ) {
        OOFEM_ERROR("dimension mismatch in a[%d]->add(b[%d])", this->giveSize(), src.giveSize());
    }

#  endif

    for ( std::size_t i = 0; i < this->size(); ++i ) {
        (*this) [ i ] -= src [ i ];
    }
}


void FloatArray :: beMaxOf(const FloatArray &a, const FloatArray &b)
{
    std::size_t n = a.size();

    if ( a.size() == 0 ) {
        *this = b;
        return;
    } else if ( b.size() == 0 ) {
        *this = a;
        return;
    }

#  ifndef NDEBUG
    if ( n != b.size() ) {
        OOFEM_ERROR("dimension mismatch in beMaxOf(a[%d],b[%d])", n, b.giveSize());
    }

#  endif

    FAST_RESIZE(n);

    for ( std::size_t i = 0; i < n; i++ ) {
        (*this) [ i ] = max( a [ i ], b [ i ] );
    }
}


void FloatArray :: beMinOf(const FloatArray &a, const FloatArray &b)
{
    std::size_t n = a.size();

    if ( a.size() == 0 ) {
        *this = b;
        return;
    } else if ( b.size() == 0 ) {
        *this = a;
        return;
    }

#  ifndef NDEBUG
    if ( n != b.size() ) {
        OOFEM_ERROR("dimension mismatch in beMinOf(a[%d],b[%d])", n, b.giveSize());
    }

#  endif

    FAST_RESIZE(n);
    for ( std::size_t i = 0; i < n; i++ ) {
        (*this) [ i ] = min( a [ i ], b [ i ] );
    }
}


void FloatArray :: beDifferenceOf(const FloatArray &a, const FloatArray &b)
{
#ifndef NDEBUG
    if ( a.size() != b.size() ) {
        OOFEM_ERROR("size mismatch (%d : %d)", a.giveSize(), b.giveSize());
    }

#endif
#if 0
    FAST_RESIZE(a.giveSize());
    for ( int i = 0; i < this->giveSize(); ++i ) {
        (*this) [ i ] = a [ i ] - b [ i ];
    }
#else
    this->values.reserve(a.giveSize());
    this->values.resize(0);
    for ( std::size_t i = 0; i < a.size(); ++i ) {
        this->values.push_back( a[i] - b[i] );
    }

#endif
}

void FloatArray :: beDifferenceOf(const FloatArray &a, const FloatArray &b, std::size_t n)
{
#ifndef NDEBUG
    if ( a.size() < n || b.size() < n ) {
        OOFEM_ERROR("wrong size ", a.giveSize(), b.giveSize());
    }

#endif
    FAST_RESIZE(n);
    for ( std::size_t i = 0; i < n; ++i ) {
        (*this) [ i ] = a [ i ] - b [ i ];
    }
}


void FloatArray :: beSubArrayOf(const FloatArray &src, const IntArray &indx)
//
// Returns a subVector of the receiver constructed from an index array.
// this->at(i) = src.at(indx->at(i))
//
{
#ifndef NDEBUG
    if ( indx.maximum() > src.giveSize() || indx.minimum() < 1 ) {
        OOFEM_ERROR("index points outside of source");
    }
#endif

    std::size_t n = indx.size();
    FAST_RESIZE(n);
    for ( std::size_t i = 1; i <= n; i++ ) {
        this->at(i) = src.at( indx.at(i) );
    }
}


void FloatArray :: addSubVector(const FloatArray &src, std::size_t si)
{
    std::size_t reqSize, n = src.size();

    si--;
    reqSize = si + n;
    if ( this->size() < reqSize ) {
        this->resizeWithValues(reqSize);
    }

    for (std::size_t i = 0; i < n; i++ ) {
        (*this) [si + i] += src [ i ];
    }
}


void FloatArray :: beVectorProductOf(const FloatArray &v1, const FloatArray &v2)
{
#  ifndef NDEBUG
    // check proper bounds
    if ( ( v1.giveSize() != 3 ) || ( v2.giveSize() != 3 ) ) {
        OOFEM_ERROR("size mismatch, size is not equal to 3");
    }
#  endif

    FAST_RESIZE(3);

    this->at(1) = v1.at(2) * v2.at(3) - v1.at(3) * v2.at(2);
    this->at(2) = v1.at(3) * v2.at(1) - v1.at(1) * v2.at(3);
    this->at(3) = v1.at(1) * v2.at(2) - v1.at(2) * v2.at(1);
}

int FloatArray :: giveIndexMinElem()
{
    std::size_t index = 1;
    if ( !this->giveSize() ) {
        return -1;
    }
    double val = (*this) [ 0 ];
    for (std::size_t i = 1; i < this->size(); i++ ) {
        if ( val > (*this) [ i ] ) {
            val = (*this) [ i ];
            index = i + 1;
        }
    }
    return index;
}

int FloatArray :: giveIndexMaxElem()
{
    std::size_t index = 1;
    if ( !this->giveSize() ) {
        return -1;
    }
    double val = (*this) [ 0 ];
    for (std::size_t i = 1; i < this->size(); i++ ) {
        if ( val < (*this) [ i ] ) {
            val = (*this) [ i ];
            index = i + 1;
        }
    }
    return index;
}

double FloatArray :: dotProduct(const FloatArray &x) const
{
#  ifndef NDEBUG
    if ( this->size() != x.size() ) {
        OOFEM_ERROR("dimension mismatch in a[%d]->dotProduct(b[%d])", this->giveSize(), x.giveSize());
    }

#  endif

    return std::inner_product(this->begin(), this->end(), x.begin(), 0.);
}


double FloatArray :: dotProduct(const FloatArray &x, std::size_t size) const
{
#  ifndef NDEBUG
    if ( size > this->size() || size > x.size() ) {
        OOFEM_ERROR("dimension mismatch in a[%d]->dotProduct(b[%d])", this->giveSize(), x.giveSize());
    }

#  endif

    return std::inner_product(this->begin(), this->begin()+size, x.begin(), 0.);
}


double FloatArray :: distance(const FloatArray &x) const
{
    return sqrt( this->distance_square(x) );
}

double FloatArray :: distance(const FloatArray &iP1, const FloatArray &iP2, double &oXi, double &oXiUnbounded) const
{
    return sqrt( distance_square(iP1, iP2, oXi, oXiUnbounded) );
}

double FloatArray :: distance_square(const FloatArray &iP1, const FloatArray &iP2, double &oXi, double &oXiUnbounded) const
{
    const double l2 = iP1.distance_square(iP2);

    if ( l2 > 0.0 ) {

        const double s = (dot(*this, iP2) - dot(*this, iP1) ) - ( dot(iP1, iP2) - dot(iP1, iP1) );

        if ( s < 0.0 ) {
            // X is closest to P1
            oXi = 0.0;
            oXiUnbounded = s/l2;
            return this->distance_square(iP1);
        } else {
            if ( s > l2 ) {
                // X is closest to P2
                oXi = 1.0;
                oXiUnbounded = s/l2;
                return this->distance_square(iP2);
            } else {
                oXi = s / l2;
                oXiUnbounded = s/l2;
                const FloatArray q = ( 1.0 - oXi ) * iP1 + oXi * iP2;
                return this->distance_square(q);
            }
        }
    } else {
        // If the points P1 and P2 coincide,
        // we can compute the distance to any
        // of these points.
        oXi = 0.5;
        oXiUnbounded = 0.5;
        return this->distance_square(iP1);
    }
}


double FloatArray :: distance_square(const FloatArray &from) const
// returns distance between receiver and from from
// computed using generalized pythagorean formulae
{
    double dist = 0.;
    std::size_t s = min(this->size(), from.size());
    for (std::size_t i = 1; i <= s; ++i ) {
        double dx = this->at(i) - from.at(i);
        dist += dx * dx;
    }

    return dist;
}


void FloatArray :: assemble(const FloatArray &fe, const IntArray &loc)
// Assembles the array fe (typically, the load vector of a finite
// element) to the receiver, using loc as location array.
{
    std::size_t n = fe.size();
#  ifndef NDEBUG
    if ( n != loc.size() ) {
        OOFEM_ERROR("dimensions of 'fe' (%d) and 'loc' (%d) mismatch", fe.giveSize(), loc.giveSize() );
    }

#  endif

    for (std::size_t i = 1; i <= n; i++ ) {
        int ii = loc.at(i);
        if ( ii ) { // if non 0 coefficient,
            this->at(ii) += fe.at(i);
        }
    }
}


void FloatArray :: assembleSquared(const FloatArray &fe, const IntArray &loc)
// Assembles the array fe (typically, the load vector of a finite
// element) to the receiver, using loc as location array.
{
    std::size_t n = fe.size();
#  ifndef NDEBUG
    if ( n != loc.size() ) {
        OOFEM_ERROR("dimensions of 'fe' (%d) and 'loc' (%d) mismatch", fe.giveSize(), loc.giveSize() );
    }

#  endif

    for (std::size_t i = 1; i <= n; i++ ) {
        int ii = loc.at(i);
        if ( ii ) { // if non 0 coefficient,
            this->at(ii) += fe.at(i) * fe.at(i);
        }
    }
}


void FloatArray :: checkSizeTowards(const IntArray &loc)
// Expands the receiver if loc points to coefficients beyond the size of
// the receiver.
{
    int high = 0;
    for ( int p : loc ) {
        high = max( high, p );
    }

    if ( high > this->giveSize() ) {   // receiver must be expanded
        this->values.resize(high);
    }
}


void FloatArray :: reserve(int s)
{
    this->values.reserve(s);
    this->values.clear();
}


void FloatArray :: resizeWithValues(std::size_t n, std::size_t allocChunk)
{
#ifndef NDEBUG
    if ( allocChunk < 0 ) {
        OOFEM_FATAL("allocChunk must be non-negative; %d", allocChunk);
    }

#endif

    if ( allocChunk > 0 && this->values.capacity() < n ) {
        this->values.reserve(n + allocChunk);
    }

    this->values.resize(n);
}

void FloatArray :: resize(int n)
{
    this->values.resize(n);
    ///@todo Change to this (faster) 
    //this->values.assign(n, 0.);
}


void FloatArray :: hardResize(int n)
{
    this->values.assign(n, 0.);
    this->values.shrink_to_fit();
}


bool FloatArray :: containsOnlyZeroes() const
{
    for ( double x : *this ) {
        if ( x != 0. ) {
            return false;
        }
    }

    return true;
}


void FloatArray :: zero()
{
    std::fill(this->begin(), this->end(), 0.);
}


void FloatArray :: append(const FloatArray &a)
{
    this->values.insert(this->end(), a.begin(), a.end());
}


void FloatArray :: append(double a)
{
    this->values.push_back(a);
}


void FloatArray :: beProductOf(const FloatMatrix &aMatrix, const FloatArray &anArray)
// Stores the product of aMatrix * anArray in to receiver
{
    std::size_t nColumns = aMatrix.giveNumberOfColumns();
    std::size_t nRows = aMatrix.giveNumberOfRows();

    FAST_RESIZE(nRows);

#  ifndef NDEBUG
    if ( aMatrix.giveNumberOfColumns() != anArray.giveSize() ) {
        OOFEM_ERROR("dimension mismatch (%d, %d) . (%d)", 
                    aMatrix.giveNumberOfRows(), aMatrix.giveNumberOfColumns(), anArray.giveSize());
    }
#  endif

#ifdef __LAPACK_MODULE
    double alpha = 1., beta = 0.;
    int inc = 1;
    dgemv_("n", & nRows, & nColumns, & alpha, aMatrix.givePointer(), & nRows, anArray.givePointer(), & inc, & beta, this->givePointer(), & inc, nColumns, nColumns, nRows);
#else
    for (std::size_t i = 1; i <= nRows; i++ ) {
        double sum = 0.;
        for (std::size_t j = 1; j <= nColumns; j++ ) {
            sum += aMatrix.at(i, j) * anArray.at(j);
        }

        this->at(i) = sum;
    }
#endif
}


void FloatArray :: beTProductOf(const FloatMatrix &aMatrix, const FloatArray &anArray)
// Stores the product of aMatrix^T * anArray in to receiver
{
    std::size_t nRows = aMatrix.giveNumberOfRows();
    std::size_t nColumns = aMatrix.giveNumberOfColumns();

#  ifndef NDEBUG
    if ( aMatrix.giveNumberOfRows() != anArray.giveSize() ) {
        OOFEM_ERROR( "dimension mismatch, matrix rows = %d, array size = %d", aMatrix.giveNumberOfRows(), anArray.giveSize() );
    }

#  endif
    FAST_RESIZE(nColumns);

#ifdef __LAPACK_MODULE
    double alpha = 1., beta = 0.;
    int inc = 1;
    dgemv_("t", & nRows, & nColumns, & alpha, aMatrix.givePointer(), & nRows, anArray.givePointer(), & inc, & beta, this->givePointer(), & inc, nColumns, nColumns, nRows);
#else
    for (std::size_t i = 1; i <= nColumns; i++ ) {
        double sum = 0.;
        for (std::size_t j = 1; j <= nRows; j++ ) {
            sum += aMatrix.at(j, i) * anArray.at(j);
        }

        this->at(i) = sum;
    }
#endif
}


void FloatArray :: negated()
{
    for ( double &x: *this ) {
        x = -x;
    }
}


void FloatArray :: printYourself() const
// Prints the receiver on screen.
{
    printf("FloatArray of size : %d \n", this->giveSize());
    for ( double x: *this ) {
        printf( "%10.3e  ", x );
    }

    printf("\n");
}


void FloatArray :: printYourself(const std::string &name) const
// Prints the receiver on screen.
{
    printf("%s (%d): \n", name.c_str(), this->giveSize());
    for ( double x: *this ) {
        printf( "%10.3e  ", x );
    }

    printf("\n");
}

void FloatArray :: printYourselfToFile(const std::string filename, const bool showDimensions) const
// Prints the receiver to file.
{
    std :: ofstream arrayfile (filename);
    if (arrayfile.is_open()) {
        if (showDimensions)
            arrayfile << "FloatArray of size : " << this->giveSize() << "\n";
        arrayfile << std::scientific << std::right << std::setprecision(3);
        for ( double x: *this ) {
            arrayfile << std::setw(10) << x << "\t";
        }
        arrayfile.close();
    } else {
        OOFEM_ERROR("Failed to write to file");
    }
}

void FloatArray :: pY() const
// Prints the receiver on screen with higher accuracy than printYourself.
{
    printf("[");
    for ( double x: *this ) {
        printf( "%20.14e; ", x );
    }

    printf("];\n");
}


void FloatArray :: rotatedWith(FloatMatrix &r, char mode)
// Returns the receiver 'a' rotated according the change-of-base matrix r.
// If mode = 't', the method performs the operation  a = r(transp) * a .
// If mode = 'n', the method performs the operation  a = r * a .
{
    FloatArray rta;

    if ( mode == 't' ) {
        rta.beTProductOf(r, * this);
    } else if ( mode == 'n' ) {
        rta.beProductOf(r, * this);
    } else {
        OOFEM_ERROR("unsupported mode");
    }

    * this = rta;
}


void FloatArray :: times(double factor)
// Multiplies every coefficient of the receiver by factor.
{
    //std::transform(this->begin(), this->end(), this->begin(), std::bind2nd(std::multiplies<double>(), factor));
    for ( double &x : *this ) {
        x *= factor;
    }
}


double FloatArray :: normalize()
{
    double norm = this->computeNorm();
    if ( norm < 1.e-80 ) {
        OOFEM_ERROR("cannot norm receiver, norm is too small");
    }

    this->times(1. / norm);
    return norm;
}


double FloatArray :: computeNorm() const
{
    return sqrt( this->computeSquaredNorm() );
}


double FloatArray :: computeSquaredNorm() const
{
    return std::inner_product(this->begin(), this->end(), this->begin(), 0.);
}


double FloatArray :: sum() const
{
    return std::accumulate(this->begin(), this->end(), 0.);
}


double FloatArray :: product() const
{
    return std::accumulate(this->begin(), this->end(), 1.0, [](double a, double b) { return a*b; });
}


void FloatArray :: copySubVector(const FloatArray &src, int si)
{
    si--;
    this->resizeWithValues(si + src.giveSize());
    std :: copy( src.begin(), src.end(), this->begin() + si);
}


contextIOResultType FloatArray :: storeYourself(DataStream &stream) const
// writes receiver's binary image into stream
// use id to distinguish some instances
// return value >0 success
//              =0 file i/o error
{
    // write size
    std::size_t size = this->size();
    if ( !stream.write(size) ) {
        return CIO_IOERR;
    }

    // write raw data
    if ( size ) {
        if ( !stream.write(this->givePointer(), size) ) {
            return CIO_IOERR;
        }
    }

    // return result back
    return CIO_OK;
}

contextIOResultType FloatArray :: restoreYourself(DataStream &stream)
// reads receiver from stream
// warning - overwrites existing data!
// returns 0 if file i/o error
//        -1 if id od class id is not correct
{
    // read size
    std::size_t size;
    if ( !stream.read(size) ) {
        return CIO_IOERR;
    }

    this->values.resize(size);

    // read raw data
    if ( size ) {
        if ( !stream.read(this->givePointer(), size) ) {
            return CIO_IOERR;
        }
    }

    // return result back
    return CIO_OK;
}


int FloatArray :: givePackSize(DataStream &buff) const
{
    return buff.givePackSizeOfSizet(1) + buff.givePackSizeOfDouble(this->giveSize());
}

// IML compat

FloatArray &FloatArray :: operator = ( const double & val )
{
    std::fill(this->begin(), this->begin(), val);
    return * this;
}

FloatArray &operator *= ( FloatArray & x, const double & a )
{
    x.times(a);
    return x;
}

FloatArray operator *( const double & a, const FloatArray & x )
{
    FloatArray result;
    result.beScaled(a, x);
    return result;
}

FloatArray operator *( const FloatArray & x, const double & a )
{
    FloatArray result;
    result.beScaled(a, x);
    return result;
}

FloatArray operator + ( const FloatArray & x, const FloatArray & y )
{
    FloatArray result(x);
    result.add(y);
    return result;
}

FloatArray operator - ( const FloatArray & x, const FloatArray & y )
{
    FloatArray result;
    result.beDifferenceOf(x, y);
    return result;
}

FloatArray &operator += ( FloatArray & x, const FloatArray & y )
{
    x.add(y);
    return x;
}

FloatArray &operator -= ( FloatArray & x, const FloatArray & y )
{
    x.subtract(y);
    return x;
}

double dot(const FloatArray &x, const FloatArray &y)
{
    return x.dotProduct(y);
}

double distance(const FloatArray &x, const FloatArray &y)
{
    return x.distance(y);
}

double distance_square(const FloatArray &x, const FloatArray &y)
{
    return x.distance_square(y);
}

double norm(const FloatArray &x)
{
    return x.computeNorm();
}

double norm_square(const FloatArray &x)
{
    return x.computeSquaredNorm();
}

bool isfinite(const FloatArray &x)
{
    return x.isFinite();
}

bool iszero(const FloatArray &x)
{
    return x.containsOnlyZeroes();
}

double sum(const FloatArray & x)
{
    return x.sum();
}

double product(const FloatArray & x)
{
    return x.product();
}


// End of IML compat

void FloatArray :: beVectorForm(const FloatMatrix &aMatrix)
{
    // Rewrites the matrix on vector form, order: 11, 22, 33, 23, 13, 12, 32, 31, 21
#  ifndef NDEBUG
    if (  aMatrix.giveNumberOfColumns() != 3 || aMatrix.giveNumberOfRows() != 3 ) {
        OOFEM_ERROR("matrix dimension is not 3x3");
    }

#  endif
    *this = {
        aMatrix.at(1, 1),
        aMatrix.at(2, 2),
        aMatrix.at(3, 3),
        aMatrix.at(2, 3),
        aMatrix.at(1, 3),
        aMatrix.at(1, 2),
        aMatrix.at(3, 2),
        aMatrix.at(3, 1),
        aMatrix.at(2, 1)
    };
}

void FloatArray :: beSymVectorFormOfStrain(const FloatMatrix &aMatrix)
{
    // Revrites a symmetric strain matrix on reduced vector form, order: 11, 22, 33, 23, 13, 12
    // shear components are multiplied with a factor 2
#  ifndef NDEBUG
    if (  aMatrix.giveNumberOfColumns() != 3 || aMatrix.giveNumberOfRows() != 3 ) {
        OOFEM_ERROR("matrix dimension is not 3x3");
    }
#  endif

    *this = {
        aMatrix.at(1, 1),
        aMatrix.at(2, 2),
        aMatrix.at(3, 3),
        ( aMatrix.at(2, 3) + aMatrix.at(3, 2) ),
        ( aMatrix.at(1, 3) + aMatrix.at(3, 1) ),
        ( aMatrix.at(1, 2) + aMatrix.at(2, 1) )
    };
}


void FloatArray :: beSymVectorForm(const FloatMatrix &aMatrix)
{
    // Revrites the  matrix on vector form (symmetrized matrix used), order: 11, 22, 33, 23, 13, 12
#  ifndef NDEBUG
    if (  aMatrix.giveNumberOfColumns() != 3 || aMatrix.giveNumberOfRows() != 3 ) {
        OOFEM_ERROR("matrix dimension is not 3x3");
    }

#  endif

    *this = {
        aMatrix.at(1, 1),
        aMatrix.at(2, 2),
        aMatrix.at(3, 3),
        0.5 * ( aMatrix.at(2, 3) + aMatrix.at(3, 2) ),
        0.5 * ( aMatrix.at(1, 3) + aMatrix.at(3, 1) ),
        0.5 * ( aMatrix.at(1, 2) + aMatrix.at(2, 1) )
    };
}

void FloatArray :: changeComponentOrder()
{
    // OOFEM:           11, 22, 33, 23, 13, 12, 32, 31, 21
    // UMAT:            11, 22, 33, 12, 13, 23, 32, 21, 31

    if ( this->giveSize() == 6 ) {
        std :: swap( this->at(4), this->at(6) );
    } else if ( this->giveSize() == 9 )    {
        // OOFEM:       11, 22, 33, 23, 13, 12, 32, 31, 21
        // UMAT:        11, 22, 33, 12, 13, 23, 32, 21, 31
        const int abq2oo [ 9 ] = {
            1,  2,  3,  6,  5,  4,  7,  9,  8
        };

        FloatArray tmp(9);
        for ( int i = 0; i < 9; i++ ) {
            tmp(i) = this->at(abq2oo [ i ]);
        }

        * this = tmp;
    }
}


std :: ostream &operator << ( std :: ostream & out, const FloatArray & x )
{
    out << x.giveSize();
    for ( double xi : x ) {
        out << " " << xi;
    }
    return out;
}

//void FloatArray :: beReducedVectorFormOfStrain(const FloatMatrix &aMatrix)
//{
//    // Revrites the  matrix on vector form (symmetrized matrix used), order: 11, 22, 33, 23, 13, 12
//#  ifndef NDEBUG
//    if (  aMatrix.giveNumberOfColumns() !=3 || aMatrix.giveNumberOfColumns() !=3) {
//        OOFEM_ERROR("matrix dimension is not 3x3");
//    }
//
//#  endif
//
//    this->values.resize(6);
//    this->at(1) = aMatrix.at(1,1); this->at(2) = aMatrix.at(2,2); this->at(3) = aMatrix.at(3,3);
//    this->at(4) = ( aMatrix.at(2,3) + aMatrix.at(3,2) ); // Shear strains multiplied with a factor of 2
//    this->at(5) = ( aMatrix.at(1,3) + aMatrix.at(3,1) );
//    this->at(6) = ( aMatrix.at(1,2) + aMatrix.at(2,1) );
//}

void FloatArray :: power(const double exponent)
{
    for ( double &x : *this ) {
        x = pow(x, exponent);
    }
}



void FloatArray :: beColumnOf(const FloatMatrix &mat, int col)
{
    ///@todo This duplicates the "copyColumn" from FloatMatrix.
    mat.copyColumn(*this, col);
}

void FloatArray :: beRowOf(const FloatMatrix &mat, std::size_t row)
{
    std::size_t nRows = mat.giveRowSize();
    std::size_t nColumns = mat.giveColSize();

#  ifndef NDEBUG
    if (row>nRows) {
        OOFEM_ERROR( "dimension mismatch, matrix rows = %d, wanted row = %d", nRows, row)    
    }
#  endif    
    
    FAST_RESIZE(nColumns);
    for ( std::size_t i = 1; i <= nColumns; i++ ) {
        this->at(i) = mat.at(row,i);
    }
}
    
    
} // end namespace oofem
