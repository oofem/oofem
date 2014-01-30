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
/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
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
#if __cplusplus > 199711L
#include <memory>
#endif

#ifdef __PARALLEL_MODE
 #include "combuff.h"
#endif

#define ALLOC(size) new double[size];

#define RESIZE(n) \
    { \
        size = n; \
        if ( n > allocatedSize ) { \
            allocatedSize = n; \
            if ( values ) { delete[] values; } \
            values = ALLOC(size); \
        } \
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
FloatArray :: FloatArray() : size(0), allocatedSize(0), values(NULL) {}


FloatArray :: FloatArray(int n) :
    size(n),
    allocatedSize(n)
{
    if ( size ) {
        values = ALLOC(size);
        memset( values, 0, size * sizeof( double ) );
#ifdef DEBUG
        if ( !values ) {
            OOFEM_FATAL2("FloatArray :: FloatArray - Failed in allocating %d doubles", n);
        }

#endif
    } else {
        values = NULL;
    }
}

FloatArray :: FloatArray(const FloatArray &src) :
    size(src.size),
    allocatedSize(src.size)
{
    // copy constructor
    if ( size ) {
        values = ALLOC(size);
#ifdef DEBUG
        if ( !values ) {
            OOFEM_FATAL2("FloatArray :: FloatArray - Failed in allocating %d doubles", size);
        }
#endif
        memcpy( this->values, src.values, size * sizeof( double ) );
    } else {
        values = NULL;
    }
}

#if __cplusplus > 199711L
FloatArray :: FloatArray(std :: initializer_list< double >list)
{
    this->size = this->allocatedSize = (int)list.size();
    if ( this->size ) {
        this->values = ALLOC(this->size);
        std :: uninitialized_copy(list.begin(), list.end(), this->values);
    } else {
        this->values = NULL;
    }
}


FloatArray &FloatArray :: operator=(std :: initializer_list< double >list)
{
    RESIZE( (int)list.size() );
    std :: uninitialized_copy(list.begin(), list.end(), this->values);
    return * this;
}

#endif

FloatArray :: ~FloatArray()
{
    // Note! It is actually OK to delete NULL !
    if ( values ) {
        delete[] values;
    }
}

FloatArray &
FloatArray :: operator=(const FloatArray &src)
{
    if ( this != & src ) { // beware of s=s;
        RESIZE(src.size);
        memcpy( this->values, src.values, size * sizeof( double ) );
    }

    return * this;
}

#ifdef DEBUG
double &
FloatArray :: operator()(int i)
{
    if ( i >= size ) {
        OOFEM_ERROR2("FloatArray :: operator() : array error on index : %d <= 0 \n", i);
    }
    return values [ i ];
}

const double &
FloatArray :: operator()(int i) const
{
    if ( i >= size ) {
        OOFEM_ERROR2("FloatArray :: operator() : array error on index : %d <= 0 \n", i);
    }
    return values [ i ];
}
#endif


void
FloatArray :: setValues(int n, ...)
{
    va_list vl;
    va_start(vl, n);
    RESIZE(n);
    for ( int i = 0; i < n; i++ ) {
        this->values [ i ] = va_arg(vl, double);
    }
    va_end(vl);
}


void
FloatArray :: beScaled(double s, const FloatArray &b)
{
    RESIZE(b.size);
    for ( int i = 0; i < this->size; ++i ) {
        this->values [ i ] = s * b.values [ i ];
    }
}


void FloatArray :: add(const FloatArray &b)
// Performs the operation a=a+b, where a stands for the receiver. If the
// receiver's size is 0, adjusts its size to that of b. Returns the
// receiver.
{
    if ( b.size == 0 ) {
        return;
    }

    if ( !size ) {
        * this = b;
        return;
    }

#  ifdef DEBUG
    if ( size != b.size ) {
        OOFEM_ERROR3("FloatArray :: add :  dimension mismatch in a[%d]->add(b[%d])\n", size, b.size);
    }

#  endif

#ifdef __LAPACK_MODULE
    int inc = 1;
    double s = 1.;
    daxpy_(& size, & s, b.values, & inc, this->values, & inc, b.size, this->size);
#else
    for ( int i = 0; i < this->size; i++ ) {
        this->values [ i ] += b.values [ i ];
    }
#endif
}


void FloatArray :: add(double offset)
{
    for ( int i = 0; i < this->size; i++ ) {
        this->values [ i ] += offset;
    }
}


void FloatArray :: add(double factor, const FloatArray &b)
// Performs the operation a=a+factor*b, where a stands for the receiver. If the
// receiver's size is 0, adjusts its size to that of b.
{
    if ( this->size == 0 ) {
        this->beScaled(factor, b);
        return;
    }

#  ifdef DEBUG
    if ( size != b.size ) {
        OOFEM_ERROR3("FloatArray :: add :  dimension mismatch in a[%d]->add(b[%d])\n", size, b.size);
    }

#  endif

#ifdef __LAPACK_MODULE
    int inc = 1;
    daxpy_(& size, & factor, b.values, & inc, this->values, & inc, b.size, this->size);
#else
    for ( int i = 0; i < this->size; ++i ) {
        this->values [ i ] += factor * b.values [ i ];
    }
#endif
}


void FloatArray :: plusProduct(const FloatMatrix &b, const FloatArray &s, double dV)
// Performs the operation a += b^T . s * dV
{
    if ( this->size == 0 ) {
        this->beTProductOf(b, s);
        this->times(dV);
        return;
    }

#  ifdef DEBUG
    if ( size != b.giveNumberOfColumns() ) {
        OOFEM_ERROR3( "FloatArray :: plusProduct :  dimension mismatch in a[%d] and b[%d, *]\n", size, b.giveNumberOfColumns() );
    }
#  endif

#ifdef __LAPACK_MODULE
    int nRows = b.giveNumberOfRows();
    int nColumns = b.giveNumberOfColumns();
    double beta = 1.;
    int inc = 1;
    dgemv_("t", & nRows, & nColumns, & dV, b.givePointer(), & nRows, s.values, & inc, & beta, this->values, & inc, nColumns, nColumns, nRows);
#else
    for ( int i = 1; i <= b.giveNumberOfColumns(); i++ ) {
        double sum = 0.;
        for ( int j = 1; j <= b.giveNumberOfRows(); j++ ) {
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
    if ( src.giveSize() == 0 ) {
        return;
    }

    if ( !size ) {
        RESIZE(src.size);
        for ( int i = 0; i < this->size; ++i ) {
            this->values [ i ] = -src.values [ i ];
        }

        return;
    }

#  ifdef DEBUG
    if ( this->size != src.size ) {
        OOFEM_ERROR3("FloatArray dimension mismatch in a[%d]->add(b[%d])\n", size, src.size);
    }

#  endif

    for ( int i = 0; i < this->size; ++i ) {
        this->values [ i ] -= src.values [ i ];
    }
}


void FloatArray :: beMaxOf(const FloatArray &a, const FloatArray &b)
{
    int n = a.giveSize();

#  ifdef DEBUG
    if ( n != b.size ) {
        OOFEM_ERROR3("FloatArray dimension mismatch in beMaxOf(a[%d],b[%d])\n", n, b.size);
    }

#  endif

    RESIZE(n);
    for ( int i = 0; i < n; i++ ) {
        this->values [ i ] = max( a(i), b(i) );
    }
}


void FloatArray :: beMinOf(const FloatArray &a, const FloatArray &b)
{
    int n = a.giveSize();

#  ifdef DEBUG
    if ( n != b.size ) {
        OOFEM_ERROR3("FloatArray dimension mismatch in beMinOf(a[%d],b[%d])\n", n, b.size);
    }

#  endif

    RESIZE(n);
    for ( int i = 0; i < n; i++ ) {
        this->values [ i ] = min( a(i), b(i) );
    }
}


void FloatArray :: beDifferenceOf(const FloatArray &a, const FloatArray &b)
{
#ifdef DEBUG
    if ( a.size != b.size ) {
        OOFEM_ERROR3("FloatArray :: beDifferenceOf - size mismatch (%d : %d)", a.size, b.size);
    }

#endif
    RESIZE(a.size);
    for ( int i = 0; i < this->size; ++i ) {
        this->values [ i ] = a.values [ i ] - b.values [ i ];
    }
}

void FloatArray :: beDifferenceOf(const FloatArray &a, const FloatArray &b, int n)
{
#ifdef DEBUG
    if ( a.size < n || b.size < n ) {
        OOFEM_ERROR3("FloatArray :: beDifferenceOf - wrong size ", a.size, b.size);
    }

#endif
    RESIZE(n);
    for ( int i = 0; i < n; ++i ) {
        this->values [ i ] = a.values [ i ] - b.values [ i ];
    }
}


void FloatArray :: beSubArrayOf(const FloatArray &src, const IntArray &indx)
//
// Returns a subVector of the receiver constructed from an index array.
// this->at(i) = src.at(indx->at(i))
//
{
#ifdef DEBUG
    if ( indx.maximum() > src.giveSize() || indx.minimum() < 1 ) {
        OOFEM_ERROR("FloatArray :: beSubArrayOf - index points outside of source");
    }
#endif

    int n = indx.giveSize();
    RESIZE(n);
    for ( int i = 1; i <= n; i++ ) {
        this->at(i) = src.at( indx.at(i) );
    }
}


void FloatArray :: addSubVector(const FloatArray &src, int si)
{
    int reqSize, n = src.giveSize();

    si--;
    reqSize = si + n;
    if ( this->giveSize() < reqSize ) {
        this->resizeWithValues(reqSize);
    }

    for ( int i = 1; i <= n; i++ ) {
        this->at(si + i) += src.at(i);
    }
}


void FloatArray :: beVectorProductOf(const FloatArray &v1, const FloatArray &v2)
{
#  if DEBUG
    // check proper bounds
    if ( ( v1.giveSize() != 3 ) || ( v2.giveSize() != 3 ) ) {
        OOFEM_ERROR(" FloatArray::VectorProduct : size mismatch, size is not equal to 3");
    }
#  endif

    RESIZE(3);

    this->at(1) = v1.at(2) * v2.at(3) - v1.at(3) * v2.at(2);
    this->at(2) = v1.at(3) * v2.at(1) - v1.at(1) * v2.at(3);
    this->at(3) = v1.at(1) * v2.at(2) - v1.at(2) * v2.at(1);
}

int FloatArray :: giveIndexMinElem()
{
    int index = 1;
    if ( !this->size ) {
        return -1;
    }
    double val = this->values [ 0 ];
    for ( int i = 1; i < this->size; i++ ) {
        if ( val > this->values [ i ] ) {
            val = this->values [ i ];
            index = i + 1;
        }
    }
    return index;
}

int FloatArray :: giveIndexMaxElem()
{
    int index = 1;
    if ( !this->size ) {
        return -1;
    }
    double val = this->values [ 0 ];
    for ( int i = 1; i < this->size; i++ ) {
        if ( val < this->values [ i ] ) {
            val = this->values [ i ];
            index = i + 1;
        }
    }
    return index;
}

double FloatArray :: dotProduct(const FloatArray &x) const
{
#  ifdef DEBUG
    if ( this->size != x.size ) {
        OOFEM_ERROR3("FloatArray :: dotProduct :  dimension mismatch in a[%d]->dotProduct(b[%d])\n", this->size, x.size);
    }

#  endif

    double dp = 0;
    for ( int i = 0; i < this->size; i++ ) {
        dp += this->values [ i ] * x.values [ i ];
    }

    return dp;
}


double FloatArray :: dotProduct(const FloatArray &x, int size) const
{
#  ifdef DEBUG
    if ( size > this->size || size > x.size ) {
        OOFEM_ERROR3("FloatArray :: dotProduct :  dimension mismatch in a[%d]->dotProduct(b[%d])\n", this->size, x.size);
    }

#  endif

    double dp = 0.;
    for ( int i = 0; i < size; i++ ) {
        dp += this->values [ i ] * x.values [ i ];
    }

    return dp;
}


double FloatArray :: distance(const FloatArray &x) const
{
    return sqrt( this->distance_square(x) );
}

double FloatArray :: distance(const FloatArray &iP1, const FloatArray &iP2, double &oXi) const
{
    double dist = 0.0;

    // Vector from start P1 to point X
    FloatArray u;
    u.beDifferenceOf(* this, iP1);

    // Line tangent vector
    FloatArray t;
    t.beDifferenceOf(iP2, iP1);
    double l = norm(t);

    if ( l > 0.0 ) {
        t.normalize();
        double s = dot(u, t);

        if ( s < 0.0 ) {
            // X is closest to P1
            dist = this->distance(iP1);
            oXi = 0.0;
            return dist;
        } else {
            if ( s > l ) {
                // X is closest to P2
                dist = this->distance(iP2);
                oXi = 1.0;
                return dist;
            } else {
                oXi = s / l;
                FloatArray q = ( 1.0 - oXi ) * iP1 + oXi * iP2;
                dist = this->distance(q);
                return dist;
            }
        }
    } else {
        // If the points P1 and P2 coincide,
        // we can compute the distance to any
        // of these points.
        dist = this->distance(iP1);
        oXi = 0.5;
        return dist;
    }
}


double FloatArray :: distance_square(const FloatArray &from) const
// returns distance between receiver and from from
// computed using generalized pythagorean formulae
{
    double *p1, *p2, dx, dist = 0.;

    p1 = this->values;
    p2 = from.values;
    int s = min(size, from.size);
    while ( s-- ) {
        dx = ( * p1 ) - ( * p2 );
        dist += dx * dx;
        p1++;
        p2++;
    }

    return dist;
}


void FloatArray :: assemble(const FloatArray &fe, const IntArray &loc)
// Assembles the array fe (typically, the load vector of a finite
// element) to the receiver, using loc as location array.
{
    int n = fe.giveSize();
#  ifdef DEBUG
    if ( n != loc.giveSize() ) {
        OOFEM_ERROR3( "FloatArray::assemble : dimensions of 'fe' (%d) and 'loc' (%d) mismatch", fe.giveSize(), loc.giveSize() );
    }

#  endif

    for ( int i = 1; i <= n; i++ ) {
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
    int n = fe.giveSize();
#  ifdef DEBUG
    if ( n != loc.giveSize() ) {
        OOFEM_ERROR3( "FloatArray::assemble : dimensions of 'fe' (%d) and 'loc' (%d) mismatch", fe.giveSize(), loc.giveSize() );
    }

#  endif

    for ( int i = 1; i <= n; i++ ) {
        int ii = loc.at(i);
        if ( ii ) { // if non 0 coefficient,
            this->at(ii) += fe.at(i) * fe.at(i);
        }
    }
}


#ifdef DEBUG
double &FloatArray :: at(int i)
// Returns the i-th coefficient of the receiver. Slow but safe.
{
    this->checkBounds(i);
    return values [ i - 1 ];
}

double FloatArray :: at(int i) const
// Returns the i-th coefficient of the receiver. Slow but safe.
{
    this->checkBounds(i);
    return values [ i - 1 ];
}
#endif


#ifdef DEBUG
void FloatArray :: checkBounds(int i) const
// Checks that the receiver's size is not smaller than 'i'.
{
    if ( i <= 0 ) {
        OOFEM_ERROR2("FloatArray :: checkBounds : array error on index : %d <= 0 \n", i);
    }

    if ( i > size ) {
        OOFEM_ERROR3("FloatArray :: checkBounds : array error on index : %d > %d \n", i, size);
    }
}
#endif


void FloatArray :: checkSizeTowards(const IntArray &loc)
// Expands the receiver if loc points to coefficients beyond the size of
// the receiver.
{
    int n, high;

    high = 0;
    n = loc.giveSize();
    for ( int i = 1; i <= n; i++ ) {
        high = max( high, ( loc.at(i) ) );
    }

    if ( high > size ) {   // receiver must be expanded
        RESIZE(high);
    }
}


void FloatArray :: resizeWithValues(int n, int allocChunk)
{
#ifdef DEBUG
    if ( allocChunk < 0 ) {
        OOFEM_FATAL2("FloatArray :: resizeWithValues - allocChunk must be non-negative; %d", allocChunk);
    }

#endif

    if ( n <= allocatedSize ) {
        size = n;
        return;
    }
    allocatedSize = n + allocChunk;

    double *newValues = ALLOC(allocatedSize);
#ifdef DEBUG
    if ( !newValues ) {
        OOFEM_FATAL2("FloatArray :: resizeWithValues - Failed in allocating %d doubles", n + allocChunk);
    }
#endif
    memcpy( newValues, values, size * sizeof( double ) );
    memset( & newValues [ size ], 0, ( allocatedSize - size ) * sizeof( double ) );

    if ( values ) {
        delete[] values;
    }
    values = newValues;
    size = n;
}

void FloatArray :: resize(int n)
{
#if 0
    // Check if code is calling resize instead of resizeWithValues with minumum false-positives:
    if ( n != 0 && size != 0 && this->computeSquaredNorm() > 0 ) {
        this->printYourself();
        printf("check\n");
    }
#endif
#if 1
    ///@todo To be removed:
    this->resizeWithValues(n);
    return;

#endif

    size = n;
    if ( n <= allocatedSize ) {
        memset( this->values, 0, this->size * sizeof( double ) );
        return;
    }
    allocatedSize = n;

    if ( values ) {
        delete[] values;
    }
    values = ALLOC(allocatedSize);
    memset(this->values, 0, allocatedSize);
#ifdef DEBUG
    if ( !values ) {
        OOFEM_FATAL2("FloatArray :: simple - Failed in allocating %d doubles", n);
    }
#endif
}


void FloatArray :: hardResize(int n)
// Reallocates the receiver with new size.
{
    allocatedSize = n;
    double *newValues = ALLOC(allocatedSize);
#ifdef DEBUG
    if ( !newValues ) {
        OOFEM_FATAL2("FloatArray :: hardResize - Failed in allocating %d doubles", n);
    }
#endif
    memcpy( newValues, values, size * sizeof( double ) );
    memset( & newValues [ size ], 0, ( allocatedSize - size ) * sizeof( double ) );

    if ( values ) {
        delete[] values;
    }
    values = newValues;
    size = n;
}


bool FloatArray :: containsOnlyZeroes() const
{
    for ( int i = 0; i < this->size; i++ ) {
        if ( this->values [ i ] != 0. ) {
            return false;
        }
    }

    return true;
}


void FloatArray :: zero()
{
    memset( this->values, 0, this->size * sizeof( double ) );
}


void FloatArray :: beProductOf(const FloatMatrix &aMatrix, const FloatArray &anArray)
// Stores the product of aMatrix * anArray in to receiver
{
    int nColumns = aMatrix.giveNumberOfColumns();
    int nRows = aMatrix.giveNumberOfRows();

    RESIZE(nRows);

#  ifdef DEBUG
    if ( aMatrix.giveNumberOfColumns() != anArray.giveSize() ) {
        OOFEM_ERROR("FloatArray :: beProductOf : dimension mismatch");
    }
#  endif

#ifdef __LAPACK_MODULE
    double alpha = 1., beta = 0.;
    int inc = 1;
    dgemv_("n", & nRows, & nColumns, & alpha, aMatrix.givePointer(), & nRows, anArray.values, & inc, & beta, this->values, & inc, nColumns, nColumns, nRows);
#else
    for ( int i = 1; i <= nRows; i++ ) {
        double sum = 0.;
        for ( int j = 1; j <= nColumns; j++ ) {
            sum += aMatrix.at(i, j) * anArray.at(j);
        }

        this->at(i) = sum;
    }
#endif
}


void FloatArray :: beTProductOf(const FloatMatrix &aMatrix, const FloatArray &anArray)
// Stores the product of aMatrix^T * anArray in to receiver
{
    int nRows = aMatrix.giveNumberOfRows();
    int nColumns = aMatrix.giveNumberOfColumns();

#  ifdef DEBUG
    if ( aMatrix.giveNumberOfRows() != anArray.giveSize() ) {
        OOFEM_ERROR3("FloatArray :: beTProductOf : dimension mismatch, matrix rows = %d, array size = %d", aMatrix.giveNumberOfRows(), anArray.giveSize());
    }

#  endif
    RESIZE(nColumns);

#ifdef __LAPACK_MODULE
    double alpha = 1., beta = 0.;
    int inc = 1;
    dgemv_("t", & nRows, & nColumns, & alpha, aMatrix.givePointer(), & nRows, anArray.values, & inc, & beta, this->values, & inc, nColumns, nColumns, nRows);
#else
    for ( int i = 1; i <= nColumns; i++ ) {
        double sum = 0.;
        for ( int j = 1; j <= nRows; j++ ) {
            sum += aMatrix.at(j, i) * anArray.at(j);
        }

        this->at(i) = sum;
    }
#endif
}


void FloatArray :: negated()
{
    for ( int i = 0; i < this->size; i++ ) {
        this->values [ i ] = -this->values [ i ];
    }
}


void FloatArray :: printYourself() const
// Prints the receiver on screen.
{
    printf("FloatArray of size : %d \n", size);
    for ( int i = 1; i <= size; ++i ) {
        printf( "%10.3e  ", this->at(i) );
    }

    printf("\n");
}


void FloatArray :: pY() const
// Prints the receiver on screen with higher accuracy than printYourself.
{
    printf("[");
    for ( int i = 1; i <= size; ++i ) {
        printf( "%20.14e; ", this->at(i) );
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
        OOFEM_ERROR("FloatArray :: rotatedWith: unsupported mode");
    }

    * this = rta;
}


void FloatArray :: times(double factor)
// Multiplies every coefficient of the receiver by factor.
{
    for ( int i = 0; i < this->size; i++ ) {
        this->values [ i ] *= factor;
    }
}


double FloatArray :: normalize()
{
    double norm = this->computeNorm();
    if ( norm < 1.e-80 ) {
        OOFEM_ERROR("FloatArray::normalize : cannot norm receiver, norm is too small");
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
    int i;
    double *p, norm2 = 0.;

    p = this->values;
    i = this->size;
    while ( i-- ) {
        norm2 += ( * p ) * ( * p );
        p++;
    }

    return norm2;
}


double FloatArray :: sum() const
{
    int i;
    double *p, sum = 0;

    p = this->values;
    i = this->size;
    while ( i-- ) {
        sum += * ( p++ );
    }

    return sum;
}

void FloatArray :: copySubVector(const FloatArray &src, int si)
{
    si--;

    int reqSize = si + src.size;
    this->resizeWithValues(reqSize);

    memcpy( values + si, src.values, src.size * sizeof( double ) );
}


contextIOResultType FloatArray :: storeYourself(DataStream *stream, ContextMode mode)
// writes receiver's binary image into stream
// use id to distinguish some instances
// return value >0 success
//              =0 file i/o error
{
    // write size
    if ( !stream->write(& size, 1) ) {
        return CIO_IOERR;
    }

    // write raw data
    if ( size ) {
        if ( !stream->write(values, size) ) {
            return CIO_IOERR;
        }
    }

    // return result back
    return CIO_OK;
}

contextIOResultType FloatArray :: restoreYourself(DataStream *stream, ContextMode mode)
// reads receiver from stream
// warning - overwrites existing data!
// returns 0 if file i/o error
//        -1 if id od class id is not correct
{
    // read size
    if ( !stream->read(& size, 1) ) {
        return CIO_IOERR;
    }

    if ( size > allocatedSize ) {
        if ( values ) {
            delete[] values;
        }

        values = ALLOC(size);
#ifdef DEBUG
        if ( !values ) {
            OOFEM_FATAL2("FloatArray :: restoreYourself - Failed in allocating %d doubles", size);
        }

#endif
        allocatedSize = size;
    }

    // read raw data
    if ( size ) {
        if ( !stream->read(values, size) ) {
            return CIO_IOERR;
        }
    }

    // return result back
    return CIO_OK;
}


#ifdef __PARALLEL_MODE
int FloatArray :: packToCommBuffer(CommunicationBuffer &buff) const
{
    int result = 1;
    // pack size
    result &= buff.packInt(size);
    // pack data
    result &= buff.packArray(this->values, size);

    return result;
}

int FloatArray :: unpackFromCommBuffer(CommunicationBuffer &buff)
{
    int newSize, result = 1;
    // unpack size
    result &= buff.unpackInt(newSize);
    // resize yourself
    RESIZE(newSize);
    result &= buff.unpackArray(this->values, newSize);

    return result;
}

int FloatArray :: givePackSize(CommunicationBuffer &buff) const
{
    return buff.givePackSize(MPI_INT, 1) + buff.givePackSize(MPI_DOUBLE, this->size);
}

#endif

// IML compat

FloatArray &FloatArray :: operator=(const double &val)
{
    for ( int i = 0; i < size; i++ ) {
        values [ i ] = val;
    }

    return * this;
}

FloatArray &operator*=(FloatArray &x, const double &a)
{
    int N = x.giveSize();
    for ( int i = 0; i < N; i++ ) {
        x(i) *= a;
    }

    return x;
}


FloatArray operator*(const double &a, const FloatArray &x)
{
    int N = x.giveSize();
    FloatArray result(N);
    for ( int i = 0; i < N; i++ ) {
        result(i) = x(i) * a;
    }

    return result;
}

FloatArray operator*(const FloatArray &x, const double &a)
{
    // This is the other commutative case of vector*scalar.
    // It should be just defined to be
    // "return operator*(a,x);"
    // but some compilers (e.g. Turbo C++ v.3.0) have trouble
    // determining the proper template match.  For the moment,
    // we'll just duplicate the code in the scalar * vector
    // case above.

    int N = x.giveSize();
    FloatArray result(N);
    for ( int i = 0; i < N; i++ ) {
        result(i) = x(i) * a;
    }

    return result;
}

FloatArray operator+(const FloatArray &x, const FloatArray &y)
{
    int N = x.giveSize();
#ifdef DEBUG
    if ( N != y.giveSize() ) {
        OOFEM_ERROR("loatArray operator+ : incompatible vector lengths");
    }
#endif

    FloatArray result(N);
    for ( int i = 0; i < N; i++ ) {
        result(i) = x(i) + y(i);
    }

    return result;
}

FloatArray operator-(const FloatArray &x, const FloatArray &y)
{
    int N = x.giveSize();
#ifdef DEBUG
    if ( N != y.giveSize() ) {
        OOFEM_ERROR("FloatArray operator- : incompatible vector lengths");
    }
#endif

    FloatArray result(N);
    for ( int i = 0; i < N; i++ ) {
        result(i) = x(i) - y(i);
    }

    return result;
}


FloatArray &operator+=(FloatArray &x, const FloatArray &y)
{
    int N = x.giveSize();
#ifdef DEBUG
    if ( N != y.giveSize() ) {
        OOFEM_ERROR("FloatArray& operator+= : incompatible vector lengths");
    }
#endif

    for ( int i = 0; i < N; i++ ) {
        x(i) += y(i);
    }

    return x;
}


FloatArray &operator-=(FloatArray &x, const FloatArray &y)
{
    int N = x.giveSize();
#ifdef DEBUG
    if ( N != y.giveSize() ) {
        OOFEM_ERROR("FloatArray& operator-= : incompatible vector lengths");
    }
#endif

    for ( int i = 0; i < N; i++ ) {
        x(i) -= y(i);
    }

    return x;
}


double dot(const FloatArray &x, const FloatArray &y)
{
    //  Check for compatible dimensions:
#ifdef DEBUG
    if ( x.giveSize() != y.giveSize() ) {
        OOFEM_ERROR("dot : incompatible dimensions");
    }
#endif

    double temp =  0;
    for ( int i = 0; i < x.giveSize(); i++ ) {
        temp += x(i) * y(i);
    }

    return temp;
}

double norm(const FloatArray &x)
{
    double temp = oofem :: dot(x, x);
    return sqrt(temp);
}

// End of IML compat

void FloatArray :: beVectorForm(const FloatMatrix &aMatrix)
{
    // Rewrites the matrix on vector form, order: 11, 22, 33, 23, 13, 12, 32, 31, 21
#  ifdef DEBUG
    if (  aMatrix.giveNumberOfColumns() != 3 || aMatrix.giveNumberOfColumns() != 3 ) {
        OOFEM_ERROR("FloatArray :: beFullVectorForm : matrix dimension is not 3x3");
    }

#  endif
    RESIZE(9);
    this->at(1) = aMatrix.at(1, 1);
    this->at(2) = aMatrix.at(2, 2);
    this->at(3) = aMatrix.at(3, 3);
    this->at(4) = aMatrix.at(2, 3);
    this->at(5) = aMatrix.at(1, 3);
    this->at(6) = aMatrix.at(1, 2);
    this->at(7) = aMatrix.at(3, 2);
    this->at(8) = aMatrix.at(3, 1);
    this->at(9) = aMatrix.at(2, 1);
}

void FloatArray :: beSymVectorFormOfStrain(const FloatMatrix &aMatrix)
{
    // Revrites a symmetric strain matrix on reduced vector form, order: 11, 22, 33, 23, 13, 12
    // shear components are multiplied with a factor 2
#  ifdef DEBUG
    if (  aMatrix.giveNumberOfColumns() != 3 || aMatrix.giveNumberOfColumns() != 3 ) {
        OOFEM_ERROR("FloatArray :: beReducedVectorFormOfStrain : matrix dimension is not 3x3");
    }
#  endif

    RESIZE(6);
    this->at(1) = aMatrix.at(1, 1);
    this->at(2) = aMatrix.at(2, 2);
    this->at(3) = aMatrix.at(3, 3);
    this->at(4) = ( aMatrix.at(2, 3) + aMatrix.at(3, 2) );
    this->at(5) = ( aMatrix.at(1, 3) + aMatrix.at(3, 1) );
    this->at(6) = ( aMatrix.at(1, 2) + aMatrix.at(2, 1) );
}


void FloatArray :: beSymVectorForm(const FloatMatrix &aMatrix)
{
    // Revrites the  matrix on vector form (symmetrized matrix used), order: 11, 22, 33, 23, 13, 12
#  ifdef DEBUG
    if (  aMatrix.giveNumberOfColumns() != 3 || aMatrix.giveNumberOfColumns() != 3 ) {
        OOFEM_ERROR("FloatArray :: beReducedVectorForm : matrix dimension is not 3x3");
    }

#  endif

    RESIZE(6);
    this->at(1) = aMatrix.at(1, 1);
    this->at(2) = aMatrix.at(2, 2);
    this->at(3) = aMatrix.at(3, 3);
    this->at(4) = 0.5 * ( aMatrix.at(2, 3) + aMatrix.at(3, 2) );
    this->at(5) = 0.5 * ( aMatrix.at(1, 3) + aMatrix.at(3, 1) );
    this->at(6) = 0.5 * ( aMatrix.at(1, 2) + aMatrix.at(2, 1) );
}

void FloatArray :: changeComponentOrder()
{
    // OOFEM: 			11, 22, 33, 23, 13, 12, 32, 31, 21
	// UMAT:			11, 22, 33, 12, 13, 23, 32, 21, 31

	if(this->giveSize() == 6) {
		std::swap(this->at(4), this->at(6));
	}
	else if( this->giveSize() == 9 ) {
	    // OOFEM: 			11, 22, 33, 23, 13, 12, 32, 31, 21
		// UMAT:			11, 22, 33, 12, 13, 23, 32, 21, 31
		const int abq2oo[9] = {  1,  2,  3,  6,  5,  4,  7,  9,  8};

		FloatArray tmp(9);
		for(int i = 1; i <= 9; i++) {
				tmp.at(i) = this->at( abq2oo[i-1]);
		}

		*this = tmp;
    }
}


//void FloatArray :: beColumnOf(const FloatMatrix &mat, int col)
//{
//#  ifdef DEBUG
//    if (  col > mat.giveNumberOfColumns() ) {
//        OOFEM_ERROR3("FloatArray :: beColumnOf: column index (%d) exceeds number of columns in input matrix (%d)", col, mat.giveNumberOfColumns() );
//    }
//#  endif
//
//    int nRows = mat.giveNumberOfRows();
//    this->resize(nRows);
//    for ( int i = 1; i <= nRows; i++ ) {
//        this->at(i) = mat.at(i,col);
//    }
//
//}

std :: ostream &operator<<(std :: ostream &out, const FloatArray &x)
{
    out << x.size;
    for ( int i = 0; i < x.size; ++i ) {
        out << " " << x(i);
    }
    return out;
}

//void FloatArray :: beReducedVectorFormOfStrain(const FloatMatrix &aMatrix)
//{
//    // Revrites the  matrix on vector form (symmetrized matrix used), order: 11, 22, 33, 23, 13, 12
//#  ifdef DEBUG
//    if (  aMatrix.giveNumberOfColumns() !=3 || aMatrix.giveNumberOfColumns() !=3) {
//        OOFEM_ERROR("FloatArray :: beReducedVectorForm : matrix dimension is not 3x3");
//    }
//
//#  endif
//
//    this->resize(6);
//    this->at(1) = aMatrix.at(1,1); this->at(2) = aMatrix.at(2,2); this->at(3) = aMatrix.at(3,3);
//    this->at(4) = ( aMatrix.at(2,3) + aMatrix.at(3,2) ); // Shear strains multiplied with a factor of 2
//    this->at(5) = ( aMatrix.at(1,3) + aMatrix.at(3,1) );
//    this->at(6) = ( aMatrix.at(1,2) + aMatrix.at(2,1) );
//}


void FloatArray :: beColumnOf(const FloatMatrix &mat, int col)
{
#  ifdef DEBUG
    if ( col > mat.giveNumberOfColumns() ) {
        OOFEM_ERROR3( "FloatArray :: beColumnOf: column index (%d) exceeds number of columns in input matrix (%d)", col, mat.giveNumberOfColumns() );
    }
#  endif

    int nRows = mat.giveNumberOfRows();
    this->resize(nRows);
    for ( int i = 1; i <= nRows; i++ ) {
        this->at(i) = mat.at(i, col);
    }
}
} // end namespace oofem
