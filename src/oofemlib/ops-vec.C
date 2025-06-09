#include"floatmatrix.h"
#include"floatarray.h"
#include"intarray.h"
#ifndef __MATH_INTERNAL
    #define __MATH_INTERNAL
#endif
#include"ops-vec.h"

#include <numeric>
#include <cmath>
#include <fstream>
#include <ios>
#include <iomanip>



namespace oofem::vec{

void printYourself(const FloatArray& a)
// Prints the receiver on screen.
{
    printf("FloatArray of size : %d \n", a.giveSize());
    for ( double x: a ) {
        printf( "%10.3e  ", x );
    }

    printf("\n");
}


void printYourself(const FloatArray& a, const std::string &name)
// Prints the receiver on screen.
{
    printf("%s (%d): \n", name.c_str(), a.giveSize());
    for ( double x: a ) {
        printf( "%10.3e  ", x );
    }

    printf("\n");
}

void printYourselfToFile(const FloatArray& a, const std::string filename, const bool showDimensions)
// Prints the receiver to file.
{
    std :: ofstream arrayfile (filename);
    if (arrayfile.is_open()) {
        if (showDimensions)
            arrayfile << "FloatArray of size : " << a.giveSize() << "\n";
        arrayfile << std::scientific << std::right << std::setprecision(3);
        for ( double x: a ) {
            arrayfile << std::setw(10) << x << "\t";
        }
        arrayfile.close();
    } else {
        OOFEM_ERROR("Failed to write to file");
    }
}

void pY(const FloatArray& a)
// Prints the receiver on screen with higher accuracy than printYourself.
{
    printf("[");
    for ( double x: a ) {
        printf( "%20.14e; ", x );
    }

    printf("];\n");
}


void rotatedWith(FloatArray& a, FloatMatrix &r, char mode)
// Returns the receiver 'a' rotated according the change-of-base matrix r.
// If mode = 't', the method performs the operation  a = r(transp) * a .
// If mode = 'n', the method performs the operation  a = r * a .
{
    FloatArray rta;

    if ( mode == 't' ) {
        rta.beTProductOf(r, a);
    } else if ( mode == 'n' ) {
        rta.beProductOf(r, a);
    } else {
        OOFEM_ERROR("unsupported mode");
    }

    a = rta;
}


double distance(const FloatArray& pos, const FloatArray &iP1, const FloatArray &iP2, double &oXi, double &oXiUnbounded)
{
    return sqrt( distance_square(pos, iP1, iP2, oXi, oXiUnbounded) );
}

double distance_square(const FloatArray& pos, const FloatArray &iP1, const FloatArray &iP2, double &oXi, double &oXiUnbounded)
{
    const double l2 = vec::distance_square(iP1,iP2);

    if ( l2 > 0.0 ) {

        const double s = (vec::dot(pos, iP2) - vec::dot(pos, iP1) ) - ( vec::dot(iP1, iP2) - vec::dot(iP1, iP1) );

        if ( s < 0.0 ) {
            // X is closest to P1
            oXi = 0.0;
            oXiUnbounded = s/l2;
            return vec::distance_square(pos,iP1);
        } else {
            if ( s > l2 ) {
                // X is closest to P2
                oXi = 1.0;
                oXiUnbounded = s/l2;
                return vec::distance_square(pos,iP2);
            } else {
                oXi = s / l2;
                oXiUnbounded = s/l2;
                const FloatArray q = ( 1.0 - oXi ) * iP1 + oXi * iP2;
                return vec::distance_square(pos,q);
            }
        }
    } else {
        // If the points P1 and P2 coincide,
        // we can compute the distance to any
        // of these points.
        oXi = 0.5;
        oXiUnbounded = 0.5;
        return vec::distance_square(pos,iP1);
    }
}



FloatArray beVectorForm( const FloatMatrix &aMatrix)
{
    // Rewrites the matrix on vector form, order: 11, 22, 33, 23, 13, 12, 32, 31, 21
#  ifndef NDEBUG
    if (  aMatrix.giveNumberOfColumns() != 3 || aMatrix.giveNumberOfRows() != 3 ) {
        OOFEM_ERROR("matrix dimension is not 3x3");
    }

#  endif
    FloatArray ret;
    ret = Vec9(
        aMatrix.at(1, 1),
        aMatrix.at(2, 2),
        aMatrix.at(3, 3),
        aMatrix.at(2, 3),
        aMatrix.at(1, 3),
        aMatrix.at(1, 2),
        aMatrix.at(3, 2),
        aMatrix.at(3, 1),
        aMatrix.at(2, 1)
    );
    return ret;
}

FloatArray beSymVectorFormOfStrain(const FloatMatrix &aMatrix)
{
    // Revrites a symmetric strain matrix on reduced vector form, order: 11, 22, 33, 23, 13, 12
    // shear components are multiplied with a factor 2
#  ifndef NDEBUG
    if (  aMatrix.giveNumberOfColumns() != 3 || aMatrix.giveNumberOfRows() != 3 ) {
        OOFEM_ERROR("matrix dimension is not 3x3");
    }
#  endif
    FloatArray ret;
    ret = Vec6(
        aMatrix.at(1, 1),
        aMatrix.at(2, 2),
        aMatrix.at(3, 3),
        ( aMatrix.at(2, 3) + aMatrix.at(3, 2) ),
        ( aMatrix.at(1, 3) + aMatrix.at(3, 1) ),
        ( aMatrix.at(1, 2) + aMatrix.at(2, 1) )
    );
    return ret;
}


FloatArray beSymVectorForm(const FloatMatrix &aMatrix)
{
    // Revrites the  matrix on vector form (symmetrized matrix used), order: 11, 22, 33, 23, 13, 12
#  ifndef NDEBUG
    if (  aMatrix.giveNumberOfColumns() != 3 || aMatrix.giveNumberOfRows() != 3 ) {
        OOFEM_ERROR("matrix dimension is not 3x3");
    }

#  endif
    FloatArray ret;
    ret = Vec6(
        aMatrix.at(1, 1),
        aMatrix.at(2, 2),
        aMatrix.at(3, 3),
        0.5 * ( aMatrix.at(2, 3) + aMatrix.at(3, 2) ),
        0.5 * ( aMatrix.at(1, 3) + aMatrix.at(3, 1) ),
        0.5 * ( aMatrix.at(1, 2) + aMatrix.at(2, 1) )
    );
    return ret;
}

void changeComponentOrder(FloatArray& a)
{
    // OOFEM:           11, 22, 33, 23, 13, 12, 32, 31, 21
    // UMAT:            11, 22, 33, 12, 13, 23, 32, 21, 31

    if ( a.giveSize() == 6 ) {
        std :: swap( a.at(4), a.at(6) );
    } else if ( a.giveSize() == 9 )    {
        // OOFEM:       11, 22, 33, 23, 13, 12, 32, 31, 21
        // UMAT:        11, 22, 33, 12, 13, 23, 32, 21, 31
        const int abq2oo [ 9 ] = {
            1,  2,  3,  6,  5,  4,  7,  9,  8
        };

        FloatArray tmp(9);
        for ( int i = 0; i < 9; i++ ) {
            tmp(i) = a.at(abq2oo [ i ]);
        }

        a = tmp;
    }
}



void assemble(FloatArray& dst, const FloatArray &fe, const IntArray &loc)
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
            dst.at(ii) += fe.at(i);
        }
    }
}


void assembleSquared(FloatArray& dst, const FloatArray &fe, const IntArray &loc)
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
            dst.at(ii) += fe.at(i) * fe.at(i);
        }
    }
}





double dot(const FloatArray &x, const FloatArray &y)
{
    return x.dotProduct(y);
}

double distance(const FloatArray &x, const FloatArray &y)
{
    return std::sqrt(vec::distance_square(x,y)); // x.distance(y);
}

double distance_square(const FloatArray &x, const FloatArray &y)
{
    double dist = 0.;
    std::size_t s = std::min(x.size(), y.size());
    for (std::size_t i = 1; i <= s; ++i ) {
        double dx = x.at(i) - y.at(i);
        dist += dx * dx;
    }

    return dist;
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
    return x.isAllFinite();
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


}; /* namespace oofem::vec */

