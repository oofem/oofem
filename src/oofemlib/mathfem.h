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

/*
 * This file contains some functions used in the finite element application.
 * ref : Lippman p 104
 */
#ifndef mathfem_h
#define mathfem_h

#include "error.h"
#include "compiler.h"

#ifndef __MAKEDEPEND
 #include <cmath>
 #include <cfloat> // For _isnan
#endif
#include "error.h"

namespace oofem {

class FloatArray;

#ifndef HAVE_M_PI
 #define M_PI 3.1415926535897932384626433832795029L
#endif
#ifndef HAVE_M_LN2
 #define M_LN2 0.6931471805599453094172321214581766L
#endif

/// Returns smaller value from two given decimals
inline int min(int i, int j)
{ return ( i <= j ? i : j ); }

/// Returns smaller value from two given long decimals
inline long min(long i, long j)
{ return ( i <= j ? i : j ); }

/// Returns smaller value from two given floats
inline double min(double i, double j)
{ return ( i <= j ? i : j ); }

/// Returns bigger value form two given decimals
inline int max(int i, int j)
{ return ( i >= j ? i : j ); }

/// Returns the clamped value of a between upper and lower
inline double clamp(int a, int lower, int upper)
{ return ( a <= lower ? lower : ( a >= upper ? upper : a) ); }

/// Returns bigger value form two given long decimals
inline long max(long i, long j)
{ return ( i >= j ? i : j ); }

/// Returns bigger value form two given floats
inline double max(double i, double j)
{ return ( i >= j ? i : j ); }

/// Returns the clamped value of a between upper and lower
inline double clamp(double a, double lower, double upper)
{ return ( a <= lower ? lower : ( a >= upper ? upper : a) ); }

/// Returns the signum of given value (if value is < 0 returns -1, otherwise returns 1)
inline double sgn(double i)
{ return ( i < 0. ? -1. : 1. ); }

#ifndef HAVE_ISNAN
#ifdef _MSC_VER
/// Returns true is x is NaN
inline bool isnan(double x) { return _isnan(x); }
#endif
#endif

#ifndef HAVE_NEAREST
/// Returns the nearest integer
inline int nearest(double x) { return (int)floor( x + 0.5 ); }
#endif

#ifndef HAVE_CBRT
/// Returns the cubic root of x.
inline double cbrt(double x) { return sgn(x)*pow(fabs(x),1.0/3.0); }
#endif

/// Returns the positive part of given float
inline double macbra(double x) { return ( x >= 0 ? x : 0 ); }
/// Returns the negative part of given float
inline double negbra(double x) { return ( x <= 0 ? x : 0 ); }

/**
 * Solves cubic equation for real roots.
 * The coefficients a to d gives the equation @f$ a x^3 + b x^2 + c x + d = 0@f$.
 * @param a Coefficient
 * @param b Coefficient
 * @param c Coefficient
 * @param d Coefficient
 * @param r1 First root
 * @param r2 Second root
 * @param r3 Third root
 * @param num Number of roots resolved (only first num roots are valid).
 */
void cubic(double a, double b, double c, double d, double *r1, double *r2, double *r3, int *num);

/**
 * Solves cubic equation for real roots, assuming that if cubic polynomial given then the only possibility
 * is that only three real roots exists. But also accepts cubic coefficient degenerated to
 * quadratic or linear equation.
 * This is used by algorithms for computing principal strain/stresses to
 * overcome rounding errors.
 * The coefficients a to d gives the equation @f$ a x^3 + b x^2 + c x + d = 0@f$.
 * @param a Coefficient
 * @param b Coefficient
 * @param c Coefficient
 * @param d Coefficient
 * @param r1 First root.
 * @param r2 Second root.
 * @param r3 Third root.
 * @param num Number of roots resolved (only first num roots are valid).
 */
void cubic3r(double a, double b, double c, double d, double *r1, double *r2, double *r3, int *num);

/**
 * Returns iperm of val, in specific rank
 */
int iperm(int val, int rank);


#define MATHFEM_C 0.38196601
#define MATHFEM_R ( 1 - MATHFEM_C )
#define MATHFEM_BRENT_MAXITER 100

template< class T >class mem_fun
{
    double ( T :: *pmf )( double );
    T *ptr;
public:
    mem_fun( T *o, double( T :: *p )( double ) ) : pmf(p), ptr(o) { }
    double operator()(double x) const { return ( ptr->*pmf )(x); }
};


class c_fun
{
    double ( *func )(double);
public:
    c_fun( double( * p )( double ) ) : func(p) { }
    double operator()(double x) const { return ( * func )( x ); }
};



/**
 * Minimize function of one variable using golden section search
 *
 * golden section search routine for finding the minimum of given function represented by functor f.
 * Input parameters:
 * ax, bx, cx -> three x-coordinates bracketing the minima (ax < bx < cx and f(bx) < f(ax) and f(bx) < f(cx))
 * tol - tolerance
 * Output parameters:
 * xmin coordinate of minima
 * return value - the minimum found
 *
 * Done according to Scientific Computation WS 2001/2002 by Gaston Gonnet
 * http://linneus20.ethz.ch:8080/wsrscript.html
 */
template< class T >double gss(double ax, double bx, double cx, const T &f,
                              double tol, double &xmin)
{
    int ii = 0;
    double f1, f2, x0, x1, x2, x3;

    x0 = ax;
    x3 = cx;

    // initialization x0-x1 made the smaller segment
    if ( fabs(cx - bx) > fabs(bx - ax) ) {
        x1 = bx;
        x2 = bx + MATHFEM_C * ( cx - bx );
    } else {
        x2 = bx;
        x1 = bx - MATHFEM_C * ( bx - ax );
    }

    f1 = f(x1);
    f2 = f(x2);

    // iteration loop
    while ( fabs(x3 - x0) > tol * ( fabs(x1) + fabs(x2) ) ) {
        if ( f2 < f1 ) {
            // minimum bracketed by (x1,x2,x3) triplet
            x0 = x1;
            x1 = x2;
            x2 = MATHFEM_R * x1 + MATHFEM_C * x3; // x2=x1+C*(x3-x1)

            f1 = f2;
            f2 = f(x2);
        } else {
            // minimum bracketed by (x0,x1,x2) triplet
            x3 = x2;
            x2 = x1;
            x1 = MATHFEM_R * x2 + MATHFEM_C * x0; // x1 = x2+c*(x0-x2)

            f2 = f1;
            f1 = f(x1);
        }

        ii++;
    }

    //printf ("gss: convergence reached in %d iterations\n", ii);
    if ( f1 < f2 ) {
        xmin = x1;
        return f1;
    } else {
        xmin = x2;
        return f2;
    }
}

template< class T >double brent(double ax, double bx, double cx, const T &f,
                                double tol, double &xmin)
{
    int ii;
    double x_left = ax, x_right = cx;
    double x, x_midpoint, v, w, u, tol1, tol2, p, q, r, e_tmp, d = 0.0, fx, fv, fw, fu;
    double e = 0.0;

    x = v = w = bx;
    fx = fv = fw = f(x);

    for ( ii = 1; ii <= MATHFEM_BRENT_MAXITER; ii++ ) {
        x_midpoint = 0.5 * ( x_left + x_right );

        // check for convergence here
        tol1 = tol * fabs(x) + 1.0e-10;
        tol2 = 2.0 * tol1;
        if ( fabs(x - x_midpoint) <= ( tol2 - 0.5 * ( x_right - x_left ) ) ) {
            //printf ("brent: convergence in %d iterations\n", ii);
            xmin = x;
            return fx;
        }

        if ( fabs(e) > tol1 ) {
            /* fit parabola */
            r = ( x - w ) * ( fx - fv );
            q = ( x - v ) * ( fx - fw );
            p = ( x - v ) * q - ( x - w ) * r;
            q = 2.0 * ( q - r );

            if ( q > 0 ) {
                p = -p;
            } else {
                q = -q;
            }

            e_tmp = e;
            e = d;

            if ( fabs(p) < fabs(0.5 * q * e_tmp) && p < q * ( x - x_left ) && p < q * ( x_right - x ) ) {
                d = p / q;
                u = x + d;
                if ( ( u - x_left ) < tol2 || ( x_right - u ) < tol2 ) {
                    d = ( x < x_midpoint ) ? tol1 : -tol1;
                }
            } else {
                e = ( x < x_midpoint ) ? x_right - x : -( x - x_left );
                d = MATHFEM_C * e;
            }
        } else {
            e = ( x < x_midpoint ) ? x_right - x : -( x - x_left );
            d = MATHFEM_C * e;
        }

        if ( fabs(d) >= tol1 ) {
            u = x + d;
        } else {
            u = x + ( ( d > 0 ) ? tol1 : -tol1 );
        }

        fu = f(u);

        if ( fu <= fx ) {
            if ( u >= x ) {
                x_left = x;
            } else {
                x_right = x;
            }

            v = w;
            w = x;
            x = u;
            fv = fw;
            fw = fx;
            fx = fu;
        } else {
            if ( u < x ) {
                x_left = u;
            } else {
                x_right = u;
            }

            if ( fu <= fw || w == x ) {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else if ( fu <= fv || v == x || v == w ) {
                v = u;
                fv = fu;
            }
        }
    }

    // too many iterations
    OOFEM_WARNING("brent : too many iterations\n");
    xmin = x;
    return fx;
}

/**
 * Least-square fit of 2nd degree polynomial @f$ y = a_0 + a_1 x + a_2 x^2 @f$.
 * @param x X-values.
 * @param y Y-values.
 * @param a Computed coefficients.
 */
void ls2fit(const FloatArray &x, const FloatArray &y, FloatArray &a);

} // end namespace oofem
#endif // mathfem_h
