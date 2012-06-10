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


#include "gjacobi.h"
#include "mathfem.h"
#include "flotmtrx.h"
#include "nmstatus.h"

namespace oofem {
GJacobi :: GJacobi(int i, Domain *d, EngngModel *m) :
    NumericalMethod(i, d, m) {
    //
    // constructor
    //
    //a      = NULL ;
    //b      = NULL ;
    //eigv   = NULL ;  // only pointer to caller data, not ownership
    //x      = NULL ;  // only pointer to caller data, not ownership

    nsmax  = 15;     // default maximum number of sweeps allowed
    rtol   = 10.E-12; // convergence tolerance
    n = 0;
    solved = 0;
}

GJacobi ::  ~GJacobi() { }


#define GJacobi_ZERO_CHECK_TOL 1.e-40

NM_Status
GJacobi :: solve(FloatMatrix *a, FloatMatrix *b, FloatArray *eigv, FloatMatrix *x)
//void  GJacobi :: solveYourselfAt (TimeStep* tNow)
//
// this function solve the generalized eigenproblem using the Generalized
// jacobi iteration
//
//
{
    int i, j, k, nsweep, nr, jj;
    double eps, eptola, eptolb, akk, ajj, ab, check, sqch, d1, d2, den, ca, cg, ak, bk, xj, xk;
    double aj, bj, tol, dif, epsa, epsb, bb;
    int jm1, kp1, km1, jp1;


    // first check whether Amatrix is defined
    if ( !a ) {
        _error("solveYourselfAt: unknown A matrix");
    }

    // and whether Bmatrix
    if ( !b ) {
        _error("solveYourselfAt: unknown Bmatrx");
    }

    if ( ( a->giveNumberOfRows() != b->giveNumberOfRows() ) ||
        ( !a->isSquare() ) || ( !b->isSquare() ) ) {
        _error("solveYourselfAt: A matrix, B mtrix -> size mismatch");
    }

    int n = a->giveNumberOfRows();
    //
    // Check output  arrays
    //
    if ( !eigv ) {
        _error("solveYourselfAt: unknown eigv array");
    }

    if ( !x ) {
        _error("solveYourselfAt: unknown x    mtrx ");
    }

    if ( eigv->giveSize() != n ) {
        _error("solveYourselfAt: eigv size mismatch");
    }

    if ( ( !x->isSquare() ) || ( x->giveNumberOfRows() != n ) ) {
        _error("solveYourselfAt: x size mismatch");
    }

    //
    // Create temporary arrays
    //
    FloatArray *d = new FloatArray(n);
    //
    // Initialize EigenValue and EigenVector Matrices
    //
    for ( i = 1; i <= n; i++ ) {
        //      if((a->at(i,i) <= 0. ) && (b->at(i,i) <= 0.))
        //        _error ("solveYourselfAt: Matrices are not positive definite");
        d->at(i) = a->at(i, i) / b->at(i, i);
        eigv->at(i) = d->at(i);
    }

    for ( i = 1; i <= n; i++ ) {
        for ( j = 1; j <= n; j++ ) {
            x->at(i, j) = 0.0;
        }

        x->at(i, i) = 1.0;
    }

    if ( n == 1 ) {
        return NM_Success;
    }

    //
    // Initialize sweep counter and begin iteration
    //
    nsweep = 0;
    nr = n - 1;

    do {
        nsweep++;
# ifdef DETAILED_REPORT
        OOFEM_LOG_DEBUG("*GJacobi*: sweep number %4d\n", nsweep);
#endif
        //
        // check if present off-diagonal element is large enough to require zeroing
        //
        eps = __OOFEM_POW(0.01, ( double ) nsweep);
        eps *= eps;
        for ( j = 1; j <= nr; j++ ) {
            jj = j + 1;
            for ( k = jj; k <= n; k++ ) {
                eptola = ( a->at(j, k) * a->at(j, k) ) / ( a->at(j, j) * a->at(k, k) );
                eptolb = ( b->at(j, k) * b->at(j, k) ) / ( b->at(j, j) * b->at(k, k) );
                if ( ( eptola  < eps ) && ( eptolb < eps ) ) {
                    continue;
                }

                //
                // if zeroing is required, calculate the rotation matrix elements ca and cg
                //
                akk = a->at(k, k) * b->at(j, k) - b->at(k, k) * a->at(j, k);
                ajj = a->at(j, j) * b->at(j, k) - b->at(j, j) * a->at(j, k);
                ab  = a->at(j, j) * b->at(k, k) - a->at(k, k) * b->at(j, j);
                check = ( ab * ab + 4.0 * akk * ajj ) / 4.0;
                if ( fabs(check) < GJacobi_ZERO_CHECK_TOL ) {
                    check = fabs(check);
                } else if ( check < 0.0 ) {
                    _error("solveYourselfAt: Matrices are not positive definite");
                }

                sqch = sqrt(check);
                d1 = ab / 2. + sqch;
                d2 = ab / 2. - sqch;
                den = d1;
                if ( fabs(d2) > fabs(d1) ) {
                    den = d2;
                }

                if ( den != 0.0 ) {                  // strange !
                    ca = akk / den;
                    cg = -ajj / den;
                } else {
                    ca = 0.0;
                    cg = -a->at(j, k) / a->at(k, k);
                }

                //
                // perform the generalized rotation to zero
                //
                if ( ( n - 2 ) != 0 ) {
                    jp1 = j + 1;
                    jm1 = j - 1;
                    kp1 = k + 1;
                    km1 = k - 1;
                    if ( ( jm1 - 1 ) >= 0 ) {
                        for ( i = 1; i <= jm1; i++ ) {
                            aj = a->at(i, j);
                            bj = b->at(i, j);
                            ak = a->at(i, k);
                            bk = b->at(i, k);
                            a->at(i, j) = aj + cg * ak;
                            b->at(i, j) = bj + cg * bk;
                            a->at(i, k) = ak + ca * aj;
                            b->at(i, k) = bk + ca * bj;
                        }
                    }

                    if ( ( kp1 - n ) <= 0 ) {
                        for ( i = kp1; i <= n; i++ ) { // label 140
                            aj = a->at(j, i);
                            bj = b->at(j, i);
                            ak = a->at(k, i);
                            bk = b->at(k, i);
                            a->at(j, i) = aj + cg * ak;
                            b->at(j, i) = bj + cg * bk;
                            a->at(k, i) = ak + ca * aj;
                            b->at(k, i) = bk + ca * bj;
                        }
                    }

                    if ( ( jp1 - km1 ) <= 0 ) { // label 160
                        for ( i = jp1; i <= km1; i++ ) {
                            aj = a->at(j, i);
                            bj = b->at(j, i);
                            ak = a->at(i, k);
                            bk = b->at(i, k);
                            a->at(j, i) = aj + cg * ak;
                            b->at(j, i) = bj + cg * bk;
                            a->at(i, k) = ak + ca * aj;
                            b->at(i, k) = bk + ca * bj;
                        }
                    }
                }                           // label 190

                ak = a->at(k, k);
                bk = b->at(k, k);
                a->at(k, k) = ak + 2.0 *ca *a->at(j, k) + ca *ca *a->at(j, j);
                b->at(k, k) = bk + 2.0 *ca *b->at(j, k) + ca *ca *b->at(j, j);
                a->at(j, j) = a->at(j, j) + 2.0 *cg *a->at(j, k) + cg * cg * ak;
                b->at(j, j) = b->at(j, j) + 2.0 *cg *b->at(j, k) + cg * cg * bk;
                a->at(j, k) = 0.0;
                b->at(j, k) = 0.0;
                //
                // update the eigenvector matrix after each rotation
                //
                for ( i = 1; i <= n; i++ ) {
                    xj = x->at(i, j);
                    xk = x->at(i, k);
                    x->at(i, j) = xj + cg * xk;
                    x->at(i, k) = xk + ca * xj;
                }                        // label 200

            }
        }                                // label 210

        //
        // update the eigenvalues after each sweep
        //
#ifdef DETAILED_REPORT
        OOFEM_LOG_DEBUG("GJacobi: a,b after sweep\n");
        a->printYourself();
        b->printYourself();
#endif
        for ( i = 1; i <= n; i++ ) {
            // in original uncommented
            //      if ((a->at(i,i) <= 0.) || (b->at(i,i) <= 0.))
            //        error ("solveYourselfAt: Matrices are not positive definite");
            eigv->at(i) = a->at(i, i) / b->at(i, i);
        }                                          // label 220

# ifdef DETAILED_REPORT
        OOFEM_LOG_DEBUG("GJacobi: current eigenvalues are:\n");
        eigv->printYourself();
        OOFEM_LOG_DEBUG("GJacobi: current eigenvectors are:\n");
        x->printYourself();
# endif
        //
        // check for convergence
        //
        for ( i = 1; i <= n; i++ ) {       // label 230
            tol = rtol * d->at(i);
            dif = ( eigv->at(i) - d->at(i) );
            if ( fabs(dif) > tol ) {
                goto label280;
            }
        }                                 // label 240

        //
        // check all off-diagonal elements to see if another sweep is required
        //
        eps = rtol * rtol;
        for ( j = 1; j <= nr; j++ ) {
            jj = j + 1;
            for ( k = jj; k <= n; k++ ) {
                epsa = ( a->at(j, k) * a->at(j, k) ) / ( a->at(j, j) * a->at(k, k) );
                epsb = ( b->at(j, k) * b->at(j, k) ) / ( b->at(j, j) * b->at(k, k) );
                if ( ( epsa < eps ) && ( epsb < eps ) ) {
                    continue;
                }

                goto label280;
            }
        }                                 // label 250

        //
        // fill out bottom triangle of resultant matrices and scale eigenvectors
        //
        break;
        // goto label255 ;
        //
        // update d matrix and start new sweep, if allowed
        //
label280:
        for ( i = 1; i <= n; i++ ) {
            d->at(i) = eigv->at(i);
        }
    } while ( nsweep < nsmax );

    // label255:
    for ( i = 1; i <= n; i++ ) {
        for ( j = 1; j <= n; j++ ) {
            a->at(j, i) = a->at(i, j);
            b->at(j, i) = b->at(i, j);
        }                               // label 260

    }

    for ( j = 1; j <= n; j++ ) {
        bb = sqrt( fabs( b->at(j, j) ) );
        for ( k = 1; k <= n; k++ ) {
            x->at(k, j) /= bb;
        }
    }                                  // label 270

    solved = 1;
    delete d;
    return NM_Success;
}

IRResultType
GJacobi :: initializeFrom(InputRecord *ir)
{
    return IRRT_OK;
}
} // end namespace oofem
