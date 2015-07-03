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

//#define DETAILED_REPORT

#include "subspaceit.h"
#include "engngm.h"
#include "domain.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "gjacobi.h"

namespace oofem {
SubspaceIteration :: SubspaceIteration(Domain *d, EngngModel *m) :
    SparseGeneralEigenValueSystemNM(d, m)
{
    //
    // constructor
    //
    //a      = NULL ;
    //b      = NULL ;
    //_eigv  = NULL ;   // not ownership
    //_r     = NULL ;   // not ownership

    //nroot  = 0 ;
    nc     = 0;
    //rtol   = 10.E-6 ;   // convergence tolerance
    n = 0;
    nitem = 40; // max number of iterations
    solved = 0;
}


SubspaceIteration :: ~SubspaceIteration()
{ }

NM_Status
SubspaceIteration :: solve(SparseMtrx &a, SparseMtrx &b, FloatArray &_eigv, FloatMatrix &_r, double rtol, int nroot)
//
// this function solve the generalized eigenproblem using the Generalized
// jacobi iteration
//
//
{
    FILE *outStream;
    FloatArray temp, w, d, tt, rtolv, eigv;
    FloatMatrix r;
    int nn, nc1, i, j, k, ij = 0, nite, is;
    double rt, art, brt, eigvt, dif;
    FloatMatrix ar, br, vec;

    GJacobi mtd(domain, engngModel);
    outStream = domain->giveEngngModel()->giveOutputStream();
    nc = min(2 * nroot, nroot + 8);
    //
    // check matrix size
    //
    if ( a.giveNumberOfColumns() != b.giveNumberOfColumns() ) {
        OOFEM_ERROR("matrices size mismatch");
    }

    // check matrix for factorization support
    if ( !a.canBeFactorized() ) {
        OOFEM_ERROR("a matrix not support factorization");
    }

    //
    // check for wery small problem
    //
    nn = a.giveNumberOfColumns();
    if ( nc > nn ) {
        nc = nn;
    }

    ar.resize(nc, nc);
    ar.zero();
    br.resize(nc, nc);
    br.zero();

    //
    // creation of initial iteration vectors
    //
    nc1 = nc - 1;

    w.resize(nn);
    w.zero();
    d.resize(nc);
    d.zero();
    tt.resize(nn);
    tt.zero();
    rtolv.resize(nc);
    rtolv.zero();
    vec.resize(nc, nc);
    vec.zero();                   // eigen vectors of reduced problem

    _r.resize(nn, nroot);
    _eigv.resize(nroot);

    //
    // create work arrays
    //
    r.resize(nn, nc);
    r.zero();
    eigv.resize(nc);
    eigv.zero();

    FloatArray h(nn);
    for ( i = 1; i <= nn; i++ ) {
        h.at(i) = 1.0;
        w.at(i) = b.at(i, i) / a.at(i, i);
    }

    b.times(h, tt);
    for ( i = 1; i <= nn; i++ ) {
        r.at(i, 1) = tt.at(i);
    }

    for ( j = 2; j <= nc; j++ ) {
        rt = 0.0;
        for ( i = 1; i <= nn; i++ ) {
            if ( fabs( w.at(i) ) >= rt ) {
                rt = fabs( w.at(i) );
                ij = i;
            }
        }

        tt.at(j) = ij;
        w.at(ij) = 0.;
        for ( i = 1; i <= nn; i++ ) {
            if ( i == ij ) {
                h.at(i) = 1.0;
            } else {
                h.at(i) = 0.0;
            }
        }

        b.times(h, tt);
        for ( i = 1; i <= nn; i++ ) {
            r.at(i, j) = tt.at(i);
        }
    } // (r = z)

# ifdef DETAILED_REPORT
    printf("SubspaceIteration :: solveYourselfAt: Degrees of freedom invoked by initial vectors :\n");
    tt.printYourself();
    printf("SubspaceIteration :: solveYourselfAt: initial vectors for iteration:\n");
    r.printYourself();
# endif

    //ish = 0;
    a.factorized();
    //
    // start of iteration loop
    //
    for ( nite = 0; ; ++nite ) {               // label 100
# ifdef DETAILED_REPORT
        printf("SubspaceIteration :: solveYourselfAt: Iteration loop no. %d\n", nite);
# endif
        //
        // compute projection ar and br of matrices a , b
        //
        for ( j = 1; j <= nc; j++ ) {
            for ( k = 1; k <= nn; k++ ) {
                tt.at(k) = r.at(k, j);
            }

            //a. forwardReductionWith(&tt) -> diagonalScalingWith (&tt)
            //  -> backSubstitutionWithtt) ;
            a.backSubstitutionWith(tt);

            for ( i = j; i <= nc; i++ ) {
                art = 0.;
                for ( k = 1; k <= nn; k++ ) {
                    art += r.at(k, i) * tt.at(k);
                }

                ar.at(j, i) = art;
            }

            for ( k = 1; k <= nn; k++ ) {
                r.at(k, j) = tt.at(k);                 // (r = xbar)
            }
        }

        ar.symmetrized();        // label 110
#ifdef DETAILED_REPORT
        printf("SubspaceIteration :: solveYourselfAt: Printing projection matrix ar\n");
        ar.printYourself();
#endif
        //
        for ( j = 1; j <= nc; j++ ) {
            for ( k = 1; k <= nn; k++ ) {
                tt.at(k) = r.at(k, j);
            }

            b.times(tt, temp);
            for ( i = j; i <= nc; i++ ) {
                brt = 0.;
                for ( k = 1; k <= nn; k++ ) {
                    brt += r.at(k, i) * temp.at(k);
                }

                br.at(j, i) = brt;
            }                   // label 180

            for ( k = 1; k <= nn; k++ ) {
                r.at(k, j) = temp.at(k);             // (r=zbar)
            }
        }                       // label 160

        br.symmetrized();
#ifdef DETAILED_REPORT
        printf("SubspaceIteration :: solveYourselfAt: Printing projection matrix br\n");
        br.printYourself();
#endif

        //
        // solution of reduced eigenvalue problem
        //
        mtd.solve(ar, br, eigv, vec);

        /// START EXPERIMENTAL
        // solve the reduced problem by Inverse iteration
        /*
         * {
         * FloatMatrix x(nc,nc), z(nc,nc), zz(nc,nc), arinv;
         * FloatArray  w(nc), ww(nc), tt(nc), t(nc);
         * double c;
         * int ii, i,j,k,ac;
         *
         *
         * //  initial setting
         * for (i=1;i<=nc;i++){
         * ww.at(i)=1.0;
         * }
         *
         * for (i=1;i<=nc;i++)
         * for (j=1; j<=nc;j++)
         *   z.at(i,j)=1.0;
         *
         * arinv.beInverseOf (ar);
         *
         * for (i=0;i<nitem;i++) {
         *
         * //  copy zz=z
         * zz = z;
         *
         * // solve matrix equation K.X = M.X
         * x.beProductOf (arinv, z);
         * //  evaluation of Rayleigh quotients
         * for (j=1;j<=nc;j++){
         *   w.at(j) = 0.0;
         *   for (k = 1; k<= nc; k++) w.at(j) += zz.at(k,j)*x.at(k,j);
         * }
         *
         * z.beProductOf (br, x);
         *
         * for (j=1;j<=nc;j++){
         *   c = 0;
         *   for (k = 1; k<= nc; k++) c+= z.at(k,j)*x.at(k,j);
         *   w.at(j) /= c;
         * }
         *
         * //  check convergence
         * ac=0;
         * for (j=1;j<=nc;j++){
         *   if (fabs((ww.at(j)-w.at(j))/w.at(j))< rtol)  ac++;
         *   ww.at(j)=w.at(j);
         * }
         *
         * printf ("\n iterace cislo  %d   %d",i,ac);
         * w.printYourself();
         *
         * //  Gramm-Schmidt ortogonalization
         * for (j=1;j<=nc;j++) {
         *   for (k = 1; k<= nc; k++) tt.at(k) = x.at(k,j) ;
         *   t.beProductOf(br,tt) ;
         *   for (ii=1;ii<j;ii++) {
         *     c = 0.0;
         *     for (k = 1; k<= nc; k++) c += x.at(k,ii)*t.at(k);
         *     for (k = 1; k<= nc; k++) x.at(k,j) -= x.at(k,ii)*c;
         *   }
         *   for (k = 1; k<= nc; k++) tt.at(k) = x.at(k,j) ;
         *   t.beProductOf (br,tt) ;
         *   c = 0.0;
         *   for (k = 1; k<= nc; k++) c += x.at(k,j)*t.at(k);
         *   for (k = 1; k<= nc; k++) x.at(k,j) /= sqrt(c);
         * }
         *
         * if (ac>nroot){
         *   break;
         * }
         *
         * //  compute new approximation of Z
         * z.beProductOf (br,x);
         * }
         *
         * eigv = w;
         * vec = x;
         * }
         * /// END EXPERIMANTAL
         */


        //
        // sorting eigenvalues according to their values
        //
        do {
            is = 0; // label 350
            for ( i = 1; i <= nc1; i++ ) {
                if ( fabs( eigv.at(i + 1) ) < fabs( eigv.at(i) ) ) {
                    is++;
                    eigvt = eigv.at(i + 1);
                    eigv.at(i + 1) = eigv.at(i);
                    eigv.at(i)   = eigvt;
                    for ( k = 1; k <= nc; k++ ) {
                        rt = vec.at(k, i + 1);
                        vec.at(k, i + 1) = vec.at(k, i);
                        vec.at(k, i)   = rt;
                    }
                }
            }                   // label 360
        } while ( is != 0 );

# ifdef DETAILED_REPORT
        printf("SubspaceIteration :: solveYourselfAt: current eigen values of reduced problem \n");
        eigv.printYourself();
        printf("SubspaceIteration :: solveYourselfAt: current eigen vectors of reduced problem \n");
        vec.printYourself();
# endif
        //
        // compute eigenvectors
        //
        for ( i = 1; i <= nn; i++ ) { // label 375
            for ( j = 1; j <= nc; j++ ) {
                tt.at(j) = r.at(i, j);
            }

            for ( k = 1; k <= nc; k++ ) {
                rt = 0.;
                for ( j = 1; j <= nc; j++ ) {
                    rt += tt.at(j) * vec.at(j, k);
                }

                r.at(i, k) = rt;
            }
        }                       // label 420   (r = z)

        //
        // convergency check
        //
        for ( i = 1; i <= nc; i++ ) {
            dif = ( eigv.at(i) - d.at(i) );
            rtolv.at(i) = fabs( dif / eigv.at(i) );
        }

# ifdef DETAILED_REPORT
        printf("SubspaceIteration :: solveYourselfAt: Reached precision of eigenvalues:\n");
        rtolv.printYourself();
# endif
        for ( i = 1; i <= nroot; i++ ) {
            if ( rtolv.at(i) > rtol ) {
                goto label400;
            }
        }

        fprintf(outStream,
                "SubspaceIteration :: solveYourselfAt: Convergence reached for RTOL=%20.15f",
                rtol);
        break;
label400:
        if ( nite >= nitem ) {
            fprintf(outStream, " SubspaceIteration :: solveYourselfAt: Convergence not reached in %d iteration - using current values", nitem);
            break;
        }

        for ( i = 1; i <= nc; i++ ) {
            d.at(i) = eigv.at(i);                     // label 410 and 440
        }

        continue;
    }


    // compute eigenvectors
    for ( j = 1; j <= nc; j++ ) {
        for ( k = 1; k <= nn; k++ ) {
            tt.at(k) = r.at(k, j);
        }

        a.backSubstitutionWith(tt);
        for ( k = 1; k <= nn; k++ ) {
            r.at(k, j) = tt.at(k);                   // (r = xbar)
        }
    }

    // one cad add a normalization of eigen-vectors here

    for ( i = 1; i <= nroot; i++ ) {
        _eigv.at(i) = eigv.at(i);
        for ( j = 1; j <= nn; j++ ) {
            _r.at(j, i) = r.at(j, i);
        }
    }

    solved = 1;
    return NM_Success;
}
} // end namespace oofem
