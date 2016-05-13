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
#include "sparselinsystemnm.h"
#include "classfactory.h"

namespace oofem {
SubspaceIteration :: SubspaceIteration(Domain *d, EngngModel *m) :
    SparseGeneralEigenValueSystemNM(d, m)
{
    nitem = 40; // max number of iterations
}


SubspaceIteration :: ~SubspaceIteration()
{ }

NM_Status
SubspaceIteration :: solve(SparseMtrx &a, SparseMtrx &b, FloatArray &_eigv, FloatMatrix &_r, double rtol, int nroot)
//
// this function solve the generalized eigenproblem using the Generalized
// jacobi iteration
//
{
    if ( a.giveNumberOfColumns() != b.giveNumberOfColumns() ) {
        OOFEM_ERROR("matrices size mismatch");
    }

    FloatArray temp, w, d, tt, f, rtolv, eigv;
    FloatMatrix r;
    int nn, nc1, ij = 0, is;
    double rt, art, brt, eigvt;
    FloatMatrix ar, br, vec;
    SparseLinearSystemNM *solver = GiveClassFactory().createSparseLinSolver(ST_Direct, domain, engngModel);

    GJacobi mtd(domain, engngModel);
    int nc = min(2 * nroot, nroot + 8);
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

    FloatMatrix tmp;
    a.toFloatMatrix(tmp);
    tmp.writeCSV("tmp.txt");
    b.toFloatMatrix(tmp);
    tmp.writeCSV("tmpb.txt");

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
    for ( int i = 1; i <= nn; i++ ) {
        h.at(i) = 1.0;
        w.at(i) = b.at(i, i) / a.at(i, i);
    }

    b.times(h, tt);
    r.setColumn(tt, 1);

    for ( int j = 2; j <= nc; j++ ) {
        rt = 0.0;
        for ( int i = 1; i <= nn; i++ ) {
            if ( fabs( w.at(i) ) >= rt ) {
                rt = fabs( w.at(i) );
                ij = i;
            }
        }

        tt.at(j) = ij;
        w.at(ij) = 0.;
        for ( int i = 1; i <= nn; i++ ) {
            if ( i == ij ) {
                h.at(i) = 1.0;
            } else {
                h.at(i) = 0.0;
            }
        }

        b.times(h, tt);
        r.setColumn(tt, j);
    } // (r = z)

# ifdef DETAILED_REPORT
    OOFEM_LOG_INFO("SubspaceIteration :: solveYourselfAt: Degrees of freedom invoked by initial vectors :\n");
    tt.printYourself();
    OOFEM_LOG_INFO("SubspaceIteration :: solveYourselfAt: initial vectors for iteration:\n");
    r.printYourself();
# endif

    //ish = 0;
    a.factorized();
    //
    // start of iteration loop
    //
    for ( int nite = 0; ; ++nite ) {               // label 100
# ifdef DETAILED_REPORT
        printf("SubspaceIteration :: solveYourselfAt: Iteration loop no. %d\n", nite);
# endif
        //
        // compute projection ar and br of matrices a , b
        //
        for ( int j = 1; j <= nc; j++ ) {
            f.beColumnOf(r, j);

            solver->solve(a, f, tt);

            for ( int i = j; i <= nc; i++ ) {
                art = 0.;
                for ( int k = 1; k <= nn; k++ ) {
                    art += r.at(k, i) * tt.at(k);
                }

                ar.at(j, i) = art;
            }

            r.setColumn(tt, j);            // (r = xbar)
        }

        ar.symmetrized();        // label 110
#ifdef DETAILED_REPORT
        OOFEM_LOG_INFO("SubspaceIteration :: solveYourselfAt: Printing projection matrix ar\n");
        ar.printYourself();
#endif
        //
        for ( int j = 1; j <= nc; j++ ) {
            tt.beColumnOf(r, j);

            b.times(tt, temp);
            for ( int i = j; i <= nc; i++ ) {
                brt = 0.;
                for ( int k = 1; k <= nn; k++ ) {
                    brt += r.at(k, i) * temp.at(k);
                }

                br.at(j, i) = brt;
            }                   // label 180

            r.setColumn(temp, j);        // (r=zbar)
        }                       // label 160

        br.symmetrized();
#ifdef DETAILED_REPORT
        OOFEM_LOG_INFO("SubspaceIteration :: solveYourselfAt: Printing projection matrix br\n");
        br.printYourself();
#endif

        //
        // solution of reduced eigenvalue problem
        //
        mtd.solve(ar, br, eigv, vec);

        // START EXPERIMENTAL
#if 0
        // solve the reduced problem by Inverse iteration
        {
            FloatMatrix x(nc,nc), z(nc,nc), zz(nc,nc), arinv;
            FloatArray  w(nc), ww(nc), tt(nc), t(nc);
            double c;

            //  initial setting
            for ( int i = 1;i <= nc; i++ ) {
                ww.at(i)=1.0;
            }
            
            
            for ( int i = 1;i <= nc; i++ )
                for ( int j = 1; j <= nc;j++ )
                    z.at(i,j)=1.0;
            
            arinv.beInverseOf (ar);
            
            for ( int i = 0;i < nitem; i++ ) {
                //  copy zz=z
                zz = z;
                
                // solve matrix equation K.X = M.X
                x.beProductOf(arinv, z);
                //  evaluation of Rayleigh quotients
                for ( int j = 1;j <= nc; j++ ) {
                    w.at(j) = 0.0;
                    for (k = 1; k<= nc; k++) w.at(j) += zz.at(k,j) * x.at(k,j);
                }

                z.beProductOf (br, x);

                for ( int j = 1;j <= nc; j++ ) {
                    c = 0;
                    for ( int k = 1; k<= nc; k++ ) c += z.at(k,j) * x.at(k,j);
                    w.at(j) /= c;
                }

                //  check convergence
                int ac = 0;
                for ( int j = 1;j <= nc; j++ ) {
                    if (fabs((ww.at(j)-w.at(j))/w.at(j))< rtol)  ac++;
                    ww.at(j) = w.at(j);
                }

                //printf ("\n iterace cislo  %d   %d",i,ac);
                //w.printYourself();

                //  Gramm-Schmidt ortogonalization
                for ( int j = 1;j <= nc;j++ ) {
                    for ( int k = 1; k<= nc; k++ ) tt.at(k) = x.at(k,j);
                    t.beProductOf(br,tt) ;
                    for ( int ii = 1;ii < j; ii++ ) {
                        c = 0.0;
                        for ( int k = 1; k<= nc; k++ ) c += x.at(k,ii) * t.at(k);
                        for ( int k = 1; k<= nc; k++ ) x.at(k,j) -= x.at(k,ii) * c;
                    }
                    for ( int k = 1; k<= nc; k++) tt.at(k) = x.at(k,j);
                    t.beProductOf(br, tt);
                    c = 0.0;
                    for ( int k = 1; k<= nc; k++) c += x.at(k,j)*t.at(k);
                    for ( int k = 1; k<= nc; k++) x.at(k,j) /= sqrt(c);
                }

                if ( ac > nroot ) {
                    break;
                }

                //  compute new approximation of Z
                z.beProductOf(br,x);
            }
            
            eigv = w;
            vec = x;
        }
#endif


        //
        // sorting eigenvalues according to their values
        //
        do {
            is = 0; // label 350
            for ( int i = 1; i <= nc1; i++ ) {
                if ( fabs( eigv.at(i + 1) ) < fabs( eigv.at(i) ) ) {
                    is++;
                    eigvt = eigv.at(i + 1);
                    eigv.at(i + 1) = eigv.at(i);
                    eigv.at(i)   = eigvt;
                    for ( int k = 1; k <= nc; k++ ) {
                        rt = vec.at(k, i + 1);
                        vec.at(k, i + 1) = vec.at(k, i);
                        vec.at(k, i)   = rt;
                    }
                }
            }                   // label 360
        } while ( is != 0 );

# ifdef DETAILED_REPORT
        OOFEM_LOG_INFO("SubspaceIteration :: solveYourselfAt: current eigen values of reduced problem \n");
        eigv.printYourself();
        OOFEM_LOG_INFO("SubspaceIteration :: solveYourselfAt: current eigen vectors of reduced problem \n");
        vec.printYourself();
# endif
        //
        // compute eigenvectors
        //
        for ( int i = 1; i <= nn; i++ ) { // label 375
            for ( int j = 1; j <= nc; j++ ) {
                tt.at(j) = r.at(i, j);
            }

            for ( int k = 1; k <= nc; k++ ) {
                rt = 0.;
                for ( int j = 1; j <= nc; j++ ) {
                    rt += tt.at(j) * vec.at(j, k);
                }

                r.at(i, k) = rt;
            }
        }                       // label 420   (r = z)

        //
        // convergency check
        //
        for ( int i = 1; i <= nc; i++ ) {
            double dif = ( eigv.at(i) - d.at(i) );
            rtolv.at(i) = fabs( dif / eigv.at(i) );
        }

# ifdef DETAILED_REPORT
        OOFEM_LOG_INFO("SubspaceIteration :: solveYourselfAt: Reached precision of eigenvalues:\n");
        rtolv.printYourself();
# endif
        for ( int i = 1; i <= nroot; i++ ) {
            if ( rtolv.at(i) > rtol ) {
                goto label400;
            }
        }

        OOFEM_LOG_INFO("SubspaceIteration :: solveYourselfAt: Convergence reached for RTOL=%20.15f", rtol);
        break;
label400:
        if ( nite >= nitem ) {
            OOFEM_WARNING("SubspaceIteration :: solveYourselfAt: Convergence not reached in %d iteration - using current values", nitem);
            break;
        }

        d = eigv;                     // label 410 and 440

        continue;
    }


    // compute eigenvectors
    for ( int j = 1; j <= nc; j++ ) {
        tt.beColumnOf(r, j);

        a.backSubstitutionWith(tt);
        r.setColumn(tt, j);                          // r = xbar
    }

    // one cad add a normalization of eigen-vectors here

    for ( int i = 1; i <= nroot; i++ ) {
        _eigv.at(i) = eigv.at(i);
        for ( int j = 1; j <= nn; j++ ) {
            _r.at(j, i) = r.at(j, i);
        }
    }

    return NM_Success;
}
} // end namespace oofem
