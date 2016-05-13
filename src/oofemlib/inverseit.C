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

#include "inverseit.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "sparsemtrx.h"
#include "mathfem.h"
#include "sparselinsystemnm.h"
#include "classfactory.h"

namespace oofem {

InverseIteration :: InverseIteration(Domain *d, EngngModel *m) :
    SparseGeneralEigenValueSystemNM(d, m)
{
    nitem = 100; // max number of iterations
}


InverseIteration :: ~InverseIteration() { }

NM_Status
InverseIteration :: solve(SparseMtrx &a, SparseMtrx &b, FloatArray &_eigv, FloatMatrix &_r, double rtol, int nroot)
{
    if ( a.giveNumberOfColumns() != b.giveNumberOfColumns() ) {
        OOFEM_ERROR("matrices size mismatch");
    }

    SparseLinearSystemNM *solver = GiveClassFactory().createSparseLinSolver(ST_Direct, domain, engngModel);
   
    int nn = a.giveNumberOfColumns();
    int nc = min(2 * nroot, nroot + 8);
    nc = min(nc, nn);

    FloatArray w(nc), ww(nc), t;
    std :: vector< FloatArray > z(nc, nn), zz(nc, nn), x(nc, nn);

    /*  initial setting  */
    ww.add(1.0);
    for ( int j = 0; j < nc; j++ ) {
        z[j].add(1.0);
    }

    for ( int i = 0; i < nitem; i++ ) {
        /*  copy zz=z  */
        for ( int j = 0; j < nc; j++ ) {
            zz[j] = z[j];
        }

        /*  solve matrix equation K.X = M.X  */
        for ( int j = 0; j < nc; j++ ) {
            solver->solve(a, z[j], x[j]);
        }

        /*  evaluation of Rayleigh quotients  */
        for ( int j = 0; j < nc; j++ ) {
            w.at(j + 1) = zz[j].dotProduct(x[j]);
        }

        //mmc_sky (m,x,z,adrm,n,niv);
        for ( int j = 0; j < nc; j++ ) {
            b.times(x[j], z[j]);
        }

        for ( int j = 0; j < nc; j++ ) {
            w.at(j + 1) /= z[j].dotProduct(x[j]);
        }

        /*  check convergence  */
        int ac = 0;
        for ( int j = 1; j <= nc; j++ ) {
            if ( fabs( ww.at(j) - w.at(j) ) <= fabs( w.at(j) * rtol ) ) {
                ac++;
            }

            ww.at(j) = w.at(j);
        }

        //printf ("\n iteration  %d   %d",i,ac);
        //w.printYourself();

        /*  Gramm-Schmidt ortogonalization   */
        for ( int j = 0; j < nc; j++ ) {
            if ( j != 0 ) {
                b.times(x[j], t);
            }

            for ( int ii = 0; ii < j; ii++ ) {
                x[j].add( -x[ii].dotProduct(t), x[ii] );
            }

            b.times(x[j], t);
            x[j].times( 1.0 / sqrt( x[j].dotProduct(t) ) );
        }

        if ( ac > nroot ) {
            break;
        }

        /*  compute new approximation of Z  */
        for ( int j = 0; j < nc; j++ ) {
            b.times(x[j], z[j]);
        }
    }

    // copy results
    _eigv.resize(nroot);
    _r.resize(nn, nroot);
    int i;
    for ( i = 1; i <= nroot; i++ ) {
        _eigv.at(i) = w.at(i);
        _r.setColumn(x[i - 1], i);
    }

    if ( i < nitem ) {
        OOFEM_LOG_INFO("InverseIt info: convergence reached in %d iterations\n", i);
    } else {
        OOFEM_WARNING("convergence not reached after %d iterations\n", i);
    }

    return NM_Success;
}
} // end namespace oofem
