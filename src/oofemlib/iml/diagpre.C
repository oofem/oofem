/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*             ********   ***                                 SparseLib++    */
/*          *******  **  ***       ***      ***               v. 1.5c        */
/*           *****      ***     ******** ********                            */
/*            *****    ***     ******** ********              R. Pozo        */
/*       **  *******  ***   **   ***      ***                 K. Remington   */
/*        ********   ********                                 A. Lumsdaine   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                     SparseLib++ : Sparse Matrix Library                   */
/*                                                                           */
/*               National Institute of Standards and Technology              */
/*                        University of Notre Dame                           */
/*              Authors: R. Pozo, K. Remington, A. Lumsdaine                 */
/*                                                                           */
/*                                 NOTICE                                    */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that the above notice appear in all copies and supporting        */
/* documentation.                                                            */
/*                                                                           */
/* Neither the Institutions (National Institute of Standards and Technology, */
/* University of Notre Dame) nor the Authors make any representations about  */
/* the suitability of this software for any purpose.  This software is       */
/* provided ``as is'' without expressed or implied warranty.                 */
/*                                                                           */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "diagpre.h"
#include "sparsemtrx.h"

namespace oofem {
DiagPreconditioner :: DiagPreconditioner(const SparseMtrx &C, InputRecord &attributes) : Preconditioner(C, attributes),
    diag( C.giveNumberOfRows() )
{ }


void
DiagPreconditioner :: init(const SparseMtrx &C)
{
    int n = C.giveNumberOfRows();

    diag.resize( n );
    diag.zero();

    /* Find the diagonal elements */
    for ( int i = 1; i <= n; i++ ) {
        double d = C.at(i, i);
        if ( d  == 0 ) {
            OOFEM_ERROR("failed, zero diagonal detected in equation %d", i);
        }

        diag.at(i) = 1. / d;
    }
}


void
DiagPreconditioner :: solve(const FloatArray &x, FloatArray &y) const
{
    y.resize( x.giveSize() );
    for ( int i = 0; i < x.giveSize(); i++ ) {
        y[i] = x[i] * diag[i];
    }
}


void
DiagPreconditioner :: trans_solve(const FloatArray &x, FloatArray &y) const
{
    y.resize( x.giveSize() );

    for ( int i = 0; i < x.giveSize(); i++ ) {
        y[i] = x[i] * diag[i];
    }
}
} // end namespace oofem
