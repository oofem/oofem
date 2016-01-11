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

#include "sparselinsystemnm.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "sparsemtrx.h"

namespace oofem {
SparseLinearSystemNM :: SparseLinearSystemNM(Domain *d, EngngModel *m) : NumericalMethod(d, m)
{ }

SparseLinearSystemNM :: ~SparseLinearSystemNM()
{ }

NM_Status SparseLinearSystemNM :: solve(SparseMtrx &A, FloatMatrix &B, FloatMatrix &X)
{
    NM_Status status = NM_None;
    int ncol = A.giveNumberOfRows();
    int nrhs = B.giveNumberOfColumns();
    if ( A.giveNumberOfRows() != B.giveNumberOfRows() ) {
        OOFEM_ERROR("A and B matrix mismatch");
    }
    FloatArray bi(ncol), xi(ncol);
    X.resize(ncol, nrhs);
    for ( int i = 1; i <= nrhs; ++i ) {
        B.copyColumn(bi, i);
        status &= this->solve(A, bi, xi);
        if ( status & NM_NoSuccess ) {
            return NM_NoSuccess;
        }
        X.setColumn(xi, i);
    }
    return status;
}
} // end namespace oofem
