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

#include "sparselinsystemnm.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "sparsemtrx.h"

namespace oofem {
SparseLinearSystemNM :: SparseLinearSystemNM(int i, Domain *d, EngngModel *m) : NumericalMethod(i, d, m)
{ }

SparseLinearSystemNM :: ~SparseLinearSystemNM()
{ }

NM_Status SparseLinearSystemNM :: solve(SparseMtrx *A, FloatMatrix &B, FloatMatrix &X)
{
    NM_Status status;
    int ncol = A->giveNumberOfRows();
    int nrhs = B.giveNumberOfColumns();
    if (A->giveNumberOfRows() != B.giveNumberOfRows()) {
        OOFEM_ERROR("SparseLinearSystemNM :: solve - A and B matrix mismatch");
    }
    FloatArray bi(ncol), xi(ncol);
    X.resize(ncol,nrhs);
    for (int i = 1; i <= nrhs; ++i ) {
        B.copyColumn(bi, i);
        status = this->solve(A, &bi, &xi);
        if (status == NM_NoSuccess) {
            return NM_NoSuccess;
        }
        X.setColumn(xi, i);
    }
    return status;
}

} // end namespace oofem
