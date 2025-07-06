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

#include "dyncomprow.h"
#include "ilucomprowprecond.h"
#include "verbose.h"

namespace oofem {
CompRow_ILUPreconditioner ::
CompRow_ILUPreconditioner(const SparseMtrx &A, InputRecord &attributes) : Preconditioner(A, attributes)
{ }

void
CompRow_ILUPreconditioner :: initializeFrom(InputRecord &ir)
{
    Preconditioner :: initializeFrom(ir);
    this->drop_tol = 1.e-8;
    IR_GIVE_OPTIONAL_FIELD(ir, this->drop_tol, _IFT_CompRow_ILUPrecond_droptol);

    part_fill = 5;
    IR_GIVE_OPTIONAL_FIELD(ir, part_fill, _IFT_CompRow_ILUPrecond_partfill);
}


void
CompRow_ILUPreconditioner :: init(const SparseMtrx &A)
{
    if ( dynamic_cast< const DynCompRow *>(& A) ) {
        this->A = static_cast< const DynCompRow &>(A);
        this->A.ILUPYourself(part_fill, drop_tol);
    } else {
        OOFEM_ERROR("unsupported sparse matrix type");
    }
}


void
CompRow_ILUPreconditioner :: solve(const FloatArray &x, FloatArray &y) const
{
    A.ILUPsolve(x, y);
}


void
CompRow_ILUPreconditioner :: trans_solve(const FloatArray &x, FloatArray &y) const
{
    A.ILUPtrans_solve(x, y);
}
} // end namespace oofem
