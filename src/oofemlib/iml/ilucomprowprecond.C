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

#include "dyncomprow.h"
#include "ilucomprowprecond.h"
#include "verbose.h"

#ifdef TIME_REPORT
 #include "clock.h"
#endif

namespace oofem {
CompRow_ILUPreconditioner ::
CompRow_ILUPreconditioner(const SparseMtrx &A, InputRecord &attributes) : Preconditioner(A, attributes)
{ }

IRResultType
CompRow_ILUPreconditioner :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    Preconditioner :: initializeFrom(ir);

    this->drop_tol = 1.e-8;
    IR_GIVE_OPTIONAL_FIELD(ir, this->drop_tol, IFT_CompRowPrecond_droptol, "droptol"); // Macro

    part_fill = 5;
    IR_GIVE_OPTIONAL_FIELD(ir, part_fill, IFT_CompRowPrecond_partfill, "partfill"); // Macro

    return IRRT_OK;
}


void
CompRow_ILUPreconditioner :: init(const SparseMtrx &A)
{
    if ( A.giveType() == SMT_DynCompRow ) {
        this->A =  ( * ( ( DynCompRow * ) & A ) );
        ( this->A ).ILUPYourself(part_fill, drop_tol);
    } else {
        OOFEM_ERROR("CompRow_ILUPreconditioner::init : unsupported sparse matrix type");
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
