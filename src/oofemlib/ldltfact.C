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

#include "ldltfact.h"

namespace oofem {
LDLTFactorization :: LDLTFactorization(int i, Domain *d, EngngModel *m) :
    SparseLinearSystemNM(i, d, m) {
    //
    // constructor
    //
}

LDLTFactorization ::  ~LDLTFactorization() {
    //
    // destructor
    //
}

NM_Status
LDLTFactorization :: solve(SparseMtrx *A, FloatArray *b, FloatArray *x)
{
    int size;

    // first check whether Lhs is defined
    if ( !A ) {
        _error("solveYourselfAt: unknown Lhs");
    }

    // and whether Rhs
    if ( !b ) {
        _error("solveYourselfAt: unknown Rhs");
    }

    // and whether previous Solution exist
    if ( !x ) {
        _error("solveYourselfAt: unknown solution array");
    }

    if ( ( size = x->giveSize() ) != b->giveSize() ) {
        _error("solveYourselfAt: size mismatch");
    }

    // check whether Lhs supports factorization
    if ( !A->canBeFactorized() ) {
        _error("solveYourselfAt: Lhs not support factorization");
    }

    for ( int i = 1; i <= size; i++ ) {
        x->at(i) = b->at(i);
    }

    // solving
    A->factorized()->backSubstitutionWith(* x);

    return NM_Success;
}

IRResultType
LDLTFactorization :: initializeFrom(InputRecord *ir)
//
//
//
{
    return IRRT_OK;
}
} // end namespace oofem
