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

#ifndef voidprecond_h
#define voidprecond_h

#include "floatarray.h"
#include "sparsemtrx.h"
#include "precond.h"

namespace oofem {
/**
 * Class implementing void preconditioner.
 */
class OOFEM_EXPORT VoidPreconditioner : public Preconditioner
{
public:
    /// Constructor. Creates the empty preconditioner.
    VoidPreconditioner(const SparseMtrx & a, InputRecord & attributes);
    /// Constructor. The user should call initializeFrom and init services in this given order to ensure consistency.
    VoidPreconditioner();
    /// Destructor
    virtual ~VoidPreconditioner(void) { }

    void init(const SparseMtrx &a) override { }
    void solve(const FloatArray &rhs, FloatArray &solution) const override { solution = rhs; }
    void trans_solve(const FloatArray &rhs, FloatArray &solution) const override { solution = rhs; }

    const char *giveClassName() const override { return "VoidPreconditioner"; }
};
} // end namespace oofem
#endif // voidprecond_h
