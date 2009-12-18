/* $Header: /home/cvs/bp/oofem/oofemlib/src/precond.h,v 1.7 2003/04/06 14:08:25 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   ****************************
//   *** CLASS PRECONDITIONER ***
//   ****************************


#ifndef precond_h
#define precond_h

#include "flotarry.h"
#include "sparsemtrx.h"
#include "inputrecord.h"
#ifndef __MAKEDEPEND
#include <string.h>
#endif

namespace oofem {

/**
 * Abstract class for IML++ compatible preconditioner.
 * Each preconditioner provides solve() and transpose_solve() functionality,
 * so that they can be used interchangeably in the same base iterative method code.
 *
 * Preconditioner matrix M is typically used to compute M^{-1}x or (M^T)^{-1}x during the
 * course of a basic iterartion, and thus can be seen as taking some input vector
 * and return a corresponding vector.
 */
class Preconditioner
{
public:
    /**
     * Constructor. Initializes the the receiver (constructs the precontioning matrix M) of given matrix.
     * Calls virtual init service.
     * @param a sparse matrix to be preconditioned
     * @param attributes attributes of receiver
     */
    Preconditioner(const SparseMtrx &a, InputRecord &attributes);
    /**
     * Constructor. The user should call initializeFrom and init services in this given order to ensure consistency.
     */
    Preconditioner() { }
    /// Destructor
    virtual ~Preconditioner(void) { };

    /**
     * Initializes the receiver (constructs the precontioning matrix M) of given matrix.
     * Virtual service, to be implemented by derived classes. Should be called after initializeFrom service.
     * @param a sparse matrix to be preconditioned
     */
    virtual void init(const SparseMtrx &a) { };

    /// Solves the linear system
    FloatArray solve(const FloatArray &x) const { FloatArray y;
                                                  this->solve(x, y);
                                                  return y; }
    /// Solves transposed system
    FloatArray trans_solve(const FloatArray &x) const { FloatArray y;
                                                        this->trans_solve(x, y);
                                                        return y; }
    /// Solves the linear system
    virtual void solve(const FloatArray &x, FloatArray &y) const = 0;
    /// Solves transposed system
    virtual void trans_solve(const FloatArray &x, FloatArray &y) const = 0;
    /// returns the preconditioner name
    virtual const char *giveClassName() const { return "Preconditioner"; }
    /// Initializes receiver from given record. Empty implementation.
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
};

} // end namespace oofem
#endif // precond_h
