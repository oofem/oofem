/* $Header: /home/cvs/bp/oofem/oofemlib/src/ldltfact.h,v 1.10 2003/04/06 14:08:24 bp Exp $ */
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


//   ******************************************
//   *** CLASS LDL^T Factorization          ***
//   ******************************************


#ifndef dsssolver_h
#define dsssolver_h

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif
#include "sparselinsystemnm.h"
#include "sparsemtrx.h"
#include "flotarry.h"

namespace oofem {

class Domain;
class EngngModel;
class FloatMatrix;

/**
 * Implements the solution of linear system of equation in the form Ax=b using direct factorization method.
 * Can work with any sparse matrix implementation. Hovewer, the sparse matrix implementation have to support
 * its factorization (canBeFactorized method).
 */
class DSSSolver : public SparseLinearSystemNM
{
private:

public:
    /// Constructor - creates new instance of LDLTFactorization, with number i, belonging to domain d and Engngmodel m.
    DSSSolver(int i, Domain *d, EngngModel *m);
    /// Destructor
    ~DSSSolver();    // destructor

    /**
     * Solves the given linear system by LDL^T factorization.
     * Implementation rely on factorization support provided by mapped sparse matrix.
     * It calls Lhs->factorized()->backSubstitutionWith(*solutionArray). Sets solved flag to 1 if o.k.
     * @param A coefficient matrix
     * @param b right hand side
     * @param x solution array
     * @return NM_Status value
     * @param tNow time step
     */
    NM_Status solve(SparseMtrx *A, FloatArray *b, FloatArray *x);
    /// Initializes receiver from given record. Empty implementation.
    IRResultType initializeFrom(InputRecord *ir);

    // identification
    /// Returns "LDLTFactorization" - class name of the receiver.
    const char *giveClassName() const { return "LDLTFactorization"; }
    /// Returns LDLTFactorizationClass - classType id of receiver.
    classType giveClassID() const { return LDLTFactorizationClass; }
    LinSystSolverType giveLinSystSolverType() const { return ST_DSS; }
};

} // end namespace oofem
#endif // dsssolver_h









