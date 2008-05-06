/* $Header: /home/cvs/bp/oofem/oofemlib/src/sparsegeneigenvalsystemnm.h,v 1.3 2003/04/06 14:08:26 bp Exp $ */
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


//   *********************************************
//   *** CLASS SparseGeneralEigenValueSystemNM ***
//   *********************************************


#ifndef sparsegeneigenvalsystemnm_h
#define sparsegeneigenvalsystemnm_h


#include "nummet.h"

#include "classtype.h"
#include "nmstatus.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

class EngngModel;
class SparseMtrx;
class FloatArray;


/**
 * This base class is an abstraction for all numerical methods solving sparse
 * linear system of equations. The purpose of this class is to declare
 * the general interface to all numerical methods solving this kind of
 * problem. This interface allows to use any suitable
 * instance of the Numerical method class to the solve problem,
 * and leave the  whole engineering model code,
 * including mapping, unchanged, because all instances of this class
 * provide the common interface.
 */
class SparseGeneralEigenValueSystemNM : public NumericalMethod
{
protected:
public:
    /// Constructor
    SparseGeneralEigenValueSystemNM(int i, Domain *d, EngngModel *m);
    /// Destructor
    ~SparseGeneralEigenValueSystemNM();

    // identification
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "SparseGeneralEigenValueSystemNM"; }
    /** Returns classType id of receiver.
     * @see FEMComponent::giveClassID
     */
    classType giveClassID() const { return SparseGeneralEigenValueSystemNMClass; }

    /**
     * Solves the given sparse generalized eigen value system of equations Ax = o^2 Bx.
     * @param A coefficient matrix
     * @param B coefficient matrix
     * @param x eigen vector(s)
     * @param o eigen value(s)
     * @param rtol tolerance
     * @param nroot number of required eigenvalues
     * @return NM_Status value
     */
    virtual NM_Status solve(SparseMtrx *A, SparseMtrx *B, FloatArray *ev, FloatMatrix *x, double rtol, int nroot) = 0;


public:
};

#endif // sparsegeneigenvalsystemnm_h









