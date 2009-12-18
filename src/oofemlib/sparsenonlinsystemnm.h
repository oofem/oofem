/* $Header: /home/cvs/bp/oofem/oofemlib/src/sparsenonlinsystemnm.h,v 1.6 2003/04/06 14:08:26 bp Exp $ */
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


//   *************************************
//   *** CLASS SparseNonLinearSystemNM ***
//   *************************************


#ifndef sparsenonlinsystemnm_h
#define sparsenonlinsystemnm_h


#include "nummet.h"

#include "equationid.h"
#include "nmstatus.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#endif

namespace oofem {

class EngngModel;
class SparseMtrx;
class FloatArray;


/**
 * This base class is an abstraction for all numerical methods solving sparse
 * nonlinear system of equations. The purpose of this class is to declare
 * the general interface to all numerical methods solving this kind of
 * problem. This interface allows to use any suitable
 * instance of the Numerical method class to the solve problem,
 * and leave the  whole engineering model code,
 * including mapping, unchanged, because all instances of this class
 * provide the common interface.
 */
class SparseNonLinearSystemNM : public NumericalMethod
{
public:
    /**
     * The following parameter allows to specify how the reference load vector
     * is obtained from given totalLoadVector and initialLoadVector.
     * The initialLoadVector desribes the part of loading which does not scale.
     * If refLoadInputMode is rlm_total (default) then the reference incremental load vector is defined as
     * totalLoadVector assembled at given time.
     * If refLoadInputMode is rlm_inceremental then the reference load vector is
     * obtained as incremental load vector at given time.
     */
    enum referenceLoadInputModeType { rlm_total=0, rlm_inceremental=1 };


protected:
    EquationID ut;
public:
    /// Constructor
    SparseNonLinearSystemNM(int i, Domain *d, EngngModel *m, EquationID ut);
    /// Destructor
    ~SparseNonLinearSystemNM();

    // identification
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "SparseNonLinearSystemNM"; }
    /** Returns classType id of receiver.
     * @see FEMComponent::giveClassID
     */
    classType giveClassID() const { return SparseNonLinearSystemNMClass; }

    /**
     * Solves the given sparse linear system of equations g(x,l)=l*R+R0-F(x); dx=K^{-1}g+ dl K^{-1}R.
     * Total load vector not passed, it is defined as l*R+R0, where l is scale factor
     * @param K coefficient matrix (K = dF/dx; stiffness matrix)
     * @param R  reference incremental Rhs (incremental load)
     * @param R0 initial Rhs (initial load)
     * @param Rr linearization of K*rri, where rri is increment of prescribed displacements
     * @param r  total solution (total displacement)
     * @param dr increment of solution (incremental displacaments)
     * @param l  Rhs scale factor (load level)
     * @param F  InternalRhs (real internal forces)
     * @param nite - number of iterations needed
     * @param rlm - reference load mode
     * @return NM_Status value
     */
    virtual NM_Status solve(SparseMtrx *K, FloatArray *R, FloatArray *R0,
                            FloatArray *Rr, FloatArray *r, FloatArray *dr, FloatArray *F,
                            double &l, referenceLoadInputModeType rlm,
                            int &nite, TimeStep *) = 0;

    virtual double giveCurrentStepLength() = 0;
    virtual void   setStepLength(double l) = 0;
    /**
     * Prints status mesage of receiver to output stream.
     * Prints the message corresponding to last solve.
     */
    virtual void   printState(FILE *outputStream) { }

public:
};

} // end namespace oofem
#endif // sparsenonlinsystemnm_h
