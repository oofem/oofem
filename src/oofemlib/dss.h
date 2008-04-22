/* $Header: /home/cvs/bp/oofem/oofemlib/src/symcompcol.h,v 1.3 2003/04/06 14:08:26 bp Exp $ */
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



#ifndef dss_h
#define dss_h

#include "sparsemtrx.h"
#include "intarray.h"


class DSSolver;
class SparseMatrixF;

/**
 * Interface to Direct Sparse Solver written by R.Vonracek.
 * This class represent the sparse matrix interface to DSS library. It allows to build internal structure,
 * assemble the DSS sparse matrix, and to factorize and back substitution operations.
 */
class DSSMatrix : public SparseMtrx
{
public:
    /// possible storage schemes and factorization types
    enum dssType { sym_LDL, sym_LL, unsym_LU };

protected:
    /// poiner  to SparseMatrixF class rep
    SparseMatrixF *_sm;
    /// pointer to DSSolver class (representation of Direct Sparse Solver in DSS lib)
    DSSolver *_dss;
    /// Flag indicating whether factorized.
    int isFactorized;
    /// type of storage & factorization
    dssType _type;

public:
    /** Constructor. Before any operation an internal profile must be built.
     *  @see builInternalStructure
     */
    DSSMatrix(dssType _t, int n);
    /** Constructor. Before any operation an internal profile must be built.
     *  @see builInternalStructure
     */
    DSSMatrix(dssType _t);
    /// Copy constructor
    DSSMatrix(const DSSMatrix &S);
    /// Destructor
    ~DSSMatrix();
    /** Returns {\bf newly allocated} copy of receiver. Programmer must take
     *  care about proper deallocation of allocated space.
     *  @return newly allocated copy of receiver */
    SparseMtrx *GiveCopy() const;
    /** Evaluates a product of receiver with vector.
     *  @param x array to be multiplied with receiver
     *  @param answer result of product of receiver and x parameter
     */
    void times(const FloatArray &x, FloatArray &answer) const;
    /** Multiplies receiver by scalar value.
     *  @param x value to multiply receiver
     */
    virtual void times(double x);
    /// Builds internal structure of receiver
    int buildInternalStructure(EngngModel *, int, EquationID);
    /** Assembles receiver from local element contributions.
     *  @param loc location array. The values corresponding to zero loc array value are not assembled.
     *  @param mat contribution to be assembled using loc array.
     */
    int assemble(const IntArray &loc, const FloatMatrix &mat);
    /** Assembles receiver from local element contributions.
     *  @param rloc row location array. The values corresponding to zero loc array value are not assembled.
     *  @param cloc column location array. The values corresponding to zero loc array value are not assembled.
     *  @param mat contribution to be assembled using loc array.
     */
    int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat);

    /// Determines, whether receiver can be factorized.
    int canBeFactorized() const { return 1; }
    /**
     * Returns the receiver factorized. \f$L^T D L\f$ form is used.
     * @return pointer to the receiver
     */
    virtual SparseMtrx *factorized();
    /**
     * Computes the solution of linear system \f$A x = \f$. A is receiver.
     */
    void solve(FloatArray *b, FloatArray *x);

    /// Zeroes the receiver.
    virtual SparseMtrx *zero();

    /// Returns coefficient at position (i,j). 1-based element access
    virtual double &at(int i, int j);
    /// Returns coefficient at position (i,j). 1-based element access
    virtual double at(int i, int j) const;

    SparseMtrxType  giveType() const { return SMT_SymCompCol; }
    int isAntisymmetric() const { return 0; }

    virtual void toFloatMatrix(FloatMatrix &answer) const { }
    virtual void printYourself() const { }


protected:


    /***********************************/
    /*  General access function (slow) */
    /***********************************/
    /// implements 0-based acess
    double operator()(int i, int j) const;
    /// implements 0-based acess
    double &operator()(int i, int j);

#ifdef IML_COMPAT
    /***********************************/
    /*  Matrix/Vector multiply         */
    /***********************************/

    FloatArray operator*(const FloatArray &x) const
    { FloatArray answer;
      this->times(x, answer);
      return answer; }
    FloatArray trans_mult(const FloatArray &x) const
    { FloatArray answer;
      this->times(x, answer);
      return answer; }

#endif
};


#endif // dss_h

