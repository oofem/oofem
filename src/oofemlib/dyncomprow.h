/* $Header: /home/cvs/bp/oofem/oofemlib/src/dyncomprow.h,v 1.3 2003/04/06 14:08:24 bp Exp $ */
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

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*        Compressed row sparse matrix  (O-based, Fortran)               */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


#ifndef dyncomprow_h
#define dyncomprow_h

#include "sparsemtrx.h"
#include "intarray.h"



// alloc chunk for columns
#define DynCompRow_CHUNK 8


/**
 * Implementation of sparse matrix stored in compressed column storage.
 * Designed to allow simple dynamic runtime grow of receiver.
 */
class DynCompRow : public SparseMtrx
{
protected:

    FloatArray **rows_;   // data values per column
    IntArray **colind_;   // row_ind per column
    IntArray diag_rowptr_; // pointers to the diagonal elements; needed only for ILU

    int base_;              // index base: offset of first element

public:
    /** Constructor. Before any operation an internal profile must be built.
     * @see builInternalStructure
     */
    DynCompRow(int n);
    /** Constructor. Before any operation an internal profile must be built.
     * @see builInternalStructure
     */
    DynCompRow();
    /// Copy constructor
    DynCompRow(const DynCompRow &S);
    /// Assignment operator
    DynCompRow &operator=(const DynCompRow &C);
    /// Destructor
    ~DynCompRow();
    /** Returns {\bf newly allocated} copy of receiver. Programmer must take
     * care about proper deallocation of allocated space.
     * @return newly allocated copy of receiver */
    SparseMtrx *GiveCopy() const;
    /** Evaluates a product of receiver with vector.
     * @param x array to be multiplied with receiver
     * @param answer result of product of receiver and x parameter
     */
    void times(const FloatArray &x, FloatArray &answer) const;
    /** Multiplies receiver by scalar value.
     * @param x value to multiply receiver
     */
    virtual void times(double x);
    /// Builds internal structure of receiver
    int buildInternalStructure(EngngModel *, int, EquationID, const UnknownNumberingScheme&);
    /** Assembles receiver from local element contributions.
     * @param loc location array. The values corresponding to zero loc array value are not assembled.
     * @param mat contribution to be assembled using loc array.
     */
    int assemble(const IntArray &loc, const FloatMatrix &mat);
    /** Assembles receiver from local element contributions.
     * @param rloc row location array. The values corresponding to zero loc array value are not assembled.
     * @param cloc column location array. The values corresponding to zero loc array value are not assembled.
     * @param mat contribution to be assembled using loc array.
     */
    int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat);

    /// Determines, whether receiver can be factorized.
    int canBeFactorized() const { return 0; }
    /// Zeroes the receiver.
    virtual SparseMtrx *zero();
    SparseMtrxType  giveType() const { return SMT_DynCompRow; }
    int isAntisymmetric() const { return 1; }
    virtual void printStatistics() const;

    /// Returns coefficient at position (i,j). 1-based element access
    virtual double &at(int i, int j);
    /// Returns coefficient at position (i,j). 1-based element access
    virtual double at(int i, int j) const;
    virtual void toFloatMatrix(FloatMatrix &answer) const;
    /// Prints receiver to stdout. Works only for relatively small matrices.
    virtual void printYourself() const;

    /** Performs LU factorization on yourself; modifies receiver
     * This routine computes the L and U factors of the ILU(p).
     */
    void ILUPYourself(int part_fill = 5, double drop_tol = 1.e-8);
    void ILUPsolve(const FloatArray &x, FloatArray &y) const;
    void ILUPtrans_solve(const FloatArray &x, FloatArray &y) const;

    /*******************************/
    /*  Access and info functions  */
    /*******************************/
    /// Returns col indx for i-th  row
    const IntArray *col_ind(int i) const { return colind_ [ i ]; }
    /// Returns row values
    const FloatArray *row(int i) const { return rows_ [ i ]; }

protected:

    /***********************************/
    /*  General access function (slow) */
    /***********************************/
    /// implements 0-based acess
    double operator()(int i, int j) const;
    /// implements 0-based acess
    double &operator()(int i, int j);

    /// returns the column index of given column at given row, else returns zero.
    int giveColIndx(int row, int col) const;
    /// insert column entry into row, preserving order of column indexes, returns the index of new row.
    int insertColInRow(int row, int col);

#ifdef IML_COMPAT
    /***********************************/
    /*  Matrix/Vector multiply         */
    /***********************************/

    FloatArray operator*(const FloatArray &x) const;
    FloatArray trans_mult(const FloatArray &x) const;

#endif

    void          checkSizeTowards(IntArray &);
    void          checkSizeTowards(const IntArray &rloc, const IntArray &cloc);
    void          growTo(int);

    void qsortRow(IntArray &ind, IntArray &ir, FloatArray &val, int l, int r);
    int qsortRowPartition(IntArray &ind, IntArray &ir, FloatArray &val, int l, int r);
};


#endif // dyncomprow_h

