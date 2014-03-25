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

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*        Compressed row sparse matrix  (O-based, Fortran)               */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


#ifndef dyncomprow_h
#define dyncomprow_h

#include "sparsemtrx.h"
#include "intarray.h"

namespace oofem {
// alloc chunk for columns
#define DynCompRow_CHUNK 8


/**
 * Implementation of sparse matrix stored in compressed column storage.
 * Designed to allow simple dynamic runtime grow of receiver.
 */
class OOFEM_EXPORT DynCompRow : public SparseMtrx
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
    DynCompRow(const DynCompRow & S);
    /// Assignment operator
    DynCompRow &operator = ( const DynCompRow & C );
    /// Destructor
    virtual ~DynCompRow();

    // Overloaded methods:
    SparseMtrx *GiveCopy() const;
    void times(const FloatArray &x, FloatArray &answer) const;
    void timesT(const FloatArray &x, FloatArray &answer) const;
    virtual void times(double x);
    int buildInternalStructure(EngngModel *, int, EquationID, const UnknownNumberingScheme &);
    int assemble(const IntArray &loc, const FloatMatrix &mat);
    int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat);
    bool canBeFactorized() const { return false; }
    void zero();
    SparseMtrxType  giveType() const { return SMT_DynCompRow; }
    bool isAsymmetric() const { return true; }
    void printStatistics() const;
    double &at(int i, int j);
    double at(int i, int j) const;

    /** Performs LU factorization on yourself; modifies receiver
     * This routine computes the L and U factors of the ILU(p).
     */
    void ILUPYourself(int part_fill = 5, double drop_tol = 1.e-8);
    void ILUPsolve(const FloatArray &x, FloatArray &y) const;
    void ILUPtrans_solve(const FloatArray &x, FloatArray &y) const;

    /*******************************/
    /*  Access and info functions  */
    /*******************************/
    /// Returns col index for i-th  row
    const IntArray *col_ind(int i) const { return colind_ [ i ]; }
    /// Returns row values
    const FloatArray *row(int i) const { return rows_ [ i ]; }

protected:

    /***********************************/
    /*  General access function (slow) */
    /***********************************/
    /// implements 0-based access
    double operator() (int i, int j) const;
    /// implements 0-based access
    double &operator() (int i, int j);

    /// returns the column index of given column at given row, else returns zero.
    int giveColIndx(int row, int col) const;
    /// insert column entry into row, preserving order of column indexes, returns the index of new row.
    int insertColInRow(int row, int col);

    void checkSizeTowards(IntArray &);
    void checkSizeTowards(const IntArray &rloc, const IntArray &cloc);
    void growTo(int);

    void qsortRow(IntArray &ind, IntArray &ir, FloatArray &val, int l, int r);
    int qsortRowPartition(IntArray &ind, IntArray &ir, FloatArray &val, int l, int r);
};
} // end namespace oofem
#endif // dyncomprow_h
