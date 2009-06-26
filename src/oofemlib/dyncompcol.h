/* $Header: /home/cvs/bp/oofem/oofemlib/src/dyncompcol.h,v 1.3 2003/04/06 14:08:24 bp Exp $ */
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
/*        Compressed column sparse matrix  (O-based, Fortran)            */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


#ifndef dyncompcol_h
#define dyncompcol_h

#include "sparsemtrx.h"
#include "intarray.h"



// alloc chunk for columns
#define DynCompCol_CHUNK 8
// turns on the stl set container, used to store column data, usually slower
//#define  DynCompCol_USE_STL_SETS


#ifdef DynCompCol_USE_STL_SETS
#ifndef __MAKEDEPEND
#include <map>
#endif
#endif
/**
 * Implementation of sparse matrix stored in compressed column storage.
 * Designed to allow simple dynamic runtime grow of receiver.
 */
class DynCompCol : public SparseMtrx
{
protected:

#ifndef DynCompCol_USE_STL_SETS
    FloatArray * * columns_; // data values per column
    IntArray **rowind_;   // row_ind per column
#else
    std :: map< int, double > ** columns;
#endif

    int base_;              // index base: offset of first element

public:
    /** Constructor. Before any operation an internal profile must be built.
     * @see builInternalStructure
     */
    DynCompCol(int n);
    /** Constructor. Before any operation an internal profile must be built.
     * @see builInternalStructure
     */
    DynCompCol();
    /// Copy constructor
    DynCompCol(const DynCompCol &S);
    /// Assignment operator
    DynCompCol &operator=(const DynCompCol &C);
    /// Destructor
    ~DynCompCol();
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
    SparseMtrxType  giveType() const { return SMT_DynCompCol; }
    int isAntisymmetric() const { return 1; }

    /// Returns coefficient at position (i,j). 1-based element access
    virtual double &at(int i, int j);
    /// Returns coefficient at position (i,j). 1-based element access
    virtual double at(int i, int j) const;
    virtual void toFloatMatrix(FloatMatrix &answer) const;
    /// Prints receiver to stdout. Works only for relatively small matrices.
    virtual void printYourself() const;
    virtual void printStatistics() const;

    /*******************************/
    /*  Access and info functions  */
    /*******************************/
#ifndef DynCompCol_USE_STL_SETS
    /// Returns row indx for i-th  column
    const IntArray *row_ind(int i) const { return rowind_ [ i ]; }
    /// Returns column values
    const FloatArray *column(int i) const { return columns_ [ i ]; }
#else
    /// Returns column values
    std :: map< int, double > *column(int i) const { return columns [ i ]; }
#endif

protected:

    /***********************************/
    /*  General access function (slow) */
    /***********************************/
    /// implements 0-based acess
    double operator()(int i, int j) const;
    /// implements 0-based acess
    double &operator()(int i, int j);

#ifndef DynCompCol_USE_STL_SETS
    /// returns the row index of given row at given column, else returns zero.
    int giveRowIndx(int col, int row) const;
    /// insert row entry into column, preserving order of row indexes, returns the index of new row.
    int insertRowInColumn(int col, int row);
#endif
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
};


#endif // dyncompcol_h

