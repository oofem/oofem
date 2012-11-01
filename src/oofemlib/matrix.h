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

/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

#ifndef matrix_h
#define matrix_h

namespace oofem {

/**
 * Abstract root class for matrices.
 * The matrix is characterized by its size.
 */
class Matrix
{
protected:
    /// Number of rows.
    int nRows;
    /// Number of columns.
    int nColumns;

public:
    /// Constructor, creates zero sized matrix.
    Matrix() { nRows = nColumns = 0; }
    /**
     * Constructor, creates zero matrix, with given size.
     * @param n Number of rows.
     * @param m Number of rows.
     */
    Matrix(int n, int m)  {
        nRows = n;
        nColumns = m;
    }
    /// Destructor.
    virtual ~Matrix() { }

    /**
     * Checks size of receiver towards requested bounds.
     * Current implementation will call exit(1), if positions are outside bounds.
     * @param i Required number of rows.
     * @param j Required number of columns.
     */
    void checkBounds(int i, int j) const;
    /// Returns number of rows of receiver.
    int giveNumberOfRows() const { return nRows; }
    /// Returns number of columns of receiver.
    int giveNumberOfColumns() const { return nColumns; }
    /// Returns nonzero if receiver is square matrix.
    bool isSquare() const { return nRows == nColumns; }
    /// Tests for empty matrix.
    bool isNotEmpty() const { return nRows > 0 && nColumns > 0; }
    /// Prints receiver on stdout. Useful mainly for debugging.
    virtual void printYourself() const;
};
} // end namespace oofem
#endif // matrix_h
