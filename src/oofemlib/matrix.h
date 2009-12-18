/* $Header: /home/cvs/bp/oofem/oofemlib/src/matrix.h,v 1.5 2003/04/06 14:08:25 bp Exp $ */
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

/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

//   ********************
//   *** CLASS MATRIX ***
//   ********************


#ifndef matrix_h
#define matrix_h

namespace oofem {

/* class FloatMatrix ; class IntArray ; */

/** Abstract root class  for matrices. Matrix is characterized by its size.
 */
class Matrix
{
    /*
     * This abstract class is the superclass of the class that implement
     * matrices (FloatMatrix, PolynomialMatrix,...).
     * DESCRIPTION :
     * A matrix is characterized by its number of rows and columns.
     * TASKS :
     * Its tasks are defined in the subclasses.
     */


protected:
    /// Number of rows
    int nRows;
    /// Number of columns
    int nColumns;

public:
    /// Constructor, creates zero sized matrix
    Matrix()             { nRows = nColumns = 0; }      // constructors
    /** Constructor. Creates zero matrix, with given size.
     * @param n number of rows
     * @param m number of rows
     */
    Matrix(int n, int m)  { nRows = n;
                            nColumns = m; }
    /// Destructor.
    virtual ~Matrix()    { }                            // destructor

    /** checks size of receiver towards requested bounds.
     * Current implementation will call exit(1), if dimension
     * mismatch found.
     * @param i required number of rows
     * @param j required number of clumns
     */
    void          checkBounds(int i, int j) const;
    /// Returns number of rows of receiver
    int           giveNumberOfRows()    const { return nRows; }
    /// Returns number of columns of receiver
    int           giveNumberOfColumns() const { return nColumns; }
    /// Returns nonzero if receiver is square matrix.
    int           isSquare()            const { return ( !nRows - nColumns ); }
    /// Tests for empty matrix.
    int           isNotEmpty()             const { return ( ( ( nRows == 0 ) && ( nColumns == 0 ) ) ? 0 : 1 ); }
    /// Prints receiver on stdin. Usefull mainly for debugging
    virtual void  printYourself() const { }
};

} // end namespace oofem
#endif // matrix_h
