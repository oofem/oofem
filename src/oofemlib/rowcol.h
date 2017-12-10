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

#ifndef rowcol_h
#define rowcol_h

#include "oofemcfg.h"

#include "vector"

namespace oofem {
class IntArray;
class FloatArray;

/**
 * This class implements a segment of a unsymmetric matrix stored in
 * segmented form (skyline).
 * A row-column segment i contains the following items :
 * - the i-th row of the lower half of the matrix ('row')
 * - the i-th column of the upper part of the matrix ('column')
 * - the i-th diagonal coefficient of the matrix ('diag').
 * Since the profile of the matrix is supposed to be symmetric (but not
 * its coefficients), the row and the column have the same length.
 * The row stores coefficients from left to right, the column from up to
 * down ; both exclude the diagonal coefficient. Both start at the first
 * non-zero coefficient, whose position is 'start'.
 *
 * Tasks:
 * - storing and returning coefficients. For segment k:
 *   - any coefficient A(i,k) of the upper part of the matrix (@f$ i < k @f$) is accessed through method 'atU(i)'. It belongs to the column of the segment ;
 *   - any coefficient A(k,i) of the lower part of the matrix (@f$ i < k @f$) is accessed through method 'atL(i)'. It belongs to the row of the segment ;
 *   - the coefficient A(k,k) is obtained through method 'atDiag()' ;
 * - enlarging itself in order to accommodate more coefficients (method 'growTo')
 * - resetting to zero all of its coefficients (method 'reinitialized').
 */
class OOFEM_NO_EXPORT RowColumn
{
protected:
    int number;
    int start;
    std::vector<double> row;
    std::vector<double> column;
    double diag;

public:
    RowColumn(int n);
    RowColumn(int n, int start);

#ifdef DEBUG
    double &atU(int);
    double &atL(int);
    double atU(int) const;
    double atL(int) const;
#else
    double &atU(int i) { return column [ i - start ]; }
    double &atL(int i) { return row [ i - start ]; }
    double atU(int i) const { return column [ i - start ]; }
    double atL(int i) const { return row [ i - start ]; }
#endif
    double &atDiag() { return diag; }
    double atDiag() const { return diag; }
    void checkBounds(int)  const;
    void checkSizeTowards(const IntArray &);
    double dot(const FloatArray &, char, int, int)  const;
    int giveStart()  const { return start; }
    void growTo(int);
    void zero();
    void printYourself() const;
    int giveSize() const { return 1 + 2 * ( number - start ); }
};
} // end namespace oofem
#endif // rowcol_h
