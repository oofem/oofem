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

/*
 * The original idea for this class comes from
 * Dubois-Pelerin, Y.: "Object-Oriented  Finite Elements: Programming concepts and Implementation",
 * PhD Thesis, EPFL, Lausanne, 1992.
 */

#ifndef flotmtrx_h
#define flotmtrx_h

#include "oofemcfg.h"
#include "contextioresulttype.h"
#include "contextmode.h"

#include <vector>
#include <iosfwd>
#include <initializer_list>
#include <algorithm>

#ifdef BOOST_PYTHON
namespace boost {
namespace python {
namespace api {
class object;
};
};
};
#endif

namespace oofem {
class FloatArray;
class IntArray;
class DataStream;

/**
 * Implementation of matrix containing floating point numbers. FloatMatrix can grow and shrink
 * to requested size if required. Implementation includes many computing and manipulation methods.
 * Rows and Columns are indexed starting from 1.
 *
 * The matrix stores its nRows*nColumns coefficients in 'values'. The coefficients
 * are stored column by column, like in Fortran.
 *
 * Tasks:
 * - Storing and retrieving a coefficient (method 'at') ;
 * - Performing standard operations : addition, multiplication, transposition,
 *   inversion, lumping, rotation, etc ;
 * - Introduced allocatedSize variable to allow dynamic rescaling of matrix
 *   size possibly without memory reallocation. At startup matrix occupies space
 *   given by allocatedSpace = nRows*nColums. Then there can be
 *   further request for resizing matrix to smaller dimension
 *   then we only change nRows and nColumns variables , but allocatedSize
 *   variable remain untouched - expecting possible matrix grow and then re-using
 *   previously allocated space.
 *   If further request for growing then is necessary memory reallocation.
 *   This process is controlled in resize member function.
 * 
 * @author Mikael Öhman
 * @author Jim Brouzoulis 
 * @author many others (please add yourselves) 
 */
class OOFEM_EXPORT FloatMatrix
{
protected:
    /// Number of rows.
    int nRows;
    /// Number of columns.
    int nColumns;
    /// Values of matrix stored column wise.
    std :: vector< double >values;

public:
    /**
     * Creates matrix of given size.
     * @param n Number of rows.
     * @param m Requested number of columns.
     */
    FloatMatrix(int n, int m) : nRows(n), nColumns(m), values(n * m) {}
    /// Creates zero sized matrix.
    FloatMatrix() : nRows(0), nColumns(0), values() {}
    /**
     * Constructor. Creates float matrix from float vector. Vector may be stored row wise
     * or column wise, depending on second parameter.
     * @param vector Float vector from which matrix is constructed
     * @param transpose If false (default) then a matrix of size (vector->giveSize(),1)
     * will be created and initialized, if true then a matrix of size (1,vector->giveSize())
     * will be created.
     */
    FloatMatrix(const FloatArray &vector, bool transpose = false);
    /// Copy constructor.
    FloatMatrix(const FloatMatrix &mat) : nRows(mat.nRows), nColumns(mat.nColumns), values(mat.values) {}
    /// Copy constructor.
    FloatMatrix(FloatMatrix && mat) : nRows(mat.nRows), nColumns(mat.nColumns), values( std :: move(mat.values) ) {}
    /// Initializer list constructor.
    FloatMatrix(std :: initializer_list< std :: initializer_list< double > >mat);
    /// Assignment operator.
    FloatMatrix &operator=(std :: initializer_list< std :: initializer_list< double > >mat);
    /// Assignment operator.
    FloatMatrix &operator=(std :: initializer_list< FloatArray >mat);
    /// Assignment operator, adjusts size of the receiver if necessary.
    FloatMatrix &operator=(const FloatMatrix &mat) {
        nRows = mat.nRows;
        nColumns = mat.nColumns;
        values = mat.values;
        return * this;
    }
    FloatMatrix &operator=(FloatMatrix && mat) {
        nRows = std :: move(mat.nRows);
        nColumns = std :: move(mat.nColumns);
        values = std :: move(mat.values);
        return * this;
    }
    /// Destructor.
    ~FloatMatrix() {}

    /**
     * Checks size of receiver towards requested bounds.
     * Current implementation will call exit(1), if positions are outside bounds.
     * @param i Required number of rows.
     * @param j Required number of columns.
     */
    void checkBounds(int i, int j) const;
    /// Returns number of rows of receiver.
    inline int giveNumberOfRows() const { return nRows; }
    /// Returns number of columns of receiver.
    inline int giveNumberOfColumns() const { return nColumns; }
    /// Returns nonzero if receiver is square matrix.
    inline bool isSquare() const { return nRows == nColumns; }
    /// Tests for empty matrix.
    inline bool isNotEmpty() const { return nRows > 0 && nColumns > 0; }

    /// Returns true if no element is NAN or infinite
    bool isFinite() const;

    /**
     * Coefficient access function. Returns value of coefficient at given
     * position of the receiver. Implements 1-based indexing.
     * @param i Row position of coefficient.
     * @param j Column position of coefficient.
     */
#ifdef DEBUG
    double at(int i, int j) const;
#else
    inline double at(int i, int j) const { return values [ ( j - 1 ) * nRows + i - 1 ]; }
#endif
    /**
     * Coefficient access function. Returns value of coefficient at given
     * position of the receiver. Implements 1-based indexing.
     * @param i Row position of coefficient.
     * @param j Column position of coefficient.
     */
#ifdef DEBUG
    double &at(int i, int j);
#else
    inline double &at(int i, int j) { return values [ ( j - 1 ) * nRows + i - 1 ]; }
#endif

    /**
     * Coefficient access function. Returns l-value of coefficient at given
     * position of the receiver. Implements 0-based indexing.
     * @param i Row position of coefficient.
     * @param j Column position of coefficient.
     */
#ifdef DEBUG
    double &operator()(int i, int j);
#else
    inline double &operator()(int i, int j) { return values [ j * nRows + i ]; }
#endif
    /**
     * Coefficient access function. Implements 0-based indexing.
     * @param i Row position of coefficient.
     * @param j Column position of coefficient.
     */
#ifdef DEBUG
    double operator()(int i, int j) const;
#else
    inline double operator()(int i, int j) const { return values [ j * nRows + i ]; }
#endif
    /**
     * Assembles the contribution using localization array into receiver. The receiver must
     * have dimensions large enough to localize contribution.
     * @param src Source to be assembled.
     * @param loc Localization indices.
     */
    void assemble(const FloatMatrix &src, const IntArray &loc);
    /**
     * Assembles the contribution using localization array into receiver. The receiver must
     * have dimensions large enough to localize contribution.
     * @param src Source to be assembled.
     * @param rowind Row localization indices.
     * @param colind Column localization indices.
     */
    void assemble(const FloatMatrix &src, const IntArray &rowind, const IntArray &colind);
    /**
     * Assembles the contribution using localization array into receiver. The receiver must
     * have dimensions large enough to localize contribution.
     * @param src Source to be assembled.
     * @param rowind Row localization indices.
     * @param colind Column localization indices.
     */
    void assemble(const FloatMatrix &src, const int *rowind, const int *colind);

    /**
     * Computes the Frobenius norm of the receiver.
     * The Frobenius norm is defined as the square root of the sum of the absolute squares of its elements.
     * @return Frobenius norm.
     */
    double computeFrobeniusNorm() const;
    /**
     * Computes the operator norm of the receiver.
     * @param p Norm type, '1' for 1 norm, '2' for 2 norm.
     * @return Norm of receiver.
     */
    double computeNorm(char p) const;
    /**
     * Computes the conditioning of the receiver. From 0 to 1, where 0 is singular and 1 best.
     * The receiver must be square.
     * Works identically as MATLAB/Octaves rcond().
     * @param p Norm type, '1' for 1 norm, '2' for 2 norm.
     * @return Conditioning of receiver.
     */
    double computeReciprocalCondition(char p = '1') const;
    /*
     * Computes the eigenvalues of a symmetric matrix.
     * The receiver must be square and symmetric.
     * @param lambda Eigenvalues.
     * @param v Eigenvectors (stored column wise).
     * @param neigs If set, only neigs largest eigenvalues are computed.
     * @return True if successful.
     */
    //bool computeEigenValuesSymmetric(FloatArray &lambda, FloatMatrix &v, int neigs = 0) const;
    /**
     * Modifies receiver to be a diagonal matrix with the components specified in diag.
     * @return Determinant of receiver.
     */
    void beDiagonal(const FloatArray &diag);
    /**
     * Returns the trace of the receiver. The receiver should be a square matrix.
     * @return Trace of receiver.
     */
    double giveTrace() const;
    /**
     * Returns the trace (sum of diagonal components) of the receiver. The receiver should be a square matrix.
     * @return Trace of receiver.
     */
    double giveDeterminant() const;

    /// Zeroes all coefficient of receiver.
    void zero();
    /// Sets receiver to unity matrix.
    void beUnitMatrix();
    /// Sets receiver to the inverse of scaling matrix P multiplied by the deviatoric projector ID.
    void bePinvID();
    /**
     * Assigns to the receiver the transposition of parameter.
     * Grows or shrinks if necessary.
     */
    void beTranspositionOf(const FloatMatrix &src);
    /**
     * Assigns to the receiver product of @f$ a \cdot b @f$ .
     * Grows or shrinks if necessary.
     */
    void beProductOf(const FloatMatrix &a, const FloatMatrix &b);
    /**
     * Adds to the receiver product of @f$ a \cdot b @f$ .
     */
    void addProductOf(const FloatMatrix &a, const FloatMatrix &b);
    /**
     * Adds to the receiver product of @f$ a^{\mathrm{T}} \cdot b @f$ .
     */
    void addTProductOf(const FloatMatrix &a, const FloatMatrix &b);
    /**
     * Assigns to the receiver product of @f$ a^{\mathrm{T}} \cdot b @f$ .
     * Grows or shrinks if necessary.
     */
    void beTProductOf(const FloatMatrix &a, const FloatMatrix &b);
    /**
     * Assigns to the receiver product of @f$ a \cdot b^{\mathrm{T}}@f$ .
     * Grows or shrinks if necessary.
     */
    void beProductTOf(const FloatMatrix &a, const FloatMatrix &b);
    /**
     * Assigns to the receiver the dyadic product @f$ v_1 \cdot v_2^{\mathrm{T}} @f$ .
     * Grows or shrinks if necessary.
     */
    void beDyadicProductOf(const FloatArray &vec1, const FloatArray &vec2);
    /**
     * Assigns the receiver to be a repeated diagonal matrix.
     * @param n Vector with components which will appear in respective diagonal.
     * @param nsd Number of spatial dimensions
     */
    void beNMatrixOf(const FloatArray &n, int nsd);
    /**
     * Makes receiver the local coordinate for the given normal.
     * Implemented for 2D and 3D.
     * @param normal Normal (normalized).
     */
    void beLocalCoordSys(const FloatArray &normal);
    /**
     * Adds the given matrix as sub-matrix to receiver. The sub-matrix values will be added to receivers
     * corresponding receiver's values at positions (ri...ri+src.nrows, ci....ci+src.ncolumns).
     * The size of receiver will be adjusted, if necessary.
     * @param src Sub-matrix to be inserted.
     * @param sr Determines the row position (receiver's 1-based index) of first src value to be added.
     * @param sc Determines the column position (receiver's 1-based index) of first src value to be added.
     */
    void setSubMatrix(const FloatMatrix &src, int sr, int sc);
    /**
     * Adds the given matrix as sub-matrix to receiver.
     * The same as setSubMatrix but transposes the source matrix.
     * @param src Sub-matrix to be inserted transposed.
     * @param sr Determines the row position (receiver's 1-based index) of first src value to be added.
     * @param sc Determines the column position (receiver's 1-based index) of first src value to be added.
     */
    void setTSubMatrix(const FloatMatrix &src, int sr, int sc);
    /**
     * Assigns to the receiver the sub-matrix of another matrix.
     * @param src Matrix from which sub-matrix is taken
     * @param topRow Index of top row, where sub matrix row index starts
     * @param bottomRow Index of bottom row, where sub matrix ends (including this row).
     * @param topCol Index of top column of sub-matrix.
     * @param bottomCol index of bottom column of sub-matrix.
     */
    void beSubMatrixOf(const FloatMatrix &src, int topRow, int bottomRow, int topCol, int bottomCol);
    /**
     * Modifies receiver to be a sub-matrix of another matrix.
     * @param src Matrix from which sub-matrix is taken
     * @param indxRow Array of the row indices to be extracted
     * @param indxCol Array of the column indices to be extracted
     */
    void beSubMatrixOf(const FloatMatrix &src, const IntArray &indxRow, const IntArray &indxCol);
    /**
     * Adds given vector to receiver starting at given position.
     * @param src Source matrix.
     * @param sr Starting row position.
     * @param sc Starting column position.
     */
    void addSubVectorRow(const FloatArray &src, int sr, int sc);
    /**
     * Adds given vector to receiver starting at given position.
     * @param src Source matrix.
     * @param sr Starting row position.
     * @param sc Starting column position.
     */
    void addSubVectorCol(const FloatArray &src, int sr, int sc);

    /**
     * Copy (set) given vector to receiver row sr, starting at column sc.
     * @param src Source matrix.
     * @param sr Starting row position.
     * @param sc Starting column position.
     */
    void copySubVectorRow(const FloatArray &src, int sr, int sc);
    /**
     * Sets the values of the matrix in specified column. If matrix size is zero, the size is adjusted.
     * @param src Array to set at column c.
     * @param c Column to set.
     */
    void setColumn(const FloatArray &src, int c);
    /**
     * Fetches the values from the specified column. Output array is resized to fit column.
     * @param dest Array to copy values to.
     * @param c Column to copy.
     */
    void copyColumn(FloatArray &dest, int c) const;
    /**
     * Modifies receiver to become inverse of given parameter. Size of receiver will be adjusted.
     * @param src Matrix to be inverted.
     */
    void beInverseOf(const FloatMatrix &src);
    /**
     * Solves the  system of linear equations @f$ K\cdot a = b @f$ . Uses Gaussian elimination with pivoting directly on receiver.
     * @param b RHS of linear system.
     * @param answer Solution of linear equations.
     * @param transpose Solves for the transpose of K.
     * @return False if K is singular, otherwise true.
     */
    bool solveForRhs(const FloatArray &b, FloatArray &answer, bool transpose = false);
    /**
     * Solves the  system of linear equations @f$ K\cdot A = B @f$ . Uses Gaussian elimination with pivoting directly on receiver.
     * @param B RHS of linear system.
     * @param answer Solution of linear equations, each column corresponding to columns in B.
     * @param transpose Solves for the transpose of K.
     */
    void solveForRhs(const FloatMatrix &B, FloatMatrix &answer, bool transpose = false);
    /**
     * Adds to the receiver the product @f$ a^{\mathrm{T}}\cdot b \mathrm{d}V @f$. If the receiver has zero size, it is expanded.
     * Assumes that receiver and product @f$ a^{\mathrm{T}}\cdot b \mathrm{d}V @f$ are symmetric matrices. Computes only the
     * upper half of receiver.
     * @param a Matrix a in equation.
     * @param b Matrix b in equation.
     * @param dV Scaling factor.
     */
    void plusProductSymmUpper(const FloatMatrix &a, const FloatMatrix &b, double dV);
    /**
     * Adds to the receiver the dyadic product @f$ a \otimes a \mathrm{d}V @f$. If the receiver has zero size, it is expanded.
     * Computes only the upper half of receiver.
     * @param a Array a in equation.
     * @param dV Scaling factor.
     */
    void plusDyadSymmUpper(const FloatArray &a, double dV);
    /**
     * Adds to the receiver the product @f$a^{\mathrm{T}} \cdot b \mathrm{d}V@f$. If the receiver has zero size, it is expanded.
     * @param a Matrix a in equation.
     * @param b Matrix b in equation.
     * @param dV Scaling factor.
     */
    void plusProductUnsym(const FloatMatrix &a, const FloatMatrix &b, double dV);
    /**
     * Adds to the receiver the product @f$a \otimes b \mathrm{d}V@f$. If the receiver has zero size, it is expanded.
     * @param a Array a in equation.
     * @param b Array b in equation.
     * @param dV Scaling factor.
     */
    void plusDyadUnsym(const FloatArray &a, const FloatArray &b, double dV);

    /**
     * Adds matrix to the receiver.
     * If receiver has zero size, size is accordingly adjusted.
     * @param a Matrix to be added.
     */
    void add(const FloatMatrix &a);
    /**
     * Adds matrix to the receiver scaled by s.
     * If receiver has zero size, size is accordingly adjusted.
     * @param s Scaling factor.
     * @param a Matrix to be added.
     */
    void add(double s, const FloatMatrix &a);
    /**
     * Subtracts matrix from the receiver.
     * @param a Matrix to be subtracted.
     */
    void subtract(const FloatMatrix &a);
    /**
     * Multiplies receiver by factor f.
     * @param f Factor to multiply by.
     */
    void times(double f);
    /**
     * Changes sign of receiver values.
     */
    void negated();
    /**
     * Assigns to receiver one column or one row matrix containing vector.
     * @param vector Source vector.
     * @param transposed If false then (vector->giveSize(),1) FloatMatrix is assigned
     * else (1,vector->giveSize()) FloatMatrix is assigned.
     */
    void initFromVector(const FloatArray &vector, bool transposed);
    /**
     * Initializes the lower half of the receiver according to the upper half.
     */
    void symmetrized();
    /**
     * Returns the receiver 'a' transformed using give transformation matrix r.
     * The method performs the operation  @f$ a = r^{\mathrm{T}} \cdot a \cdot r@f$ .
     * @param r Transformation matrix.
     * @param mode If set to 't' then the transpose of the rotation matrix is used instead.
     */
    void rotatedWith(const FloatMatrix &r, char mode = 'n');
    /**
     * Checks size of receiver towards requested bounds.
     * If dimension mismatch, size is adjusted accordingly.
     * Content is always zeroed.
     * @param rows New number of rows.
     * @param cols New number of columns.
     */
    void resize(int rows, int cols);
    /**
     * Checks size of receiver towards requested bounds.
     * If dimension mismatch, size is adjusted accordingly.
     * Note: New coefficients are initialized to zero, old are kept.
     */
    void resizeWithData(int, int);
    /**
     * Resizing that enforces reallocation of memory.
     * Data is zeroed.
     * @param r Number of rows.
     * @param c Number of columns.
     */
    void hardResize(int r, int c);
    /// Sets size of receiver to be an empty matrix. It will have zero rows and zero columns size.
    void clear() {
        this->nRows = 0;
        this->nColumns = 0;
    }
    /**
     * Computes eigenvalues and eigenvectors of receiver (must be symmetric)
     * The receiver is preserved.
     * @param eval Requested eigenvalues.
     * @param v Requested eigenvectors (stored colum wise).
     * @param nf Number of significant figures.
     * @return True if ok,otherwise false.
     */
    bool jaco_(FloatArray &eval, FloatMatrix &v, int nf);

    /// Prints matrix to stdout. Useful for debugging.
    void printYourself() const;
    /**
     * Print receiver on stdout with custom name.
     * @param name Display name of reciever.
     */
    void printYourself(const std::string &name) const;

    /// Higher accuracy than printYourself.
    void pY() const;

    /**
     * Writes receiver as CSV (comma seperated values)
     * @param name Filename
     */
    void writeCSV(const std :: string &name) const;

    /**
     * Exposes the internal values of the matrix. Should typically not be used outside of matrix classes.
     * @return Pointer to the values of the matrix.
     */
    inline const double *givePointer() const { return values.data(); }
    inline double *givePointer() { return values.data(); }

    /**
     * Reciever will be a 3x3 matrix formed from a vector with either 9 or 6 components.
     * Order of matrix components in vector: 11, 22, 33, 23, 13, 12, 32, 31, 21
     * If size(aArray) = 6, a symmetric matrix will be created.
     * @param aArray Array to transform.
     */
    void beMatrixFormOfStress(const FloatArray &aArray);
    void beMatrixForm(const FloatArray &aArray);

    /**
     * Swaps the indices in the 6x6 matrix such that it converts between OOFEM's
     * and Abaqus' way of writing matrices. Currently used to convert the 6x6 Jacobian
     * from Abaqus UMAT to OOFEM.
     */
    void changeComponentOrder();

    // Overloaded methods:
    contextIOResultType storeYourself(DataStream &stream) const;
    contextIOResultType restoreYourself(DataStream &stream);
    int givePackSize(DataStream &buff) const;

    friend std :: ostream &operator<<(std :: ostream &out, const FloatMatrix &r);


#ifdef BOOST_PYTHON
    void __setitem__(boost :: python :: api :: object t, double val);
    double __getitem__(boost :: python :: api :: object t);
    void beCopyOf(const FloatMatrix &src) { this->operator=(src); }
#endif
};
} // end namespace oofem
#endif // flotmtrx_h
