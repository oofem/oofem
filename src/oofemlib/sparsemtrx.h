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

#ifndef sparsemtrx_h
#define sparsemtrx_h

#include "oofemcfg.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "error.h"
#include "sparsemtrxtype.h"

namespace oofem {
class EngngModel;
class TimeStep;
class UnknownNumberingScheme;

/**
 * Base class for all matrices stored in sparse format. Basically sparse matrix
 * contains contribution of local element matrices. Localization of local element
 * matrix into global (structural) matrix is determined using element code numbers.
 * Basic methods include
 * - Building internal structure of sparse matrix (according to code numbers of elements)
 * - Assembling of local element matrices
 * - Multiplication by array
 * - Possible factorization and back substitution.
 */
class OOFEM_EXPORT SparseMtrx
{
public:
    typedef long SparseMtrxVersionType;

protected:
    /// Number of rows.
    int nRows;
    /// Number of columns.
    int nColumns;

    /**
     * Allows to track if receiver changes.
     * Any change on receiver should increment this variable.
     * The factorization should not change version (nrsolver direct control relies on that).
     * This version info is used for example by preconditioners associated to
     * particular matrix; the preconditioner initialization can be demanding
     * and this versioning allows to reuse initialized preconditioner for same
     * matrix, if there is no change;
     */
    SparseMtrxVersionType version;

public:
    /**
     * Constructor, creates (n,m) sparse matrix. Due to sparsity character of matrix,
     * not all coefficient are physically stored (in general, zero members are omitted).
     */
    SparseMtrx(int n, int m) : nRows(n), nColumns(m), version(0) { }
    /// Constructor
    SparseMtrx() : nRows(0), nColumns(0), version(0) { }
    /// Destructor
    virtual ~SparseMtrx() { }

    /// Return receiver version.
    SparseMtrxVersionType giveVersion() { return this->version; }

    /**
     * Checks size of receiver towards requested bounds.
     * Current implementation will call exit(1), if positions are outside bounds.
     * @param i Required number of rows.
     * @param j Required number of columns.
     */
    void checkBounds(int i, int j) const {
        if ( i <= 0 ) {
            OOFEM_ERROR("matrix error on rows : %d <= 0", i);
        } else if ( j <= 0 ) {
            OOFEM_ERROR("matrix error on columns : %d <= 0", j);
        } else if ( i > nRows ) {
            OOFEM_ERROR("matrix error on rows : %d => %d", i, nRows);
        } else if ( j > nColumns ) {
            OOFEM_ERROR("matrix error on columns : %d => %d", j, nColumns);
        }
    }
    /// Returns number of rows of receiver.
    int giveNumberOfRows() const { return nRows; }
    /// Returns number of columns of receiver.
    int giveNumberOfColumns() const { return nColumns; }
    /// Returns nonzero if receiver is square matrix.
    bool isSquare() const { return nRows == nColumns; }
    /// Tests for empty matrix.
    bool isNotEmpty() const { return nRows > 0 && nColumns > 0; }

    /**
     * Returns a <em>newly allocated</em> copy of receiver. Programmer must take
     * care about proper deallocation of allocated space.
     * @return Newly allocated copy of receiver.
     */
    virtual SparseMtrx *GiveCopy() const {
        OOFEM_ERROR("Not implemented");
        return NULL;
    }

    /**
     * Evaluates @f$ y = A \cdot x @f$
     * @param x Array to be multiplied with receiver.
     * @param answer y.
     */
    virtual void times(const FloatArray &x, FloatArray &answer) const { OOFEM_ERROR("Not implemented"); }
    /**
     * Evaluates @f$ y = A^{\mathrm{T}} \cdot x @f$
     * @param x Array to be multiplied with transpose of the receiver.
     * @param answer y.
     */
    virtual void timesT(const FloatArray &x, FloatArray &answer) const { OOFEM_ERROR("Not implemented"); }
    /**
     * Evaluates @f$ C = A^{\mathrm{T}} \cdot B @f$
     * @param B Array to be multiplied with receiver.
     * @param answer C.
     */
    virtual void times(const FloatMatrix &B, FloatMatrix &answer) const { OOFEM_ERROR("Not implemented"); }
    /**
     * Evaluates @f$ C = A^{\mathrm{T}} \cdot B @f$
     * @param B Matrix to be multiplied with receiver.
     * @param answer C.
     */
    virtual void timesT(const FloatMatrix &B, FloatMatrix &answer) const { OOFEM_ERROR("Not implemented"); }
    /**
     * Multiplies receiver by scalar value.
     * @param x Value to multiply receiver.
     */
    virtual void times(double x) { OOFEM_ERROR("Not implemented"); }
    /**
     * Adds x * m.
     * @param x Value to multiply m by.
     * @param m Matrix to add (should of the same matrix type).
     */
    virtual void add(double x, SparseMtrx &m) { OOFEM_ERROR("Not implemented"); }
    /**
     * Adds x * m (treats m as a diagonal matrix, stored as an array)
     * @param x Value to multiply m by.
     * @param m Matrix to add.
     */
    virtual void addDiagonal(double x, FloatArray &m) {
        for ( int i = 1; i <= m.giveSize(); ++i ) {
            this->at(i, i) += x * m.at(i);
        }
    }

    /**
     * Builds internal structure of receiver. This method determines the internal profile
     * of sparse matrix, allocates necessary space for storing nonzero coefficients and
     * initializes receiver. In general, the profile of sparse matrix is determined
     * using one (or more) loop over local code numbers of elements.
     * This method must be called before any operation, like assembly, zeroing,
     * or multiplication.
     * @param eModel Pointer to corresponding engineering model.
     * @param di Domain index specify which domain to use.
     * @param s Determines unknown numbering scheme.
     * @return Zero iff successful.
     */
    virtual int buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme &s) = 0;
    /**
     * Build internal structure of receiver.
     * @see buildInternalStructure
     * @param eModel Pointer to corresponding engineering model.
     * @param di Domain index specify which domain to use.
     * @param r_s Determines unknown numbering scheme for the rows.
     * @param c_s Determines unknown numbering scheme for the columns.
     * @return Zero iff successful.
     */
    virtual int buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme &r_s,
                                       const UnknownNumberingScheme &c_s) {
        OOFEM_ERROR("Not implemented");
        return 0;
    }
    /**
     * Assembles sparse matrix from contribution of local elements. This method for
     * each element adds its contribution to itself. Mapping between local element
     * contribution and its global position is given by local code numbers of element.
     * @param loc Location array. The values corresponding to zero loc array value are not assembled.
     * @param mat Contribution to be assembled using loc array.
     * @return Zero iff successful.
     */
    virtual int assemble(const IntArray &loc, const FloatMatrix &mat) = 0;
    /**
     * Assembles sparse matrix from contribution of local elements. This method for
     * each element adds its contribution to itself. Mapping between local element
     * contribution and its global position is given by row and column local code numbers.
     * @param rloc Row location array. The values corresponding to zero loc array value are not assembled.
     * @param cloc Column location array. The values corresponding to zero loc array value are not assembled.
     * @param mat Contribution to be assembled using rloc and cloc arrays. The rloc position determines the row, the
     * cloc position determines the corresponding column.
     * @return Zero iff successful.
     */
    virtual int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat) = 0;

    /// Starts assembling the elements.
    virtual int assembleBegin() { return 1; }
    /// Returns when assemble is completed.
    virtual int assembleEnd() { return 1; }

    /// Determines, whether receiver can be factorized.
    virtual bool canBeFactorized() const = 0;
    /**
     * Returns the receiver factorized. @f$ L^{\mathrm{T}} \cdot D \cdot L @f$ form is used.
     * @return pointer to the receiver
     */
    virtual SparseMtrx *factorized() { return NULL; }
    /**
     * Computes the solution of linear system @f$ A\cdot x = y @f$ where A is receiver.
     * Solution vector x overwrites the right hand side vector y.
     * Receiver must be in factorized form.
     * @param y Right hand side on input, solution on output.
     * @return Pointer to y array.
     */
    virtual FloatArray *backSubstitutionWith(FloatArray &y) const { return NULL; }
    /// Zeroes the receiver.
    virtual void zero() = 0;

    /// Returns the norm of receiver.
    virtual double computeNorm() const {
        OOFEM_ERROR("Not implemented");
        return 0.0;
    }

    virtual SparseMtrx *giveSubMatrix(const IntArray &rows, const IntArray &cols) 
    { OOFEM_ERROR("Not implemented"); return NULL; }
    
    /// Returns coefficient at position (i,j).
    virtual double &at(int i, int j) = 0;
    /// Returns coefficient at position (i,j).
    virtual double at(int i, int j) const = 0;
    /// Checks whether memory is allocated at position (i,j).
    virtual bool isAllocatedAt(int i, int j) const { return false; }
    /// Converts receiving sparse matrix to a dense float matrix.
    virtual void toFloatMatrix(FloatMatrix &answer) const { OOFEM_ERROR("Not implemented"); }
    /// Prints the receiver statistics (one-line) to stdout.
    virtual void printStatistics() const {  OOFEM_LOG_INFO("Not implemented file %s, line %d\n", __FILE__, __LINE__); }
    /// Prints receiver to stdout. Works only for relatively small matrices.
    virtual void printYourself() const { OOFEM_LOG_INFO("Not implemented file %s, line %d\n", __FILE__, __LINE__); }
    /// Helpful for debugging, writes the matrix to given file.
    virtual void writeToFile(const char *fname) const { OOFEM_LOG_INFO("Not implemented file %s, line %d\n", __FILE__, __LINE__); }
    /// Sparse matrix type identification
    virtual SparseMtrxType giveType() const = 0;
    /// Returns true if asymmetric
    virtual bool isAsymmetric() const = 0;

    virtual const char *giveClassName() const = 0;
    /// Error printing helper.
    std :: string errorInfo(const char *func) const { return std :: string(giveClassName()) + func; }
    ///@name IML compatibility
    //@{
    /// IML compatibility, @f$ A \cdot x@f$
    FloatArray operator *( const FloatArray & x ) const
    {
        FloatArray answer;
        this->times(x, answer);
        return answer;
    }
    /// IML compatibility, @f$ A^{\mathrm{T}} \cdot x@f$
    FloatArray trans_mult(const FloatArray &x) const
    {
        FloatArray answer;
        this->timesT(x, answer);
        return answer;
    }
    //@}
};
} // end namespace oofem
#endif // sparsemtrx_h
