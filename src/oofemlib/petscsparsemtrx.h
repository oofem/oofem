/* $Header: /home/cvs/bp/oofem/oofemlib/src/Attic/petscsparsemtrx.h,v 1.1.2.1 2004/04/05 15:19:43 bp Exp $ */
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
#ifndef petscsparsemtrx_h
#define petscsparsemtrx_h

#ifdef __PETSC_MODULE

 #include "sparsemtrx.h"

 #ifndef __MAKEDEPEND
  #include "petscksp.h"
 #endif

namespace oofem {
/**
 * This class provides an sparse matrix interface to PETSc Matrices
 */
class PetscSparseMtrx : public SparseMtrx
{
protected:
    Mat mtrx;
    bool symmFlag;
    MatType mType;
    int leqs;
    int geqs;
    int di;
    EquationID ut;

    EngngModel *emodel;

public:
 #ifdef __PARALLEL_MODE
    PetscSparseMtrx(int n, int m) : SparseMtrx(n, m)
    {
        di = 0;
        leqs = n;
        geqs = n;
        symmFlag = false;
        mtrx = NULL;
    }
 #else
    PetscSparseMtrx(int n, int m) : SparseMtrx(n, m)
    {
        di = 0;
        leqs = n;
        geqs = n;
        symmFlag = false;
        mtrx = NULL;
    }
 #endif
    PetscSparseMtrx() : SparseMtrx() {
        di = 0;
        mtrx = NULL;
    }
    ~PetscSparseMtrx() { MatDestroy(mtrx); }

    /** Returns {\bf newly allocated} copy of receiver. Programmer must take
     * care about proper deallocation of allocated space.
     * @return newly allocated copy of receiver */
    virtual SparseMtrx *GiveCopy() const;

    /**
     * Evaluates @f$ y = A\cdot x @f$
     * @param x Array to be multiplied with receiver.
     * @param answer y.
     */
    virtual void times(const FloatArray &x, FloatArray &answer) const;
    /**
     * Evaluates @f$ y = A^{\rm T}\cdot x @f$
     * @param x Array to be multiplied with transpose of the receiver.
     * @param answer y.
     */
    virtual void timesT(const FloatArray &x, FloatArray &answer) const;
    /**
     * Evaluates @f$ C = A^{\rm T}\cdot B @f$
     * @param B array to be multiplied with receiver.
     * @param answer C.
     */
    virtual void times(const FloatMatrix &B, FloatMatrix &answer) const;
    /**
     * Evaluates @f$ C = A^{\rm T}\cdot B @f$
     * @param x matrix to be multiplied with receiver.
     * @param answer C.
     */
    virtual void timesT(const FloatMatrix &B, FloatMatrix &answer) const;
    /**
     * Multiplies receiver by scalar value.
     * @param x value to multiply receiver
     */
    virtual void times(double x);
    /**
     * Builds internal structure of receiver. This method determines the internal profile
     * of sparse matrix, allocates necessary space for storing nonzero coefficients and
     * initializes receiver. In general, the profile of sparse matrix is determined
     * using one (or more) loop over local code numbers of elements.
     * This method must be called before any operation, like assembly, zeroing,
     * or multiplication.
     * @param eModel pointer to corresponding engineering model
     * @param di domain index specify which domain to use
     */
    virtual int buildInternalStructure(EngngModel * eModel, int di, EquationID, const UnknownNumberingScheme & s);
    /**
     * Builds internal structure.
     * Same as buildInternalStructure except with different numbering schemes for rows and columns
     * @param eModel pointer to corresponding engineering model.
     * @param di domain index to specify which domain to use.
     * @param eid equation id to use to numbering scheme.
     * @param r_s numbering scheme for rows.
     * @param c_s numbering scheme for columns.
     */
    virtual int buildInternalStructure(EngngModel *eModel, int di, EquationID eid, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);
    // virtual int assemble (FloatMatrix*, IntArray*) = 0;
    /**
     * Assembles sparse matrix from contribution of local elements. This method for
     * each element adds its contribution to itself. Mapping between local element
     * contribution and its global position is given by local code numbers of element.
     * @param loc location array. The values corresponding to zero loc array value are not assembled.
     * @param mat contribution to be assembled using loc array.
     */
    virtual int assemble(const IntArray &loc, const FloatMatrix &mat);
    /**
     * Assembles sparse matrix from contribution of local elements. This method for
     * each element adds its contribution to itself. Mapping between local element
     * contribution and its global position is given by row and column local code numbers.
     * @param rloc row location array. The values corresponding to zero loc array value are not assembled.
     * @param cloc column location array. The values corresponding to zero loc array value are not assembled.
     * @param mat contribution to be assembled using rloc and cloc arrays. The rloc position determines the row, the
     * cloc determines the corresponding column.
     */
    virtual int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat);

    int assembleBegin();
    int assembleEnd();

    /// Determines, whether receiver can be factorized.
    virtual int canBeFactorized() const { return 0; }
    /**
     * Returns the receiver factorized. \f$L^T D L\f$ form is used.
     * @return pointer to the receiver
     */
    virtual SparseMtrx *factorized() { return NULL; }
    /**
     * Computes the solution of linear system \f$A x = y\f$. A is receiver.
     * solution vector x overwrites the right hand side vector y.
     * Receiver must be in factorized form.
     * @param y right hand side on input, solution on output.
     * @return pointer to y array
     * @see factorized method
     */
    virtual FloatArray *backSubstitutionWith(FloatArray &y) const { return NULL; }
    /// Zeroes the receiver.
    virtual SparseMtrx *zero();

    /// Returns coefficient at position (i,j).
    virtual double &at(int i, int j);
    /// Returns coefficient at position (i,j).
    virtual double at(int i, int j) const;
    virtual void toFloatMatrix(FloatMatrix &answer) const;
    /// Prints the receiver statistics (one-line) to stdout.
    virtual void printStatistics() const;
    /// Prints receiver to stdout. Works only for relatively small matrices.
    virtual void printYourself() const;

    /// Sparse matrix type identification
    virtual SparseMtrxType  giveType() const { return SMT_PetscMtrx; }
    /// Returns nonzero if anti-symmetric
    virtual int isAntisymmetric() const { if ( symmFlag ) { return 0; } else { return 1; } }

    Mat *giveMtrx() { return & this->mtrx; }
    int       giveSymmetryFlag() const { return symmFlag; }
    int       setOption(MatOption op, PetscTruth flag) { return MatSetOption(this->mtrx, op, flag); }
    EquationID       giveEquationID() { return ut; }

 #ifdef IML_COMPAT
    // /***********************************/
    //  /*  Matrix/Vector multiply         */
    //  /***********************************/

    virtual FloatArray operator*(const FloatArray &x) const
    {
        FloatArray answer;
        this->times(x, answer);
        return answer;
    }
    virtual FloatArray trans_mult(const FloatArray &x) const
    {
        FloatArray answer;
        this->timesT(x, answer);
        return answer;
    }

 #endif

    int giveLeqs() { return leqs; }
    int giveDomainIndex() { return di; }
};
} // end namespace oofem
#endif
#endif // petscsparsemtrx_h
