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

#ifndef skyline_h
#define skyline_h

#include "sparsemtrx.h"

namespace oofem {
/**
 * Class implementing sparse matrix stored in skyline form. This class
 * assumes symmetric form of matrix to be stored, i.e., only upper half
 * part is therefore stored. Coefficients are stored in one dimensional
 * float array, containing for each column the values from the diagonal
 * to the last nonzero entry.
 * This requires to remember the addresses of diagonal members in such
 * array (adr attribute)
 *
 * Attribute 'isFactorized' is True if the skyline is already in U(T).D.U
 *         factorized form, else it is False.
 * Attribute adr  is array of diagonal members, it's size is size+1 (adr[0]=neni, adr[1]=1)
 * Attribute mtrx is double pointer to skyline stored in a array form
 *         (but we start from index 1)
 *
 * Tasks:
 * - building its internal storage structure (method 'buildInternalStructure')
 * - store and localize local matrices (method 'localize')
 * - performing standard operations : multiplication by array (method 'times')
 * - possible factorization and backSubstitution (recognized by nonzero result of
 *   canBeFactorized) (methods 'factorize' and 'backSubstitution')
 * - setting all coefficients to zero (method 'zero')
 *
 * @see SparseMtrx class
 */
class OOFEM_EXPORT Skyline : public SparseMtrx
{
protected:
    /// Total number of nonzero coefficients stored.
    int nwk;
    /// Integer array holding addresses of diagonal members.
    IntArray adr;
    /// Vector of stored coefficients.
    double *mtrx;
    /// Flag indicating whether factorized.
    int isFactorized;

public:
    /**
     * Constructor. Before any operation an internal profile must be built.
     * @see builInternalStructure
     */
    Skyline(int n);
    /**
     * Constructor. Before any operation an internal profile must be built.
     * @see builInternalStructure
     */
    Skyline();
    /// Destructor
    virtual ~Skyline();

    virtual SparseMtrx *GiveCopy() const;

    virtual void times(const FloatArray &x, FloatArray &answer) const;
    virtual void timesT(const FloatArray &x, FloatArray &answer) const { this->times(x, answer); }
    virtual void times(double x);
    virtual void add(double x, SparseMtrx &m);
    virtual int buildInternalStructure(EngngModel *, int, const UnknownNumberingScheme &);

    //virtual Skyline *giveSubMatrix(Skyline &mat, IntArray &rows, IntArray &cols);
    //virtual Skyline *beSubMatrixOf(const Skyline &mat, IntArray &rows, IntArray &cols);
    //virtual SparseMtrx *beSubMatrixOf(const SparseMtrx &mat, IntArray &rows, IntArray &cols);
    virtual SparseMtrx *giveSubMatrix(const IntArray &rows, const IntArray &cols); 
    /**
     * Allocates and builds internal structure according to given
     * array holding addresses of diagonal members values (adr).
     */
    int setInternalStructure(IntArray &a);

    virtual int assemble(const IntArray &loc, const FloatMatrix &mat);
    virtual int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat);

    virtual bool canBeFactorized() const { return true; }
    virtual SparseMtrx *factorized();
    virtual FloatArray *backSubstitutionWith(FloatArray &) const;
    virtual void zero();
    /**
     * Splits the receiver to LDLT form,
     * and computes the rigid body motions.
     * @param r matrix containing the rigid body motions base vectors.
     * @param nse number of rigid body motions
     * @param se array containing indexes of singular eqs.
     * @param limit determines linear dependence or independence
     * @param tc type of calculation
     * = 1 = Decomposition to LDL form
     * = 2 = ? konstruovani baze prostoru ? Ker A
     * = 3 = Decomposition to LDL form and ? konstruovani baze prostoru ? Ker A
     */
    void rbmodes(FloatMatrix &r, int &nse, IntArray &se, double limit, int tc);
    /**
     * Solves the singular system of equations, the receiver should be factorized
     * using rbmodes service.
     * @param x solution vector
     * @param y right hand side
     * @param nse number of rigid body motions
     * @param limit determines the linear dependence or independence
     * @param se indexes of singular equations
     */
    void ldl_feti_sky(FloatArray &x, FloatArray &y, int nse, double limit, IntArray &se);
    /// Returns coefficient at position (i,j).
    virtual double &at(int, int);
    /// Returns coefficient at position (i,j).
    virtual double at(int i, int j) const;
    /// Returns 0 if the memory is not allocated at position (i,j).
    virtual bool isAllocatedAt(int i, int j) const;
    int giveNumberOfNonZeros() const { return this->nwk; }
    virtual void toFloatMatrix(FloatMatrix &answer) const;
    /// Prints receiver to stdout.
    virtual void printYourself() const;
    virtual void writeToFile(const char *fname) const;
    int giveAllocatedSize() { return nwk; }

    virtual SparseMtrxType giveType() const { return SMT_Skyline; }
    virtual bool isAsymmetric() const { return false; }

    virtual const char *giveClassName() const { return "Skyline"; }

protected:
    Skyline(int, int, double *, const IntArray &);
};
} // end namespace oofem
#endif // skyline_h
