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

#define _IFT_Skyline_Name "skyline"

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
 * Attribute mtrx is matrix values to skyline stored in a array form
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
    /// Vector of stored coefficients.
    FloatArray mtrx;
    /// Integer array holding addresses of diagonal members.
    IntArray adr;
    /// Flag indicating whether factorized.
    int isFactorized;

public:
    /**
     * Constructor. Before any operation an internal profile must be built.
     * @see builInternalStructure
     */
    Skyline(int n=0);
    Skyline(const Skyline &s);
    Skyline(int n, FloatArray mtrx, IntArray adr);
    /// Destructor
    virtual ~Skyline() {}

    std::unique_ptr<SparseMtrx> clone() const override;

    void times(const FloatArray &x, FloatArray &answer) const override;
    void timesT(const FloatArray &x, FloatArray &answer) const override { this->times(x, answer); }
    void times(double x) override;
    void add(double x, SparseMtrx &m) override;
    int buildInternalStructure(EngngModel *, int, const UnknownNumberingScheme &) override;

    std::unique_ptr<SparseMtrx> giveSubMatrix(const IntArray &rows, const IntArray &cols) override; 
    /**
     * Allocates and builds internal structure according to given
     * array holding addresses of diagonal members values (adr).
     */
    int setInternalStructure(IntArray a);

    int assemble(const IntArray &loc, const FloatMatrix &mat) override;
    int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat) override;

    bool canBeFactorized() const override { return true; }
    SparseMtrx *factorized() override;
    FloatArray *backSubstitutionWith(FloatArray &) const override;
    void zero() override;
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

    double &at(int, int) override;
    double at(int i, int j) const override;
    bool isAllocatedAt(int i, int j) const override;
    int giveNumberOfNonZeros() const { return this->mtrx.giveSize(); }
    void toFloatMatrix(FloatMatrix &answer) const override;
    void printYourself() const override;
    void writeToFile(const char *fname) const override;
    int giveAllocatedSize() { return this->mtrx.giveSize(); }

    SparseMtrxType giveType() const override { return SMT_Skyline; }
    bool isAsymmetric() const override { return false; }

    const char *giveClassName() const override { return "Skyline"; }
};
} // end namespace oofem
#endif // skyline_h
