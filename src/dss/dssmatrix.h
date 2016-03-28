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

#ifndef dss_h
#define dss_h

#include "sparsemtrx.h"
#include "intarray.h"
#include "floatarray.h"

#include <memory>

/* DSS module lives in global namespace, not in oofem namespace */
class DSSolver;
struct SparseMatrixF;

namespace oofem {
/**
 * Interface to Direct Sparse Solver written by R.Vonracek.
 * This class represent the sparse matrix interface to DSS library. It allows to build internal structure,
 * assemble the DSS sparse matrix, and to factorize and back substitution operations.
 */
class OOFEM_EXPORT DSSMatrix : public SparseMtrx
{
public:
    /// Possible storage schemes and factorization types
    enum dssType { sym_LDL, sym_LL, unsym_LU };

protected:
    /// Pointer to SparseMatrixF class rep
    std :: unique_ptr< SparseMatrixF > _sm;
    /// pointer to DSSolver class (representation of Direct Sparse Solver in DSS lib)
    std :: unique_ptr< DSSolver > _dss;
    /// Flag indicating whether factorized.
    bool isFactorized;
    /// type of storage & factorization
    dssType _type;

    /// implements 0-based access
    double operator()(int i, int j) const;
    /// implements 0-based access
    double &operator()(int i, int j);

public:
    /**
     * Constructor. Before any operation an internal profile must be built.
     * @param t Storage type
     * @param n Size of row and columns of square matrix
     * @see buildInternalStructure
     */
    DSSMatrix(dssType t, int n);
    /**
     * Constructor. Before any operation an internal profile must be built.
     * @param t Storage type
     * @see buildInternalStructure
     */
    DSSMatrix(dssType t);
    /// Copy constructor
    DSSMatrix(const DSSMatrix &S);
    /// Destructor
    virtual ~DSSMatrix();

    // Overloaded methods
    virtual SparseMtrx *GiveCopy() const;
    virtual void times(const FloatArray &x, FloatArray &answer) const;
    virtual void times(double x);
    virtual int buildInternalStructure(EngngModel *, int, const UnknownNumberingScheme & s);
    virtual int assemble(const IntArray &loc, const FloatMatrix &mat);
    virtual int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat);
    virtual bool canBeFactorized() const { return true; }
    virtual SparseMtrx *factorized();
    virtual void solve(FloatArray &b, FloatArray &x);
    virtual void zero();
    virtual double &at(int i, int j);
    virtual double at(int i, int j) const;
    virtual SparseMtrxType giveType() const { return SMT_SymCompCol; }
    virtual bool isAsymmetric() const { return false; }

    virtual const char *giveClassName() const { return "DSSMatrix"; }
};

class DSSMatrixLDL : public DSSMatrix
{
public:
    DSSMatrixLDL() : DSSMatrix(sym_LDL) {}
};

class DSSMatrixLL : public DSSMatrix
{
public:
    DSSMatrixLL() : DSSMatrix(sym_LL) {}
};

class DSSMatrixLU : public DSSMatrix
{
public:
    DSSMatrixLU() : DSSMatrix(unsym_LU) {}
};

} // end namespace oofem

#endif // dss_h

