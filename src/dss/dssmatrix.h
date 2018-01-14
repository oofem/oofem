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
#include "DSSolver.h"
#include "SparseMatrixF.h"

#include <memory>

#define _IFT_DSSMatrix_Name "dss"

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


    /// Pointer ro row indices
    std:: unique_ptr<unsigned long[]> rowind_;
    /// Pointer to columns
    std:: unique_ptr<unsigned long[]> colptr_;
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
    DSSMatrix(dssType t, int n=0);
    /// Copy constructor
    DSSMatrix(const DSSMatrix &S);
    /// Destructor
    virtual ~DSSMatrix() {}

    void times(const FloatArray &x, FloatArray &answer) const override;
    void times(double x) override;
    int buildInternalStructure(EngngModel *, int, const UnknownNumberingScheme & s) override;
    int assemble(const IntArray &loc, const FloatMatrix &mat) override;
    int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat) override;
    bool canBeFactorized() const override { return true; }
    SparseMtrx *factorized() override;
    void solve(FloatArray &b, FloatArray &x);
    void zero() override;
    double &at(int i, int j) override;
    double at(int i, int j) const override;
    SparseMtrxType giveType() const override { return SMT_SymCompCol; }
    bool isAsymmetric() const override { return false; }

    const char *giveClassName() const override { return "DSSMatrix"; }
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

