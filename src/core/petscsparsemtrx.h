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
#ifndef petscsparsemtrx_h
#define petscsparsemtrx_h

#include "sparsemtrx.h"

#include <petscksp.h>

#define _IFT_PetscSparseMtrx_Name "petsc"

namespace oofem {
/**
 * This class provides an sparse matrix interface to PETSc Matrices
 */
class OOFEM_EXPORT PetscSparseMtrx : public SparseMtrx
{
protected:
    Mat mtrx;
    bool symmFlag;
    MatType mType;
    int leqs;
    int geqs;
    int blocksize;
    int di;
    EngngModel *emodel;

    /// Linear solver context.
    KSP ksp;
    /// Flag if context initialized.
    bool kspInit;
    /// Flag if matrix has changed since last solve.
    bool newValues;

    /// Context or scattering/collecting parallel PETSc vectors
    IS localIS, globalIS;

public:
    PetscSparseMtrx(int n=0, int m=0);

    virtual ~PetscSparseMtrx();

    std::unique_ptr<SparseMtrx> clone() const override;
    void times(const FloatArray &x, FloatArray &answer) const override;
    void timesT(const FloatArray &x, FloatArray &answer) const override;
    void times(const FloatMatrix &B, FloatMatrix &answer) const override;
    void timesT(const FloatMatrix &B, FloatMatrix &answer) const override;
    void times(double x) override;
    void add(double x, SparseMtrx &m) override;
    void addDiagonal(double x, FloatArray &m) override;
    int buildInternalStructure(EngngModel *eModel, int n, int m, const IntArray &I, const IntArray &J) override;
    int buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme &s) override;
    int buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) override;
    int assemble(const IntArray &loc, const FloatMatrix &mat) override;
    int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat) override;
    int assembleBegin() override;
    int assembleEnd() override;
    std::unique_ptr<SparseMtrx> giveSubMatrix(const IntArray &rows, const IntArray &cols) override;
    bool canBeFactorized() const override { return false; }
    SparseMtrx *factorized() override { return NULL; }
    FloatArray *backSubstitutionWith(FloatArray &y) const override { return NULL; }
    void zero() override;
    double computeNorm() const override;
    double &at(int i, int j) override;
    double at(int i, int j) const override;
    void toFloatMatrix(FloatMatrix &answer) const override;
    void printStatistics() const override;
    void printYourself() const override;
    void printMatlab() const;
    SparseMtrxType giveType() const override;
    bool isAsymmetric() const override;
    void writeToFile(const char *fname) const override;
    const char *giveClassName() const override { return "PetscSparseMtrx"; }

    /// Creates a global vector that fits the instance of this matrix.
    void createVecGlobal(Vec *answer) const;
    /// Scatters global vector to natural one.
    int scatterG2L(Vec src, FloatArray &dest) const;
    /**
     * Scatters and gathers vector in local form to global (parallel) one.
     * Only local entries are processed.
     */
    int scatterL2G(const FloatArray &src, Vec dest) const;

    // Internals (should be documented)
    Mat *giveMtrx();
    bool giveSymmetryFlag() const;
    int setOption(MatOption op, PetscBool flag);
    int giveLeqs();
    int giveDomainIndex() const;

    friend class PetscSolver;
};
} // end namespace oofem
#endif
