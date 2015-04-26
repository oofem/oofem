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
    PetscSparseMtrx(int n, int m);
    PetscSparseMtrx();

    virtual ~PetscSparseMtrx();

    // Overloaded methods:
    virtual SparseMtrx *GiveCopy() const;
    virtual void times(const FloatArray &x, FloatArray &answer) const;
    virtual void timesT(const FloatArray &x, FloatArray &answer) const;
    virtual void times(const FloatMatrix &B, FloatMatrix &answer) const;
    virtual void timesT(const FloatMatrix &B, FloatMatrix &answer) const;
    virtual void times(double x);
    virtual void add(double x, SparseMtrx &m);
    virtual void addDiagonal(double x, FloatArray &m);
    virtual int buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme &s);
    virtual int buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);
    virtual int assemble(const IntArray &loc, const FloatMatrix &mat);
    virtual int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat);
    virtual int assembleBegin();
    virtual int assembleEnd();
    virtual bool canBeFactorized() const { return false; }
    virtual SparseMtrx *factorized() { return NULL; }
    virtual FloatArray *backSubstitutionWith(FloatArray &y) const { return NULL; }
    virtual void zero();
    virtual double computeNorm() const;
    virtual double &at(int i, int j);
    virtual double at(int i, int j) const;
    virtual void toFloatMatrix(FloatMatrix &answer) const;
    virtual void printStatistics() const;
    virtual void printYourself() const;
    void printMatlab() const;
    virtual SparseMtrxType  giveType() const;
    virtual bool isAsymmetric() const;
    virtual void writeToFile(const char *fname) const;
    virtual const char *giveClassName() const { return "PetscSparseMtrx"; }

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
