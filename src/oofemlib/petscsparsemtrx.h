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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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
  #include <petscksp.h>
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

    /// Linear solver context.
    KSP ksp;
    /// Flag if context initialized.
    bool kspInit;
    /// Flag if matrix has changed since last solve.
    bool newValues;

public:
    PetscSparseMtrx(int n, int m) : SparseMtrx(n, m),
        mtrx(NULL), symmFlag(false), leqs(0), geqs(0), di(0), kspInit(false), newValues(true) {}
    PetscSparseMtrx() : SparseMtrx(),
        mtrx(NULL), symmFlag(false), leqs(0), geqs(0), di(0), kspInit(false), newValues(true)  {}

    ~PetscSparseMtrx() {
        MatDestroy(this->mtrx);
        if ( this->kspInit ) {
            KSPDestroy(this->ksp);
        }
    }

    // Overloaded methods:
    virtual SparseMtrx *GiveCopy() const;
    virtual void times(const FloatArray &x, FloatArray &answer) const;
    virtual void timesT(const FloatArray &x, FloatArray &answer) const;
    virtual void times(const FloatMatrix &B, FloatMatrix &answer) const;
    virtual void timesT(const FloatMatrix &B, FloatMatrix &answer) const;
    virtual void times(double x);
    virtual int buildInternalStructure(EngngModel *eModel, int di, EquationID, const UnknownNumberingScheme & s);
    virtual int buildInternalStructure(EngngModel *eModel, int di, EquationID eid, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);
    virtual int assemble(const IntArray &loc, const FloatMatrix &mat);
    virtual int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat);
    int assembleBegin();
    int assembleEnd();
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
    virtual SparseMtrxType  giveType() const { return SMT_PetscMtrx; }
    virtual bool isAsymmetric() const { return !symmFlag; }

    // Internals (should be documented)
    Mat *giveMtrx() { return & this->mtrx; }
    bool giveSymmetryFlag() const { return symmFlag; }
    int setOption(MatOption op, PetscTruth flag) { return MatSetOption(this->mtrx, op, flag); }
    EquationID giveEquationID() { return ut; }
    int giveLeqs() { return leqs; }
    int giveDomainIndex() { return di; }

    friend class PetscSolver;
    friend class PETScSNESNM;
};
} // end namespace oofem
#endif
#endif // petscsparsemtrx_h
