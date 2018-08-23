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

#ifndef spoolessparsemtrx_h
#define spoolessparsemtrx_h

#include "sparsemtrx.h"

extern "C" {
#include <spooles/misc.h>
#include <spooles/FrontMtx.h>
#include <spooles/SymbFac.h>
};

#define _IFT_SpoolesSparseMtrx_Name "spooles"

namespace oofem {
/**
 * This class provides an sparse matrix interface to SPOOLES InpMtrx
 */
class OOFEM_EXPORT SpoolesSparseMtrx : public SparseMtrx
{
protected:
    InpMtx *mtrx;
    int type;
    int sflag;

public:
    SpoolesSparseMtrx(int n=0, int m=0, int _sflag=SPOOLES_SYMMETRIC, int _type=SPOOLES_REAL) : SparseMtrx(n, m),
        mtrx(nullptr),
        type(_type),
        sflag(_sflag)
    { }
    virtual ~SpoolesSparseMtrx() {
        if ( mtrx ) {
            InpMtx_free(mtrx);
        }
    }

    void times(const FloatArray &x, FloatArray &answer) const override;
    void timesT(const FloatArray &x, FloatArray &answer) const override;
    void times(double x) override;
    int buildInternalStructure(EngngModel *eModel, int di, const UnknownNumberingScheme &s) override;
    int assemble(const IntArray &loc, const FloatMatrix &mat) override;
    int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat) override;
    bool canBeFactorized() const override { return false; }
    SparseMtrx *factorized() override { return nullptr; }
    FloatArray *backSubstitutionWith(FloatArray &y) const override { return nullptr; }
    void zero() override;
    double &at(int i, int j) override;
    double at(int i, int j) const override;
    void printStatistics() const override;
    void printYourself() const override;
    SparseMtrxType  giveType() const override { return SMT_SpoolesMtrx; }
    bool isAsymmetric() const override { return this->type == SPOOLES_NONSYMMETRIC; }

    // Exposed internals
    InpMtx *giveInpMtrx() { return this->mtrx; }
    int giveValueType() const { return type; }
    int giveSymmetryFlag() const { return sflag; }
};
} // end namespace oofem
#endif
