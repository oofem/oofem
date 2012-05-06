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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef spoolessparsemtrx_h
#define spoolessparsemtrx_h

#ifdef __SPOOLES_MODULE

 #include "sparsemtrx.h"
 #include "spoolesinterface.h"

namespace oofem {
/**
 * This class provides an sparse matrix interface to SPOOLES InpMtrx
 */
class SpoolesSparseMtrx : public SparseMtrx
{
protected:
    InpMtx *mtrx;
    int type;
    int nent;
    int sflag;

public:
    SpoolesSparseMtrx(int _type, int _nent, int _sflag, int n, int m) : SparseMtrx(n, m)
    {
        type = _type;
        nent = _nent;
        sflag = _sflag;
        mtrx = NULL;
    }
    SpoolesSparseMtrx() : SparseMtrx() {
        type = SPOOLES_REAL;
        nent = 0;
        sflag = SPOOLES_SYMMETRIC;
        mtrx = NULL;
    }
    ~SpoolesSparseMtrx() { if ( mtrx ) { InpMtx_free(mtrx); } }

    // Overloaded methods
    virtual SparseMtrx *GiveCopy() const;
    virtual void times(const FloatArray &x, FloatArray &answer) const;
    virtual void timesT(const FloatArray &x, FloatArray &answer) const;
    virtual void times(double x);
    virtual int buildInternalStructure(EngngModel *eModel, int di, EquationID ut, const UnknownNumberingScheme &s);
    virtual int assemble(const IntArray &loc, const FloatMatrix &mat);
    virtual int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat);
    virtual bool canBeFactorized() const { return false; }
    virtual SparseMtrx *factorized() { return NULL; }
    virtual FloatArray *backSubstitutionWith(FloatArray &y) const { return NULL; }
    virtual void zero();
    virtual double &at(int i, int j);
    virtual double at(int i, int j) const;
    virtual void printStatistics() const;
    virtual void printYourself() const;
    virtual SparseMtrxType  giveType() const { return SMT_SpoolesMtrx; }
    virtual bool isAsymmetric() const { return this->type == SPOOLES_NONSYMMETRIC; }

    // Exposed internals
    InpMtx *giveInpMtrx() { return this->mtrx; }
    int giveValueType() const { return type; }
    int giveSymmetryFlag() const { return sflag; }

};
} // end namespace oofem
#endif
#endif // spoolessparsemtrx_h
