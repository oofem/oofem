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

#ifndef skylineu_h
#define skylineu_h

#include "sparsemtrx.h"
#include "rowcol.h"

namespace oofem {
/// "zero" pivot for SkylineUnsym class
#define SkylineUnsym_TINY_PIVOT 1.e-30

/**
 * This class implements a nonsymmetric matrix stored in a compacted
 * (skyline) form. Its shape is symmetric, but not its coefficients.
 * The Skyline contains 'size' row-column segments, each of them of any size;
 * these are stored in the attribute 'rowColumns'.
 */
class OOFEM_EXPORT SkylineUnsym : public SparseMtrx
{
protected:
    /// Row column segments
    RowColumn **rowColumns;
    /// Size of receiver
    int size;
    /// Factorization flag
    int isFactorized;

public:
    /**
     * Constructor. Before any operation an internal profile must be built.
     * @param n Size of matrix
     * @see buildInternalStructure
     */
    SkylineUnsym(int n);
    /**
     * Constructor. Before any operation an internal profile must be built.
     * @see buildInternalStructure
     */
    SkylineUnsym();
    /// Destructor
    virtual ~SkylineUnsym();

    // Overloaded methods:
    virtual SparseMtrx *GiveCopy() const;
    virtual void times(const FloatArray &x, FloatArray &answer) const;
    virtual void timesT(const FloatArray &x, FloatArray &answer) const;
    virtual void times(double x);
    virtual int buildInternalStructure(EngngModel *, int, const UnknownNumberingScheme &s);
    int setInternalStructure(IntArray &a);
    virtual int assemble(const IntArray &loc, const FloatMatrix &mat);
    virtual int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat);
    virtual bool canBeFactorized() const { return true; }
    virtual SparseMtrx *factorized();
    virtual FloatArray *backSubstitutionWith(FloatArray &) const;
    virtual void zero();
    virtual double &at(int i, int j);
    virtual double at(int i, int j) const;
    virtual void toFloatMatrix(FloatMatrix &answer) const;
    virtual void printYourself() const;
    virtual void printStatistics() const;
    virtual void writeToFile(const char *fname) const;
    virtual SparseMtrxType giveType() const { return SMT_SkylineU; }
    virtual bool isAsymmetric() const { return true; }
    virtual const char *giveClassName() const { return "SkylineU"; }

protected:
    void checkSizeTowards(const IntArray &);
    void checkSizeTowards(const IntArray &rloc, const IntArray &cloc);
    RowColumn *giveRowColumn(int j) const;
    void growTo(int);

    SkylineUnsym(RowColumn * *, int, int);
};
} // end namespace oofem
#endif // skylineu_h
