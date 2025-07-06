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

#define _IFT_SkylineUnsym_Name "skylineu"


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
    std::vector<RowColumn> rowColumns;
    /// Factorization flag
    int isFactorized;

public:
    /**
     * Constructor. Before any operation an internal profile must be built.
     * @param n Size of matrix
     * @see buildInternalStructure
     */
    SkylineUnsym(int n=0);
    SkylineUnsym(const SkylineUnsym &s);
    /// Destructor
    virtual ~SkylineUnsym() {}

    // Overloaded methods:
    std::unique_ptr<SparseMtrx> clone() const override;
    void times(const FloatArray &x, FloatArray &answer) const override;
    void timesT(const FloatArray &x, FloatArray &answer) const override;
    void times(double x) override;
    int buildInternalStructure(EngngModel *, int, const UnknownNumberingScheme &s) override;
    int setInternalStructure(IntArray &a);
    int assemble(const IntArray &loc, const FloatMatrix &mat) override;
    int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat) override;
    bool canBeFactorized() const override { return true; }
    SparseMtrx *factorized() override;
    FloatArray *backSubstitutionWith(FloatArray &) const override;
    void zero() override;
    double &at(int i, int j) override;
    double at(int i, int j) const override;
    void toFloatMatrix(FloatMatrix &answer) const override;
    void printYourself() const override;
    void printStatistics() const override;
    void writeToFile(const char *fname) const override;
    SparseMtrxType giveType() const override { return SMT_SkylineU; }
    bool isAsymmetric() const override { return true; }
    const char *giveClassName() const override { return "SkylineU"; }

protected:
    void checkSizeTowards(const IntArray &);
    void checkSizeTowards(const IntArray &rloc, const IntArray &cloc);
    void growTo(int);

    SkylineUnsym(int n, std::vector<RowColumn> data, bool isFact);
};
} // end namespace oofem
#endif // skylineu_h
