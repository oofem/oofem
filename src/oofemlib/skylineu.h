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

//   ***********************************
//   *** CLASS SKYLINE (UNSYMMETRIC) ***
//   ***********************************

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
class SkylineUnsym : public SparseMtrx
{
protected:
    /// Row column segments
    RowColumn **rowColumns;
    /// Size of receiver
    int size;
    /// Factorization flag
    int isFactorized;

public:
    /** Constructor. Before any operation an internal profile must be built.
     * @param n Size of matrix
     * @see buildInternalStructure
     */
    SkylineUnsym(int n);
    /** Constructor. Before any operation an internal profile must be built.
     * @see buildInternalStructure
     */
    SkylineUnsym();
    /// Destructor
    ~SkylineUnsym();

    // Overloaded methods:
    SparseMtrx *GiveCopy() const;
    void times(const FloatArray &x, FloatArray &answer) const;
    void timesT(const FloatArray &x, FloatArray &answer) const;
    virtual void times(double x);
    int buildInternalStructure(EngngModel *, int, EquationID, const UnknownNumberingScheme & s);
    int setInternalStructure(IntArray *a);
    int assemble(const IntArray &loc, const FloatMatrix &mat);
    int assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat);
    bool canBeFactorized() const { return true; }
    SparseMtrx *factorized();
    FloatArray *backSubstitutionWith(FloatArray &) const;
    void zero();
    double &at(int i, int j);
    double at(int i, int j) const;
    void toFloatMatrix(FloatMatrix &answer) const;
    virtual void printYourself() const;
    virtual void printStatistics() const;
    virtual void writeToFile(const char* fname) const;
    SparseMtrxType giveType() const { return SMT_SkylineU; }
    bool isAsymmetric() const { return true; }

protected:
    /*
     *    FloatMatrix*  AsFloatMatrix () ;
     *    double&       at (int i,int j) ;
     *    void          assemble (FloatMatrix*,IntArray*) ;
     *    void          assemble (FloatMatrix*,IntArray*,IntArray*) ;
     *    FloatArray*   backSubstitutionWith (FloatArray*) ;
     *    void          carveYourselfFor (Domain*) ;
     *    int           computeNumberNegativeEigenvalue();
     *    void          checkSizeTowards (IntArray*) ;
     *    Skyline*      diagonalScalingWith (FloatArray*) ;
     *    Skyline*      factorized () ;
     *    Skyline*      forwardReductionWith (FloatArray*) ;
     */
    void checkSizeTowards(const IntArray &);
    void checkSizeTowards(const IntArray &rloc, const IntArray &cloc);
    RowColumn *giveRowColumn(int j) const;
    void growTo(int);
    /*
     * void          printYourself () ;
     * Skyline*      reinitialized () ;
     * SkylineSym*   giveSymmetricPart();
     * int       computeNumberNegativePivotsOfSymPart();
     */

    SkylineUnsym(RowColumn **, int, int);
};
} // end namespace oofem
#endif // skylineu_h
