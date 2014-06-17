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

#ifndef PRESCRIBEDGRADIENTBCWEAK_H_
#define PRESCRIBEDGRADIENTBCWEAK_H_

#include "prescribedgradientbc.h"

#define _IFT_PrescribedGradientBCWeak_Name   "prescribedgradientbcweak"
#define _IFT_PrescribedGradientBCWeak_TractionInterpOrder   "tractioninterporder"
#define _IFT_PrescribedGradientBCWeak_NumTractionNodesAtIntersections   "numnodesatintersections"
#define _IFT_PrescribedGradientBCWeak_NumTractionNodeSpacing   "tractionnodespacing"
#define _IFT_PrescribedGradientBCWeak_TractionOnGammaPlus   "periodic"

namespace oofem {

class TractionElement;

/*
 * Imposes a prescribed gradient weakly on the boundary
 * with an independent traction discretization.
 *
 * @author Erik Svenning
 * @date April 17, 2014
 */
class PrescribedGradientBCWeak : public PrescribedGradientBC {
public:
    PrescribedGradientBCWeak(int n, Domain * d);
    virtual ~PrescribedGradientBCWeak();

    virtual int giveNumberOfInternalDofManagers() {return mpTractionNodes.size();}
    virtual DofManager *giveInternalDofManager(int i);

    virtual bcType giveType() const { return UnknownBT; }

    /**
     * Initializes receiver according to object description stored in input record.
     * The input record contains two fields;
     * - devGradient \#columns { d_11 d_22 ... d_21 ... } (required)
     * - pressure p (required)
     * The gradient should be in Voigt notation (only the deviatoric part will be used)
     */
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual void postInitialize();

    virtual void scale(double s);

    virtual void assembleVector(FloatArray &answer, TimeStep *tStep,
                                CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, FloatArray *eNorm = NULL);

    virtual void assemble(SparseMtrx *answer, TimeStep *tStep,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);

    virtual void giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                                    const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);

    virtual const char *giveClassName() const { return "PrescribedGradientBCWeak"; }
    virtual const char *giveInputRecordName() const { return _IFT_PrescribedGradientBCWeak_Name; }

protected:

    // Options

    /// Order of interpolation for traction (0->piecewise constant, 1->piecewise linear)
    int mTractionInterpOrder;

    /**
     * If traction nodes should be inserted where cracks intersect
     * the RVE boundary.
     * 0 -> do not insert node.
     * 1 -> insert node.
     * 2 -> insert duplicated node.
     */
    int mNumTractionNodesAtIntersections;

    /**
     * Use every (mTractionNodeSpacing) displacement nodes when
     * constructing the traction element mesh.
     */
    int mTractionNodeSpacing;

    /**
     * true -> the traction lives only on gammaPlus, so that we
     * get strong periodicity as a special case.
     * false -> the traction lives everywhere on gamma, so
     * that we get Dirichlet as a special case
     */
    bool mTractionLivesOnGammaPlus;

    /// Lower corner of domain (assuming a rectangular RVE)
    FloatArray mLC;

    /// Upper corner of domain (assuming a rectangular RVE)
    FloatArray mUC;


    /// DOF-managers for the independent traction discretization
    std::vector<Node*> mpTractionNodes;

    /// Elements for the independent traction discretization
    std::vector<TractionElement*> mpTractionElements;


    /**
     * Map from a traction element to displacement elements it
     * interacts with everywhere on gamma.
     */
    std :: unordered_map< int, std :: vector< int > > mMapTractionElDispElGamma;

    /**
     * Map from a traction element to displacement elements it
     * interacts with on gammaPlus.
     */
    std :: unordered_map< int, std :: vector< int > > mMapTractionElDispElGammaPlus;

    /**
     * Map from a traction element to displacement elements it
     * interacts with on gammaMinus.
     */
    std :: unordered_map< int, std :: vector< int > > mMapTractionElDispElGammaMinus;

    /**
     * Map from a nodes global number to its local index
     * in the traction element.
     */
    std :: unordered_map< int, int > mNodeTractionElLocalInd;



    void createTractionMesh();
};

class TractionElement {
public:
    TractionElement() {}
    virtual ~TractionElement() {}
    std::vector<int> mTractionNodeInd;

    FloatArray mStartCoord;
    FloatArray mEndCoord;
};


} /* namespace oofem */

#endif /* PRESCRIBEDGRADIENTBCWEAK_H_ */
