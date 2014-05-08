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
#define _IFT_PrescribedGradientBCWeak_DuplicateCornerNodes   "duplicatecornernodes"

namespace oofem {

class TractionElement;
class Line;
class IntegrationRule;

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

    virtual int giveNumberOfInternalDofManagers();
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

    virtual void assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, FloatArray *eNorm = NULL);

    virtual void assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);

    virtual void giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, EquationID eid, CharType type,
                               		const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);

    virtual void giveTractionLocationArrays(int iTracElInd, IntArray &rows, EquationID eid, CharType type,
                                    const UnknownNumberingScheme &s);

    virtual void giveDisplacementLocationArrays(int iTracElInd, IntArray &rows, EquationID eid, CharType type,
                                    const UnknownNumberingScheme &s);

    virtual const char *giveClassName() const { return "PrescribedGradientBCWeak"; }
    virtual const char *giveInputRecordName() const { return _IFT_PrescribedGradientBCWeak_Name; }

    /**
     * Computes the homogenized, macroscopic, field (stress).
     * @param sigma Output quantity (typically stress).
     * @param eid Equation ID to which sigma belongs.
     * @param tStep Active time step.
     */
    void computeField(FloatArray &sigma, EquationID eid, TimeStep *tStep);

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


    /**
     * 0 -> Do not duplicate corner traction nodes
     * 1 -> Duplicate corner traction nodes
     */
    bool mDuplicateCornerNodes;

    /// Lower corner of domain (assuming a rectangular RVE)
    FloatArray mLC;

    /// Upper corner of domain (assuming a rectangular RVE)
    FloatArray mUC;


    /// DOF-managers for the independent traction discretization
    std::vector<Node*> mpTractionNodes;

    /// Lock displacements in one node i periodic
    Node *mpDisplacementLock;
    double mDispLockScaling;

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
     * Map from a nodes global number to its local (zero based) index
     * in the traction element. One map for each traction element!
     * First index: Traction el index.
     * Second index: Global node number.
     */
    std :: unordered_map< int, std :: unordered_map< int, int > > mNodeTractionElLocalInd;

    /**
     * Map from a traction element to displacement node and
     * crack intersection coordinates inside the element.
     * This map is used when creating the integration rule.
     */
    std :: unordered_map<int, std::vector<FloatArray> > mTractionElInteriorCoordinates;


    std :: unordered_map<int, std::vector<int> > mTracElDispNodes;

    void createTractionMesh();

    void integrateTangent(FloatMatrix &oTangent, size_t iTracElInd);

    IntegrationRule *createNewIntegrationRule(int iTracElInd);

    void computeNTraction(FloatArray &oN, const double &iXi, const TractionElement &iEl) const;

    void giveTractionUnknows(FloatArray &oTracUnknowns, ValueModeType mode, TimeStep *tStep, int iTracElInd);
    void giveDisplacementUnknows(FloatArray &oDispUnknowns, ValueModeType mode, TimeStep *tStep, int iTracElInd);

    double domainSize();

    bool pointIsOnGammaPlus(const FloatArray &iPos) const;
    void giveMirroredPointOnGammaMinus(FloatArray &oPosMinus, const FloatArray &iPosPlus) const;
};

class TractionElement {
public:
    TractionElement() {}
    virtual ~TractionElement() {}

    void computeN_Constant(FloatArray &oN, const double &iXi) const {oN = {1.0};}
    void computeN_Linear(FloatArray &oN, const double &iXi) const {oN = {0.5*(1.0-iXi), 0.5*(1.0+iXi)};}

    std::vector<int> mTractionNodeInd;

    // Interior segments used for Gaussian quadrature
    std::vector<Line> mInteriorSegments;

    FloatArray mStartCoord;
    FloatArray mEndCoord;
};

class ArcPosSortFunction {
public:
    ArcPosSortFunction(const FloatArray &iStartPos):mStartPos(iStartPos) {};
    ~ArcPosSortFunction() {};

    bool operator()(const FloatArray &iVec1, const FloatArray &iVec2) const
    {
                return mStartPos.distance_square(iVec1) < mStartPos.distance_square(iVec2);
    }

private:
    const FloatArray mStartPos;
};


} /* namespace oofem */

#endif /* PRESCRIBEDGRADIENTBCWEAK_H_ */
