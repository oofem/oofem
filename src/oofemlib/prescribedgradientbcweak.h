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

#include "geometry.h"
#include "gausspoint.h"

#define _IFT_PrescribedGradientBCWeak_Name   "prescribedgradientbcweak"
#define _IFT_PrescribedGradientBCWeak_TractionInterpOrder   "tractioninterporder"
#define _IFT_PrescribedGradientBCWeak_NumTractionNodesAtIntersections   "numnodesatintersections"
#define _IFT_PrescribedGradientBCWeak_NumTractionNodeSpacing   "tractionnodespacing"
#define _IFT_PrescribedGradientBCWeak_DuplicateCornerNodes   "duplicatecornernodes"
#define _IFT_PrescribedGradientBCWeak_TangDistPadding   "tangdistpadding"

namespace oofem {

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

    virtual void assembleVector(FloatArray &answer, TimeStep *tStep,
                                CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, FloatArray *eNorm = NULL);

    virtual void assemble(SparseMtrx *answer, TimeStep *tStep,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);

    virtual void giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                               		const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);

    virtual void giveTractionLocationArray(IntArray &rows,
                                    const UnknownNumberingScheme &s);

    virtual void giveTractionLocationArrays(int iTracElInd, IntArray &rows, CharType type,
                                    const UnknownNumberingScheme &s);

    virtual void giveDisplacementLocationArrays(int iTracElInd, IntArray &rows, CharType type,
                                    const UnknownNumberingScheme &s);

    virtual const char *giveClassName() const { return "PrescribedGradientBCWeak"; }
    virtual const char *giveInputRecordName() const { return _IFT_PrescribedGradientBCWeak_Name; }

    /**
     * Computes the homogenized, macroscopic, field (stress).
     * @param sigma Output quantity (typically stress).
     * @param eid Equation ID to which sigma belongs.
     * @param tStep Active time step.
     */
    void computeField(FloatArray &sigma, TimeStep *tStep);

    /// Routines for postprocessing
    size_t giveNumberOfTractionElements() const {return mpTractionElements.size();}
    void giveTractionElCoord(size_t iElInd, FloatArray &oStartCoord, FloatArray &oEndCoord) const {oStartCoord = mpTractionElements[iElInd]->mStartCoord; oEndCoord = mpTractionElements[iElInd]->mEndCoord;}
    void giveTractionElNormal(size_t iElInd, FloatArray &oNormal, FloatArray &oTangent) const;

    void giveTraction(size_t iElInd, FloatArray &oStartTraction, FloatArray &oEndTraction, ValueModeType mode, TimeStep *tStep);

    // TODO: Consider moving this function to Domain.
    void computeDomainBoundingBox(Domain &iDomain, FloatArray &oLC, FloatArray &oUC);

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
    bool mMeshIsPeriodic;


    /**
     * 0 -> Do not duplicate corner traction nodes
     * 1 -> Duplicate corner traction nodes
     */
    bool mDuplicateCornerNodes;

    /**
     * Parameter for creation of traction mesh
     */
    double mTangDistPadding;

    /// Lower corner of domain (assuming a rectangular RVE)
    FloatArray mLC;

    /// Upper corner of domain (assuming a rectangular RVE)
    FloatArray mUC;


    /// DOF-managers for the independent traction discretization
    std::vector<Node*> mpTractionNodes;
    std::vector<Node*> mpTractionMasterNodes;

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
     * Map from a traction element to displacement node and
     * crack intersection coordinates inside the element.
     * This map is used when creating the integration rule.
     */
    std :: unordered_map<int, std::vector<FloatArray> > mTractionElInteriorCoordinates;


    std :: unordered_map<int, std::vector<int> > mTracElDispNodes;

    void createTractionMesh(bool iEnforceCornerPeriodicity);
    void buildMaps(const std::vector< std::pair<FloatArray, bool> > &iBndNodeCoordsFull);
    void createTractionElements(const std::vector<FloatArray> &iTractionNodeCoord, const double &iNodeDistTol);

    void integrateTangent(FloatMatrix &oTangent, size_t iTracElInd);

    void assembleTangentGPContribution(FloatMatrix &oTangent, size_t iTracElInd, GaussPoint &iGP, const FloatArray &iBndCoord, std :: unordered_map<int, IntArray > &iGlobalNodeIndToPosInLocalLocArray, const double &iScaleFactor);

    IntegrationRule *createNewIntegrationRule(int iTracElInd);

    void computeNTraction(FloatArray &oN, const double &iXi, const TractionElement &iEl) const;

    void giveTractionUnknows(FloatArray &oTracUnknowns, ValueModeType mode, TimeStep *tStep, int iTracElInd);
    void giveDisplacementUnknows(FloatArray &oDispUnknowns, ValueModeType mode, TimeStep *tStep, int iTracElInd);

    double domainSize();

    bool pointIsOnGammaPlus(const FloatArray &iPos) const;
    void giveMirroredPointOnGammaMinus(FloatArray &oPosMinus, const FloatArray &iPosPlus) const;
    void giveMirroredPointOnGammaPlus(FloatArray &oPosPlus, const FloatArray &iPosMinus) const;
    bool pointIsMapapble(const FloatArray &iPos) const;

    virtual void giveBoundaryCoordVector(FloatArray &oX, const FloatArray &iPos) const = 0;
    virtual void checkIfCorner(bool &oIsCorner, bool &oDuplicatable, const FloatArray &iPos, const double &iNodeDistTol) const = 0;
    virtual bool boundaryPointIsOnActiveBoundary(const FloatArray &iPos) const = 0;
};

class ArcPosSortFunction {
public:
    ArcPosSortFunction(const FloatArray &iStartPos):mStartPos(iStartPos) {}
    ~ArcPosSortFunction() {}

    bool operator()(const FloatArray &iVec1, const FloatArray &iVec2) const
    {
                return mStartPos.distance_square(iVec1) < mStartPos.distance_square(iVec2);
    }

private:
    const FloatArray mStartPos;
};

template <class T>
class ArcPosSortFunction2 {
public:
	ArcPosSortFunction2(const FloatArray &iLC, const FloatArray &iUC, const double &iTol):
	mLC(iLC),
	mUC(iUC),
	mTol(iTol)
	{

	}

	~ArcPosSortFunction2() {}

    bool operator()(const std::pair<FloatArray, T> &iVec1, const std::pair<FloatArray,int> &iVec2) const
    {
                return calcArcPos(iVec1.first) < calcArcPos(iVec2.first);
    }

	double calcArcPos(const FloatArray &iPos) const
	{
		double Lx = mUC[0] - mLC[0];
		double Ly = mUC[1] - mLC[1];

		if( iPos[1] < (mLC[1]+mTol) ){
			// Edge 1
			double dist = iPos.distance(mLC);
			return dist;
		}

		if( iPos[0] > (mUC[0]-mTol) ){
			// Edge 2
			const FloatArray &x = {mUC[0], mLC[1]};
			double dist = Lx + iPos.distance(x);
			return dist;
		}

		if( iPos[1] > (mUC[1]-mTol) ){
			// Edge 3
			double dist = Lx + Ly + iPos.distance(mUC);
			return dist;
		}

		if( iPos[0] < (mLC[0]+mTol) ){
			// Edge 4
			const FloatArray &x = {mLC[0], mUC[1]};
			double dist = Lx + Ly + Lx + iPos.distance(x);
			return dist;
		}

		OOFEM_ERROR("Could not compute distance.")
		return 0.0;
	}

private:
	const FloatArray mLC;
	const FloatArray mUC;
	const double mTol;
};


class EdgeTracker {
public:
	EdgeTracker()
	{
		mNodePair.first 	= -1;
		mNodePair.second 	= -1;
		mStoredSecondIndex = -1;
		mIsFirstPoint = true;
	}

	~EdgeTracker() {}

	void addPoint(int iPointIndex, bool iIsDuplicated, bool iPeriodic)
	{
		if( iIsDuplicated || (iPeriodic && mIsFirstPoint) ) {
			// Will be added in first position (i.e. master node)

			if( mNodePair.first == -1 && mNodePair.second == -1 ) {
				mNodePair.first = iPointIndex;
			}
			else {
//				printf("Pair with slave %d and master %d\n", mNodePair.second, mNodePair.first);
				mStoredSecondIndex = mNodePair.first;
//				OOFEM_ERROR("The node pair has already been initialized.")
			}

			mIsFirstPoint = false;
		}
		else {
			// Will be added in second position (i.e. slave node)

			// Add to second position if we have a pair where
			// the first position has been set. Otherwise,
			// store it for later use.

			if( mNodePair.first != -1 ) {
				mNodePair.second = iPointIndex;

				// We prefer to have the highest index as slave
				if( mNodePair.first > mNodePair.second ) {
					std::swap( mNodePair.first, mNodePair.second );
				}

				mMapSlaveToMaster[mNodePair.second] = mNodePair.first;

//				printf("Creating pair with slave %d and master %d\n", mNodePair.second, mNodePair.first);

				mNodePair.first 	= -1;
				mNodePair.second 	= -1;
			}
			else {
				mStoredSecondIndex = iPointIndex;
			}

		}
	}

	bool giveMasterIndex(int &oMasterIndex, int iSlaveIndex) const
	{
		auto it = mMapSlaveToMaster.find(iSlaveIndex);
		if( it == mMapSlaveToMaster.end() ) {
			return false;
		}

		oMasterIndex = it->second;
		return true;
	}

	void addStoredSecondIndex()
	{
		addPoint(mStoredSecondIndex, false, false);
	}

	std::pair<int,int> mNodePair;
	std::unordered_map<int, int> mMapSlaveToMaster;

	int mStoredSecondIndex;
	bool mIsFirstPoint;
};


} /* namespace oofem */

#endif /* PRESCRIBEDGRADIENTBCWEAK_H_ */
