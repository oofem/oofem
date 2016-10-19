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

#include "prescribedgradienthomogenization.h"
#include "activebc.h"
#include "geometry.h"
#include "dofiditem.h"

#include <unordered_map>
#include <memory>
#include "node.h"

#define _IFT_PrescribedGradientBCWeak_Name   "prescribedgradientbcweak"
#define _IFT_PrescribedGradientBCWeak_TractionInterpOrder   "tractioninterporder"
#define _IFT_PrescribedGradientBCWeak_NumTractionNodesAtIntersections   "numnodesatintersections"
#define _IFT_PrescribedGradientBCWeak_NumTractionNodeSpacing   "tractionnodespacing"
#define _IFT_PrescribedGradientBCWeak_DuplicateCornerNodes   "duplicatecornernodes"
#define _IFT_PrescribedGradientBCWeak_TangDistPadding   "tangdistpadding"
#define _IFT_PrescribedGradientBCWeak_TracDofScaling   "tracdofscaling"
#define _IFT_PrescribedGradientBCWeak_PeriodicityNormal   "periodicitynormal"
#define _IFT_PrescribedGradientBCWeak_MirrorFunction   "mirrorfunction"

namespace oofem {
class IntegrationRule;
class Node;
class GaussPoint;

class TracSegArray
{
public:
	TracSegArray() {}
    virtual ~TracSegArray() {}

    void printYourself() {
    	printf("\nTracSegArray segments:\n");
    	for(auto &l: mInteriorSegments) {
    		printf("\n");
    		l.giveVertex(1).printYourself();
    		l.giveVertex(2).printYourself();
    	}
    }

    double giveLength() {
    	double l = 0.0;
    	for(Line &line : mInteriorSegments) {
    		l += line.giveLength();
    	}

    	return l;
    }

    void giveTractionLocationArray(IntArray &rows, CharType type, const UnknownNumberingScheme &s);

    void setupIntegrationRuleOnEl();

    std :: vector< Line >mInteriorSegments;

    // Interior segments used for Gaussian quadrature
    std :: vector< Line > mInteriorSegmentsFine;

    std :: vector< FloatArray > mInteriorSegmentsPointsFine;



    std :: unique_ptr< Node > mFirstNode;

    std :: unique_ptr< IntegrationRule > mIntRule;
};

/**
 * Imposes a prescribed gradient weakly on the boundary
 * with an independent traction discretization.
 *
 * @author Erik Svenning
 * @date April 17, 2014
 */
class PrescribedGradientBCWeak : public ActiveBoundaryCondition, public PrescribedGradientHomogenization
{
public:
    PrescribedGradientBCWeak(int n, Domain *d);
    virtual ~PrescribedGradientBCWeak();

    void clear();

    virtual double domainSize() {return PrescribedGradientHomogenization::domainSize(this->giveDomain(), this->giveSetNumber());}

    virtual int giveNumberOfInternalDofManagers();
    virtual DofManager *giveInternalDofManager(int i);

    virtual bcType giveType() const { return UnknownBT; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual void postInitialize();

    virtual void computeField(FloatArray &sigma, TimeStep *tStep);
    virtual void computeTangent(FloatMatrix &E, TimeStep *tStep);

    virtual void assembleVector(FloatArray &answer, TimeStep *tStep,
                                CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, FloatArray *eNorm = NULL);

    void computeExtForceElContrib(FloatArray &oContrib, TracSegArray &iEl, int iDim, TimeStep *tStep);
    void computeIntForceGPContrib(FloatArray &oContrib_disp, IntArray &oDisp_loc_array, FloatArray &oContrib_trac, IntArray &oTrac_loc_array,TracSegArray &iEl, GaussPoint &iGP, int iDim, TimeStep *tStep, const FloatArray &iBndCoord, const double &iScaleFac, ValueModeType mode, CharType type, const UnknownNumberingScheme &s);


    virtual void assemble(SparseMtrx &answer, TimeStep *tStep,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);

    virtual void assembleExtraDisplock(SparseMtrx &answer, TimeStep *tStep,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);

    virtual void assembleGPContrib(SparseMtrx &answer, TimeStep *tStep,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, TracSegArray &iEl, GaussPoint &iGP);

    virtual void giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                                    const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);

    virtual void giveTractionLocationArray(IntArray &rows,
                                           const UnknownNumberingScheme &s);

    virtual void giveDisplacementLocationArray(IntArray &rows,
                                           const UnknownNumberingScheme &s);

    void compute_x_times_N_1(FloatMatrix &o_x_times_N);
    void compute_x_times_N_2(FloatMatrix &o_x_times_N);


//    virtual void giveTractionLocationArrays(int iTracElInd, IntArray &rows, CharType type,
//                                            const UnknownNumberingScheme &s);
//
//    virtual void giveDisplacementLocationArrays(int iTracElInd, IntArray &rows, CharType type,
//                                                const UnknownNumberingScheme &s);

    virtual const char *giveClassName() const { return "PrescribedGradientBCWeak"; }
    virtual const char *giveInputRecordName() const { return _IFT_PrescribedGradientBCWeak_Name; }

    // Routines for postprocessing
    size_t giveNumberOfTractionElements() const { return mpTracElNew.size(); }
    void giveTractionElCoord(size_t iElInd, FloatArray &oStartCoord, FloatArray &oEndCoord) const { oStartCoord = mpTracElNew [ iElInd ]->mInteriorSegments[0].giveVertex(1); oEndCoord = mpTracElNew [ iElInd ]->mInteriorSegments.back().giveVertex(2); }
    void giveTractionElNormal(size_t iElInd, FloatArray &oNormal, FloatArray &oTangent) const;
    void giveTractionElArcPos(size_t iElInd, double &oXiStart, double &oXiEnd) const;
    void giveBoundaries(IntArray &oBoundaries);

    void giveTraction(size_t iElInd, FloatArray &oStartTraction, FloatArray &oEndTraction, ValueModeType mode, TimeStep *tStep);

    // TODO: Consider moving this function to Domain.
    void computeDomainBoundingBox(Domain &iDomain, FloatArray &oLC, FloatArray &oUC);


    const IntArray &giveTracDofIDs() const {return mTractionDofIDs;}
    const IntArray &giveDispLockDofIDs() const {return mDispLockDofIDs;}
    const IntArray &giveRegularDispDofIDs() const {return mRegularDispDofIDs;}


    // Functions mainly for testing
    void setPeriodicityNormal(const FloatArray &iPeriodicityNormal) {mPeriodicityNormal = iPeriodicityNormal; };
    void setDomainSize(double iDomainSize) {mDomainSize = std::move(iDomainSize);};
    void setLowerCorner(FloatArray iLC) {mLC = std::move(iLC);};
    void setUpperCorner(FloatArray iUC) {mUC = std::move(iUC);};

    void setMirrorFunction(int iMirrorFunction) {mMirrorFunction = iMirrorFunction;};

protected:

    const IntArray mTractionDofIDs;
    const IntArray mDispLockDofIDs;
    const IntArray mRegularDispDofIDs;


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

    double mTracDofScaling;

    /// Lower corner of domain (assuming a rectangular RVE)
    FloatArray mLC;

    /// Upper corner of domain (assuming a rectangular RVE)
    FloatArray mUC;


    /// Lock displacements in one node if periodic
    Node *mpDisplacementLock;
    int mLockNodeInd;
    double mDispLockScaling;

    int mSpringNodeInd1;
    int mSpringNodeInd2;
    int mSpringNodeInd3;
    double mSpringPenaltyStiffness;


    /// Elements for the independent traction discretization
    std :: vector< TracSegArray * > mpTracElNew;


    /**
     * Periodicity direction.
     */
    FloatArray mPeriodicityNormal;

    double mDomainSize;

    /**
     * Mirror function (i.e. mapping between gamma^+ and gamma^-).
     * 0 -> Standard periodicity, 1 -> Shifted stacking, 2 -> Rotation
     */
    int mMirrorFunction;

public:
    void recomputeTractionMesh();

    void giveMirroredPointOnGammaMinus(FloatArray &oPosMinus, const FloatArray &iPosPlus) const;
    void giveMirroredPointOnGammaPlus(FloatArray &oPosPlus, const FloatArray &iPosMinus) const;

protected:
    void createTractionMesh(bool iEnforceCornerPeriodicity, int iNumSides);

    void splitSegments(std :: vector< TracSegArray * > &ioElArray);

    bool damageExceedsTolerance(Element *el);

    void assembleTangentGPContributionNew(FloatMatrix &oTangent, TracSegArray &iEl, GaussPoint &iGP, const double &iScaleFactor, const FloatArray &iBndCoord);

    bool pointIsOnGammaPlus(const FloatArray &iPos) const;

    virtual void giveBoundaryCoordVector(FloatArray &oX, const FloatArray &iPos) const = 0;
    virtual void checkIfCorner(bool &oIsCorner, bool &oDuplicatable, const FloatArray &iPos, const double &iNodeDistTol) const = 0;
    virtual bool boundaryPointIsOnActiveBoundary(const FloatArray &iPos) const = 0;

    int giveSideIndex(const FloatArray &iPos) const;


    void findHoleCoord(std::vector<FloatArray> &oHoleCoordUnsorted, std::vector<FloatArray> &oAllCoordUnsorted);
    void findCrackBndIntersecCoord(std::vector<FloatArray> &oHoleCoordUnsorted);
    void findPeriodicityCoord(std::vector<FloatArray> &oHoleCoordUnsorted);

    void removeClosePoints(std::vector<FloatArray> &ioCoords, const double &iAbsTol);
    void removeSegOverHoles(TracSegArray &ioTSeg, const double &iAbsTol);
};

class ArcPosSortFunction
{
public:
    ArcPosSortFunction(const FloatArray &iStartPos) : mStartPos(iStartPos) {}
    ~ArcPosSortFunction() {}

    bool operator()(const FloatArray &iVec1, const FloatArray &iVec2) const
    {
        return mStartPos.distance_square(iVec1) < mStartPos.distance_square(iVec2);
    }

private:
    const FloatArray mStartPos;
};

template< class T >
class ArcPosSortFunction3
{
public:
    ArcPosSortFunction3(const FloatArray &iLC, const FloatArray &iUC, const double &iTol, int iSideInd) :
        mLC(iLC),
        mUC(iUC),
        mTol(iTol),
        mSideInd(iSideInd)
    {}

    ~ArcPosSortFunction3() {}

    bool operator()(const std :: pair< FloatArray, T > &iVec1, const std :: pair< FloatArray, int > &iVec2) const
    {
        return calcArcPos(iVec1.first) < calcArcPos(iVec2.first);
    }

    double calcArcPos(const FloatArray &iPos) const
    {
        double Lx = mUC [ 0 ] - mLC [ 0 ];
        double Ly = mUC [ 1 ] - mLC [ 1 ];

        if ( mSideInd == 0 ) {
            const FloatArray &x = { mUC [ 0 ], mLC [ 1 ] };
            double dist = Lx + iPos.distance(x);
            return dist;
        }

        if ( mSideInd == 1 ) {
            double dist = Lx + Ly + iPos.distance(mUC);
            return dist;
        }

        if ( mSideInd == 2 ) {
            const FloatArray &x = { mLC [ 0 ], mUC [ 1 ] };
            double dist = Lx + Ly + Lx + iPos.distance(x);
            return dist;
        }

        if ( mSideInd == 3 ) {
            double dist = iPos.distance(mLC);
            return dist;
        }

        OOFEM_ERROR("Could not compute distance.")
        return 0.0;
    }

private:
    const FloatArray mLC;
    const FloatArray mUC;
    const double mTol;

    // 0->x=L, 1->y=L, 2->x=0, 3->y=0
    const int mSideInd;
};


class ArcPosSortFunction4
{
public:
    ArcPosSortFunction4(const FloatArray &iLC, const FloatArray &iUC, const double &iRelTol) :
        mLC(iLC),
        mUC(iUC),
        mRelTol(iRelTol)
    {}

    ~ArcPosSortFunction4() {}

    bool operator()(const FloatArray &iVec1, const FloatArray &iVec2) const
    {
        return calcArcPos(iVec1) < calcArcPos(iVec2);
    }

    double calcArcPos(const FloatArray &iPos) const
    {
        double Lx = mUC [ 0 ] - mLC [ 0 ];
        double Ly = mUC [ 1 ] - mLC [ 1 ];


        int sideInd = -1;

        if( iPos[0] > Lx - Lx*mRelTol ) {
        	sideInd = 0;
        }

        if( iPos[1] > Ly - Ly*mRelTol ) {
        	sideInd = 1;
        }

        if( iPos[0] < Lx*mRelTol ) {
        	sideInd = 2;
        }

        if( iPos[1] < Ly*mRelTol ) {
        	sideInd = 3;
        }

        if ( sideInd == 0 ) {
            const FloatArray &x = { mUC [ 0 ], mLC [ 1 ] };
            double dist = Lx + iPos.distance(x);
            return dist;
        }

        if ( sideInd == 1 ) {
            double dist = Lx + Ly + iPos.distance(mUC);
            return dist;
        }

        if ( sideInd == 2 ) {
            const FloatArray &x = { mLC [ 0 ], mUC [ 1 ] };
            double dist = Lx + Ly + Lx + iPos.distance(x);
            return dist;
        }

        if ( sideInd == 3 ) {
            double dist = iPos.distance(mLC);
            return dist;
        }

        OOFEM_ERROR("Could not compute distance.")
        return 0.0;
    }

private:
    const FloatArray mLC;
    const FloatArray mUC;
    const double mRelTol;

};

} /* namespace oofem */

#endif /* PRESCRIBEDGRADIENTBCWEAK_H_ */
