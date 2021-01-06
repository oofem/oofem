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

#include "prescribedgradientbcweak.h"
#include "classfactory.h"
#include "node.h"
#include "masterdof.h"
#include "element.h"
#include "feinterpol.h"
#include "feinterpol2d.h"
#include "gausspoint.h"
#include "sparsemtrx.h"
#include "xfem/xfemelementinterface.h"
#include "xfem/integrationrules/discsegintegrationrule.h"
#include "spatiallocalizer.h"
#include "geometry.h"
#include "dynamicinputrecord.h"
#include "timestep.h"
#include "function.h"
#include "engngm.h"
#include "mathfem.h"
#include "sparselinsystemnm.h"
#include "unknownnumberingscheme.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/EngineeringModels/staticstructural.h"

#include "timer.h"

#include "xfem/XFEMDebugTools.h"

#include <cmath>
#include <string>
#include <sstream>
#include <iterator>

#ifdef _OPENMP
#include <omp.h>
#endif


namespace oofem {

void TracSegArray::giveTractionLocationArray(IntArray &rows, CharType type, const UnknownNumberingScheme &s)
{
    mFirstNode->giveLocationArray({Trac_u, Trac_v}, rows, s);
}

void TracSegArray :: setupIntegrationRuleOnEl()
{
    mIntRule = std::make_unique<DiscontinuousSegmentIntegrationRule>(1, nullptr, mInteriorSegmentsFine);

    int numPointsPerSeg = 1;
    mIntRule->SetUpPointsOnLine(numPointsPerSeg, _PlaneStrain);
}


PrescribedGradientBCWeak :: PrescribedGradientBCWeak(int n, Domain *d) :
    ActiveBoundaryCondition(n, d),
    PrescribedGradientHomogenization(),
    mTractionDofIDs( {Trac_u, Trac_v} ),
    mDispLockDofIDs( {LMP_u, LMP_v} ),
    mRegularDispDofIDs( {D_u, D_v} ),
    mTractionInterpOrder(0),
    mNumTractionNodesAtIntersections(1),
    mTractionNodeSpacing(1),
    mMeshIsPeriodic(false),
    mDuplicateCornerNodes(false),
    mTangDistPadding(0.0),
    mTracDofScaling(1.0e0),
    mLockNodeInd(0),
    mDispLockScaling(1.0),
    mSpringNodeInd1(-1),
    mSpringNodeInd2(-1),
    mSpringPenaltyStiffness(1.0e-3),
    mPeriodicityNormal({0.0, 1.0}),
    mDomainSize(0.0),
    mMirrorFunction(0)
{
    if ( d ) {
        // Compute bounding box of the domain
        computeDomainBoundingBox(* d, mLC, mUC);
    }
}


PrescribedGradientBCWeak :: ~PrescribedGradientBCWeak()
{
    clear();
}


void PrescribedGradientBCWeak :: clear()
{
    mpTracElNew.clear();
    mpDisplacementLock = nullptr;
}

//#define DAMAGE_TEST

int PrescribedGradientBCWeak :: giveNumberOfInternalDofManagers()
{
    int numDMan = mpTracElNew.size();

    if ( mpDisplacementLock ) {
        numDMan++;
    }

    return numDMan;
}


DofManager *PrescribedGradientBCWeak :: giveInternalDofManager(int i)
{
    if ( i - 1 < int( mpTracElNew.size() ) ) {
        return mpTracElNew [ i - 1 ].mFirstNode.get();
    } else {
        OOFEM_ERROR("return mpDisplacementLock")
        return mpDisplacementLock.get();
    }
}


void PrescribedGradientBCWeak :: initializeFrom(InputRecord &ir)
{
    ActiveBoundaryCondition :: initializeFrom(ir);
    PrescribedGradientHomogenization :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, mTractionInterpOrder, _IFT_PrescribedGradientBCWeak_TractionInterpOrder);
//    printf("mTractionInterpOrder: %d\n", mTractionInterpOrder);

    IR_GIVE_FIELD(ir, mNumTractionNodesAtIntersections, _IFT_PrescribedGradientBCWeak_NumTractionNodesAtIntersections);
//    printf("mNumTractionNodesAtIntersections: %d\n", mNumTractionNodesAtIntersections);

    if ( mNumTractionNodesAtIntersections > 1 && mTractionInterpOrder == 0 ) {
        OOFEM_ERROR("mNumTractionNodesAtIntersections > 1 is not allowed if mTractionInterpOrder == 0.")
    }

    IR_GIVE_FIELD(ir, mTractionNodeSpacing, _IFT_PrescribedGradientBCWeak_NumTractionNodeSpacing);
//    printf("mTractionNodeSpacing: %d\n", mTractionNodeSpacing);

    int duplicateCorners = 0;
    IR_GIVE_FIELD(ir, duplicateCorners, _IFT_PrescribedGradientBCWeak_DuplicateCornerNodes);
//    printf("duplicateCorners: %d\n", duplicateCorners);

    if ( duplicateCorners == 1 ) {
        mDuplicateCornerNodes = true;
    } else {
        mDuplicateCornerNodes = false;
    }

    mTangDistPadding = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, mTangDistPadding, _IFT_PrescribedGradientBCWeak_TangDistPadding);
//    printf("mTangDistPadding: %e\n", mTangDistPadding);

    IR_GIVE_OPTIONAL_FIELD(ir, mTracDofScaling, _IFT_PrescribedGradientBCWeak_TracDofScaling);
//    printf("mTracDofScaling: %e\n", mTracDofScaling );

    IR_GIVE_OPTIONAL_FIELD(ir, mPeriodicityNormal, _IFT_PrescribedGradientBCWeak_PeriodicityNormal);
    mPeriodicityNormal.normalize();
//    printf("mPeriodicityNormal: "); mPeriodicityNormal.printYourself();


    IR_GIVE_OPTIONAL_FIELD(ir, mMirrorFunction, _IFT_PrescribedGradientBCWeak_MirrorFunction);
//    printf("mMirrorFunction: %d\n", mMirrorFunction );

    if ( mMirrorFunction == 0 ) {
        mPeriodicityNormal = {0.0, 1.0};
    }
}


void PrescribedGradientBCWeak :: giveInputRecord(DynamicInputRecord &input)
{
    ActiveBoundaryCondition :: giveInputRecord(input);
    PrescribedGradientHomogenization :: giveInputRecord(input);

    input.setField(mTractionInterpOrder, _IFT_PrescribedGradientBCWeak_TractionInterpOrder);
    input.setField(mNumTractionNodesAtIntersections, _IFT_PrescribedGradientBCWeak_NumTractionNodesAtIntersections);
    input.setField(mTractionNodeSpacing, _IFT_PrescribedGradientBCWeak_NumTractionNodeSpacing);

    if ( mDuplicateCornerNodes ) {
        input.setField(1, _IFT_PrescribedGradientBCWeak_DuplicateCornerNodes);
    } else {
        input.setField(0, _IFT_PrescribedGradientBCWeak_DuplicateCornerNodes);
    }

    input.setField(mTangDistPadding, _IFT_PrescribedGradientBCWeak_TangDistPadding);
    input.setField(mTracDofScaling, _IFT_PrescribedGradientBCWeak_TracDofScaling);

    if ( mMirrorFunction > 0 ) {
        input.setField(mMirrorFunction, _IFT_PrescribedGradientBCWeak_MirrorFunction);
        input.setField(mPeriodicityNormal, _IFT_PrescribedGradientBCWeak_PeriodicityNormal);
    }
}

void PrescribedGradientBCWeak :: postInitialize()
{}

void PrescribedGradientBCWeak :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                                CharType type, ValueModeType mode,
                                                const UnknownNumberingScheme &s, 
                                                FloatArray *eNorm,
                                                void* lock)
{
    int dim = domain->giveNumberOfSpatialDimensions();

    if ( type == ExternalForcesVector ) {
        // The external force vector is given by
        // f_ext = int N^trac H . (x - x_c)


        for ( auto &el : mpTracElNew ) {

            FloatArray contrib;
            computeExtForceElContrib(contrib, el, dim, tStep);

            IntArray rows;
            el.giveTractionLocationArray(rows, type, s);
#ifdef _OPENMP
            if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
            answer.assemble(contrib, rows);
#ifdef _OPENMP
            if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif

        }

    } else if ( type == InternalForcesVector ) {

        for ( auto &el : mpTracElNew ) {

            for ( auto &gp: *el.mIntRule ) {

                // Contribution on gamma_plus
                FloatArray contrib_disp, contrib_trac;
                IntArray disp_loc_array, trac_loc_array;
                computeIntForceGPContrib(contrib_disp, disp_loc_array, contrib_trac, trac_loc_array, el, *gp, dim, tStep, gp->giveGlobalCoordinates(), 1.0, mode, type, s);
#ifdef _OPENMP
                if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
                answer.assemble(contrib_disp, disp_loc_array);
                answer.assemble(contrib_trac, trac_loc_array);
#ifdef _OPENMP
                if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif

                // Contribution on gamma_minus
                contrib_disp.clear(); contrib_trac.clear();
                disp_loc_array.clear(); trac_loc_array.clear();
                FloatArray xMinus;
                this->giveMirroredPointOnGammaMinus(xMinus, gp->giveGlobalCoordinates());
                computeIntForceGPContrib(contrib_disp, disp_loc_array, contrib_trac, trac_loc_array, el, *gp, dim, tStep, xMinus, -1.0, mode, type, s);
#ifdef _OPENMP
                if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
                answer.assemble(contrib_disp, disp_loc_array);
                answer.assemble(contrib_trac, trac_loc_array);
#ifdef _OPENMP
                if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif

            }
    	}

        if ( mpDisplacementLock ) {
            IntArray dispLockRows;
            mpDisplacementLock->giveLocationArray(giveDispLockDofIDs(), dispLockRows, s);

            FloatArray fe_dispLock;

            int lockNodePlaceInArray = domain->giveDofManPlaceInArray(mLockNodeInd);
            FloatArray nodeUnknowns;
            domain->giveDofManager(lockNodePlaceInArray)->giveUnknownVector(nodeUnknowns,this->giveRegularDispDofIDs(), mode, tStep);

            for ( int i = 0; i < domain->giveNumberOfSpatialDimensions(); i++ ) {
                fe_dispLock.push_back(nodeUnknowns [ i ]);
            }

            fe_dispLock.times(mDispLockScaling);
#ifdef _OPENMP
            if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
            answer.assemble(fe_dispLock, dispLockRows);
#ifdef _OPENMP
            if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif
        }

    }
}

void PrescribedGradientBCWeak :: computeExtForceElContrib(FloatArray &oContrib, TracSegArray &iEl, int iDim, TimeStep *tStep)
{

    oContrib.clear();
    FloatArray contrib_gp;

    for ( auto &gp: *iEl.mIntRule ) {

        // Fetch global coordinate x
        const FloatArray &x = gp->giveGlobalCoordinates();

        // Compute H.[x]
        FloatArray temp;
        giveBoundaryCoordVector(temp, x);

        FloatArray Hx;
        FloatMatrix grad2D = mGradient;
        grad2D.resizeWithData(2,2);
        Hx.beProductOf(grad2D, temp);


        // For now, assume piecewise constant approx
        FloatArray Ntrac = FloatArray { 1.0*mTracDofScaling };

        // N-matrix
        FloatMatrix Nmat;
        Nmat.beNMatrixOf(Ntrac, iDim);


        // Assemble contribution to the global vector directly
        contrib_gp.beTProductOf(Nmat, Hx);
        double detJ = 0.5 * iEl.giveLength();
        double loadLevel = this->giveTimeFunction()->evaluateAtTime(tStep->giveTargetTime());
        contrib_gp.times(-detJ * gp->giveWeight()*loadLevel);

        oContrib.add(contrib_gp);
    }

}

void PrescribedGradientBCWeak :: computeIntForceGPContrib(FloatArray &oContrib_disp, IntArray &oDisp_loc_array, FloatArray &oContrib_trac, IntArray &oTrac_loc_array,TracSegArray &iEl, GaussPoint &iGP, int iDim, TimeStep *tStep, const FloatArray &iBndCoord, const double &iScaleFac, ValueModeType mode, CharType type, const UnknownNumberingScheme &s)
{

    SpatialLocalizer *localizer = domain->giveSpatialLocalizer();

    FloatMatrix contrib;
    assembleTangentGPContributionNew(contrib, iEl, iGP, iScaleFac, iBndCoord);

    // Compute vector of traction unknowns
    FloatArray tracUnknowns;
    iEl.mFirstNode->giveUnknownVector(tracUnknowns, giveTracDofIDs(), mode, tStep);

    iEl.giveTractionLocationArray(oTrac_loc_array, type, s);

    FloatArray dispElLocCoord, closestPoint;
    Element *dispEl = localizer->giveElementClosestToPoint(dispElLocCoord, closestPoint, iBndCoord );

    // Compute vector of displacement unknowns
    FloatArray dispUnknowns;
    int numDMan = dispEl->giveNumberOfDofManagers();
    for ( int i = 1; i <= numDMan; i++ ) {
        FloatArray nodeUnknowns;
        DofManager *dMan = dispEl->giveDofManager(i);

        IntArray dispIDs = giveRegularDispDofIDs();
        if ( domain->hasXfemManager() ) {
            XfemManager *xMan = domain->giveXfemManager();
            dispIDs.followedBy(xMan->giveEnrichedDofIDs(*dMan));
        }

        dMan->giveUnknownVector(nodeUnknowns, dispIDs,mode, tStep);
        dispUnknowns.append(nodeUnknowns);

    }

    dispEl->giveLocationArray(oDisp_loc_array, s);


    oContrib_disp.beTProductOf(contrib, tracUnknowns);
    oContrib_disp.negated();

    oContrib_trac.beProductOf(contrib, dispUnknowns);
    oContrib_trac.negated();
}

void PrescribedGradientBCWeak :: assemble( SparseMtrx &answer,
                                        TimeStep *tStep,
                                        CharType type,
                                        const UnknownNumberingScheme &r_s,
                                        const UnknownNumberingScheme &c_s,
                                        double scale,
                                        void* lock)
{
    std::vector<FloatArray> gpCoordArray;

    if ( type == TangentStiffnessMatrix || type == SecantStiffnessMatrix || type == ElasticStiffnessMatrix ) {

        for ( auto &el : mpTracElNew ) {

            for ( auto &gp: *el.mIntRule ) {

                gpCoordArray.push_back( gp->giveGlobalCoordinates() );
                assembleGPContrib(answer, tStep, type, r_s, c_s, el, *gp, scale);
            }
        }

        if ( mpDisplacementLock ) {
            int nsd = domain->giveNumberOfSpatialDimensions();
            FloatMatrix KeDispLock(nsd, nsd);
            KeDispLock.beUnitMatrix();
            KeDispLock.times(mDispLockScaling);

            int placeInArray = domain->giveDofManPlaceInArray(mLockNodeInd);
            DofManager *node = domain->giveDofManager(placeInArray);

            IntArray lockRows, lockCols, nodeRows, nodeCols;
            mpDisplacementLock->giveLocationArray(giveDispLockDofIDs(), lockRows, r_s);
            node->giveLocationArray(giveRegularDispDofIDs(), nodeRows, r_s);

            mpDisplacementLock->giveLocationArray(giveDispLockDofIDs(), lockCols, c_s);
            node->giveLocationArray(giveRegularDispDofIDs(), nodeCols, c_s);
#ifdef _OPENMP
            if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
            answer.assemble(lockRows, nodeCols, KeDispLock);
            answer.assemble(nodeRows, lockCols, KeDispLock);
#ifdef _OPENMP
            if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif

            FloatMatrix KZero( lockRows.giveSize(), lockCols.giveSize() );
            KZero.zero();
            KZero.times(scale);
#ifdef _OPENMP
            if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
            answer.assemble(lockRows, lockCols, KZero);
#ifdef _OPENMP
            if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif
        }


        int nsd = domain->giveNumberOfSpatialDimensions();
        FloatMatrix KeDispLock(nsd, nsd);
        KeDispLock.beUnitMatrix();
        KeDispLock.times(mSpringPenaltyStiffness);

//        printf("mSpringNodeInd1: %d\n", mSpringNodeInd1);
        int placeInArray = domain->giveDofManPlaceInArray(mSpringNodeInd1);
        DofManager *node1 = domain->giveDofManager(placeInArray);

        IntArray nodeRows, nodeCols;
        node1->giveLocationArray(giveRegularDispDofIDs(), nodeRows, r_s);
        node1->giveLocationArray(giveRegularDispDofIDs(), nodeCols, c_s);
        KeDispLock.times(scale);
#ifdef _OPENMP
        if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
        answer.assemble(nodeRows, nodeCols, KeDispLock);
#ifdef _OPENMP
        if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif


    } else {
        printf("Skipping assembly in PrescribedGradientBCWeak::assemble().\n");
    }

//    std :: string fileName("TracGpCoord.vtk");
//    XFEMDebugTools :: WritePointsToVTK(fileName, gpCoordArray);
}

void PrescribedGradientBCWeak :: assembleExtraDisplock(SparseMtrx &answer, TimeStep *tStep,
                       CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)

{
//    printf("Entering PrescribedGradientBCWeak :: assembleExtraDisplock.\n");

#if 1

    int nsd = domain->giveNumberOfSpatialDimensions();
    FloatMatrix KeDispLock(nsd, nsd);
    KeDispLock.zero();

    // Lock in y-direction to get rid of rigid body rotation.
    double scaling = 1.0e0;
    KeDispLock.at(2,2) = mSpringPenaltyStiffness*scaling;

    IntArray nodeRows, nodeCols;

//    printf("mSpringNodeInd2: %d\n", mSpringNodeInd2);
    int placeInArray = domain->giveDofManPlaceInArray(mSpringNodeInd2);
    DofManager *node1 = domain->giveDofManager(placeInArray);

    node1->giveLocationArray(giveRegularDispDofIDs(), nodeRows, r_s);
    node1->giveLocationArray(giveRegularDispDofIDs(), nodeCols, c_s);

    answer.assemble(nodeRows, nodeCols, KeDispLock);


#else
    int nsd = domain->giveNumberOfSpatialDimensions();
    FloatMatrix KeDispLock(2*nsd, 2*nsd);
    KeDispLock.zero();

    // Lock in y-direction to get rid of rigid body rotation.
    double scaling = 1.0e0;
    KeDispLock.at(2,2)         =  mSpringPenaltyStiffness*scaling;
    KeDispLock.at(2,nsd+1)     = -mSpringPenaltyStiffness*scaling;
    KeDispLock.at(nsd+1,2)     = -mSpringPenaltyStiffness*scaling;
    KeDispLock.at(nsd+1,nsd+1) =  mSpringPenaltyStiffness*scaling;
//    KeDispLock.times(mSpringPenaltyStiffness);

    IntArray nodeRows, nodeCols;


    int placeInArray = domain->giveDofManPlaceInArray(mSpringNodeInd2);
    DofManager *node1 = domain->giveDofManager(placeInArray);

    IntArray nodeRows1, nodeCols1;
    node1->giveLocationArray(giveRegularDispDofIDs(), nodeRows1, r_s);
    node1->giveLocationArray(giveRegularDispDofIDs(), nodeCols1, c_s);

    nodeRows.followedBy(nodeRows1);
    nodeCols.followedBy(nodeCols1);


    placeInArray = domain->giveDofManPlaceInArray(mSpringNodeInd3);
    DofManager *node2 = domain->giveDofManager(placeInArray);

    IntArray nodeRows2, nodeCols2;
    node2->giveLocationArray(giveRegularDispDofIDs(), nodeRows2, r_s);
    node2->giveLocationArray(giveRegularDispDofIDs(), nodeCols2, c_s);

    nodeRows.followedBy(nodeRows2);
    nodeCols.followedBy(nodeCols2);

    answer.assemble(nodeRows, nodeCols, KeDispLock);
#endif

}

void PrescribedGradientBCWeak :: assembleGPContrib(SparseMtrx &answer, TimeStep *tStep,
                                                    CharType type, const UnknownNumberingScheme &r_s, 
                                                    const UnknownNumberingScheme &c_s, TracSegArray &iEl, 
                                                    GaussPoint &iGP, double k, void* lock)
{

    SpatialLocalizer *localizer = domain->giveSpatialLocalizer();

    ///////////////
    // Gamma_plus
    FloatMatrix contrib;
    assembleTangentGPContributionNew(contrib, iEl, iGP, -1.0, iGP.giveGlobalCoordinates());

    // Compute vector of traction unknowns
    FloatArray tracUnknowns;
    iEl.mFirstNode->giveUnknownVector(tracUnknowns, giveTracDofIDs(), VM_Total, tStep);

    IntArray trac_rows;
    iEl.giveTractionLocationArray(trac_rows, type, r_s);


    FloatArray dispElLocCoord, closestPoint;
    Element *dispEl = localizer->giveElementClosestToPoint(dispElLocCoord, closestPoint, iGP.giveGlobalCoordinates() );

    IntArray disp_cols;
    dispEl->giveLocationArray(disp_cols, c_s);

    contrib.times(k);
#ifdef _OPENMP
    if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
    answer.assemble(trac_rows, disp_cols, contrib);
#ifdef _OPENMP
    if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif

    FloatMatrix contribT;
    contribT.beTranspositionOf(contrib);
#ifdef _OPENMP
    if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
    answer.assemble(disp_cols, trac_rows, contribT);
#ifdef _OPENMP
    if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif
    ///////////////
    // Gamma_minus
    contrib.clear();
    FloatArray xMinus;
    this->giveMirroredPointOnGammaMinus(xMinus, iGP.giveGlobalCoordinates() );
    assembleTangentGPContributionNew(contrib, iEl, iGP, 1.0, xMinus);

    // Compute vector of traction unknowns
    tracUnknowns.clear();
    iEl.mFirstNode->giveUnknownVector(tracUnknowns, giveTracDofIDs(), VM_Total, tStep);

    trac_rows.clear();
    iEl.giveTractionLocationArray(trac_rows, type, r_s);


    dispElLocCoord.clear(); closestPoint.clear();
    dispEl = localizer->giveElementClosestToPoint(dispElLocCoord, closestPoint, xMinus );

    disp_cols.clear();
    dispEl->giveLocationArray(disp_cols, c_s);
    contrib.times(k);
#ifdef _OPENMP
    if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
    answer.assemble(trac_rows, disp_cols, contrib);
#ifdef _OPENMP
    if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif
    contribT.clear();
    contribT.beTranspositionOf(contrib);
#ifdef _OPENMP
    if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
    answer.assemble(disp_cols, trac_rows, contribT);
#ifdef _OPENMP
    if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif

    // Assemble zeros on diagonal (required by PETSc solver)
    FloatMatrix KZero(1,1);
    KZero.zero();
#ifdef _OPENMP
    if (lock) omp_set_lock(static_cast<omp_lock_t*>(lock));
#endif
    for ( int i :  trac_rows) {
        answer.assemble(IntArray({i}), IntArray({i}), KZero);
    }
#ifdef _OPENMP
    if (lock) omp_unset_lock(static_cast<omp_lock_t*>(lock));
#endif
}

void PrescribedGradientBCWeak :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                                                    const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{

}

void PrescribedGradientBCWeak :: giveTractionLocationArray(IntArray &rows,
                                                           const UnknownNumberingScheme &s)
{
    // Used for the condensation when computing the macroscopic tangent.

    rows.clear();

    // Loop over traction elements
    for ( auto &el : mpTracElNew ) {

        IntArray trac_loc_r;
        el.mFirstNode->giveLocationArray(giveTracDofIDs(), trac_loc_r, s);

        rows.followedBy(trac_loc_r);
    }

#if 0
    rows.clear();

    // Loop over traction elements
    for ( size_t tracElInd = 0; tracElInd < mpTractionElements.size(); tracElInd++ ) {
        IntArray tracElRows, trac_loc_r;

        const TractionElement &tEl = mpTractionElements [ tracElInd ];
        for ( int tracNodeInd : tEl.mTractionNodeInd ) {
            Node *tNode = mpTractionNodes [ tracNodeInd ];
            tNode->giveLocationArray(giveTracDofIDs(), trac_loc_r, s);
            tracElRows.followedBy(trac_loc_r);
        }

        rows.followedBy(tracElRows);
    }

    if ( mpDisplacementLock ) {
        IntArray dispLock_r;
        mpDisplacementLock->giveLocationArray(giveDispLockDofIDs(), dispLock_r, s);

        rows.followedBy(dispLock_r);
    }
#endif
}

void PrescribedGradientBCWeak :: giveDisplacementLocationArray(IntArray &rows, const UnknownNumberingScheme &s)
{
    // Used for the condensation when computing the macroscopic tangent.
}

void PrescribedGradientBCWeak :: compute_x_times_N_1(FloatMatrix &o_x_times_N)
{
    //double rve_size = this->domainSize();
    const int dim = domain->giveNumberOfSpatialDimensions();

    IntArray loc_t;
    EModelDefaultEquationNumbering fnum;
    giveTractionLocationArray(loc_t, fnum);

    int num_t_eq = loc_t.giveSize();

    //o_x_times_N.resize(3, num_t_eq);
    o_x_times_N.resize(num_t_eq,3);

    IntArray cols = {1,2,3};

    int trac_el_ind = 1;

    for ( auto &el : mpTracElNew ) {

        IntArray rows = {trac_el_ind*2-1,trac_el_ind*2};

        for ( auto &gp : *el.mIntRule ) {

            FloatMatrix contrib(2,3);

            // For now, assume piecewise constant approx
            FloatArray Ntrac = FloatArray { 1.0*mTracDofScaling };

            // N-matrix
            FloatMatrix Nmat;
            Nmat.beNMatrixOf(Ntrac, dim);
//            printf("Nmat: "); Nmat.printYourself();

            // Fetch global coordinate x
            const FloatArray &x = gp->giveGlobalCoordinates();

//            // Compute vector of traction unknowns
//            FloatArray tracUnknowns;
//            el.mFirstNode->giveUnknownVector(tracUnknowns, giveTracDofIDs(), VM_Total, tStep);


//            FloatArray traction;
//            traction.beProductOf(Nmat, tracUnknowns);

            FloatArray tmp;
            giveBoundaryCoordVector(tmp, x);

            FloatMatrix coord_mat(2,3);
            coord_mat.at(1,1) = tmp.at(1);
            coord_mat.at(2,2) = tmp.at(2);

//            coord_mat.at(3,2) = tmp.at(1);

//            coord_mat.at(3,1) = tmp.at(2);
//            coord_mat.at(4,2) = tmp.at(1);

//            coord_mat.at(4,1) = tmp.at(2);
//            coord_mat.at(3,2) = tmp.at(1);

            coord_mat.at(1,3) = 0.5*tmp.at(2);
            coord_mat.at(2,3) = 0.5*tmp.at(1);

//            coord_mat.at(4,2) = 0.5*tmp.at(2);
//            coord_mat.at(4,1) = 0.5*tmp.at(1);


//            FloatMatrix contrib;
//            contrib.beDyadicProductOf(traction, tmp);

            contrib.beProductOf(Nmat, coord_mat);

//            contrib = coord_mat;

            double detJ = 0.5 * el.giveLength();
            contrib.times( detJ * gp->giveWeight() );

//            printf("\n\ncontrib: "); contrib.printYourself();
//            printf("rows: "); rows.printYourself();
//            printf("cols: "); cols.printYourself();

            o_x_times_N.assemble(contrib, rows, cols);
        }

        trac_el_ind++;
    }
}

void PrescribedGradientBCWeak :: compute_x_times_N_2(FloatMatrix &o_x_times_N)
{
    //double rve_size = this->domainSize();
    const int dim = domain->giveNumberOfSpatialDimensions();

    IntArray loc_t;
    EModelDefaultEquationNumbering fnum;
    giveTractionLocationArray(loc_t, fnum);

    int num_t_eq = loc_t.giveSize();

    //o_x_times_N.resize(4, num_t_eq);
    o_x_times_N.resize(3, num_t_eq);

    IntArray rows = {1,2,3};

    int trac_el_ind = 1;

    for ( auto &el : mpTracElNew ) {

    	IntArray cols = {trac_el_ind*2-1,trac_el_ind*2};

        for ( auto &gp : *el.mIntRule ) {

        	FloatMatrix contrib(4,2);

            // For now, assume piecewise constant approx
            FloatArray Ntrac = FloatArray { 1.0*mTracDofScaling };

            // N-matrix
            FloatMatrix Nmat;
            Nmat.beNMatrixOf(Ntrac, dim);
//            printf("Nmat: "); Nmat.printYourself();

            // Fetch global coordinate x
            const FloatArray &x = gp->giveGlobalCoordinates();

//            // Compute vector of traction unknowns
//            FloatArray tracUnknowns;
//            el.mFirstNode->giveUnknownVector(tracUnknowns, giveTracDofIDs(), VM_Total, tStep);


//            FloatArray traction;
//            traction.beProductOf(Nmat, tracUnknowns);

            FloatArray tmp;
            giveBoundaryCoordVector(tmp, x);

            FloatMatrix coord_mat(4,2);
            coord_mat.at(1,1) = tmp.at(1);
            coord_mat.at(2,2) = tmp.at(2);

//            coord_mat.at(3,2) = tmp.at(1);
//            coord_mat.at(3,1) = tmp.at(2);

//            coord_mat.at(3,1) = tmp.at(2);
//            coord_mat.at(4,2) = tmp.at(1);

            coord_mat.at(3,1) = 0.5*tmp.at(2);
            coord_mat.at(3,2) = 0.5*tmp.at(1);
//
//            coord_mat.at(4,2) = 0.5*tmp.at(2);
//            coord_mat.at(4,1) = 0.5*tmp.at(1);


//            FloatMatrix contrib;
//            contrib.beDyadicProductOf(traction, tmp);

            contrib.beProductOf(coord_mat, Nmat);

//            contrib = coord_mat;

            double detJ = 0.5 * el.giveLength();
            contrib.times( detJ * gp->giveWeight() );

//            printf("\n\ncontrib: "); contrib.printYourself();
//            printf("rows: "); rows.printYourself();
//            printf("cols: "); cols.printYourself();

            o_x_times_N.assemble(contrib, rows, cols);
        }

        trac_el_ind++;
    }
}


void PrescribedGradientBCWeak :: computeField(FloatArray &sigma, TimeStep *tStep)
{
    double Lx = mUC[0] - mLC[0];
    double Ly = mUC[1] - mLC[1];
    double dSize = Lx*Ly;
//    printf("dSize: %e\n", dSize);

    const int dim = domain->giveNumberOfSpatialDimensions();
    FloatMatrix stressMatrix(dim, dim);

    for ( auto &el : mpTracElNew ) {

        for ( auto &gp : *el.mIntRule ) {

            // For now, assume piecewise constant approx
            FloatArray Ntrac = FloatArray { 1.0*mTracDofScaling };

            // N-matrix
            FloatMatrix Nmat;
            Nmat.beNMatrixOf(Ntrac, dim);

            // Fetch global coordinate x
            const FloatArray &x = gp->giveGlobalCoordinates();

            // Compute vector of traction unknowns
            FloatArray tracUnknowns;
            el.mFirstNode->giveUnknownVector(tracUnknowns, giveTracDofIDs(), VM_Total, tStep);


            FloatArray traction;
            traction.beProductOf(Nmat, tracUnknowns);

            FloatArray tmp;
            giveBoundaryCoordVector(tmp, x);

            FloatMatrix contrib;
            contrib.beDyadicProductOf(traction, tmp);

            double detJ = 0.5 * el.giveLength();
            contrib.times( detJ * gp->giveWeight() );

            for ( int m = 0; m < dim; m++ ) {
                for ( int n = 0; n < dim; n++ ) {
                    stressMatrix(m, n) += contrib(m, n);
                }
            }

        }

    }

    if ( dim == 2 ) {
        sigma = {
            stressMatrix(0, 0), stressMatrix(1, 1), 0.0, 0.0, 0.0, 0.5*(stressMatrix(0, 1) + stressMatrix(1, 0))
        };
    } else {
        sigma.beVectorForm(stressMatrix);
    }

    sigma.times(1.0 / dSize);

}

//#define TIME_INFO

void PrescribedGradientBCWeak :: computeTangent(FloatMatrix& E, TimeStep* tStep)
{
#ifdef TIME_INFO
    static double tot_time = 0.0;
    Timer timer;
    timer.startTimer();

    static double assemble_time = 0.0;
    Timer assemble_timer;

#endif

    E.resize(9,9);


    // Extract the relevant submatrices from the RVE problem.
    // At equilibrium, the RVE problem has the structure
    //
    // [ S C; C^T 0 ]*[a_u a_t] = [0; f]
    //
    // where a_u and a_t denote displacement and traction dofs, respectively.
    // We need to extract S and C.

    EngngModel *rve = this->giveDomain()->giveEngngModel();
    ///@todo Get this from engineering model
    std :: unique_ptr< SparseLinearSystemNM > solver(
        classFactory.createSparseLinSolver( ST_Petsc, this->domain, this->domain->giveEngngModel() ) ); // = rve->giveLinearSolver();
    bool symmetric_matrix = false;
    SparseMtrxType stype = solver->giveRecommendedMatrix(symmetric_matrix);
//    double rve_size = this->domainSize();
    double Lx = mUC[0] - mLC[0];
    double Ly = mUC[1] - mLC[1];
    double rve_size = Lx*Ly;


    EModelDefaultEquationNumbering fnum;
    std :: unique_ptr< SparseMtrx > Kmicro( classFactory.createSparseMtrx(stype) );
    if ( !Kmicro ) {
        OOFEM_ERROR("Couldn't create sparse matrix of type %d\n", stype);
    }

#ifdef __SM_MODULE
    StaticStructural *rveStatStruct = dynamic_cast<StaticStructural*>(rve);
    if ( rveStatStruct ) {
        //printf("Successfully casted rve to StaticStructural.\n");

        if ( rveStatStruct->stiffnessMatrix ) {
            Kmicro = rveStatStruct->stiffnessMatrix->clone();
        }
    }
#endif


#ifdef TIME_INFO
    assemble_timer.startTimer();
#endif

    if ( Kmicro->giveNumberOfColumns() == 0 ) {
        //printf("Rebuilding stiffness matrix.\n");
        Kmicro->buildInternalStructure(rve, this->domain->giveNumber(), fnum);
        rve->assemble(*Kmicro, tStep, TangentAssembler(TangentStiffness), fnum, this->domain);
    }
//    else {
//        printf("Using existing stiffness matrix.\n");
//    }

    assembleExtraDisplock(*Kmicro, tStep, TangentStiffnessMatrix, fnum, fnum);
#ifdef TIME_INFO
    assemble_timer.stopTimer();
    assemble_time += assemble_timer.getUtime();
    printf("Assembly time for RVE tangent: %e\n", assemble_time);
#endif

    // Fetch displacement and traction location arrays
    IntArray loc_u, loc_t;

    giveTractionLocationArray(loc_t, fnum);

    int neq = Kmicro->giveNumberOfRows();
    loc_u.resize(neq - loc_t.giveSize());
    int k = 1;
    for ( int i = 1; i <= neq; i++ ) {
        if ( !loc_t.contains(i) ) {
            loc_u.at(k) = i;
            k++;
        }
    }

    // Fetch the submatrices
    std :: unique_ptr< SparseMtrx > S = Kmicro->giveSubMatrix(loc_u, loc_u);
    // NOTE: Kus is actually a dense matrix, but we have to make it a dense matrix first
    std :: unique_ptr< SparseMtrx > C = Kmicro->giveSubMatrix(loc_u, loc_t);
    FloatMatrix Cd;
    C->toFloatMatrix(Cd);

    // Let Sm = C. Solve for m.
    FloatMatrix m;
    solver->solve(*S, Cd, m);

    // Compute G := C^T m. (Which could, formally, be written as G = C^T S^-1 C.)
    FloatMatrix G;
    G.beTProductOf(Cd, m);

    // Compute D := \int x \otimes N
    FloatMatrix D;
    compute_x_times_N_1(D);
//    FloatMatrix D;
//    compute_x_times_N_1(DT2);
//    D.beTranspositionOf(DT);

    // Let Gp = D. Solve for p. (Which could, formally, be written as p = G^-1 D = (C^T S^-1 C)^-1 D.)
    // Need to make G sparse;
    std :: unique_ptr< SparseMtrx > Gs( classFactory.createSparseMtrx(stype) );
    if ( !Gs ) {
        OOFEM_ERROR("Couldn't create sparse matrix of type %d\n", stype);
    }

    int num_eq_G = G.giveNumberOfRows();
    IntArray loc_G;
    loc_G.enumerate(num_eq_G);

    Gs->buildInternalStructure(rve, num_eq_G, num_eq_G, loc_G, loc_G);
    Gs->assemble(loc_G, loc_G, G);

    Gs->assembleBegin();
    Gs->assembleEnd();

//    Gs->writeToFile("Gs.txt");

    FloatMatrix p;
    solver->solve(*Gs, D, p);

    // Compute d_sigma_depsilon = D^T p.
    FloatMatrix Ered;
    Ered.beTProductOf(D,p);
//    rve_size = 25.0;
    Ered.times(1.0/rve_size);
//    printf("rve_size: %e\n", rve_size);
//    Ered.printYourself();

    IntArray indx = {1,2,6,9}; // to make it independednt of sm module
    //StructuralMaterial :: giveVoigtVectorMask(indx, _PlaneStress);
    

//    FloatMatrix EredT;
//    EredT.beTranspositionOf(Ered);
//    Ered.add(EredT);
//    Ered.times(0.5);

//    Ered.at(1,3) = 0.0;
//    Ered.at(2,3) = 0.0;
//    Ered.at(3,1) = 0.0;
//    Ered.at(3,2) = 0.0;

    E.assemble(Ered, indx, indx);
//    E.printYourself();

#ifdef TIME_INFO
    timer.stopTimer();
    tot_time += timer.getUtime();
//    printf("Total time for RVE tangent: %e\n", tot_time);
#endif

}

void PrescribedGradientBCWeak :: giveTractionElNormal(size_t iElInd, FloatArray &oNormal, FloatArray &oTangent) const
{
    FloatArray xS, xE;
    giveTractionElCoord(iElInd, xS, xE);

    oTangent.beDifferenceOf(xE, xS);
    oTangent.normalize();

    oNormal = {
        oTangent [ 1 ], -oTangent [ 0 ]
    };
}

void PrescribedGradientBCWeak :: giveTractionElArcPos(size_t iElInd, double &oXiStart, double &oXiEnd) const
{
    FloatArray xS, xE;
    giveTractionElCoord(iElInd, xS, xE);

    FloatArray xC;
    xC.beScaled(0.5, xS);
    xC.add(0.5, xE);
    int sideIndex = giveSideIndex(xC);

    const double nodeDistTol = 1.0e-15;
    ArcPosSortFunction3< bool >sortFunc(mLC, mUC, nodeDistTol, sideIndex);

    oXiStart = sortFunc.calcArcPos(xS);
    oXiEnd = sortFunc.calcArcPos(xE);
}

void PrescribedGradientBCWeak :: giveBoundaries(IntArray &oBoundaries)
{
    Set *setPointer = this->giveDomain()->giveSet(this->set);
    oBoundaries = setPointer->giveBoundaryList();
}

void PrescribedGradientBCWeak :: giveTraction(size_t iElInd, FloatArray &oStartTraction, FloatArray &oEndTraction, ValueModeType mode, TimeStep *tStep)
{
    // For now, assuming piecewise constant traction
    mpTracElNew[iElInd].mFirstNode->giveUnknownVector(oStartTraction, giveTracDofIDs(), mode, tStep);
    oStartTraction.times(mTracDofScaling);
    mpTracElNew[iElInd].mFirstNode->giveUnknownVector(oEndTraction, giveTracDofIDs(), mode, tStep);
    oEndTraction.times(mTracDofScaling);
}

void PrescribedGradientBCWeak :: recomputeTractionMesh()
{
    clear();
    postInitialize();
}

void PrescribedGradientBCWeak :: createTractionMesh(bool iEnforceCornerPeriodicity, int iNumSides)
{
    bool split_at_holes = true;

    const double l_s = mUC[0] - mLC[0];
    const double minPointDist = 1.0e-4*l_s;

    // Find holes intersecting the RVE boundary so that these can be excluded
    std::vector<FloatArray> holeCoordUnsorted, allCoordUnsorted;
    findHoleCoord(holeCoordUnsorted, allCoordUnsorted);

    // Add corner points
    holeCoordUnsorted.push_back( {mUC[0], mLC[1]} );
    allCoordUnsorted.push_back( {mUC[0], mLC[1]} );

    holeCoordUnsorted.push_back( {mUC[0], mUC[1]} );
    allCoordUnsorted.push_back( {mUC[0], mUC[1]} );

    holeCoordUnsorted.push_back( {mLC[0], mUC[1]} );
    allCoordUnsorted.push_back( {mLC[0], mUC[1]} );


    // Add crack-boundary intersections
    findCrackBndIntersecCoord(holeCoordUnsorted);

    // Add periodicity points
    findPeriodicityCoord(holeCoordUnsorted);


    // Sort arrays in terms of arc length along the RVE boundary
    std :: sort( holeCoordUnsorted.begin(), holeCoordUnsorted.end(), ArcPosSortFunction4( mLC, mUC, 1.0e-4 ) );
    std :: sort( allCoordUnsorted.begin(), allCoordUnsorted.end(), ArcPosSortFunction4( mLC, mUC, 1.0e-4 ) );


    // Add refinement points
    int pointsPassed = 0;
    for ( const auto &x : allCoordUnsorted ) {

        if ( pointsPassed >= mTractionNodeSpacing ) {
            holeCoordUnsorted.push_back(x);
            pointsPassed = 0;
        }

        pointsPassed++;
    }

    // Sort again
    std :: sort( holeCoordUnsorted.begin(), holeCoordUnsorted.end(), ArcPosSortFunction4( mLC, mUC, 1.0e-4 ) );
    std :: sort( allCoordUnsorted.begin(), allCoordUnsorted.end(), ArcPosSortFunction4( mLC, mUC, 1.0e-4 ) );

    // Remove points that are too close to each other
    removeClosePoints(holeCoordUnsorted, minPointDist);
    removeClosePoints(allCoordUnsorted, minPointDist);

    // Create two arrays of segments, where each array represents the coarsest possible traction
    // mesh on one side of the RVE
    ArcPosSortFunction4 arcPosFunc( mLC, mUC, 1.0e-4 );

    std :: vector< TracSegArray > tracElNew0, tracElNew1;
    tracElNew0.emplace_back();
    //tracElNew1.emplace_back();

    for (size_t i = 1; i < holeCoordUnsorted.size(); i++) {

        FloatArray xS = holeCoordUnsorted[i-1];
        xS.resizeWithValues(2);
        FloatArray xE = holeCoordUnsorted[i];
        xE.resizeWithValues(2);
        const FloatArray xC = {0.5*(xS[0]+xE[0]), 0.5*(xS[1]+xE[1])};

        if ( arcPosFunc.calcArcPos(xC) < 2.*l_s ) {
            tracElNew0[0].mInteriorSegments.emplace_back(xS, xE);
        } else {
            tracElNew1[0].mInteriorSegments.emplace_back(xS, xE);
        }
    }

    // Remove segments located in holes
    removeSegOverHoles(tracElNew0[0], 1.0e-4);
    removeSegOverHoles(tracElNew1[0], 1.0e-4);

    if ( split_at_holes ) {
        splitSegments(tracElNew0);
        splitSegments(tracElNew1);
    }

    // Identify additional points that can be used to refine the traction mesh


    //////////////////////////////////////////////////
    // Create traction dofs
    int numNodes = domain->giveNumberOfDofManagers();
    int totNodesCreated = 0;


    // For each side (0 and 1), loop over over elements

    // We may always create the first node on the element.
    // For the linear approximation, it may need to be a slave node,
    // depending on which element it is.

    // For now, consider only piecewise constant approximations. Then,
    // we can always create on node on each element.

    // RVE side at x=L
    for ( auto &el : tracElNew0 ) {

        //////////////////////////////////////////////////////
        // Create first node
        totNodesCreated++;

        el.mFirstNode = std::make_unique<Node>(numNodes + 1, domain);
        el.mFirstNode->setGlobalNumber(numNodes + 1);
        for ( auto &dofId: giveTracDofIDs() ) {
        	el.mFirstNode->appendDof( new MasterDof(el.mFirstNode.get(), ( DofIDItem ) dofId) );
        }

        el.mFirstNode->setCoordinates( el.mInteriorSegments[0].giveVertex(1) );
        numNodes++;
    }


    // RVE side at y=L
    for ( auto &el : tracElNew1 ) {

        //////////////////////////////////////////////////////
        // Create first node
        totNodesCreated++;

        el.mFirstNode = std::make_unique<Node>(numNodes + 1, domain);
        el.mFirstNode->setGlobalNumber(numNodes + 1);
        for ( auto &dofId: giveTracDofIDs() ) {
        	el.mFirstNode->appendDof( new MasterDof(el.mFirstNode.get(), ( DofIDItem ) dofId) );
        }

        el.mFirstNode->setCoordinates( el.mInteriorSegments[0].giveVertex(1) );
        numNodes++;
    }

    if ( mMeshIsPeriodic && false ) {
        // Lock displacement in one node if we use periodic BCs

        int numNodes = domain->giveNumberOfDofManagers();
        mpDisplacementLock = std::make_unique<Node>(numNodes + 1, domain);
        mLockNodeInd = domain->giveElement(1)->giveNode(1)->giveGlobalNumber();


        for ( auto &dofid: giveDispLockDofIDs() ) {
            mpDisplacementLock->appendDof( new MasterDof(mpDisplacementLock.get(), ( DofIDItem ) dofid) );
        }
    }


    // Nodes to lock in order to prevent rigid body motion
    SpatialLocalizer *localizer = domain->giveSpatialLocalizer();
    FloatArray x1 = mLC;
//    printf("x1: "); x1.printYourself();
    double maxDist = 1.0e10;
    Node *node1 = localizer->giveNodeClosestToPoint(x1, maxDist);
    mSpringNodeInd1 = node1->giveGlobalNumber();
//    printf("mSpringNodeInd1: %d\n", mSpringNodeInd1 );

    FloatArray x2 = {mUC.at(1), mLC.at(2)};
//    FloatArray x2 = {mUC.at(1), mUC.at(2)};
//    printf("x2: "); x2.printYourself();
    Node *node2 = localizer->giveNodeClosestToPoint(x2, maxDist);
    mSpringNodeInd2 = node2->giveGlobalNumber();
//    printf("mSpringNodeInd2: %d\n", mSpringNodeInd2 );


    FloatArray x3 = {mLC.at(1), mUC.at(2)};
//    FloatArray x2 = {mUC.at(1), mUC.at(2)};
//    printf("x3: "); x3.printYourself();
    Node *node3 = localizer->giveNodeClosestToPoint(x3, maxDist);
    mSpringNodeInd3 = node3->giveGlobalNumber();
//    printf("mSpringNodeInd3: %d\n", mSpringNodeInd3 );


    mpTracElNew.reserve(mpTracElNew.size() + tracElNew0.size() + tracElNew1.size());
    std::move(tracElNew0.begin(), tracElNew0.end(), std::back_inserter(mpTracElNew));
    std::move(tracElNew1.begin(), tracElNew1.end(), std::back_inserter(mpTracElNew));
    tracElNew0.clear();
    tracElNew1.clear();


    ////////////
    // Segment arrays for Gauss quadrature
    size_t i = 0;

    for ( auto & el : mpTracElNew ) {

        const FloatArray &xS = el.mInteriorSegments[0].giveVertex(1);
        const double arcPosXS = arcPosFunc.calcArcPos(xS);

        const FloatArray &xE = el.mInteriorSegments.back().giveVertex(2);
        const double arcPosXE = arcPosFunc.calcArcPos(xE);

        while (i < allCoordUnsorted.size()) {

            FloatArray x = allCoordUnsorted[i];
            x.resizeWithValues(2);
            const double arcPosX = arcPosFunc.calcArcPos(x);

            if ( arcPosX > (arcPosXS+minPointDist) && arcPosX < (arcPosXE-minPointDist) ) {
                el.mInteriorSegmentsPointsFine.push_back(std::move(x));
            }

            if ( arcPosX > arcPosXE ) {
                break;
            }

            i++;
        }
    }

    // Now we have the necessary points on each traction element.
    // The next step is to create splitted segments.
    for ( auto & el : mpTracElNew ) {

        i = 0;

        for ( auto &line : el.mInteriorSegments ) {
            FloatArray xS = line.giveVertex(1);
            xS.resizeWithValues(2);
            const double arcPosXS = arcPosFunc.calcArcPos(xS);

            FloatArray xE = line.giveVertex(2);
            xE.resizeWithValues(2);
            const double arcPosXE = arcPosFunc.calcArcPos(xE);

            if ( el.mInteriorSegmentsPointsFine.size() == 0 ) {
                Line newLine(xS, xE);
                el.mInteriorSegmentsFine.push_back(newLine);
            } else {
                while ( i < el.mInteriorSegmentsPointsFine.size() ) {

                    const FloatArray &x = el.mInteriorSegmentsPointsFine[i];
                    const double arcPosX = arcPosFunc.calcArcPos(x);

                    if ( arcPosX < arcPosXS ) {
                        OOFEM_ERROR("Error in PrescribedGradientBCWeak :: createTractionMesh.")
                    }

                    if ( arcPosX < arcPosXE ) {
                        // Split from start pos to x
                        Line newLine(xS, x);
                        el.mInteriorSegmentsFine.push_back(newLine);

                        xS = x;
                    } else {
                        // Split from x to end pos
                        Line newLine(xS, xE);
                        el.mInteriorSegmentsFine.push_back(newLine);

                        break;
                    }

                    if ( i == (el.mInteriorSegmentsPointsFine.size()-1) ) {
                        // Split from x to end pos
                        Line newLine(xS, xE);
                        el.mInteriorSegmentsFine.push_back(newLine);
                    }

                    i++;
                }
            }
        }
    }

    // Create integration rules
    for ( auto & el : mpTracElNew ) {
        el.setupIntegrationRuleOnEl();
    }

    // Write discontinuity points to debug vtk
    std :: vector< FloatArray > discPoints;
    for ( auto & el : mpTracElNew ) {

        discPoints.push_back( el.mInteriorSegments[0].giveVertex(1) );
        discPoints.push_back( el.mInteriorSegments.back().giveVertex(2) );
    }

//    std :: string fileName("DiscontPoints.vtk");
//    XFEMDebugTools :: WritePointsToVTK(fileName, discPoints);

}

void PrescribedGradientBCWeak :: splitSegments(std :: vector< TracSegArray > &ioElArray)
{
    std :: vector< TracSegArray > newArray;

    for ( auto &el : ioElArray ) {
        for ( auto &line : el.mInteriorSegments ) {
            TracSegArray newEl;
            newEl.mInteriorSegments.push_back(line);
            newArray.push_back(std::move(newEl));
        }
    }

    ioElArray = std::move(newArray);
}

bool PrescribedGradientBCWeak :: damageExceedsTolerance(Element *el)
{
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();

    double maxDamage = 0.0, damageTol = 0.01;
    for ( auto &gp: *el->giveDefaultIntegrationRulePtr() ) {
        FloatArray damage;
        el->giveIPValue(damage, gp, IST_DamageScalar, tStep);
//                    printf("damage: "); damage.printYourself();
        maxDamage = std::max(maxDamage, damage(0));
    }

//    printf("maxDamage: %e\n", maxDamage);

    return maxDamage >= damageTol;
}


void PrescribedGradientBCWeak :: assembleTangentGPContributionNew(FloatMatrix &oTangent, TracSegArray &iEl, GaussPoint &iGP, const double &iScaleFactor, const FloatArray &iBndCoord)
{
    int dim = domain->giveNumberOfSpatialDimensions();
    double detJ = 0.5 * iEl.giveLength();

    //////////////////////////////////
    // Compute traction N-matrix
    // For now, assume piecewise constant approx
    FloatArray Ntrac = FloatArray { 1.0*mTracDofScaling };
    FloatMatrix NtracMat;
    NtracMat.beNMatrixOf(Ntrac, dim);

    //////////////////////////////////
    // Compute displacement N-matrix
    // Identify the displacement element
    // we are currently standing in
    // and compute local coordinates on
    // the displacement element
    SpatialLocalizer *localizer = domain->giveSpatialLocalizer();
    FloatArray dispElLocCoord, closestPoint;
    Element *dispEl = localizer->giveElementClosestToPoint(dispElLocCoord, closestPoint, iBndCoord );

    // Compute basis functions
    XfemElementInterface *xfemElInt = dynamic_cast< XfemElementInterface * >( dispEl );
    FloatMatrix NdispMat;

    if ( xfemElInt && domain->hasXfemManager() ) {
        // If the element is an XFEM element, we use the XfemElementInterface to compute the N-matrix
        // of the enriched element.
        xfemElInt->XfemElementInterface_createEnrNmatrixAt(NdispMat, dispElLocCoord, * dispEl, false);
    } else {
        // Otherwise, use the usual N-matrix.
        FloatArray N;

        int dim = dispEl->giveSpatialDimension();

        dispEl->giveInterpolation()->evalN( N, dispElLocCoord, FEIElementGeometryWrapper(dispEl) );

        NdispMat.beNMatrixOf(N, dim);
    }

    oTangent.beTProductOf(NtracMat, NdispMat);
    oTangent.times( iScaleFactor * detJ * iGP.giveWeight() );
}

bool PrescribedGradientBCWeak :: pointIsOnGammaPlus(const FloatArray &iPos) const
{
    const double distTol = 1.0e-12;

    if ( iPos [ 0 ] > mUC [ 0 ] - distTol ) {
        return true;
    }

    if ( iPos [ 1 ] > mUC [ 1 ] - distTol ) {
        return true;
    }

    return false;
}

void PrescribedGradientBCWeak :: giveMirroredPointOnGammaMinus(FloatArray &oPosMinus, const FloatArray &iPosPlus) const
{
//#if 0
    if ( mMirrorFunction == 0 ) {
        oPosMinus = iPosPlus;
        const double distTol = 1.0e-12;

//      if ( distance(iPosPlus, mUC) < distTol ) {
//          printf("iPosPlus: %.12e %.12e\n", iPosPlus [ 0 ], iPosPlus [ 1 ]);
//          OOFEM_ERROR("Unmappable point.")
//      }

        if ( iPosPlus [ 0 ] > mUC [ 0 ] - distTol ) {
            oPosMinus [ 0 ] = mLC [ 0 ];
            return;
        } else if ( iPosPlus [ 1 ] > mUC [ 1 ] - distTol ) {
            oPosMinus [ 1 ] = mLC [ 1 ];
            return;
        }

        iPosPlus.printYourself();
        OOFEM_ERROR("Mapping failed.")

        //    printf("iPosPlus: "); iPosPlus.printYourself();
        //    printf("oPosMinus: "); oPosMinus.printYourself();
    //#else
    } else {

#if 1

        const double distTol = 1.0e-12;
//        bool mappingPerformed = false;

        FloatArray n = mPeriodicityNormal;
        FloatArray t = {n(1),-n(0)};
        t.normalize();

        double l_s = mUC[0] - mLC[0];

        // Compute angle
        double alpha = 0.0, a = 0.0;
        if ( fabs(t(0)) > 1.0e-6 && fabs(t(1)) > 1.0e-6 ) {
            alpha = atan(t(1)/t(0));

            if ( alpha > 45.0*M_PI/180.0 ) {
                a = l_s/tan(alpha);
            } else {
                a = l_s*tan(alpha);
            }
        } else {
            // 90 degrees or 0 degrees
//            alpha = 1.57079632679490e+00;
            a = 0.0;
        }

        if ( alpha > 45.0*M_PI/180.0 ) {

            // alpha > 45 degrees

            if ( iPosPlus [ 0 ] > mUC [ 0 ] - distTol ) {
                // Gamma_1_plus
                oPosMinus = {0.0, iPosPlus[1]};
                return;
            }

            if ( iPosPlus [ 1 ] > mUC [ 1 ] - distTol ) {
                // Gamma_2_plus

                if ( iPosPlus[0] < a ) {
                    oPosMinus = {l_s - a + iPosPlus[0], 0.0};
                    return;
                } else {
                    oPosMinus = {iPosPlus[0] - a, 0.0};
                    return;
                }

            }

        } else {

            // alpha <= 45 degrees

            if ( iPosPlus [ 0 ] > mUC [ 0 ] - distTol ) {
                // Gamma_1_plus
                if ( iPosPlus[1] < a ) {
                    oPosMinus = {0.0, l_s - a + iPosPlus[1]};
                    return;
                } else {
                    oPosMinus = {0.0, iPosPlus[1] - a};
                    return;
                }

            }

            if ( iPosPlus [ 1 ] > mUC [ 1 ] - distTol ) {
                // Gamma_2_plus

                oPosMinus = {iPosPlus[0], 0.0};
                return;
            }

        }

#else

#endif
    }
    //#endif
}

void PrescribedGradientBCWeak :: giveMirroredPointOnGammaPlus(FloatArray &oPosPlus, const FloatArray &iPosMinus) const
{
//#if 0
    if ( mMirrorFunction == 0 ) {

        const double l_box = mUC(0) - mLC(0);

        oPosPlus = iPosMinus;
//        const double distTol = 1.0e-16;
        const double distTol = l_box*1.0e-10;

//        if ( distance(iPosMinus, mLC) < distTol ) {
//            printf("iPosMinus: %.12e %.12e\n", iPosMinus [ 0 ], iPosMinus [ 1 ]);
//            OOFEM_ERROR("Unmappable point.")
//        }

        if ( iPosMinus [ 0 ] < mLC [ 0 ] + distTol ) {
            oPosPlus [ 0 ] = mUC [ 0 ];
            return;
        } else if ( iPosMinus [ 1 ] < mLC [ 1 ] + distTol ) {
            oPosPlus [ 1 ] = mUC [ 1 ];
            return;
        }

        iPosMinus.printYourself();
        OOFEM_ERROR("Mapping failed.")
//#else
    } else {

#if 1

        const double distTol = 1.0e-12;

        FloatArray n = mPeriodicityNormal;
        FloatArray t = {n(1),-n(0)};
        t.normalize();

        double l_s = mUC[0] - mLC[0];

        // Compute angle
        double alpha = 0.0, a = 0.0;
        if ( fabs(t(0)) > 1.0e-6 && fabs(t(1)) > 1.0e-6 ) {
            alpha = atan(t(1)/t(0));

            if ( alpha > 45.0*M_PI/180.0 ) {
                a = l_s/tan(alpha);
            } else {
                a = l_s*tan(alpha);
            }
        } else {
            // 90 degrees
            a = 0.0;
        }

//	printf("t(1)/t(0): %e\n", t(1)/t(0));
//	printf("a: %e\n", a);

        if ( alpha > 45.0*M_PI/180.0 ) {

            // alpha > 45 degrees

            if ( iPosMinus [ 0 ] < mLC [ 0 ] + distTol ) {
                // Gamma_1_minus
                oPosPlus = {l_s, iPosMinus[1]};
                return;
            }

            if ( iPosMinus [ 1 ] < mLC [ 1 ] + distTol ) {
                // Gamma_2_minus

                if ( iPosMinus[0] < l_s - a) {
                    oPosPlus = {iPosMinus[0] + a, l_s};
                    return;
                }
                else {
                    oPosPlus = {iPosMinus[0] - (l_s - a), l_s};
                    return;
                }

            }

        } else {
            // alpha <= 45 degrees

            if ( iPosMinus [ 0 ] < mLC [ 0 ] + distTol ) {
                // Gamma_1_minus

                if ( iPosMinus[1] < l_s - a ) {
                    oPosPlus = {l_s, iPosMinus[1] + a};
                    return;
                } else {
                    oPosPlus = {l_s, iPosMinus[1] - (l_s - a) };
                    return;
                }
            }

            if ( iPosMinus [ 1 ] < mLC [ 1 ] + distTol ) {
                // Gamma_2_minus

                oPosPlus = {iPosMinus[0], l_s};
                return;

            }

        }
#endif
    }
//#endif
}


void PrescribedGradientBCWeak :: computeDomainBoundingBox(Domain &iDomain, FloatArray &oLC, FloatArray &oUC)
{
    // Compute LC and UC by assuming a rectangular domain.
    int numNodes = iDomain.giveNumberOfDofManagers();
    int nsd = iDomain.giveNumberOfSpatialDimensions();

    FloatArray lc = iDomain.giveDofManager(1)->giveCoordinates();
    FloatArray uc = iDomain.giveDofManager(1)->giveCoordinates();

    for ( int i = 1; i <= numNodes; i++ ) {
        DofManager *dMan = iDomain.giveDofManager(i);
        const auto &coord = dMan->giveCoordinates();

        for ( int j = 0; j < nsd; j++ ) {
            if ( coord [ j ] < lc [ j ] ) {
                lc [ j ] = coord [ j ];
            }

            if ( coord [ j ] > uc [ j ] ) {
                uc [ j ] = coord [ j ];
            }
        }
    }

    oLC = std::move(lc);
    oUC = std::move(uc);
}

int PrescribedGradientBCWeak :: giveSideIndex(const FloatArray &iPos) const
{
    const double distTol = 1.0e-12;

    if ( iPos [ 0 ] > mUC [ 0 ] - distTol ) {
        return 0;
    }

    if ( iPos [ 1 ] > mUC [ 1 ] - distTol ) {
        return 1;
    }

    if ( iPos [ 0 ] < mLC [ 0 ] + distTol ) {
        return 2;
    }

    if ( iPos [ 1 ] < mLC [ 1 ] + distTol ) {
        return 3;
    }

    OOFEM_ERROR("Could not identify side index.")

    return -1;
}

void PrescribedGradientBCWeak :: findHoleCoord(std::vector<FloatArray> &oHoleCoordUnsorted, std::vector<FloatArray> &oAllCoordUnsorted)
{
    Set *setPointer = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = setPointer->giveBoundaryList();

    // Loop over boundary nodes and check how many times they occur:
    // 1 -> at the edge of an inclusion, therefore must be retained
    // 2 -> connected to two segments, optional to keep

    std::unordered_map<int,int> map_bnd_node_ind_to_num_occurences;
    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {

    	int elIndex = boundaries.at(pos * 2 - 1);
        Element *e = this->giveDomain()->giveElement( elIndex );
        int boundary = boundaries.at(pos * 2);
        const auto &bNodes = e->giveInterpolation()->boundaryGiveNodes(boundary);
        DofManager *startNode   = e->giveDofManager(bNodes [ 0 ]);
        int startNodeInd = startNode->giveNumber();
        DofManager *endNode     = e->giveDofManager(bNodes [ 1 ]);
        int endNodeInd = endNode->giveNumber();


        auto res = map_bnd_node_ind_to_num_occurences.find(startNodeInd);
        if ( res != map_bnd_node_ind_to_num_occurences.end() ) {
            map_bnd_node_ind_to_num_occurences[startNodeInd]++;
        } else {
            map_bnd_node_ind_to_num_occurences[startNodeInd] = 1;
        }

        res = map_bnd_node_ind_to_num_occurences.find(endNodeInd);
        if ( res != map_bnd_node_ind_to_num_occurences.end() ) {
            map_bnd_node_ind_to_num_occurences[endNodeInd]++;
        } else {
            map_bnd_node_ind_to_num_occurences[endNodeInd] = 1;
        }

    }


    for ( auto it = map_bnd_node_ind_to_num_occurences.begin(); it != map_bnd_node_ind_to_num_occurences.end(); ++it ) {

        bool mandatory_to_keep = false;
        if ( it->second == 1 ) {
            mandatory_to_keep = true;
        }

        DofManager *bndNode = domain->giveDofManager(it->first);
        const auto &x = bndNode->giveCoordinates();
        FloatArray xPlus = x;

        if ( !boundaryPointIsOnActiveBoundary(x) ) {
            giveMirroredPointOnGammaPlus(xPlus, x);
        }

        if ( mandatory_to_keep ) {
            oHoleCoordUnsorted.push_back(xPlus);
        }

        oAllCoordUnsorted.push_back(xPlus);
    }
}


void PrescribedGradientBCWeak :: findCrackBndIntersecCoord(std::vector<FloatArray> &oHoleCoordUnsorted)
{
    Set *setPointer = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = setPointer->giveBoundaryList();

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);

        const auto &bNodes = e->giveInterpolation()->boundaryGiveNodes(boundary);

        // Add the start and end nodes of the segment
        DofManager *startNode = e->giveDofManager(bNodes [ 0 ]);
        const auto &xS = startNode->giveCoordinates();

        DofManager *endNode = e->giveDofManager(bNodes [ 1 ]);
        const auto &xE = endNode->giveCoordinates();

        FloatArray xC;
        xC.beScaled(0.5, xS);
        xC.add(0.5, xE);


        // Add the points where cracks intersect the boundary
        XfemElementInterface *xfemElInt = dynamic_cast< XfemElementInterface * >( e );

        if ( xfemElInt && domain->hasXfemManager() ) {
            std :: vector< Line >segments;
            std :: vector< FloatArray >intersecPoints;
            xfemElInt->partitionEdgeSegment(boundary, segments, intersecPoints, mTangDistPadding);

            for ( auto x : intersecPoints ) {

                if ( pointIsOnGammaPlus(x) ) {
                    oHoleCoordUnsorted.push_back(x);
                } else {
                    FloatArray xPlus;
                    giveMirroredPointOnGammaPlus(xPlus, x);
                    oHoleCoordUnsorted.push_back( std::move(xPlus) );
                }
            }
        }
    }

    // Also add traction nodes where cohesive zone elements intersect the boundary
    // TODO: Check later, this should already be covered by the new approach for
    //       identifying inclusions.
}

void PrescribedGradientBCWeak :: findPeriodicityCoord(std::vector<FloatArray> &oHoleCoordUnsorted)
{
    const double l_s = mUC[0] - mLC[0];

    FloatArray n = mPeriodicityNormal;
    FloatArray t = {n(1),-n(0)};

    if ( mMirrorFunction == 1 || mMirrorFunction == 2 ) {

        if ( fabs(n(1)) <= fabs(n(0)) ) {
            // a <= l_s/2
            double a = 0.5*l_s*( 1.0 +  n(1)/n(0) );

            FloatArray p1 = {2.0*a, 0.0};
            FloatArray p1Plus;
            giveMirroredPointOnGammaPlus(p1Plus, p1);
            oHoleCoordUnsorted.push_back(std::move(p1Plus));
        } else {
            // a > l_s/2
            double c = l_s - 0.5*l_s*( 1.0 + t(1)/t(0) );

            FloatArray p1 = {l_s, l_s-2.0*c};
            oHoleCoordUnsorted.push_back(std::move(p1));

            FloatArray p2 = {0, 2.0*c};
            FloatArray p2Plus;
            giveMirroredPointOnGammaPlus(p2Plus, p2);
            oHoleCoordUnsorted.push_back(std::move(p2Plus));
        }
    }
}


void PrescribedGradientBCWeak :: removeClosePoints(std::vector<FloatArray> &ioCoords, const double &iAbsTol)
{
    if ( ioCoords.size() == 0 ) {
        return;
    }

    const double tol2 = iAbsTol*iAbsTol;

    std::vector<FloatArray> tmp = { ioCoords[0] };
    size_t j = 0;

    for ( size_t i = 1; i < ioCoords.size(); i++ ) {
        if ( distance_square(ioCoords[i], tmp[j]) > tol2 ) {
            tmp.push_back(ioCoords[i]);
            j++;
        }
    }

    ioCoords = std::move(tmp);
}


void PrescribedGradientBCWeak :: removeSegOverHoles(TracSegArray &ioTSeg, const double &iAbsTol)
{
    // Idea: 	Loop over segments and check if the center point of each
    //			segment is located on an element. If not, the segment is
    //			located over a hole and needs to be removed.

    std :: vector< Line > tmp;

    SpatialLocalizer *localizer = domain->giveSpatialLocalizer();
    const double tol2 = iAbsTol*iAbsTol;

    for ( auto &l : ioTSeg.mInteriorSegments ) {
        const auto &xS = l.giveVertex(1);
        const auto &xE = l.giveVertex(2);
        FloatArray xPlus = {0.5*(xS[0]+xE[0]), 0.5*(xS[1]+xE[1])};

        FloatArray lcoordsPlus, closestPlus;
        localizer->giveElementClosestToPoint(lcoordsPlus, closestPlus, xPlus);

        FloatArray xMinus;
        giveMirroredPointOnGammaMinus(xMinus, xPlus);
        FloatArray lcoordsMinus, closestMinus;
        localizer->giveElementClosestToPoint(lcoordsMinus, closestMinus, xMinus);

        if ( !(distance_square(xPlus, closestPlus) > tol2 || distance_square(xMinus, closestMinus) > tol2) ) {
            tmp.push_back(l);
        }
    }

    ioTSeg.mInteriorSegments = std::move(tmp);
}

} /* namespace oofem */
