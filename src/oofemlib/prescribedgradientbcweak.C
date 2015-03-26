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

#include "xfem/XFEMDebugTools.h"

#include <cmath>
#include <string>
#include <sstream>

namespace oofem {
PrescribedGradientBCWeak :: PrescribedGradientBCWeak(int n, Domain *d) :
    ActiveBoundaryCondition(n, d),
    PrescribedGradientHomogenization(),
    mTractionInterpOrder(0),
    mNumTractionNodesAtIntersections(1),
    mTractionNodeSpacing(1),
    mMeshIsPeriodic(false),
    mDuplicateCornerNodes(false),
    mTangDistPadding(0.0),
    mTracDofScaling(1.0e0),
    mpDisplacementLock(NULL),
    mLockNodeInd(0),
    mDispLockScaling(1.0)
{
    // Compute bounding box of the domain
    computeDomainBoundingBox(* d, mLC, mUC);
}

PrescribedGradientBCWeak :: ~PrescribedGradientBCWeak()
{
    clear();
}

void PrescribedGradientBCWeak :: clear()
{
    mpTractionNodes.clear();

    for ( size_t i = 0; i < mpTractionMasterNodes.size(); i++ ) {
        if ( mpTractionMasterNodes [ i ] != NULL ) {
            delete mpTractionMasterNodes [ i ];
            mpTractionMasterNodes [ i ] = NULL;
        }
    }
    mpTractionMasterNodes.clear();

    for ( size_t i = 0; i < mpTractionElements.size(); i++ ) {
        if ( mpTractionElements [ i ] != NULL ) {
            delete mpTractionElements [ i ];
            mpTractionElements [ i ] = NULL;
        }
    }
    mpTractionElements.clear();

    if ( mpDisplacementLock != NULL ) {
        delete mpDisplacementLock;
        mpDisplacementLock = NULL;
    }

    mMapTractionElDispElGamma.clear();
    mTractionElInteriorCoordinates.clear();
    mTracElDispNodes.clear();
}

//#define DAMAGE_TEST

int PrescribedGradientBCWeak :: giveNumberOfInternalDofManagers()
{
    int numDMan = mpTractionMasterNodes.size();

    if ( mpDisplacementLock != NULL ) {
        numDMan++;
    }

    return numDMan;
}

DofManager *PrescribedGradientBCWeak :: giveInternalDofManager(int i)
{
    if ( i - 1 < int( mpTractionMasterNodes.size() ) ) {
        return mpTractionMasterNodes [ i - 1 ];
    } else   {
        return mpDisplacementLock;
    }
}

IRResultType PrescribedGradientBCWeak :: initializeFrom(InputRecord *ir)
{
    IRResultType result;
    result = ActiveBoundaryCondition :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }
    result = PrescribedGradientHomogenization :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    IR_GIVE_FIELD(ir, mTractionInterpOrder, _IFT_PrescribedGradientBCWeak_TractionInterpOrder);
    printf("mTractionInterpOrder: %d\n", mTractionInterpOrder);

    IR_GIVE_FIELD(ir, mNumTractionNodesAtIntersections, _IFT_PrescribedGradientBCWeak_NumTractionNodesAtIntersections);
    printf("mNumTractionNodesAtIntersections: %d\n", mNumTractionNodesAtIntersections);

    if ( mNumTractionNodesAtIntersections > 1 && mTractionInterpOrder == 0 ) {
        OOFEM_ERROR("mNumTractionNodesAtIntersections > 1 is not allowed if mTractionInterpOrder == 0.")
    }

    IR_GIVE_FIELD(ir, mTractionNodeSpacing, _IFT_PrescribedGradientBCWeak_NumTractionNodeSpacing);
    printf("mTractionNodeSpacing: %d\n", mTractionNodeSpacing);

    int duplicateCorners = 0;
    IR_GIVE_FIELD(ir, duplicateCorners, _IFT_PrescribedGradientBCWeak_DuplicateCornerNodes);
    printf("duplicateCorners: %d\n", duplicateCorners);

    if ( duplicateCorners == 1 ) {
        mDuplicateCornerNodes = true;
    } else   {
        mDuplicateCornerNodes = false;
    }

    mTangDistPadding = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, mTangDistPadding, _IFT_PrescribedGradientBCWeak_TangDistPadding);
    printf("mTangDistPadding: %e\n", mTangDistPadding);

    IR_GIVE_OPTIONAL_FIELD(ir, mTracDofScaling, _IFT_PrescribedGradientBCWeak_TracDofScaling);
    printf("mTracDofScaling: %e\n", mTracDofScaling );

    return IRRT_OK;
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
    } else   {
        input.setField(0, _IFT_PrescribedGradientBCWeak_DuplicateCornerNodes);
    }

    input.setField(mTangDistPadding, _IFT_PrescribedGradientBCWeak_TangDistPadding);
    input.setField(mTracDofScaling, _IFT_PrescribedGradientBCWeak_TracDofScaling);
}

void PrescribedGradientBCWeak :: postInitialize()
{}

void PrescribedGradientBCWeak :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                                CharType type, ValueModeType mode,
                                                const UnknownNumberingScheme &s, FloatArray *eNorm)
{
    int dim = domain->giveNumberOfSpatialDimensions();

    if ( type == ExternalForcesVector ) {
        // The external force vector is given by
        // f_ext = int N^trac H . (x - x_c)

        for ( size_t i = 0; i < mpTractionElements.size(); i++ ) {
            // Compute f_ext
            FloatArray fExt;

            if ( mTractionInterpOrder == 0 ) {
                fExt.resize(1 * dim);
            } else if ( mTractionInterpOrder == 1 )    {
                fExt.resize(2 * dim);
            }
            fExt.zero();


            const TractionElement *el = mpTractionElements [ i ];
            IntegrationRule *ir = createNewIntegrationRule(i);

            for ( GaussPoint *gp: *ir ) {
                const FloatArray &locCoordsOnLine = gp->giveNaturalCoordinates();

                // Compute N^trac
                FloatArray N, Ntrac;
                computeNTraction(Ntrac, locCoordsOnLine [ 0 ], * el);
                el->computeN_Linear(N, locCoordsOnLine [ 0 ]);

                // Compute x
                FloatArray x( el->mStartCoord.giveSize() );
                x.zero();

                x.add(N [ 0 ], el->mStartCoord);
                x.add(N [ 1 ], el->mEndCoord);

                // Compute H.(x - x_c) or H.[x]
                FloatArray temp;
                giveBoundaryCoordVector(temp, x);

                FloatArray Hx;
                Hx.beProductOf(mGradient, temp);


                // N-matrix
                FloatMatrix Nmat;
                Nmat.beNMatrixOf(Ntrac, dim);


                // Add contribution to fExt
                FloatArray contrib;
                contrib.beTProductOf(Nmat, Hx);
                double detJ = 0.5 * el->mStartCoord.distance(el->mEndCoord);
                fExt.add(detJ * gp->giveWeight(), contrib);
            }

            // Fetch location arrays for current traction element
            IntArray rows;
            giveTractionLocationArrays(i, rows, type, s);

            fExt.negated();

            double loadLevel = this->giveTimeFunction()->evaluateAtTime(tStep->giveTargetTime());
            fExt.times(loadLevel);

            // Assemble
            answer.assemble(fExt, rows);

            delete ir;
        }
    } else if ( type == InternalForcesVector ) {
        FloatArray fe_trac, fe_disp;

        for ( size_t i = 0; i < mpTractionElements.size(); i++ ) {
            // Use the tangent to reduce duplication of code
            FloatMatrix Ke;
            this->integrateTangent(Ke, i);

            // Compute vector of traction unknowns
            FloatArray tracUnknowns;
            giveTractionUnknows(tracUnknowns, mode, tStep, i);

            fe_disp.beTProductOf(Ke, tracUnknowns);
            fe_disp.negated();


            // Compute vector of displacement unknowns
            FloatArray dispUnknowns;
            giveDisplacementUnknows(dispUnknowns, mode, tStep, i);

#ifdef DEBUG
            if ( !dispUnknowns.isFinite() ) {
                OOFEM_ERROR("!dispUnknowns.isFinite()")
            }
#endif

            fe_trac.beProductOf(Ke, dispUnknowns);
            fe_trac.negated();


            // Fetch location arrays
            IntArray tracRrows;
            giveTractionLocationArrays(i, tracRrows, type, s);

            IntArray dispRrows;
            giveDisplacementLocationArrays(i, dispRrows, type, s);

#ifdef DEBUG
            if ( !fe_trac.isFinite() ) {
                OOFEM_ERROR("!fe_trac.isFinite()")
            }

            if ( !fe_disp.isFinite() ) {
                OOFEM_ERROR("!fe_trac.isFinite()")
            }
#endif
            // Assemble
            answer.assemble(fe_trac, tracRrows);
            answer.assemble(fe_disp, dispRrows);
        }

        if ( mpDisplacementLock != NULL ) {
            IntArray dispLockRows;
            mpDisplacementLock->giveCompleteLocationArray(dispLockRows, s);

            FloatArray fe_dispLock;

            int lockNodePlaceInArray = domain->giveDofManPlaceInArray(mLockNodeInd);
            FloatArray nodeUnknowns;
            domain->giveDofManager(lockNodePlaceInArray)->giveCompleteUnknownVector(nodeUnknowns, mode, tStep);

            for ( int i = 0; i < domain->giveNumberOfSpatialDimensions(); i++ ) {
                fe_dispLock.push_back(nodeUnknowns [ i ]);
            }

            fe_dispLock.times(mDispLockScaling);

            answer.assemble(fe_dispLock, dispLockRows);
        }
    }
}


void PrescribedGradientBCWeak :: assemble(SparseMtrx &answer, TimeStep *tStep,
                                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    if ( type == TangentStiffnessMatrix || type == SecantStiffnessMatrix || type == ElasticStiffnessMatrix ) {
        FloatMatrix Ke, KeT;
        IntArray tracRows, tracCols, dispRows, dispCols;

        for ( size_t i = 0; i < mpTractionElements.size(); i++ ) {
            // Rows and columns for displacement and traction contributions
            giveTractionLocationArrays(i, tracRows, type, r_s);
            giveDisplacementLocationArrays(i, dispRows, type, r_s);

            giveTractionLocationArrays(i, tracCols, type, c_s);
            giveDisplacementLocationArrays(i, dispCols, type, c_s);

            this->integrateTangent(Ke, i);
            Ke.negated();
            KeT.beTranspositionOf(Ke);

            answer.assemble(tracRows, dispCols, Ke);
            answer.assemble(dispRows, tracCols, KeT);


            FloatMatrix KZero( tracRows.giveSize(), tracCols.giveSize() );
            KZero.zero();
            answer.assemble(tracRows, tracCols, KZero);
        }


        if ( mpDisplacementLock != NULL ) {
            int nsd = domain->giveNumberOfSpatialDimensions();
            FloatMatrix KeDispLock(nsd, nsd);
            KeDispLock.beUnitMatrix();
            KeDispLock.times(mDispLockScaling);

            int placeInArray = domain->giveDofManPlaceInArray(mLockNodeInd);
            DofManager *node = domain->giveDofManager(placeInArray);

            IntArray lockRows, lockCols, nodeRows, nodeCols;
            mpDisplacementLock->giveCompleteLocationArray(lockRows, r_s);
            node->giveCompleteLocationArray(nodeRows, r_s);

            // TODO: Can be done nicer by prescribing which dof IDs to lock.
            IntArray nodeRowsRed;
            for ( int m = 0; m < nsd; m++ ) {
                nodeRowsRed.followedBy(nodeRows [ m ]);
            }

            mpDisplacementLock->giveCompleteLocationArray(lockCols, c_s);
            node->giveCompleteLocationArray(nodeCols, c_s);

            IntArray nodeColsRed;
            for ( int m = 0; m < nsd; m++ ) {
                nodeColsRed.followedBy(nodeCols [ m ]);
            }

            answer.assemble(lockRows, nodeColsRed, KeDispLock);
            answer.assemble(nodeRowsRed, lockCols, KeDispLock);

            FloatMatrix KZero( lockRows.giveSize(), lockCols.giveSize() );
            KZero.zero();
            answer.assemble(lockRows, lockCols, KZero);
        }
    } else {
        printf("Skipping assembly in PrescribedGradientBCWeak::assemble().\n");
    }
}

void PrescribedGradientBCWeak :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                                                    const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    rows.clear();
    cols.clear();

    // Loop over traction elements
    for ( size_t tracElInd = 0; tracElInd < mpTractionElements.size(); tracElInd++ ) {
        IntArray tracElCols, tracElRows, trac_loc_c, trac_loc_r;

        const TractionElement &tEl = * ( mpTractionElements [ tracElInd ] );
        for ( int tracNodeInd : tEl.mTractionNodeInd ) {
            Node *tNode = mpTractionNodes [ tracNodeInd ];
            tNode->giveCompleteLocationArray(trac_loc_r, r_s);
            tracElRows.followedBy(trac_loc_r);

            tNode->giveCompleteLocationArray(trac_loc_c, c_s);
            tracElCols.followedBy(trac_loc_c);
        }


        // Fetch displacement elements that interact
        // with the current traction element.
        const std :: vector< int > &dispElInds = mMapTractionElDispElGamma [ tracElInd ];

        IntArray dispElCols, dispElRows, disp_loc_c, disp_loc_r;

        for ( int dispElInd : dispElInds ) {
            Element *e = this->giveDomain()->giveElement(dispElInd);
            e->giveLocationArray(disp_loc_r, r_s);
            dispElRows.followedBy(disp_loc_r);

            e->giveLocationArray(disp_loc_c, c_s);
            dispElCols.followedBy(disp_loc_c);
        }

        rows.push_back(dispElRows);
        cols.push_back(tracElCols);

        rows.push_back(tracElRows);
        cols.push_back(dispElCols);
    }


    if ( mpDisplacementLock != NULL ) {
        IntArray dispLock_r, dispLock_c;
        mpDisplacementLock->giveCompleteLocationArray(dispLock_r, r_s);
        mpDisplacementLock->giveCompleteLocationArray(dispLock_c, c_s);

        int nodeInd = 1;
        DofManager *node = domain->giveDofManager(nodeInd);
        IntArray node_r, node_c;
        node->giveCompleteLocationArray(node_r, r_s);
        node->giveCompleteLocationArray(node_c, c_s);

        rows.push_back(dispLock_r);
        cols.push_back(node_c);

        rows.push_back(node_r);
        cols.push_back(dispLock_c);
    }
}

void PrescribedGradientBCWeak :: giveTractionLocationArray(IntArray &rows,
                                                           const UnknownNumberingScheme &s)
{
    rows.clear();

    // Loop over traction elements
    for ( size_t tracElInd = 0; tracElInd < mpTractionElements.size(); tracElInd++ ) {
        IntArray tracElRows, trac_loc_r;

        const TractionElement &tEl = * ( mpTractionElements [ tracElInd ] );
        for ( int tracNodeInd : tEl.mTractionNodeInd ) {
            Node *tNode = mpTractionNodes [ tracNodeInd ];
            tNode->giveCompleteLocationArray(trac_loc_r, s);
            tracElRows.followedBy(trac_loc_r);
        }

        rows.followedBy(tracElRows);
    }


    if ( mpDisplacementLock != NULL ) {
        IntArray dispLock_r;
        mpDisplacementLock->giveCompleteLocationArray(dispLock_r, s);

        rows.followedBy(dispLock_r);
    }
}

void PrescribedGradientBCWeak :: giveTractionLocationArrays(int iTracElInd, IntArray &rows, CharType type,
                                                            const UnknownNumberingScheme &s)
{
    rows.clear();

    IntArray tracElRows, trac_loc_r;

    const TractionElement &tEl = * ( mpTractionElements [ iTracElInd ] );
    for ( int tracNodeInd : tEl.mTractionNodeInd ) {
        Node *tNode = mpTractionNodes [ tracNodeInd ];
        tNode->giveCompleteLocationArray(trac_loc_r, s);
        tracElRows.followedBy(trac_loc_r);
    }

    rows = tracElRows;
}

void PrescribedGradientBCWeak :: giveDisplacementLocationArrays(int iTracElInd, IntArray &rows, CharType type,
                                                                const UnknownNumberingScheme &s)
{
    rows.clear();

    for ( int nodeInd : mTracElDispNodes [ iTracElInd ] ) {
        IntArray nodeLocationArray;
        domain->giveDofManager(nodeInd)->giveCompleteLocationArray(nodeLocationArray, s);
        rows.followedBy(nodeLocationArray);
    }
}

void PrescribedGradientBCWeak :: computeField(FloatArray &sigma, TimeStep *tStep)
{
    double dSize = domainSize(this->giveDomain(), this->giveSetNumber());

    const int dim = domain->giveNumberOfSpatialDimensions();
    FloatMatrix stressMatrix(dim, dim);


    for ( size_t i = 0; i < mpTractionElements.size(); i++ ) {
        const TractionElement *el = mpTractionElements [ i ];
        IntegrationRule *ir = createNewIntegrationRule(i);

        for ( GaussPoint *gp: *ir ) {
            const FloatArray &locCoordsOnLine = gp->giveNaturalCoordinates();

            // Compute N^trac
            FloatArray N, Ntrac;
            computeNTraction(Ntrac, locCoordsOnLine [ 0 ], * el);
            el->computeN_Linear(N, locCoordsOnLine [ 0 ]);

            // Compute x
            FloatArray x( el->mStartCoord.giveSize() );
            x.zero();

            x.add(N [ 0 ], el->mStartCoord);
            x.add(N [ 1 ], el->mEndCoord);

            // N-matrix
            FloatMatrix Nmat;
            Nmat.beNMatrixOf(Ntrac, dim);

            // Interpolate traction
            FloatArray tracUnknowns;
            giveTractionUnknows(tracUnknowns, VM_Total, tStep, i);

            FloatArray traction;
            traction.beProductOf(Nmat, tracUnknowns);

            FloatArray tmp;
            giveBoundaryCoordVector(tmp, x);

            FloatMatrix contrib;
            contrib.beDyadicProductOf(traction, tmp);

            double detJ = 0.5 * el->mStartCoord.distance(el->mEndCoord);
            contrib.times( detJ * gp->giveWeight() );

            for ( int m = 0; m < dim; m++ ) {
                for ( int n = 0; n < dim; n++ ) {
                    stressMatrix(m, n) += contrib(m, n);
                }
            }
        }


        delete ir;
    }


    if ( dim == 2 ) {
        sigma = {
            stressMatrix(0, 0), stressMatrix(1, 1), stressMatrix(0, 1), stressMatrix(1, 0)
        };
    } else   {
        sigma.beVectorForm(stressMatrix);
    }

    sigma.times(1.0 / dSize);

#if 0
    FloatArray tx, ty;
    for ( Node *node : mpTractionNodes ) {
        FloatArray tNode;
        node->giveCompleteUnknownVector(tNode, VM_Total, tStep);
        tx.push_back(tNode [ 0 ]);
        ty.push_back(tNode [ 1 ]);
    }

    printf("\n\n\n");
    printf("tx: ");
    tx.printYourself();
    printf("\n\n\n");
    printf("ty: ");
    ty.printYourself();
#endif
}

void PrescribedGradientBCWeak :: computeTangent(FloatMatrix& E, TimeStep* tStep)
{
    OOFEM_ERROR("Not implemented yet.");
}

void PrescribedGradientBCWeak :: giveTractionElNormal(size_t iElInd, FloatArray &oNormal, FloatArray &oTangent) const
{
    const FloatArray &xS = mpTractionElements [ iElInd ]->mStartCoord;
    const FloatArray &xE = mpTractionElements [ iElInd ]->mEndCoord;

    oTangent.beDifferenceOf(xE, xS);
    oTangent.normalize();

    oNormal = {
        oTangent [ 1 ], -oTangent [ 0 ]
    };
}

void PrescribedGradientBCWeak :: giveTractionElArcPos(size_t iElInd, double &oXiStart, double &oXiEnd) const
{
    const FloatArray &xS = mpTractionElements [ iElInd ]->mStartCoord;
    const FloatArray &xE = mpTractionElements [ iElInd ]->mEndCoord;

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
    mpTractionNodes [ mpTractionElements [ iElInd ]->mTractionNodeInd [ 0 ] ]->giveCompleteUnknownVector(oStartTraction, mode, tStep);

    if ( mpTractionElements [ iElInd ]->mTractionNodeInd.size() < 2 ) {
        mpTractionNodes [ mpTractionElements [ iElInd ]->mTractionNodeInd [ 0 ] ]->giveCompleteUnknownVector(oEndTraction, mode, tStep);
    } else   {
        mpTractionNodes [ mpTractionElements [ iElInd ]->mTractionNodeInd [ 1 ] ]->giveCompleteUnknownVector(oEndTraction, mode, tStep);
    }
}

void PrescribedGradientBCWeak :: recomputeTractionMesh()
{
    printf("Recomputing traction mesh.\n");
    clear();
    postInitialize();
}

void PrescribedGradientBCWeak :: createTractionMesh(bool iEnforceCornerPeriodicity, int iNumSides)
{
    const double nodeDistTol = 1.0e-15;
    const double meshTol = 1.0e-8; // Minimum distance between traction nodes

    /**
     * first:   coordinates
     * second:  bool telling if the point must be included in the
     *          traction mesh (e.g. a corner node),
     *          or if it can be omitted.
     */
    // Side 1: x = L, side 1: y = L
    std :: vector< std :: vector< std :: pair< FloatArray, bool > > >bndNodeCoords;
    std :: vector< std :: pair< FloatArray, bool > >emptyVec;

    for ( int i = 0; i < iNumSides; i++ ) {
        bndNodeCoords.push_back(emptyVec);
    }

    Set *setPointer = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = setPointer->giveBoundaryList();
    IntArray bNodes;

    // Loop over all boundary segments twice:
    // first add mesh points...
    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);

        e->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);

        // Add the start and end nodes of the segment
        DofManager *startNode   = e->giveDofManager(bNodes [ 0 ]);
        const FloatArray &xS    = * ( startNode->giveCoordinates() );

        DofManager *endNode     = e->giveDofManager(bNodes [ 1 ]);
        const FloatArray &xE    = * ( endNode->giveCoordinates() );

        FloatArray xC;
        xC.beScaled(0.5, xS);
        xC.add(0.5, xE);

        const double meshTol2 = meshTol * meshTol;

        if ( boundaryPointIsOnActiveBoundary(xC) ) {
            int sideInd = giveSideIndex(xC);

            if ( !closePointExists(bndNodeCoords [ sideInd ], xS, meshTol2) ) {
                std :: pair< FloatArray, bool >nodeCoord = {
                    xS, false
                };
                bndNodeCoords [ sideInd ].push_back(nodeCoord);
            }

            if ( !closePointExists(bndNodeCoords [ sideInd ], xE, meshTol2) ) {
                std :: pair< FloatArray, bool >nodeCoord = {
                    xE, false
                };
                bndNodeCoords [ sideInd ].push_back(nodeCoord);
            }
        } // if pointIsOnGammaPlus
        else {
            if ( pointIsMapapble(xS) ) {
                if ( pointIsMapapble(xS) ) {
                    FloatArray xSPlus;
                    giveMirroredPointOnGammaPlus(xSPlus, xS);

                    int sideIndS = giveSideIndex(xSPlus);

                    if ( !closePointExists(bndNodeCoords [ sideIndS ], xSPlus, meshTol2) ) {
                        std :: pair< FloatArray, bool >nodeCoord = {
                            xSPlus, false
                        };
                        bndNodeCoords [ sideIndS ].push_back(nodeCoord);
                    }
                }



                if ( pointIsMapapble(xE) ) {
                    FloatArray xEPlus;
                    giveMirroredPointOnGammaPlus(xEPlus, xE);

                    int sideIndE = giveSideIndex(xEPlus);

                    if ( !closePointExists(bndNodeCoords [ sideIndE ], xEPlus, meshTol2) ) {
                        std :: pair< FloatArray, bool >nodeCoord = {
                            xEPlus, false
                        };
                        bndNodeCoords [ sideIndE ].push_back(nodeCoord);
                    }
                }
            }
        }
    }

    // ...and then add points where cracks intersect the boundary
    // (by doing this in two steps with two loops, we avoid getting
    // in trouble if cracks intersect the boundary close to domain corners.
    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);

        e->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);

        // Add the start and end nodes of the segment
        DofManager *startNode   = e->giveDofManager(bNodes [ 0 ]);
        const FloatArray &xS    = * ( startNode->giveCoordinates() );

        DofManager *endNode     = e->giveDofManager(bNodes [ 1 ]);
        const FloatArray &xE    = * ( endNode->giveCoordinates() );

        FloatArray xC;
        xC.beScaled(0.5, xS);
        xC.add(0.5, xE);


        // Add traction nodes where cracks intersect the boundary if desired
        XfemElementInterface *xfemElInt = dynamic_cast< XfemElementInterface * >( e );

        if ( xfemElInt != NULL && domain->hasXfemManager() ) {
            std :: vector< Line >segments;
            std :: vector< FloatArray >intersecPoints;
            xfemElInt->partitionEdgeSegment(boundary, segments, intersecPoints, mTangDistPadding);

            for ( size_t i = 0; i < intersecPoints.size(); i++ ) {
                if ( boundaryPointIsOnActiveBoundary(intersecPoints [ i ]) ) {
                    int sideInd = giveSideIndex(intersecPoints [ i ]);

                    int numClosePoints = 0;
                    const double meshTol2 = meshTol * meshTol;
                    for ( auto &bndPos : bndNodeCoords [ sideInd ] ) {
                        if ( bndPos.first.distance_square(intersecPoints [ i ]) < meshTol2 ) {

                            if(numClosePoints < mNumTractionNodesAtIntersections) {
                                bndPos.second = true;
                            }

                            numClosePoints++;
                        }
                    }

                    if ( numClosePoints < mNumTractionNodesAtIntersections ) {
                        for ( int j = 0; j < mNumTractionNodesAtIntersections - numClosePoints; j++ ) {
                            std :: pair< FloatArray, bool >nodeCoord = {
                                intersecPoints [ i ], true
                            };
                            bndNodeCoords [ sideInd ].push_back(nodeCoord);
                        }
                    }
                } else   {
                    if ( mMeshIsPeriodic ) {

                        if ( pointIsMapapble(intersecPoints [ i ]) ) {
                            FloatArray xPlus;
                            giveMirroredPointOnGammaPlus(xPlus, intersecPoints [ i ]);

                            int sideInd = giveSideIndex(xPlus);

                            int numClosePoints = 0;
                            const double meshTol2 = meshTol * meshTol;
                            for ( auto &bndPos : bndNodeCoords [ sideInd ] ) {
                                if ( bndPos.first.distance_square(xPlus) < meshTol2 ) {

                                    if(numClosePoints < mNumTractionNodesAtIntersections) {
                                        bndPos.second = true;
                                    }

                                    numClosePoints++;
                                }
                            }

                            if ( numClosePoints < mNumTractionNodesAtIntersections ) {
                                for ( int j = 0; j < mNumTractionNodesAtIntersections - numClosePoints; j++ ) {
                                    std :: pair< FloatArray, bool >nodeCoord = {
                                        xPlus, true
                                    };
                                    bndNodeCoords [ sideInd ].push_back(nodeCoord);
                                }
                            }
                        }

                    }
                }
            }
        }




    }



    SpatialLocalizer *localizer = domain->giveSpatialLocalizer();

    // Also add traction nodes where cohesive zone elements intersect the boundary
    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);

        e->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);

        // Add the start and end nodes of the segment
        DofManager *startNode   = e->giveDofManager(bNodes [ 0 ]);
        const FloatArray &xS    = * ( startNode->giveCoordinates() );

        DofManager *endNode     = e->giveDofManager(bNodes [ 1 ]);
        const FloatArray &xE    = * ( endNode->giveCoordinates() );

        double radius = xS.distance(xE)*1.0e-2;

        std :: set< int >elListS;
        localizer->giveAllElementsWithNodesWithinBox(elListS, xS, radius);

        std :: vector< FloatArray >intersecPoints;


        // Also add traction nodes where cohesive zone elements intersect the boundary
        for(int elInd : elListS) {

            Element *el = domain->giveElement(elInd);

            if( strcmp(el->giveClassName(),"IntElLine1" ) == 0 || strcmp(el->giveClassName(),"IntElLine2" ) == 0 )
            {
#ifdef DAMAGE_TEST
                if(damageExceedsTolerance(el)) {
//                    OOFEM_ERROR("Damage exceeds tolerance.")
//                    printf("Damage exceeds tolerance. Adding traction node.\n");
                    intersecPoints.push_back(xS);
                    break;
                }

#else
                intersecPoints.push_back(xS);
                break;

#endif
            }

        }


        std :: set< int >elListE;
        localizer->giveAllElementsWithNodesWithinBox(elListE, xE, radius);



        for(int elInd : elListE) {

            Element *el = domain->giveElement(elInd);

            if( strcmp(el->giveClassName(),"IntElLine1" ) == 0 || strcmp(el->giveClassName(),"IntElLine2" ) == 0 )
            {

#ifdef DAMAGE_TEST
                if(damageExceedsTolerance(el)) {
//                    OOFEM_ERROR("Damage exceeds tolerance.")
//                    printf("Damage exceeds tolerance. Adding traction node.\n");
                    intersecPoints.push_back(xE);
                    break;
                }
#else
                intersecPoints.push_back(xE);
                break;
#endif


            }

        }



        for ( size_t i = 0; i < intersecPoints.size(); i++ ) {
            if ( boundaryPointIsOnActiveBoundary(intersecPoints [ i ]) ) {
                int sideInd = giveSideIndex(intersecPoints [ i ]);

                int numClosePoints = 0;
                const double meshTol2 = meshTol * meshTol;
                for ( auto &bndPos : bndNodeCoords [ sideInd ] ) {
                    if ( bndPos.first.distance_square(intersecPoints [ i ]) < meshTol2 ) {

                        if(numClosePoints < mNumTractionNodesAtIntersections) {
                            bndPos.second = true;
                        }

                        numClosePoints++;
                    }
                }

                if ( numClosePoints < mNumTractionNodesAtIntersections ) {
                    for ( int j = 0; j < mNumTractionNodesAtIntersections - numClosePoints; j++ ) {
                        std :: pair< FloatArray, bool >nodeCoord = {
                            intersecPoints [ i ], true
                        };

                        bndNodeCoords [ sideInd ].push_back(nodeCoord);
                    }
                }
            } else   {
                if ( mMeshIsPeriodic ) {

                    if ( pointIsMapapble(intersecPoints [ i ]) ) {

                        FloatArray xPlus;
                        giveMirroredPointOnGammaPlus(xPlus, intersecPoints [ i ]);

                        int sideInd = giveSideIndex(xPlus);

                        int numClosePoints = 0;
                        const double meshTol2 = meshTol * meshTol;
                        for ( auto &bndPos : bndNodeCoords [ sideInd ] ) {
                            if ( bndPos.first.distance_square(xPlus) < meshTol2 ) {

                                if(numClosePoints < mNumTractionNodesAtIntersections) {
                                    bndPos.second = true;
                                }
                                numClosePoints++;
                            }
                        }

                        if ( numClosePoints < mNumTractionNodesAtIntersections ) {
                            for ( int j = 0; j < mNumTractionNodesAtIntersections - numClosePoints; j++ ) {
                                std :: pair< FloatArray, bool >nodeCoord = {
                                    xPlus, true
                                };
                                bndNodeCoords [ sideInd ].push_back(nodeCoord);
                            }
                        }
                    }

                }
            }
        }

    }

    // Sort boundary nodes
    for ( size_t arrayInd = 0; arrayInd < bndNodeCoords.size(); arrayInd++ ) {
        std :: sort( bndNodeCoords [ arrayInd ].begin(), bndNodeCoords [ arrayInd ].end(), ArcPosSortFunction3< bool >( mLC, mUC, nodeDistTol, int( arrayInd ) ) );

        // Make sure that the first and last point on each side are retained.
        bndNodeCoords [ arrayInd ] [ 0 ].second = true;
        bndNodeCoords [ arrayInd ].back().second = true;


#if 0
        printf("\n\ncoordArray: ");
        for ( auto pos : coordArray ) {
            pos.first.printYourself();
        }
#endif
    }


    // Create traction dofs
    int nsd = domain->giveNumberOfSpatialDimensions();
    std :: vector< int >dofIds;
    for ( int j = 0; j < nsd; j++ ) {
        dofIds.push_back( this->domain->giveNextFreeDofID() );
    }

    std :: vector< FloatArray >tractionNodeCoord;
    int numNodes = domain->giveNumberOfDofManagers();
    int numPointsPassed = 0;

    int totNodesCreated = 0;
    for ( auto coordArray: bndNodeCoords ) {
        int startNodeInd = mpTractionNodes.size();

        //////////////////////////////////////////////////////
        // Create first node
        totNodesCreated++;
        numPointsPassed++;
        Node *firstNode = new Node(numNodes + 1, domain);
        firstNode->setGlobalNumber(numNodes + 1);
        for ( auto &dofId: dofIds ) {
            firstNode->appendDof( new MasterDof(firstNode, ( DofIDItem ) dofId) );
        }


        firstNode->setCoordinates(coordArray [ 0 ].first);
        mpTractionNodes.push_back(firstNode);
        mpTractionMasterNodes.push_back(firstNode);

        tractionNodeCoord.push_back(coordArray [ 0 ].first);

        numNodes++;
        //////////////////////////////////////////////////////


        std :: vector< FloatArray >coordsToKeep;
        for ( size_t i = 0; i < coordArray.size(); i++ ) {
            numPointsPassed++;

            if ( coordArray [ i ].second || numPointsPassed >= mTractionNodeSpacing ) {
                numPointsPassed = 0;
                coordsToKeep.push_back(coordArray [ i ].first);
            }
        }


        for ( size_t i = 1; i < coordsToKeep.size(); i++ ) {
            // Create the second node if desired
            numPointsPassed++;

            if ( !( ( i == ( coordsToKeep.size() - 1 ) ) && mTractionInterpOrder == 0 ) ) {
                totNodesCreated++;

                numPointsPassed = 0;

                bool createSlaveNode = false;

                int masterInd = 0;
                if ( ( ( i == coordsToKeep.size() - 1 ) && mTractionInterpOrder == 1 ) || ( ( i == coordsToKeep.size() - 2 ) && mTractionInterpOrder == 0 ) ) {
                    //                    printf("Creating slave node for i: %lu\n", i);
                    createSlaveNode = true;
                    masterInd = startNodeInd;
                }



                if ( !iEnforceCornerPeriodicity ) {
                    createSlaveNode = false;
                }

                if ( mTractionInterpOrder == 0 ) {
                    createSlaveNode = false;
                }

                if ( createSlaveNode ) {
                    Node *masterNode = mpTractionNodes [ masterInd ];
                    //                      printf("Creating a slave of %d with coord: ",masterInd ); coordsToKeep[i].printYourself();

                    mpTractionNodes.push_back(masterNode);
                    tractionNodeCoord.push_back(coordsToKeep [ i ]);
                } else   {
                    Node *node = new Node(numNodes + 1, domain);
                    node->setGlobalNumber(numNodes + 1);
                    for ( auto &dofid: dofIds ) {
                        node->appendDof( new MasterDof(node, ( DofIDItem ) dofid) );
                    }

                    //                      printf("Creating master node with coord: "); coordsToKeep[i].first.printYourself();

                    node->setCoordinates(coordsToKeep [ i ]);
                    mpTractionNodes.push_back(node);
                    mpTractionMasterNodes.push_back(node);

                    tractionNodeCoord.push_back(coordsToKeep [ i ]);

                    numNodes++;
                }
            }




            // Create traction elements
            if ( mTractionInterpOrder == 0 ) {
                // Piecewise constant traction
                // (Not stable in terms of the LBB condition,
                //  but interesting for comparison.)


                if ( i == coordsToKeep.size() - 1 ) {
                    TractionElement *tractionEl = new TractionElement();

                    tractionEl->mTractionNodeInd.push_back(mpTractionNodes.size() - 1);
                    tractionEl->mStartCoord = tractionNodeCoord [ mpTractionNodes.size() - 1 ];

                    tractionEl->mEndCoord = coordsToKeep [ i ];
                    mpTractionElements.push_back(tractionEl);
                } else   {
                    TractionElement *tractionEl = new TractionElement();

                    tractionEl->mTractionNodeInd.push_back(mpTractionNodes.size() - 2);
                    tractionEl->mStartCoord = tractionNodeCoord [ mpTractionNodes.size() - 2 ];

                    tractionEl->mEndCoord = tractionNodeCoord [ mpTractionNodes.size() - 1 ];
                    mpTractionElements.push_back(tractionEl);
                }
            } else if ( mTractionInterpOrder == 1 )      {
                // Piecewise linear traction
                TractionElement *tractionEl = new TractionElement();
                tractionEl->mStartCoord = tractionNodeCoord [ mpTractionNodes.size() - 2 ];
                tractionEl->mTractionNodeInd.push_back(mpTractionNodes.size() - 2);

                tractionEl->mEndCoord = tractionNodeCoord [ mpTractionNodes.size() - 1 ];
                tractionEl->mTractionNodeInd.push_back(mpTractionNodes.size() - 1);

                if ( tractionEl->mStartCoord.distance(tractionEl->mEndCoord) > nodeDistTol ) {
                    mpTractionElements.push_back(tractionEl);
                } else   {
                    delete tractionEl;
                }
            }
        }
    }
#if 0
    printf("bndNodeCoords: \n");
    for(auto coordArray : bndNodeCoords) {
        for ( auto x :  coordArray ) {
            x.first.printYourself();
        }
    }
#endif


    // Construct maps necessary for assembly
    std :: vector< std :: pair< FloatArray, bool > >allBndNodeCoords;
    for ( auto &coordArray : bndNodeCoords ) {
        allBndNodeCoords.insert( allBndNodeCoords.end(), coordArray.begin(), coordArray.end() );
    }

    buildMaps(allBndNodeCoords);

    // Write traction nodes to debug vtk
    std :: vector< FloatArray >nodeCoord;
    for ( Node *node : mpTractionNodes ) {
        nodeCoord.push_back( * ( node->giveCoordinates() ) );
    }

    std :: string fileName("TractionNodeCoord.vtk");
    XFEMDebugTools :: WritePointsToVTK(fileName, nodeCoord);


    if ( mMeshIsPeriodic ) {
        // Lock displacement in one node if we use periodic BCs

        int numNodes = domain->giveNumberOfDofManagers();
        mpDisplacementLock = new Node(numNodes + 1, domain);
        mLockNodeInd = domain->giveElement(1)->giveNode(1)->giveGlobalNumber();


        int nsd = domain->giveNumberOfSpatialDimensions();
        std :: vector< int >dofIds;
        for ( int j = 0; j < nsd; j++ ) {
            dofIds.push_back( this->domain->giveNextFreeDofID() );
        }

        for ( auto &dofid: dofIds ) {
            mpDisplacementLock->appendDof( new MasterDof(mpDisplacementLock, ( DofIDItem ) dofid) );
        }
    }
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

void PrescribedGradientBCWeak :: buildMaps(const std :: vector< std :: pair< FloatArray, bool > > &iBndNodeCoordsFull)
{
    // Create map from a traction element to displacement elements it
    // interacts with everywhere on gamma.
    SpatialLocalizer *localizer = domain->giveSpatialLocalizer();
    mMapTractionElDispElGamma.clear();
    for ( size_t i = 0; i < mpTractionElements.size(); i++ ) {
        // Elements interacting on Gamma plus
        FloatArray xS_plus = mpTractionElements [ i ]->mStartCoord;
        FloatArray xE_plus = mpTractionElements [ i ]->mEndCoord;
        FloatArray xC_plus;
        xC_plus.beScaled(0.5, xS_plus);
        xC_plus.add(0.5, xE_plus);

        double elLength_plus = xS_plus.distance(xE_plus);
        std :: set< int >elList_plus;
        // TODO: What if an element is cut by two cracks, so that the
        // traction element becomes shorter than the displacement element?
        // Make sure that the search radius is never smaller than the
        // largest displacement element length along the boundary.
        localizer->giveAllElementsWithNodesWithinBox(elList_plus, xC_plus, 0.51 * elLength_plus);

        if ( elList_plus.empty() ) {
            FloatArray lCoords, closestPoint;
            Element *el = localizer->giveElementClosestToPoint(lCoords, closestPoint, xC_plus);
            int elPlaceInArray = domain->giveElementPlaceInArray( el->giveGlobalNumber() );
            elList_plus.insert(elPlaceInArray);
        }

        std :: vector< int >displacementElements, displacementElements_plus;
        for ( int elNum: elList_plus ) {
            // Check if the traction element and the displacement element intersect
            // Intersection occurs if at least one displacement element node is
            // on the traction element.
            Element *el = domain->giveElement(elNum);

            Line line_plus(xS_plus, xE_plus);
            if ( line_plus.intersects(el) ) {
                displacementElements.push_back(elNum);
                displacementElements_plus.push_back(elNum);
            }
        }

        if ( mMeshIsPeriodic ) {
            //mMapTractionElDispElGamma[i] = displacementElements;

            if ( xS_plus.distance(mUC) < 1.0e-12 ) {
                // Perturb in direction of xE
                FloatArray t;
                t.beDifferenceOf(xE_plus, xS_plus);
                xS_plus.add(1.0e-6, t);
                //printf("xS_plus: %.12e %.12e\n", xS_plus[0], xS_plus[1]);
            }

            if ( xE_plus.distance(mUC) < 1.0e-12 ) {
                // Perturb in direction of xS
                FloatArray t;
                t.beDifferenceOf(xS_plus, xE_plus);
                xE_plus.add(1.0e-6, t);
                //printf("xE_plus: %.12e %.12e\n", xE_plus[0], xE_plus[1]);
            }

            // Elements interacting on Gamma minus
            FloatArray xS_minus;
            giveMirroredPointOnGammaMinus(xS_minus, xS_plus);
            FloatArray xE_minus;
            giveMirroredPointOnGammaMinus(xE_minus, xE_plus);
            FloatArray xC_minus;
            xC_minus.beScaled(0.5, xS_minus);
            xC_minus.add(0.5, xE_minus);

            double elLength_minus = xS_minus.distance(xE_minus);
            std :: set< int >elList_minus;
            // TODO: What if an element is cut by two cracks, so that the
            // traction element becomes shorter than the displacement element?
            // Make sure that the search radius is never smaller than the
            // largest displacement element length along the boundary.
            localizer->giveAllElementsWithNodesWithinBox(elList_minus, xC_minus, 0.51 * elLength_minus);

            if ( elList_minus.empty() ) {
                FloatArray lCoords, closestPoint;
                Element *el = localizer->giveElementClosestToPoint(lCoords, closestPoint, xC_minus);
                int elPlaceInArray = domain->giveElementPlaceInArray( el->giveGlobalNumber() );
                elList_minus.insert(elPlaceInArray);
            }


            for ( int elNum: elList_minus ) {
                // Check if the traction element and the displacement element intersect
                // Intersection occurs if at least one displacement element node is
                // on the traction element.
                Element *el = domain->giveElement(elNum);

                Line line_minus(xS_minus, xE_minus);
                if ( line_minus.intersects(el) ) {
                    displacementElements.push_back(elNum);
                }
            }
        }



        mMapTractionElDispElGamma [ i ] = displacementElements;
    }


    // Create map from a traction element to displacement element and crack
    // intersection coordinates in the traction element.
    mTractionElInteriorCoordinates.clear();
    const double distTol = 1.0e-12;

    for ( size_t i = 0; i < mpTractionElements.size(); i++ ) {
        const FloatArray &xS = mpTractionElements [ i ]->mStartCoord;
        const FloatArray &xE = mpTractionElements [ i ]->mEndCoord;

        PolygonLine pl;
        pl.insertVertexBack(xS);
        pl.insertVertexBack(xE);

        for ( auto x :  iBndNodeCoordsFull ) {
            double distN;
            pl.computeNormalSignDist(distN, x.first);

            double distT, arcPos;
            pl.computeTangentialSignDist(distT, x.first, arcPos);

            if ( fabs(distN) < distTol && distT > -distTol && distT < ( 1.0 + distTol ) ) {
                mTractionElInteriorCoordinates [ i ].push_back(x.first);
            }
        }

        // Sort based on distance to start point
        std :: sort( mTractionElInteriorCoordinates [ i ].begin(), mTractionElInteriorCoordinates [ i ].end(), ArcPosSortFunction(xS) );
    }


    // Create map from a nodes global number to its local index
    // in the traction element.
    mTracElDispNodes.clear();
    for ( size_t i = 0; i < mpTractionElements.size(); i++ ) {
        std :: vector< int >tracElDispNodes;

        for ( int elIndex : mMapTractionElDispElGamma [ i ] ) {
            const IntArray &elNodes = domain->giveElement(elIndex)->giveDofManArray();

            for ( int nodeInd : elNodes ) {
                tracElDispNodes.push_back(nodeInd);
            }
        }


        // We now have an array of displacement nodes that affect the traction element.
        // It may contain duplicates and it is not sorted in any particular order.

        // Sort the array
        std :: sort( tracElDispNodes.begin(), tracElDispNodes.end() );

        // Remove duplicates
        tracElDispNodes.erase( std :: unique( tracElDispNodes.begin(), tracElDispNodes.end() ), tracElDispNodes.end() );

        mTracElDispNodes [ i ] = tracElDispNodes;
    }
}

void PrescribedGradientBCWeak :: integrateTangent(FloatMatrix &oTangent, size_t iTracElInd)
{
    // Compute the tangent stiffness contribution
    // K = int (N^trac)^T . N^disp dGamma

    int dim = domain->giveNumberOfSpatialDimensions();

    IntegrationRule *ir = createNewIntegrationRule(iTracElInd);

    // Number of rows and colums
    int numRows = 0;
    if ( mTractionInterpOrder == 0 ) {
        numRows = dim;
    } else if ( mTractionInterpOrder == 1 ) {
        numRows = 2 * dim;
    }

    int numCols = 0;
    std :: unordered_map< int, IntArray >globalNodeIndToPosInLocalLocArray;
    for ( auto nodeInd : mTracElDispNodes [ iTracElInd ] ) {
        DofManager *dMan = domain->giveDofManager(nodeInd);

        for ( int i = 0; i < dMan->giveNumberOfDofs(); i++ ) {
            globalNodeIndToPosInLocalLocArray [ nodeInd ].followedBy(numCols + i + 1);
        }

        numCols += dMan->giveNumberOfDofs();
    }

    oTangent.resize(numRows, numCols);


    for ( GaussPoint *gp: *ir ) {
        /*
         * Contribution from Gamma if Dirichlet,
         * or Gamma_plus is Periodic
         */

        // Fetch GP coordinates
        const FloatArray &globalCoord = gp->giveGlobalCoordinates();

        assembleTangentGPContribution(oTangent, iTracElInd, * gp, globalCoord, globalNodeIndToPosInLocalLocArray, 1.0);

        /*
         * Contribution from Gamma_minus if periodic.
         */
        if ( mMeshIsPeriodic ) {
            FloatArray globalCoord_minus;
            giveMirroredPointOnGammaMinus(globalCoord_minus, globalCoord);
            assembleTangentGPContribution(oTangent, iTracElInd, * gp, globalCoord_minus, globalNodeIndToPosInLocalLocArray, -1.0);
        }
    }

    delete ir;
}

void PrescribedGradientBCWeak :: assembleTangentGPContribution(FloatMatrix &oTangent, size_t iTracElInd, GaussPoint &iGP, const FloatArray &iBndCoord, std :: unordered_map< int, IntArray > &iGlobalNodeIndToPosInLocalLocArray, const double &iScaleFactor)
{
    int dim = domain->giveNumberOfSpatialDimensions();

    const TractionElement &tEl = * ( mpTractionElements [ iTracElInd ] );
    double detJ = 0.5 * tEl.mStartCoord.distance(tEl.mEndCoord);
    const FloatArray &locCoordsOnLine = iGP.giveNaturalCoordinates();

    //////////////////////////////////
    // Compute traction N-matrix
    FloatArray N, Ntrac;
    computeNTraction(Ntrac, locCoordsOnLine [ 0 ], tEl);
    tEl.computeN_Linear(N, locCoordsOnLine [ 0 ]);


    FloatMatrix NtracMat;
    NtracMat.beNMatrixOf(Ntrac, dim);


    //////////////////////////////////
    // Compute displacement N-matrix

    // Identify the displacement element
    // we are currently standing in
    // and compute local coordinates on
    // the displacement element
    SpatialLocalizer *localizer = domain->giveSpatialLocalizer();
    FloatArray dispElLocCoord_minus, closestPoint_minus;
    Element *dispEl_minus = localizer->giveElementClosestToPoint(dispElLocCoord_minus, closestPoint_minus, iBndCoord);

    // Compute basis functions
    XfemElementInterface *xfemElInt_minus = dynamic_cast< XfemElementInterface * >( dispEl_minus );
    FloatMatrix NdispMat_minus;

    if ( xfemElInt_minus != NULL && domain->hasXfemManager() ) {
        // If the element is an XFEM element, we use the XfemElementInterface to compute the N-matrix
        // of the enriched element.
        xfemElInt_minus->XfemElementInterface_createEnrNmatrixAt(NdispMat_minus, dispElLocCoord_minus, * dispEl_minus, false);
    } else   {
        // Otherwise, use the usual N-matrix.
        const int numNodes = dispEl_minus->giveNumberOfDofManagers();
        FloatArray N(numNodes);

        const int dim = dispEl_minus->giveSpatialDimension();

        NdispMat_minus.resize(dim, dim * numNodes);
        NdispMat_minus.zero();
        dispEl_minus->giveInterpolation()->evalN( N, dispElLocCoord_minus, FEIElementGeometryWrapper(dispEl_minus) );

        NdispMat_minus.beNMatrixOf(N, dim);

    }

    FloatMatrix contrib;
    contrib.beTProductOf(NtracMat, NdispMat_minus);
    contrib.times( iScaleFactor * detJ * iGP.giveWeight() );


    // Create local location arrays
    IntArray rows, cols;
    for ( int i = 1; i <= ( mTractionInterpOrder + 1 ) * dim; i++ ) {
        rows.followedBy(i);
    }

    const IntArray &dispElNodes_minus = dispEl_minus->giveDofManArray();
    for ( int nodeInd : dispElNodes_minus ) {
        cols.followedBy(iGlobalNodeIndToPosInLocalLocArray [ nodeInd ]);
    }

    if(contrib.giveNumberOfRows() == rows.giveSize() && contrib.giveNumberOfColumns() == cols.giveSize()) {
        oTangent.assemble(contrib, rows, cols);
    }
    else {
        contrib.resize(rows.giveSize(), cols.giveSize());
        contrib.zero();
        oTangent.assemble(contrib, rows, cols);
        printf("Warning in PrescribedGradientBCWeak :: assembleTangentGPContribution: rows.giveSize(): %d cols.giveSize(): %d\n", rows.giveSize(), cols.giveSize() );
    }
}

IntegrationRule *PrescribedGradientBCWeak :: createNewIntegrationRule(int iTracElInd)
{
    std :: vector< FloatArray >tracGpCoord;

    // Create integration rule
    // -> need segments defined by displacement nodes, traction nodes and cracks.
    int numSeg = int( mTractionElInteriorCoordinates [ iTracElInd ].size() ) - 1;
    std :: vector< Line >segments;


    for ( int segIndex = 0; segIndex < numSeg; segIndex++ ) {
        const FloatArray &xS = mTractionElInteriorCoordinates [ iTracElInd ] [ segIndex ];
        const FloatArray &xE = mTractionElInteriorCoordinates [ iTracElInd ] [ segIndex + 1 ];

        if ( xS.distance(xE) > 1.0e-12 ) {
            segments.push_back( Line(xS, xE) );
        }
    }

    Element *dummyEl = NULL;
    const FloatArray &xS = mTractionElInteriorCoordinates [ iTracElInd ] [ 0 ];
    const FloatArray &xE = mTractionElInteriorCoordinates [ iTracElInd ].back();

    IntegrationRule *ir = new DiscontinuousSegmentIntegrationRule(1, dummyEl, segments, xS, xE);

    // Take material mode from first element
    MaterialMode matMode = domain->giveElement(1)->giveMaterialMode();

    int numPointsPerSeg = 3;
    ir->SetUpPointsOnLine(numPointsPerSeg, matMode);

    for ( GaussPoint *gp : *ir ) {
        tracGpCoord.push_back( gp->giveGlobalCoordinates() );
    }

    std :: stringstream str3;
    str3 << "tracGpCoordEl" << iTracElInd << ".vtk";
    std :: string name3 = str3.str();

    XFEMDebugTools :: WritePointsToVTK(name3, tracGpCoord);

    return ir;
}

void PrescribedGradientBCWeak :: computeNTraction(FloatArray &oN, const double &iXi, const TractionElement &iEl) const
{
    if ( mTractionInterpOrder == 0 ) {
        iEl.computeN_Constant(oN, iXi);
    } else if ( mTractionInterpOrder == 1 ) {
        iEl.computeN_Linear(oN, iXi);
    }

    oN.times(mTracDofScaling);
}

void PrescribedGradientBCWeak :: giveTractionUnknows(FloatArray &oTracUnknowns, ValueModeType mode, TimeStep *tStep, int iTracElInd)
{
    oTracUnknowns.clear();

    const std :: vector< int > &tNodeInd = mpTractionElements [ iTracElInd ]->mTractionNodeInd;

    for ( int nodeInd : tNodeInd ) {
        FloatArray nodeUnknowns;
        mpTractionNodes [ nodeInd ]->giveCompleteUnknownVector(nodeUnknowns, mode, tStep);
        oTracUnknowns.append(nodeUnknowns);
    }
}

void PrescribedGradientBCWeak :: giveDisplacementUnknows(FloatArray &oDispUnknowns, ValueModeType mode, TimeStep *tStep, int iTracElInd)
{
    oDispUnknowns.clear();

    for ( int nodeInd : mTracElDispNodes [ iTracElInd ] ) {
        FloatArray nodeUnknowns;
        domain->giveDofManager(nodeInd)->giveCompleteUnknownVector(nodeUnknowns, mode, tStep);
        oDispUnknowns.append(nodeUnknowns);
    }
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
    oPosMinus = iPosPlus;
    const double distTol = 1.0e-12;

    if ( iPosPlus.distance(mUC) < distTol ) {
        printf("iPosPlus: %.12e %.12e\n", iPosPlus [ 0 ], iPosPlus [ 1 ]);
        OOFEM_ERROR("Unmappable point.")
    }

    double mappingPerformed = false;

    if ( iPosPlus [ 0 ] > mUC [ 0 ] - distTol ) {
        oPosMinus [ 0 ] = mLC [ 0 ];
        mappingPerformed = true;
    }

    if ( iPosPlus [ 1 ] > mUC [ 1 ] - distTol ) {
        oPosMinus [ 1 ] = mLC [ 1 ];
        mappingPerformed = true;
    }

    if ( !mappingPerformed ) {
        iPosPlus.printYourself();
        OOFEM_ERROR("Mapping failed.")
    }

    //    printf("iPosPlus: "); iPosPlus.printYourself();
    //    printf("oPosMinus: "); oPosMinus.printYourself();
}

void PrescribedGradientBCWeak :: giveMirroredPointOnGammaPlus(FloatArray &oPosPlus, const FloatArray &iPosMinus) const
{
    oPosPlus = iPosMinus;
    const double distTol = 1.0e-16;

    if ( iPosMinus.distance(mLC) < distTol ) {
        printf("iPosMinus: %.12e %.12e\n", iPosMinus [ 0 ], iPosMinus [ 1 ]);
        OOFEM_ERROR("Unmappable point.")
    }

    double mappingPerformed = false;

    if ( iPosMinus [ 0 ] < mLC [ 0 ] + distTol ) {
        oPosPlus [ 0 ] = mUC [ 0 ];
        mappingPerformed = true;
    }

    if ( iPosMinus [ 1 ] < mLC [ 1 ] + distTol ) {
        oPosPlus [ 1 ] = mUC [ 1 ];
        mappingPerformed = true;
    }

    if ( !mappingPerformed ) {
        iPosMinus.printYourself();
        OOFEM_ERROR("Mapping failed.")
    }
}

bool PrescribedGradientBCWeak :: pointIsMapapble(const FloatArray &iPos) const
{
    const double distTol = 1.0e-13;
    return !( iPos.distance(mLC) < distTol );
}

void PrescribedGradientBCWeak :: computeDomainBoundingBox(Domain &iDomain, FloatArray &oLC, FloatArray &oUC)
{
    // Compute LC and UC by assuming a rectangular domain.
    int numNodes = iDomain.giveNumberOfDofManagers();
    int nsd = iDomain.giveNumberOfSpatialDimensions();

    oLC = * ( iDomain.giveDofManager(1)->giveCoordinates() );
    oUC = * ( iDomain.giveDofManager(1)->giveCoordinates() );

    for ( int i = 1; i <= numNodes; i++ ) {
        DofManager *dMan = iDomain.giveDofManager(i);
        const FloatArray &coord = * ( dMan->giveCoordinates() );
        bool nodeIsLC = true;
        bool nodeIsUC = true;
        for ( int j = 0; j < nsd; j++ ) {
            if ( coord [ j ] > oLC [ j ] ) {
                nodeIsLC = false;
            }

            if ( coord [ j ] < oUC [ j ] ) {
                nodeIsUC = false;
            }
        }

        if ( nodeIsLC ) {
            oLC = coord;
        }

        if ( nodeIsUC ) {
            oUC = coord;
        }
    }
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

bool PrescribedGradientBCWeak :: closePointExists(const std :: vector< std :: pair< FloatArray, bool > > &iCoordArray, const FloatArray &iPos, const double &iMeshTol2) const
{
    for ( auto bndPos : iCoordArray ) {
        if ( bndPos.first.distance_square(iPos) < iMeshTol2 ) {
            return true;
        }
    }
    return false;
}
} /* namespace oofem */
