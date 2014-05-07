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

#include "xfem/XFEMDebugTools.h"

#include <cmath>
#include <string>
#include <sstream>

namespace oofem {
REGISTER_BoundaryCondition(PrescribedGradientBCWeak);

PrescribedGradientBCWeak::PrescribedGradientBCWeak(int n, Domain * d):
PrescribedGradientBC(n, d),
mTractionInterpOrder(0),
mNumTractionNodesAtIntersections(1),
mTractionNodeSpacing(1),
mTractionLivesOnGammaPlus(false),
mDuplicateCornerNodes(false)
{

    // Compute LC and UC by assuming a rectangular domain.
    int numNodes = d->giveNumberOfDofManagers();
    int nsd = d->giveNumberOfSpatialDimensions();

    mLC = *(d->giveDofManager(1)->giveCoordinates());
    mUC = *(d->giveDofManager(1)->giveCoordinates());

    for(int i = 1; i <= numNodes; i++) {
        DofManager *dMan = d->giveDofManager(i);
        const FloatArray &coord = *(dMan->giveCoordinates());
        bool nodeIsLC = true;
        bool nodeIsUC = true;
        for(int j = 0; j < nsd; j++) {

            if( coord[j] > mLC[j] ) {
                nodeIsLC = false;
            }

            if( coord[j] < mUC[j] ) {
                nodeIsUC = false;
            }

        }

        if(nodeIsLC) {
            mLC = coord;
        }

        if(nodeIsUC) {
            mUC = coord;
        }

    }
}

PrescribedGradientBCWeak::~PrescribedGradientBCWeak() {

    for(size_t i = 0; i < mpTractionNodes.size(); i++) {
        if(mpTractionNodes[i] != NULL) {
            delete mpTractionNodes[i];
            mpTractionNodes[i] = NULL;
        }
    }
    mpTractionNodes.clear();

    for(size_t i = 0; i < mpTractionElements.size(); i++) {
        if(mpTractionElements[i] != NULL) {
            delete mpTractionElements[i];
            mpTractionElements[i] = NULL;
        }
    }
    mpTractionElements.clear();
}

DofManager* PrescribedGradientBCWeak::giveInternalDofManager(int i)
{
    return mpTractionNodes[i-1];
}

IRResultType PrescribedGradientBCWeak :: initializeFrom(InputRecord *ir)
{
    IRResultType result;
    PrescribedGradientBC :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, mTractionInterpOrder, _IFT_PrescribedGradientBCWeak_TractionInterpOrder);
    printf("mTractionInterpOrder: %d\n", mTractionInterpOrder);

    IR_GIVE_FIELD(ir, mNumTractionNodesAtIntersections, _IFT_PrescribedGradientBCWeak_NumTractionNodesAtIntersections);
    printf("mNumTractionNodesAtIntersections: %d\n", mNumTractionNodesAtIntersections);

    if( mNumTractionNodesAtIntersections > 1 && mTractionInterpOrder == 0 ) {
        OOFEM_ERROR("mNumTractionNodesAtIntersections > 1 is not allowed if mTractionInterpOrder == 0.")
    }

    IR_GIVE_FIELD(ir, mTractionNodeSpacing, _IFT_PrescribedGradientBCWeak_NumTractionNodeSpacing);
    printf("mTractionNodeSpacing: %d\n", mTractionNodeSpacing);

    int periodic = 0;
    IR_GIVE_FIELD(ir, periodic, _IFT_PrescribedGradientBCWeak_TractionOnGammaPlus);
    printf("periodic: %d\n", periodic);

    if(periodic == 1) {
        mTractionLivesOnGammaPlus = true;
    }
    else {
        mTractionLivesOnGammaPlus = false;
    }

    int duplicateCorners = 0;
    IR_GIVE_FIELD(ir, duplicateCorners, _IFT_PrescribedGradientBCWeak_DuplicateCornerNodes);
    printf("duplicateCorners: %d\n", duplicateCorners);

    if(duplicateCorners == 1) {
        mDuplicateCornerNodes = true;
    }
    else {
        mDuplicateCornerNodes = false;
    }

    return IRRT_OK;
}

void PrescribedGradientBCWeak :: giveInputRecord(DynamicInputRecord &input)
{

}

void PrescribedGradientBCWeak::postInitialize()
{
    createTractionMesh();
}


void PrescribedGradientBCWeak::scale(double s)
{

}

void PrescribedGradientBCWeak::assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                            CharType type, ValueModeType mode,
                            const UnknownNumberingScheme &s, FloatArray *eNorm)
{
    int dim = domain->giveNumberOfSpatialDimensions();

    if ( type == ExternalForcesVector ) {
        // The external force vector is given by
        // f_ext = int N^trac H . (x - x_c)

        for(size_t i = 0; i < mpTractionElements.size(); i++) {

            // Compute f_ext
            FloatArray fExt;

            if( mTractionInterpOrder == 0) {
                fExt.resize(1*dim);
            }
            else if( mTractionInterpOrder == 1 ) {
                fExt.resize(2*dim);
            }
            fExt.zero();


            const TractionElement *el = mpTractionElements[i];
            IntegrationRule *ir = createNewIntegrationRule(i);

            for ( GaussPoint *gp: *ir ) {
                const FloatArray &locCoordsOnLine = * gp->giveLocalCoordinates();

                // Compute N^trac
                FloatArray N, Ntrac;
                computeNTraction(Ntrac, locCoordsOnLine[0], *el);
                el->computeN_Linear(N, locCoordsOnLine[0]);

                // Compute x
                FloatArray x(el->mStartCoord.giveSize());
                x.zero();

                x.add( N[0], el->mStartCoord );
                x.add( N[1], el->mEndCoord );

                // Compute H.(x - x_c)
                FloatArray temp;

                if(mGradient.giveNumberOfRows() == 2) {
                    temp = {x[0]-mCenterCoord[0], x[1]-mCenterCoord[1]};
                }
                else {
                    temp.beDifferenceOf(x, mCenterCoord);
                }

                FloatArray Hx;
                Hx.beProductOf(mGradient, temp);


                // N-matrix
                FloatMatrix Nmat;
                Nmat.beNMatrixOf(Ntrac, dim);


                // Add contribution to fExt
                FloatArray contrib;
                contrib.beTProductOf(Nmat, Hx);
                double detJ = 0.5*el->mStartCoord.distance(el->mEndCoord);
                fExt.add( detJ*gp->giveWeight(), contrib );
            }

            // Fetch location arrays for current traction element
            IntArray rows;
            giveTractionLocationArrays( i, rows, eid, type, s);

            fExt.negated();

            // Assemble
            answer.assemble(fExt, rows);

            delete ir;
        }


    } else if ( type == InternalForcesVector ) {

        FloatArray fe_trac, fe_disp;

        for(size_t i = 0; i < mpTractionElements.size(); i++) {

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

            fe_trac.beProductOf(Ke, dispUnknowns);
            fe_trac.negated();


            // Fetch location arrays
            IntArray tracRrows;
            giveTractionLocationArrays( i, tracRrows, eid, type, s);

            IntArray dispRrows;
            giveDisplacementLocationArrays( i, dispRrows, eid, type, s);


            // Assemble
            answer.assemble(fe_trac, tracRrows);
            answer.assemble(fe_disp, dispRrows);

        }
    }


}


void PrescribedGradientBCWeak::assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
                      CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    if ( eid == EID_MomentumBalance_ConservationEquation ) {
        // TODO: Check if this applies here
        eid = EID_MomentumBalance;
    }

    if ( eid != EID_MomentumBalance ) {
        printf("returning.");
        return;
    }

    if ( type == TangentStiffnessMatrix || type == SecantStiffnessMatrix || type == StiffnessMatrix || type == ElasticStiffnessMatrix ) {

        FloatMatrix Ke, KeT;
        IntArray tracRows, tracCols, dispRows, dispCols;

        for(size_t i = 0; i < mpTractionElements.size(); i++) {

            // Rows and columns for displacement and traction contributions
            giveTractionLocationArrays( i, tracRows, eid, type, r_s);
            giveDisplacementLocationArrays( i, dispRows, eid, type, r_s);

            giveTractionLocationArrays( i, tracCols, eid, type, c_s);
            giveDisplacementLocationArrays( i, dispCols, eid, type, c_s);

            this->integrateTangent(Ke, i);
            Ke.negated();
            KeT.beTranspositionOf(Ke);

            answer->assemble(tracRows, dispCols, Ke);
            answer->assemble(dispRows, tracCols, KeT);
        }
    }
    else {
        printf("Skipping assembly in PrescribedGradientBCWeak::assemble().\n");
    }

}

void PrescribedGradientBCWeak::giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, EquationID eid, CharType type,
                                const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    rows.clear();
    cols.clear();

    if(!mTractionLivesOnGammaPlus) {
        // Loop over traction elements
        for(size_t tracElInd = 0; tracElInd < mpTractionElements.size(); tracElInd++) {

            IntArray tracElCols, tracElRows, trac_loc_c, trac_loc_r;

            const TractionElement &tEl = *(mpTractionElements[tracElInd]);
            for(int tracNodeInd : tEl.mTractionNodeInd) {
                Node *tNode = mpTractionNodes[tracNodeInd];
                tNode->giveCompleteLocationArray(trac_loc_r, r_s);
                tracElRows.followedBy(trac_loc_r);

                tNode->giveCompleteLocationArray(trac_loc_c, c_s);
                tracElCols.followedBy(trac_loc_c);
            }


            // Fetch displacement elements that interact
            // with the current traction element.
            const std::vector<int> &dispElInds = mMapTractionElDispElGamma[tracElInd];

            IntArray dispElCols, dispElRows, disp_loc_c, disp_loc_r;

            for( int dispElInd : dispElInds ) {
                Element *e = this->giveDomain()->giveElement( dispElInd );
                e->giveLocationArray(disp_loc_r, eid, r_s);
                dispElRows.followedBy(disp_loc_r);

                e->giveLocationArray(disp_loc_c, eid, c_s);
                dispElCols.followedBy(disp_loc_c);
            }

            rows.push_back(dispElRows);
            cols.push_back(tracElCols);

            rows.push_back(tracElRows);
            cols.push_back(dispElCols);
        }
    }
    else {
        OOFEM_ERROR("Not implemented.");
    }

}

void PrescribedGradientBCWeak::giveLocationArrays(int iTracElInd, std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, EquationID eid, CharType type,
                                const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    rows.clear();
    cols.clear();

    if(!mTractionLivesOnGammaPlus) {

        IntArray tracElCols, tracElRows, trac_loc_c, trac_loc_r;

        const TractionElement &tEl = *(mpTractionElements[iTracElInd]);
        for(int tracNodeInd : tEl.mTractionNodeInd) {
            Node *tNode = mpTractionNodes[tracNodeInd];
            tNode->giveCompleteLocationArray(trac_loc_r, r_s);
            tracElRows.followedBy(trac_loc_r);

            tNode->giveCompleteLocationArray(trac_loc_c, c_s);
            tracElCols.followedBy(trac_loc_c);
        }


        // Fetch displacement elements that interact
        // with the current traction element.
        const std::vector<int> &dispElInds = mMapTractionElDispElGamma[iTracElInd];

        IntArray dispElCols, dispElRows, disp_loc_c, disp_loc_r;

        for( int dispElInd : dispElInds ) {
            Element *e = this->giveDomain()->giveElement( dispElInd );
            e->giveLocationArray(disp_loc_r, eid, r_s);
            dispElRows.followedBy(disp_loc_r);

            e->giveLocationArray(disp_loc_c, eid, c_s);
            dispElCols.followedBy(disp_loc_c);
        }

        rows.push_back(dispElRows);
        cols.push_back(tracElCols);

        rows.push_back(tracElRows);
        cols.push_back(dispElCols);
    }
    else {
        OOFEM_ERROR("Not implemented.");
    }
}

void PrescribedGradientBCWeak::giveTractionLocationArrays(int iTracElInd, IntArray &rows, EquationID eid, CharType type,
                                const UnknownNumberingScheme &s)
{
    rows.clear();

    if(!mTractionLivesOnGammaPlus) {

        IntArray tracElRows, trac_loc_c, trac_loc_r;

        const TractionElement &tEl = *(mpTractionElements[iTracElInd]);
        for(int tracNodeInd : tEl.mTractionNodeInd) {
            Node *tNode = mpTractionNodes[tracNodeInd];
            tNode->giveCompleteLocationArray(trac_loc_r, s);
            tracElRows.followedBy(trac_loc_r);
        }

        rows = tracElRows;
    }
    else {
        OOFEM_ERROR("Not implemented.");
    }
}

void PrescribedGradientBCWeak::giveDisplacementLocationArrays(int iTracElInd, IntArray &rows, EquationID eid, CharType type,
                                const UnknownNumberingScheme &s)
{
    rows.clear();

    for(int nodeInd : mTracElDispNodes[iTracElInd]) {
        IntArray nodeLocationArray;
        domain->giveDofManager(nodeInd)->giveCompleteLocationArray(nodeLocationArray, s);
        rows.followedBy(nodeLocationArray);
    }

}

void PrescribedGradientBCWeak::computeField(FloatArray &sigma, EquationID eid, TimeStep *tStep)
{
    double dSize = domainSize();

    const int dim = domain->giveNumberOfSpatialDimensions();
    FloatMatrix stressMatrix(dim,dim);


    for(size_t i = 0; i < mpTractionElements.size(); i++) {

        const TractionElement *el = mpTractionElements[i];
        IntegrationRule *ir = createNewIntegrationRule(i);

        for ( GaussPoint *gp: *ir ) {
            const FloatArray &locCoordsOnLine = * gp->giveLocalCoordinates();

            // Compute N^trac
            FloatArray N, Ntrac;
            computeNTraction(Ntrac, locCoordsOnLine[0], *el);
            el->computeN_Linear(N, locCoordsOnLine[0]);

            // Compute x
            FloatArray x(el->mStartCoord.giveSize());
            x.zero();

            x.add( N[0], el->mStartCoord );
            x.add( N[1], el->mEndCoord );


            // N-matrix
            FloatMatrix Nmat;
            Nmat.beNMatrixOf(Ntrac, dim);

            // Interpolate traction
            FloatArray tracUnknowns;
            giveTractionUnknows(tracUnknowns, VM_Total, tStep, i);

            FloatArray traction;
            traction.beProductOf(Nmat, tracUnknowns);
            //printf("traction: "); traction.printYourself();

            FloatArray tmp;
            tmp.beDifferenceOf(x, mCenterCoord);

            FloatMatrix contrib;
            contrib.beDyadicProductOf(traction, tmp);

            double detJ = 0.5*el->mStartCoord.distance(el->mEndCoord);
            contrib.times( detJ * gp->giveWeight() );

            for(int m = 0; m < dim; m++) {
                for(int n = 0; n < dim; n++) {
                    stressMatrix(m,n) += contrib(m,n);
                }
            }
        }


        delete ir;
    }


    if(dim == 2) {
        sigma = { stressMatrix(0,0), stressMatrix(1,1), stressMatrix(0,1), stressMatrix(1,0)};
    }
    else {
        sigma.beVectorForm(stressMatrix);
    }

    sigma.times(1.0/dSize);

#if 0
    FloatArray tx,ty;
    for( Node *node : mpTractionNodes ) {
        FloatArray tNode;
        node->giveCompleteUnknownVector(tNode, VM_Total, tStep);
        tx.push_back(tNode[0]);
        ty.push_back(tNode[1]);
    }

    printf("\n\n\n");
    printf("tx: "); tx.printYourself();
    printf("\n\n\n");
    printf("ty: "); ty.printYourself();
#endif
}

void PrescribedGradientBCWeak::createTractionMesh()
{
    const double nodeDistTol = 1.0e-15;

    std::vector<FloatArray> bndNodeCoords, bndNodeCoordsFull;
    std::vector<bool> mandatoryToKeep;

    if(!mTractionLivesOnGammaPlus) {
        //printf("Creating Dirichlet tractions.\n");

        Set *setPointer = this->giveDomain()->giveSet(this->set);
        const IntArray &boundaries = setPointer->giveBoundaryList();
        IntArray bNodes;

        bool duplicateAtCorner = mDuplicateCornerNodes;

        // Loop over all boundary segments
        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
            int boundary = boundaries.at(pos * 2);

            e->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);

            // Add the start node of the segment
            DofManager *startNode = e->giveDofManager( bNodes[0] );
            bndNodeCoords.push_back( *(startNode->giveCoordinates()) );
            bndNodeCoordsFull.push_back( *(startNode->giveCoordinates()) );

            bool isCorner = false;
#if 1
            FloatArray cornerPos = mLC;
            if( startNode->giveCoordinates()->distance(cornerPos) < nodeDistTol ) {
                isCorner = true;
            }

            cornerPos = {mUC[0], mLC[1]};
            if( startNode->giveCoordinates()->distance( cornerPos ) < nodeDistTol ) {
                isCorner = true;
            }

            cornerPos = {mUC[0], mUC[1]};
            if( startNode->giveCoordinates()->distance( cornerPos ) < nodeDistTol ) {
                isCorner = true;
            }

            cornerPos = {mLC[0], mUC[1]};
            if( startNode->giveCoordinates()->distance( cornerPos ) < nodeDistTol ) {
                isCorner = true;
            }
#endif


            if(isCorner) {

                if(duplicateAtCorner) {
                    printf("Found corner: "); startNode->giveCoordinates()->printYourself();
                    printf("Duplicating...\n");
                    bndNodeCoords.push_back( *(startNode->giveCoordinates()) );
                    bndNodeCoordsFull.push_back( *(startNode->giveCoordinates()) );
                    mandatoryToKeep.push_back(true);
                }
            }

            if(isCorner) {
                mandatoryToKeep.push_back(true);
            }
            else {
                mandatoryToKeep.push_back(false);
            }

            // Add traction nodes where cracks intersect the boundary if desired
            XfemElementInterface *xfemElInt = dynamic_cast<XfemElementInterface*> (e);

            if(xfemElInt != NULL && domain->hasXfemManager() ) {
                std::vector<Line> segments;
                std::vector<FloatArray> intersecPoints;
                xfemElInt->partitionEdgeSegment( boundary, segments, intersecPoints );

                for(size_t i = 0; i < intersecPoints.size(); i++) {
                    for(int j = 0; j < mNumTractionNodesAtIntersections; j++) {
                        bndNodeCoords.push_back( intersecPoints[i] );
                        mandatoryToKeep.push_back(true);
                    }

                    bndNodeCoordsFull.push_back( intersecPoints[i] );

                }
            }

        }

        std::vector<FloatArray> bndNodeCoordsToKeep;
        int numPointsPassed = 0;
        for(size_t i = 0; i < bndNodeCoords.size(); i++) {
            numPointsPassed++;
            if( mandatoryToKeep[i] || numPointsPassed >= mTractionNodeSpacing ) {
                bndNodeCoordsToKeep.push_back(bndNodeCoords[i]);
                //printf("Keeping node %lu located at: ", i); bndNodeCoords[i].printYourself();
                numPointsPassed = 0;
            }

        }

        bndNodeCoords = bndNodeCoordsToKeep;


#if 0
        printf("bndNodeCoords: \n");
        for(FloatArray x :  bndNodeCoords) {
            x.printYourself();
        }
#endif

        // Create traction dofs
        int nsd = domain->giveNumberOfSpatialDimensions();
        std::vector<int> dofIds;
        for(int j = 0; j < nsd; j++) {
            dofIds.push_back( this->domain->giveNextFreeDofID() );
        }
        for(size_t i = 0; i < bndNodeCoords.size(); i++) {

            int numNodes = domain->giveNumberOfDofManagers();
            Node *node = new Node(numNodes+1, domain);
            for(int j = 0; j < nsd; j++) {
                node->appendDof( new MasterDof( j + 1, node, ( DofIDItem ) ( dofIds[j] ) ) );
            }

            node->setCoordinates(bndNodeCoords[i]);
            mpTractionNodes.push_back(node);
        }

        // Create traction elements
        if(mTractionInterpOrder == 0) {
            for(size_t i = 0; i < mpTractionNodes.size(); i++) {

                TractionElement *tractionEl = new TractionElement();

                tractionEl->mTractionNodeInd.push_back( i );
                tractionEl->mStartCoord = *(mpTractionNodes[i]->giveCoordinates());

                if( i < mpTractionNodes.size()-1) {
                    tractionEl->mEndCoord = *(mpTractionNodes[i+1]->giveCoordinates());
                }
                else {
                    tractionEl->mEndCoord = *(mpTractionNodes[0]->giveCoordinates());
                }

                mpTractionElements.push_back(tractionEl);
            }
        }
        else if(mTractionInterpOrder == 1) {
            for(size_t i = 0; i < mpTractionNodes.size(); i++) {

                TractionElement *tractionEl = new TractionElement();

                if( i < mpTractionNodes.size()-1) {
                    tractionEl->mStartCoord = *(mpTractionNodes[i]->giveCoordinates());
                    tractionEl->mTractionNodeInd.push_back( i );

                    tractionEl->mEndCoord = *(mpTractionNodes[i+1]->giveCoordinates());
                    tractionEl->mTractionNodeInd.push_back( i+1 );
                }
                else {
#if 1
                    tractionEl->mStartCoord = *(mpTractionNodes[0]->giveCoordinates());
                    tractionEl->mTractionNodeInd.push_back( 0 );

                    tractionEl->mEndCoord = *(mpTractionNodes[i]->giveCoordinates());
                    tractionEl->mTractionNodeInd.push_back( i );
#else
                    tractionEl->mStartCoord = *(mpTractionNodes[i]->giveCoordinates());
                    tractionEl->mTractionNodeInd.push_back( i );

                    tractionEl->mEndCoord = *(mpTractionNodes[0]->giveCoordinates());
                    tractionEl->mTractionNodeInd.push_back( 0 );
#endif
                }


                // Take care of duplicated nodes at crack intersections
                if( tractionEl->mStartCoord.distance( tractionEl->mEndCoord ) > nodeDistTol ) {
                    mpTractionElements.push_back(tractionEl);
                }
                else {
                    printf("Skipping traction el with start coord: "); tractionEl->mStartCoord.printYourself();
                    delete tractionEl;
                }

            }
        }


        // Create map from a traction element to displacement elements it
        // interacts with everywhere on gamma.
        SpatialLocalizer *localizer = domain->giveSpatialLocalizer();
        mMapTractionElDispElGamma.clear();
        for(size_t i = 0; i < mpTractionElements.size(); i++) {
            const FloatArray &xS = mpTractionElements[i]->mStartCoord;
            const FloatArray &xE = mpTractionElements[i]->mEndCoord;
            FloatArray xC;
            xC.beScaled(0.5, xS);
            xC.add(0.5, xE);

            double elLength = xS.distance(xE);
            std :: set< int >elList;
            // TODO: What if an element is cut by two cracks, so that the
            // traction element becomes shorter than the displacement element?
            // Make sure that the search radius is never smaller than the
            // largest displacement element length along the boundary.
            localizer->giveAllElementsWithNodesWithinBox(elList, xC, 0.51*elLength );

            std :: vector< int > displacementElements;
            for ( int elNum: elList ) {

                // Check if the traction element and the displacement element intersect
                // Intersection occurs if at least one displacement element node is
                // on the traction element.
                Element *el = domain->giveElement(elNum);

                Line line(xS, xE);
                if( line.intersects(el) ) {
                    displacementElements.push_back(elNum);
                }
            }

            mMapTractionElDispElGamma[i] = displacementElements;
        }
    }
    else{
        printf("Creating periodic tractions.\n");
        OOFEM_ERROR("Not implemented.");
    }


    // Create map from a traction element to displacement element and crack
    // intersection coordinates in the traction element.
    mTractionElInteriorCoordinates.clear();
    const double distTol = 1.0e-12;

    for(size_t i = 0; i < mpTractionElements.size(); i++) {

        const FloatArray &xS = mpTractionElements[i]->mStartCoord;
        const FloatArray &xE = mpTractionElements[i]->mEndCoord;

        PolygonLine pl;
        pl.insertVertexBack(xS);
        pl.insertVertexBack(xE);

        for(FloatArray x :  bndNodeCoordsFull) {

            double distN;
            pl.computeNormalSignDist(distN, x);

            double distT, arcPos;
            pl.computeTangentialSignDist(distT, x, arcPos);

            if( fabs(distN) < distTol && distT > -distTol && distT < (1.0+distTol) ) {
//                printf("Found point inside: "); x.printYourself();
                mTractionElInteriorCoordinates[i].push_back(x);
            }

        }

        // Sort based on distance to start point
        std::sort(mTractionElInteriorCoordinates[i].begin(), mTractionElInteriorCoordinates[i].end(), ArcPosSortFunction(xS) );

    }


    // Create map from a nodes global number to its local index
    // in the traction element.
    mNodeTractionElLocalInd.clear();
    mTracElDispNodes.clear();
    for(size_t i = 0; i < mpTractionElements.size(); i++) {

        std::vector<int> tracElDispNodes;

        for(int elIndex : mMapTractionElDispElGamma[i]) {
            const IntArray &elNodes = domain->giveElement(elIndex)->giveDofManArray();

            for(int nodeInd : elNodes) {
                tracElDispNodes.push_back(nodeInd);
            }
        }

//        printf("tracElDispNodes: ");
//        for(int ind : tracElDispNodes) {
//            printf("%d ", ind);
//        }
//        printf("\n");

        // We now an array of displacement nodes that affect the traction element.
        // It may contain duplicates and it is not sorted in any particular order.

        // Sort the array
        std::sort(tracElDispNodes.begin(), tracElDispNodes.end());

//        printf("tracElDispNodes sorted: ");
//        for(int ind : tracElDispNodes) {
//            printf("%d ", ind);
//        }
//        printf("\n");

        // Remove duplicates
        tracElDispNodes.erase( std::unique( tracElDispNodes.begin(), tracElDispNodes.end() ), tracElDispNodes.end() );

        //printf("tracElDispNodes unique: ");
        for(int j = 0; j < int(tracElDispNodes.size()); j++) {
            //printf("%d ", tracElDispNodes[j]);
            mNodeTractionElLocalInd[i][tracElDispNodes[j]] = j;
        }
        //printf("\n");

        mTracElDispNodes[i] = tracElDispNodes;
    }


    // Write traction nodes to debug vtk
    std :: vector<FloatArray> nodeCoord;
    for( Node *node : mpTractionNodes ) {
        nodeCoord.push_back( *(node->giveCoordinates()) );
    }

    std :: string fileName("TractionNodeCoord.vtk");
    XFEMDebugTools::WritePointsToVTK(fileName, nodeCoord);

}

void PrescribedGradientBCWeak::integrateTangent(FloatMatrix &oTangent, size_t iTracElInd)
{
    // Compute the tangent stiffness contribution
    // K = int (N^trac)^T . N^disp dGamma

    //printf("Entering PrescribedGradientBCWeak::integrateTangent.\n");
    const TractionElement &tEl = *(mpTractionElements[iTracElInd]);

    int dim = domain->giveNumberOfSpatialDimensions();

    IntegrationRule *ir = createNewIntegrationRule(iTracElInd);

    // Number of rows and colums
    int numRows = 0;
    if( mTractionInterpOrder == 0) {
        numRows = dim;
    }
    else if( mTractionInterpOrder == 1 ) {
        numRows = 2*dim;
    }

    int numCols = 0;
    std :: unordered_map<int, IntArray > globalNodeIndToPosInLocalLocArray;
    for(auto nodeInd : mTracElDispNodes[iTracElInd]) {
        DofManager *dMan = domain->giveDofManager(nodeInd);

        for(int i = 0; i < dMan->giveNumberOfDofs(); i++) {
            globalNodeIndToPosInLocalLocArray[nodeInd].followedBy( numCols+i+1 );
        }

        numCols += dMan->giveNumberOfDofs();
    }

    oTangent.resize(numRows, numCols);

    SpatialLocalizer *localizer = domain->giveSpatialLocalizer();

    for ( GaussPoint *gp: *ir ) {
        // Fetch local coordinates on traction element
        const FloatArray &locCoordsOnLine = * gp->giveLocalCoordinates();
        const FloatArray &globalCoord = * gp->giveCoordinates();

        //////////////////////////////////
        // Compute traction N-matrix
        FloatArray N, Ntrac;
        computeNTraction(Ntrac, locCoordsOnLine[0], tEl);
        tEl.computeN_Linear(N, locCoordsOnLine[0]);


        FloatMatrix NtracMat;
        NtracMat.beNMatrixOf(Ntrac, dim);


        //////////////////////////////////
        // Compute displacement N-matrix

        // Identify the displacement element
        // we are currently standing in
        // and compute local coordinates on
        // the displacement element
        FloatArray dispElLocCoord, closestPoint;
        Element *dispEl = localizer->giveElementClosestToPoint(dispElLocCoord, closestPoint, globalCoord);

        // Compute basis functions
        XfemElementInterface *xfemElInt = dynamic_cast<XfemElementInterface*> (dispEl);
        FloatMatrix NdispMat;

        if(xfemElInt != NULL && domain->hasXfemManager()) {
            //printf("Computing enriched N-matrix.\n");
            xfemElInt->XfemElementInterface_createEnrNmatrixAt(NdispMat, dispElLocCoord, *dispEl);
        }
        else {
            OOFEM_ERROR("Unable to compute N-matrix.")
        }


        FloatMatrix contrib;
        contrib.beTProductOf(NtracMat, NdispMat);
        double detJ = 0.5*tEl.mStartCoord.distance(tEl.mEndCoord);
        //printf("detJ: %e\n", detJ);
        //printf("gp->giveWeight(): %e\n", gp->giveWeight());
        contrib.times( detJ * gp->giveWeight() );


        // Create local location arrays
        IntArray rows, cols;
        for(int i = 1; i <= (mTractionInterpOrder+1)*dim; i++) {
            rows.followedBy(i);
        }

        const IntArray &dispElNodes = dispEl->giveDofManArray();
        for(int nodeInd : dispElNodes) {
            cols.followedBy( globalNodeIndToPosInLocalLocArray[nodeInd] );
        }

        oTangent.assemble(contrib, rows, cols);
    }

    delete ir;
}

IntegrationRule *PrescribedGradientBCWeak::createNewIntegrationRule(int iTracElInd)
{
    std::vector<FloatArray> tracGpCoord;

    // Create integration rule
    // -> need segments defined by displacement nodes, traction nodes and cracks.
    int numSeg = int(mTractionElInteriorCoordinates[iTracElInd].size())-1;
    std :: vector< Line > segments;

    for(int segIndex = 0; segIndex < numSeg; segIndex++) {
        const FloatArray &xS = mTractionElInteriorCoordinates[iTracElInd][segIndex];
        const FloatArray &xE = mTractionElInteriorCoordinates[iTracElInd][segIndex+1];
        segments.push_back( Line(xS, xE) );
    }

    Element *dummyEl = NULL;
    const FloatArray &xS = mTractionElInteriorCoordinates[iTracElInd][0];
    const FloatArray &xE = mTractionElInteriorCoordinates[iTracElInd].back();

    IntegrationRule *ir = new DiscontinuousSegmentIntegrationRule(1, dummyEl, segments, xS, xE);

    // Take material mode from first element
    MaterialMode matMode = domain->giveElement(1)->giveMaterialMode();

    int numPointsPerSeg = 3;
    ir->SetUpPointsOnLine(numPointsPerSeg, matMode);

    for(GaussPoint *gp : *ir) {
        tracGpCoord.push_back( *(gp->giveCoordinates()) );
    }

    std :: stringstream str3;
    str3 << "tracGpCoordEl" << iTracElInd << ".vtk";
    std :: string name3 = str3.str();

    XFEMDebugTools::WritePointsToVTK(name3, tracGpCoord);

    return ir;
}

void PrescribedGradientBCWeak::computeNTraction(FloatArray &oN, const double &iXi, const TractionElement &iEl) const
{
    if( mTractionInterpOrder == 0) {
        iEl.computeN_Constant(oN, iXi);
    }
    else if( mTractionInterpOrder == 1 ) {
        iEl.computeN_Linear(oN, iXi);
    }
}

void PrescribedGradientBCWeak::giveTractionUnknows(FloatArray &oTracUnknowns, ValueModeType mode, TimeStep *tStep, int iTracElInd)
{
    oTracUnknowns.clear();

    const std::vector<int> &tNodeInd = mpTractionElements[iTracElInd]->mTractionNodeInd;

    for(int nodeInd : tNodeInd) {
        FloatArray nodeUnknowns;
        mpTractionNodes[nodeInd]->giveCompleteUnknownVector(nodeUnknowns, mode, tStep);
        oTracUnknowns.append(nodeUnknowns);
    }

}

void PrescribedGradientBCWeak::giveDisplacementUnknows(FloatArray &oDispUnknowns, ValueModeType mode, TimeStep *tStep, int iTracElInd)
{
    oDispUnknowns.clear();

    for(int nodeInd : mTracElDispNodes[iTracElInd]) {
        FloatArray nodeUnknowns;
        domain->giveDofManager(nodeInd)->giveCompleteUnknownVector(nodeUnknowns, mode, tStep);
        oDispUnknowns.append(nodeUnknowns);
    }
}

double PrescribedGradientBCWeak :: domainSize()
{
    int nsd = this->domain->giveNumberOfSpatialDimensions();
    double domain_size = 0.0;
    // This requires the boundary to be consistent and ordered correctly.
    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);
        FEInterpolation *fei = e->giveInterpolation();
        domain_size += fei->evalNXIntegral( boundary, FEIElementGeometryWrapper(e) );
    }
    return fabs(domain_size / nsd);
}

} /* namespace oofem */
