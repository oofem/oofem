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

#include <cmath>

namespace oofem {
REGISTER_BoundaryCondition(PrescribedGradientBCWeak);

PrescribedGradientBCWeak::PrescribedGradientBCWeak(int n, Domain * d):
PrescribedGradientBC(n, d),
mTractionInterpOrder(0),
mNumTractionNodesAtIntersections(1),
mTractionNodeSpacing(1),
mTractionLivesOnGammaPlus(false)
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

//    printf("mLC: "); mLC.printYourself();
//    printf("mUC: "); mUC.printYourself();

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

}


void PrescribedGradientBCWeak::assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
                      CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{

}

void PrescribedGradientBCWeak::giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, EquationID eid, CharType type,
                                const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{

}

void PrescribedGradientBCWeak::createTractionMesh()
{
    printf("Entering PrescribedGradientBCWeak::createTractionMesh().\n");

    std::vector<FloatArray> bndNodeCoords;

    if(!mTractionLivesOnGammaPlus) {
        printf("Creating Dirichlet tractions.\n");

        Set *setPointer = this->giveDomain()->giveSet(this->set);
        const IntArray &boundaries = setPointer->giveBoundaryList();
        IntArray bNodes;

        // Loop over all boundary segments
        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
            int boundary = boundaries.at(pos * 2);

            e->giveInterpolation()->boundaryGiveNodes(bNodes, boundary);
//            printf("bNodes: "); bNodes.printYourself();

            // Add the start node of the segment
            DofManager *startNode = e->giveDofManager( bNodes[0] );
//            printf("startNode: "); startNode->giveCoordinates()->printYourself();
            bndNodeCoords.push_back( *(startNode->giveCoordinates()) );

            // Add traction nodes where cracks intersect the boundary if desired
            XfemElementInterface *xfemElInt = dynamic_cast<XfemElementInterface*> (e);

            if(xfemElInt != NULL && domain->hasXfemManager() ) {
                std::vector<Line> segments;
                std::vector<FloatArray> intersecPoints;
                xfemElInt->partitionEdgeSegment( boundary, segments, intersecPoints );

                for(size_t i = 0; i < intersecPoints.size(); i++) {
//                    printf("i: %lu intersecPoints[i]: ", i); intersecPoints[i].printYourself();

                    for(int j = 0; j < mNumTractionNodesAtIntersections; j++) {
                        bndNodeCoords.push_back( intersecPoints[i] );
                    }
                }
            }

        }

        // TODO: Sort bndNodeCoords and pick every n:th node.

        // Create traction dofs
        int nsd = domain->giveNumberOfSpatialDimensions();
//        printf("nsd: %d\n", nsd);
        std::vector<int> dofIds;
        for(int j = 0; j < nsd; j++) {
            dofIds.push_back( this->domain->giveNextFreeDofID() );
        }
        for(size_t i = 0; i < bndNodeCoords.size(); i++) {
//            printf("i: %d bndNodeCoords[i]: ", i); bndNodeCoords[i].printYourself();

            int numNodes = domain->giveNumberOfDofManagers();
            Node *node = new Node(numNodes+1, domain);
            for(int j = 0; j < nsd; j++) {
//                printf("dofIds[j]: %d\n", dofIds[j]);
                node->appendDof( new MasterDof( j + 1, node, ( DofIDItem ) ( dofIds[j] ) ) );
            }

            node->setCoordinates(bndNodeCoords[i]);
            mpTractionNodes.push_back(node);
        }

        // Create traction elements
        if(mTractionInterpOrder == 0) {
            for(size_t i = 0; i < mpTractionNodes.size(); i++) {
                TractionElement *tractionEl = new TractionElement();
//                printf("i: %d\n", i );
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
            localizer->giveAllElementsWithNodesWithinBox(elList, xC, 0.51*elLength );

//            printf("Close displacement elements: ");
            std :: vector< int > displacementElements;
            for ( int elNum: elList ) {
//                printf("elNum: %d\n", elNum);

                // Check if the traction element and the displacement element intersect
                // Intersection occurs if at least one displacement element node is
                // on the traction element.
                Element *el = domain->giveElement(elNum);

                Line line(xS, xE);
                if( line.intersects(el) ) {
//                    printf("line.intersects(el).\n");
                    displacementElements.push_back(i);
                }
            }

            mMapTractionElDispElGamma[i] = displacementElements;
        }
    }
    else{
        printf("Creating periodic tractions.\n");
    }
}

} /* namespace oofem */
