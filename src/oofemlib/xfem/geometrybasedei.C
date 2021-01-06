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

#include "xfem/geometrybasedei.h"
#include "xfemmanager.h"

#include "spatiallocalizer.h"
#include "classfactory.h"
#include "element.h"
#include "node.h"
#include "mathfem.h"
#include "feinterpol.h"
#include "gausspoint.h"
#include "dynamicinputrecord.h"
#include "dynamicdatareader.h"
#include "geometry.h"

#include "xfem/enrichmentfunction.h"
#include "xfem/propagationlaw.h"
#include "xfem/enrichmentfronts/enrichmentfrontdonothing.h"

#include "engngm.h"
#include "timestep.h"

#include <string>
#include <algorithm>
#include <set>
#include <memory>

namespace oofem {
//REGISTER_EnrichmentItem(GeometryBasedEI)

GeometryBasedEI :: GeometryBasedEI(int n, XfemManager *xm, Domain *aDomain) :
    EnrichmentItem(n, xm, aDomain)
{}

GeometryBasedEI :: ~GeometryBasedEI()
{}

int GeometryBasedEI :: instanciateYourself(DataReader &dr)
{
    std :: string name;

    // Instantiate enrichment function
    {
        auto &mir = dr.giveInputRecord(DataReader :: IR_enrichFuncRec, 1);
        mir.giveRecordKeywordField(name);

        mpEnrichmentFunc = classFactory.createEnrichmentFunction( name.c_str(), 1, this->giveDomain() );
        if ( mpEnrichmentFunc ) {
            mpEnrichmentFunc->initializeFrom(mir);
        } else {
            OOFEM_ERROR( "failed to create enrichment function (%s)", name.c_str() );
        }
    }

    // Instantiate geometry
    {
        auto &mir = dr.giveInputRecord(DataReader :: IR_geoRec, 1);
        mir.giveRecordKeywordField(name);
        mpBasicGeometry = classFactory.createGeometry( name.c_str() );
        if ( !mpBasicGeometry ) {
            OOFEM_ERROR( "unknown geometry domain (%s)", name.c_str() );
        }

        mpBasicGeometry->initializeFrom(mir);
    }

    // Instantiate EnrichmentFront
    if ( mEnrFrontIndex == 0 ) {
        mpEnrichmentFrontStart = std::make_unique<EnrFrontDoNothing>();
        mpEnrichmentFrontEnd = std::make_unique<EnrFrontDoNothing>();
    } else {
        std :: string enrFrontNameStart, enrFrontNameEnd;

        auto &enrFrontStartIr = dr.giveInputRecord(DataReader :: IR_enrichFrontRec, mEnrFrontIndex);
        enrFrontStartIr.giveRecordKeywordField(enrFrontNameStart);

        mpEnrichmentFrontStart = classFactory.createEnrichmentFront( enrFrontNameStart.c_str() );
        if ( mpEnrichmentFrontStart ) {
            mpEnrichmentFrontStart->initializeFrom(enrFrontStartIr);
        } else {
            OOFEM_ERROR( "Failed to create enrichment front (%s)", enrFrontNameStart.c_str() );
        }

        auto &enrFrontEndIr = dr.giveInputRecord(DataReader :: IR_enrichFrontRec, mEnrFrontIndex);
        enrFrontEndIr.giveRecordKeywordField(enrFrontNameEnd);

        mpEnrichmentFrontEnd = classFactory.createEnrichmentFront( enrFrontNameEnd.c_str() );
        if ( mpEnrichmentFrontEnd ) {
            mpEnrichmentFrontEnd->initializeFrom(enrFrontEndIr);
        } else {
            OOFEM_ERROR( "Failed to create enrichment front (%s)", enrFrontNameEnd.c_str() );
        }
    }


    // Instantiate PropagationLaw
    if ( mPropLawIndex == 0 ) {
        mpPropagationLaw = std::make_unique<PLDoNothing>();
    } else {
        std :: string propLawName;

        auto &propLawir = dr.giveInputRecord(DataReader :: IR_propagationLawRec, mPropLawIndex);
        propLawir.giveRecordKeywordField(propLawName);

        mpPropagationLaw = classFactory.createPropagationLaw( propLawName.c_str() );
        if ( mpPropagationLaw ) {
            mpPropagationLaw->initializeFrom(propLawir);
        } else {
            OOFEM_ERROR( "Failed to create propagation law (%s)", propLawName.c_str() );
        }
    }

    // Set start of the enrichment dof pool for the given EI
    int xDofPoolAllocSize = this->giveDofPoolSize();
    this->startOfDofIdPool = this->giveDomain()->giveNextFreeDofID(xDofPoolAllocSize);
    this->endOfDofIdPool = this->startOfDofIdPool + xDofPoolAllocSize - 1;


    
    //    mpEnrichmentDomain->CallNodeEnrMarkerUpdate(* this, * xMan);
    // this->updateNodeEnrMarker(* xMan); // moved to postInitialize

//    writeVtkDebug();

    return 1;
}

void GeometryBasedEI :: postInitialize()
{
    XfemManager *xMan = this->giveDomain()->giveXfemManager();
    this->updateNodeEnrMarker(* xMan);
}

void GeometryBasedEI :: updateDofIdPool()
{
    // Set start of the enrichment dof pool for the given EI
    int xDofPoolAllocSize = this->giveDofPoolSize();
    this->startOfDofIdPool = this->giveDomain()->giveNextFreeDofID(xDofPoolAllocSize);
    this->endOfDofIdPool = this->startOfDofIdPool + xDofPoolAllocSize - 1;

//    printf("startOfDofIdPool: %d\n", startOfDofIdPool);
//    printf("endOfDofIdPool: %d\n", endOfDofIdPool);

    XfemManager *xMan = this->giveDomain()->giveXfemManager();
    //    mpEnrichmentDomain->CallNodeEnrMarkerUpdate(* this, * xMan);
    this->updateNodeEnrMarker(* xMan);
}

void GeometryBasedEI :: appendInputRecords(DynamicDataReader &oDR)
{
    auto eiRec = std::unique_ptr<DynamicInputRecord>();
    FEMComponent :: giveInputRecord(* eiRec);

    eiRec->setField(mEnrFrontIndex,                     _IFT_EnrichmentItem_front);
    eiRec->setField(mPropLawIndex,                      _IFT_EnrichmentItem_propagationlaw);

    if ( mInheritBoundaryConditions ) {
        eiRec->setField(_IFT_EnrichmentItem_inheritbc);
    }
    if ( mInheritOrderedBoundaryConditions ) {
        eiRec->setField(_IFT_EnrichmentItem_inheritorderedbc);
    }

    oDR.insertInputRecord(DataReader :: IR_enrichItemRec, std::move(eiRec));

    // Enrichment function
    auto efRec = std::unique_ptr<DynamicInputRecord>();
    mpEnrichmentFunc->giveInputRecord(* efRec);
    oDR.insertInputRecord(DataReader :: IR_enrichFuncRec, std::move(efRec));

    // Geometry
    auto geoRec = std::unique_ptr<DynamicInputRecord>();
    mpBasicGeometry->giveInputRecord(* geoRec);
    oDR.insertInputRecord(DataReader :: IR_geoRec, std::move(geoRec));


    // Enrichment front
    if ( mEnrFrontIndex != 0 ) {
        auto efrRecStart = std::unique_ptr<DynamicInputRecord>();
        mpEnrichmentFrontStart->giveInputRecord(* efrRecStart);
        oDR.insertInputRecord(DataReader :: IR_enrichFrontRec, std::move(efrRecStart));

        auto efrRecEnd = std::unique_ptr<DynamicInputRecord>();
        mpEnrichmentFrontEnd->giveInputRecord(* efrRecEnd);
        oDR.insertInputRecord(DataReader :: IR_enrichFrontRec, std::move(efrRecEnd));
    }

    // Propagation law
    if ( mPropLawIndex != 0 ) {
        auto plRec = std::unique_ptr<DynamicInputRecord>();
        this->mpPropagationLaw->giveInputRecord(* plRec);
        oDR.insertInputRecord(DataReader :: IR_propagationLawRec, std::move(plRec));
    }
}

void GeometryBasedEI :: updateGeometry()
{
    // Update enrichments ...
    XfemManager *xMan = this->giveDomain()->giveXfemManager();

    this->updateNodeEnrMarker(* xMan);
    // ... and create new dofs if necessary.
    createEnrichedDofs();
}

void GeometryBasedEI :: updateNodeEnrMarker(XfemManager &ixFemMan)
{
    updateLevelSets(ixFemMan);

    Domain *domain = giveDomain();
    SpatialLocalizer *localizer = domain->giveSpatialLocalizer();

    mNodeEnrMarkerMap.clear();
    TipInfo tipInfoStart, tipInfoEnd;
    bool foundTips = mpBasicGeometry->giveTips(tipInfoStart, tipInfoEnd);


    FloatArray center;
    double radius = 0.0;
    giveBoundingSphere(center, radius);


    IntArray elList;
    localizer->giveAllElementsWithNodesWithinBox(elList, center, radius);

    // Loop over elements and use the level sets to mark nodes belonging to completely cut elements.
    for ( int elNum: elList ) {
        Element *el = domain->giveElement(elNum);
        int nElNodes = el->giveNumberOfNodes();

        double minSignPhi  = 1, maxSignPhi         = -1;
        double minPhi = std :: numeric_limits< double > :: max();
        double maxPhi = std :: numeric_limits< double > :: min();

        FloatArray elCenter(2);
        elCenter.zero();

        for ( int elNodeInd = 1; elNodeInd <= nElNodes; elNodeInd++ ) {
            int nGlob = el->giveNode(elNodeInd)->giveGlobalNumber();

            double levelSetNormalNode = 0.0;
            if ( evalLevelSetNormalInNode( levelSetNormalNode, nGlob, el->giveNode(elNodeInd)->giveCoordinates() ) ) {
                minSignPhi = std :: min( sgn(minSignPhi), sgn(levelSetNormalNode) );
                maxSignPhi = std :: max( sgn(maxSignPhi), sgn(levelSetNormalNode) );

                minPhi = std :: min(minPhi, levelSetNormalNode);
                maxPhi = std :: max(maxPhi, levelSetNormalNode);
            }

            elCenter.at(1) += el->giveDofManager(elNodeInd)->giveCoordinate(1) / double ( nElNodes );
            elCenter.at(2) += el->giveDofManager(elNodeInd)->giveCoordinate(2) / double ( nElNodes );
        }


        int numEdgeIntersec = 0;

        if ( minPhi * maxPhi < mLevelSetTol ) { // If the level set function changes sign within the element.
            // Count the number of element edges intersected by the interface
            //int numEdges = nElNodes; // TODO: Is this assumption always true?
            int numEdges = el->giveInterpolation()->giveNumberOfEdges(); //JIM

            for ( int edgeIndex = 1; edgeIndex <= numEdges; edgeIndex++ ) {
                const auto &bNodes = el->giveInterpolation()->boundaryGiveNodes(edgeIndex);

                int niLoc = bNodes.at(1);
                int niGlob = el->giveNode(niLoc)->giveGlobalNumber();
                const auto &nodePosI = el->giveNode(niLoc)->giveCoordinates();
                int njLoc = bNodes.at(2);
                int njGlob = el->giveNode(njLoc)->giveGlobalNumber();
                const auto &nodePosJ = el->giveNode(njLoc)->giveCoordinates();

                double levelSetNormalNodeI = 0.0;
                double levelSetNormalNodeJ = 0.0;
                if ( evalLevelSetNormalInNode(levelSetNormalNodeI, niGlob, nodePosI) && evalLevelSetNormalInNode(levelSetNormalNodeJ, njGlob, nodePosJ) ) {
                    if ( levelSetNormalNodeI * levelSetNormalNodeJ < mLevelSetTol ) {
                        double xi = calcXiZeroLevel(levelSetNormalNodeI, levelSetNormalNodeJ);

                        // Compute the exact value of the tangential level set
                        // from the discretized geometry instead of interpolating.
                        double tangDist = 0.0, arcPos = 0.0;
                        const auto &posI = el->giveDofManager(niLoc)->giveCoordinates();
                        const auto &posJ = el->giveDofManager(njLoc)->giveCoordinates();
                        FloatArray pos;
                        pos.add(0.5 * ( 1.0 - xi ), posI);
                        pos.add(0.5 * ( 1.0 + xi ), posJ);
                        pos.resizeWithValues(2);

                        mpBasicGeometry->computeTangentialSignDist(tangDist, pos, arcPos);

                        double gamma = tangDist;

                        if ( gamma > 0.0 ) {
                            numEdgeIntersec++;
                        }
                    }
                }
            }


            if ( numEdgeIntersec >= 1 ) {
                // If we captured a cut element.
                for ( int elNodeInd = 1; elNodeInd <= nElNodes; elNodeInd++ ) {
                    int nGlob = el->giveNode(elNodeInd)->giveGlobalNumber();

                    auto res = mNodeEnrMarkerMap.find(nGlob);
                    if ( res == mNodeEnrMarkerMap.end() ) {
                        mNodeEnrMarkerMap [ nGlob ] = NodeEnr_BULK;
                    }
                }
            }
        }
    }

    // Mark tip nodes for special treatment.
    if(foundTips) {
		XfemManager *xMan = this->giveDomain()->giveXfemManager();
		mpEnrichmentFrontStart->MarkNodesAsFront(mNodeEnrMarkerMap, * xMan, mLevelSetNormalDirMap, mLevelSetTangDirMap, tipInfoStart);
		mpEnrichmentFrontEnd->MarkNodesAsFront(mNodeEnrMarkerMap, * xMan, mLevelSetNormalDirMap, mLevelSetTangDirMap, tipInfoEnd);
    }
}

void GeometryBasedEI :: updateLevelSets(XfemManager &ixFemMan)
{
    mLevelSetNormalDirMap.clear();
    mLevelSetTangDirMap.clear();

    FloatArray center;
    double radius = 0.0;
    giveBoundingSphere(center, radius);

    Domain *domain = giveDomain();
    SpatialLocalizer *localizer = domain->giveSpatialLocalizer();

    std :: list< int >nodeList;
    localizer->giveAllNodesWithinBox(nodeList, center, radius);

    for ( int nodeNum: nodeList ) {
        Node *node = ixFemMan.giveDomain()->giveNode(nodeNum);

        // Extract node coord
        FloatArray pos( node->giveCoordinates() );
        pos.resizeWithValues(2);

        // Calc normal sign dist
        double phi = 0.0;
        mpBasicGeometry->computeNormalSignDist(phi, pos);
        mLevelSetNormalDirMap [ nodeNum ] = phi;

        // Calc tangential sign dist
        double gamma = 0.0, arcPos = -1.0;
        mpBasicGeometry->computeTangentialSignDist(gamma, pos, arcPos);
        mLevelSetTangDirMap [ nodeNum ] = gamma;
    }

    mLevelSetsNeedUpdate = false;
}

void GeometryBasedEI :: evaluateEnrFuncInNode(std :: vector< double > &oEnrFunc, const Node &iNode) const
{
    double levelSetGP = 0.0;
    const FloatArray &globalCoord = iNode.giveCoordinates();
    int nodeInd = iNode.giveNumber();
    this->evalLevelSetNormalInNode(levelSetGP, nodeInd, globalCoord);

    //    const int dim = this->giveDomain()->giveNumberOfSpatialDimensions();
    //    FloatArray gradLevelSetGP(dim);

    double tangDist = 0.0, minDistArcPos = 0.0;
    mpBasicGeometry->computeTangentialSignDist(tangDist, globalCoord, minDistArcPos);

    FloatArray edGlobalCoord, localTangDir;

    mpBasicGeometry->giveGlobalCoordinates(edGlobalCoord, minDistArcPos);
    mpBasicGeometry->giveTangent(localTangDir, minDistArcPos);


    const EfInput efInput(globalCoord, levelSetGP, nodeInd, edGlobalCoord, minDistArcPos, localTangDir);


    if ( nodeInd == -1 ) {
        // Bulk enrichment
        oEnrFunc.resize(1, 0.0);
        mpEnrichmentFunc->evaluateEnrFuncAt(oEnrFunc [ 0 ], globalCoord, levelSetGP);
    } else {
        auto res = mNodeEnrMarkerMap.find(nodeInd);
        if ( res != mNodeEnrMarkerMap.end() ) {
            switch ( res->second ) {
            case NodeEnr_NONE:
                break;
            case NodeEnr_BULK:
                // Bulk enrichment
                oEnrFunc.resize(1, 0.0);
                mpEnrichmentFunc->evaluateEnrFuncAt(oEnrFunc [ 0 ], globalCoord, levelSetGP);
                break;
            case NodeEnr_START_TIP:
                mpEnrichmentFrontStart->evaluateEnrFuncAt(oEnrFunc, efInput);
                break;
            case NodeEnr_END_TIP:
                mpEnrichmentFrontEnd->evaluateEnrFuncAt(oEnrFunc, efInput);
                break;
            case NodeEnr_START_AND_END_TIP:
                mpEnrichmentFrontStart->evaluateEnrFuncAt(oEnrFunc, efInput);
                mpEnrichmentFrontEnd->evaluateEnrFuncAt(oEnrFunc, efInput);
                break;
            }
        } else {
            printf("In EnrichmentItem :: evaluateEnrFuncAt: mNodeEnrMarkerMap not found for iNodeInd %d\n", nodeInd);
        }
    }
}

void GeometryBasedEI :: evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iGlobalCoord, const FloatArray &iLocalCoord, int iNodeInd, const Element &iEl) const
{
    FloatArray N;
    iEl.giveInterpolation()->evalN( N, iLocalCoord, FEIElementGeometryWrapper(& iEl) );

    evaluateEnrFuncAt( oEnrFunc, iGlobalCoord, iLocalCoord, iNodeInd, iEl, N, iEl.giveDofManArray() );
}

void GeometryBasedEI :: evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iGlobalCoord, const FloatArray &iLocalCoord, int iNodeInd, const Element &iEl, const FloatArray &iN, const IntArray &iElNodes) const
{
    double levelSetGP = 0.0;
    evalLevelSetNormal(levelSetGP, iGlobalCoord, iN, iElNodes);

    //    const int dim = this->giveDomain()->giveNumberOfSpatialDimensions();
    //    FloatArray gradLevelSetGP(dim);

    double tangDist = 0.0, minDistArcPos = 0.0;
    const FloatArray globalCoord = {
        iGlobalCoord [ 0 ], iGlobalCoord [ 1 ]
    };
    mpBasicGeometry->computeTangentialSignDist(tangDist, globalCoord, minDistArcPos);

    FloatArray edGlobalCoord, localTangDir;

    mpBasicGeometry->giveGlobalCoordinates(edGlobalCoord, minDistArcPos);
    mpBasicGeometry->giveTangent(localTangDir, minDistArcPos);


    const EfInput efInput(iGlobalCoord, levelSetGP, iNodeInd, edGlobalCoord, minDistArcPos, localTangDir);


    if ( iNodeInd == -1 ) {
        // Bulk enrichment
        oEnrFunc.resize(1, 0.0);
        mpEnrichmentFunc->evaluateEnrFuncAt(oEnrFunc [ 0 ], iGlobalCoord, levelSetGP);
    } else {
        auto res = mNodeEnrMarkerMap.find(iNodeInd);
        if ( res != mNodeEnrMarkerMap.end() ) {
            switch ( res->second ) {
            case NodeEnr_NONE:
                break;
            case NodeEnr_BULK:
                // Bulk enrichment
                oEnrFunc.resize(1, 0.0);
                mpEnrichmentFunc->evaluateEnrFuncAt(oEnrFunc [ 0 ], iGlobalCoord, levelSetGP);
                break;
            case NodeEnr_START_TIP:
                mpEnrichmentFrontStart->evaluateEnrFuncAt(oEnrFunc, efInput);
                break;
            case NodeEnr_END_TIP:
                mpEnrichmentFrontEnd->evaluateEnrFuncAt(oEnrFunc, efInput);
                break;
            case NodeEnr_START_AND_END_TIP:
                mpEnrichmentFrontStart->evaluateEnrFuncAt(oEnrFunc, efInput);
                mpEnrichmentFrontEnd->evaluateEnrFuncAt(oEnrFunc, efInput);
                break;
            }
        } else {
            printf("In EnrichmentItem :: evaluateEnrFuncAt: mNodeEnrMarkerMap not found for iNodeInd %d\n", iNodeInd);
        }
    }
}

void GeometryBasedEI :: evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const FloatArray &iGlobalCoord, const FloatArray &iLocalCoord, int iNodeInd, const Element &iEl) const
{
    FloatMatrix dNdx;
    FloatArray N;
    FEInterpolation *interp = iEl.giveInterpolation();
    const FEIElementGeometryWrapper geomWrapper(& iEl);
    interp->evaldNdx(dNdx, iLocalCoord, geomWrapper);
    interp->evalN(N, iLocalCoord, geomWrapper);

    evaluateEnrFuncDerivAt( oEnrFuncDeriv, iGlobalCoord, iLocalCoord, iNodeInd, iEl, N, dNdx, iEl.giveDofManArray() );
}

void GeometryBasedEI :: evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const FloatArray &iGlobalCoord, const FloatArray &iLocalCoord, int iNodeInd, const Element &iEl, const FloatArray &iN, const FloatMatrix &idNdX, const IntArray &iElNodes) const
{
    auto res = mNodeEnrMarkerMap.find(iNodeInd);
    if ( res != mNodeEnrMarkerMap.end() ) {
        double levelSetGP = 0.0;
        evalLevelSetNormal(levelSetGP, iGlobalCoord, iN, iElNodes);

        const int dim = 2;
        FloatArray gradLevelSetGP(dim);
        evalGradLevelSetNormal(gradLevelSetGP, iGlobalCoord, idNdX, iElNodes);

        double tangDist = 0.0, minDistArcPos = 0.0;
        mpBasicGeometry->computeTangentialSignDist(tangDist, iGlobalCoord, minDistArcPos);

        FloatArray edGlobalCoord, localTangDir;


        mpBasicGeometry->giveGlobalCoordinates(edGlobalCoord, minDistArcPos);
        mpBasicGeometry->giveTangent(localTangDir, minDistArcPos);

        const EfInput efInput(iGlobalCoord, levelSetGP, iNodeInd, edGlobalCoord, minDistArcPos, localTangDir);

        switch ( res->second ) {
        case NodeEnr_NONE:
            break;
        case NodeEnr_BULK:
            oEnrFuncDeriv.resize(1);
            mpEnrichmentFunc->evaluateEnrFuncDerivAt(oEnrFuncDeriv [ 0 ], iGlobalCoord, levelSetGP, gradLevelSetGP);
            break;
        case NodeEnr_START_TIP:
            mpEnrichmentFrontStart->evaluateEnrFuncDerivAt(oEnrFuncDeriv, efInput, gradLevelSetGP);
            break;
        case NodeEnr_END_TIP:
            mpEnrichmentFrontEnd->evaluateEnrFuncDerivAt(oEnrFuncDeriv, efInput, gradLevelSetGP);
            break;
        case NodeEnr_START_AND_END_TIP:
            mpEnrichmentFrontStart->evaluateEnrFuncDerivAt(oEnrFuncDeriv, efInput, gradLevelSetGP);
            mpEnrichmentFrontEnd->evaluateEnrFuncDerivAt(oEnrFuncDeriv, efInput, gradLevelSetGP);
            break;
        }
    } else {
        printf("In EnrichmentItem :: evaluateEnrFuncDerivAt: mNodeEnrMarkerMap not found for iNodeInd %d\n", iNodeInd);
    }
}

void GeometryBasedEI :: evaluateEnrFuncJumps(std :: vector< double > &oEnrFuncJumps, int iNodeInd, GaussPoint &iGP, bool iGPLivesOnCurrentCrack) const
{
    double normalSignDist = 0.0;

    SpatialLocalizer *localizer = this->giveDomain()->giveSpatialLocalizer();

    Element *el = localizer->giveElementContainingPoint( iGP.giveGlobalCoordinates() );
    if ( el != NULL ) {
        FloatArray N;
        FEInterpolation *interp = el->giveInterpolation();
        interp->evalN( N, iGP.giveNaturalCoordinates(), FEIElementGeometryWrapper(el) );
        //el->computeLocalCoordinates(locCoord, iGlobalCoord);

        evalLevelSetNormal( normalSignDist, iGP.giveGlobalCoordinates(), N, el->giveDofManArray() );
    }

    //    this->interpLevelSet(normalSignDist, iGP.giveGlobalCoordinates() );

    auto res = mNodeEnrMarkerMap.find(iNodeInd);
    if ( res != mNodeEnrMarkerMap.end() ) {
        switch ( res->second ) {
        case NodeEnr_NONE:
            break;
        case NodeEnr_BULK:
            oEnrFuncJumps.resize(1);
            mpEnrichmentFunc->giveJump(oEnrFuncJumps);
            break;
        case NodeEnr_START_TIP:
            mpEnrichmentFrontStart->evaluateEnrFuncJumps(oEnrFuncJumps, iGP, iNodeInd, iGPLivesOnCurrentCrack, normalSignDist);
            break;
        case NodeEnr_END_TIP:
            mpEnrichmentFrontEnd->evaluateEnrFuncJumps(oEnrFuncJumps, iGP, iNodeInd, iGPLivesOnCurrentCrack, normalSignDist);
            break;
        case NodeEnr_START_AND_END_TIP:
            mpEnrichmentFrontStart->evaluateEnrFuncJumps(oEnrFuncJumps, iGP, iNodeInd, iGPLivesOnCurrentCrack, normalSignDist);
            mpEnrichmentFrontEnd->evaluateEnrFuncJumps(oEnrFuncJumps, iGP, iNodeInd, iGPLivesOnCurrentCrack, normalSignDist);
            break;
        }
    } else {
        printf("In EnrichmentItem :: evaluateEnrFuncDerivAt: evaluateEnrFuncJumps not found for iNodeInd %d\n", iNodeInd);
    }
}

void GeometryBasedEI :: computeIntersectionPoints(std :: vector< FloatArray > &oIntersectionPoints, std :: vector< int > &oIntersectedEdgeInd, Element *element, std :: vector< double > &oMinDistArcPos) const
{
    if ( isElementEnriched(element) ) {
        // Use the level set functions to compute intersection points

        // Loop over element edges; an edge is intersected if the
        // node values of the level set functions have different signs

        //      int numEdges = element->giveNumberOfBoundarySides();
        //int numEdges = element->giveNumberOfNodes(); // TODO: Is this assumption always true?
        int numEdges = element->giveInterpolation()->giveNumberOfEdges();

        for ( int edgeIndex = 1; edgeIndex <= numEdges; edgeIndex++ ) {
            const auto &bNodes = element->giveInterpolation()->boundaryGiveNodes(edgeIndex);

            int nsLoc = bNodes.at(1);
            int nsGlob = element->giveNode(nsLoc)->giveGlobalNumber();
            int neLoc = bNodes.at(2);
            int neGlob = element->giveNode(neLoc)->giveGlobalNumber();

            double phiS = 1.0;
            bool foundPhiS = evalLevelSetNormalInNode( phiS, nsGlob, element->giveNode(nsLoc)->giveCoordinates() );

            double phiE = 1.0;
            bool foundPhiE = evalLevelSetNormalInNode( phiE, neGlob, element->giveNode(neLoc)->giveCoordinates() );

            const auto &xS = element->giveNode(nsLoc)->giveCoordinates();
            const auto &xE = element->giveNode(neLoc)->giveCoordinates();
            const double edgeLength2 = distance_square(xS, xE);
            const double gammaRelTol = 1.0e-2;

            if ( ( foundPhiS && foundPhiE ) && phiS * phiE < mLevelSetRelTol * mLevelSetRelTol * edgeLength2 ) {
                // Intersection detected

                double xi = calcXiZeroLevel(phiS, phiE);

                // Compute the exact value of the tangential level set
                // from the discretized geometry instead of interpolating.
                double tangDist = 0.0, arcPos = 0.0;
                const auto &posI = element->giveDofManager(nsLoc)->giveCoordinates();
                const auto &posJ = element->giveDofManager(neLoc)->giveCoordinates();
                FloatArray pos;
                pos.add(0.5 * ( 1.0 - xi ), posI);
                pos.add(0.5 * ( 1.0 + xi ), posJ);
                pos.resizeWithValues(2);
                mpBasicGeometry->computeTangentialSignDist(tangDist, pos, arcPos);
                double gamma = tangDist;


                // If we are inside in tangential direction
                if ( gamma > -gammaRelTol * sqrt(edgeLength2) ) {
                    if ( fabs(phiS - phiE) < mLevelSetTol ) {
                        // If the crack is parallel to the edge.

                        FloatArray ps( element->giveDofManager(nsLoc)->giveCoordinates() );
                        ps.resizeWithValues(2);
                        FloatArray pe( element->giveDofManager(neLoc)->giveCoordinates() );
                        pe.resizeWithValues(2);

                        // Check that the intersection points have not already been identified.
                        // This may happen if the crack intersects the element exactly at a node,
                        // so that intersection is detected for both element edges in that node.

                        bool alreadyFound = false;

                        int numPointsOld = oIntersectionPoints.size();
                        for ( int k = 1; k <= numPointsOld; k++ ) {
                            double dist = distance(ps, oIntersectionPoints [ k - 1 ]);

                            if ( dist < mLevelSetTol ) {
                                alreadyFound = true;
                                break;
                            }
                        }

                        if ( !alreadyFound ) {
                            oIntersectionPoints.push_back(ps);

                            double arcPos = 0.0, tangDist = 0.0;
                            mpBasicGeometry->computeTangentialSignDist(tangDist, ps, arcPos);
                            oMinDistArcPos.push_back(arcPos);

                            oIntersectedEdgeInd.push_back(edgeIndex);
                        }

                        alreadyFound = false;

                        numPointsOld = oIntersectionPoints.size();
                        for ( int k = 1; k <= numPointsOld; k++ ) {
                            double dist = distance(pe, oIntersectionPoints [ k - 1 ]);

                            if ( dist < mLevelSetTol ) {
                                alreadyFound = true;
                                break;
                            }
                        }

                        if ( !alreadyFound ) {
                            oIntersectionPoints.push_back(pe);

                            double arcPos = 0.0, tangDist = 0.0;
                            mpBasicGeometry->computeTangentialSignDist(tangDist, pe, arcPos);
                            oMinDistArcPos.push_back(arcPos);

                            oIntersectedEdgeInd.push_back(edgeIndex);
                        }
                    } else {
                        const auto &ps = element->giveDofManager(nsLoc)->giveCoordinates();
                        const auto &pe = element->giveDofManager(neLoc)->giveCoordinates();

                        FloatArray p(2);

                        for ( int i = 1; i <= 2; i++ ) {
                            ( p.at(i) ) = 0.5 * ( 1.0 - xi ) * ( ( ps.at(i) ) ) + 0.5 * ( 1.0 + xi ) * ( ( pe.at(i) ) );
                        }


                        // Check that the intersection point has not already been identified.
                        // This may happen if the crack intersects the element exactly at a node,
                        // so that intersection is detected for both element edges in that node.

                        bool alreadyFound = false;


                        int numPointsOld = oIntersectionPoints.size();
                        for ( int k = 1; k <= numPointsOld; k++ ) {
                            double dist = distance(p, oIntersectionPoints [ k - 1 ]);

                            if ( dist < mLevelSetTol ) {
                                alreadyFound = true;
                                break;
                            }
                        }

                        if ( !alreadyFound ) {
                            oIntersectionPoints.push_back(p);

                            double arcPos = 0.0, tangDist = 0.0;
                            mpBasicGeometry->computeTangentialSignDist(tangDist, p, arcPos);
                            oMinDistArcPos.push_back(arcPos);

                            oIntersectedEdgeInd.push_back(edgeIndex);
                        }
                    }
                }
            }
        }
    }
}

void GeometryBasedEI :: computeIntersectionPoints(std :: vector< FloatArray > &oIntersectionPoints, std :: vector< int > &oIntersectedEdgeInd, Element *element, const Triangle &iTri, std :: vector< double > &oMinDistArcPos) const
{
    // Use the level set functions to compute intersection points

    // Loop over element edges; an edge is intersected if the
    // node values of the level set functions have different signs

    const int numEdges = 3;

    for ( int edgeIndex = 1; edgeIndex <= numEdges; edgeIndex++ ) {
        FloatArray xS, xE;

        // Global coordinates of vertices
        switch ( edgeIndex ) {
        case 1:
            xS = iTri.giveVertex(1);
            xE = iTri.giveVertex(2);
            break;
        case 2:
            xS = iTri.giveVertex(2);
            xE = iTri.giveVertex(3);
            break;

        case 3:
            xS = iTri.giveVertex(3);
            xE = iTri.giveVertex(1);
            break;
        default:
            break;
        }

        // Local coordinates of vertices
        FloatArray xiS;
        element->computeLocalCoordinates(xiS, xS);
        FloatArray xiE;
        element->computeLocalCoordinates(xiE, xE);

        const IntArray &elNodes = element->giveDofManArray();
        FloatArray Ns, Ne;
        FEInterpolation *interp = element->giveInterpolation();

        interp->evalN( Ns, xiS, FEIElementGeometryWrapper(element) );
        interp->evalN( Ne, xiE, FEIElementGeometryWrapper(element) );


        double phiS         = 0.0, phiE     = 0.0;
        double gammaS       = 0.0, gammaE   = 0.0;

        bool levelSetDefinedInAllNodes = true;
        for ( int i = 1; i <= Ns.giveSize(); i++ ) {
            double phiSNode = 0.0;
            if ( evalLevelSetNormalInNode( phiSNode, elNodes [ i - 1 ], element->giveNode(i)->giveCoordinates() ) ) {
                phiS += Ns.at(i) * phiSNode;
            } else   {
                levelSetDefinedInAllNodes = false;
            }

            double gammaSNode = 0.0;
            if ( evalLevelSetTangInNode( gammaSNode, elNodes [ i - 1 ], element->giveNode(i)->giveCoordinates() ) ) {
                gammaS += Ns.at(i) * gammaSNode;
            } else   {
                levelSetDefinedInAllNodes = false;
            }

            double phiENode = 0.0;
            if ( evalLevelSetNormalInNode( phiENode, elNodes [ i - 1 ], element->giveNode(i)->giveCoordinates() ) ) {
                phiE += Ne.at(i) * phiENode;
            } else   {
                levelSetDefinedInAllNodes = false;
            }

            double gammaENode = 0.0;
            if ( evalLevelSetTangInNode( gammaENode, elNodes [ i - 1 ], element->giveNode(i)->giveCoordinates() ) ) {
                gammaE += Ne.at(i) * gammaENode;
            } else   {
                levelSetDefinedInAllNodes = false;
            }
        }

        if ( phiS * phiE < mLevelSetTol2 && levelSetDefinedInAllNodes ) {
            // Intersection detected

            double xi = calcXiZeroLevel(phiS, phiE);
            double gamma = 0.5 * ( 1.0 - xi ) * gammaS + 0.5 * ( 1.0 + xi ) * gammaE;

            // If we are inside in tangential direction
            if ( gamma > 0.0 ) {
                if ( fabs(phiS - phiE) < mLevelSetTol ) {
                    // If the crack is parallel to the edge.

                    FloatArray ps(xS);
                    ps.resizeWithValues(2);
                    FloatArray pe(xE);
                    pe.resizeWithValues(2);

                    // Check that the intersection points have not already been identified.
                    // This may happen if the crack intersects the element exactly at a node,
                    // so that intersection is detected for both element edges in that node.

                    bool alreadyFound = false;

                    int numPointsOld = oIntersectionPoints.size();
                    for ( int k = 1; k <= numPointsOld; k++ ) {
                        double dist = distance(ps, oIntersectionPoints [ k - 1 ]);

                        if ( dist < mLevelSetTol ) {
                            alreadyFound = true;
                            break;
                        }
                    }

                    if ( !alreadyFound ) {
                        oIntersectionPoints.push_back(ps);

                        double arcPos = 0.0, tangDist = 0.0;
                        mpBasicGeometry->computeTangentialSignDist(tangDist, ps, arcPos);
                        oMinDistArcPos.push_back(arcPos);

                        oIntersectedEdgeInd.push_back(edgeIndex);
                    }

                    alreadyFound = false;

                    numPointsOld = oIntersectionPoints.size();
                    for ( int k = 1; k <= numPointsOld; k++ ) {
                        double dist = distance(pe, oIntersectionPoints [ k - 1 ]);

                        if ( dist < mLevelSetTol ) {
                            alreadyFound = true;
                            break;
                        }
                    }

                    if ( !alreadyFound ) {
                        oIntersectionPoints.push_back(pe);

                        double arcPos = 0.0, tangDist = 0.0;
                        mpBasicGeometry->computeTangentialSignDist(tangDist, pe, arcPos);
                        oMinDistArcPos.push_back(arcPos);

                        oIntersectedEdgeInd.push_back(edgeIndex);
                    }
                } else {
                    FloatArray ps(xS);
                    FloatArray pe(xE);

                    int nDim = std :: min( ps.giveSize(), pe.giveSize() );
                    FloatArray p;
                    p.resize(nDim);

                    for ( int i = 1; i <= nDim; i++ ) {
                        ( p.at(i) ) = 0.5 * ( 1.0 - xi ) * ( ( ps.at(i) ) ) + 0.5 * ( 1.0 + xi ) * ( ( pe.at(i) ) );
                    }


                    // Check that the intersection point has not already been identified.
                    // This may happen if the crack intersects the element exactly at a node,
                    // so that intersection is detected for both element edges in that node.

                    bool alreadyFound = false;


                    int numPointsOld = oIntersectionPoints.size();
                    for ( int k = 1; k <= numPointsOld; k++ ) {
                        double dist = distance(p, oIntersectionPoints [ k - 1 ]);

                        if ( dist < mLevelSetTol ) {
                            alreadyFound = true;
                            break;
                        }
                    }

                    if ( !alreadyFound ) {
                        oIntersectionPoints.push_back(p);

                        double arcPos = 0.0, tangDist = 0.0;
                        p.resizeWithValues(2);
                        mpBasicGeometry->computeTangentialSignDist(tangDist, p, arcPos);
                        oMinDistArcPos.push_back(arcPos);

                        oIntersectedEdgeInd.push_back(edgeIndex);
                    }
                }
            }
        }
    }
}

void GeometryBasedEI :: writeVtkDebug() const
{
    // For debugging only
	int tStepInd = 0;
	TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep(false);
	if(tStep != NULL) {
		tStepInd = tStep->giveNumber();
	}
    this->mpBasicGeometry->printVTK(tStepInd, number);
}

void GeometryBasedEI :: giveSubPolygon(std :: vector< FloatArray > &oPoints, const double &iXiStart, const double &iXiEnd) const
{
    mpBasicGeometry->giveSubPolygon(oPoints, iXiStart, iXiEnd);
}

void GeometryBasedEI :: propagateFronts(bool &oFrontsHavePropagated)
{
    oFrontsHavePropagated = false;

    TipPropagation tipPropStart;
    if ( mpPropagationLaw->propagateInterface(* giveDomain(), * mpEnrichmentFrontStart, tipPropStart) ) {
        //        mpEnrichmentDomain->propagateTip(tipPropStart);
        // TODO: Generalize
        // Propagate start point
        FloatArray pos( mpBasicGeometry->giveVertex(1) );
        pos.add(tipPropStart.mPropagationLength, tipPropStart.mPropagationDir);
        mpBasicGeometry->insertVertexFront(pos);

        oFrontsHavePropagated = true;
    }

    TipPropagation tipPropEnd;
    if ( mpPropagationLaw->propagateInterface(* giveDomain(), * mpEnrichmentFrontEnd, tipPropEnd) ) {
        //        mpEnrichmentDomain->propagateTip(tipPropEnd);
        // TODO: Generalize
        // Propagate end point
        FloatArray pos( mpBasicGeometry->giveVertex( mpBasicGeometry->giveNrVertices() ) );
        pos.add(tipPropEnd.mPropagationLength, tipPropEnd.mPropagationDir);
        mpBasicGeometry->insertVertexBack(pos);

        oFrontsHavePropagated = true;
    }

#if 0
    // For debugging only
    if ( mpEnrichmentDomain->getVtkDebug() ) {
        int tStepInd = this->domain->giveEngngModel()->giveCurrentStep()->giveNumber();

        EnrichmentDomain_BG *enrDomBG = dynamic_cast< EnrichmentDomain_BG * >( mpEnrichmentDomain );

        if ( enrDomBG != NULL ) {
            PolygonLine *pl = dynamic_cast< PolygonLine * >( enrDomBG->bg );
            if ( pl != NULL ) {
                pl->printVTK(tStepInd, number);
            }
        }
    }
#endif
    updateGeometry();

//    if( domain->giveEngngModel()->giveProblemScale() == macroScale ) {
//   	writeVtkDebug();
//    }
}

bool GeometryBasedEI :: giveElementTipCoord(FloatArray &oCoord, double &oArcPos,  Element &iEl, const FloatArray &iElCenter) const
{
    TipInfo tipInfoStart, tipInfoEnd;
    mpBasicGeometry->giveTips(tipInfoStart, tipInfoEnd);


    std :: vector< TipInfo >tipInfos = {
        tipInfoStart, tipInfoEnd
    };

    double minDist2 = std :: numeric_limits< double > :: max();
    size_t minIndex = 0;
    bool foundTip = false;

    for ( size_t i = 0; i < tipInfos.size(); i++ ) {
        double d2 = distance_square(tipInfos [ i ].mGlobalCoord, iElCenter);
        if ( d2 < minDist2 ) {
            minDist2 = d2;
            minIndex = i;
            foundTip = true;
        }
    }

    if ( !foundTip ) {
        return false;
    } else {
        // Check if the tip point is inside the element
        const FloatArray &globCoord = tipInfos [ minIndex ].mGlobalCoord;
        FloatArray locCoord;
        if ( !iEl.computeLocalCoordinates(locCoord, globCoord) ) {
            return false;
        }

        oCoord = tipInfos [ minIndex ].mGlobalCoord;
        oArcPos = tipInfos [ minIndex ].mArcPos;
        return true;
    }
}

void GeometryBasedEI :: giveBoundingSphere(FloatArray &oCenter, double &oRadius)
{
    // Compute bounding sphere from enrichment domain ...
    mpBasicGeometry->giveBoundingSphere(oCenter, oRadius);

    // ... increase the radius to cover the support of
    //     the enrichment front ...
    oRadius += max( mpEnrichmentFrontStart->giveSupportRadius(), mpEnrichmentFrontEnd->giveSupportRadius() );


    if ( domain->giveNumberOfSpatialDimensions() == 2 ) {
        // Compute mean area if applicable
        double meanArea = domain->giveArea() / domain->giveNumberOfElements();
        double meanLength = sqrt(meanArea);
        //printf("meanLength: %e\n", meanLength );

        oRadius = max(oRadius, meanLength);
    }

    // ... and make sure that all nodes of partly cut elements are included.
    oRadius *= 2.0;     // TODO: Compute a better estimate based on maximum element size. /ES
}
} /* namespace oofem */
