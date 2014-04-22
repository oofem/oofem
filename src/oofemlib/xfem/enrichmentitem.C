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

#include "xfemmanager.h"
#include "floatmatrix.h"
#include "enrichmentitem.h"
#include "element.h"
#include "enrichmentfunction.h"
#include "enrichmentdomain.h"
#include "cltypes.h"
#include "connectivitytable.h"
#include "classfactory.h"
#include "fracturemanager.h"
#include "mathfem.h"
#include "feinterpol.h"
#include "masterdof.h"
#include "propagationlaw.h"
#include "dynamicinputrecord.h"
#include "dynamicdatareader.h"
#include "XFEMDebugTools.h"
#include "spatiallocalizer.h"
#include "gausspoint.h"
#include <algorithm>
#include <limits>
#include <sstream>

#include "enrichmentfronts/enrichmentfront.h"
#include "enrichmentfronts/enrichmentfrontdonothing.h"

namespace oofem {
const double EnrichmentItem :: mLevelSetTol = 1.0e-12;

REGISTER_EnrichmentItem(Inclusion)


EnrichmentItem :: EnrichmentItem(int n, XfemManager *xMan, Domain *aDomain) : FEMComponent(n, aDomain),
    mpEnrichmentDomain(NULL),
    mpEnrichmentFunc(NULL),
    mpEnrichmentFront(NULL),
    mEnrFrontIndex(0),
    mpPropagationLaw(NULL),
    mPropLawIndex(0),
    mInheritBoundaryConditions(false),
    startOfDofIdPool(-1),
    endOfDofIdPool(-1),
    mpEnrichesDofsWithIdArray(),
    mLevelSetsNeedUpdate(true),
    mLevelSetTol2(1.0e-12)
{}

EnrichmentItem :: ~EnrichmentItem()
{
    if ( mpEnrichmentDomain != NULL ) {
        delete mpEnrichmentDomain;
        mpEnrichmentDomain = NULL;
    }

    if ( mpEnrichmentFunc != NULL ) {
        delete mpEnrichmentFunc;
        mpEnrichmentFunc = NULL;
    }

    if ( mpEnrichmentFront != NULL ) {
        delete mpEnrichmentFront;
        mpEnrichmentFront = NULL;
    }

    if ( mpPropagationLaw != NULL ) {
        delete mpPropagationLaw;
        mpPropagationLaw = NULL;
    }
}

IRResultType EnrichmentItem :: initializeFrom(InputRecord *ir)
{
    IRResultType result; // Required by IR_GIVE_FIELD macro

    mEnrFrontIndex = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, mEnrFrontIndex, _IFT_EnrichmentItem_front);


    mPropLawIndex = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, mPropLawIndex, _IFT_EnrichmentItem_propagationlaw);

    if ( ir->hasField(_IFT_EnrichmentItem_inheritbc) ) {
        mInheritBoundaryConditions = true;
    }

    return IRRT_OK;
}

void EnrichmentItem :: appendInputRecords(DynamicDataReader &oDR)
{
    DynamicInputRecord *eiRec = new DynamicInputRecord();
    FEMComponent :: giveInputRecord(* eiRec);

    eiRec->setField(mEnrFrontIndex,                     _IFT_EnrichmentItem_front);
    eiRec->setField(mPropLawIndex,                      _IFT_EnrichmentItem_propagationlaw);

    if ( mInheritBoundaryConditions ) {
        eiRec->setField(_IFT_EnrichmentItem_inheritbc);
    }

    oDR.insertInputRecord(DataReader :: IR_enrichItemRec, eiRec);


    // Enrichment function
    DynamicInputRecord *efRec = new DynamicInputRecord();
    mpEnrichmentFunc->giveInputRecord(* efRec);
    oDR.insertInputRecord(DataReader :: IR_enrichFuncRec, efRec);


    // Enrichment domain
    DynamicInputRecord *edRec = new DynamicInputRecord();
    mpEnrichmentDomain->giveInputRecord(* edRec);
    oDR.insertInputRecord(DataReader :: IR_geoRec, edRec);


    // Enrichment front
    if ( mEnrFrontIndex != 0 ) {
        DynamicInputRecord *efrRec = new DynamicInputRecord();
        mpEnrichmentFront->giveInputRecord(* efrRec);
        oDR.insertInputRecord(DataReader :: IR_enrichFrontRec, efrRec);
    }

    // Propagation law
    if ( mPropLawIndex != 0 ) {
        DynamicInputRecord *plRec = new DynamicInputRecord();
        this->mpPropagationLaw->giveInputRecord(* plRec);
        oDR.insertInputRecord(DataReader :: IR_propagationLawRec, plRec);
    }
}


int EnrichmentItem :: instanciateYourself(DataReader *dr)
{
    IRResultType result; // Required by IR_GIVE_FIELD macro
    std :: string name;

    // Instantiate enrichment function
    InputRecord *mir = dr->giveInputRecord(DataReader :: IR_enrichFuncRec, 1);
    result = mir->giveRecordKeywordField(name);

    if ( result != IRRT_OK ) {
        mir->report_error(this->giveClassName(), __func__, "", result, __FILE__, __LINE__);
    }

    mpEnrichmentFunc = classFactory.createEnrichmentFunction( name.c_str(), 1, this->giveDomain() );
    if ( mpEnrichmentFunc != NULL ) {
        mpEnrichmentFunc->initializeFrom(mir);
    } else {
        OOFEM_ERROR( "failed to create enrichment function (%s)", name.c_str() );
    }


    // Instantiate enrichment domain
    mir = dr->giveInputRecord(DataReader :: IR_geoRec, 1);
    result = mir->giveRecordKeywordField(name);
    if ( result != IRRT_OK ) {
        mir->report_error(this->giveClassName(), __func__, "", result, __FILE__, __LINE__);
    }

    mpEnrichmentDomain = classFactory.createEnrichmentDomain( name.c_str() );
    if ( mpEnrichmentDomain == NULL ) {
        OOFEM_ERROR( "unknown enrichment domain (%s)", name.c_str() );
    }

    if ( giveDomain()->giveXfemManager()->giveVtkDebug() ) {
        mpEnrichmentDomain->setVtkDebug(true);
    }

    mpEnrichmentDomain->initializeFrom(mir);


    // Instantiate EnrichmentFront
    if ( mEnrFrontIndex == 0 ) {
        mpEnrichmentFront = new EnrFrontDoNothing();
    } else {
        std :: string enrFrontName;

        InputRecord *enrFrontir = dr->giveInputRecord(DataReader :: IR_enrichFrontRec, mEnrFrontIndex);
        result = enrFrontir->giveRecordKeywordField(enrFrontName);

        mpEnrichmentFront = classFactory.createEnrichmentFront( enrFrontName.c_str() );
        if ( mpEnrichmentFront != NULL ) {
            mpEnrichmentFront->initializeFrom(enrFrontir);
        } else {
            OOFEM_ERROR( "Failed to create enrichment front (%s)", enrFrontName.c_str() );
        }
    }


    // Instantiate PropagationLaw
    if ( mPropLawIndex == 0 ) {
        mpPropagationLaw = new PLDoNothing();
    } else {
        std :: string propLawName;

        InputRecord *propLawir = dr->giveInputRecord(DataReader :: IR_propagationLawRec, mPropLawIndex);
        result = propLawir->giveRecordKeywordField(propLawName);

        mpPropagationLaw = classFactory.createPropagationLaw( propLawName.c_str() );
        if ( mpPropagationLaw != NULL ) {
            mpPropagationLaw->initializeFrom(propLawir);
        } else {
            OOFEM_ERROR( "Failed to create propagation law (%s)", propLawName.c_str() );
        }
    }

    // Set start of the enrichment dof pool for the given EI
    int xDofPoolAllocSize = this->giveEnrichesDofsWithIdArray()->giveSize() * this->giveNumberOfEnrDofs() * 1 + 5; // TODO: overload for Crack
    this->startOfDofIdPool = this->giveDomain()->giveNextFreeDofID(xDofPoolAllocSize);
    this->endOfDofIdPool = this->startOfDofIdPool + xDofPoolAllocSize - 1;


    XfemManager *xMan = this->giveDomain()->giveXfemManager();
    mpEnrichmentDomain->CallNodeEnrMarkerUpdate(* this, * xMan);


    // For debugging only
    if ( mpEnrichmentDomain->getVtkDebug() ) {
        int tStepInd = 0; //this->domain->giveEngngModel()->giveCurrentStep()->giveNumber();

        EnrichmentDomain_BG *enrDomBG = dynamic_cast< EnrichmentDomain_BG * >( mpEnrichmentDomain );

        if ( enrDomBG != NULL ) {
            PolygonLine *pl = dynamic_cast< PolygonLine * >( enrDomBG->bg );
            if ( pl != NULL ) {
                pl->printVTK(tStepInd, number);
            }
        }
    }


    return 1;
}

int
EnrichmentItem :: giveNumberOfEnrDofs() const
{
    int numEnrDofs = mpEnrichmentFunc->giveNumberOfDofs();

    if ( mpEnrichmentFront != NULL ) {
        numEnrDofs = max( numEnrDofs, mpEnrichmentFront->giveMaxNumEnrichments() );
    }

    return numEnrDofs;
}

bool EnrichmentItem :: isElementEnriched(const Element *element) const
{
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        if ( this->isDofManEnriched( * ( element->giveDofManager(i) ) ) ) {
            return true;
        }
    }

    return false;
}

int EnrichmentItem :: giveNumDofManEnrichments(const DofManager &iDMan) const
{
    int nodeInd     = iDMan.giveGlobalNumber();
    auto res = mNodeEnrMarkerMap.find(nodeInd);

    if ( res != mNodeEnrMarkerMap.end() ) {
        if ( res->second == 1 ) {
            // Bulk enrichment
            return 1;
        } else {
            // Front enrichment
            return mpEnrichmentFront->giveNumEnrichments(iDMan);
        }
    } else   {
        return 0;
    }
}

int EnrichmentItem :: giveNumEnrichedDofs(const DofManager &iDMan) const
{
    int numEnrDofs = 0;

    int startId = giveStartOfDofIdPool();
    int endId = this->giveEndOfDofIdPool();

    for( Dof *dof: iDMan ) {
        int dofId = dof->giveDofID();

        // If the dof belongs to this enrichment item
        if ( dofId >= startId && dofId <= endId ) {
            numEnrDofs++;
        }
    }

    return numEnrDofs;
}

bool EnrichmentItem :: isMaterialModified(GaussPoint &iGP, Element &iEl, CrossSection * &opCS) const
{
    return false;
}

void EnrichmentItem :: updateGeometry()
{
    // Update enrichments ...
    XfemManager *xMan = this->giveDomain()->giveXfemManager();
    mpEnrichmentDomain->CallNodeEnrMarkerUpdate(* this, * xMan);

    // ... and create new dofs if necessary.
    createEnrichedDofs();
}

void EnrichmentItem :: propagateFronts()
{
    // Propagate interfaces
    mpPropagationLaw->propagateInterfaces(* giveDomain(), * mpEnrichmentDomain);

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

    updateGeometry();
}

void
EnrichmentItem :: computeEnrichedDofManDofIdArray(IntArray &oDofIdArray, DofManager &iDMan)
{
    // Gives an array containing the dofId's that
    // are candidated for enrichment. At the moment,
    // regular dofs are considered as candidates. In
    // the future, we may also consider enriching
    // enriched dofs from other enrichment items.
    const IntArray *enrichesDofsWithIdArray = this->giveEnrichesDofsWithIdArray();

    // Number of candidates for enrichment
    int numEnrCand = enrichesDofsWithIdArray->giveSize();

    // Number of active enrichment functions
    int numEnrFunc = this->giveNumDofManEnrichments(iDMan);

    // Go through the list of dofs that the EI supports
    // and compare with the available dofs in the dofMan.
    int count = 0;

    for ( int i = 1; i <= numEnrFunc; i++ ) {
        for ( int j = 1; j <= numEnrCand; j++ ) {
            if ( iDMan.hasDofID( ( DofIDItem ) enrichesDofsWithIdArray->at(j) ) ) {
                count++;
            }
        }
    }

    oDofIdArray.resize(count);
    for ( int i = 1; i <= count; i++ ) {
        oDofIdArray.at(i) = this->giveStartOfDofIdPool() + i - 1;
    }
}

void
EnrichmentItem :: giveEIDofIdArray(IntArray &answer) const
{
    // Returns an array containing the dof Id's of the new enrichment dofs pertinent to the ei.
    // Note: the dof managers may not support these dofs/all potential dof id's
    const IntArray *enrichesDofsWithIdArray = this->giveEnrichesDofsWithIdArray();
    int eiEnrSize = enrichesDofsWithIdArray->giveSize();

    answer.resize(eiEnrSize);
    for ( int i = 1; i <= eiEnrSize; i++ ) {
        answer.at(i) = this->giveStartOfDofIdPool() + i - 1;
    }
}

void EnrichmentItem :: evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iPos, const double &iLevelSet, int iNodeInd) const
{
    if ( iNodeInd == -1 ) {
        // Bulk enrichment
        oEnrFunc.resize(1, 0.0);
        mpEnrichmentFunc->evaluateEnrFuncAt(oEnrFunc [ 0 ], iPos, iLevelSet, mpEnrichmentDomain);
    } else {
        auto res = mNodeEnrMarkerMap.find(iNodeInd);
        if ( res != mNodeEnrMarkerMap.end() ) {
            if ( res->second == 1 ) {
                // Bulk enrichment
                oEnrFunc.resize(1, 0.0);
                mpEnrichmentFunc->evaluateEnrFuncAt(oEnrFunc [ 0 ], iPos, iLevelSet, mpEnrichmentDomain);
            } else   {
                // Front enrichment
                mpEnrichmentFront->evaluateEnrFuncAt(oEnrFunc, iPos, iLevelSet, iNodeInd);
            }
        } else   {
            printf("In EnrichmentItem :: evaluateEnrFuncAt: mNodeEnrMarkerMap not found for iNodeInd %d\n", iNodeInd);
        }
    }
}




void EnrichmentItem :: evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet, int iNodeInd) const
{
    auto res = mNodeEnrMarkerMap.find(iNodeInd);
    if ( res != mNodeEnrMarkerMap.end() ) {
        if ( res->second == 1 ) {
            // Bulk enrichment
            oEnrFuncDeriv.resize(1);
            mpEnrichmentFunc->evaluateEnrFuncDerivAt(oEnrFuncDeriv [ 0 ], iPos, iLevelSet, iGradLevelSet, mpEnrichmentDomain);
        } else   {
            // Front enrichment
            mpEnrichmentFront->evaluateEnrFuncDerivAt(oEnrFuncDeriv, iPos, iLevelSet, iGradLevelSet, iNodeInd);
        }
    } else   {
        printf("In EnrichmentItem :: evaluateEnrFuncDerivAt: mNodeEnrMarkerMap not found for iNodeInd %d\n", iNodeInd);
    }
}

void EnrichmentItem :: evaluateEnrFuncJumps(std :: vector< double > &oEnrFuncJumps, int iNodeInd, GaussPoint &iGP) const
{
    auto res = mNodeEnrMarkerMap.find(iNodeInd);
    if ( res != mNodeEnrMarkerMap.end() ) {
        if ( res->second == 1 ) {
            // Bulk enrichment
            oEnrFuncJumps.resize(1);
            mpEnrichmentFunc->giveJump(oEnrFuncJumps);
        } else   {
            // Front enrichment
            mpEnrichmentFront->evaluateEnrFuncJumps(oEnrFuncJumps, iGP, iNodeInd);
        }
    } else   {
        printf("In EnrichmentItem :: evaluateEnrFuncDerivAt: evaluateEnrFuncJumps not found for iNodeInd %d\n", iNodeInd);
    }
}

bool EnrichmentItem :: evalLevelSetNormalInNode(double &oLevelSet, int iNodeInd) const
{
    auto res = mLevelSetNormalDirMap.find(iNodeInd);
    if ( res != mLevelSetNormalDirMap.end() ) {
        oLevelSet = res->second;
        return true;
    } else   {
        oLevelSet = 0.0;
        return false;
    }
}

bool EnrichmentItem :: evalLevelSetTangInNode(double &oLevelSet, int iNodeInd) const
{
    auto res = mLevelSetTangDirMap.find(iNodeInd);
    if ( res != mLevelSetTangDirMap.end() ) {
        oLevelSet = res->second;
        return true;
    } else   {
        oLevelSet = 0.0;
        return false;
    }
}

bool EnrichmentItem :: evalNodeEnrMarkerInNode(double &oNodeEnrMarker, int iNodeInd) const
{
    auto res = mNodeEnrMarkerMap.find(iNodeInd);
    if ( res != mNodeEnrMarkerMap.end() ) {
        oNodeEnrMarker = res->second;
        return true;
    } else   {
        oNodeEnrMarker = 0.0;
        return false;
    }
}

bool EnrichmentItem :: levelSetChangesSignInEl(const IntArray &iElNodes) const
{
    double maxLevelSet = 0.0, minLevelSet = 0.0;
    double levelSetNode = 0.0;

    if ( evalLevelSetNormalInNode( levelSetNode, iElNodes.at(1) ) ) {
        maxLevelSet = levelSetNode;
        minLevelSet = levelSetNode;
    }

    for ( int j = 2; j < iElNodes.giveSize(); j++ ) {
        if ( evalLevelSetNormalInNode( levelSetNode, iElNodes.at(j) ) ) {
            maxLevelSet = std :: max(maxLevelSet, levelSetNode);
            minLevelSet = std :: min(minLevelSet, levelSetNode);
        }
    }

    if ( maxLevelSet * minLevelSet < 0.0 ) {
        return true;
    }

    return false;
}

void EnrichmentItem :: updateLevelSets(XfemManager &ixFemMan)
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
        const FloatArray &pos( * node->giveCoordinates() );

        // Calc normal sign dist
        double phi = 0.0;
        mpEnrichmentDomain->computeNormalSignDist(phi, pos);
        mLevelSetNormalDirMap [ nodeNum ] = phi;

        // Calc tangential sign dist
        double gamma = 0.0, arcPos = -1.0;
        mpEnrichmentDomain->computeTangentialSignDist(gamma, pos, arcPos);
        mLevelSetTangDirMap [ nodeNum ] = gamma;
    }

    mLevelSetsNeedUpdate = false;
}




void EnrichmentItem :: updateNodeEnrMarker(XfemManager &ixFemMan, const EnrichmentDomain_BG &iEnrichmentDomain_BG)
{
    updateLevelSets(ixFemMan);

    Domain *d = ixFemMan.giveDomain();
    mNodeEnrMarkerMap.clear();
    std :: vector< TipInfo >tipInfoArray;
    mpEnrichmentDomain->giveTipInfos(tipInfoArray);


    FloatArray center;
    double radius = 0.0;
    giveBoundingSphere(center, radius);

    Domain *domain = giveDomain();
    SpatialLocalizer *localizer = domain->giveSpatialLocalizer();

    std :: set< int >elList;
    localizer->giveAllElementsWithNodesWithinBox(elList, center, radius);

    // Loop over elements and use the level sets to mark nodes belonging to completely cut elements.
    for ( int elNum: elList ) {
        Element *el = d->giveElement(elNum);
        int nElNodes = el->giveNumberOfNodes();

        double minSignPhi  = 1, maxSignPhi         = -1;
        double minPhi = std :: numeric_limits< double > :: max();
        double maxPhi = std :: numeric_limits< double > :: min();

        FloatArray elCenter(2);
        elCenter.zero();

        for ( int elNodeInd = 1; elNodeInd <= nElNodes; elNodeInd++ ) {
            int nGlob = el->giveNode(elNodeInd)->giveGlobalNumber();

            double levelSetNormalNode = 0.0;
            if ( evalLevelSetNormalInNode(levelSetNormalNode, nGlob) ) {
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
            int numEdges = nElNodes; // TODO: Is this assumption always true?

            for ( int edgeIndex = 1; edgeIndex <= numEdges; edgeIndex++ ) {
                IntArray bNodes;
                el->giveInterpolation()->boundaryGiveNodes(bNodes, edgeIndex);

                int niLoc = bNodes.at(1);
                int niGlob = el->giveNode(niLoc)->giveGlobalNumber();
                int njLoc = bNodes.at( bNodes.giveSize() );
                int njGlob = el->giveNode(njLoc)->giveGlobalNumber();

                double levelSetNormalNodeI = 0.0;
                double levelSetNormalNodeJ = 0.0;
                if ( evalLevelSetNormalInNode(levelSetNormalNodeI, niGlob) && evalLevelSetNormalInNode(levelSetNormalNodeJ, njGlob) ) {
                    if ( levelSetNormalNodeI * levelSetNormalNodeJ < mLevelSetTol ) {
                        double xi = calcXiZeroLevel(levelSetNormalNodeI, levelSetNormalNodeJ);

                        // Compute the exact value of the tangential level set
                        // from the discretized geometry instead of interpolating.
                        double tangDist = 0.0, arcPos = 0.0;
                        const FloatArray &posI = * ( el->giveDofManager(niLoc)->giveCoordinates() );
                        const FloatArray &posJ = * ( el->giveDofManager(njLoc)->giveCoordinates() );
                        FloatArray pos;
                        pos.add(0.5 * ( 1.0 - xi ), posI);
                        pos.add(0.5 * ( 1.0 + xi ), posJ);
                        mpEnrichmentDomain->computeTangentialSignDist(tangDist, pos, arcPos);
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
                        mNodeEnrMarkerMap [ nGlob ] = 1;
                    }
                }
            }
        }
    }

    // Mark tip nodes for special treatment.
    XfemManager *xMan = this->giveDomain()->giveXfemManager();
    mpEnrichmentFront->MarkNodesAsFront(mNodeEnrMarkerMap, * xMan, mLevelSetNormalDirMap, mLevelSetTangDirMap, tipInfoArray);
}

void EnrichmentItem :: updateNodeEnrMarker(XfemManager &ixFemMan, const DofManList &iDofManList)
{
    updateLevelSets(ixFemMan);

    mNodeEnrMarkerMap.clear();

    //printf("\n The following nodes are enriched ");
    // Loop over nodes in the DofManList and mark nodes as enriched.
    const std :: vector< int > &dofList = iDofManList.giveDofManList();
    for ( int i = 0; i < int ( dofList.size() ); i++ ) {
        mNodeEnrMarkerMap [ dofList [ i ] ] = 1;
        //  printf(" %i", dofList [ i ]);
    }

    // Set level set fields to zero
    //mLevelSetSurfaceNormalDir.resize(nNodes, 0.0); // New /JB
}

void EnrichmentItem :: updateNodeEnrMarker(XfemManager &ixFemMan, const WholeDomain &iWholeDomain)
{
    // Mark all nodes for enrichment
    Domain *d = ixFemMan.giveDomain();
    int nNodes = d->giveNumberOfDofManagers();

    mNodeEnrMarkerMap.clear();
    for ( int i = 1; i <= nNodes; i++ ) {
        mNodeEnrMarkerMap [ i ] = 1;
    }
}

void EnrichmentItem :: createEnrichedDofs()
{
    // Creates new dofs due to the enrichment and appends them to the dof managers

    int nrDofMan = this->giveDomain()->giveNumberOfDofManagers();
    IntArray dofIdArray;

    int bcIndex = -1;
    int icIndex = -1;

    // Create new dofs
    for ( int i = 1; i <= nrDofMan; i++ ) {
        DofManager *dMan = this->giveDomain()->giveDofManager(i);

        if ( isDofManEnriched(* dMan) ) {
            int dofsPassed = dMan->giveNumberOfDofs();
            //printf("dofMan %i is enriched \n", dMan->giveNumber());
            computeEnrichedDofManDofIdArray(dofIdArray, * dMan);
            for ( int m = 1; m <= dofIdArray.giveSize(); m++ ) {
                dofsPassed++;

                if ( !dMan->hasDofID( ( DofIDItem ) ( dofIdArray.at(m) ) ) ) {
                    if ( mInheritBoundaryConditions ) {
                        // Check if the other dofs in the dof manager have
                        // Dirichlet BCs. If so, let the new enriched dof
                        // inherit the same BC.
                        bool foundBC = false;
                        for( Dof *dof: *dMan ) {
                            if ( dof->giveBcId() > 0 ) {
                                foundBC = true;
                                bcIndex = dof->giveBcId();
                                break;
                            }
                        }

                        if ( foundBC ) {
                            // Append dof with BC
                            dMan->appendDof( new MasterDof( dofsPassed, dMan, bcIndex, icIndex, ( DofIDItem ) ( dofIdArray.at(m) ) ) );
                        } else   {
                            // No BC found, append enriched dof without BC
                            dMan->appendDof( new MasterDof( dofsPassed, dMan, ( DofIDItem ) ( dofIdArray.at(m) ) ) );
                        }
                    } else   {
                        // Append enriched dof without BC
                        dMan->appendDof( new MasterDof( dofsPassed, dMan, ( DofIDItem ) ( dofIdArray.at(m) ) ) );
                    }
                }
            }
        }
    }

    // Remove old dofs
    int poolStart       = giveStartOfDofIdPool();
    int poolEnd         = giveEndOfDofIdPool();

    for ( int i = 1; i <= nrDofMan; i++ ) {
        DofManager *dMan = this->giveDomain()->giveDofManager(i);

        computeEnrichedDofManDofIdArray(dofIdArray, * dMan);
        std :: vector< DofIDItem >dofsToRemove;
        for ( Dof *dof: *dMan ) {
            DofIDItem dofID = dof->giveDofID();

            if ( dofID >= DofIDItem(poolStart) && dofID <= DofIDItem(poolEnd) ) {
                bool dofIsInIdArray = false;
                for ( int k = 1; k <= dofIdArray.giveSize(); k++ ) {
                    if ( dofID == DofIDItem( dofIdArray.at(k) ) ) {
                        dofIsInIdArray = true;
                        break;
                    }
                }

                if ( !dofIsInIdArray ) {
                    dofsToRemove.push_back(dofID);
                }
            }
        }

        for ( size_t j = 0; j < dofsToRemove.size(); j++ ) {
            dMan->removeDof(dofsToRemove [ j ]);
        }
    }
}

void EnrichmentItem :: computeIntersectionPoints(std :: vector< FloatArray > &oIntersectionPoints, std :: vector< int > &oIntersectedEdgeInd, Element *element, std :: vector< double > &oMinDistArcPos) const
{
    if ( isElementEnriched(element) ) {
        // Use the level set functions to compute intersection points

        // Loop over element edges; an edge is intersected if the
        // node values of the level set functions have different signs

        //		int numEdges = element->giveNumberOfBoundarySides();
        int numEdges = element->giveNumberOfNodes(); // TODO: Is this assumption always true?

        for ( int edgeIndex = 1; edgeIndex <= numEdges; edgeIndex++ ) {
            IntArray bNodes;
            element->giveInterpolation()->boundaryGiveNodes(bNodes, edgeIndex);

            int nsLoc = bNodes.at(1);
            int nsGlob = element->giveNode(nsLoc)->giveGlobalNumber();
            int neLoc = bNodes.at( bNodes.giveSize() );
            int neGlob = element->giveNode(neLoc)->giveGlobalNumber();

            double phiS = 1.0;
            bool foundPhiS = evalLevelSetNormalInNode(phiS, nsGlob);

            double phiE = 1.0;
            bool foundPhiE = evalLevelSetNormalInNode(phiE, neGlob);

            if ( (foundPhiS && foundPhiE) && phiS * phiE < mLevelSetTol2 ) {
                // Intersection detected

                double xi = calcXiZeroLevel(phiS, phiE);

                // Compute the exact value of the tangential level set
                // from the discretized geometry instead of interpolating.
                double tangDist = 0.0, arcPos = 0.0;
                const FloatArray &posI = * ( element->giveDofManager(nsLoc)->giveCoordinates() );
                const FloatArray &posJ = * ( element->giveDofManager(neLoc)->giveCoordinates() );
                FloatArray pos;
                pos.add(0.5 * ( 1.0 - xi ), posI);
                pos.add(0.5 * ( 1.0 + xi ), posJ);
                mpEnrichmentDomain->computeTangentialSignDist(tangDist, pos, arcPos);
                double gamma = tangDist;


                // If we are inside in tangential direction
                if ( gamma > 0.0 ) {
                    if ( fabs(phiS - phiE) < mLevelSetTol ) {
                        // If the crack is parallel to the edge.

                        FloatArray ps( * ( element->giveDofManager(nsLoc)->giveCoordinates() ) );
                        FloatArray pe( * ( element->giveDofManager(neLoc)->giveCoordinates() ) );

                        // Check that the intersection points have not already been identified.
                        // This may happen if the crack intersects the element exactly at a node,
                        // so that intersection is detected for both element edges in that node.

                        bool alreadyFound = false;

                        int numPointsOld = oIntersectionPoints.size();
                        for ( int k = 1; k <= numPointsOld; k++ ) {
                            double dist = ps.distance(oIntersectionPoints [ k - 1 ]);

                            if ( dist < mLevelSetTol ) {
                                alreadyFound = true;
                                break;
                            }
                        }

                        if ( !alreadyFound ) {
                            oIntersectionPoints.push_back(ps);

                            double arcPos = 0.0, tangDist = 0.0;
                            mpEnrichmentDomain->computeTangentialSignDist(tangDist, ps, arcPos);
                            oMinDistArcPos.push_back(arcPos);

                            oIntersectedEdgeInd.push_back(edgeIndex);
                        }

                        alreadyFound = false;

                        numPointsOld = oIntersectionPoints.size();
                        for ( int k = 1; k <= numPointsOld; k++ ) {
                            double dist = pe.distance(oIntersectionPoints [ k - 1 ]);

                            if ( dist < mLevelSetTol ) {
                                alreadyFound = true;
                                break;
                            }
                        }

                        if ( !alreadyFound ) {
                            oIntersectionPoints.push_back(pe);

                            double arcPos = 0.0, tangDist = 0.0;
                            mpEnrichmentDomain->computeTangentialSignDist(tangDist, pe, arcPos);
                            oMinDistArcPos.push_back(arcPos);

                            oIntersectedEdgeInd.push_back(edgeIndex);
                        }
                    } else {
                        FloatArray ps( * ( element->giveDofManager(nsLoc)->giveCoordinates() ) );
                        FloatArray pe( * ( element->giveDofManager(neLoc)->giveCoordinates() ) );

                        int nDim = ps.giveSize();
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
                            double dist = p.distance(oIntersectionPoints [ k - 1 ]);

                            if ( dist < mLevelSetTol ) {
                                alreadyFound = true;
                                break;
                            }
                        }

                        if ( !alreadyFound ) {
                            oIntersectionPoints.push_back(p);

                            double arcPos = 0.0, tangDist = 0.0;
                            mpEnrichmentDomain->computeTangentialSignDist(tangDist, p, arcPos);
                            oMinDistArcPos.push_back(arcPos);

                            oIntersectedEdgeInd.push_back(edgeIndex);
                        }
                    }
                }
            }
        }
    }
}

void EnrichmentItem :: computeIntersectionPoints(std :: vector< FloatArray > &oIntersectionPoints, std :: vector< int > &oIntersectedEdgeInd, Element *element, const Triangle &iTri, std :: vector< double > &oMinDistArcPos) const
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

        for ( int i = 1; i <= Ns.giveSize(); i++ ) {
            double phiSNode = 0.0;
            if ( evalLevelSetNormalInNode(phiSNode, elNodes [ i - 1 ]) ) {
                phiS += Ns.at(i) * phiSNode;
            }

            double gammaSNode = 0.0;
            if ( evalLevelSetTangInNode(gammaSNode, elNodes [ i - 1 ]) ) {
                gammaS += Ns.at(i) * gammaSNode;
            }

            double phiENode = 0.0;
            if ( evalLevelSetNormalInNode(phiENode, elNodes [ i - 1 ]) ) {
                phiE += Ne.at(i) * phiENode;
            }

            double gammaENode = 0.0;
            if ( evalLevelSetTangInNode(gammaENode, elNodes [ i - 1 ]) ) {
                gammaE += Ne.at(i) * gammaENode;
            }
        }

        if ( phiS * phiE < mLevelSetTol2 ) {
            // Intersection detected

            double xi = calcXiZeroLevel(phiS, phiE);
            double gamma = 0.5 * ( 1.0 - xi ) * gammaS + 0.5 * ( 1.0 + xi ) * gammaE;

            // If we are inside in tangential direction
            if ( gamma > 0.0 ) {
                if ( fabs(phiS - phiE) < mLevelSetTol ) {
                    // If the crack is parallel to the edge.

                    FloatArray ps(xS);
                    FloatArray pe(xE);

                    // Check that the intersection points have not already been identified.
                    // This may happen if the crack intersects the element exactly at a node,
                    // so that intersection is detected for both element edges in that node.

                    bool alreadyFound = false;

                    int numPointsOld = oIntersectionPoints.size();
                    for ( int k = 1; k <= numPointsOld; k++ ) {
                        double dist = ps.distance(oIntersectionPoints [ k - 1 ]);

                        if ( dist < mLevelSetTol ) {
                            alreadyFound = true;
                            break;
                        }
                    }

                    if ( !alreadyFound ) {
                        oIntersectionPoints.push_back(ps);

                        double arcPos = 0.0, tangDist = 0.0;
                        mpEnrichmentDomain->computeTangentialSignDist(tangDist, ps, arcPos);
                        oMinDistArcPos.push_back(arcPos);

                        oIntersectedEdgeInd.push_back(edgeIndex);
                    }

                    alreadyFound = false;

                    numPointsOld = oIntersectionPoints.size();
                    for ( int k = 1; k <= numPointsOld; k++ ) {
                        double dist = pe.distance(oIntersectionPoints [ k - 1 ]);

                        if ( dist < mLevelSetTol ) {
                            alreadyFound = true;
                            break;
                        }
                    }

                    if ( !alreadyFound ) {
                        oIntersectionPoints.push_back(pe);

                        double arcPos = 0.0, tangDist = 0.0;
                        mpEnrichmentDomain->computeTangentialSignDist(tangDist, pe, arcPos);
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
                        double dist = p.distance(oIntersectionPoints [ k - 1 ]);

                        if ( dist < mLevelSetTol ) {
                            alreadyFound = true;
                            break;
                        }
                    }

                    if ( !alreadyFound ) {
                        oIntersectionPoints.push_back(p);

                        double arcPos = 0.0, tangDist = 0.0;
                        mpEnrichmentDomain->computeTangentialSignDist(tangDist, p, arcPos);
                        oMinDistArcPos.push_back(arcPos);

                        oIntersectedEdgeInd.push_back(edgeIndex);
                    }
                }
            }
        }
    }
}


bool EnrichmentItem :: giveElementTipCoord(FloatArray &oCoord, double &oArcPos, int iElIndex, const FloatArray &iElCenter) const
{
    std :: vector< TipInfo >tipInfos;
    mpEnrichmentDomain->giveTipInfos(tipInfos);

    double minDist2 = std :: numeric_limits< double > :: max();
    size_t minIndex = 0;
    bool foundTip = false;

    for ( size_t i = 0; i < tipInfos.size(); i++ ) {
        if ( tipInfos [ i ].mGlobalCoord.distance_square(iElCenter) < minDist2 ) {
            minDist2 = tipInfos [ i ].mGlobalCoord.distance_square(iElCenter);
            minIndex = i;
            foundTip = true;
        }
    }

    if ( !foundTip ) {
        return false;
    } else   {
        oCoord = tipInfos [ minIndex ].mGlobalCoord;
        oArcPos = tipInfos [ minIndex ].mArcPos;
        return true;
    }
}

bool EnrichmentItem :: giveElementTipCoord(FloatArray &oCoord, double &oArcPos, int iElIndex, const Triangle &iTri, const FloatArray &iElCenter) const
{
    std :: vector< TipInfo >tipInfos;
    mpEnrichmentDomain->giveTipInfos(tipInfos);

    double minDist2 = std :: numeric_limits< double > :: max();
    size_t minIndex = 0;
    bool foundTip = false;

    for ( size_t i = 0; i < tipInfos.size(); i++ ) {
        if ( tipInfos [ i ].mGlobalCoord.distance_square(iElCenter) < minDist2 ) {
            if ( iTri.pointIsInTriangle(tipInfos [ i ].mGlobalCoord) ) {
                minDist2 = tipInfos [ i ].mGlobalCoord.distance_square(iElCenter);
                minIndex = i;
                foundTip = true;
            }
        }
    }

    if ( !foundTip ) {
        return false;
    } else   {
        oCoord = tipInfos [ minIndex ].mGlobalCoord;
        oArcPos = tipInfos [ minIndex ].mArcPos;
        return true;
    }
}

double EnrichmentItem :: calcXiZeroLevel(const double &iQ1, const double &iQ2)
{
    double xi = 0.0;

    if ( fabs(iQ1 - iQ2) > mLevelSetTol ) {
        xi = ( iQ1 + iQ2 ) / ( iQ1 - iQ2 );
    }

    if ( xi < -1.0 ) {
        xi = -1.0;
    }

    if ( xi > 1.0 ) {
        xi = 1.0;
    }

    return xi;
}

void EnrichmentItem :: calcPolarCoord(double &oR, double &oTheta, const FloatArray &iOrigin, const FloatArray &iPos, const FloatArray &iN, const FloatArray &iT)
{
    FloatArray q = {
        iPos.at(1) - iOrigin.at(1), iPos.at(2) - iOrigin.at(2)
    };

    const double tol = 1.0e-12;

    // Compute polar coordinates
    oR = iOrigin.distance(iPos);

    if ( oR > tol ) {
        q.times(1.0 / oR);
    }

    if ( q.dotProduct(iN) > 0.0 ) {
        oTheta =  acos( q.dotProduct(iT) );
    } else {
        oTheta = -acos( q.dotProduct(iT) );
    }
}

void EnrichmentItem :: giveSubPolygon(std :: vector< FloatArray > &oPoints, const double &iXiStart, const double &iXiEnd) const
{
    mpEnrichmentDomain->giveSubPolygon(oPoints, iXiStart, iXiEnd);
}

void EnrichmentItem :: callGnuplotExportModule(GnuplotExportModule &iExpMod)
{
    //iExpMod.outputXFEM(*this);
}

void EnrichmentItem :: giveBoundingSphere(FloatArray &oCenter, double &oRadius)
{
    // Compute bounding sphere from enrichment domain ...
    mpEnrichmentDomain->giveBoundingSphere(oCenter, oRadius);

    // ... increase the radius to cover the support of
    //	   the enrichment front ...
    oRadius += mpEnrichmentFront->giveSupportRadius();

    // ... and make sure that all nodes of partly cut elements are included.
    oRadius *= 2.0;     // TODO: Compute a better estimate based on maximum element size. /ES
}


Inclusion :: Inclusion(int n, XfemManager *xm, Domain *aDomain) :
    EnrichmentItem(n, xm, aDomain),
    mpCrossSection(NULL)
{
    mpEnrichesDofsWithIdArray = {
        D_u, D_v, D_w
    };
}

Inclusion :: ~Inclusion()
{
    if ( mpCrossSection != NULL ) {
        mpCrossSection = NULL;
    }
}

bool Inclusion :: isMaterialModified(GaussPoint &iGP, Element &iEl, CrossSection * &opCS) const
{
    // Check if the point is located inside the inclusion

    FloatArray N;
    FEInterpolation *interp = iEl.giveInterpolation();
    interp->evalN( N, * iGP.giveCoordinates(), FEIElementGeometryWrapper(& iEl) );

    const IntArray &elNodes = iEl.giveDofManArray();

    double levelSetGP = 0.0;
    this->interpLevelSet(levelSetGP, N, elNodes);

    if ( levelSetGP < 0.0 ) {
        opCS = mpCrossSection;
        return true;
    }

    return false;
}

IRResultType Inclusion :: initializeFrom(InputRecord *ir)
{
    EnrichmentItem :: initializeFrom(ir);
    IRResultType result;
    int crossSectionIndex = 0;
    IR_GIVE_FIELD(ir, crossSectionIndex, _IFT_Inclusion_CrossSection);
    mpCrossSection = this->giveDomain()->giveCrossSection(crossSectionIndex);

    return IRRT_OK;
}
} // end namespace oofem
