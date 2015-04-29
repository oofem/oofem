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
#include "xfemtolerances.h"
#include "spatiallocalizer.h"
#include "gausspoint.h"
#include "enrichmentfronts/enrichmentfront.h"
#include "enrichmentfronts/enrichmentfrontdonothing.h"
#include "engngm.h"

#include <algorithm>
#include <limits>
#include <sstream>
#include <string>


namespace oofem {
const double EnrichmentItem :: mLevelSetTol = 1.0e-12;
const double EnrichmentItem :: mLevelSetRelTol = 1.0e-3;


EnrichmentItem :: EnrichmentItem(int n, XfemManager *xMan, Domain *aDomain) : FEMComponent(n, aDomain),
    mpEnrichmentFunc(NULL),
    mpEnrichmentFrontStart(NULL),
    mpEnrichmentFrontEnd(NULL),
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
    delete mpEnrichmentFunc;
    delete mpEnrichmentFrontStart;
    delete mpEnrichmentFrontEnd;
    delete mpPropagationLaw;
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


int
EnrichmentItem :: giveDofPoolSize() const
{
    return this->giveEnrichesDofsWithIdArray()->giveSize() * this->giveNumberOfEnrDofs();
}

int
EnrichmentItem :: giveNumberOfEnrDofs() const
{
    int numEnrDofs = mpEnrichmentFunc->giveNumberOfDofs();

    if ( mpEnrichmentFrontStart != NULL ) {
        numEnrDofs = max( numEnrDofs, mpEnrichmentFrontStart->giveMaxNumEnrichments() );
    }

    if ( mpEnrichmentFrontEnd != NULL ) {
        numEnrDofs = max( numEnrDofs, mpEnrichmentFrontEnd->giveMaxNumEnrichments() );
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
        switch ( res->second ) {
        case NodeEnr_NONE:
            return 0;

            break;
        case NodeEnr_BULK:
            return 1;

            break;
        case NodeEnr_START_TIP:
            return mpEnrichmentFrontStart->giveNumEnrichments(iDMan);

            break;
        case NodeEnr_END_TIP:
            return mpEnrichmentFrontEnd->giveNumEnrichments(iDMan);

            break;
        case NodeEnr_START_AND_END_TIP:
            return mpEnrichmentFrontStart->giveNumEnrichments(iDMan) + mpEnrichmentFrontEnd->giveNumEnrichments(iDMan);

            break;
        }
    }

    return 0;
}

bool EnrichmentItem :: isMaterialModified(GaussPoint &iGP, Element &iEl, CrossSection * &opCS) const
{
    return false;
}

bool EnrichmentItem :: hasPropagatingFronts() const
{
    if ( mpPropagationLaw == NULL ) {
        return false;
    }

    return mpPropagationLaw->hasPropagation();
}

void
EnrichmentItem :: computeEnrichedDofManDofIdArray(IntArray &oDofIdArray, DofManager &iDMan)
{
    // Gives an array containing the dofId's that
    // are candidates for enrichment. At the moment,
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

bool EnrichmentItem :: evalLevelSetNormalInNode(double &oLevelSet, int iNodeInd, const FloatArray &iGlobalCoord) const
{
    auto res = mLevelSetNormalDirMap.find(iNodeInd);
    if ( res != mLevelSetNormalDirMap.end() ) {
        oLevelSet = res->second;
        return true;
    } else {
        oLevelSet = 0.0;
        return false;
    }
}

bool EnrichmentItem :: evalLevelSetTangInNode(double &oLevelSet, int iNodeInd, const FloatArray &iGlobalCoord) const
{
    auto res = mLevelSetTangDirMap.find(iNodeInd);
    if ( res != mLevelSetTangDirMap.end() ) {
        oLevelSet = res->second;
        return true;
    } else {
        oLevelSet = 0.0;
        return false;
    }
}

bool EnrichmentItem :: evalNodeEnrMarkerInNode(double &oNodeEnrMarker, int iNodeInd) const
{
    auto res = mNodeEnrMarkerMap.find(iNodeInd);
    if ( res != mNodeEnrMarkerMap.end() ) {
        oNodeEnrMarker = double( res->second );
        return true;
    } else {
        oNodeEnrMarker = 0.0;
        return false;
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
            //printf("dofMan %i is enriched \n", dMan->giveNumber());
            computeEnrichedDofManDofIdArray(dofIdArray, * dMan);
            for ( auto &dofid: dofIdArray ) {
                if ( !dMan->hasDofID( ( DofIDItem ) ( dofid ) ) ) {
                    if ( mInheritBoundaryConditions ) {
                        // Check if the other dofs in the dof manager have
                        // Dirichlet BCs. If so, let the new enriched dof
                        // inherit the same BC.
                        bool foundBC = false;
                        for ( Dof *dof: *dMan ) {
                            if ( dof->giveBcId() > 0 ) {
                                foundBC = true;
                                bcIndex = dof->giveBcId();
                                break;
                            }
                        }

                        if ( foundBC ) {
                            // Append dof with BC
                            dMan->appendDof( new MasterDof(dMan, bcIndex, icIndex, ( DofIDItem ) dofid) );
                        } else {
                            // No BC found, append enriched dof without BC
                            dMan->appendDof( new MasterDof(dMan, ( DofIDItem ) dofid) );
                        }
                    } else {
                        // Append enriched dof without BC
                        dMan->appendDof( new MasterDof(dMan, ( DofIDItem ) dofid) );
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

void EnrichmentItem :: calcPolarCoord(double &oR, double &oTheta, const FloatArray &iOrigin, const FloatArray &iPos, const FloatArray &iN, const FloatArray &iT, const EfInput &iEfInput, bool iFlipTangent)
{
    FloatArray q = {
        iPos.at(1) - iOrigin.at(1), iPos.at(2) - iOrigin.at(2)
    };

    const double tol = 1.0e-20;

    // Compute polar coordinates
    oR = iOrigin.distance(iPos);

    if ( oR > tol ) {
        q.times(1.0 / oR);
    }

    const double pi = M_PI;

    //    if( q.dotProduct(iT) > 0.0 ) {
    //        oTheta = asin( q.dotProduct(iN) );
    //    } else {
    //        if ( q.dotProduct(iN) > 0.0 ) {
    //            oTheta = pi - asin( q.dotProduct(iN) );
    //        } else {
    //            oTheta = -pi - asin( q.dotProduct(iN) );
    //        }
    //    }


    const double tol_q = 1.0e-3;
    double phi = iEfInput.mLevelSet;

    if ( iFlipTangent ) {
        phi *= -1.0;
    }

    double phi_r = 0.0;
    if ( oR > tol ) {
        phi_r = fabs(phi / oR);
    }

    if ( phi_r > 1.0 - XfemTolerances :: giveApproxZero() ) {
        phi_r = 1.0 - XfemTolerances :: giveApproxZero();
    }

    if ( iEfInput.mArcPos < tol_q || iEfInput.mArcPos > ( 1.0 - tol_q ) ) {
        double q_dot_n = q.dotProduct(iN);
        if ( q_dot_n > 1.0 - XfemTolerances :: giveApproxZero() ) {
            q_dot_n = 1.0 - XfemTolerances :: giveApproxZero();
        }

        oTheta = asin(q_dot_n);
    } else {
        if ( phi > 0.0 ) {
            oTheta = pi - asin( fabs(phi_r) );
        } else {
            oTheta = -pi + asin( fabs(phi_r) );
        }
    }
}

void EnrichmentItem :: callGnuplotExportModule(GnuplotExportModule &iExpMod, TimeStep *tStep)
{
    //iExpMod.outputXFEM(*this);
}

void EnrichmentItem :: setEnrichmentFrontStart(EnrichmentFront *ipEnrichmentFrontStart)
{
    delete mpEnrichmentFrontStart;
    mpEnrichmentFrontStart = ipEnrichmentFrontStart;
}

void EnrichmentItem :: setEnrichmentFrontEnd(EnrichmentFront *ipEnrichmentFrontEnd)
{
    delete mpEnrichmentFrontEnd;
    mpEnrichmentFrontEnd = ipEnrichmentFrontEnd;
}

bool EnrichmentItem :: tipIsTouchingEI(const TipInfo &iTipInfo)
{
    double tol = 1.0e-9;
    SpatialLocalizer *localizer = giveDomain()->giveSpatialLocalizer();

    Element *tipEl = localizer->giveElementContainingPoint(iTipInfo.mGlobalCoord);
    if ( tipEl != NULL ) {
        // Check if the candidate tip is located on the current crack
        FloatArray N;
        FloatArray locCoord;
        tipEl->computeLocalCoordinates(locCoord, iTipInfo.mGlobalCoord);
        FEInterpolation *interp = tipEl->giveInterpolation();
        interp->evalN( N, locCoord, FEIElementGeometryWrapper(tipEl) );

        double normalSignDist;
        evalLevelSetNormal( normalSignDist, iTipInfo.mGlobalCoord, N, tipEl->giveDofManArray() );

        double tangSignDist;
        evalLevelSetTangential( tangSignDist, iTipInfo.mGlobalCoord, N, tipEl->giveDofManArray() );

        if ( fabs(normalSignDist) < tol && tangSignDist > tol ) {
            return true;
        }
    }

    return false;
}
} // end namespace oofem
