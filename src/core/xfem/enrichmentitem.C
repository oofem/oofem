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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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
    mEnrFrontIndex(0),
    mPropLawIndex(0),
    mInheritBoundaryConditions(false),
    mInheritOrderedBoundaryConditions(false),
    startOfDofIdPool(-1),
    endOfDofIdPool(-1),
    mpEnrichesDofsWithIdArray(),
    mLevelSetsNeedUpdate(true),
    mLevelSetTol2(1.0e-12)
{}

EnrichmentItem :: ~EnrichmentItem()
{
}

void EnrichmentItem :: initializeFrom(InputRecord &ir)
{
    thisIr=ir.clone();
    mEnrFrontIndex = ir.giveGroupCount(_IFT_EnrichmentItem_front,"EnrichmentFront",/*optional*/true);
    mPropLawIndex = ir.hasChild(_IFT_EnrichmentItem_propagationlaw,"PropagationLaw",/*optional*/true);

    if ( ir.hasField(_IFT_EnrichmentItem_inheritbc) ) {
        mInheritBoundaryConditions = true;
    }
    if ( ir.hasField(_IFT_EnrichmentItem_inheritorderedbc) ) {
        mInheritOrderedBoundaryConditions = true;
    }
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

    if ( mpEnrichmentFrontStart ) {
        numEnrDofs = max( numEnrDofs, mpEnrichmentFrontStart->giveMaxNumEnrichments() );
    }

    if ( mpEnrichmentFrontEnd ) {
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
    int nodeInd = iDMan.giveGlobalNumber();
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
    if ( mpPropagationLaw == nullptr ) {
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
    int eiEnrSize = enrichesDofsWithIdArray->giveSize(); // Not necessarily true: for example, we may add 8 enrichments (4 in x and 4 in y) at a crack tip where we enrich D_u and D_v. /ES

    answer.resize(eiEnrSize);
    for ( int i = 1; i <= eiEnrSize; i++ ) {
        answer.at(i) = this->giveStartOfDofIdPool() + i - 1;
    }
}

void EnrichmentItem :: givePotentialEIDofIdArray(IntArray &answer) const {
    answer.clear();
    for ( int i = this->giveStartOfDofIdPool(); i <= this->giveEndOfDofIdPool(); i++ ) {
        answer.followedBy(i);
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
    IntArray EnrDofIdArray;

    mEIDofIdArray.clear();

    //int bcIndex = -1;
    int icIndex = -1;

    // Create new dofs
    for ( int i = 1; i <= nrDofMan; i++ ) {
        DofManager *dMan = this->giveDomain()->giveDofManager(i);

        if ( isDofManEnriched(* dMan) ) {
            //printf("dofMan %i is enriched \n", dMan->giveNumber());
            computeEnrichedDofManDofIdArray(EnrDofIdArray, * dMan);

            // Collect boundary condition ID of existing dofs
            IntArray bcIndexArray;
            for ( Dof *dof: *dMan ) {
                bcIndexArray.followedBy(dof->giveBcId());
            }

            bool foundBC = false;
            IntArray nonZeroBC;
            if ( !bcIndexArray.containsOnlyZeroes() ) {
                // BC is found on dofs  
                foundBC = true;
                nonZeroBC.findNonzeros(bcIndexArray);
            }

            int iDof(1);
            for ( auto &dofid: EnrDofIdArray ) {
                if ( !dMan->hasDofID( ( DofIDItem ) ( dofid ) ) ) {
                    if ( mInheritBoundaryConditions || mInheritOrderedBoundaryConditions ) {

                        if ( foundBC ) {
                            // Append dof with BC
                            if ( mInheritOrderedBoundaryConditions ) {
                                ///TODO: add choise of inheriting only specific BC. 
                                // Assume order type of new dofs are the same as original 
                                dMan->appendDof( new MasterDof(dMan, bcIndexArray.at(iDof), icIndex, ( DofIDItem ) dofid) );
                            } else {
                                // Append enriched dofs with same BC 
                                dMan->appendDof( new MasterDof(dMan, bcIndexArray.at(nonZeroBC.at(1)), icIndex, ( DofIDItem ) dofid) );
                            }
                        } else {
                            // No BC found, append enriched dof without BC
                            dMan->appendDof( new MasterDof(dMan, ( DofIDItem ) dofid) );
                        }
                    } else {
                        // Append enriched dof without BC
                        dMan->appendDof( new MasterDof(dMan, ( DofIDItem ) dofid) );
                    }
                }
                iDof++;
            }
        }
    }

    // Remove old dofs
    int poolStart       = giveStartOfDofIdPool();
    int poolEnd         = giveEndOfDofIdPool();

    for ( int i = 1; i <= nrDofMan; i++ ) {
        DofManager *dMan = this->giveDomain()->giveDofManager(i);

        computeEnrichedDofManDofIdArray(EnrDofIdArray, * dMan);
        std :: vector< DofIDItem >dofsToRemove;
        for ( auto &dof: *dMan ) {
            DofIDItem dofID = dof->giveDofID();

            if ( dofID >= DofIDItem(poolStart) && dofID <= DofIDItem(poolEnd) ) {
                bool dofIsInIdArray = false;
                for ( int k = 1; k <= EnrDofIdArray.giveSize(); k++ ) {
                    if ( dofID == DofIDItem( EnrDofIdArray.at(k) ) ) {
                        dofIsInIdArray = true;
                        break;
                    }
                }

                if ( !dofIsInIdArray ) {
                    dofsToRemove.push_back(dofID);
                }


                if(mEIDofIdArray.findFirstIndexOf(dofID) == 0 && dofIsInIdArray) {
                	mEIDofIdArray.followedBy(dofID);
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
    FloatArray q = Vec2(
        iPos.at(1) - iOrigin.at(1), iPos.at(2) - iOrigin.at(2)
    );

    const double tol = 1.0e-20;

    // Compute polar coordinates
    oR = distance(iOrigin, iPos);

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

void EnrichmentItem :: setPropagationLaw(std::unique_ptr<PropagationLaw> ipPropagationLaw)
{
    mpPropagationLaw = std::move(ipPropagationLaw);
    mPropLawIndex = 1;
}

void EnrichmentItem :: callGnuplotExportModule(GnuplotExportModule &iExpMod, TimeStep *tStep)
{
    //iExpMod.outputXFEM(*this);
}

void EnrichmentItem :: setEnrichmentFrontStart(std::unique_ptr<EnrichmentFront> ipEnrichmentFrontStart, bool iDeleteOld)
{
    if( ! iDeleteOld ) {
        mpEnrichmentFrontStart.release(); ///@note This is bad bad code
    }

    mpEnrichmentFrontStart = std::move(ipEnrichmentFrontStart);
}

void EnrichmentItem :: setEnrichmentFrontEnd(std::unique_ptr<EnrichmentFront> ipEnrichmentFrontEnd, bool iDeleteOld)
{
    if ( !iDeleteOld ) {
        mpEnrichmentFrontEnd.release(); ///@note This is bad bad code
    }

    mpEnrichmentFrontEnd = std::move(ipEnrichmentFrontEnd);
}

bool EnrichmentItem :: tipIsTouchingEI(const TipInfo &iTipInfo)
{
    double tol = 1.0e-9;
    SpatialLocalizer *localizer = giveDomain()->giveSpatialLocalizer();

    Element *tipEl = localizer->giveElementContainingPoint(iTipInfo.mGlobalCoord);
    if ( tipEl ) {
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
