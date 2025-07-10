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

#include "xfem/listbasedei.h"
#include "xfemmanager.h"
#include "domain.h"
#include "classfactory.h"
#include "engngm.h"
#include "timestep.h"
#include "xfem/propagationlaw.h"

#include <string>
#include <algorithm>

namespace oofem {
//REGISTER_EnrichmentItem(ListBasedEI)

ListBasedEI :: ListBasedEI(int n, XfemManager *xm, Domain *aDomain) :
    EnrichmentItem(n, xm, aDomain)
{}

ListBasedEI :: ~ListBasedEI()
{}

void ListBasedEI :: updateGeometry()
{
    // Update enrichments ...
    XfemManager *xMan = this->giveDomain()->giveXfemManager();
    //    mpEnrichmentDomain->CallNodeEnrMarkerUpdate(* this, * xMan);

    this->updateNodeEnrMarker(* xMan);
    // ... and create new dofs if necessary.
    createEnrichedDofs();
}

void ListBasedEI :: propagateFronts(bool &oFrontsHavePropagated)
{
    oFrontsHavePropagated = false;
    
    TipPropagation tipProp;
    if ( mpPropagationLaw->propagateInterface(* giveDomain(), * mpEnrichmentFrontStart, tipProp) ) {
        // TODO: Generalize
        // Propagate front
        
        for ( int i = 1; i <= tipProp.mPropagationDofManNumbers.giveSize(); i++ ) {
            //std::list< int > :: iterator p;
            std :: vector< int > :: iterator p;
            p = std :: find( this->dofManList.begin(), this->dofManList.end(), tipProp.mPropagationDofManNumbers.at(i) );
            if ( p == this->dofManList.end() ) {          // if new node
                this->dofManList.push_back( tipProp.mPropagationDofManNumbers.at(i) );
            }
        }

        std :: sort( dofManList.begin(), this->dofManList.end() );

        oFrontsHavePropagated = true;
    }
    
    this->updateGeometry();
    
    
}

void 
ListBasedEI :: initiateFronts(bool &oFrontsHavePropagated, IntArray &initiateDofMans) 
{    
    oFrontsHavePropagated = false;
    
    printf("\n Enrichment %i - The following nodes have been initiated: ",this->giveNumber());
    for ( int i = 1; i <= initiateDofMans.giveSize(); i++ ) {
        //std::list< int > :: iterator p;
        std :: vector< int > :: iterator p;
        p = std :: find( this->dofManList.begin(), this->dofManList.end(), initiateDofMans.at(i) );
        if ( p == this->dofManList.end() ) {          // if new node
            printf(" %i", initiateDofMans.at(i) ); 
            this->dofManList.push_back( initiateDofMans.at(i) );
        }
    }
    printf(" \n");

    std :: sort( dofManList.begin(), this->dofManList.end() );
    
    oFrontsHavePropagated = true;
    
    this->updateGeometry();
    
}

void ListBasedEI :: updateNodeEnrMarker(XfemManager &ixFemMan)
{
    mNodeEnrMarkerMap.clear();
    TipInfo tipInfo;
    //TipInfo *tipInfo = mpEnrichmentFrontStart->giveTipInfo();
    //tipInfo->mTipDofManNumbers.clear();
    //this->mpEnrichmentFrontStart->giveTipInfo().mTipDofManNumbers.clear();
    
    // Update the dofManList and the enrichment boundaries.
    // For now let all nodes in dofManList be the boundary
    
    // Loop over nodes in the DofManList and mark nodes as enriched.
    bool printed = false;
    for ( auto &dman : dofManList ) {
        mNodeEnrMarkerMap [ dman ] = NodeEnr_BULK;
        if ( !printed ) {
            printf("\n Enrichment %i - The following nodes are enriched:",this->giveNumber());
            printed = true;
        }
        printf(" %i", dman ); 
        //TODO change this so only the boundaries are added to tipInfo, now all nodes are added   
        tipInfo.mTipDofManNumbers.insertSorted(dman);
        //this->mpEnrichmentFrontStart->giveTipInfo().mTipDofManNumbers.insertSorted(dman);
    }
    printf(" \n");
    //tipInfo.mTipDofManNumbers.printYourself("tip dofManNumbers (updateNodeEnrMarker)");
    //this->mpEnrichmentFrontStart->giveTipInfo().mTipDofManNumbers.printYourself("tip dofManNumbers (updateNodeEnrMarker)");
    
    
    // Update tipInfo
    mpEnrichmentFrontStart->setTipInfo( tipInfo );

    // Set level set fields to zero
    //mLevelSetSurfaceNormalDir.resize(nNodes, 0.0); // New /JB
}

bool ListBasedEI :: giveElementTipCoord(FloatArray &oCoord, double &oArcPos,  Element &iEl, const FloatArray &iElCenter) const
{
    // No tips can be defined for list based enrichment items.
    // Hence, we return false.
    return false;
}
} /* namespace oofem */
