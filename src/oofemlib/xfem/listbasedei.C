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

#include "classfactory.h"

#include <string>

#include "enrichmentdomain.h"

namespace oofem {

//REGISTER_EnrichmentItem(ListBasedEI)

ListBasedEI::ListBasedEI(int n, XfemManager *xm, Domain *aDomain) :
EnrichmentItem(n, xm, aDomain)
{
    // TODO Auto-generated constructor stub

}

ListBasedEI::~ListBasedEI() {
    // TODO Auto-generated destructor stub
}

void ListBasedEI::updateGeometry()
{
    // Update enrichments ...
    XfemManager *xMan = this->giveDomain()->giveXfemManager();
//    mpEnrichmentDomain->CallNodeEnrMarkerUpdate(* this, * xMan);

    this->updateNodeEnrMarker(*xMan);
    // ... and create new dofs if necessary.
    createEnrichedDofs();

}

void ListBasedEI::updateNodeEnrMarker(XfemManager &ixFemMan)
{
    //    updateLevelSets(ixFemMan);
    DofManList *dManList = dynamic_cast<DofManList*>(mpEnrichmentDomain);

    mNodeEnrMarkerMap.clear();

    //printf("\n The following nodes are enriched ");
    // Loop over nodes in the DofManList and mark nodes as enriched.
    const std :: vector< int > &dofList = dManList->giveDofManList();
    for ( int i = 0; i < int ( dofList.size() ); i++ ) {
        mNodeEnrMarkerMap [ dofList [ i ] ] = NodeEnr_BULK;
        //  printf(" %i", dofList [ i ]);
    }

    // Set level set fields to zero
    //mLevelSetSurfaceNormalDir.resize(nNodes, 0.0); // New /JB
}


} /* namespace oofem */
