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

#include "enrichmentfront.h"
#include "xfem/tipinfo.h"
#include "domain.h"
#include "xfem/xfemmanager.h"
#include "spatiallocalizer.h"
#include "element.h"

namespace oofem {
void EnrichmentFront :: addTipIndexToNode(int iNodeInd, int iTipInd)
{
    // If the node is already enriched by the tip,
    // append the new index to the list.
    for ( size_t i = 0; i < mNodeTipIndices.size(); i++ ) {
        if ( mNodeTipIndices [ i ].first == iNodeInd ) {
            for ( size_t j = 0; j < mNodeTipIndices [ i ].second.size(); j++ ) {
                if ( mNodeTipIndices [ i ].second [ j ] == iTipInd ) {
                    // If the index is already present, we do not
                    // need to do anything
                    return;
                }
            }

            mNodeTipIndices [ i ].second.push_back(iTipInd);
            return;
        }
    }

    // If not, create a new pair
    std :: vector< int >tipIndices;
    tipIndices.push_back(iTipInd);
    std :: pair< int, std :: vector< int > >nodeTipInd = make_pair(iNodeInd, tipIndices);
    mNodeTipIndices.push_back(nodeTipInd);
}

void EnrichmentFront :: giveNodeTipIndices(int iNodeInd, std :: vector< int > &oTipIndices) const
{
    for ( size_t i = 0; i < mNodeTipIndices.size(); i++ ) {
        if ( mNodeTipIndices [ i ].first == iNodeInd ) {
            oTipIndices = mNodeTipIndices [ i ].second;
            return;
        }
    }
}

void EnrichmentFront :: MarkTipElementNodesAsFront(std :: unordered_map< int, int > &ioNodeEnrMarkerMap, XfemManager &ixFemMan,  const std :: unordered_map< int, double > &iLevelSetNormalDirMap, const std :: unordered_map< int, double > &iLevelSetTangDirMap, const std :: vector< TipInfo > &iTipInfo)
{
    mTipInfo = iTipInfo;
    mNodeTipIndices.clear();

    Domain &d = * ( ixFemMan.giveDomain() );

    for ( size_t tipInd = 0; tipInd < iTipInfo.size(); tipInd++ ) {
        Element *el = d.giveSpatialLocalizer()->giveElementContainingPoint(iTipInfo [ tipInd ].mGlobalCoord);

        if ( el != NULL ) {
            const IntArray &elNodes = el->giveDofManArray();

            for ( int i = 1; i <= elNodes.giveSize(); i++ ) {
                ioNodeEnrMarkerMap [ elNodes.at(i) ] = 2;
                addTipIndexToNode(elNodes.at(i), tipInd);
            }
        }
    }
}

} // end namespace oofem
