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
void EnrichmentFront :: MarkTipElementNodesAsFront(std :: unordered_map< int, NodeEnrichmentType > &ioNodeEnrMarkerMap, XfemManager &ixFemMan,  const std :: unordered_map< int, double > &iLevelSetNormalDirMap, const std :: unordered_map< int, double > &iLevelSetTangDirMap, const TipInfo &iTipInfo)
{
    mTipInfo = iTipInfo;

    Domain &d = * ( ixFemMan.giveDomain() );

    Element *el = d.giveSpatialLocalizer()->giveElementContainingPoint(mTipInfo.mGlobalCoord);

    if ( el != NULL ) {
        const IntArray &elNodes = el->giveDofManArray();

        for ( int i = 1; i <= elNodes.giveSize(); i++ ) {
            if ( ioNodeEnrMarkerMap [ elNodes.at(i) ] == NodeEnr_START_TIP || ioNodeEnrMarkerMap [ elNodes.at(i) ] == NodeEnr_END_TIP ) {
                ioNodeEnrMarkerMap [ elNodes.at(i) ] = NodeEnr_START_AND_END_TIP;
            } else   {
                if ( mTipInfo.mTipIndex == 0 ) {
                    ioNodeEnrMarkerMap [ elNodes.at(i) ] = NodeEnr_START_TIP;
                }

                if ( mTipInfo.mTipIndex == 1 ) {
                    ioNodeEnrMarkerMap [ elNodes.at(i) ] = NodeEnr_END_TIP;
                }
            }
        }
    }
}

void EnrichmentFront :: computeCrackTangent(FloatArray &oTangent, FloatArray &oNormal, bool &oFlipTangent, const EfInput &iEfInput) const
{
    oTangent = iEfInput.mLocalTangDir;

    if ( oTangent.dotProduct(mTipInfo.mTangDir) < 0.0 ) {
        oTangent.times(-1.0);
        oFlipTangent = true;
    } else   {
        oFlipTangent = false;
    }

    oNormal = {
        -oTangent.at(2), oTangent.at(1)
    };
}
} // end namespace oofem
