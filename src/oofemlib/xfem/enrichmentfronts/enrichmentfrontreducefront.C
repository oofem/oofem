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

#include "enrichmentfrontreducefront.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
#include "xfem/xfemmanager.h"
#include "domain.h"
#include "connectivitytable.h"
#include "spatiallocalizer.h"
#include "element.h"

namespace oofem {
REGISTER_EnrichmentFront(EnrFrontReduceFront)

void EnrFrontReduceFront :: MarkNodesAsFront(std :: unordered_map< int, int > &ioNodeEnrMarkerMap, XfemManager &ixFemMan, const std :: unordered_map< int, double > &iLevelSetNormalDirMap, const std :: unordered_map< int, double > &iLevelSetTangDirMap, const std :: vector< TipInfo > &iTipInfo)
{
    // Remove nodes touched by the crack tip
    Domain &d = * ( ixFemMan.giveDomain() );

    for ( size_t tipInd = 0; tipInd < iTipInfo.size(); tipInd++ ) {
        //      printf("iTipInfo[tipInd].mElIndex: %d\n", iTipInfo[tipInd].mElIndex );

        Element *el = d.giveSpatialLocalizer()->giveElementContainingPoint(iTipInfo [ tipInd ].mGlobalCoord);

        const IntArray &elNodes = el->giveDofManArray();

        for ( int i = 1; i <= elNodes.giveSize(); i++ ) {
            ioNodeEnrMarkerMap [ elNodes.at(i) ] = 0;
        }
    }
}

void EnrFrontReduceFront :: giveInputRecord(DynamicInputRecord &input)
{
    int number = 1;
    input.setRecordKeywordField(this->giveInputRecordName(), number);
}
} // end namespace oofem
