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

#include "mathfem.h"
#include "enrichmentdomain.h"
#include "enrichmentitem.h"
#include "element.h"
#include "dofmanager.h"
#include "connectivitytable.h"
#include "classfactory.h"
#include "enrichmentfunction.h"
#include "xfemmanager.h"
#include "dynamicinputrecord.h"
#include "set.h"

#include <algorithm>

#include <cmath>

namespace oofem {
REGISTER_EnrichmentDomain(DofManList)
REGISTER_EnrichmentDomain(WholeDomain)
REGISTER_EnrichmentDomain(EDBGCircle)

REGISTER_EnrichmentDomain(EDCrack)


EnrichmentDomain :: EnrichmentDomain() :
    mDebugVTK(false)
{ }

void EnrichmentDomain_BG :: giveInputRecord(DynamicInputRecord &input)
{
    input.setRecordKeywordField(this->giveInputRecordName(), 1);

    bg->giveInputRecord(input);
}

void EnrichmentDomain_BG :: giveBoundingSphere(FloatArray &oCenter, double &oRadius)
{
    int nVert = bg->giveNrVertices();
    oCenter = {
        0.0, 0.0
    };
    oRadius = 0.0;

    if ( nVert > 0 ) {
        for ( int i = 1; i <= nVert; i++ ) {
            oCenter.add( bg->giveVertex(i) );
        }

        oCenter.times( 1.0 / double( nVert ) );

        for ( int i = 1; i <= nVert; i++ ) {
            oRadius = std :: max( oRadius, oCenter.distance( bg->giveVertex(i) ) );
        }
    }
}

void EDBGCircle :: giveBoundingSphere(FloatArray &oCenter, double &oRadius)
{
    Circle *circle = dynamic_cast< Circle * >( bg );

    if ( circle == NULL ) {
        OOFEM_ERROR("In EDBGCircle::giveBoundingSphere(): Failed to cast to Circle.")
    }

    oCenter = bg->giveVertex(1);
    oRadius = circle->giveRadius();
}


IRResultType EDCrack :: initializeFrom(InputRecord *ir)
{
    IRResultType result = bg->initializeFrom(ir);

    return result;
}

void
DofManList :: addDofManagers(IntArray &dofManNumbers)
{
    for ( int i = 1; i <= dofManNumbers.giveSize(); i++ ) {
        //std::list< int > :: iterator p;
        std :: vector< int > :: iterator p;
        p = std :: find( this->dofManList.begin(), this->dofManList.end(), dofManNumbers.at(i) );
        if ( p == this->dofManList.end() ) {          // if new node
            this->dofManList.push_back( dofManNumbers.at(i) );
        }
    }

    std :: sort( dofManList.begin(), this->dofManList.end() );
}

bool EDCrack :: giveTipInfos(TipInfo &oStartTipInfo, TipInfo &oEndTipInfo) const
{
    int nVert = bg->giveNrVertices();
    if ( nVert > 1 ) {
        // Start tip
        TipInfo info1;
        const FloatArray &p1S = ( bg->giveVertex(1) );
        const FloatArray &p2S = ( bg->giveVertex(2) );

        // Tip position
        info1.mGlobalCoord = p1S;

        // Tip tangent
        info1.mTangDir.beDifferenceOf(p1S, p2S);
        info1.mTangDir.normalize();

        // Tip normal
        info1.mNormalDir = {
            -info1.mTangDir.at(2), info1.mTangDir.at(1)
        };

        info1.mTipIndex = 0;
        info1.mArcPos = 0.0;

        oStartTipInfo = info1;

        // End tip
        TipInfo info2;
        const FloatArray &p1E = ( bg->giveVertex(nVert - 1) );
        const FloatArray &p2E = ( bg->giveVertex(nVert) );

        // Tip position
        info2.mGlobalCoord = p2E;

        // Tip tangent
        info2.mTangDir.beDifferenceOf(p2E, p1E);
        info2.mTangDir.normalize();

        // Tip normal
        info2.mNormalDir = {
            -info2.mTangDir.at(2), info2.mTangDir.at(1)
        };

        info2.mTipIndex = 1;
        info2.mArcPos = 1.0;

        oEndTipInfo = info2;

        return true;
    }

    return false;
}

bool EDCrack :: propagateTip(const TipPropagation &iTipProp) {
        if ( iTipProp.mTipIndex == 0 ) {
            // Propagate start point
            FloatArray pos( bg->giveVertex(1) );
            pos.add(iTipProp.mPropagationLength, iTipProp.mPropagationDir);
            bg->insertVertexFront(pos);
        } else if ( iTipProp.mTipIndex == 1 ) {
            // Propagate end point
            FloatArray pos( bg->giveVertex( bg->giveNrVertices() ) );
            pos.add(iTipProp.mPropagationLength, iTipProp.mPropagationDir);
            bg->insertVertexBack(pos);
        }

    return true;
}

void EDCrack :: cropPolygon(const double &iArcPosStart, const double &iArcPosEnd)
{
    PolygonLine *pl = dynamic_cast<PolygonLine*> (bg);

    if(pl == NULL) {
        OOFEM_ERROR("Failed to cast bg to PolygonLine.")
    }

    std::vector<FloatArray> points;
    pl->giveSubPolygon(points, iArcPosStart, iArcPosEnd);
    pl->setVertices(points);

    const double tol2 = 1.0e-18;
    pl->removeDuplicatePoints(tol2);

}

void DofManList :: computeSurfaceNormalSignDist(double &oDist, const FloatArray &iPoint) const
{
    oDist = iPoint.at(3) - this->xi;     // will only work for plane el
}

IRResultType DofManList :: initializeFrom(InputRecord *ir)
{
    IRResultType result;     // Required by IR_GIVE_FIELD macro

    IntArray idList;
    IR_GIVE_FIELD(ir, idList, _IFT_DofManList_list);
    for ( int i = 1; i <= idList.giveSize(); i++ ) {
        this->dofManList.push_back( idList.at(i) );
    }

    std :: sort( dofManList.begin(), this->dofManList.end() );
    //IR_GIVE_FIELD(ir, this->xi, _IFT_DofManList_DelaminationLevel);

    return IRRT_OK;
}

int
DofManList::instanciateYourself( Domain *d )
{

    // Set the nodes based on a given node set
    if(this->setNumber > 0) {

        Set *set = d->giveSet( this->setNumber );
        const IntArray &nodes = set->giveNodeList( );
        //nodes.printYourself();
        for(int i = 1; i <= nodes.giveSize( ); i++) {
            this->dofManList.push_back( nodes.at( i ) );
        }
    }

    std::sort( dofManList.begin( ), this->dofManList.end( ) );

    return 1;
}

void DofManList :: giveInputRecord(DynamicInputRecord &input)
{
    input.setRecordKeywordField(this->giveInputRecordName(), 1);

    IntArray idList;
    idList.resize( dofManList.size() );
    for ( size_t i = 0; i < dofManList.size(); i++ ) {
        idList.at(i + 1) = dofManList [ i ];
    }

    input.setField(idList, _IFT_DofManList_list);
}

void DofManList :: giveBoundingSphere(FloatArray &oCenter, double &oRadius)
{
    // TODO: Compute tighter bounds. /ES
    oCenter = {
        0.0, 0.0
    };
    oRadius = std :: numeric_limits< double > :: max();
}


void WholeDomain :: giveInputRecord(DynamicInputRecord &input)
{
    input.setRecordKeywordField(this->giveInputRecordName(), 1);
}

void WholeDomain :: giveBoundingSphere(FloatArray &oCenter, double &oRadius)
{
    // TODO: Compute tighter bounds. /ES
    oCenter = {
        0.0, 0.0
    };
    oRadius = std :: numeric_limits< double > :: max();
}
} // end namespace oofem
