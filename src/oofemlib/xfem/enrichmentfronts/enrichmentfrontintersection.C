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

#include "enrichmentfrontintersection.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
#include "xfem/xfemmanager.h"
#include "domain.h"
#include "connectivitytable.h"
#include "spatiallocalizer.h"
#include "element.h"
#include "gausspoint.h"

namespace oofem {
REGISTER_EnrichmentFront(EnrFrontIntersection)

EnrFrontIntersection::EnrFrontIntersection() {

}

EnrFrontIntersection::~EnrFrontIntersection() {

}

void EnrFrontIntersection::MarkNodesAsFront(std :: unordered_map< int, int > &ioNodeEnrMarkerMap, XfemManager &ixFemMan,  const std :: unordered_map< int, double > &iLevelSetNormalDirMap, const std :: unordered_map< int, double > &iLevelSetTangDirMap, const std :: vector< TipInfo > &iTipInfo)
{
    MarkTipElementNodesAsFront(ioNodeEnrMarkerMap, ixFemMan, iLevelSetNormalDirMap, iLevelSetTangDirMap, iTipInfo);
}

int EnrFrontIntersection :: giveNumEnrichments(const DofManager &iDMan) const
{
    std :: vector< int >tipIndices;
    int nodeInd = iDMan.giveGlobalNumber();
    giveNodeTipIndices(nodeInd, tipIndices);

    return 1 * tipIndices.size();
}

void EnrFrontIntersection :: evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iPos, const double &iLevelSet, int iNodeInd) const
{
    oEnrFunc.clear();

    std :: vector< int >tipIndices;
    giveNodeTipIndices(iNodeInd, tipIndices);

    for ( size_t i = 0; i < tipIndices.size(); i++ ) {
        int tipInd = tipIndices [ i ];
        FloatArray xTip = {
            mTipInfo [ tipInd ].mGlobalCoord.at(1), mTipInfo [ tipInd ].mGlobalCoord.at(2)
        };

        FloatArray pos = {
            iPos.at(1), iPos.at(2)
        };

        // Crack tip tangent and normal
        // TODO: We may wish to have a user defined tangential direction here!
        const FloatArray &t = mTipInfo [ tipInd ].mTangDir;
        const FloatArray &n = mTipInfo [ tipInd ].mNormalDir;

        FloatArray tipToPos = { iPos(0)-xTip(0), iPos(1)-xTip(1) };

        // Heaviside in normal direction
        double Hn = 0.0;
        if( tipToPos.dotProduct(n) > 0.0 ) {
            Hn = 1.0;
        }

        // Heaviside in tangential direction
        double Ht = 0.0;
        if( tipToPos.dotProduct(t) < 0.0 ) {
            Ht = 1.0;
        }

        oEnrFunc.push_back(Hn*Ht);
    }
}

void EnrFrontIntersection :: evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet, int iNodeInd) const
{
    oEnrFuncDeriv.clear();

    std :: vector< int >tipIndices;
    giveNodeTipIndices(iNodeInd, tipIndices);

    for ( size_t i = 0; i < tipIndices.size(); i++ ) {
        FloatArray enrFuncDeriv = {0.0, 0.0};
        oEnrFuncDeriv.push_back(enrFuncDeriv);
    }

}

void EnrFrontIntersection :: evaluateEnrFuncJumps(std :: vector< double > &oEnrFuncJumps, GaussPoint &iGP, int iNodeInd) const
{
    // TODO: Implement
    printf("Warning:EnrFrontIntersection :: evaluateEnrFuncJumps() is not implemented yet.\n ");

    oEnrFuncJumps.clear();

    std :: vector< int >tipIndices;
    giveNodeTipIndices(iNodeInd, tipIndices);

    for ( size_t i = 0; i < tipIndices.size(); i++ ) {
        int tipInd = tipIndices [ i ];
        const FloatArray &xTip = mTipInfo [ tipInd ].mGlobalCoord;
        const FloatArray &gpCoord = *(iGP.giveCoordinates());

        double radius = gpCoord.distance(xTip);

        std :: vector< double > jumps;
        jumps.push_back(0.0);

        oEnrFuncJumps.insert( oEnrFuncJumps.end(), jumps.begin(), jumps.end() );
    }
}

IRResultType EnrFrontIntersection :: initializeFrom(InputRecord *ir)
{
    IRResultType result; // Required by IR_GIVE_FIELD macro
    IR_GIVE_FIELD(ir, mTangent, _IFT_EnrFrontIntersection_Tangent);
    printf("mTangent: "); mTangent.printYourself();

    return IRRT_OK;
}

void EnrFrontIntersection :: giveInputRecord(DynamicInputRecord &input)
{
    int number = 1;
    input.setRecordKeywordField(this->giveInputRecordName(), number);

    input.setField(mTangent, _IFT_EnrFrontIntersection_Tangent);
}


} /* namespace oofem */
