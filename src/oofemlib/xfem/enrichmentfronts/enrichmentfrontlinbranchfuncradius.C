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


#include "enrichmentfrontlinbranchfuncradius.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
#include "xfem/xfemmanager.h"
#include "domain.h"
#include "connectivitytable.h"
#include "gausspoint.h"

namespace oofem {
REGISTER_EnrichmentFront(EnrFrontLinearBranchFuncRadius)

EnrFrontLinearBranchFuncRadius :: EnrFrontLinearBranchFuncRadius() :
    mEnrichmentRadius(0.0)
{
    mpBranchFunc = new LinElBranchFunction();
}

EnrFrontLinearBranchFuncRadius :: ~EnrFrontLinearBranchFuncRadius()
{
    if ( mpBranchFunc != NULL ) {
        delete mpBranchFunc;
        mpBranchFunc = NULL;
    }
}

void EnrFrontLinearBranchFuncRadius :: MarkNodesAsFront(std :: unordered_map< int, int > &ioNodeEnrMarkerMap, XfemManager &ixFemMan, const std :: unordered_map< int, double > &iLevelSetNormalDirMap, const std :: unordered_map< int, double > &iLevelSetTangDirMap, const std :: vector< TipInfo > &iTipInfo)
{
    // Enrich all nodes within a prescribed radius around the crack tips.
    // TODO: If performance turns out to be an issue, we may wish
    // to put the nodes in a Kd tree (or similar) to speed up searching.
    // For now, loop over all nodes.

    mTipInfo = iTipInfo;
    mNodeTipIndices.clear();

    Domain *d = ixFemMan.giveDomain();
    int nNodes = d->giveNumberOfDofManagers();

    for ( int i = 1; i <= nNodes; i++ ) {
        DofManager *dMan = d->giveDofManager(i);
        const FloatArray &nodePos = * ( dMan->giveCoordinates() );

        for ( int j = 0; j < int ( iTipInfo.size() ); j++ ) {
            double radius2 = iTipInfo [ j ].mGlobalCoord.distance_square(nodePos);

            if ( radius2 < mEnrichmentRadius * mEnrichmentRadius ) {
                ioNodeEnrMarkerMap [ i ] = 2;
                addTipIndexToNode(i, j);
            }
        }
    }
}

int EnrFrontLinearBranchFuncRadius :: giveNumEnrichments(const DofManager &iDMan) const
{
    std :: vector< int >tipIndices;
    int nodeInd = iDMan.giveGlobalNumber();
    giveNodeTipIndices(nodeInd, tipIndices);

    return 4 * tipIndices.size();
}

void EnrFrontLinearBranchFuncRadius :: evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iPos, const double &iLevelSet, int iNodeInd) const
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
        const FloatArray &t = mTipInfo [ tipInd ].mTangDir;
        const FloatArray &n = mTipInfo [ tipInd ].mNormalDir;

        double r = 0.0, theta = 0.0;
        EnrichmentItem :: calcPolarCoord(r, theta, xTip, pos, n, t);

        mpBranchFunc->evaluateEnrFuncAt(oEnrFunc, r, theta);
    }
}

void EnrFrontLinearBranchFuncRadius :: evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet, int iNodeInd) const
{
    oEnrFuncDeriv.clear();

    std :: vector< int >tipIndices;
    giveNodeTipIndices(iNodeInd, tipIndices);

    for ( size_t i = 0; i < tipIndices.size(); i++ ) {
        int tipInd = tipIndices [ i ];
        const FloatArray &xTip = mTipInfo [ tipInd ].mGlobalCoord;

        // Crack tip tangent and normal
        const FloatArray &t = mTipInfo [ tipInd ].mTangDir;
        const FloatArray &n = mTipInfo [ tipInd ].mNormalDir;

        double r = 0.0, theta = 0.0;
        EnrichmentItem :: calcPolarCoord(r, theta, xTip, iPos, n, t);


        size_t sizeStart = oEnrFuncDeriv.size();
        mpBranchFunc->evaluateEnrFuncDerivAt(oEnrFuncDeriv, r, theta);

        /**
         * Transform to global coordinates.
         */
        FloatMatrix E;
        E.resize(2, 2);
        E.setColumn(t, 1);
        E.setColumn(n, 2);


        for ( size_t j = sizeStart; j < oEnrFuncDeriv.size(); j++ ) {
            FloatArray enrFuncDerivGlob;
            enrFuncDerivGlob.beProductOf(E, oEnrFuncDeriv [ j ]);
            oEnrFuncDeriv [ j ] = enrFuncDerivGlob;
        }
    }
}

void EnrFrontLinearBranchFuncRadius :: evaluateEnrFuncJumps(std :: vector< double > &oEnrFuncJumps, GaussPoint &iGP, int iNodeInd) const
{
	oEnrFuncJumps.clear();

    std :: vector< int >tipIndices;
    giveNodeTipIndices(iNodeInd, tipIndices);

    for ( size_t i = 0; i < tipIndices.size(); i++ ) {
        int tipInd = tipIndices [ i ];
        const FloatArray &xTip = mTipInfo [ tipInd ].mGlobalCoord;
        const FloatArray &gpCoord = *(iGP.giveCoordinates());
        double radius = gpCoord.distance(xTip);

        std :: vector< double > jumps;
        mpBranchFunc->giveJump(jumps, radius);

        oEnrFuncJumps.insert( oEnrFuncJumps.end(), jumps.begin(), jumps.end() );
    }
}

IRResultType EnrFrontLinearBranchFuncRadius :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    IR_GIVE_FIELD(ir, mEnrichmentRadius, _IFT_EnrFrontLinearBranchFuncRadius_Radius);

    return IRRT_OK;
}

void EnrFrontLinearBranchFuncRadius :: giveInputRecord(DynamicInputRecord &input)
{
    int number = 1;
    input.setRecordKeywordField(this->giveInputRecordName(), number);

    input.setField(mEnrichmentRadius, _IFT_EnrFrontLinearBranchFuncRadius_Radius);
}
} // end namespace oofem