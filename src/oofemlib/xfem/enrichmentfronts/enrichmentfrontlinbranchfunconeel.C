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

#include "enrichmentfrontlinbranchfunconeel.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
#include "xfem/xfemmanager.h"
#include "domain.h"
#include "connectivitytable.h"
#include "spatiallocalizer.h"
#include "element.h"
#include "gausspoint.h"

namespace oofem {
REGISTER_EnrichmentFront(EnrFrontLinearBranchFuncOneEl)


EnrFrontLinearBranchFuncOneEl :: EnrFrontLinearBranchFuncOneEl()
{
    mpBranchFunc = new LinElBranchFunction();
}

EnrFrontLinearBranchFuncOneEl :: ~EnrFrontLinearBranchFuncOneEl()
{
    if ( mpBranchFunc != NULL ) {
        delete mpBranchFunc;
        mpBranchFunc = NULL;
    }
}


void EnrFrontLinearBranchFuncOneEl :: MarkNodesAsFront(std :: unordered_map< int, NodeEnrichmentType > &ioNodeEnrMarkerMap, XfemManager &ixFemMan,  const std :: unordered_map< int, double > &iLevelSetNormalDirMap, const std :: unordered_map< int, double > &iLevelSetTangDirMap, const TipInfo &iTipInfo)
{
    MarkTipElementNodesAsFront(ioNodeEnrMarkerMap, ixFemMan, iLevelSetNormalDirMap, iLevelSetTangDirMap, iTipInfo);
}

int EnrFrontLinearBranchFuncOneEl :: giveNumEnrichments(const DofManager &iDMan) const
{
    return 4;
}

void EnrFrontLinearBranchFuncOneEl :: evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iPos, const double &iLevelSet, int iNodeInd) const
{
    FloatArray xTip = { mTipInfo.mGlobalCoord.at(1), mTipInfo.mGlobalCoord.at(2) };

    FloatArray pos = { iPos.at(1), iPos.at(2) };

    // Crack tip tangent and normal
    const FloatArray &t = mTipInfo.mTangDir;
    const FloatArray &n = mTipInfo.mNormalDir;

    double r = 0.0, theta = 0.0;
    EnrichmentItem :: calcPolarCoord(r, theta, xTip, pos, n, t);

    mpBranchFunc->evaluateEnrFuncAt(oEnrFunc, r, theta);
}

void EnrFrontLinearBranchFuncOneEl :: evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet, int iNodeInd) const
{
    const FloatArray &xTip = mTipInfo.mGlobalCoord;

    // Crack tip tangent and normal
    const FloatArray &t = mTipInfo.mTangDir;
    const FloatArray &n = mTipInfo.mNormalDir;

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

void EnrFrontLinearBranchFuncOneEl :: evaluateEnrFuncJumps(std :: vector< double > &oEnrFuncJumps, GaussPoint &iGP, int iNodeInd, bool iGPLivesOnCurrentCrack, const double &iNormalSignDist) const
{
	const FloatArray &xTip = mTipInfo.mGlobalCoord;
	const FloatArray &gpCoord = *(iGP.giveCoordinates());

	double radius = gpCoord.distance(xTip);

	std :: vector< double > jumps;
	mpBranchFunc->giveJump(jumps, radius);

	oEnrFuncJumps.insert( oEnrFuncJumps.end(), jumps.begin(), jumps.end() );
}

IRResultType EnrFrontLinearBranchFuncOneEl :: initializeFrom(InputRecord *ir)
{
    return IRRT_OK;
}

void EnrFrontLinearBranchFuncOneEl :: giveInputRecord(DynamicInputRecord &input)
{
    int number = 1;
    input.setRecordKeywordField(this->giveInputRecordName(), number);
}
} // end namespace oofem
