/*
 * enrichmentfrontcohesivebranchfunconeel.C
 *
 *  Created on: Nov 28, 2014
 *      Author: svennine
 */

#include "enrichmentfrontcohesivebranchfunconeel.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"
#include "xfem/xfemmanager.h"
#include "domain.h"
#include "connectivitytable.h"
#include "spatiallocalizer.h"
#include "element.h"
#include "gausspoint.h"

namespace oofem {
REGISTER_EnrichmentFront(EnrFrontCohesiveBranchFuncOneEl)

EnrFrontCohesiveBranchFuncOneEl::EnrFrontCohesiveBranchFuncOneEl()
{
    mpBranchFunc = new CohesiveBranchFunction();
}

EnrFrontCohesiveBranchFuncOneEl::~EnrFrontCohesiveBranchFuncOneEl()
{
    delete mpBranchFunc;
}


void EnrFrontCohesiveBranchFuncOneEl :: MarkNodesAsFront(std :: unordered_map< int, NodeEnrichmentType > &ioNodeEnrMarkerMap, XfemManager &ixFemMan,  const std :: unordered_map< int, double > &iLevelSetNormalDirMap, const std :: unordered_map< int, double > &iLevelSetTangDirMap, const TipInfo &iTipInfo)
{
    MarkTipElementNodesAsFront(ioNodeEnrMarkerMap, ixFemMan, iLevelSetNormalDirMap, iLevelSetTangDirMap, iTipInfo);
}

int EnrFrontCohesiveBranchFuncOneEl :: giveNumEnrichments(const DofManager &iDMan) const
{
    return 1;
}

void EnrFrontCohesiveBranchFuncOneEl :: evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const EfInput &iEfInput) const
{
    FloatArray xTip = {
        mTipInfo.mGlobalCoord.at(1), mTipInfo.mGlobalCoord.at(2)
    };

    FloatArray pos = {
        iEfInput.mPos.at(1), iEfInput.mPos.at(2)
    };

    // Crack tangent and normal
    FloatArray t, n;
    bool flipTangent = false;
    computeCrackTangent(t, n, flipTangent, iEfInput);

    double r = 0.0, theta = 0.0;
    EnrichmentItem :: calcPolarCoord(r, theta, xTip, pos, n, t, iEfInput, flipTangent);

    mpBranchFunc->evaluateEnrFuncAt(oEnrFunc, r, theta);

#ifdef DEBUG
    for ( double val:oEnrFunc ) {
        if ( !std :: isfinite(val) ) {
            printf("r: %e theta: %e\n", r, theta);
            OOFEM_ERROR("!std::isfinite(val)")
        }
    }
#endif
}

void EnrFrontCohesiveBranchFuncOneEl :: evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const EfInput &iEfInput, const FloatArray &iGradLevelSet) const
{
    const FloatArray &xTip = mTipInfo.mGlobalCoord;

    // Crack tangent and normal
    FloatArray t, n;
    bool flipTangent = false;
    computeCrackTangent(t, n, flipTangent, iEfInput);

    double r = 0.0, theta = 0.0;
    EnrichmentItem :: calcPolarCoord(r, theta, xTip, iEfInput.mPos, n, t, iEfInput, flipTangent);


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

void EnrFrontCohesiveBranchFuncOneEl :: evaluateEnrFuncJumps(std :: vector< double > &oEnrFuncJumps, GaussPoint &iGP, int iNodeInd, bool iGPLivesOnCurrentCrack, const double &iNormalSignDist) const
{
    const FloatArray &xTip = mTipInfo.mGlobalCoord;
    const FloatArray &gpCoord = iGP.giveGlobalCoordinates();

    double radius = gpCoord.distance(xTip);

    std :: vector< double >jumps;
    mpBranchFunc->giveJump(jumps, radius);

    oEnrFuncJumps.insert( oEnrFuncJumps.end(), jumps.begin(), jumps.end() );
}

IRResultType EnrFrontCohesiveBranchFuncOneEl :: initializeFrom(InputRecord *ir)
{
    return IRRT_OK;
}

void EnrFrontCohesiveBranchFuncOneEl :: giveInputRecord(DynamicInputRecord &input)
{
    int number = 1;
    input.setRecordKeywordField(this->giveInputRecordName(), number);
}

} /* namespace oofem */
