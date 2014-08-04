/*
 * prescribedgradientbcneumann.C
 *
 *  Created on: Mar 5, 2014
 *      Author: svennine
 */

#include "prescribedgradientbcneumann.h"
#include "classfactory.h"
#include "node.h"
#include "masterdof.h"
#include "element.h"
#include "feinterpol.h"
#include "feinterpol2d.h"
#include "gausspoint.h"
#include "sparsemtrx.h"
#include "xfem/xfemelementinterface.h"
#include "xfem/integrationrules/discsegintegrationrule.h"

#include <cmath>

namespace oofem {
REGISTER_BoundaryCondition(PrescribedGradientBCNeumann);

PrescribedGradientBCNeumann::PrescribedGradientBCNeumann(int n, Domain * d):
PrescribedGradientBCWeakPeriodic(n, d)
{

}

PrescribedGradientBCNeumann::~PrescribedGradientBCNeumann()
{

}

IRResultType PrescribedGradientBCNeumann :: initializeFrom(InputRecord *ir)
{
    IRResultType result = PrescribedGradientBC :: initializeFrom(ir);

    mTractionInterpOrder = 1;
    mNumTractionNodesAtIntersections = 0;
    mTractionNodeSpacing = std::numeric_limits<int>::max();
    mDuplicateCornerNodes = true;
    mMeshIsPeriodic = true;

    return result;
}

} /* namespace oofem */
