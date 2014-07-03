/*
 * prescribedgradientbcweakperiodic.C
 *
 *  Created on: May 19, 2014
 *      Author: svennine
 */

#include "prescribedgradientbcweakperiodic.h"

#include "classfactory.h"

namespace oofem {
REGISTER_BoundaryCondition(PrescribedGradientBCWeakPeriodic);

PrescribedGradientBCWeakPeriodic::PrescribedGradientBCWeakPeriodic(int n, Domain * d):
PrescribedGradientBCWeak(n,d)
{
	// TODO Auto-generated constructor stub

}

PrescribedGradientBCWeakPeriodic::~PrescribedGradientBCWeakPeriodic() {
	// TODO Auto-generated destructor stub
}

IRResultType PrescribedGradientBCWeakPeriodic :: initializeFrom(InputRecord *ir)
{
	PrescribedGradientBCWeak :: initializeFrom(ir);
	mMeshIsPeriodic = true;

    return IRRT_OK;
}

void PrescribedGradientBCWeakPeriodic :: postInitialize()
{
	bool enforceCornerPeriodicity = true;
    createTractionMesh(enforceCornerPeriodicity);
}

void PrescribedGradientBCWeakPeriodic :: giveBoundaryCoordVector(FloatArray &oX, const FloatArray &iPos) const
{
    FloatArray xMinus;
    giveMirroredPointOnGammaMinus(xMinus, iPos);

    oX = {iPos[0]-xMinus[0], iPos[1]-xMinus[1]};
}

void PrescribedGradientBCWeakPeriodic :: checkIfCorner(bool &oIsCorner, bool &oDuplicatable, const FloatArray &iPos, const double &iNodeDistTol) const
{
	oIsCorner = false;
	oDuplicatable = false;

    FloatArray cornerPos = mLC;
    if( iPos.distance(cornerPos) < iNodeDistTol ) {
    	oIsCorner = true;
    }

    cornerPos = {mUC[0], mLC[1]};
    if( iPos.distance( cornerPos ) < iNodeDistTol ) {
    	oIsCorner = true;
    }

    cornerPos = {mUC[0], mUC[1]};
    if( iPos.distance( cornerPos ) < iNodeDistTol ) {
    	oIsCorner = true;
    	oDuplicatable = true;
    }

    cornerPos = {mLC[0], mUC[1]};
    if( iPos.distance( cornerPos ) < iNodeDistTol ) {
    	oIsCorner = true;
    }
}

} /* namespace oofem */