/*
 * prescribedgradientbcweakdirichlet.C
 *
 *  Created on: May 19, 2014
 *      Author: svennine
 */

#include "prescribedgradientbcweakdirichlet.h"

#include "classfactory.h"

namespace oofem {
REGISTER_BoundaryCondition(PrescribedGradientBCWeakDirichlet);

PrescribedGradientBCWeakDirichlet :: PrescribedGradientBCWeakDirichlet(int n, Domain *d) :
    PrescribedGradientBCWeak(n, d)
{
    // TODO Auto-generated constructor stub
}

PrescribedGradientBCWeakDirichlet :: ~PrescribedGradientBCWeakDirichlet()
{
    // TODO Auto-generated destructor stub
}

IRResultType PrescribedGradientBCWeakDirichlet :: initializeFrom(InputRecord *ir)
{
    mMeshIsPeriodic = false;

    return PrescribedGradientBCWeak :: initializeFrom(ir);
}

void PrescribedGradientBCWeakDirichlet :: postInitialize()
{
    bool enforceCornerPeriodicity = false;
    int numSides = 4;
    createTractionMesh(enforceCornerPeriodicity, numSides);
}

void PrescribedGradientBCWeakDirichlet :: giveBoundaryCoordVector(FloatArray &oX, const FloatArray &iPos) const
{
    oX = {
        iPos [ 0 ] - mCenterCoord [ 0 ], iPos [ 1 ] - mCenterCoord [ 1 ]
    };
}

void PrescribedGradientBCWeakDirichlet :: checkIfCorner(bool &oIsCorner, bool &oDuplicatable, const FloatArray &iPos, const double &iNodeDistTol) const
{
    oIsCorner = false;
    oDuplicatable = false;

    FloatArray cornerPos = mLC;
    if ( iPos.distance(cornerPos) < iNodeDistTol ) {
        oIsCorner = true;
        oDuplicatable = true;
    }

    cornerPos = {
        mUC [ 0 ], mLC [ 1 ]
    };
    if ( iPos.distance(cornerPos) < iNodeDistTol ) {
        oIsCorner = true;
        oDuplicatable = true;
    }

    cornerPos = {
        mUC [ 0 ], mUC [ 1 ]
    };
    if ( iPos.distance(cornerPos) < iNodeDistTol ) {
        oIsCorner = true;
        oDuplicatable = true;
    }

    cornerPos = {
        mLC [ 0 ], mUC [ 1 ]
    };
    if ( iPos.distance(cornerPos) < iNodeDistTol ) {
        oIsCorner = true;
        oDuplicatable = true;
    }
}
} /* namespace oofem */
