/*
 * prescribedgradientbcweakperiodic.h
 *
 *  Created on: May 19, 2014
 *      Author: svennine
 */

#ifndef PRESCRIBEDGRADIENTBCWEAKPERIODIC_H_
#define PRESCRIBEDGRADIENTBCWEAKPERIODIC_H_

#include "prescribedgradientbcweak.h"

#define _IFT_PrescribedGradientBCWeakPeriodic_Name   "prescribedgradientbcweakperiodic"

namespace oofem {
class PrescribedGradientBCWeakPeriodic : public PrescribedGradientBCWeak
{
public:
    PrescribedGradientBCWeakPeriodic(int n, Domain *d);
    virtual ~PrescribedGradientBCWeakPeriodic();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void postInitialize();

    virtual const char *giveClassName() const { return "PrescribedGradientBCWeakPeriodic"; }
    virtual const char *giveInputRecordName() const { return _IFT_PrescribedGradientBCWeakPeriodic_Name; }

protected:
    virtual void giveBoundaryCoordVector(FloatArray &oX, const FloatArray &iPos) const;
    virtual void checkIfCorner(bool &oIsCorner, bool &oDuplicatable, const FloatArray &iPos, const double &iNodeDistTol) const;

    virtual bool boundaryPointIsOnActiveBoundary(const FloatArray &iPos) const { return pointIsOnGammaPlus(iPos); }
};
} /* namespace oofem */

#endif /* PRESCRIBEDGRADIENTBCWEAKPERIODIC_H_ */
