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

    void initializeFrom(InputRecord &ir) override;

    void postInitialize() override;

    const char *giveClassName() const override { return "PrescribedGradientBCWeakPeriodic"; }
    const char *giveInputRecordName() const override { return _IFT_PrescribedGradientBCWeakPeriodic_Name; }

protected:
    void giveBoundaryCoordVector(FloatArray &oX, const FloatArray &iPos) const override;
    void checkIfCorner(bool &oIsCorner, bool &oDuplicatable, const FloatArray &iPos, const double &iNodeDistTol) const override;

    bool boundaryPointIsOnActiveBoundary(const FloatArray &iPos) const override { return pointIsOnGammaPlus(iPos); }
};
} /* namespace oofem */

#endif /* PRESCRIBEDGRADIENTBCWEAKPERIODIC_H_ */
