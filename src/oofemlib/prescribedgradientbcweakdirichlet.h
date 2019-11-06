/*
 * prescribedgradientbcweakdirichlet.h
 *
 *  Created on: May 19, 2014
 *      Author: svennine
 */

#ifndef PRESCRIBEDGRADIENTBCWEAKDIRICHLET_H_
#define PRESCRIBEDGRADIENTBCWEAKDIRICHLET_H_

#include "prescribedgradientbcweak.h"

#define _IFT_PrescribedGradientBCWeakDirichlet_Name   "prescribedgradientbcweakdirichlet"

namespace oofem {
class PrescribedGradientBCWeakDirichlet : public PrescribedGradientBCWeak
{
public:
    PrescribedGradientBCWeakDirichlet(int n, Domain *d);
    virtual ~PrescribedGradientBCWeakDirichlet();

    void initializeFrom(InputRecord &ir) override;

    void postInitialize() override;

    const char *giveClassName() const override { return "PrescribedGradientBCWeakDirichlet"; }
    const char *giveInputRecordName() const override { return _IFT_PrescribedGradientBCWeakDirichlet_Name; }

protected:
    void giveBoundaryCoordVector(FloatArray &oX, const FloatArray &iPos) const override;
    void checkIfCorner(bool &oIsCorner, bool &oDuplicatable, const FloatArray &iPos, const double &iNodeDistTol) const override;
    bool boundaryPointIsOnActiveBoundary(const FloatArray &iPos) const override { return true; }
};
} /* namespace oofem */

#endif /* PRESCRIBEDGRADIENTBCWEAKDIRICHLET_H_ */
