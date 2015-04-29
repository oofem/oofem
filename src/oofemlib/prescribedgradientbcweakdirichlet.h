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

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void postInitialize();

    virtual const char *giveClassName() const { return "PrescribedGradientBCWeakDirichlet"; }
    virtual const char *giveInputRecordName() const { return _IFT_PrescribedGradientBCWeakDirichlet_Name; }

protected:
    virtual void giveBoundaryCoordVector(FloatArray &oX, const FloatArray &iPos) const;
    virtual void checkIfCorner(bool &oIsCorner, bool &oDuplicatable, const FloatArray &iPos, const double &iNodeDistTol) const;
    virtual bool boundaryPointIsOnActiveBoundary(const FloatArray &iPos) const { return true; }
};
} /* namespace oofem */

#endif /* PRESCRIBEDGRADIENTBCWEAKDIRICHLET_H_ */
