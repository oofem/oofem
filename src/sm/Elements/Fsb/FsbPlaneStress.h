#ifndef FsbPlaneStress_H
#define FsbPlaneStress_H

#include "../sm/Elements/structuralelement.h"
#include "dofmanager.h"
#include "feinterpol.h"
#include "FsbLinearStatic2D.h"

///@name Input fields for FsbPlaneStress
//@{
#define _IFT_FsbPlaneStress_Name "FsbPlaneStress"
#define _IFT_FsbPlaneStress_dofstocondense "dofstocondense"
#define _IFT_FsbPlaneStress_refnode "refnode"
#define _IFT_FsbPlaneStress_refangle "refangle"
#define _IFT_FsbPlaneStress_zaxis "zaxis"
//@}

namespace oofem {

class FsbPlaneStress : public FsbLinearStatic2D
{
protected:
    /// Geometry interpolator only.
    static FEI2dQuadLin interpolator;

    virtual MaterialMode giveMaterialMode() { return _PlaneStress; }
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li = 1, int ui = ALL_STRAINS);

public:
    FsbPlaneStress(int n, Domain *d);
    virtual ~FsbPlaneStress() {}

    virtual FEInterpolation *giveInterpolation() const;
    virtual const char *giveClassName() const { return "FsbPlaneStress"; }
    virtual const char *giveInputRecordName() const { return _IFT_FsbPlaneStress_Name; }

    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const { answer = { D_u, D_v }; }

    virtual void initializeFrom(InputRecord &ir);
    virtual void updateLocalNumbering(EntityRenumberingFunctor &f);
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
};

} // end namespace oofem

#endif // FsbPlaneStress_H
