#ifndef FsbPlate_H
#define FsbPlate_H

#include "../sm/Elements/structuralelement.h"
#include "dofmanager.h"
#include "feinterpol.h"
#include "FsbLinearStatic2D.h"

///@name Input fields for FsbPlate
//@{
#define _IFT_FsbPlate_Name "FsbPlate"
#define _IFT_FsbPlate_dofstocondense "dofstocondense"
#define _IFT_FsbPlate_refnode "refnode"
#define _IFT_FsbPlate_refangle "refangle"
#define _IFT_FsbPlate_zaxis "zaxis"
//@}

namespace oofem {

class FsbPlate : public FsbLinearStatic2D
{
protected:
    /// Geometry interpolator only.
    static FEI2dQuadLin interpolator;

    virtual MaterialMode giveMaterialMode() { return _2dPlate; }
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li = 1, int ui = ALL_STRAINS);
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);

public:
    FsbPlate(int n, Domain *d);
    virtual ~FsbPlate() {}

    virtual FEInterpolation *giveInterpolation() const;
    virtual const char *giveClassName() const { return "FsbPlate"; }
    virtual const char *giveInputRecordName() const { return _IFT_FsbPlate_Name; }

    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const { answer = {D_w, R_u, R_v}; }

    virtual void initializeFrom(InputRecord &ir);
    virtual void updateLocalNumbering(EntityRenumberingFunctor &f);
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
};

} // end namespace oofem

#endif // FsbPlate_H
