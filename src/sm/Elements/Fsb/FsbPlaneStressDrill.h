#ifndef FsbPlaneStressDrill_H
#define FsbPlaneStressDrill_H

#include "../sm/Elements/structuralelement.h"
#include "dofmanager.h"
#include "feinterpol.h"
#include "FsbPlaneStress.h"
#include "FsbInterpolatorPlaneStress.h"

///@name Input fields for FsbPlaneStressDrill
//@{
#define _IFT_FsbPlaneStressDrill_Name "FsbPlaneStressDrill"
#define _IFT_FsbPlaneStressDrill_dofstocondense "dofstocondense"
#define _IFT_FsbPlaneStressDrill_refnode "refnode"
#define _IFT_FsbPlaneStressDrill_refangle "refangle"
#define _IFT_FsbPlaneStressDrill_zaxis "zaxis"
//@}

namespace oofem {

class FsbPlaneStressDrill : public FsbLinearStatic2D
{
protected:

    /// Geometry interpolator only.
    static FsbInterpolatorPlaneStress interpolator;

    virtual MaterialMode giveMaterialMode() { return _PlaneStress; }
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li = 1, int ui = ALL_STRAINS);

public:
    FsbPlaneStressDrill(int n, Domain *d);
    virtual ~FsbPlaneStressDrill() {}

    virtual FEInterpolation *giveInterpolation() const;
    virtual const char *giveClassName() const { return "FsbPlaneStressDrill"; }
    virtual const char *giveInputRecordName() const { return _IFT_FsbPlaneStressDrill_Name; }

    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const { answer = { D_u, D_v, R_w }; }
    int computeNumberOfGlobalDofs() {
        return this->computeNumberOfDofs();
    }

    virtual void initializeFrom(InputRecord &ir);
    virtual void updateLocalNumbering(EntityRenumberingFunctor &f);
    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer);
    virtual int computeLoadGToLRotationMtrx(FloatMatrix &answer);
};

} // end namespace oofem

#endif // FsbPlaneStressDrill_H
