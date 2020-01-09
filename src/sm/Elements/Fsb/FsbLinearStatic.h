#ifndef FsbLinearStatic_H
#define FsbLinearStatic_H

#include "../sm/Elements/structuralelement.h"
#include "dofmanager.h"
#include "feinterpol.h"

namespace oofem {

class FEI2dQuadLin;

class FsbLinearStatic : public StructuralElement
{
protected:
    int numberOfDofs;
    int referenceNode;
    FloatMatrix lcsMatrix;

    virtual void computeGaussPoints();
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li = 1, int ui = ALL_STRAINS) = 0;
//    This method is not used in computations as long as OOFEG
//    macro is not defined
//    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);

public:
    FsbLinearStatic(int n, Domain *d);
    virtual ~FsbLinearStatic() {}

    virtual FEInterpolation *giveInterpolation() const = 0;
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const = 0;
    virtual int giveNumberOfDofManagers() const { return numberOfDofMans; }
    virtual int computeNumberOfDofs() {
        return numberOfDofs*numberOfDofMans;
    }

    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    virtual bool computeGtoLRotationMatrix(FloatMatrix &answer) = 0;
    virtual void computeLCS();

    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    virtual double computeVolumeAround(GaussPoint *gp) = 0;
};

} // end namespace oofem

#endif // FsbLinearStatic_H
