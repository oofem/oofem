#ifndef FSBLINEARSTATIC2D_H
#define FSBLINEARSTATIC2D_H

#include "FsbLinearStatic.h"
#include "dofmanager.h"
#include "feinterpol.h"

namespace oofem {

class FsbLinearStatic2D : public FsbLinearStatic
{
public:
    FsbLinearStatic2D(int n, Domain *d);
    virtual ~FsbLinearStatic2D() {}

    virtual double computeVolumeAround(GaussPoint *gp);
    virtual const char *giveClassName() const = 0;
};

} // end namespace oofem

#endif // FSBLINEARSTATIC2D_H
