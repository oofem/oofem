#include "FsbLinearStatic2D.h"
#include "Materials/structuralmaterial.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "engngm.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "fei2dquadlin.h"
#include "classfactory.h"
#include "elementinternaldofman.h"
#include "../sm/CrossSections/structuralcrosssection.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {

FsbLinearStatic2D::FsbLinearStatic2D(int n, Domain *d):
    FsbLinearStatic(n, d)
{

}

double FsbLinearStatic2D::computeVolumeAround(GaussPoint *gp)
{
    // Computes the volume element dV associated with the given gp.

    double weight = gp->giveWeight();
    const FloatArray &lCoords = gp->giveNaturalCoordinates(); // local/natural coords of the gp (parent domain)
    double detJ = fabs( this->giveInterpolation()->giveTransformationJacobian( lCoords, FEIElementGeometryWrapper(this) ) );
    double thickness = this->giveCrossSection()->give(CS_Thickness, gp); // the cross section keeps track of the thickness

    return detJ * thickness * weight; // dV
}

} // end namespace oofem
