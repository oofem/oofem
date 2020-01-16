#include "FsbPlaneStress.h"
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
REGISTER_Element(FsbPlaneStress);

FEI2dQuadLin FsbPlaneStress :: interpolator(1, 2);

FsbPlaneStress::FsbPlaneStress(int n, Domain *d):
    FsbLinearStatic2D(n, d)
{
    numberOfDofs = 2;
    numberOfDofMans = 4;
}

FEInterpolation *FsbPlaneStress::giveInterpolation() const
{
    return & interpolator;
}

void FsbPlaneStress::initializeFrom(InputRecord &ir)
{
    numberOfGaussPoints = 4;
    StructuralElement :: initializeFrom(ir);

    if ( numberOfGaussPoints != 1 && numberOfGaussPoints != 4 && numberOfGaussPoints != 9 && numberOfGaussPoints != 16 && numberOfGaussPoints != 25 ) {
        numberOfGaussPoints = 4;
        OOFEM_WARNING("Number of Gauss points enforced to 4");
    }
}

void FsbPlaneStress::updateLocalNumbering(EntityRenumberingFunctor &f)
{
    StructuralElement :: updateLocalNumbering(f);
    if ( this->referenceNode ) {
        this->referenceNode = f(this->referenceNode, ERS_DofManager);
    }
}

void FsbPlaneStress::computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
{
    FloatMatrix ndx;
    interpolator.evaldNdx(ndx, gp->giveNaturalCoordinates(),
        FEIElementGeometryWrapper(this));
    int colN = ndx.giveNumberOfRows();
    answer.resize(3, 2*colN);

    for ( int i = 0; i < colN; i++ ) {
        answer.at(1, 2*i + 1) = ndx.at(i + 1, 1);
        answer.at(1, 2*i + 2) = 0;

        answer.at(2, 2*i + 1) = 0;
        answer.at(2, 2*i + 2) = ndx.at(i + 1, 2);

        answer.at(3, 2*i + 1) = ndx.at(i + 1, 2);
        answer.at(3, 2*i + 2) = ndx.at(i + 1, 1);
    }
}

bool FsbPlaneStress::computeGtoLRotationMatrix(FloatMatrix &answer)
{
    // From basiclsr element
    this->computeLCS();
    int dofNb = this->computeNumberOfDofs();
    answer.resize(dofNb, dofNb);
    answer.zero();

    for ( int i = 0; i < numberOfDofMans; i++ ) { // Loops over nodes
        // In each node, transform global c.s. {D_u, D_v} into local c.s.
        answer.setSubMatrix(this->lcsMatrix, 1 + i * 2, 1 + i * 2);     // Displacements
    }

    return true;
}

} // end namespace oofem
