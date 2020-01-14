#include "FsbPlate.h"
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
REGISTER_Element(FsbPlate);

FEI2dQuadLin FsbPlate :: interpolator(1, 2);

FsbPlate::FsbPlate(int n, Domain *d):
    FsbLinearStatic2D(n, d)
{
    numberOfDofs = 3;
    numberOfDofMans = 4;
}

FEInterpolation *FsbPlate::giveInterpolation() const
{
    return & interpolator;
}

IRResultType FsbPlate::initializeFrom(InputRecord *ir)
{
    // REIMPLEMENT if neccessary!
    numberOfGaussPoints = 4;
    IRResultType result = StructuralElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    if ( numberOfGaussPoints != 1 && numberOfGaussPoints != 4 && numberOfGaussPoints != 9 && numberOfGaussPoints != 16 && numberOfGaussPoints != 25 ) {
        numberOfGaussPoints = 4;
        OOFEM_WARNING("Number of Gauss points enforced to 4");
    }

    return IRRT_OK;
}

void FsbPlate::updateLocalNumbering(EntityRenumberingFunctor &f)
{
    // REIMPLEMENT if neccessary!
    StructuralElement :: updateLocalNumbering(f);
    if ( this->referenceNode ) {
        this->referenceNode = f(this->referenceNode, ERS_DofManager);
    }
}

// This method might not work properly because interpolator.evalN(N, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this)); is in comments
void FsbPlate::computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
{
    FloatMatrix ndx;
    interpolator.evaldNdx(ndx, gp->giveNaturalCoordinates(),
        FEIElementGeometryWrapper(this));
    //int colN = ndx.giveNumberOfRows();
    answer.resize(3, 12);

    FloatMatrix cosMat(4, 4);
    FloatMatrix sinMat(4, 4);
    FloatMatrix lMat(4, 4);
    FloatArray a(4);
    FloatArray b(4);
    FloatArray c(4);
    FloatArray d(4);
    FloatArray e(4);

    int i1[4] = {1, 2, 3, 4};
    int j1[4] = {2, 3, 4, 1};

    for (int it = 0; it < 4; it++)
    {
        double x = this->giveNode(j1[it])->giveCoordinates()->at(1) -
                this->giveNode(i1[it])->giveCoordinates()->at(1);

        double y = this->giveNode(j1[it])->giveCoordinates()->at(2) -
                this->giveNode(i1[it])->giveCoordinates()->at(2);

        lMat.at(i1[it], j1[it]) = sqrt(x*x + y*y);

        cosMat.at(i1[it], j1[it]) = y/lMat.at(i1[it], j1[it]);
        sinMat.at(i1[it], j1[it]) = -x/lMat.at(i1[it], j1[it]);

        double lTemp = 1/lMat.at(i1[it], j1[it]);

        c.at(i1[it]) = -x*x*lTemp;
        e.at(i1[it]) = -y*y*lTemp;

        lTemp /= lMat.at(i1[it], j1[it]);

        a.at(i1[it]) = -0.75*x*y*lTemp;
        b.at(i1[it]) = (0.5*y*y - 0.25*x*x)*lTemp;
        d.at(i1[it]) = (0.5*x*x - 0.25*y*y)*lTemp;
    }

    FloatMatrix dH(12, 4); // [dvx/dx, dvy/dx, dvx/dy, dvy/dy]
    FloatMatrix dNdx;
    // The following line is in comments because an error appears saying that N is not defined.
	//interpolator.evalN(N, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this));
    interpolator.evaldNdx(dNdx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this));

    for (int i = 1; i <= 2; i++)
    {
        dH.at(1, 2*i-1)  = 1.5*(c.at(5)*dNdx.at(5, i) - c.at(8)*dNdx.at(8, i));
        dH.at(2, 2*i-1)  = a.at(5)*dNdx.at(5, i) + a.at(8)*dNdx.at(8, i);
        dH.at(3, 2*i-1)  = -dNdx.at(1, i) - b.at(5)*dNdx.at(5, i) - b.at(8)*dNdx.at(8, i);
        dH.at(4, 2*i-1)  = 1.5*(c.at(6)*dNdx.at(6, i) - c.at(5)*dNdx.at(5, i));
        dH.at(5, 2*i-1)  = a.at(6)*dNdx.at(6, i) + a.at(5)*dNdx.at(5, i);
        dH.at(6, 2*i-1)  = -dNdx.at(2, i) - b.at(6)*dNdx.at(6, i) - b.at(5)*dNdx.at(5, i);
        dH.at(7, 2*i-1)  = 1.5*(c.at(7)*dNdx.at(7, i) - c.at(6)*dNdx.at(6, i));
        dH.at(8, 2*i-1)  = a.at(7)*dNdx.at(7, i) + a.at(6)*dNdx.at(6, i);
        dH.at(9, 2*i-1)  = -dNdx.at(3, i) - b.at(7)*dNdx.at(7, i) - b.at(6)*dNdx.at(6, i);
        dH.at(10, 2*i-1) = 1.5*(c.at(8)*dNdx.at(8, i) - c.at(7)*dNdx.at(7, i));
        dH.at(11, 2*i-1) = a.at(8)*dNdx.at(8, i) + a.at(7)*dNdx.at(7, i);
        dH.at(12, 2*i-1) = -dNdx.at(4, i) - b.at(8)*dNdx.at(8, i) - b.at(7)*dNdx.at(7, i);

        dH.at(1, 2*i)  = 1.5*(e.at(5)*dNdx.at(5, i) - e.at(8)*dNdx.at(8, i));
        dH.at(2, 2*i)  = dNdx.at(1, i) + d.at(5)*dNdx.at(5, i) + d.at(8)*dNdx.at(8, i);
        dH.at(3, 2*i)  = -dH.at(2, 1);
        dH.at(4, 2*i)  = 1.5*(e.at(6)*dNdx.at(6, i) - e.at(5)*dNdx.at(5, i));
        dH.at(5, 2*i)  = dNdx.at(2, i) + d.at(6)*dNdx.at(6, i) + d.at(5)*dNdx.at(5, i);
        dH.at(6, 2*i)  = -dH.at(5, 1);
        dH.at(7, 2*i)  = 1.5*(e.at(7)*dNdx.at(7, i) - e.at(6)*dNdx.at(6, i));
        dH.at(8, 2*i)  = dNdx.at(3, i) + d.at(7)*dNdx.at(7, i) + d.at(6)*dNdx.at(6, i);
        dH.at(9, 2*i)  = -dH.at(8, 1);
        dH.at(10, 2*i) = 1.5*(e.at(8)*dNdx.at(8, i) - e.at(7)*dNdx.at(7, i));
        dH.at(11, 2*i) = dNdx.at(4, i) + d.at(8)*dNdx.at(8, i) + d.at(7)*dNdx.at(7, i);
        dH.at(12, 2*i) = -dH.at(11, 1);
    }

    for (int i = 1; i <= 12; i++)
    {
        answer.at(1, i) = -dH.at(i, 1);
        answer.at(2, i) = -dH.at(i, 4);
        answer.at(3, i) = -dH.at(i, 2) - dH.at(i, 3);
    }
}

void FsbPlate::computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
{
    // a, b, c, d, e - should be initialized in some common space
    /*
    H.at(1, 1)  = 1.5*(c.at(5)*N.at(5) - c.at(8)*N.at(8));
    H.at(2, 1)  = a.at(5)*N.at(5) + a.at(8)*N.at(8);
    H.at(3, 1)  = -N.at(1) - b.at(5)*N.at(5) - b.at(8)*N.at(8);
    H.at(4, 1)  = 1.5*(c.at(6)*N.at(6) - c.at(5)*N.at(5));
    H.at(5, 1)  = a.at(6)*N.at(6) + a.at(5)*N.at(5);
    H.at(6, 1)  = -N.at(2) - b.at(6)*N.at(6) - b.at(5)*N.at(5);
    H.at(7, 1)  = 1.5*(c.at(7)*N.at(7) - c.at(6)*N.at(6));
    H.at(8, 1)  = a.at(7)*N.at(7) + a.at(6)*N.at(6);
    H.at(9, 1)  = -N.at(3) - b.at(7)*N.at(7) - b.at(6)*N.at(6);
    H.at(10, 1) = 1.5*(c.at(8)*N.at(8) - c.at(7)*N.at(7));
    H.at(11, 1) = a.at(8)*N.at(8) + a.at(7)*N.at(7);
    H.at(12, 1) = -N.at(4) - b.at(8)*N.at(8) - b.at(7)*N.at(7);

    H.at(1, 2)  = 1.5*(e.at(5)*N.at(5) - e.at(8)*N.at(8));
    H.at(2, 2)  = N.at(1) + d.at(5)*N.at(5) + d.at(8)*N.at(8);
    H.at(3, 2)  = -H.at(2, 1);
    H.at(4, 2)  = 1.5*(e.at(6)*N.at(6) - e.at(5)*N.at(5));
    H.at(5, 2)  = N.at(2) + d.at(6)*N.at(6) + d.at(5)*N.at(5);
    H.at(6, 2)  = -H.at(5, 1);
    H.at(7, 2)  = 1.5*(e.at(7)*N.at(7) - e.at(6)*N.at(6));
    H.at(8, 2)  = N.at(3) + d.at(7)*N.at(7) + d.at(6)*N.at(6);
    H.at(9, 2)  = -H.at(8, 1);
    H.at(10, 2) = 1.5*(e.at(8)*N.at(8) - e.at(7)*N.at(7));
    H.at(11, 2) = N.at(4) + d.at(8)*N.at(8) + d.at(7)*N.at(7);
    H.at(12, 2) = -H.at(11, 1);
    */
}

bool FsbPlate::computeGtoLRotationMatrix(FloatMatrix &answer)
{
    // From basiclsr element
    this->computeLCS();
    answer.resize(8, 8);
    answer.zero();

    int nodeNum = dofManArray.giveSize();
    for ( int i = 0; i < nodeNum; i++ ) { // Loops over nodes
        // In each node, transform global c.s. {D_u, D_v} into local c.s.
        answer.setSubMatrix(this->lcsMatrix, 1 + i * 2, 1 + i * 2);     // Displacements
    }

    return true;
}

} // end namespace oofem
