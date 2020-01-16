#include "FsbPlaneStressDrill.h"
#include "Materials/structuralmaterial.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "engngm.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "FsbInterpolatorPlaneStress.h"
#include "classfactory.h"
#include "elementinternaldofman.h"
#include "../sm/CrossSections/structuralcrosssection.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(FsbPlaneStressDrill);

FsbInterpolatorPlaneStress FsbPlaneStressDrill :: interpolator(1, 2);

FsbPlaneStressDrill::FsbPlaneStressDrill(int n, Domain *d):
    FsbLinearStatic2D(n, d)
{
    numberOfDofs = 3;
    numberOfDofMans = 4;
}

FEInterpolation *FsbPlaneStressDrill::giveInterpolation() const
{
    return & interpolator;
}

void FsbPlaneStressDrill::initializeFrom(InputRecord &ir)
{
    numberOfGaussPoints = 4;
    StructuralElement :: initializeFrom(ir);

    if ( numberOfGaussPoints != 1 && numberOfGaussPoints != 4 && numberOfGaussPoints != 9 && numberOfGaussPoints != 16 && numberOfGaussPoints != 25 ) {
        numberOfGaussPoints = 4;
        OOFEM_WARNING("Number of Gauss points enforced to 4");
    }
}

void FsbPlaneStressDrill::updateLocalNumbering(EntityRenumberingFunctor &f)
{
    StructuralElement :: updateLocalNumbering(f);
    if ( this->referenceNode ) {
        this->referenceNode = f(this->referenceNode, ERS_DofManager);
    }
}

void FsbPlaneStressDrill::computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
{
    FloatMatrix ndx;
    interpolator.evaldNdx(ndx, gp->giveNaturalCoordinates(),
        FEIElementGeometryWrapper(this));
    //int colN = ndx.giveNumberOfRows();
    answer.resize(3, 12); // size should be 3, 3*colN? // 3*4 = 12

    FloatMatrix cosMat(4, 4);
    FloatMatrix sinMat(4, 4);
    FloatMatrix lMat(4, 4);

    int i1[4] = {1, 2, 3, 4};
    int j1[4] = {2, 3, 4, 1};

    for (int it = 0; it < 4; it++)
    {
        double x = this->giveNode(j1[it])->giveCoordinates().at(1) -
                this->giveNode(i1[it])->giveCoordinates().at(1);

        double y = this->giveNode(j1[it])->giveCoordinates().at(2) -
                this->giveNode(i1[it])->giveCoordinates().at(2);

        lMat.at(i1[it], j1[it]) = sqrt(x*x + y*y);

        cosMat.at(i1[it], j1[it]) = y/lMat.at(i1[it], j1[it]);
        sinMat.at(i1[it], j1[it]) = -x/lMat.at(i1[it], j1[it]);
    }

    int i[4] = {1, 2, 3, 4};
    int m[4] = {5, 6, 7, 8};
    int l[4] = {8, 5, 6, 7};
    int j[4] = {4, 1, 2, 3};
    int k[4] = {2, 3, 4, 1};

    for ( int it = 0; it < 4; it++ ) {
        // N matrix part
        answer.at(1, 3*it + 1) = ndx.at(it + 1, 1);
        answer.at(1, 3*it + 2) = 0;

        answer.at(2, 3*it + 1) = 0;
        answer.at(2, 3*it + 2) = ndx.at(it + 1, 2);

        answer.at(3, 3*it + 1) = ndx.at(it + 1, 2);
        answer.at(3, 3*it + 2) = ndx.at(it + 1, 1);

        // G matrix part
        answer.at(1, 3*it + 3) =
            0.125*(lMat.at(i[it], j[it]) * cosMat.at(i[it], j[it]) * ndx.at(l[it],1) -
                   lMat.at(i[it], k[it]) * cosMat.at(i[it], k[it]) * ndx.at(m[it],2)
                  );

        answer.at(2, 3*it + 3) =
            0.125*(lMat.at(i[it], j[it]) * sinMat.at(i[it], j[it]) * ndx.at(l[it],1) -
                   lMat.at(i[it], k[it]) * sinMat.at(i[it], k[it]) * ndx.at(m[it],2)
                  );

        answer.at(3, 3*it + 3) =
            0.125*((lMat.at(i[it], j[it]) * cosMat.at(i[it], j[it]) * ndx.at(l[it],2) -
                    lMat.at(i[it], k[it]) * cosMat.at(i[it], k[it]) * ndx.at(m[it],2)
                   ) +
                   (lMat.at(i[it], j[it]) * sinMat.at(i[it], j[it]) * ndx.at(l[it],1) -
                    lMat.at(i[it], k[it]) * sinMat.at(i[it], k[it]) * ndx.at(m[it],1)
                   ));
    }

    int abcd = 0;
}

bool FsbPlaneStressDrill::computeGtoLRotationMatrix(FloatMatrix &answer)
{
    this->computeLCS();
    answer.resize(12, 12);
    answer.zero();

    for ( int i = 0; i <= 3; i++ ) { // Loops over nodes
        // In each node, transform global c.s. {D_u, D_v,D_w} into local c.s.
        answer.setSubMatrix(this->lcsMatrix, 1 + i * 3, 1 + i * 3);     // Displacements
    }

    return true;
}

int FsbPlaneStressDrill::computeLoadGToLRotationMtrx(FloatMatrix &answer)
// Returns the rotation matrix of the receiver of the size [5,6]
// f(local) = T * f(global)
{
    this->computeLCS();

    answer.resize(3, 3);
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        answer.at(1, i) = lcsMatrix.at(1, i);
        answer.at(2, i) = lcsMatrix.at(2, i);
        answer.at(3, i) = lcsMatrix.at(3, i);
    }

    return 1;
}

} // end namespace oofem
