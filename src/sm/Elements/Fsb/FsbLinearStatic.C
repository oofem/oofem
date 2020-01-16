#include "FsbLinearStatic.h"
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

FsbLinearStatic::FsbLinearStatic(int n, Domain *d) :
    StructuralElement(n, d),
    referenceNode(0)
{
    numberOfDofs = 0;
}

//void FsbPlaneStress::computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
//{
//    FloatArray h(4);
//    interp.evalN(h, iLocCoord,
//        FEIElementGeometryWrapper(this));
//    answer = h;
//}

void FsbLinearStatic::computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    // partialy copied from structuralelement
    double dV;
    FloatMatrix d, bj, dbj;
    bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);

    answer.clear();
    for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
        this->computeBmatrixAt(gp, bj);
        this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);
        dV = this->computeVolumeAround(gp);
        dbj.beProductOf(d, bj);
        if ( matStiffSymmFlag ) {
            answer.plusProductSymmUpper(bj, dbj, dV);
        } else {
            answer.plusProductUnsym(bj, dbj, dV);
        }
    }

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}

void FsbLinearStatic::computeLCS()
{
    // From basiclsr element
    int num = 2;
    lcsMatrix.resize(num, num); // Note! G -> L transformation matrix
    FloatArray e[3];

    // compute e1' = [N2-N1]  and  help = [N4-N1]
	// Pointer operators removed because the project could not build. Test later to check if the method still works properly.
    e[0].beDifferenceOf( this->giveNode(2)->giveCoordinates(), this->giveNode(1)->giveCoordinates() );
    e[1].beDifferenceOf( this->giveNode(4)->giveCoordinates(), this->giveNode(1)->giveCoordinates() );
    e[0].normalize();
    e[2].beVectorProductOf(e[0], e[1]);
    e[2].normalize();
    e[1].beVectorProductOf(e[2],e[0]);

    for ( int i = 1; i <= num; i++ ) {
        for ( int j = 1; j <= num; j++)
            this->lcsMatrix.at(j, i) = e[j].at(i);
    }
}

void FsbLinearStatic::computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    // copied from planestresselement (structural2delement.C)
    answer = this->giveStructuralCrossSection()->giveStiffnessMatrix_PlaneStress(rMode, gp, tStep);
}

void FsbLinearStatic::computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        // the gauss point is used only when methods from crosssection and/or material
        // classes are requested
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 2) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], this->numberOfGaussPoints, this);
    }
}

} // end namespace oofem
