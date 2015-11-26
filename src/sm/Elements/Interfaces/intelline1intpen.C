/*
 * intelline1intpen.C
 *
 *  Created on: Nov 20, 2015
 *      Author: svennine
 */

#include "intelline1intpen.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"

namespace oofem {
REGISTER_Element(IntElLine1IntPen);

IntElLine1IntPen::IntElLine1IntPen(int n, Domain * d) : IntElLine1(n, d) {

}

IntElLine1IntPen::~IntElLine1IntPen() {

}

#if 0
void
IntElLine1IntPen :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
	double xi_0 = 0.0;

    // Computes the stiffness matrix of the receiver K_cohesive = int_A ( N^t * dT/dj * N ) dA
    FloatMatrix N, D, DN;
    bool matStiffSymmFlag = this->giveCrossSection()->isCharacteristicMtrxSymmetric(rMode);
    answer.clear();

    FloatMatrix rotationMatGtoL;
    FloatArray u;

    // First loop over GP: compute projection of test function and traction.
    // The setting is as follows: we have an interface with quadratic interpolation and we
    // wish to project onto the space of constant functions on each element.

    // Projecting the basis functions gives a constant for each basis function.
    FloatMatrix proj_N_c1, proj_N_c2;
    proj_N_c1.clear();
    proj_N_c2.clear();

    FloatMatrix proj_DN_c1, proj_DN_c2;
    proj_DN_c1.clear();
    proj_DN_c2.clear();

    double area1 = 0., area2 = 0.;

    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {

        this->computeNmatrixAt(ip, N);

        if ( this->nlGeometry == 0 ) {
            this->giveStiffnessMatrix_Eng(D, rMode, ip, tStep);
        } else if ( this->nlGeometry == 1 ) {
            this->giveStiffnessMatrix_dTdj(D, rMode, ip, tStep);
        } else {
            OOFEM_ERROR("nlgeometry must be 0 or 1!")
        }

        this->computeTransformationMatrixAt(ip, rotationMatGtoL);
        D.rotatedWith(rotationMatGtoL, 't');                      // transform stiffness to global coord system

        this->computeNmatrixAt(ip, N);
        DN.beProductOf(D, N);


        double dA = this->computeAreaAround(ip);

        if(ip->giveNaturalCoordinate(1) < xi_0) {
            area1 += dA;

            proj_N_c1.add(dA, N);
            proj_DN_c1.add(dA, DN);
        }
        else {
            area2 += dA;

            proj_N_c2.add(dA, N);
            proj_DN_c2.add(dA, DN);
        }
    }

//    printf("area: %e\n", area);
    proj_N_c1.times(1./area1);
    proj_DN_c1.times(1./area1);

    proj_N_c2.times(1./area2);
    proj_DN_c2.times(1./area2);


    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {

        if ( this->nlGeometry == 0 ) {
            this->giveStiffnessMatrix_Eng(D, rMode, ip, tStep);
        } else if ( this->nlGeometry == 1 ) {
            this->giveStiffnessMatrix_dTdj(D, rMode, ip, tStep);
        } else {
            OOFEM_ERROR("nlgeometry must be 0 or 1!")
        }

        this->computeTransformationMatrixAt(ip, rotationMatGtoL);
        D.rotatedWith(rotationMatGtoL, 't');                      // transform stiffness to global coord system

        if(ip->giveNaturalCoordinate(1) < xi_0) {

			DN.beProductOf(D, proj_N_c1);


			double dA = this->computeAreaAround(ip);
			if ( matStiffSymmFlag ) {
				answer.plusProductSymmUpper(proj_N_c1, DN, dA);
			} else {
				answer.plusProductUnsym(proj_N_c1, DN, dA);
			}
        }
        else {
			DN.beProductOf(D, proj_N_c2);


			double dA = this->computeAreaAround(ip);
			if ( matStiffSymmFlag ) {
				answer.plusProductSymmUpper(proj_N_c2, DN, dA);
			} else {
				answer.plusProductUnsym(proj_N_c2, DN, dA);
			}
        }
    }


    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
}


void
IntElLine1IntPen :: giveInternalForcesVector(FloatArray &answer,
                                                       TimeStep *tStep, int useUpdatedGpRecord)
{
    // Computes internal forces
	// For this element we use an "interior penalty" formulation, where
	// the cohesive zone contribution is weakened, i.e. the traction and
	// test function for the cohesive zone are projected onto a reduced
	// space. The purpose of the projection is to improve the stability
	// properties of the formulation, thereby avoiding traction oscilations.

	double xi_0 = 0.0;

    FloatMatrix N;
    FloatArray u, traction, jump;

    this->computeVectorOf(VM_Total, tStep, u);
    // subtract initial displacements, if defined
    if ( initialDisplacements.giveSize() ) {
        u.subtract(initialDisplacements);
    }

    // zero answer will resize accordingly when adding first contribution
    answer.clear();



    // First loop over GP: compute projection of test function and traction.
    // The setting is as follows: we have an interface with quadratic interpolation and we
    // wish to project onto the space of constant functions on each element.


    FloatArray proj_jump_c1, proj_jump_c2;
    proj_jump_c1.clear();
    proj_jump_c2.clear();

    // Projecting the basis functions gives a constant for each basis function.
    FloatMatrix proj_N_c1, proj_N_c2;
    proj_N_c1.clear();
    proj_N_c2.clear();

    double area1 = 0., area2 = 0.;

    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {

        this->computeNmatrixAt(ip, N);
        jump.beProductOf(N, u);
//        this->computeTraction(traction, ip, jump, tStep);

        double dA = this->computeAreaAround(ip);

        if(ip->giveNaturalCoordinate(1) < xi_0) {
            area1 += dA;
            proj_jump_c1.add(dA, jump);
            proj_N_c1.add(dA, N);
        }
        else {
            area2 += dA;
            proj_jump_c2.add(dA, jump);
            proj_N_c2.add(dA, N);
        }

    }

//    printf("area1: %e\n", area1);
    proj_jump_c1.times(1./area1);
    proj_N_c1.times(1./area1);

//    printf("area2: %e\n", area2);
    proj_jump_c2.times(1./area2);
    proj_N_c2.times(1./area2);

//    printf("proj_N: "); proj_N.printYourself();


    // Second loop over GP: assemble contribution to internal forces as we are used to.
    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {
//        this->computeNmatrixAt(ip, N);
//        jump.beProductOf(N, u);
        if(ip->giveNaturalCoordinate(1) < xi_0) {

			this->computeTraction(traction, ip, proj_jump_c1, tStep);

			// compute internal cohesive forces as f = N^T*traction dA
			double dA = this->computeAreaAround(ip);
			answer.plusProduct(proj_N_c1, traction, dA);
        }
        else {

			this->computeTraction(traction, ip, proj_jump_c2, tStep);

			// compute internal cohesive forces as f = N^T*traction dA
			double dA = this->computeAreaAround(ip);
			answer.plusProduct(proj_N_c2, traction, dA);
        }
    }

}
#endif

void
IntElLine1IntPen :: computeNmatrixAt(GaussPoint *ip, FloatMatrix &answer)
{
    // Returns the modified N-matrix which multiplied with u give the spatial jump.

    FloatArray N;
    FEInterpolation *interp = this->giveInterpolation();
//    interp->evalN( N, ip->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    double xi_0 = 0.;
    if(ip->giveNaturalCoordinate(1) < xi_0 ){
    	N = {1., 0.};
    }
    else {
    	N = {0., 1.};
    }

    answer.resize(2, 8);
    answer.zero();
    answer.at(1, 1) = answer.at(2, 2) = -N.at(1);
    answer.at(1, 3) = answer.at(2, 4) = -N.at(2);

    answer.at(1, 5) = answer.at(2, 6) = N.at(1);
    answer.at(1, 7) = answer.at(2, 8) = N.at(2);
}

void
IntElLine1IntPen :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );

//        integrationRulesArray[ 0 ].reset( new LobattoIntegrationRule (1,this, 1, 2, false) );
//        integrationRulesArray [ 0 ]->SetUpPointsOnLine(2, _2dInterface);

        int numGP = 64;
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 2) );
        integrationRulesArray [ 0 ]->SetUpPointsOnLine(numGP, _2dInterface);
    }
}


} /* namespace oofem */
