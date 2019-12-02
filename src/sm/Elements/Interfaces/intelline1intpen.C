/*
 * intelline1intpen.C
 *
 *  Created on: Nov 20, 2015
 *      Author: svennine
 */

#include "intelline1intpen.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "dofmanager.h"
#include "fei2dlinelin.h"
#include "fei2dlinequad.h"
#include "feinterpol.h"
#include "sm/Materials/InterfaceMaterials/structuralinterfacematerialstatus.h"

namespace oofem {
REGISTER_Element(IntElLine1IntPen);

IntElLine1IntPen::IntElLine1IntPen(int n, Domain * d) : IntElLine1(n, d)
{
    numberOfDofMans = 6;
    numberOfGaussPoints = 4;
}

int IntElLine1IntPen :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    FloatArray N;
    FEInterpolation *interp = this->giveInterpolation();
    interp->evalN( N, lcoords, FEIElementGeometryWrapper(this) );

    answer.resize(this->giveDofManager(1)->giveCoordinates().giveSize());
    answer.zero();

    double xi_0 = 0.;
    FloatArray xiScaled = {0.};

    if ( lcoords.at(1) < xi_0 ) {
        xiScaled = {lcoords.at(1)*2. + 1.};
        interp->evalN( N, xiScaled, FEIElementGeometryWrapper(this) );

        const auto &x1 = this->giveDofManager(1)->giveCoordinates();
        answer.add(N.at(1), x1 );

        const FloatArray &x3 = this->giveDofManager(3)->giveCoordinates();
        answer.add(N.at(2), x3 );
    } else {
        xiScaled = {lcoords.at(1)*2. - 1.};
        interp->evalN( N, xiScaled, FEIElementGeometryWrapper(this) );

        const auto &x3 = this->giveDofManager(3)->giveCoordinates();
        answer.add(N.at(1), x3 );

        const auto &x2 = this->giveDofManager(2)->giveCoordinates();
        answer.add(N.at(2), x2 );
    }

    return true;
}

FloatArrayF<2>
IntElLine1IntPen :: computeCovarBaseVectorAt(IntegrationPoint *ip) const
{
    //printf("Entering IntElLine2IntPen :: computeCovarBaseVectorAt\n");

    // Since we are averaging over the whole element, always evaluate the base vectors at xi = 0.

    FloatArray xi_0 = {0.0};
    //FloatArray xi_0 = {ip->giveNaturalCoordinate(1)};

    FloatMatrix dNdxi;
    FEInterpolation *interp = this->giveInterpolation();
    //interp->evaldNdxi( dNdxi, ip->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    interp->evaldNdxi( dNdxi, xi_0, FEIElementGeometryWrapper(this) );

    FloatArrayF<2> G;

    double X1_i = 0.5 * ( this->giveNode(1)->giveCoordinate(1) + this->giveNode(4)->giveCoordinate(1) ); // (mean) point on the fictious mid surface
    double X2_i = 0.5 * ( this->giveNode(1)->giveCoordinate(2) + this->giveNode(4)->giveCoordinate(2) );
    G.at(1) += dNdxi.at(1, 1) * X1_i;
    G.at(2) += dNdxi.at(1, 1) * X2_i;

    X1_i = 0.5 * ( this->giveNode(2)->giveCoordinate(1) + this->giveNode(5)->giveCoordinate(1) ); // (mean) point on the fictious mid surface
    X2_i = 0.5 * ( this->giveNode(2)->giveCoordinate(2) + this->giveNode(5)->giveCoordinate(2) );
    G.at(1) += dNdxi.at(2, 1) * X1_i;
    G.at(2) += dNdxi.at(2, 1) * X2_i;
    return G;
}


void
IntElLine1IntPen :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
#if 1
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
    FloatMatrix proj_N;
    proj_N.clear();

    FloatMatrix proj_DN;
    proj_DN.clear();

    double area = 0.;

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
        D.rotatedWith(rotationMatGtoL, 'n');                      // transform stiffness to global coord system

        this->computeNmatrixAt(ip, N);
        DN.beProductOf(D, N);


        double dA = this->computeAreaAround(ip);
        area += dA;

        proj_N.add(dA, N);
        proj_DN.add(dA, DN);
    }

//    printf("area: %e\n", area);
    proj_N.times(1./area);
    proj_DN.times(1./area);

//    printf("proj_N: "); proj_N.printYourself();

    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {

        if ( this->nlGeometry == 0 ) {
            this->giveStiffnessMatrix_Eng(D, rMode, ip, tStep);
        } else if ( this->nlGeometry == 1 ) {
            this->giveStiffnessMatrix_dTdj(D, rMode, ip, tStep);
        } else {
            OOFEM_ERROR("nlgeometry must be 0 or 1!")
        }

        this->computeTransformationMatrixAt(ip, rotationMatGtoL);
        D.rotatedWith(rotationMatGtoL, 'n');                      // transform stiffness to global coord system

        DN.beProductOf(D, proj_N);


        double dA = this->computeAreaAround(ip);
        if ( matStiffSymmFlag ) {
            answer.plusProductSymmUpper(proj_N, DN, dA);
        } else {
            answer.plusProductUnsym(proj_N, DN, dA);
        }
    }


    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }

#else

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
    FloatMatrix proj_N;
    proj_N.clear();

    FloatMatrix proj_DN;
    proj_DN.clear();

    double area = 0.;

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
        D.rotatedWith(rotationMatGtoL, 'n');                      // transform stiffness to global coord system

        this->computeNmatrixAt(ip, N);
        DN.beProductOf(D, N);


        double dA = this->computeAreaAround(ip);
        area += dA;

        proj_N.add(dA, N);
        proj_DN.add(dA, DN);
    }

//    printf("area: %e\n", area);
    proj_N.times(1./area);
    proj_DN.times(1./area);

//    printf("proj_N: "); proj_N.printYourself();

    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {

        if ( this->nlGeometry == 0 ) {
            this->giveStiffnessMatrix_Eng(D, rMode, ip, tStep);
        } else if ( this->nlGeometry == 1 ) {
            this->giveStiffnessMatrix_dTdj(D, rMode, ip, tStep);
        } else {
            OOFEM_ERROR("nlgeometry must be 0 or 1!")
        }

//        this->computeTransformationMatrixAt(ip, rotationMatGtoL);
//        D.rotatedWith(rotationMatGtoL, 'n');                      // transform stiffness to global coord system
//
//        DN.beProductOf(D, proj_N);


        double dA = this->computeAreaAround(ip);
        if ( matStiffSymmFlag ) {
            answer.plusProductSymmUpper(proj_N, proj_DN, dA);
//            answer.plusProductSymmUpper(proj_DN, proj_N, dA);
        } else {
            answer.plusProductUnsym(proj_N, proj_DN, dA);
//            answer.plusProductUnsym(proj_DN, proj_N, dA);
        }
    }

    if ( matStiffSymmFlag ) {
        answer.symmetrized();
    }
#endif
}


void
IntElLine1IntPen :: giveInternalForcesVector(FloatArray &answer,
                                                       TimeStep *tStep, int useUpdatedGpRecord)
{
#if 1
    // Computes internal forces
    // For this element we use an "interior penalty" formulation, where
    // the cohesive zone contribution is weakened, i.e. the traction and
    // test function for the cohesive zone are projected onto a reduced
    // space. The purpose of the projection is to improve the stability
    // properties of the formulation, thereby avoiding traction oscilations.

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

    // The projection of t becomes a constant
//    FloatArray proj_t;
//    proj_t.clear();


    FloatArray proj_jump;
    proj_jump.clear();

    // Projecting the basis functions gives a constant for each basis function.
    FloatMatrix proj_N;
    proj_N.clear();

    double area = 0.;

    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {

        this->computeNmatrixAt(ip, N);
        jump.beProductOf(N, u);
//        this->computeTraction(traction, ip, jump, tStep);

        double dA = this->computeAreaAround(ip);
        area += dA;

        proj_jump.add(dA, jump);
        proj_N.add(dA, N);
    }

//    printf("area: %e\n", area);
    proj_jump.times(1./area);
    proj_N.times(1./area);

//    printf("proj_N: "); proj_N.printYourself();


    /////////////////////////////////////////////////////////
    //
    //	Find c2 such that c1 = t{c2}
    //
    /////////////////////////////////////////////////////////

//    // Initial guess
//    FloatArray c2 = proj_jump;
//    int maxIter = 10;
//    for(int iter = 0; iter < maxIter; iter++) {
//
//        this->computeTraction(traction, ip, c2, tStep);
//
//
//    }


    // Second loop over GP: assemble contribution to internal forces as we are used to.
    for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {
//        this->computeNmatrixAt(ip, N);
//        jump.beProductOf(N, u);
        this->computeTraction(traction, ip, proj_jump, tStep);

        // compute internal cohesive forces as f = N^T*traction dA
        double dA = this->computeAreaAround(ip);
        answer.plusProduct(proj_N, traction, dA);
    }

}

#else
// Computes internal forces
// For this element we use an "interior penalty" formulation, where
// the cohesive zone contribution is weakened, i.e. the traction and
// test function for the cohesive zone are projected onto a reduced
// space. The purpose of the projection is to improve the stability
// properties of the formulation, thereby avoiding traction oscilations.

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

// The projection of t becomes a constant
    FloatArray proj_t;
    proj_t.clear();


//FloatArray proj_jump;
//proj_jump.clear();

// Projecting the basis functions gives a constant for each basis function.
FloatMatrix proj_N;
proj_N.clear();

double area = 0.;

for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {

    this->computeNmatrixAt(ip, N);
    jump.beProductOf(N, u);
    this->computeTraction(traction, ip, jump, tStep);

    double dA = this->computeAreaAround(ip);
    area += dA;

//    printf("traction: "); traction.printYourself();
//    proj_jump.add(dA, jump);
    proj_t.add(dA, traction);
    proj_N.add(dA, N);
}

//    printf("area: %e\n", area);
//proj_jump.times(1./area);
proj_t.times(1./area);
//printf("proj_t: "); proj_t.printYourself();
//printf("\n\n");
proj_N.times(1./area);

//    printf("proj_N: "); proj_N.printYourself();

FloatMatrix rotationMatGtoL;

// Second loop over GP: assemble contribution to internal forces as we are used to.
for ( auto &ip: *this->giveDefaultIntegrationRulePtr() ) {
        this->computeNmatrixAt(ip, N);
//        jump.beProductOf(N, u);
//    this->computeTraction(traction, ip, proj_jump, tStep);

    // compute internal cohesive forces as f = N^T*traction dA
    double dA = this->computeAreaAround(ip);
    answer.plusProduct(proj_N, proj_t, dA);
//    answer.plusProduct(N, proj_t, dA);

    StructuralInterfaceMaterialStatus *status = static_cast< StructuralInterfaceMaterialStatus * >( ip->giveMaterialStatus() );
    FloatArray proj_t_gp = proj_t;

    this->computeTransformationMatrixAt(ip, rotationMatGtoL);
    proj_t_gp.rotatedWith(rotationMatGtoL, 'n');                      // transform to local coord system

    FloatArray proj_t_gp_3D = {proj_t_gp.at(1), 0., proj_t_gp.at(2)};
    status->letProjectedTractionBe(proj_t_gp_3D);

//    FloatArray proj_t_gp = proj_t;
//    printf("proj_t_gp: "); proj_t_gp.printYourself();
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
    FloatArray xiScaled = {0.};

    answer.resize(2, 12);
    answer.zero();

    if ( ip->giveNaturalCoordinate(1) < xi_0 ) {
        xiScaled = {ip->giveNaturalCoordinate(1)*2. + 1.};
        interp->evalN( N, xiScaled, FEIElementGeometryWrapper(this) );

        answer.at(1, 1) = answer.at(2, 2) = -N.at(1);
        //answer.at(1, 3) = answer.at(2, 4) = -N.at(2);
        answer.at(1, 5) = answer.at(2, 6) = -N.at(2);

        answer.at(1, 7) = answer.at(2, 8) = N.at(1);
        //answer.at(1, 9) = answer.at(2, 10) = N.at(2);
        answer.at(1, 11) = answer.at(2, 12) = N.at(2);
    } else {
        xiScaled = {ip->giveNaturalCoordinate(1)*2. - 1.};
        interp->evalN( N, xiScaled, FEIElementGeometryWrapper(this) );

        //answer.at(1, 1) = answer.at(2, 2) = -N.at(1);
        answer.at(1, 3) = answer.at(2, 4) = -N.at(2);
        answer.at(1, 5) = answer.at(2, 6) = -N.at(1);

        //answer.at(1, 7) = answer.at(2, 8) = N.at(1);
        answer.at(1, 9) = answer.at(2, 10) = N.at(2);
        answer.at(1, 11) = answer.at(2, 12) = N.at(1);
    }

}

void
IntElLine1IntPen :: computeGaussPoints()
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );

//        integrationRulesArray[ 0 ] = std::make_unique<LobattoIntegrationRule>(1,this, 1, 2, false);
//        integrationRulesArray [ 0 ]->SetUpPointsOnLine(2, _2dInterface);

        int numGP = 4;
        integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this, 1, 2);
        integrationRulesArray [ 0 ]->SetUpPointsOnLine(numGP, _2dInterface);
    }
}


} /* namespace oofem */
