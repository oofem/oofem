#include "tet21ghostsolid.h"
#include "fei3dtetquad.h"
#include "fei3dtetlin.h"
#include "nlstructuralelement.h"
#include "classfactory.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "structuralcrosssection.h"
#ifdef __FM_MODULE
 #include "fluiddynamicmaterial.h"
 #include "fluidcrosssection.h"
#endif
#include "load.h"
#include "element.h"
#include "dofmanager.h"

// Fredag?

namespace oofem {

#define USEUNCOUPLED 0
#define TESTNUMTAN 1

REGISTER_Element(tet21ghostsolid);

FEI3dTetQuad tet21ghostsolid :: interpolation;
FEI3dTetLin tet21ghostsolid :: interpolation_lin;

IntArray tet21ghostsolid :: momentum_ordering(30);
IntArray tet21ghostsolid :: conservation_ordering(4);
IntArray tet21ghostsolid :: ghostdisplacement_ordering(30);

tet21ghostsolid::tet21ghostsolid(int n, Domain *aDomain) : NLStructuralElement(n, aDomain)
{

    numberOfGaussPoints = 4;
    numberOfDofMans = 10;


    double nu=.25, E=1;
    Dghost.resize(6,6);
    Dghost.zero();
    Dghost.at(1,1) = 1-nu; Dghost.at(1,2) = Dghost.at(1,3) = nu;
    Dghost.at(2,2) = 1-nu; Dghost.at(2,1) = Dghost.at(2,3) = nu;
    Dghost.at(3,3) = 1-nu; Dghost.at(3,1) = Dghost.at(3,2) = nu;
    Dghost.at(4,4) = Dghost.at(5,5) = Dghost.at(6,6) = .5*(1-2*nu);
    Dghost.times(E/(1+nu)/(1-2*nu));

    conservation_ordering = { 7, 14, 21, 28};

    for ( int i = 0, j = 1; i < 10; ++i ) {
        momentum_ordering(i * 3 + 0) = j++;
        momentum_ordering(i * 3 + 1) = j++;
        momentum_ordering(i * 3 + 2) = j++;
        if ( i <= 3 ) { j++; }
        j=j+3;
    }

    for ( int i = 0, j = 1; i < 10; ++i ) {
        j=j+3;
        ghostdisplacement_ordering(i * 3 + 0) = j++;
        ghostdisplacement_ordering(i * 3 + 1) = j++;
        ghostdisplacement_ordering(i * 3 + 2) = j++;
        if ( i <= 3 ) { j++; }
    }

}

FEInterpolation *
tet21ghostsolid :: giveInterpolation() const
{
    return & interpolation;
}

void
tet21ghostsolid :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 6);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}

void
tet21ghostsolid :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
#ifdef __FM_MODULE
    FluidDynamicMaterial *fluidMaterial = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
#endif

    FloatMatrix Kf, G, Kx, D, B, Ed, EdB, dNx;
    FloatArray Nlin, dNv;

    for (int j = 0; j<iRule->giveNumberOfIntegrationPoints(); j++) {
        GaussPoint *gp = iRule->getIntegrationPoint(j);

        double detJ = fabs( ( this->interpolation.giveTransformationJacobian( * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) ) );
        double weight = gp->giveWeight();

        this->interpolation.evaldNdx( dNx, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
        this->interpolation_lin.evalN( Nlin, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

        dNv.resize(30); // dNv = [dN1/dx dN1/dy dN1/dz dN2/dx dN2/dy dN2/dz ... dN10/dz]

        for (int k = 0; k<dNx.giveNumberOfRows(); k++) {
            dNv.at(k*3+1) = dNx.at(k+1,1);
            dNv.at(k*3+2) = dNx.at(k+1,2);
            dNv.at(k*3+3) = dNx.at(k+1,3);
        }

        if (nlGeometry == 0) {

            this->computeBmatrixAt(gp, B);

            // Fluid part
            gp->setMaterialMode(_3dFlow);
#ifdef __FM_MODULE
            fluidMaterial->giveDeviatoricStiffnessMatrix(Ed, TangentStiffness, gp, tStep);
#else
            OOFEM_ERROR("Fluid module missing\n");
#endif
            gp->setMaterialMode(_3dMat);

            EdB.beProductOf(Ed, B);
            Kf.plusProductSymmUpper(B, EdB, detJ*weight);

            // Ghost solid part
            EdB.beProductOf(Dghost, B);
            Kx.plusProductSymmUpper(B, EdB, detJ*weight);

            // Incompressibility part
            G.plusDyadUnsym(dNv, Nlin, -detJ*weight);

        } else {

            this->computeBmatrixAt(gp, B);

            // Fluid part
            gp->setMaterialMode(_3dFlow);
            fluidMaterial->giveDeviatoricStiffnessMatrix(Ed, TangentStiffness, gp, tStep);
            gp->setMaterialMode(_3dMat);

            EdB.beProductOf(Ed, B);
            Kf.plusProductSymmUpper(B, EdB, detJ*weight);

            // Ghost solid part
            EdB.beProductOf(Dghost, B);
            Kx.plusProductSymmUpper(B, EdB, detJ*weight);

            // Incompressibility part
            G.plusDyadUnsym(dNv, Nlin, -detJ*weight);

        }

    }

    FloatMatrix GT;

    GT.beTranspositionOf(G);
    //GTdeltat.beTranspositionOf(G);
    //GTdeltat.times(deltat);
    Kf.symmetrized();
    Kx.symmetrized();
    //    Kf.printYourself();
    //    G.printYourself();
    //    GT.printYourself();
    //    Kx.printYourself();

    answer.resize(64, 64);
    answer.zero();
#define USEUNCOUPLED 0

#if USEUNCOUPLED == 1
    // Totaly uncoupled
    answer.assemble(Kf, momentum_ordering, momentum_ordering);
    answer.assemble(G, momentum_ordering, conservation_ordering);
    answer.assemble(GT, conservation_ordering, momentum_ordering);
    answer.assemble(Kx, ghostdisplacement_ordering, ghostdisplacement_ordering);
#else

    answer.assemble(Kf, ghostdisplacement_ordering, ghostdisplacement_ordering);
    answer.assemble(GT, conservation_ordering, ghostdisplacement_ordering);
    answer.assemble(Kf, ghostdisplacement_ordering, momentum_ordering);
    answer.assemble(G,  ghostdisplacement_ordering, conservation_ordering);
    answer.assemble(Kx, momentum_ordering, ghostdisplacement_ordering);
    answer.assemble(GT, conservation_ordering, momentum_ordering);

#endif

    if (this->giveNumber() == 364 ) {
        //answer.printYourself();
    }

}

void
tet21ghostsolid :: computeForceLoadVectorX(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
{

}

void
tet21ghostsolid :: computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep)
{

    answer.resize(64);
    answer.zero();

    if ( type != ExternalForcesVector ) {
        answer.resize(0);
        return;
    }

    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    FloatArray N, gVector, temparray(30), dNv, inc, a, a_prev, u, u_prev, vload;
    FloatMatrix dNx, G;

    load->computeComponentArrayAt(gVector, tStep, VM_Total);
    temparray.zero();

    vload.resize(4);
    vload.zero();

    this->giveUnknownData(a_prev, a, inc, tStep);
    u.beSubArrayOf(a, ghostdisplacement_ordering);
    u_prev.beSubArrayOf(a_prev, ghostdisplacement_ordering);

    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(k);
        FloatArray *lcoords = gp->giveNaturalCoordinates();

        double detJ = fabs( this->interpolation.giveTransformationJacobian( * lcoords, FEIElementGeometryWrapper(this) ) );
        double dA = detJ * gp->giveWeight();

        // Body load
        if ( gVector.giveSize() ) {
#ifdef __FM_MODULE
            double rho = mat->give('d', gp);
#else
            OOFEM_ERROR("Missing FM module");
            double rho = 1.0;
#endif
            this->interpolation.evalN( N, * lcoords, FEIElementGeometryWrapper(this) );

            for ( int j = 0; j < N.giveSize(); j++ ) {
                temparray(3 * j + 0) += N(j) * rho * gVector(0) * dA;
                temparray(3 * j + 1) += N(j) * rho * gVector(1) * dA;
                temparray(3 * j + 2) += N(j) * rho * gVector(2) * dA;
            }
        }

        // "load" from previous step
        this->interpolation.evaldNdx( dNx, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
        this->interpolation_lin.evalN( N, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

        dNv.resize(30);
        for (int k = 0; k<dNx.giveNumberOfRows(); k++) {
            dNv.at(k*3+1) = dNx.at(k+1,1);
            dNv.at(k*3+2) = dNx.at(k+1,2);
            dNv.at(k*3+3) = dNx.at(k+1,3);
        }

        G.plusDyadUnsym(N, dNv, dA);
        FloatMatrix GT;
        GT.beTranspositionOf(G);

        vload.plusProduct(GT, u_prev, -0.0 );
        //vload.printYourself();

    }

#if USEUNCOUPLED == 1
    // Totaly uncoupled
    answer.assemble(temparray, this->momentum_ordering);
#else
    answer.assemble(temparray, this->ghostdisplacement_ordering);
    answer.assemble(vload, this->conservation_ordering);
#endif



    // answer.printYourself();
}

void
tet21ghostsolid :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{

    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
#ifdef __FM_MODULE
    FluidDynamicMaterial *fluidMaterial = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
#endif

    FloatMatrix Kf, G, Kx, B, Ed, dNx;
    FloatArray Strain, Stress, Nlin, dNv, a, a_inc, a_prev, aVelocity, aPressure, aGhostDisplacement, aIncGhostDisplacement, fluidStress, epsf, dpdivv;
    FloatArray momentum, conservation, auxstress, divv;
    double pressure, epsvol=0.0;

    giveUnknownData(a_prev, a, a_inc, tStep);
    a_inc.operator =(a-a_prev);

#if TESTNUMTAN == 1
    FloatMatrix K;
    FloatArray f_int;
    double eps = 1e-3;
    K.resize(a.giveSize(), a.giveSize());
    K.zero();

    for (int i = 0; i<=a.giveSize(); i++) {
        giveUnknownData(a_prev, a, a_inc, tStep);
        if (i>0) { // Compute f_int(a) if i==0, otherwide compute f_int(a+delta)
            a.at(i) = a.at(i) + eps;
        }
        // Update a_inc
        a_inc.operator =(a-a_prev);
#endif

        momentum.resize(30);
        momentum.zero();
        conservation.resize(4);
        conservation.zero();
        auxstress.resize(30);
        auxstress.zero();

        aVelocity.beSubArrayOf(a, momentum_ordering);
        aPressure.beSubArrayOf(a, conservation_ordering);
        aGhostDisplacement.beSubArrayOf(a, ghostdisplacement_ordering);
        aIncGhostDisplacement.beSubArrayOf(a_inc, ghostdisplacement_ordering);

        for (int j = 0; j<iRule->giveNumberOfIntegrationPoints(); j++) {
            GaussPoint *gp = iRule->getIntegrationPoint(j);

            double detJ = fabs( ( this->interpolation.giveTransformationJacobian( * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) ) );
            double weight = gp->giveWeight();

            this->interpolation.evaldNdx( dNx, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
            this->interpolation_lin.evalN( Nlin, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

            dNv.resize(30);
            for (int k = 0; k<dNx.giveNumberOfRows(); k++) {
                dNv.at(k*3+1) = dNx.at(k+1,1);
                dNv.at(k*3+2) = dNx.at(k+1,2);
                dNv.at(k*3+3) = dNx.at(k+1,3);
            }

            if (nlGeometry == 0) {

                this->computeBmatrixAt(gp, B);
                FloatArray aTotal;
                aTotal.operator = (aVelocity + aIncGhostDisplacement); // Assume deltaT=1 gives that the increment is the velocity

                epsf.beProductOf(B, aTotal);
                pressure = Nlin.dotProduct(aPressure);

                // Momentum equation
                gp->setMaterialMode(_3dFlow);
                fluidMaterial->computeDeviatoricStressVector(fluidStress, epsvol, gp, epsf, pressure, tStep);
                gp->setMaterialMode(_3dMat);

                momentum.plusProduct(B, fluidStress, detJ*weight);
                momentum.add(-pressure * detJ * weight, dNv);

                // Conservation equation
                divv.beTProductOf(dNv, aTotal);
                divv.times(-detJ * weight);

                conservation.add(divv.at(1), Nlin);

                // Ghost solid part
                Strain.beProductOf(B, aGhostDisplacement);
                Stress.beProductOf(Dghost, Strain);
                auxstress.plusProduct(B, Stress, detJ * weight);

            } else {
                // We use small deformation on solid part since it is just there for keeping the mesh in ok shape.
                FloatMatrix B, BH;

                this->computeBHmatrixAt(gp, BH);
                this->computeBmatrixAt(gp, B);

                FloatArray aTotal;
                aTotal.operator = (aVelocity + aIncGhostDisplacement); // Assume deltaT=1 gives that the increment is the velocity

                epsf.beProductOf(B, aTotal);
                pressure = Nlin.dotProduct(aPressure);

                // Momentum equation -----

                // Compute fluid cauchy stress
                FloatArray fluidCauchy;
                FloatMatrix fluidCauchyMatrix;

                gp->setMaterialMode(_3dFlow);
                fluidMaterial->computeDeviatoricStressVector(fluidCauchy, epsvol, gp, epsf, pressure, tStep);
                gp->setMaterialMode(_3dMat);

                // Transform to 1st Piola-Kirshhoff
                FloatArray Fa;
                FloatMatrix F, Finv, FinvT, fluidStressMatrix;

                this->computeDeformationGradientVector(Fa, gp, tStep);
                F.beMatrixForm(Fa);
                Finv.beInverseOf(F);
                FinvT.beTranspositionOf(Finv);

                fluidCauchyMatrix.beMatrixFormOfStress(fluidCauchy);
                fluidStressMatrix.beProductOf(FinvT, fluidCauchyMatrix);
                fluidStressMatrix.times( detJ );
                fluidStress.beVectorForm(fluidStressMatrix);

                // Add to equation

                momentum.plusProduct(BH, fluidStress, detJ*weight);
                momentum.add(-pressure * detJ * weight, dNv);

                // Conservation equation -----
                divv.beTProductOf(dNv, aTotal);
                divv.times(-detJ * weight);

                conservation.add(divv.at(1), Nlin);

                // Ghost solid part -----
                Strain.beProductOf(B, aGhostDisplacement);
                Stress.beProductOf(Dghost, Strain);
                auxstress.plusProduct(B, Stress, detJ * weight);
            }

        }
        answer.resize(64);
        answer.zero();

#if USEUNCOUPLED == 1
        // Totaly uncoupled
        answer.assemble(momentum, momentum_ordering);
        answer.assemble(conservation, conservation_ordering);
        answer.assemble(auxstress, ghostdisplacement_ordering);
#else
        answer.assemble(momentum, ghostdisplacement_ordering);
        answer.assemble(conservation, conservation_ordering);
        answer.assemble(auxstress, momentum_ordering);
#endif

#if TESTNUMTAN == 1

        if (i==0) {
            f_int.operator =( answer );
        } else {
            // Compute derivative
            FloatArray dfda;

            //f_int.printYourself();
            //answer.printYourself();

            dfda.operator = (f_int - answer);
            dfda.times(-1.0/eps);

            //dfda.printYourself();

            K.setColumn(dfda, i);

        }
    }

    if (this->giveNumber() == 364) {
        //K.printYourself();
        //answer.printYourself();
    }
    answer.resize(64);
    answer.zero();

#if USEUNCOUPLED == 1
    // Totaly uncoupled
    answer.assemble(momentum, momentum_ordering);
    answer.assemble(conservation, conservation_ordering);
    answer.assemble(auxstress, ghostdisplacement_ordering);
#else
    answer.assemble(momentum, ghostdisplacement_ordering);
    answer.assemble(conservation, conservation_ordering);
    answer.assemble(auxstress, momentum_ordering);
#endif

//if (this->giveNumber() == 364) {
    // answer.printYourself();
//}

#if TESTNUMTAN == 1
    if (this->giveNumber() == 364) {
        FloatMatrix Ka, Diff;
        FloatArray ans;
        this->computeStiffnessMatrix(Ka, TangentStiffness, tStep);
        printf("Numerical tangent:\n");
        K.printYourself();
        printf("Analytical tangent:\n");
        Ka.printYourself();
        Diff = Ka;
        Diff.subtract(K);

        // Diff.printYourself();

        double m=0.0;
        for (int i=0; i<Diff.giveNumberOfRows(); i++) {
            for (int j=0; j<Diff.giveNumberOfColumns(); j++) {
                if (sqrt( Diff(i,j)*Diff(i,j) ) > m) m = sqrt( Diff(i,j)*Diff(i,j) );
            }
        }

        printf("%e\n", m);

        /*Ka.printYourself();
        ans.beProductOf(Ka, a);
        ans.printYourself();*/
    }
#endif

}

void
tet21ghostsolid :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    // Returns the mask for node number inode of this element.

    if ( inode <= 4 ) {
        answer = { V_u, V_v, V_w, D_u, D_v, D_w, P_f };
    } else {
        answer = { V_u, V_v, V_w, D_u, D_v, D_w };
    }
}

void
tet21ghostsolid :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [6x30] strain-displacement matrix {B} of the receiver, eva-
// luated at gp.
// B matrix  -  6 rows : epsilon-X, epsilon-Y, epsilon-Z, gamma-YZ, gamma-ZX, gamma-XY  :
{
    FloatMatrix dnx;

    this->interpolation.evaldNdx( dnx, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(6, 30);
    answer.zero();

    for ( int i = 1; i <= 10; i++ ) {
        answer.at(1, 3 * i - 2) = dnx.at(i, 1);
        answer.at(2, 3 * i - 1) = dnx.at(i, 2);
        answer.at(3, 3 * i - 0) = dnx.at(i, 3);

        answer.at(4, 3 * i - 1) = dnx.at(i, 3);
        answer.at(4, 3 * i - 0) = dnx.at(i, 2);

        answer.at(5, 3 * i - 2) = dnx.at(i, 3);
        answer.at(5, 3 * i - 0) = dnx.at(i, 1);

        answer.at(6, 3 * i - 2) = dnx.at(i, 2);
        answer.at(6, 3 * i - 1) = dnx.at(i, 1);
    }
    //answer.printYourself();
}

void
tet21ghostsolid :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    FloatMatrix dnx;

    this->interpolation.evaldNdx( dnx, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(9, 30);
    answer.zero();

    for ( int i = 1; i <= dnx.giveNumberOfRows(); i++ ) {
        answer.at(1, 3 * i - 2) = dnx.at(i, 1);     // du/dx
        answer.at(2, 3 * i - 1) = dnx.at(i, 2);     // dv/dy
        answer.at(3, 3 * i - 0) = dnx.at(i, 3);     // dw/dz
        answer.at(4, 3 * i - 1) = dnx.at(i, 3);     // dv/dz
        answer.at(7, 3 * i - 0) = dnx.at(i, 2);     // dw/dy
        answer.at(5, 3 * i - 2) = dnx.at(i, 3);     // du/dz
        answer.at(8, 3 * i - 0) = dnx.at(i, 1);     // dw/dx
        answer.at(6, 3 * i - 2) = dnx.at(i, 2);     // du/dy
        answer.at(9, 3 * i - 1) = dnx.at(i, 1);     // dv/dx
    }
}

void
tet21ghostsolid :: giveUnknownData(FloatArray &u_prev, FloatArray &u, FloatArray &inc, TimeStep *tStep)
{

    this->computeVectorOf( VM_Total, tStep, u);

    if (!tStep->isTheFirstStep()) {
        this->computeVectorOf( VM_Total, tStep->givePreviousStep(), u_prev);
        this->computeVectorOf( VM_Incremental, tStep, inc);
    } else {
        inc.resize( u.giveSize() );
        inc.zero();
        u_prev.operator = (inc);
    }

        for (int i = 0; i<this->giveNumberOfDofManagers(); i++) {

            this->giveDofManager(i+1)->giveUnknownVector(u_n, id, VM_Total, tStep );
            this->giveDofManager(i+1)->giveUnknownVector(inc_n, id, VM_Incremental, tStep );
            u_prev.at(3*i+1) = u_n.at(1)-inc_n.at(1);
            u_prev.at(3*i+2) = u_n.at(2)-inc_n.at(2);
            u_prev.at(3*i+3) = u_n.at(3)-inc_n.at(3);
            u.at(3*i+1) = u_n.at(1);
            u.at(3*i+2) = u_n.at(2);
            u.at(3*i+3) = u_n.at(3);
            inc.at(3*i+1) = inc_n.at(1);
            inc.at(3*i+2) = inc_n.at(2);
            inc.at(3*i+3) = inc_n.at(3);
        }
    }
}

}
