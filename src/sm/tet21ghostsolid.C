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

namespace oofem {

#define USEUNCOUPLED 0
#define USENUMTAN 1

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
    computeItransform = true;

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
tet21ghostsolid :: computeNumericStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

    FloatArray a, aPert, intF, intFPert, DintF;
    double eps=1e-6;

    answer.resize(64, 64);
    answer.zero();

    this->computeVectorOf(VM_Total, tStep, a);

    giveInternalForcesVectorGivenSolution(intF, tStep, 0, a);

    for (int i=1; i<=answer.giveNumberOfColumns(); i++) {
        aPert = a;
        aPert.at(i) = aPert.at(i) + eps;

        giveInternalForcesVectorGivenSolution(intFPert, tStep, 0, aPert);

        DintF.operator =(intF-intFPert);
        DintF.times(-1/eps);
        answer.setColumn(DintF, i);
    }


}

void
tet21ghostsolid :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

#ifdef USENUMTAN
    computeNumericStiffnessMatrix(answer, rMode, tStep);
    return;
#endif

#ifdef __FM_MODULE
    FluidDynamicMaterial *fluidMaterial = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
#endif

    FloatMatrix Kf, G, Kx, D, B, Ed, EdB, dNx;
    FloatArray Nlin, dNv;

    for ( auto &gp: *this->giveDefaultIntegrationRulePtr() ) {
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
#ifdef __FM_MODULE
            fluidMaterial->giveDeviatoricStiffnessMatrix(Ed, TangentStiffness, gp, tStep);
#else
            OOFEM_ERROR("Fluid module missing\n");
#endif

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
    Kf.symmetrized();
    Kx.symmetrized();

    answer.resize(64, 64);
    answer.zero();

    answer.assemble(Kf, ghostdisplacement_ordering, ghostdisplacement_ordering);
    answer.assemble(GT, conservation_ordering, ghostdisplacement_ordering);
    answer.assemble(Kf, ghostdisplacement_ordering, momentum_ordering);
    answer.assemble(G,  ghostdisplacement_ordering, conservation_ordering);
    answer.assemble(Kx, momentum_ordering, ghostdisplacement_ordering);
    answer.assemble(GT, conservation_ordering, momentum_ordering);

}

void
tet21ghostsolid :: computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep)
{

    // Compute displacements used to compute J
    /*FloatArray a, aGhostDisplacement;
    this->computeVectorOf(VM_Total, tStep, a);
    aGhostDisplacement.beSubArrayOf(a, ghostdisplacement_ordering); */

    answer.resize(64);
    answer.zero();

    if ( type != ExternalForcesVector ) {
        answer.resize(0);
        return;
    }

#ifdef __FM_MODULE
    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
#endif
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

            if (this->nlGeometry) {
                FloatArray Fa, temp;
                FloatMatrix F, Finv, FinvT;
                computeDeformationGradientVector(Fa, gp, tStep, u);
                F.beMatrixForm(Fa);
                Finv.beInverseOf(F);
                FinvT.beTranspositionOf(Finv);

                dA = dA * F.giveDeterminant();

                temp.beProductOf(Finv, gVector);
                gVector = temp;
            }

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

    if ( this->computeItransform ) {
        giveRowTransformationMatrix(Itransform, tStep);
    }

    FloatArray temp;

    temp.resize(64);
    temp.zero();
    temp.assemble(temparray, this->ghostdisplacement_ordering);
    temp.assemble(vload, this->conservation_ordering);

    answer.beProductOf(Itransform, temp);

}

void
tet21ghostsolid :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{

    FloatArray a;
    this->computeVectorOf(VM_Total, tStep, a);
    giveInternalForcesVectorGivenSolution(answer, tStep, useUpdatedGpRecord, a);

}

void
tet21ghostsolid :: giveInternalForcesVectorGivenSolution(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord, FloatArray &a)
{

    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
#ifdef __FM_MODULE
    FluidDynamicMaterial *fluidMaterial = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
#endif

    FloatMatrix Kf, G, Kx, B, Ed, dNx;
    FloatArray Strain, Stress, Nlin, dNv, a_inc, a_prev, aVelocity, aPressure, aGhostDisplacement, aIncGhostDisplacement, fluidStress, epsf, dpdivv;
    FloatArray momentum, conservation, auxstress, divv;
    double pressure, epsvol=0.0, velocityCoeff=1.0;

    if (!tStep->isTheFirstStep()) {
        this->computeVectorOf( VM_Total, tStep->givePreviousStep(), a_prev);
    } else {
        a_prev.resize( a.giveSize() );
        a_prev.zero();
    }
    a_inc.operator =(a-a_prev);

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
            aTotal.operator = (aVelocity + aIncGhostDisplacement*velocityCoeff); // Assume deltaT=1 gives that the increment is the velocity

            epsf.beProductOf(B, aTotal);
            pressure = Nlin.dotProduct(aPressure);

            // Momentum equation
            gp->setMaterialMode(_3dFlow);
#ifdef __FM_MODULE
            fluidMaterial->computeDeviatoricStressVector(fluidStress, epsvol, gp, epsf, pressure, tStep);
#else
            OOFEM_ERROR("Missing FM module");
#endif
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
            FloatMatrix B, BH;

            this->computeBHmatrixAt(gp, BH);
            this->computeBmatrixAt(gp, B);

            FloatArray aTotal;
            aTotal.operator = (aVelocity + aIncGhostDisplacement*velocityCoeff); // Assume deltaT=1 gives that the increment is the velocity

            epsf.beProductOf(B, aTotal);
            pressure = Nlin.dotProduct(aPressure);

            // Momentum equation -----

            // Compute fluid cauchy stress
            FloatArray fluidCauchy;
            FloatMatrix fluidCauchyMatrix;

            gp->setMaterialMode(_3dFlow);
            fluidMaterial->computeDeviatoricStressVector(fluidCauchy, epsvol, gp, epsf, pressure, tStep);
            //fluidCauchy.printYourself();
            gp->setMaterialMode(_3dMat);

            // Transform to 1st Piola-Kirshhoff
            FloatArray Fa;
            FloatMatrix F, Finv, FinvT, fluidStressMatrix;

            computeDeformationGradientVector(Fa, gp, tStep, aGhostDisplacement);
            F.beMatrixForm(Fa);
            double J=F.giveDeterminant();
            Finv.beInverseOf(F);
            FinvT.beTranspositionOf(Finv);

            fluidCauchyMatrix.beMatrixFormOfStress(fluidCauchy);
            fluidStressMatrix.beProductOf(FinvT, fluidCauchyMatrix);
            fluidStressMatrix.times( F.giveDeterminant() );
            fluidStress.beVectorForm(fluidStressMatrix);
            //fluidStress.printYourself();

            // Add to equation

            momentum.plusProduct(BH, fluidStress, detJ*weight);
            momentum.add(-pressure * detJ * weight *J, dNv);

            // Conservation equation -----
            divv.beTProductOf(dNv, aTotal);
            divv.times(-detJ * weight*J);

            conservation.add(divv.at(1), Nlin);

            // Ghost solid part -----
            Strain.beProductOf(B, aGhostDisplacement);
            Stress.beProductOf(Dghost, Strain);
            auxstress.plusProduct(B, Stress, detJ * weight);
        }

    }
    answer.resize(64);
    answer.zero();

    FloatArray temp;

    temp.resize(64);
    temp.zero();

    temp.assemble(momentum, ghostdisplacement_ordering);
    temp.assemble(conservation, conservation_ordering);
    temp.assemble(auxstress, momentum_ordering);

    if ( this->computeItransform ) {
        giveRowTransformationMatrix(Itransform, tStep);
    }

    answer.beProductOf(Itransform, temp);

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

}

void
tet21ghostsolid :: computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, FloatArray &u)
{
    // Computes the deformation gradient in the Voigt format at the Gauss point gp of
    // the receiver at time step tStep.
    // Order of components: 11, 22, 33, 23, 13, 12, 32, 31, 21 in the 3D.

    // Obtain the current displacement vector of the element and subtract initial displacements (if present)
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // Displacement gradient H = du/dX
    FloatMatrix B;
    this->computeBHmatrixAt(gp, B);
    answer.beProductOf(B, u);

    // Deformation gradient F = H + I
    MaterialMode matMode = gp->giveMaterialMode();
    if ( matMode == _3dMat || matMode == _PlaneStrain ) {
        answer.at(1) += 1.0;
        answer.at(2) += 1.0;
        answer.at(3) += 1.0;
    } else if ( matMode == _PlaneStress ) {
        answer.at(1) += 1.0;
        answer.at(2) += 1.0;
    } else if ( matMode == _1dMat ) {
        answer.at(1) += 1.0;
    } else {
        OOFEM_ERROR("MaterialMode is not supported yet (%s)", __MaterialModeToString(matMode) );
    }
}

double
tet21ghostsolid :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double determinant = fabs( this->interpolation.giveTransformationJacobian( * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    double weight      = gp->giveWeight();
    return ( this->computeVolume() );
    return ( determinant * weight );
}


int
tet21ghostsolid :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{

    if ( type == IST_Velocity ) {
        FloatArray N, a;
        FloatMatrix Nmat;

        this->computeVectorOf({V_u, V_v, V_w}, VM_Total, tStep, a );
        this->interpolation.evalN( N, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

        Nmat.resize(3, N.giveSize()*3);
        for (int i=1; i<= N.giveSize(); i ++ ) {
            Nmat.at(1,3*i-2) = N.at(i);
            Nmat.at(2,3*i-1) = N.at(i);
            Nmat.at(3,3*i-0) = N.at(i);
        }

        answer.beProductOf(Nmat, a);
        return 1;
    } else if (type == IST_Pressure) {
        FloatArray N, a;

        this->computeVectorOf({P_f}, VM_Total, tStep, a);
        this->interpolation_lin.evalN( N, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

        answer.resize(1);
        answer.at(1) = N.dotProduct(a);
        return 1;

    } else {
        MaterialMode matmode=gp->giveMaterialMode();
        gp->setMaterialMode(_3dFlow);
        int r = StructuralElement :: giveIPValue(answer, gp, type, tStep);
        gp->setMaterialMode(matmode);
        return r;
    }


}

bool
tet21ghostsolid :: giveRowTransformationMatrix(FloatMatrix &Itransform, TimeStep *tStep)
{

    // Create a transformation matrix that switch all rows/equations located in OmegaF but not on GammaInt, i.e where we do not have a no slip condition

    Itransform.resize(64, 64);
    Itransform.zero();

    int m=1;

    int row1, row2;
    IntArray row1List, row2List;

    for (int i=1; i<=this->giveNumberOfDofManagers(); i++) {

        IntArray DofIDs;
        int firstIndex = m;
        this->giveDofManDofIDMask(i, DofIDs);
        // DofIDs.printYourself();

        for (int j=1; j<=DofIDs.giveSize(); j++) {

            row1=m;
            row2=m;

            if (DofIDs.at(j) == V_u || DofIDs.at(j) == V_v || DofIDs.at(j) == V_w) {

                bool doSwitch=!this->giveDofManager(i)->giveDofWithID(DofIDs.at(j))->hasBc(tStep);

                if (doSwitch) {     // Boundary condition not set, make switch
                    row1=m;

                    for (int n=1; n<=DofIDs.giveSize(); n++) {

                        if ( (DofIDs.at(j) == V_u && DofIDs.at(n) == D_u) ||
                             (DofIDs.at(j) == V_v && DofIDs.at(n) == D_v) ||
                             (DofIDs.at(j) == V_w && DofIDs.at(n) == D_w) ) {
                            row2 = firstIndex + n - 1;
                        }

                    }
                    row1List.resizeWithValues(row1List.giveSize() + 1);
                    row1List.at(row1List.giveSize()) = row1;
                    row2List.resizeWithValues(row2List.giveSize() + 1);
                    row2List.at(row2List.giveSize()) = row2;

                }

            }
            m++;
        }

    }

    // Create identity matrix
    for (int i=1; i<=64; i++) {
        Itransform.at(i, i) = 1;
    }

    // Create tranformation matrix by switching rows
    for (int i=1; i<=row1List.giveSize(); i++) {
        Itransform.at(row1List.at(i), row1List.at(i)) = 0;
        Itransform.at(row2List.at(i), row2List.at(i)) = 0;

        Itransform.at(row1List.at(i), row2List.at(i)) = 1;
        Itransform.at(row2List.at(i), row1List.at(i)) = 1;
    }

    computeItransform = false;

    return 0;

}

}
