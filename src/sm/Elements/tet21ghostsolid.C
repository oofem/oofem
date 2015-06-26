/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2014   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "../sm/Elements/tet21ghostsolid.h"
#include "../sm/Elements/nlstructuralelement.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "fei3dtetquad.h"
#include "fei3dtetlin.h"
#include "classfactory.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "load.h"
#include "element.h"
#include "dofmanager.h"
#include "load.h"
#include "boundaryload.h"
#include "neumannmomentload.h"
#include "dof.h"

#ifdef __FM_MODULE
#include "fluiddynamicmaterial.h"
#include "fluidcrosssection.h"
#endif

namespace oofem {

#define USEUNCOUPLED 0
#define USENUMTAN 1

REGISTER_Element(tet21ghostsolid);

FEI3dTetQuad tet21ghostsolid :: interpolation;
FEI3dTetLin tet21ghostsolid :: interpolation_lin;

IntArray tet21ghostsolid :: momentum_ordering(30);
IntArray tet21ghostsolid :: conservation_ordering(4);
IntArray tet21ghostsolid :: ghostdisplacement_ordering(30);
IntArray tet21ghostsolid :: velocitydofsonside = {1, 2, 3, 8, 9, 10, 15, 16, 17, 22, 23, 24, 28, 29, 30, 34, 35, 36};
IntArray tet21ghostsolid :: displacementdofsonside = {4, 5, 6, 11, 12, 13, 18, 19, 20, 25, 26, 27, 31, 32, 33, 37, 38, 39};

tet21ghostsolid::tet21ghostsolid(int n, Domain *aDomain) : NLStructuralElement(n, aDomain), SpatialLocalizerInterface(this)
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

FEInterpolation *
tet21ghostsolid :: giveInterpolation(DofIDItem id) const
{
    if ( id == P_f ) {
        return & interpolation_lin;
    } else {
        return & interpolation;
    }
}

void
tet21ghostsolid :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 6) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}

void
tet21ghostsolid :: computeNumericStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

    FloatArray a, aPert, intF, intFPert, DintF;
    double eps = 1e-6;

    answer.resize(64, 64);
    answer.zero();

    this->computeVectorOf(VM_Total, tStep, a);

    giveInternalForcesVectorGivenSolution(intF, tStep, 0, a);

    for (int i = 1; i <= answer.giveNumberOfColumns(); i++) {
        aPert = a;
        aPert.at(i) += eps;

        giveInternalForcesVectorGivenSolution(intFPert, tStep, 0, aPert);

        DintF.beDifferenceOf(intF, intFPert);
        DintF.times(-1/eps);
        answer.setColumn(DintF, i);
    }


}

void
tet21ghostsolid :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveRealStress_3d(answer, gp, strain, tStep);
}

void
tet21ghostsolid :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

#ifdef USENUMTAN
    computeNumericStiffnessMatrix(answer, rMode, tStep);
/*    if (this->globalNumber == 292) {
        printf("Tangent: ");
        answer.printYourself();
    } */

    return;
#endif

#ifdef __FM_MODULE
    FluidDynamicMaterial *fluidMaterial = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
#endif

    FloatMatrix Kf, G, Kx, D, B, Ed, EdB, dNx;
    FloatArray Nlin, dNv;

    for ( auto &gp: *this->giveDefaultIntegrationRulePtr() ) {
        double detJ = fabs( ( this->interpolation.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) ) );
        double weight = gp->giveWeight();

        this->interpolation.evaldNdx( dNx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
        this->interpolation_lin.evalN( Nlin, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

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
    FloatArray N, gVector, temparray(30), dNv, inc, a, a_prev, u, u_prev, vload;
    FloatMatrix dNx, G;

    load->computeComponentArrayAt(gVector, tStep, VM_Total);
    temparray.zero();

    vload.resize(4);
    vload.zero();

    this->giveUnknownData(a_prev, a, inc, tStep);
    u.beSubArrayOf(a, ghostdisplacement_ordering);
    u_prev.beSubArrayOf(a_prev, ghostdisplacement_ordering);

    for ( auto &gp: *this->integrationRulesArray [ 0 ] ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();

        double detJ = fabs( this->interpolation.giveTransformationJacobian( lcoords, FEIElementGeometryWrapper(this) ) );
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
                computeDeformationGradientVectorFromDispl(Fa, gp, tStep, u);
                F.beMatrixForm(Fa);
                Finv.beInverseOf(F);
                FinvT.beTranspositionOf(Finv);

                dA = dA * F.giveDeterminant();

                temp.beProductOf(Finv, gVector);
                gVector = temp;
            }

            this->interpolation.evalN( N, lcoords, FEIElementGeometryWrapper(this) );

            for ( int j = 0; j < N.giveSize(); j++ ) {
                temparray(3 * j + 0) += N(j) * rho * gVector(0) * dA;
                temparray(3 * j + 1) += N(j) * rho * gVector(1) * dA;
                temparray(3 * j + 2) += N(j) * rho * gVector(2) * dA;
            }
        }

        // "load" from previous step
        this->interpolation.evaldNdx( dNx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
        this->interpolation_lin.evalN( N, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

        dNv.resize(30);
        for (int k = 0; k < dNx.giveNumberOfRows(); k++) {
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
        giveRowTransformationMatrix(tStep);
    }

    FloatArray temp;

    temp.resize(64);
    temp.zero();
    temp.assemble(temparray, this->ghostdisplacement_ordering);
    temp.assemble(vload, this->conservation_ordering);

    answer.beProductOf(Itransform, temp);

/*    if (this->globalNumber == 292) {
        printf("LoadVector: ");
        answer.printYourself();
    } */

}

void
tet21ghostsolid :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{

    FloatArray a;
    this->computeVectorOf(VM_Total, tStep, a);
    giveInternalForcesVectorGivenSolution(answer, tStep, useUpdatedGpRecord, a);

/*    if (this->globalNumber == 292) {
        printf("InternalForces: ");
        answer.printYourself();
        if ( a.at(4) < -1.337e-2 && a.at(4)> -1.34e-2) {
            //printf(" ");
        }
        printf("Current and previous solutions: ");
        this->computeVectorOf( VM_Total, tStep, a);
        a.printYourself();
        FloatArray a_prev;
        this->computeVectorOf( VM_Total, tStep->givePreviousStep(), a_prev);
        a_prev.printYourself();
    }*/

}

void
tet21ghostsolid :: giveInternalForcesVectorGivenSolution(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord, FloatArray &a)
{
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
    a_inc.beDifferenceOf(a, a_prev);

    aVelocity.beSubArrayOf(a, momentum_ordering);
    aPressure.beSubArrayOf(a, conservation_ordering);
    aGhostDisplacement.beSubArrayOf(a, ghostdisplacement_ordering);
    aIncGhostDisplacement.beSubArrayOf(a_inc, ghostdisplacement_ordering);

    for ( auto &gp: *this->giveDefaultIntegrationRulePtr() ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        double detJ = fabs( ( this->interpolation.giveTransformationJacobian( lcoords, FEIElementGeometryWrapper(this) ) ) );
        double weight = gp->giveWeight();

        this->interpolation.evaldNdx( dNx, lcoords, FEIElementGeometryWrapper(this) );
        this->interpolation_lin.evalN( Nlin, lcoords, FEIElementGeometryWrapper(this) );

        dNv.resize(30);
        for (int k = 0; k < dNx.giveNumberOfRows(); k++) {
            dNv.at(k*3+1) = dNx.at(k+1,1);
            dNv.at(k*3+2) = dNx.at(k+1,2);
            dNv.at(k*3+3) = dNx.at(k+1,3);
        }

        if (nlGeometry == 0) {

            this->computeBmatrixAt(gp, B);
            FloatArray aTotal = aVelocity + aIncGhostDisplacement*velocityCoeff; // Assume deltaT=1 gives that the increment is the velocity

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
            FloatMatrix BH;

            this->computeBHmatrixAt(gp, BH);
            this->computeBmatrixAt(gp, B);

            // Compute deformation gradient etc.
            FloatArray Fa, FinvTa;
            FloatMatrix F, Finv, FinvT, fluidStressMatrix;

            computeDeformationGradientVectorFromDispl(Fa, gp, tStep, aGhostDisplacement);
            F.beMatrixForm(Fa);
            double J = F.giveDeterminant();
            Finv.beInverseOf(F);
            FinvT.beTranspositionOf(Finv);
            FinvTa.beVectorForm(FinvT);

            FloatArray aTotal = aVelocity + aIncGhostDisplacement*velocityCoeff; // Assume deltaT=1 gives that the increment is the velocity

            epsf.beProductOf(B, aTotal);
            pressure = Nlin.dotProduct(aPressure);

            // Momentum equation --------------------------------------------------

            // Compute fluid cauchy stress
            FloatArray fluidCauchy;
            FloatMatrix fluidCauchyMatrix;

            gp->setMaterialMode(_3dFlow);
#ifdef __FM_MODULE
            fluidMaterial->computeDeviatoricStressVector(fluidCauchy, epsvol, gp, epsf, pressure, tStep);
#else
            OOFEM_ERROR("Missing FM module");
#endif
            gp->setMaterialMode(_3dMat);

            // Transform to 1st Piola-Kirshhoff
            fluidCauchyMatrix.beMatrixFormOfStress(fluidCauchy);
            fluidStressMatrix.beProductOf(fluidCauchyMatrix, FinvT);
            fluidStressMatrix.times( J );
            fluidStress.beVectorForm(fluidStressMatrix);

            // Add to equation
            momentum.plusProduct(BH, fluidStress, detJ*weight);

            // Pressure term in momentum equation
            FloatArray ptemp;
            ptemp.beTProductOf(BH, FinvTa);
            momentum.add(-pressure * detJ * weight *J, ptemp);

            // Conservation equation ---------------------------------------------
            divv.beTProductOf(ptemp, aTotal);  // Also include external forces
            divv.times(-detJ * weight* J);

            conservation.add(divv.at(1), Nlin);

            // Ghost solid part --------------------------------------------------
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
        giveRowTransformationMatrix(tStep);
    }

    answer.beProductOf(Itransform, temp);

    /*if (this->globalNumber == 292) {
        printf("InternalForces: ");
        answer.printYourself();
        printf("Temp: ");
        temp.printYourself();
    }*/

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

    this->interpolation.evaldNdx( dnx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

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
tet21ghostsolid :: computeDeformationGradientVectorAt(FloatArray &answer, FloatArray lcoord, TimeStep *tStep)
{

    FloatArray F, u;
    FloatMatrix dNdx, BH, Fmatrix, Finv;
    FEInterpolation *interpolation = this->giveInterpolation();

    // Fetch displacements
    this->computeVectorOf({1, 2, 3}, VM_Total, tStep, u);

    // Compute dNdx in point
    interpolation->evaldNdx(dNdx, lcoord, FEIElementGeometryWrapper (this) );

    // Compute displcement gradient BH
    BH.resize(9, dNdx.giveNumberOfRows() * 3);
    BH.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        BH.at(1, 3 * i - 2) = dNdx.at(i, 1);     // du/dx
        BH.at(2, 3 * i - 1) = dNdx.at(i, 2);     // dv/dy
        BH.at(3, 3 * i - 0) = dNdx.at(i, 3);     // dw/dz
        BH.at(4, 3 * i - 1) = dNdx.at(i, 3);     // dv/dz
        BH.at(7, 3 * i - 0) = dNdx.at(i, 2);     // dw/dy
        BH.at(5, 3 * i - 2) = dNdx.at(i, 3);     // du/dz
        BH.at(8, 3 * i - 0) = dNdx.at(i, 1);     // dw/dx
        BH.at(6, 3 * i - 2) = dNdx.at(i, 2);     // du/dy
        BH.at(9, 3 * i - 1) = dNdx.at(i, 1);     // dv/dx
    }

    // Finally, compute deformation gradient F=BH*u+I
    F.beProductOf(BH, u);
    F.at(1) += 1.0;
    F.at(2) += 1.0;
    F.at(3) += 1.0;

    answer = F;
}

void
tet21ghostsolid :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    FloatMatrix dnx;

    this->interpolation.evaldNdx( dnx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

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
        u_prev = inc;
    }

}

void
tet21ghostsolid :: computeDeformationGradientVectorFromDispl(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, FloatArray &u)
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
    // Deformation gradient F = H + I
    answer.beProductOf(B, u);
    answer.at(1) += 1.0;
    answer.at(2) += 1.0;
    answer.at(3) += 1.0;
}

void
tet21ghostsolid :: computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray u;
    this->computeVectorOf({1, 2, 3}, VM_Total, tStep, u);
    computeDeformationGradientVectorFromDispl(answer, gp, tStep, u);
}

double
tet21ghostsolid :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double determinant = fabs( this->interpolation.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
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
        this->interpolation.evalN( N, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

        Nmat.beNMatrixOf(N, 3);
        answer.beProductOf(Nmat, a);
        return 1;
    } else if (type == IST_Pressure) {
        FloatArray N, a;

        this->computeVectorOf({P_f}, VM_Total, tStep, a);
        this->interpolation_lin.evalN( N, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

        answer.resize(1);
        answer.at(1) = N.dotProduct(a);
        return 1;

    } else {
        MaterialMode matmode = gp->giveMaterialMode();
        gp->setMaterialMode(_3dFlow);
        int r = StructuralElement :: giveIPValue(answer, gp, type, tStep);
        gp->setMaterialMode(matmode);
        return r;
    }


}

bool
tet21ghostsolid :: giveRowTransformationMatrix(TimeStep *tStep)
{

    // Create a transformation matrix that switch all rows/equations located in OmegaF but not on GammaInt, i.e where we do not have a no slip condition

    Itransform.resize(64, 64);
    Itransform.zero();

    int m = 1;

    int row1, row2;
    IntArray row1List, row2List;

    for (int i = 1; i <= this->giveNumberOfDofManagers(); i++) {

        IntArray DofIDs;
        int firstIndex = m;
        this->giveDofManDofIDMask(i, DofIDs);
        // DofIDs.printYourself();

        for (int j = 1; j <= DofIDs.giveSize(); j++) {

            row1 = m;
            row2 = m;

            if (DofIDs.at(j) == V_u || DofIDs.at(j) == V_v || DofIDs.at(j) == V_w) {

                bool doSwitch = !this->giveDofManager(i)->giveDofWithID(DofIDs.at(j))->hasBc(tStep);

                if ( doSwitch ) {     // Boundary condition not set, make switch
                    row1=m;

                    for (int n = 1; n <= DofIDs.giveSize(); n++) {

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
    Itransform.beUnitMatrix();

    // Create tranformation matrix by switching rows
    for (int i = 1; i <= row1List.giveSize(); i++) {
        Itransform.at(row1List.at(i), row1List.at(i)) = 0;
        Itransform.at(row2List.at(i), row2List.at(i)) = 0;

        Itransform.at(row1List.at(i), row2List.at(i)) = 1;
        Itransform.at(row2List.at(i), row1List.at(i)) = 1;
    }

    computeItransform = false;

    return 0;

}

// Some extension Interfaces to follow:

Interface *tet21ghostsolid :: giveInterface(InterfaceType it)
{
    switch ( it ) {
    case NodalAveragingRecoveryModelInterfaceType:
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);

    case SpatialLocalizerInterfaceType:
        return static_cast< SpatialLocalizerInterface * >(this);

    case EIPrimaryUnknownMapperInterfaceType:
        return static_cast< EIPrimaryUnknownMapperInterface * >(this);

    default:
        return StructuralElement :: giveInterface(it);
        //return FMElement :: giveInterface(it);
    }
}

void tet21ghostsolid :: EIPrimaryUnknownMI_computePrimaryUnknownVectorAtLocal(ValueModeType mode,
                                                                              TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer)
{
    FloatArray n, n_lin;
    this->interpolation.evalN( n, lcoords, FEIElementGeometryWrapper(this) );
    this->interpolation_lin.evalN( n_lin, lcoords, FEIElementGeometryWrapper(this) );
    answer.resize(4);
    answer.zero();
    for ( int i = 1; i <= n.giveSize(); i++ ) {
        answer(0) += n.at(i) * this->giveNode(i)->giveDofWithID(V_u)->giveUnknown(mode, tStep);
        answer(1) += n.at(i) * this->giveNode(i)->giveDofWithID(V_v)->giveUnknown(mode, tStep);
        answer(2) += n.at(i) * this->giveNode(i)->giveDofWithID(V_w)->giveUnknown(mode, tStep);
    }

    for ( int i = 1; i <= n_lin.giveSize(); i++ ) {
        answer(3) += n_lin.at(i) * this->giveNode(i)->giveDofWithID(P_f)->giveUnknown(mode, tStep);
    }
}

double tet21ghostsolid :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray center;
    FloatArray lcoords = {0.3333333, 0.3333333, 0.3333333};
    this->computeGlobalCoordinates(center, lcoords);
    return center.distance(coords);
}

void
tet21ghostsolid :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_Pressure ) {
        answer.resize(1);
        if ( node <= 4 ) {
            answer.at(1) = this->giveNode(node)->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep);
        } else {
            IntArray eNodes;
            this->interpolation.computeLocalEdgeMapping(eNodes, node - 4);
            answer.at(1) = 0.5 * (
                        this->giveNode( eNodes.at(1) )->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) +
                        this->giveNode( eNodes.at(2) )->giveDofWithID(P_f)->giveUnknown(VM_Total, tStep) );
        }
    } else {
        answer.clear();
    }
}

void
tet21ghostsolid :: computeBoundaryLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep)
{
    answer.clear();
    if ( type != ExternalForcesVector ) {
        return;
    }

    FEInterpolation *fei = this->giveInterpolation();
    if ( !fei ) {
        OOFEM_ERROR("No interpolator available");
    }

    FloatArray n_vec, f;
    FloatMatrix n, T;
    FloatArray force;
    int nsd = fei->giveNsd();

    std :: unique_ptr< IntegrationRule > iRule( fei->giveBoundaryIntegrationRule(load->giveApproxOrder(), boundary) );

    for ( GaussPoint *gp: *iRule ) {
        FloatArray lcoords = gp->giveNaturalCoordinates();
        if ( load->giveFormulationType() == Load :: FT_Entity ) {
            load->computeValueAt(force, tStep, lcoords, mode);
        } else {
            FloatArray gcoords, elcoords;
            this->interpolation.surfaceLocal2global( gcoords, boundary, lcoords, FEIElementGeometryWrapper(this) );
            this->interpolation.global2local(elcoords, gcoords, FEIElementGeometryWrapper(this));
            NeumannMomentLoad *thisLoad = dynamic_cast<NeumannMomentLoad*> (load) ;
            if (thisLoad != NULL ) {
                FloatArray temp;
                thisLoad->computeValueAtBoundary(temp, tStep, gcoords, VM_Total, this, boundary);

                FloatArray F;
                FloatMatrix Fm, Finv;
                this->computeDeformationGradientVectorAt(F, elcoords, tStep);
                Fm.beMatrixForm(F);
                double J = Fm.giveDeterminant();
                Finv.beInverseOf(Fm);
                force.beTProductOf(Finv, temp);
                force.times(J);

            } else {
                load->computeValueAt(force, tStep, gcoords, VM_Total);
            }
        }

        ///@todo Make sure this part is correct.
        // We always want the global values in the end, so we might as well compute them here directly:
        // transform force
        if ( load->giveCoordSystMode() == Load :: CST_Global ) {
            // then just keep it in global c.s
        } else {
            ///@todo Support this...
            // transform from local boundary to element local c.s
            /*if ( this->computeLoadLSToLRotationMatrix(T, boundary, gp) ) {
             *  force.rotatedWith(T, 'n');
             * }*/
            // then to global c.s
            if ( this->computeLoadGToLRotationMtrx(T) ) {
                force.rotatedWith(T, 't');
            }
        }

        // Construct n-matrix
        fei->boundaryEvalN( n_vec, boundary, lcoords, FEIElementGeometryWrapper(this) );
        n.beNMatrixOf(n_vec, nsd);

        ///@todo Some way to ask for the thickness at a global coordinate maybe?
        double thickness = 1.0; // Should be the circumference for axisymm-elements.
        double dV = thickness * gp->giveWeight() * fei->boundaryGiveTransformationJacobian( boundary, lcoords, FEIElementGeometryWrapper(this) );
        f.plusProduct(n, force, dV);
    }

    answer.resize(39);
    answer.zero();
    answer.assemble(f, this->velocitydofsonside);
}

}
