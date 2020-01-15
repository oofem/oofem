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

#include "sm/Elements/Shells/solidshell.h"
#include <Materials/structuralms.h>
#include "classfactory.h"
#include "fei3dhexalin.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "sm/CrossSections/structuralcrosssection.h"

namespace oofem {
REGISTER_Element(SolidShell);

FEI3dHexaLin SolidShell :: interpolation;

SolidShell :: SolidShell(int n, Domain *aDomain) : LSpace(n, aDomain)
{
    numberOfDofMans  = 8;
    numberOfGaussPoints = 8;
}


void
SolidShell :: postInitialize()
{
    LSpace :: postInitialize();
 

    int numEASparam = 0;

    switch ( this->EAS_type) { 
      case 1:
        numEASparam = 1;
        break;

      case 2:
        numEASparam = 3;
        break;
    }

    this->u_k.resize(numberOfDofMans*3);
    this->u_k.zero();
    this->fE.resize(numEASparam);
    this->fE.zero();
    this->KEC.resize(numEASparam, numberOfDofMans*3);
    this->KEC.zero();
    this->invKEE.resize(numEASparam, numEASparam);
    this->invKEE.zero();
    this->alpha.resize(numEASparam);
    this->alpha.zero();
}

void
SolidShell :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ] = std::make_unique<GaussIntegrationRule>(1, this, 1, 6);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}

FEInterpolation *SolidShell :: giveInterpolation() const { return & interpolation; }


void
SolidShell :: initializeFrom(InputRecord &ir)
{
    numberOfGaussPoints = 8;
    NLStructuralElement :: initializeFrom(ir);

    // Check if EAS should be used
    this->EAS_type = 0;
    if ( ir.hasField(_IFT_SolidShell_EAS_type) ) {
        IR_GIVE_FIELD(ir, this->EAS_type,   _IFT_SolidShell_EAS_type);

    }
}



void
SolidShell :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [ 6 x (8*3) ] strain-displacement matrix {B} of the receiver, eva-
// luated at gp.
// B matrix  -  6 rows : epsilon-X, epsilon-Y, epsilon-Z, gamma-YZ, gamma-ZX, gamma-XY  :
{
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdx; 
    const FloatArray &lCoords = gp->giveNaturalCoordinates();
    interp->evaldNdx( dNdx, lCoords, FEIElementGeometryWrapper(this) );

    answer.resize(6, dNdx.giveNumberOfRows() * 3);
    answer.zero();

    // Do ANS
    // Need to evaluate strain at four points
    FloatArray A = { 0.0, -1.0, 0.0};
    FloatArray B = { 1.0,  0.0, 0.0};
    FloatArray C = { 0.0,  1.0, 0.0};
    FloatArray D = {-1.0,  0.0, 0.0};

    FloatMatrix dNdxA, dNdxB, dNdxC, dNdxD;
    interp->evaldNdx( dNdxA, A, FEIElementGeometryWrapper(this) );
    interp->evaldNdx( dNdxB, B, FEIElementGeometryWrapper(this) );
    interp->evaldNdx( dNdxC, C, FEIElementGeometryWrapper(this) );
    interp->evaldNdx( dNdxD, D, FEIElementGeometryWrapper(this) );

    FloatMatrix dNdx0; 
    interp->evaldNdx( dNdx0, {lCoords.at(1), lCoords.at(2), 0.0}, FEIElementGeometryWrapper(this) );

    double NA = 0.5 * ( 1.0 - lCoords.at(2) );
    double NC = 0.5 * ( 1.0 + lCoords.at(2) );
    double NB = 0.5 * ( 1.0 - lCoords.at(1) );
    double ND = 0.5 * ( 1.0 + lCoords.at(1) );

    // E_xz = E(5) = N_A(gp) * B(A)(5,:) + N_C(gp) * B(C)(5,:)
    // E_yz = E(4) = N_B(gp) * B(B)(4,:) + N_D(gp) * B(D)(4,:)

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, 3 * i - 2) = dNdx.at(i, 1);
        answer.at(2, 3 * i - 1) = dNdx.at(i, 2);
        answer.at(6, 3 * i - 2) = dNdx.at(i, 2);
        answer.at(6, 3 * i - 1) = dNdx.at(i, 1);

        // Evaluate thickness strain at the mid-surface
        //answer.at(3, 3 * i - 0) = dNdx.at(i, 3);
        answer.at(3, 3 * i - 0) = dNdx0.at(i, 3);

        answer.at(4, 3 * i - 1) = NB * dNdxB.at(i, 3) + ND * dNdxD.at(i, 3);
        answer.at(4, 3 * i - 0) = NB * dNdxB.at(i, 2) + ND * dNdxD.at(i, 2);

        answer.at(5, 3 * i - 2) = NA * dNdxA.at(i, 3) + NC * dNdxC.at(i, 3);
        answer.at(5, 3 * i - 0) = NA * dNdxA.at(i, 1) + NC * dNdxC.at(i, 1);
    }
}


void
SolidShell :: computeBHmatrixAt(FloatArray &lCoords, FloatMatrix &answer)
// Returns the [ 6 x (8*3) ] strain-displacement matrix {B} of the receiver, eva-
// luated at gp.
// B matrix  -  6 rows : epsilon-X, epsilon-Y, epsilon-Z, gamma-YZ, gamma-ZX, gamma-XY  :
{
#if 1
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdx; 
    interp->evaldNdx( dNdx, lCoords, FEIElementGeometryWrapper(this) );

    answer.resize(9, dNdx.giveNumberOfRows() * 3);
    answer.zero();
    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, 3 * i - 2) = dNdx.at(i, 1);     // du/dx
        answer.at(2, 3 * i - 1) = dNdx.at(i, 2);     // dv/dy
        answer.at(3, 3 * i - 0) = dNdx.at(i, 3);     // dw/dz
        answer.at(4, 3 * i - 1) = dNdx.at(i, 3);     // dv/dz
        answer.at(7, 3 * i - 0) = dNdx.at(i, 2);     // dw/dy
        answer.at(5, 3 * i - 2) = dNdx.at(i, 3);     // du/dz
        answer.at(8, 3 * i - 0) = dNdx.at(i, 1);     // dw/dx
        answer.at(6, 3 * i - 2) = dNdx.at(i, 2);     // du/dy
        answer.at(9, 3 * i - 1) = dNdx.at(i, 1);     // dv/dx
    }

#else
    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdx; 
    interp->evaldNdx( dNdx, lCoords, FEIElementGeometryWrapper(this) );

    answer.resize(9, dNdx.giveNumberOfRows() * 3);
    answer.zero();

    // Need to evaluate strain at four points
    FloatArray A = { 0.0, -1.0, 0.0};
    FloatArray B = { 1.0,  0.0, 0.0};
    FloatArray C = { 0.0,  1.0, 0.0};
    FloatArray D = {-1.0,  0.0, 0.0};

    FloatMatrix dNdxA, dNdxB, dNdxC, dNdxD;
    interp->evaldNdx( dNdxA, A, FEIElementGeometryWrapper(this) );
    interp->evaldNdx( dNdxB, B, FEIElementGeometryWrapper(this) );
    interp->evaldNdx( dNdxC, C, FEIElementGeometryWrapper(this) );
    interp->evaldNdx( dNdxD, D, FEIElementGeometryWrapper(this) );

    FloatMatrix dNdx0; 
    interp->evaldNdx( dNdx0, {lCoords.at(1), lCoords.at(2), 0.0}, FEIElementGeometryWrapper(this) );

    double NA = 0.5 * ( 1.0 - lCoords.at(2) );
    double NC = 0.5 * ( 1.0 + lCoords.at(2) );
    double NB = 0.5 * ( 1.0 - lCoords.at(1) );
    double ND = 0.5 * ( 1.0 + lCoords.at(1) );

    // E_xz = E(5) = N_A(gp) * B(A)(5,:) + N_C(gp) * B(C)(5,:)
    // E_yz = E(4) = N_B(gp) * B(B)(4,:) + N_D(gp) * B(D)(4,:)

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        answer.at(1, 3 * i - 2) = dNdx.at(i, 1);
        answer.at(2, 3 * i - 1) = dNdx.at(i, 2);
        answer.at(6, 3 * i - 2) = dNdx.at(i, 2);     // du/dy
        answer.at(9, 3 * i - 1) = dNdx.at(i, 1);     // dv/dx
        
        // Evaluate thickness strain at the mid-surface
        //answer.at(3, 3 * i - 0) = dNdx.at(i, 3);
        answer.at(3, 3 * i - 0) = dNdx0.at(i, 3);

        answer.at(4, 3 * i - 1) = NB * dNdxB.at(i, 3) + ND * dNdxD.at(i, 3);
        answer.at(7, 3 * i - 0) = NB * dNdxB.at(i, 2) + ND * dNdxD.at(i, 2);

        answer.at(5, 3 * i - 2) = NA * dNdxA.at(i, 3) + NC * dNdxC.at(i, 3);
        answer.at(8, 3 * i - 0) = NA * dNdxA.at(i, 1) + NC * dNdxC.at(i, 1);
    }
#endif
}


void
SolidShell :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the [ 6 x (8*3) ] strain-displacement matrix {B} of the receiver, eva-
// luated at gp.
// B matrix  -  6 rows : epsilon-X, epsilon-Y, epsilon-Z, gamma-YZ, gamma-ZX, gamma-XY  :
{
    this->computeBHmatrixAt(gp->giveNaturalCoordinates(), answer);
}


void
SolidShell :: computeEASBmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    const FloatArray &lCoords = gp->giveNaturalCoordinates();
    double xi   = lCoords.at(1);
    double eta  = lCoords.at(2);
    double zeta = lCoords.at(3);
    FloatMatrix M;

     switch ( this->EAS_type) { 
       case 1: 
        // alt 1
        M.resize(6,1);
        M.zero();
        M.at(3,1) = zeta;
        break;
      case 2:
        // alt 2
        M.resize(6,3);
        M.zero();
        M.at(1,1) = xi;
        M.at(2,2) = eta;
        M.at(3,3) = zeta;
        break;
    }

    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix J0mat; 
    FloatArray center = {0.0, 0.0, 0.0};

    interp->giveJacobianMatrixAt(J0mat, center, FEIElementGeometryWrapper(this) );
//     double J0 = interp->giveTransformationJacobian(center, FEIElementGeometryWrapper(this) );
    double J0 = J0mat.giveDeterminant();
    double J  = interp->giveTransformationJacobian(lCoords, FEIElementGeometryWrapper(this) );

    FloatMatrix T(6,6);
    this->computeBondTransformationMatrix(T, J0mat);
    //T.printYourself("T");
    //T.beUnitMatrix();
    answer.clear();
    answer.plusProductUnsym(T, M, J/J0);
}

 
void
SolidShell :: computeBondTransformationMatrix(FloatMatrix &answer, FloatMatrix &base)
{
    // given a bases {g1, g2, g3} compute the Bond (Voigt form) transformation matrix.
    // TODO is this with OOFEM ordering of components?
    FloatArray x, y, z;
    x.beColumnOf(base,1);
    y.beColumnOf(base,2);
    z.beColumnOf(base,3);
    answer = {
          { x(0) * x(0), x(1) * x(1), x(2) * x(2), x(1) * x(2), x(0) * x(2), x(0) * x(1) },
          { y(0) * y(0), y(1) * y(1), y(2) * y(2), y(1) * y(2), y(0) * y(2), y(0) * y(1) },
          { z(0) * z(0), z(1) * z(1), z(2) * z(2), z(1) * z(2), z(0) * z(2), z(0) * z(1) },
          { 2 * y(0) * z(0), 2 * y(1) * z(1), 2 * y(2) * z(2), y(2) * z(1) + y(1) * z(2), y(2) * z(0) + y(0) * z(2), y(1) * z(0) + y(0) * z(1) },
          { 2 * x(0) * z(0), 2 * x(1) * z(1), 2 * x(2) * z(2), x(2) * z(1) + x(1) * z(2), x(2) * z(0) + x(0) * z(2), x(1) * z(0) + x(0) * z(1) },
          { 2 * x(0) * y(0), 2 * x(1) * y(1), 2 * x(2) * y(2), x(2) * y(1) + x(1) * y(2), x(2) * y(0) + x(0) * y(2), x(1) * y(0) + x(0) * y(1) }
    };
}

void
SolidShell :: computeAlpha(FloatArray &answer, FloatArray &u)
{
    // compute alpha based on displacement update
    FloatArray deltaU;
    deltaU.beDifferenceOf(u,this->u_k);

    FloatMatrix KEE_inv = this->invKEE;
    FloatArray fE = this->fE;
    FloatMatrix KEC = this->KEC;
    FloatArray temp, deltaAlpha;
    temp.beProductOf(KEC, deltaU);
    temp.add(fE);
    deltaAlpha.beProductOf(KEE_inv,temp);

    FloatArray oldAlpha = this->alpha; // last converged values
    answer = this->alpha - deltaAlpha;

    // set current u-displcement
    this->u_k = u;
}

void
SolidShell :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    if ( this->EAS_type ) {
        FloatArray fC, fE, strain, u, vStrainC, vStrainE, vStress;
        FloatMatrix KEE, KBE, BC, BE, D;
        fE.clear();
        fC.clear();

        KEE.clear();

        this->computeVectorOf(VM_Total, tStep, u);
        FloatArray alpha;
        this->computeAlpha(alpha, u );
        StructuralCrossSection *cs = this->giveStructuralCrossSection();
        for ( auto &gp: *this->giveDefaultIntegrationRulePtr() ) {
            if ( nlGeometry == 0 ) {
                computeBmatrixAt(gp, BC);
                this->computeEASBmatrixAt(gp, BE);

                vStrainC.beProductOf(BC, u);
                vStrainE.beProductOf(BE, alpha);
                vStrainC.add(vStrainE);
                //this->computeStressVector(vStress, vStrainC, gp, tStep);
                cs->giveRealStresses(vStress, gp, vStrainC, tStep);

            } else if ( nlGeometry == 1 ) {  // First Piola-Kirchhoff stress

                this->computeBHmatrixAt(gp, BC);
                FloatArray vFC, VFE;
                vFC.beProductOf(BC, u);
                //vFE.beProductOf(BE, alpha);
                //vFC.add(VFE);
                cs->giveFirstPKStresses(vStress, gp, vFC, tStep);

            } else if ( nlGeometry == 2 ) {  // Second Piola-Kirchhoff stress

                this->computeBEmatrixAt(gp, BC, tStep);
                this->computeEASBmatrixAt(gp, BE);

                this->computeEVector(vStrainC, gp->giveNaturalCoordinates(), u);
                vStrainE.beProductOf(BE, alpha);
                vStrainC.add(vStrainE);
                cs->giveRealStresses(vStress, gp, vStrainC, tStep);
            }

            double dV  = this->computeVolumeAround(gp);

            // Compute nodal internal forces at nodes as f = B^T*Stress dV
            fC.plusProduct(BC, vStress, dV);
            fE.plusProduct(BE, vStress, dV);
        }
        this->fE = fE;

        FloatMatrix KEE_inv, KEC;
        KEE_inv = this->invKEE;
        KEC = this->KEC;

        FloatArray temp, f;
        temp.beProductOf(KEE_inv,fE);
        answer = fC;
        answer.plusProduct(KEC, temp, -1.0);
    } else {
        LSpace :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);
    }
}


void
SolidShell :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{
    if ( this->EAS_type ) {
        FloatArray fC, fE, strain, u, vStrainC, vStrainE, vStress;
        FloatMatrix KEE, KEC, KCC, KCC_geo, BC, BE, D, DBC, DBE;

        KEC.clear();
        KEE.clear();
        KCC.clear();

        StructuralCrossSection *cs = this->giveStructuralCrossSection();
        for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
             if ( nlGeometry == 0 ) {   
                this->computeBmatrixAt(gp, BC);
                this->computeEASBmatrixAt(gp, BE);
                double dV  = this->computeVolumeAround(gp);

                this->computeConstitutiveMatrixAt(D, rMode, gp, tStep);
                DBC.beProductOf(D, BC);
                DBE.beProductOf(D, BE);

                KCC.plusProductUnsym(BC, DBC, dV);
                KEC.plusProductUnsym(BE, DBC, dV);
                KEE.plusProductUnsym(BE, DBE, dV);

            } else if ( nlGeometry == 1 ) {
                this->computeBHmatrixAt(gp, BC);
                cs->giveStiffnessMatrix_dPdF(D, rMode, gp, tStep);

            } else if ( nlGeometry == 2 ) {
                this->computeBEmatrixAt(gp, BC, tStep);
                this->computeEASBmatrixAt(gp, BE);
                double dV  = this->computeVolumeAround(gp);

                this->computeConstitutiveMatrixAt(D, rMode, gp, tStep);
                DBC.beProductOf(D, BC);
                DBE.beProductOf(D, BE);

                KCC.plusProductUnsym(BC, DBC, dV);
                KEC.plusProductUnsym(BE, DBC, dV);
                KEE.plusProductUnsym(BE, DBE, dV);

                //TODO add geometric stiffness 
                this->computeGeometricStiffness(KCC_geo, gp, tStep);
                KCC.add(dV, KCC_geo);
            }
        }

        FloatMatrix KEE_inv;
        KEE_inv.beInverseOf(KEE);
        this->invKEE = KEE_inv;
        this->KEC = KEC;

        answer = KCC;

        FloatMatrix K, tempmat;
        tempmat.beProductOf(KEE_inv, KEC);
        answer.plusProductUnsym(KEC, tempmat, -1.0);
    } else {
        LSpace :: computeStiffnessMatrix(answer, rMode, tStep);
    }
}


void
SolidShell :: computeBEmatrixAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *tStep)
{
    FloatMatrix dN, F;
    FloatArray vF;
    // compute the derivatives of shape functions
    this->interpolation.evaldNdx( dN, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    FloatArray u;
    this->computeVectorOf(VM_Total, tStep, u); // solution vector

    this->computeFVector(vF, gp->giveNaturalCoordinates(), u);
    F.beMatrixForm(vF);

    answer.resize(6, 24);
    answer.zero();
#if 1
    for ( int i = 1, col = 0; i <= dN.giveNumberOfRows(); i++, col += 3 ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at(1,col + j) =  F.at(j,1) * dN.at(i, 1);
            answer.at(2,col + j) =  F.at(j,2) * dN.at(i, 2);
            answer.at(3,col + j) =  F.at(j,3) * dN.at(i, 3);
            answer.at(4,col + j) =  F.at(j,2) * dN.at(i, 3) + F.at(j,3) * dN.at(i, 2);
            answer.at(5,col + j) =  F.at(j,1) * dN.at(i, 3) + F.at(j,3) * dN.at(i, 1);
            answer.at(6,col + j) =  F.at(j,1) * dN.at(i, 2) + F.at(j,2) * dN.at(i, 1);      
        }
    }
#else

    FEInterpolation *interp = this->giveInterpolation();
    FloatMatrix dNdx; 
    const FloatArray &lCoords = gp->giveNaturalCoordinates();
    interp->evaldNdx( dNdx, lCoords, FEIElementGeometryWrapper(this) );
// Do ANS
    // Need to evaluate strain at four points
    FloatArray A = { 0.0, -1.0, 0.0};
    FloatArray B = { 1.0,  0.0, 0.0};
    FloatArray C = { 0.0,  1.0, 0.0};
    FloatArray D = {-1.0,  0.0, 0.0};

    FloatMatrix dNdxA, dNdxB, dNdxC, dNdxD;
    interp->evaldNdx( dNdxA, A, FEIElementGeometryWrapper(this) );
    interp->evaldNdx( dNdxB, B, FEIElementGeometryWrapper(this) );
    interp->evaldNdx( dNdxC, C, FEIElementGeometryWrapper(this) );
    interp->evaldNdx( dNdxD, D, FEIElementGeometryWrapper(this) );

    // Eval deformation gradients
    FloatArray vFA, vFB, vFC, vFD;
    this->computeFVector(vFA, A, u);
    this->computeFVector(vFB, B, u);
    this->computeFVector(vFC, C, u);
    this->computeFVector(vFD, D, u);
    FloatMatrix FA, FB, FC, FD;
    FA.beMatrixForm(vFA);
    FB.beMatrixForm(vFB);
    FC.beMatrixForm(vFC);
    FD.beMatrixForm(vFD);

    FloatMatrix dNdx0; 
    interp->evaldNdx( dNdx0, {lCoords.at(1), lCoords.at(2), 0.0}, FEIElementGeometryWrapper(this) );

    double NA = 0.5 * ( 1.0 - lCoords.at(2) );
    double NC = 0.5 * ( 1.0 + lCoords.at(2) );
    double NB = 0.5 * ( 1.0 - lCoords.at(1) );
    double ND = 0.5 * ( 1.0 + lCoords.at(1) );

    // E_xz = E(5) = N_A(gp) * B(A)(5,:) + N_C(gp) * B(C)(5,:)
    // E_yz = E(4) = N_B(gp) * B(B)(4,:) + N_D(gp) * B(D)(4,:)

    for ( int i = 1, col = 0; i <= dN.giveNumberOfRows(); i++, col += 3 ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at(1,col + j) =  F.at(j,1) * dN.at(i, 1);
            answer.at(2,col + j) =  F.at(j,2) * dN.at(i, 2);
            answer.at(3,col + j) =  F.at(j,3) * dN.at(i, 3);
            answer.at(6,col + j) =  F.at(j,1) * dN.at(i, 2) + F.at(j,2) * dN.at(i, 1);      

            answer.at(4,col + j) =  NB * ( FB.at(j,2) * dNdxB.at(i, 3) + FB.at(j,3) * dNdxB.at(i, 2) ) + ND * ( FD.at(j,2) * dNdxD.at(i, 3) + FD.at(j,3) * dNdxD.at(i, 2) )  ;
            answer.at(5,col + j) =  NA * ( FA.at(j,1) * dNdxA.at(i, 3) + FA.at(j,3) * dNdxA.at(i, 1) ) + NC * ( FC.at(j,1) * dNdxC.at(i, 3) + FC.at(j,3) * dNdxC.at(i, 1) );
        }
    }
#endif
}


void 
SolidShell :: computeFVector(FloatArray &answer, FloatArray &lCoords, FloatArray &ae)
{
    FloatMatrix B;
    this->computeBHmatrixAt(lCoords, B);
    answer.beProductOf(B, ae);

    answer.at(1) += 1.0;
    answer.at(2) += 1.0;
    answer.at(3) += 1.0;
}


void 
SolidShell :: computeEVector(FloatArray &answer, FloatArray &lCoords, FloatArray &u)
{
    FloatArray vF;
    this->computeFVector(vF, lCoords, u);
#if 0
    FloatMatrix F, E;
    F.beMatrixForm(vF);
    E.beTProductOf(F, F);
    E.at(1, 1) -= 1.0;
    E.at(2, 2) -= 1.0;
    E.at(3, 3) -= 1.0;
    E.times(0.5);
    answer.beSymVectorFormOfStrain(E);
#else
    FloatMatrix F, E;
    F.beMatrixForm(vF);
    E.beTProductOf(F, F);
    E.at(1, 1) -= 1.0;
    E.at(2, 2) -= 1.0;
    E.at(3, 3) -= 1.0;
    E.times(0.5);
    answer.beSymVectorFormOfStrain(E);     

    // Need to evaluate strain at four points
    FloatArray A = { 0.0, -1.0, 0.0};
    FloatArray B = { 1.0,  0.0, 0.0};
    FloatArray C = { 0.0,  1.0, 0.0};
    FloatArray D = {-1.0,  0.0, 0.0};

    // Eval deformation gradients
    FloatArray vFA, vFB, vFC, vFD;
    this->computeFVector(vFA, A, u);
    this->computeFVector(vFB, B, u);
    this->computeFVector(vFC, C, u);
    this->computeFVector(vFD, D, u);
    FloatMatrix FA, FB, FC, FD, EA, EB, EC, ED;
    FA.beMatrixForm(vFA);
    FB.beMatrixForm(vFB);
    FC.beMatrixForm(vFC);
    FD.beMatrixForm(vFD);

    EA.beTProductOf(FA,FA); EA.times(0.5);
    EB.beTProductOf(FB,FB); EB.times(0.5);
    EC.beTProductOf(FC,FC); EC.times(0.5);
    ED.beTProductOf(FD,FD); ED.times(0.5);

//     FloatMatrix dNdx0; 
//     interp->evaldNdx( dNdx0, {lCoords.at(1), lCoords.at(2), 0.0}, FEIElementGeometryWrapper(this) );

    double NA = 0.5 * ( 1.0 - lCoords.at(2) );
    double NC = 0.5 * ( 1.0 + lCoords.at(2) );
    double NB = 0.5 * ( 1.0 - lCoords.at(1) );
    double ND = 0.5 * ( 1.0 + lCoords.at(1) );

    answer.at(4) = NB * EB.at(2,3) + ND * ED.at(2,3);
    answer.at(5) = NA * EA.at(2,3) + NC * EC.at(2,3);
#endif
}


void 
SolidShell :: computeGeometricStiffness(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(24,24);
    answer.zero();

    StructuralMaterialStatus *matStat = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() );
    if ( !matStat ) { // no matstat created yet
        return; 
    }
    FloatArray vS;
    vS = matStat->giveStressVector();

    FloatMatrix S, dNdx, temp, G_ij;
    S.beMatrixFormOfStress(vS);
    this->interpolation.evaldNdx( dNdx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    temp.beProductTOf(S,dNdx);
    G_ij.beProductOf(dNdx, temp);
    FloatMatrix K_ij(3,3);
    K_ij.zero();

    for ( int i = 1; i <= dNdx.giveNumberOfRows(); i++ ) {
        for ( int j = i; j <= dNdx.giveNumberOfRows(); j++ ) {
            K_ij.at(1,1) = K_ij.at(2,2) = K_ij.at(3,3) = G_ij.at(i,j);
            answer.assemble(K_ij, { (i-1)*3 + 1, (i-1)*3 + 2, (i-1)*3 + 3 }, { (j-1)*3 + 1, (j-1)*3 + 2, (j-1)*3 + 3 } );
        }
    }
    answer.symmetrized();
}

} // end namespace oofem
