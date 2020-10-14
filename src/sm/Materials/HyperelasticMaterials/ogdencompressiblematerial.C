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

#include "ogdencompressiblematerial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"
#include "mathfem.h"

namespace oofem {
REGISTER_Material(OgdenCompressibleMaterial);

OgdenCompressibleMaterial::OgdenCompressibleMaterial(int n, Domain *d) : StructuralMaterial(n, d), BaseHyperElasticMaterial()
{ }

FloatArrayF< 9 >
OgdenCompressibleMaterial::giveFirstPKStressVector_3d(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const
// returns 9 components of the first Piola-Kirchhoff stress corresponding to the given deformation gradient
{
    //store deformation gradient into tensor
    Tensor2_3d F(vF), P;
    Tensor2_3d C, S;
    // @todo:replace C by  Tensor2sym_3d C = F^F, where ^ is leads to symmetric tensor
    C(i_3, j_3) = F(k_3, i_3) * F(k_3, j_3);
    S = this->giveDeviatoricSecondPKStressVector_3d(C);
    P(i_3, j_3) = F(i_3, k_3) * S(k_3, j_3) + this->compute_dVolumetricEnergy_dF(F)(i_3, j_3);
    auto vP = P.to_voigt_form();
    // update gp
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(vP);

    //
    return vP;
}


FloatMatrixF< 9, 9 >
OgdenCompressibleMaterial::give3dMaterialStiffnessMatrix_dPdF(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    double J;
    FloatArrayF< 9 >vF(status->giveTempFVector() );
    Tensor2_3d F(vF), P, C, iC, Sdev;
    Tensor4_3d dSdE, dPdF, iCiC;

    J = F.compute_determinant();
    //
    // @todo:replace C by  Tensor2sym_3d and C = F^F, where ^ leads to symmetric tensor
    C(i_3, j_3) = F(k_3, i_3) * F(k_3, j_3);
    // compute eigen values and eigen vectors of C
    auto [ eVals, eVecs ] = C.eigs();
    //compute inverse of the right Cauchy-Green tensor
    iC = C.compute_inverse();
    //compute symetric dyadic product
    iCiC(i_3, j_3, k_3, l_3) = ( iC(i_3, k_3) * iC(j_3, l_3) ) + ( iC(i_3, l_3) * iC(j_3, k_3) );
    for ( int i = 1; i <= N; i++ ) {
        double l1a, l2a, l3a, Ja;
        //
        Ja = pow(J, -alpha.at(i) / 3.);
        l1a = pow(eVals.at(1), alpha.at(i) / 2.);
        l2a = pow(eVals.at(2), alpha.at(i) / 2.);
        l3a = pow(eVals.at(3), alpha.at(i) / 2.);
        //
        double I1a = ( l1a + l2a + l3a ) / 3.;
        //
        Tensor2_3d powC, iSdev;
        powC = Tensor2_3d::computeTensorPowerFromEigs(eVals, eVecs, ( alpha.at(i) - 2. ) / 2.);
        iSdev(i_3, j_3) = mu.at(i) * Ja * ( powC(i_3, j_3) - I1a * iC(i_3, j_3) );
        //
        dSdE(i_3, j_3, k_3, l_3) += 2. * mu.at(i) * Ja *  ( Tensor2_3d::compute_dCm_dC_fromEigs( ( alpha.at(i) - 2. ) / 2., eVals, eVecs )(i_3, j_3, k_3, l_3) - alpha.at(i) / 6. * iC(i_3, j_3) * powC(k_3, l_3) + I1a / 2. *  iCiC(i_3, j_3, k_3, l_3) ) - alpha.at(i) / 3. *  iSdev(i_3, j_3) * iC(k_3, l_3);
        Sdev(i_3, j_3) += iSdev(i_3, j_3);
    }

    // transform the second elasticity tensor to the first elasticity tensor
    FloatMatrixF< 9, 9 >mdSdE;
    mdSdE = dSdE.to_voigt_form();
    FloatMatrix m2dSdE(mdSdE);
    Tensor2_3d I = Tensor2_3d::UnitTensor();
    dPdF(i_3, j_3, k_3, l_3) = I(i_3, k_3) * Sdev(j_3, l_3) + F(i_3, m_3) * dSdE(m_3, j_3, n_3, l_3) * F(k_3, n_3);
    // add volumetric part
    dPdF(i_3, j_3, k_3, l_3) += this->compute_d2VolumetricEnergy_dF2(F)(i_3, j_3, k_3, l_3);

    FloatMatrix mdP(dPdF.to_voigt_form() );
    return dPdF.to_voigt_form();
}





Tensor2_3d
OgdenCompressibleMaterial::giveDeviatoricSecondPKStressVector_3d(const Tensor2_3d &C) const
{
    // compute eigen values and eigen vectors of C
    auto [ eVals, eVecs ] = C.eigs();
    //
    double J = sqrt(C.compute_determinant() );
    Tensor2_3d S, iC;
    //compute inverse of the right Cauchy-Green tensor
    iC = C.compute_inverse();
    for ( int i = 1; i <= N; i++ ) {
        double l1a, l2a, l3a, Ja;
        Tensor2_3d Si;
        //
        l1a = pow(eVals.at(1), alpha.at(i) / 2.);
        l2a = pow(eVals.at(2), alpha.at(i) / 2.);
        l3a = pow(eVals.at(3), alpha.at(i) / 2.);
        Ja = pow(J, -alpha.at(i) / 3.);
        double I1a = ( l1a + l2a + l3a ) / 3.;
        S(i_3, j_3) += mu.at(i) * Ja * ( Tensor2_3d::computeTensorPowerFromEigs(eVals, eVecs, ( alpha.at(i) - 2. ) / 2.)(i_3, j_3) - I1a * iC(i_3, j_3) );
    }
    return S;
}





MaterialStatus *
OgdenCompressibleMaterial::CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(gp);
}


void
OgdenCompressibleMaterial::initializeFrom(InputRecord &ir)
{
    StructuralMaterial::initializeFrom(ir);
    BaseHyperElasticMaterial::initializeFrom(ir);

    IR_GIVE_FIELD(ir, alpha, _IFT_OgdenCompressibleMaterial_alpha);
    IR_GIVE_FIELD(ir, mu, _IFT_OgdenCompressibleMaterial_mu);

    N = alpha.giveSize();
    int M = mu.giveSize();
    if ( N != M ) {
        OOFEM_ERROR("Inconsistent size of alpha and mu");
    }
}
} // end namespace oofem
