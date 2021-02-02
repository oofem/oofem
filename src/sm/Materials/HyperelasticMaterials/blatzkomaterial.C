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
 *               Copyright (C) 1993 - 2020   Borek Patzak
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

#include "blatzkomaterial.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Material(BlatzKoMaterial);


void
BlatzKoMaterial::initializeFrom(InputRecord &ir)
{
    StructuralMaterial::initializeFrom(ir);
    IR_GIVE_FIELD(ir, mu, _IFT_BlatzKoMaterial_mu);
}


FloatMatrixF< 9, 9 >
BlatzKoMaterial::give3dMaterialStiffnessMatrix_dPdF(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
{
    //retrieve deformation from status and convert to matrix form
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    FloatArrayF< 9 >vF(status->giveTempFVector() );
    Tensor2_3d F(vF);

    double I2 = compute_I2_C_from_F(F);
    double I3 = compute_I3_C_from_F(F);

    Tensor2_3d dI2_dF = compute_dI2_C_dF(F);
    Tensor2_3d dI3_dF = compute_dI3_C_dF(F);

    Tensor4_3d A;
    A(i_3, j_3, k_3, l_3) = 0.5 * mu / I3  * ( compute_d2I2_C_dF2(F)(i_3, j_3, k_3, l_3) - 1. / I3 * ( dI2_dF(i_3, j_3) * dI3_dF(k_3, l_3) + dI3_dF(i_3, j_3) * dI2_dF(k_3, l_3) ) + ( sqrt(I3) - I2 / I3 ) * compute_d2I3_C_dF2(F)(i_3, j_3, k_3, l_3) + ( 2. * I2 / I3 / I3 - 0.5 / sqrt(I3) ) *  dI3_dF(i_3, j_3) * dI3_dF(k_3, l_3) );

    return A.to_voigt_form();
}

FloatArrayF< 9 >
BlatzKoMaterial::giveFirstPKStressVector_3d(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    // compute the first Piola-Kirchhoff
    Tensor2_3d F(vF), P;
    double I2 = compute_I2_C_from_F(F);
    double I3 = compute_I3_C_from_F(F);
    P(i_3, j_3) =  0.5 * mu * ( compute_dI2_C_dF(F)(i_3, j_3) / I3 + ( 1. / sqrt(I3) -  I2 / I3 / I3 )  * compute_dI3_C_dF(F)(i_3, j_3) );
    //
    auto vP = P.to_voigt_form();
    // update gp
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(vP);
    return vP;
}

MaterialStatus *
BlatzKoMaterial::CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(gp);
}
};
