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

#include "simopistermaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Material(SimoPisterMaterial);

SimoPisterMaterial::SimoPisterMaterial(int n, Domain *d) : StructuralMaterial(n, d), BaseHyperElasticMaterial()
{ }

FloatMatrixF< 9, 9 >
SimoPisterMaterial::give3dMaterialStiffnessMatrix_dPdF(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const
// returns the 9x9 tangent stiffness matrix - dP/dF
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );
    FloatArrayF< 9 >vF(status->giveTempFVector() );
    Tensor2_3d F(vF);
    auto [ J, cofF ] = F.compute_determinant_and_cofactor();
    Tensor4_3d A;
    A(i_3, j_3, k_3, l_3) = G * ( 0.5 * this->compute_d2I1_C_dF2(F)(i_3, j_3, k_3, l_3) + 1 / J / J * cofF(i_3, j_3) * cofF(k_3, l_3) - 1. / J * F.compute_tensor_cross_product()(i_3, j_3, k_3, l_3) )  + this->compute_d2VolumetricEnergy_dF2(F)(i_3, j_3, k_3, l_3);
    return A.to_voigt_form();
}

FloatArrayF< 9 >
SimoPisterMaterial::giveFirstPKStressVector_3d(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const
// returns 9 components of the first piola kirchhoff stress corresponding to the given deformation gradinet
{
    StructuralMaterialStatus *status = static_cast< StructuralMaterialStatus * >( this->giveStatus(gp) );

    Tensor2_3d F(vF), P;
    auto [ J, cofF ] = F.compute_determinant_and_cofactor();
    // compute the first Piola-Kirchhoff
    P(i_3, j_3) =  G * ( 0.5 * this->compute_dI1_C_dF(F)(i_3, j_3) - 1. / J * cofF(i_3, j_3) ) + this->compute_dVolumetricEnergy_dF(F)(i_3, j_3);
    auto vP = P.to_voigt_form();
    // update gp
    status->letTempFVectorBe(vF);
    status->letTempPVectorBe(vP);
    //
    return vP;
}


MaterialStatus *
SimoPisterMaterial::CreateStatus(GaussPoint *gp) const
{
    return new StructuralMaterialStatus(gp);
}


void
SimoPisterMaterial::initializeFrom(InputRecord &ir)
{
    StructuralMaterial::initializeFrom(ir);
    BaseHyperElasticMaterial::initializeFrom(ir);
    IR_GIVE_FIELD(ir, G, _IFT_SimoPisterMaterial_g);
}
} // end namespace oofem
