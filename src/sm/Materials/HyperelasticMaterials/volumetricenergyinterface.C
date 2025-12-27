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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "volumetricenergyinterface.h"
#include "tensor/tensor.h"
#include "dynamicinputrecord.h"



namespace oofem {

  
Tensor2_3d
VolumetricEnergyInterface :: compute_dVolumetricEnergy_dF(const Tensor2_3d &F) const

{
  // compute jacobian and its logarithm
  Tensor2_3d dVolumetricEnergy_dF;
  if(type == 1) {
    auto J = compute_determinant(F);
    auto lnJ = log(J);
    dVolumetricEnergy_dF(i_3,j_3) =  K * lnJ / J * compute_dJ_dF(F)(i_3,j_3);
  }
  return dVolumetricEnergy_dF;
}
  

Tensor4_3d
VolumetricEnergyInterface :: compute_d2VolumetricEnergy_dF2(const Tensor2_3d &F) const
{
  // compute jacobian and its logarithm

  Tensor4_3d d2VolumetricEnergy_dF2;
  if(type == 1) {
    auto [J, cofF] = compute_determinant_and_cofactor(F);
    auto lnJ = log(J);
    d2VolumetricEnergy_dF2(i_3,j_3,k_3,l_3) =  K * (1.-lnJ) / J / J * cofF(i_3, j_3) * cofF(k_3, l_3) + K * lnJ / J * compute_tensor_cross_product(F)(i_3,j_3,k_3,l_3);
   
  }
  return  d2VolumetricEnergy_dF2;
}
  
  


void
VolumetricEnergyInterface :: initializeFrom(InputRecord &ir)
{
  IR_GIVE_FIELD(ir, K, _IFT_VolumetricEnergyInterface_k);
  IR_GIVE_OPTIONAL_FIELD(ir, type, _IFT_VolumetricEnergyInterface_type)
}

} // end namespace oofem
