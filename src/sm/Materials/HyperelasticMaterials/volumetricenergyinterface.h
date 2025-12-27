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


#ifndef volumetricenergyinterface_h
#define volumetricenergyinterface_h

#include "tensor/tensor.h"

///@name Input fields for VolumetricEnergyInterface
//@{
#define _IFT_VolumetricEnergyInterface_k "k"
#define _IFT_VolumetricEnergyInterface_type "type"
//@}

namespace oofem {
  class InputRecord;

/**
 * This class implements Volumetric form of energy fo Compressible Hyperelastic materials.
 *
 * @author Martin Horák, nitramkaroh@seznam.cz
 * 
 * References: R.W. Ogden: Non-Linear Elastic Deformations,
 * de Souza Neto, Peric, Owen: Computational Methods for Plasticity: Theory and Applications
 *
 * Free energy is considered as:
 * @f[
 * \rho_0 \psi = Kf(J)
 */
  
class VolumetricEnergyInterface
{
  
 protected:
  double K;
  int type = 1;

 public:
  VolumetricEnergyInterface(){;}
  virtual ~VolumetricEnergyInterface(){;}

  Tensor2_3d compute_dVolumetricEnergy_dF(const Tensor2_3d &F) const;
  Tensor4_3d compute_d2VolumetricEnergy_dF2(const Tensor2_3d &F) const;
  void initializeFrom(InputRecord &ir);
};
 
} // end namespace oofem
#endif
