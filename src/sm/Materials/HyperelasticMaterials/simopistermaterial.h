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


#ifndef simopistermaterial_h
#define simopistermaterial_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "basehyperelasticmaterial.h"


///@name Input fields for SimoPisterMaterial
//@{
#define _IFT_SimoPisterMaterial_Name "simopistermat"
#define _IFT_SimoPisterMaterial_g "g"
//@}

namespace oofem {
/**
 * Free energy is considered as:
 * \f$[
 * \rho_0 \psi = U(J) + G( 0.5 * I_1 - ln J) ]\f$
 *  This form of energy corresponds to a neo-Hookean material which is extended to the compressible * range by adding an extra function depending on J.
 * @author Martin Horak, nitramkaroh@seznam.cz
 * @note Reference: article{simo1984remarks,
 * title={Remarks on rate constitutive equations for finite deformation problems: computational implications},
 * author={Simo, Juan C and Pister, Karl S},
 * journal={Computer Methods in Applied Mechanics and Engineering},
 * volume={46},
 * number={2},
 * pages={201--215},
 * year={1984},
 * publisher={Elsevier}
 * }
 *
 *
 */
class SimoPisterMaterial : public StructuralMaterial, public BaseHyperElasticMaterial
{
protected:
    double G = 0.; ///< Shear modulus.

public:
    SimoPisterMaterial(int n, Domain *d);

    void initializeFrom(InputRecord &ir) override;

    FloatMatrixF< 6, 6 >give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }

    FloatArrayF< 6 >giveRealStressVector_3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) const override { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }


    FloatMatrixF< 9, 9 >give3dMaterialStiffnessMatrix_dPdF(MatResponseMode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep) const override;

    FloatArrayF< 9 >giveFirstPKStressVector_3d(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;

    const char *giveInputRecordName() const override { return _IFT_SimoPisterMaterial_Name; }
    const char *giveClassName() const override { return "SimoPisterMaterial"; }
};
} // end namespace oofem
#endif
