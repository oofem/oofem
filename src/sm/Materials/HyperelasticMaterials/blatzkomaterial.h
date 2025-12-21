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

#ifndef blatzkomaterial_h
#define blatzkomaterial_h

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "basehyperelasticmaterial.h"


///@name Input fields for MooneyRivlinMaterial
//@{
#define _IFT_BlatzKoMaterial_Name "blatzkomat"
#define _IFT_BlatzKoMaterial_mu "mu"
//@}

namespace oofem {
/**
 * This class implements the Blatz-Ko model for porous materials.
 * The Blatz-Ko strain energy density function is useful for modeling compressible foam type rubbers and can be expressed as:
 *\f[ W = \mu\left(\frac{I_2}{I_3}+2\sqrt{I_3}-5\right)\f]
 * The model contains only one parametr, shear modulus, Poisson's ratio is fixed to 0.25
 * @author Ondrej Faltus
 * @author Martin Horak
 *
 */
class BlatzKoMaterial : public StructuralMaterial, public BaseHyperElasticMaterial
{
protected:
    /// Shear modulus
    double mu;

public:
    BlatzKoMaterial(int n, Domain *d) : StructuralMaterial(n, d), BaseHyperElasticMaterial() { };

    void initializeFrom(InputRecord &ir) override;

    FloatMatrixF< 6, 6 >give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }

    FloatArrayF< 6 >giveRealStressVector_3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) const override { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }

    FloatMatrixF< 9, 9 >give3dMaterialStiffnessMatrix_dPdF(MatResponseMode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep) const override;

    FloatArrayF< 9 >giveFirstPKStressVector_3d(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;

    const char *giveInputRecordName() const override { return _IFT_BlatzKoMaterial_Name; }
    const char *giveClassName() const override { return "BlatzKoMaterial"; }
};
} // end namespace oofem
#endif
