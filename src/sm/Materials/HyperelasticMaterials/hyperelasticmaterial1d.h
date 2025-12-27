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

#ifndef hyperelasticmaterial1d_h
#define hyperelasticmaterial1d_h

#include "sm/Materials/structuralms.h"
#include "sm/Materials/structuralmaterial.h"


///@name Input fields for MooneyRivlinMaterial
//@{
#define _IFT_HyperelasticMaterial1d_Name "hyperelasticmaterial1d"
#define _IFT_HyperelasticMaterial1d_type "mtype"
#define _IFT_HyperelasticMaterial1d_E "e"
//@}

namespace oofem {
/**
 * This class implements several 1d hyperelastic materials.
 * Including now Biot and Saint Venant-Kirchhoff materials
 *
 * @author Martin Horák, nitramkaroh@seznam.cz
 *
 */
class HyperelasticMaterial1d : public StructuralMaterial
{
protected:
  double E = 0;

  enum HyperElasticMaterialType {
        HEMT_Biot=0,
        HEMT_StVenantKirchhoff=1,
        HEMT_Unknown = 100
    };
    /// Parameter specifying the definition of used constitutive model
    HyperElasticMaterialType hyperelasticMaterialType = HEMT_Unknown;


public:
    HyperelasticMaterial1d(int n, Domain *d);

    void initializeFrom(InputRecord &ir) override;

    FloatMatrixF< 6, 6 >give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override { OOFEM_ERROR("not implemented, this material supports only 1d case"); }

    FloatArrayF< 6 >giveRealStressVector_3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) const override { OOFEM_ERROR("not implemented, this material supports only 1d case"); }
    FloatMatrixF< 9, 9 >give3dMaterialStiffnessMatrix_dPdF(MatResponseMode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep) const override{OOFEM_ERROR("not implemented, this material supports only 1d case");}

    FloatArrayF< 9 >giveFirstPKStressVector_3d(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const override{OOFEM_ERROR("not implemented, this material supports only 1d case");}


    FloatArrayF< 1 > giveFirstPKStressVector_1d(const FloatArrayF< 1 > &vF, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF< 1, 1 > give1dStressStiffnessMatrix_dPdF(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;
   

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;

    const char *giveInputRecordName() const override { return _IFT_HyperelasticMaterial1d_Name; }
    const char *giveClassName() const override { return "HyperelasticMaterial1d"; }
};
} // end namespace oofem
#endif
