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

#ifndef mooneyrivlinmaterial_h
#define mooneyrivlinmaterial_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "basehyperelasticmaterial.h"


///@name Input fields for MooneyRivlinMaterial
//@{
#define _IFT_MooneyRivlinCompressibleMaterial_Name "mooneyrivlincompressiblemat"
#define _IFT_MooneyRivlinCompressibleMaterial_c1 "c1"
#define _IFT_MooneyRivlinCompressibleMaterial_c2 "c2"
//@}

namespace oofem {
/**
 * This class implements Compressible Mooney - Rivlin material.
 *
 * @author Martin Horák, nitramkaroh@seznam.cz
 *
 * @note References: R.W. Ogden: Non-Linear Elastic Deformations,
 * de Souza Neto, Peric, Owen: Computational Methods for Plasticity: Theory and Applications
 *
 * Free energy is considered as:
 * \f$[
 * \rho_0 \psi = C_1(\bar{I}_1 - 3) + C_2(\bar{I}_2-3) + \frac{1}{2} K[ln(J)]^2
 * ]\f$
 * \f$C_1\f$, \f$C_2\f$ , and \f$K\f$ are material parameters.
 *
 * \f$ \bar{I}_1 = J^{-2/3}I_1 \f$, where \f$I_1 \f$ is the first invariant of the left Cauchy-Green deformation tensor C.
 *
 * \f$ \bar{I}_2 = J^{-4/3}I_2 \f$, where \f$I_2\f$ is the second invariant of the left Cauchy-Green deformation tensor C
 *
 * Compressible Neo-Hookean model is obtained by setting @f$C_2 = 0@f$
 */
class MooneyRivlinCompressibleMaterial : public StructuralMaterial, public BaseHyperElasticMaterial
{
protected:
    // Material parameters
    double C1;
    double C2;


public:
    MooneyRivlinCompressibleMaterial(int n, Domain *d);

    void initializeFrom(InputRecord &ir) override;

    FloatMatrixF< 6, 6 >give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }

    FloatArrayF< 6 >giveRealStressVector_3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) const override { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }
    FloatMatrixF< 9, 9 >give3dMaterialStiffnessMatrix_dPdF(MatResponseMode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep) const override;

    FloatArrayF< 9 >giveFirstPKStressVector_3d(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    const char *giveInputRecordName() const override { return _IFT_MooneyRivlinCompressibleMaterial_Name; }
    const char *giveClassName() const override { return "MooneyRivlinCompressibleMaterial"; }
};
} // end namespace oofem
#endif
