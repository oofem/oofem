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

#ifndef ogdencompressiblematerial_h
#define ogdencompressiblematerial_h

#include "../sm/Materials/structuralmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "basehyperelasticmaterial.h"

///@name Input fields for OgdenMaterial
//@{
#define _IFT_OgdenCompressibleMaterial_Name "ogdencompressiblemat"

#define _IFT_OgdenCompressibleMaterial_alpha "alpha"
#define _IFT_OgdenCompressibleMaterial_mu "mu"
//@}

namespace oofem {
/**
 * This class implements Compressible Ogden material.
 *
 * @author Martin Horák, nitramkaroh@seznam.cz
 *
 * @Note References: R.W. Ogden: Non-Linear Elastic Deformations, de Souza Neto, Peric, Owen: Computational Methods for Plasticity: Theory and Applications
 *
 * Free energy is considered as:
 * @f[
 * \rho_0 \psi = \sum_I^N \frac{mu_I}{alpha_I}(\lambda_1^{\alpha_I}+\lambda_2^{\alpha_I}+\lambda_3^{\alpha_I}-3) + \frac{1}{2} K[ln(J)]^2
 * @f]
 * @f$ \alpha_I @f$, @f$ mu_I @f$, and @f$K@f$ are material parameters.
 *
 * @f$ \lambda_1@f$, @f$ \lambda_2 @f$, and @f$ \lambda_3 @f$ are principal stretches
 *
 * Compressible Neo-Hookean model is obtained by setting @f$N = 1@f$, @f$\alpha_1 = 2@f$
 *
 * Compressible Mooney-Rivlin model is obtained by setting @f$N = 2@f$, @f$\alpha_1 = 2@f$, and @f$\alpha_2 = -2@f$.
 */
class OgdenCompressibleMaterial : public StructuralMaterial, public BaseHyperElasticMaterial
{
protected:
    /// Array of Exponents alpha
    FloatArray alpha;
    /// Array of Material parameters mu
    FloatArray mu;
    /// Number of material parameters in the arrays mu(alpha)
    int N;


public:
    OgdenCompressibleMaterial(int n, Domain *d);


    FloatMatrixF< 6, 6 >give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }

    FloatArrayF< 6 >giveRealStressVector_3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) const override { OOFEM_ERROR("not implemented, this material is designed for large strains only"); }
    FloatMatrixF< 9, 9 >give3dMaterialStiffnessMatrix_dPdF(MatResponseMode,
                                                           GaussPoint *gp,
                                                           TimeStep *tStep) const override;

    FloatArrayF< 9 >giveFirstPKStressVector_3d(const FloatArrayF< 9 > &vF, GaussPoint *gp, TimeStep *tStep) const override;

    virtual void initializeFrom(InputRecord &ir) override;
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const override;
    virtual const char *giveInputRecordName() const override { return _IFT_OgdenCompressibleMaterial_Name; }
    virtual const char *giveClassName() const override { return "OgdenCompressibleMaterial"; }

private:
    /**
     * Computes the deviatoric part of the second Piola-Kirchhoff stress tensor from the left Cauchy-Green tensor
     * @param C the left Cauchy-Green deformation tensor
     * @return the second Piola-Kirchhoff stress tensor
     **/
    Tensor2_3d giveDeviatoricSecondPKStressVector_3d(const Tensor2_3d &C) const;
};
} // end namespace oofem
#endif
