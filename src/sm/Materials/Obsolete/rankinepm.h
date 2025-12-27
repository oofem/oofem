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

#ifndef rankinepm_h
#define rankinepm_h

#include "mplasticmaterial.h"

///@name Input fields for RankinePlasticMaterial
//@{
#define _IFT_RankinePlasticMaterial_Name "rankine"
#define _IFT_RankinePlasticMaterial_ry "ry"
//@}

namespace oofem {
class Domain;

/**
 * This class implements a isotropic  plastic linear material (J2 plasticity condition is used)
 * in a finite element problem.
 */
class RankinePlasticMaterial : public MPlasticMaterial
{
protected:
    /// Yield value.
    double k = 0.;

public:
    RankinePlasticMaterial(int n, Domain * d);

    void initializeFrom(InputRecord &ir) override;

    const char *giveInputRecordName() const override { return _IFT_RankinePlasticMaterial_Name; }
    const char *giveClassName() const override { return "RankinePlasticMaterial"; }

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;

protected:

    //
    // yield(YC-like functions) and loading(LC-like functions) criteria specific section
    //

    double computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector,
                               const FloatArray &stressSpaceHardeningVars) const override;

    void computeHardeningReducedModuli(FloatMatrix &answer,
                                       GaussPoint *gp,
                                       const FloatArray &strainSpaceHardeningVariables,
                                       TimeStep *tStep) const override;
    void computeStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &stressVector,
                                     const FloatArray &stressSpaceHardeningVars) const override;
    void computeStressSpaceHardeningVarsReducedGradient(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp,
                                                        const FloatArray &stressVector,
                                                        const FloatArray &stressSpaceHardeningVars) const override;
    int hasHardening() const override { return 0; }
    void computeReducedGradientMatrix(FloatMatrix &answer, int isurf,
                                      GaussPoint *gp,
                                      const FloatArray &stressVector,
                                      const FloatArray &stressSpaceHardeningVars) const override;

    void computeStressSpaceHardeningVars(FloatArray &answer, GaussPoint *gp,
                                         const FloatArray &strainSpaceHardeningVariables) const override;
};
} // end namespace oofem
#endif // rankinepm_h
