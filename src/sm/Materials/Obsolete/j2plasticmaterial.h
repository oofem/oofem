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

#ifndef j2plasticmaterial_h
#define j2plasticmaterial_h

#include "plasticmaterial.h"

///@name Input fields for J2plasticMaterial
//@{
#define _IFT_J2plasticMaterial_Name "j2mat"
#define _IFT_J2plasticMaterial_ry "ry"
#define _IFT_J2plasticMaterial_khm "khm"
#define _IFT_J2plasticMaterial_ihm "ihm"
//@}

namespace oofem {
class Domain;

/**
 * This class implements a isotropic  plastic linear material (J2 plasticity condition is used)
 * in a finite element problem. A material
 * is an attribute of a domain. It is usually also attribute of many elements.
 */
class J2plasticMaterial : public PlasticMaterial
{
protected:
    int kinematicHardeningFlag = 0, isotropicHardeningFlag = 0;
    double kinematicModuli = 0., isotropicModuli = 0.;
    //double E = 0., nu = 0.; // isotropic material constants
    double k = 0.;

public:
    J2plasticMaterial(int n, Domain * d);

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;
    const char *giveInputRecordName() const override { return _IFT_J2plasticMaterial_Name; }
    const char *giveClassName() const override { return "J2plasticMaterial"; }

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;

protected:

    //
    // yield(YC-like functions) and loading(LC-like functions) criteria specific section
    //

    FloatArray *ComputeStressSpaceHardeningVars(GaussPoint *gp,
                                                FloatArray *strainSpaceHardeningVariables) const override;
    double computeYieldValueAt(GaussPoint *gp, FloatArray *stressVector,
                               FloatArray *stressSpaceHardeningVars) const override;
    void computeHardeningReducedModuli(FloatMatrix &answer, GaussPoint *gp,
                                       FloatArray *strainSpaceHardeningVariables,
                                       TimeStep *tStep) const override;
    FloatArray *ComputeStressGradient(GaussPoint *gp, FloatArray *stressVector,
                                      FloatArray *stressSpaceHardeningVars) const override;
    FloatArray *ComputeStressSpaceHardeningVarsReducedGradient(GaussPoint *gp,
                                                               FloatArray *stressVector,
                                                               FloatArray *stressSpaceHardeningVars) const override;
    int hasHardening() const override;
    void computeReducedGradientMatrix(FloatMatrix &answer, GaussPoint *gp,
                                      const FloatArray &stressVector,
                                      const FloatArray &stressSpaceHardeningVars) const override;
    void computeTrialStressIncrement(FloatArray &answer, GaussPoint *gp,
                                     const FloatArray &strainIncrement, TimeStep *tStep) const override;
    void compute3dElasticModuli(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep) const override;

    // auxiliary function
    double computeJ2InvariantAt(FloatArray *answer) const;
    int giveSizeOfFullHardeningVarsVector() const override;
    int giveSizeOfReducedHardeningVarsVector(GaussPoint *gp) const override;
    double giveIsotropicHardeningVar(FloatArray *stressSpaceHardeningVars) const;
    void giveStressBackVector(FloatArray &answer, const FloatArray &stressSpaceHardeningVars) const;
};
} // end namespace oofem
#endif // j2plasticmaterial_h
