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
 *               Copyright (C) 1993 - 2013   Borek Patzak
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

#ifndef j2mat_h
#define j2mat_h

#include "sm/Materials/ConcreteMaterials/mplasticmaterial2.h"

///@name Input fields for J2Mat
//@{
#define _IFT_J2Mat_Name "j22mat"
#define _IFT_J2Mat_ry "ry"
#define _IFT_J2Mat_khm "khm"
#define _IFT_J2Mat_ihm "ihm"
#define _IFT_J2Mat_rma "rma"
//@}

namespace oofem {
class Domain;

/**
 * This class implements a isotropic plastic linear material (J2 plasticity condition is used).
 * in a finite element problem.
 * Both kinematic and isotropic hardening is supported.
 */
class J2Mat : public MPlasticMaterial2
{
protected:
    int kinematicHardeningFlag = 0, isotropicHardeningFlag = 0;
    double kinematicModuli = 0., isotropicModuli = 0.;
    double k = 0.;

public:
    J2Mat(int n, Domain * d);

    void initializeFrom(InputRecord &ir) override;
    const char *giveInputRecordName() const override { return _IFT_J2Mat_Name; }
    const char *giveClassName() const override { return "J2Mat"; }

    int giveSizeOfFullHardeningVarsVector() const override;
    int giveSizeOfReducedHardeningVarsVector(GaussPoint *gp) const override;
    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

protected:
    int giveMaxNumberOfActiveYieldConds(GaussPoint *gp) const override { return 2; }

    double computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector,
                               const FloatArray &strainSpaceHardeningVars) const override;

    void computeStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp,
                                     const FloatArray &stressVector, const FloatArray &strainSpaceHardeningVars) const override;

    void computeStrainHardeningVarsIncrement(FloatArray &answer, GaussPoint *gp,
                                             const FloatArray &stress, const FloatArray &dlambda,
                                             const FloatArray &dplasticStrain, const IntArray &activeConditionMap) const override;
    void computeKGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, FloatArray &fullStressVector,
                                const FloatArray &strainSpaceHardeningVariables) const override;

    void computeReducedHardeningVarsSigmaGradient(FloatMatrix &answer, GaussPoint *gp, const IntArray &activeConditionMap,
                                                  const FloatArray &fullStressVector,
                                                  const FloatArray &strainSpaceHardeningVars,
                                                  const FloatArray &gamma) const override;
    void computeReducedHardeningVarsLamGradient(FloatMatrix &answer, GaussPoint *gp, int actSurf,
                                                const IntArray &activeConditionMap,
                                                const FloatArray &fullStressVector,
                                                const FloatArray &strainSpaceHardeningVars,
                                                const FloatArray &gamma) const override;
    int hasHardening() const override;

    void  computeReducedSSGradientMatrix(FloatMatrix &gradientMatrix,  int i, GaussPoint *gp, const FloatArray &fullStressVector,
                                         const FloatArray &strainSpaceHardeningVariables) const override;
    void  computeReducedSKGradientMatrix(FloatMatrix &gradientMatrix,  int i, GaussPoint *gp, const FloatArray &fullStressVector,
                                         const FloatArray &strainSpaceHardeningVariables) const override;

    // auxiliary function
    static double computeJ2InvariantAt(const FloatArray &stressVector);
    double giveIsotropicHardeningVar(GaussPoint *gp, const FloatArray &strainSpaceHardeningVars) const;
    void giveStressBackVector(FloatArray &answer, GaussPoint *gp, const FloatArray &strainSpaceHardeningVars) const;
};
} // end namespace oofem
#endif // j2mat_h
