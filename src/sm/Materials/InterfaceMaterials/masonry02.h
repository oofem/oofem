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

#ifndef masonry02_h
#define masonry02_h

#include "sm/Materials/ConcreteMaterials/mplasticmaterial2.h"

///@name Input fields for Masonry02
//@{
#define _IFT_Masonry02_Name "masonry02"
#define _IFT_Masonry02_ft0 "ft0"
#define _IFT_Masonry02_gfi "gfi"
#define _IFT_Masonry02_gfii "gfii"
#define _IFT_Masonry02_kn "kn"
#define _IFT_Masonry02_ks "ks"
#define _IFT_Masonry02_c0 "c0"
#define _IFT_Masonry02_tanfi0 "tanfi0"
#define _IFT_Masonry02_tanfir "tanfir"
#define _IFT_Masonry02_tanpsi "tanpsi"
#define _IFT_Masonry02_cnn "cnn"
#define _IFT_Masonry02_css "css"
#define _IFT_Masonry02_cn "cn"
#define _IFT_Masonry02_si "si"
#define _IFT_Masonry02_sp "sp"
#define _IFT_Masonry02_sm "sm"
#define _IFT_Masonry02_sr "sr"
#define _IFT_Masonry02_kp "kp"
#define _IFT_Masonry02_km "km"
#define _IFT_Masonry02_kr "kr"
#define _IFT_Masonry02_cplane "cplane"
//@}

namespace oofem {
class Domain;

/**
 * This class implements an interface masonry model based on
 * non associated multisurface plasticity.
 * Model follows the description from
 * Lourenco, P.B., Rots, J.G.: Multisurface Interface Model for Analysis of Masonry Structures
 * as published in Journal of Engineering Mechanics, vol 123, No. 7, 1997.
 */
class Masonry02 : public MPlasticMaterial2
{
protected:
    /// Tensile strength.
    double ft0 = 0.;
    /// Mode I GF.
    double gfI = 0.;
    /// Mode II GF.
    double gfII = 0.;
    /// Residual friction angle.
    double tanfir = 0.;
    /// Initial friction angle.
    double tanfi0 = 0.;
    /// Initial cohesion of joint.
    double c0 = 0.;
    /// Cap mode parameters.
    double Cnn = 0., Css = 0., Cn = 0.;
    // double fm = 0.;

    /// Elastic properties.
    double kn = 0., ks = 0.;

    /// Dilatancy angle.
    double tanpsi = 0.;

    /// Cap mode parameters.
    double sic = 0., spc = 0., smc = 0., src = 0.;
    double kp = 0., km = 0., kr = 0.;

public:
    Masonry02(int n, Domain * d);

    void initializeFrom(InputRecord &ir) override;
    bool hasMaterialModeCapability(MaterialMode mode) const override;

    const char *giveInputRecordName() const override { return _IFT_Masonry02_Name; }
    const char *giveClassName() const override { return "Masonry02"; }

    void giveStiffnessMatrix(FloatMatrix &answer,
                             MatResponseMode mode,
                             GaussPoint *gp,
                             TimeStep *tStep) const override;

    int giveSizeOfFullHardeningVarsVector() const override { return 3; }
    int giveSizeOfReducedHardeningVarsVector(GaussPoint *) const override { return 3; }

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;


protected:
    int giveMaxNumberOfActiveYieldConds(GaussPoint *gp) const override { return 2; }
    //
    // yield(YC-like functions) and loading(LC-like functions) criteria specific section
    //

    double computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector,
                               const FloatArray &stressSpaceHardeningVars) const override;

    void computeStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &stressVector,
                                     const FloatArray &stressSpaceHardeningVars) const override;
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
    int hasHardening() const override { return 1; }

    void computeReducedSSGradientMatrix(FloatMatrix &gradientMatrix,  int i, GaussPoint *gp, const FloatArray &fullStressVector,
                                        const FloatArray &strainSpaceHardeningVariables) const override;
    void computeReducedSKGradientMatrix(FloatMatrix &gradientMatrix,  int i, GaussPoint *gp, const FloatArray &fullStressVector,
                                        const FloatArray &strainSpaceHardeningVariables) const override;


    void give2dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                                GaussPoint *gp, TimeStep *tStep) const;

    void computeReducedElasticModuli(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep) const override;

    /// Cap mode related functions.
    double computeF3HardeningLaw(double k) const;
    double computeF3HardeningGradient(double k) const;
};
} // end namespace oofem
#endif // masonry02_h
