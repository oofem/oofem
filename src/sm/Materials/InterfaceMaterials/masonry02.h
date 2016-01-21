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

#include "Materials/ConcreteMaterials/mplasticmaterial2.h"

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
    double ft0;
    /// Mode I GF.
    double gfI;
    /// Mode II GF.
    double gfII;
    /// Residual friction angle.
    double tanfir;
    /// Initial friction angle.
    double tanfi0;
    /// Initial cohesion of joint.
    double c0;
    /// Cap mode parameters.
    double Cnn, Css, Cn;
    // double fm;

    /// Elastic properties.
    double kn, ks;

    /// Dilatancy angle.
    double tanpsi;

    /// Cap mode parameters.
    double sic, spc, smc, src;
    double kp, km, kr;
public:

    Masonry02(int n, Domain * d);
    virtual ~Masonry02();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual const char *giveInputRecordName() const { return _IFT_Masonry02_Name; }
    virtual const char *giveClassName() const { return "Masonry02"; }

    virtual void giveStiffnessMatrix(FloatMatrix &answer,
                                     MatResponseMode mode,
                                     GaussPoint *gp,
                                     TimeStep *tStep);

    virtual int giveSizeOfFullHardeningVarsVector() { return 3; }
    virtual int giveSizeOfReducedHardeningVarsVector(GaussPoint *) const { return 3; }

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;


protected:
    virtual int giveMaxNumberOfActiveYieldConds(GaussPoint *gp) { return 2; }
    //
    // yield(YC-like functions) and loading(LC-like functions) criteria specific section
    //

    virtual double computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector,
                                       const FloatArray &stressSpaceHardeningVars);

    virtual void computeStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &stressVector,
                                             const FloatArray &stressSpaceHardeningVars);
    virtual void computeStrainHardeningVarsIncrement(FloatArray &answer, GaussPoint *gp,
                                                     const FloatArray &stress, const FloatArray &dlambda,
                                                     const FloatArray &dplasticStrain, const IntArray &activeConditionMap);
    virtual void computeKGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, FloatArray &fullStressVector,
                                        const FloatArray &strainSpaceHardeningVariables);

    virtual void computeReducedHardeningVarsSigmaGradient(FloatMatrix &answer, GaussPoint *gp, const IntArray &activeConditionMap,
                                                          const FloatArray &fullStressVector,
                                                          const FloatArray &strainSpaceHardeningVars,
                                                          const FloatArray &gamma);
    virtual void computeReducedHardeningVarsLamGradient(FloatMatrix &answer, GaussPoint *gp, int actSurf,
                                                        const IntArray &activeConditionMap,
                                                        const FloatArray &fullStressVector,
                                                        const FloatArray &strainSpaceHardeningVars,
                                                        const FloatArray &gamma);
    virtual int hasHardening() { return 1; }

    virtual void computeReducedSSGradientMatrix(FloatMatrix &gradientMatrix,  int i, GaussPoint *gp, const FloatArray &fullStressVector,
                                                const FloatArray &strainSpaceHardeningVariables);
    virtual void computeReducedSKGradientMatrix(FloatMatrix &gradientMatrix,  int i, GaussPoint *gp, const FloatArray &fullStressVector,
                                                const FloatArray &strainSpaceHardeningVariables);


    void give2dInterfaceMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode,
                                                GaussPoint *gp, TimeStep *tStep);

    virtual void computeReducedElasticModuli(FloatMatrix &answer, GaussPoint *gp,
                                             TimeStep *tStep);


    /// Cap mode related functions.
    double computeF3HardeningLaw(double k);
    double computeF3HardeningGradient(double k);
};
} // end namespace oofem
#endif // masonry02_h
