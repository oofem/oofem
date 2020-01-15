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

#ifndef frcfcm_h
#define frcfcm_h


#include "concretefcm.h"

///@name Input fields for Fcm_
//@{
#define _IFT_FRCFCM_Name "frcfcm"
#define _IFT_FRCFCM_Vf "vf"
#define _IFT_FRCFCM_Lf "lf"
#define _IFT_FRCFCM_Df "df"
#define _IFT_FRCFCM_Ef "ef"
#define _IFT_FRCFCM_nuf "nuf"
#define _IFT_FRCFCM_Gfib "gfib"
#define _IFT_FRCFCM_kfib "kfib"
#define _IFT_FRCFCM_tau_0 "tau_0"
#define _IFT_FRCFCM_b0 "b0"
#define _IFT_FRCFCM_b1 "b1"
#define _IFT_FRCFCM_b2 "b2"
#define _IFT_FRCFCM_b3 "b3"
#define _IFT_FRCFCM_f "f"
#define _IFT_FRCFCM_M "m"
#define _IFT_FRCFCM_orientationVector "fibreorientationvector"
#define _IFT_FRCFCM_fssType "fsstype"
#define _IFT_FRCFCM_fDamType "fdamtype"
#define _IFT_FRCFCM_fiberType "fibertype"
#define _IFT_FRCFCM_gammaCrack "gammacrack"
#define _IFT_FRCFCM_computeCrackSpacing "computecrackspacing"
#define _IFT_FRCFCM_fibreActivationOpening "fibreactivationopening"
#define _IFT_FRCFCM_dw0 "dw0"
#define _IFT_FRCFCM_dw1 "dw1"

//@}

namespace oofem {
/**
 * This class manages the status of FRCFCM
 */
class FRCFCMStatus : public ConcreteFCMStatus
{
protected:
    /// Damage level of material.
    double damage = 0.;
    /// Non-equilibrated damage level of material.
    double tempDamage = 0.;

public:
    FRCFCMStatus(GaussPoint *g);

    /// Returns the last equilibrated damage level.
    double giveDamage() const { return damage; }
    /// Returns the temporary damage level.
    double giveTempDamage() const { return tempDamage; }
    /// Sets the temp damage level to given value.
    void setTempDamage(double newDamage) { tempDamage = newDamage; }

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    const char *giveClassName() const override { return "FRCFCMStatus"; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
};


/**
 * This class implements a FRCFCM material (Fiber Reinforced Concrete base on Fixed Crack Model)
 * in a finite element problem. This class provides an extension to the ConcreteFCM which
 * serves as a material model for matrix while the present class FRCFCM adds the contribution
 * of fibers. The contribution of the two constituents are defined by Vf parameter = volume of fibers.
 */
class FRCFCM : public ConcreteFCM
{
public:
    FRCFCM(int n, Domain *d);

    void initializeFrom(InputRecord &ir) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new FRCFCMStatus(gp); }

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

protected:
    /// fiber shear strength at zero slip
    double tau_0 = 0.;

    /// micromechanical parameter for fiber shear according to Sajdlova
    double b0 = 0.;
    /// micromechanical parameter for fiber shear according to Kabele
    double b1 = 0., b2 = 0., b3 = 0.;

    /// snubbing factor "f"
    double f = 0.;

    /// auxiliary parameter computed from snubbing factor "f"
    double g = 0.;

    /// volume fraction of fibers
    double Vf = 0.;

    /// fiber length
    double Lf = 0.;

    /// fiber diameter
    double Df = 0.;

    /// fiber Young's modulus
    double Ef = 0.;

    /// fiber shear modulus
    double Gfib = 0.;

    /// fiber cross-sectional shear factor
    double kfib = 0.;

    /// transitional opening
    double w_star = 0.;

    /// aux. factor
    double eta = 0.;

    /// shear strain at full fibers rupture
    double gammaCrackFail = 0.;

    /// minimum opening at which damage can start
    double minDamageOpening = 0.;

    /**
     * Exponent in the unloading-reloading constitutive law.
     * the function is defined by two points - the origin
     * and the stress at maximum cracking strain
     * sigma(w) = sig_max * w^M / w_max^M
     */
    int M = 0;

    /// orientation of fibres
    FloatArray orientationVector;

    /// crack opening at which the crossing fibers begin to be activated
    double fibreActivationOpening = 0.;

    /**
     * smooth transition of the bridging stress if fibreActivationOpening is applied
     * dw0 = distance from the fibreActivationOpening where the smooth transition starts
     * dw1 = distance from the fibreActivationOpening where the smooth transition ends
     * smoothen = flag
     */
    double dw0 = 0., dw1 = 0.;
    bool smoothen = false;

    /**
     * Type strength of the shear bond. This is activated only for short fibers (SAF or SRF)
     * once the maximum crack opening exceeds w*
     */
    enum FiberShearStrengthType { FSS_NONE, FSS_Sajdlova, FSS_Kabele, FSS_Havlasek, FSS_Unknown };
    FiberShearStrengthType fiberShearStrengthType = FSS_Unknown;

    /// Type of fibre damage which is triggered by crack shearing strain = w / u
    enum FiberDamageType { FDAM_NONE, FDAM_GammaCrackLin, FDAM_GammaCrackExp, FDAM_Unknown };
    FiberDamageType fiberDamageType = FDAM_Unknown;

    /**
     * Type fo fibers in the composite.
     * CAF = continuous aligned fibers
     * SAF = short aligned fibers
     * SRF = short random fibers
     * SRF2D = short fibers random in 2D
     */
    enum FiberType { FT_CAF, FT_SAF, FT_SRF, FT_SRF2D, FT_Unknown };
    FiberType fiberType = FT_Unknown;

    double giveCrackingModulus(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep, int i) override;
    /// evaluates the fiber bond if w > w*
    virtual double computeFiberBond(double w);
    double giveNormalCrackingStress(GaussPoint *gp, TimeStep *tStep, double eps_cr, int i) override;
    /// compute the nominal stress in fibers in the i-th crack
    virtual double computeStressInFibersInCracked(GaussPoint *gp, TimeStep *tStep, double eps_cr, int i);
    double computeEffectiveShearModulus(GaussPoint *gp, TimeStep *tStep, int i) override;
    double computeD2ModulusForCrack(GaussPoint *gp, TimeStep *tStep, int icrack) override;
    /// estimate shear modulus for a given crack plane (1, 2, 3). Uses equilibrated value of damage.
    virtual double estimateD2ModulusForCrack(GaussPoint *gp, TimeStep *tStep, int icrack);
    double maxShearStress(GaussPoint *gp, TimeStep *tStep, int i) override;
    /// evaluates temporary value of damage caused by fibre shearing
    virtual double computeTempDamage(GaussPoint *gp, TimeStep *tStep);
    /// computes crack spacing based on composition of the fibre composite
    virtual double computeCrackSpacing();
    /// compute the angle between the fibre and i-th crack normal
    virtual double computeCrackFibreAngle(GaussPoint *gp, int i);
    void checkSnapBack(GaussPoint *gp, TimeStep *tStep, int crack) override { }
    bool isStrengthExceeded(const FloatMatrix &base, GaussPoint *gp, TimeStep *tStep, int iCrack, double trialStress) override;
    double computeShearStiffnessRedistributionFactor(GaussPoint *gp, TimeStep *tStep, int ithCrackPlane, int jthCrackDirection) override;
    double computeOverallElasticStiffness(GaussPoint *gp, TimeStep *tStep) override;
    double computeOverallElasticShearModulus(GaussPoint *gp, TimeStep *tStep) override { return this->computeOverallElasticStiffness(gp, tStep) / ( 2. * ( 1. + linearElasticMaterial.givePoissonsRatio() ) ); }
};
} // end namespace oofem
#endif // frcfcm_h
