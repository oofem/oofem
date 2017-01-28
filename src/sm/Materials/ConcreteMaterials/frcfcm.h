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
    double damage;
    /// Non-equilibrated damage level of material.
    double tempDamage;

public:
    FRCFCMStatus(int n, Domain *d, GaussPoint *g);
    virtual ~FRCFCMStatus();

    /// Returns the last equilibrated damage level.
    double giveDamage() { return damage; }
    /// Returns the temporary damage level.
    double giveTempDamage() { return tempDamage; }
    /// Sets the temp damage level to given value.
    void setTempDamage(double newDamage) { tempDamage = newDamage; }
    /// Writes information into the output file.
    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    // definition
    virtual const char *giveClassName() const { return "FRCFCMStatus"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    // saves current context(state) into stream
    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
};


/**
 * This class implements a FRCFCM material (Fiber Reinforced Concrete base on Fixed Crack Model)
 * in a finite element problem. This class provides an extension to the ConcreteFCM which
 * serves as a material model for matrix while the present class FRCFCM adds the contribution
 * of fibers
 */
class FRCFCM : public ConcreteFCM
{
public:
    FRCFCM(int n, Domain *d);
    virtual ~FRCFCM() {}

    // identification and auxiliary functions
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new FRCFCMStatus(1, domain, gp); }

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

protected:

    /// fiber shear strength at zero slip
    double tau_0;

    /// micromechanical parameter for fiber shear according to Sajdlova
    double b0;
    /// micromechanical parameter for fiber shear according to Kabele
    double b1, b2, b3;

    /// snubbing factor "f"
    double f;

    /// auxiliary parameter computed from snubbing factor "f"
    double g;

    /// volume fraction of fibers
    double Vf;

    /// fiber length
    double Lf;

    /// fiber diameter
    double Df;

    /// fiber Young's modulus
    double Ef;

    /// fiber shear modulus
    double Gfib;

    /// fiber cross-sectional shear factor
    double kfib;

    /// transitional opening
    double w_star;

    /// aux. factor
    double eta;

    /// shear strain at full fibers rupture
    double gammaCrackFail;

    /// minimum opening at which damage can start
    double minDamageOpening;

    /**
     * Exponent in the unloading-reloading constitutive law.
     * the function is defined by two points - the origin
     * and the stress at maximum cracking strain
     * sigma(w) = sig_max * w^M / w_max^M
     */
    int M;

    /// orientation of fibres
    FloatArray orientationVector;

    /// crack opening at which the crossing fibers begin to be activated
    double fibreActivationOpening;

    /**
     * smooth transition of the bridging stress if fibreActivationOpening is applied
     * dw0 = distance from the fibreActivationOpening where the smooth transition starts
     * dw1 = distance from the fibreActivationOpening where the smooth transition ends
     * smoothen = flag
     */
    double dw0, dw1;
    bool smoothen;

    /**
     * Type strength of the shear bond. This is activated only for short fibers (SAF or SRF)
     * once the maximum crack opening exceeds w*
     */
    enum FiberShearStrengthType { FSS_NONE, FSS_Sajdlova, FSS_Kabele, FSS_Havlasek, FSS_Unknown };
    FiberShearStrengthType fiberShearStrengthType;

    /// Type of fibre damage which is triggered by crack shearing strain = w / u
    enum FiberDamageType { FDAM_NONE, FDAM_GammaCrackLin, FDAM_GammaCrackExp, FDAM_Unknown };
    FiberDamageType fiberDamageType;

    /**
     * Type fo fibers in the composite.
     * CAF = continuous aligned fibers
     * SAF = short aligned fibers
     * SRF = short random fibers
     */
    enum FiberType { FT_CAF, FT_SAF, FT_SRF, FT_Unknown };
    FiberType fiberType;

    /// returns stiffness in the normal direction of the i-th crack
    virtual double giveCrackingModulus(MatResponseMode rMode, GaussPoint *gp, int i);
    /// evaluates the fiber bond if w > w*
    virtual double computeFiberBond(double w);
    /// computes normal stress associated with i-th crack direction
    virtual double giveNormalCrackingStress(GaussPoint *gp, double eps_cr, int i);
    /// compute the nominal stress in fibers in the i-th crack
    virtual double computeStressInFibersInCracked(GaussPoint *gp, double eps_cr, int i);
    /// returns Geff which is necessary in the global stiffness matrix
    virtual double computeEffectiveShearModulus(GaussPoint *gp, int i);
    /// shear modulus for a given crack plane (1, 2, 3)
    virtual double computeD2ModulusForCrack(GaussPoint *gp, int icrack);
    /// estimate shear modulus for a given crack plane (1, 2, 3). Uses equilibrated value of damage.
    virtual double estimateD2ModulusForCrack(GaussPoint *gp, int icrack);
    /// computes the maximum value of the shear stress; if the shear stress exceeds this value, it is cropped
    virtual double maxShearStress(GaussPoint *gp, int i);
    /// evaluates temporary value of damage caused by fibre shearing
    virtual double computeTempDamage(GaussPoint *gp);
    /// computes crack spacing based on composition of the fibre composite
    virtual double computeCrackSpacing(void);
    /// compute the angle between the fibre and i-th crack normal
    virtual double computeCrackFibreAngle(GaussPoint *gp, int i);
    /// overrides real checking from concretefcm.C, here we assume that the fibers provide sufficient strength to prevent snapback
    virtual void checkSnapBack(GaussPoint *gp, int crack) { return; }
    /// the method from fcm is overridden to consider stress split between the matrix and fibers
    virtual bool isStrengthExceeded(const FloatMatrix &base, GaussPoint *gp, TimeStep *tStep, int iCrack, double trialStress);
    /// function calculating ratio used to split shear slips on two crack planes
    virtual double computeShearStiffnessRedistributionFactor(GaussPoint *gp, int ithCrackPlane, int jthCrackDirection);
    /// according to the volume fraction of fibers and the Young's moduli estimates the overall stiffness of the composite
    virtual double computeOverallElasticStiffness(void);
    /// from the Poisson's ratio of matrix and the overall stiffness estimates G
    virtual double computeOverallElasticShearModulus(void) { return this->computeOverallElasticStiffness() / ( 2. * ( 1. + this->giveLinearElasticMaterial()->givePoissonsRatio() ) ); }
};
} // end namespace oofem
#endif // frcfcm_h
