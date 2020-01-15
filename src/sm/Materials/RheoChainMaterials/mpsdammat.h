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

//   ***********************************************************************************
//   ***   CLASS Damage extension of MPS model for creep and shrinakge of concrete   ***
//   ***********************************************************************************

#ifndef mpsdammat_h
#define mpsdammat_h

#include "mps.h"

#define supplementary_info

///@name Input fields for MPSDamMaterial
//@{
#define _IFT_MPSDamMaterial_Name "mpsdammat"
//prediction of tensile strength and fracture energy accroding to fib MC2010
#define _IFT_MPSDamMaterial_timedepfracturing "timedepfracturing"
#define _IFT_MPSDamMaterial_fib_s "fib_s"

#define _IFT_MPSDamMaterial_isotropic "isotropic"

#define _IFT_MPSDamMaterial_maxOmega "maxomega"

#define _IFT_MPSDamMaterial_damageLaw "damlaw"
#define _IFT_MPSDamMaterial_checkSnapBack "checksnapback"
#define _IFT_MPSDamMaterial_ft "ft"
#define _IFT_MPSDamMaterial_gf "gf"

#define _IFT_MPSDamMaterial_ft28 "ft28"
#define _IFT_MPSDamMaterial_gf28 "gf28"

//@}

namespace oofem {
#define MPSDAMMAT_ITERATION_LIMIT 1.e-9

/**
 */
class MPSDamMaterialStatus : public MPSMaterialStatus
{
protected:
    /// Equilibrated stress vector in reduced form
    FloatArray effectiveStressVector;
    /// Temporary stress vector in reduced form (increments are used mainly in nonlinear analysis)
    FloatArray tempEffectiveStressVector;
    /// Scalar measure of the largest strain level ever reached in material.
    double kappa = 0.;
    /// Non-equilibrated scalar measure of the largest strain level.
    double tempKappa = 0.;
    /// Damage level of material.
    double damage = 0.;
    /// Non-equilibrated damage level of material.
    double tempDamage = 0.;

    /// Characteristic length
    double charLength = 0.;
    /// Crack orientation normalized to damage magnitude. This is useful for plotting cracks as a vector field (paraview etc.).
    FloatArray crackVector;

    /// hydration-degree dependent equivalent strain at stress peak
    double var_e0 = 0.;
    /// hydration-degree dependent fracture energy
    double var_gf = 0.;

#ifdef supplementary_info
    double crackWidth = 0.;
    double residTensileStrength = 0.;
#endif

public:
    MPSDamMaterialStatus(GaussPoint *g, int nunits);

    const FloatArray &giveViscoelasticStressVector() const override { return effectiveStressVector; }
    /// Assigns tempStressVector to given vector v.
    void letTempViscoelasticStressVectorBe(FloatArray v) { tempEffectiveStressVector = std :: move(v); }
    virtual const FloatArray &giveTempViscoelasticStressVector() const { return tempEffectiveStressVector; }

    /// Returns the last equilibrated scalar measure of the largest strain level.
    double giveKappa() const { return kappa; }
    /// Returns the temp. scalar measure of the largest strain level.
    double giveTempKappa() const { return tempKappa; }
    /// Sets the temp scalar measure of the largest strain level to given value.
    void setTempKappa(double newKappa) { tempKappa = newKappa; }
    /// Returns the last equilibrated damage level.
    double giveDamage() const { return damage; }
    /// Returns the temp. damage level.
    double giveTempDamage() const { return tempDamage; }
    /// Sets the temp damage level to given value.
    void setTempDamage(double newDamage) { tempDamage = newDamage; }

    /// Returns characteristic length stored in receiver.
    double giveCharLength() const { return charLength; }
    /// Sets characteristic length to given value.
    void setCharLength(double length) { charLength = length; }
    /// Returns crack vector stored in receiver. This is useful for plotting cracks as a vector field (paraview etc.).
    void giveCrackVector(FloatArray &answer) const;
    /// Sets crack vector to given value. This is useful for plotting cracks as a vector field (paraview etc.).
    void setCrackVector(FloatArray cv) { crackVector = cv; }

    void sete0(double e0) { var_e0 = e0; }
    void setgf(double gf) { var_gf = gf; }
    double givee0() const { return var_e0; }
    double givegf() const { return var_gf; }

    void printOutputAt(FILE *file, TimeStep *tStep) const override;
    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

#ifdef supplementary_info
    void setCrackWidth(double src) { crackWidth = src; }
    double giveCrackWidth(void) { return crackWidth; }
    void setResidualTensileStrength(double src) { residTensileStrength = src; }
    double giveResidualTensileStrength(void) { return residTensileStrength; }
#endif

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    // definition
    const char *giveClassName() const override { return "MPSDamMaterialStatus"; }

protected:
    /**
     * Abstract service allowing to perform some initialization, when damage first appear.
     * @param kappa Scalar measure of strain level.
     * @param totalStrainVector Current total strain vector.
     * @param gp Integration point.
     */
    virtual void initDamaged(double kappa, FloatArray &totalStrainVector, GaussPoint *gp) { }
};


/**
 * This class extends the material model based on MPS theory (microprestress-solidification) for concrete
 * creep and shrinkage by a simple isotropic damage model to take into account cracking in tension.
 * Most of the methods are simply copied from "idm1" and "isodamagemodel" files.
 * Since the damage is not the main feature of this model the methods were simplified to make the code brief.
 * The following simplifications are assumed:
 * - only two softening laws are implemented: linear and exponential
 * - the equivalent strain measure uses Rankine (standard) criterion (max. principal strain)
 * - the damage parameters (tensile strength, fracture energy) can be either specified or they can be computed
 *   according to recommendations from fib Model Code 2010
 * - on contrary to idm1, the crack vector specifies the direction of crack opening and not the crack orientation
 * - snapback is checked by default
 * - element length is computed using ecsMethod = ECSM_Projection
 */
class MPSDamMaterial : public MPSMaterial
{
protected:

    bool timeDepFracturing = false;
    double fib_s = 0.;
    double fib_fcm28 = 0.;
    bool isotropic = false;

    /// dummy Young's modulus
    double E = -1.;

    /// Maximum limit on omega. The purpose is elimination of a too compliant material which may cause convergence problems. Set to something like 0.99 if needed.
    double maxOmega = 0.999999;

    /// Equivalent strain at stress peak (or a similar parameter).
    //double const_e0;

    /// constant tensile strength
    double ft = 0.;

    /**
     * Determines the softening -> corresponds to the initial fracture energy. For a linear law, it is the area
     * under the stress/strain curve. For an exponential law, it is the area bounded by the elastic range
     * and a tangent to the softening part of the curve at the peak stress. For a bilinear law,
     * gf corresponds to area bounded by elasticity and the first linear softening line projected to zero stress.
     */
    double const_gf = 0.;

    /// Check possible snap back flag
    int checkSnapBack = 1;

    /** Type characterizing the formula for the damage law. For example, linear softening can be specified
     *   with fracturing strain or crack opening.
     */
    enum SofteningType { ST_Exponential_Cohesive_Crack, ST_Linear_Cohesive_Crack, ST_Disable_Damage };
    /// Parameter specifying the type of softening (damage law).
    SofteningType softType = ST_Exponential_Cohesive_Crack;

    /// Method used for evaluation of characteristic element size
    ElementCharSizeMethod ecsMethod = ECSM_Projection;

    /// 28-day value of tensile strength. Used only with "timedepfracturing"
    double ft28 = 0.;
    /// 28-day value of fracture energy. Used only with "timedepfracturing"
    double gf28 = 0.;

public:
    MPSDamMaterial(int n, Domain *d);

    bool hasMaterialModeCapability(MaterialMode mode) const override;

    const char *giveInputRecordName() const override { return _IFT_MPSDamMaterial_Name; }
    const char *giveClassName() const override { return "MPSDamMaterial"; }

    void initializeFrom(InputRecord &ir) override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    void giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) override;

    double givee0(GaussPoint *gp) const;
    double givegf(GaussPoint *gp) const;

    /**
     * Abstract service allowing to perform some initialization, when damage first appear.
     * @param kappa Scalar measure of strain level.
     * @param totalStrainVector Current total strain vector.
     * @param gp Integration point.
     */
    void initDamaged(double kappa, FloatArray &totalStrainVector, GaussPoint *gp, TimeStep *tStep);

    void initDamagedFib(GaussPoint *gp, TimeStep *tStep);

    //void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);

    /**
     * Computes the value of damage omega, based on given value of equivalent strain.
     * @param kappa Equivalent strain measure.
     * @param strain Total strain in full form.
     * @param gp Integration point.
     * @return Damage.
     */
    double computeDamage(double kappa, GaussPoint *gp) const;

    /**
     * computes the value of damage parameter omega,
     * based on a given value of equivalent strain,
     * using iterations to achieve objectivity,
     * based on the crack band concept (effective element size used)
     * @param kappa Equivalent strain measure.
     * @param gp Integration point.
     * @return Damage.
     */
    double computeDamageForCohesiveCrack(double kappa, GaussPoint *gp) const;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;


    virtual double computeTensileStrength(double equivalentTime) const;
    virtual double computeFractureEnergy(double equivalentTime) const;

    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp,
                                                    TimeStep *tStep) const override;

    FloatMatrixF<3,3> givePlaneStressStiffMtrx(MatResponseMode mode, GaussPoint *gp,
                                               TimeStep *tStep) const override;

    FloatMatrixF<4,4> givePlaneStrainStiffMtrx(MatResponseMode mode, GaussPoint *gp,
                                               TimeStep *tStep) const override;

    FloatMatrixF<1,1> give1dStressStiffMtrx(MatResponseMode mode, GaussPoint *gp,
                                            TimeStep *tStep) const override;
};
} // end namespace oofem
#endif // mpsdammat_h
