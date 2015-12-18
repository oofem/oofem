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
    double kappa;
    /// Non-equilibrated scalar measure of the largest strain level.
    double tempKappa;
    /// Damage level of material.
    double damage;
    /// Non-equilibrated damage level of material.
    double tempDamage;

    /// Characteristic length
    double charLength;
    /// Crack orientation normalized to damage magnitude. This is useful for plotting cracks as a vector field (paraview etc.).
    FloatArray crackVector;

    /// hydration-degree dependent equivalent strain at stress peak
    double var_e0;
    /// hydration-degree dependent fracture energy
    double var_gf;

#ifdef supplementary_info
    double crackWidth;
    double residTensileStrength;
#endif

public:
    MPSDamMaterialStatus(int n, Domain *d, GaussPoint *g, int nunits);
    virtual ~MPSDamMaterialStatus() { }

    virtual const FloatArray &giveViscoelasticStressVector() const { return effectiveStressVector; }
    /// Assigns tempStressVector to given vector v.
    void letTempViscoelasticStressVectorBe(FloatArray v) { tempEffectiveStressVector = std :: move(v); }
    virtual const FloatArray &giveTempViscoelasticStressVector() const { return tempEffectiveStressVector; }

    /// Returns the last equilibrated scalar measure of the largest strain level.
    double giveKappa() { return kappa; }
    /// Returns the temp. scalar measure of the largest strain level.
    double giveTempKappa() { return tempKappa; }
    /// Sets the temp scalar measure of the largest strain level to given value.
    void setTempKappa(double newKappa) { tempKappa = newKappa; }
    /// Returns the last equilibrated damage level.
    double giveDamage() { return damage; }
    /// Returns the temp. damage level.
    double giveTempDamage() { return tempDamage; }
    /// Sets the temp damage level to given value.
    void setTempDamage(double newDamage) { tempDamage = newDamage; }

    /// Returns characteristic length stored in receiver.
    double giveCharLength() { return charLength; }
    /// Sets characteristic length to given value.
    void setCharLength(double length) { charLength = length; }
    /// Returns crack vector stored in receiver. This is useful for plotting cracks as a vector field (paraview etc.).
    void giveCrackVector(FloatArray &answer);
    /// Sets crack vector to given value. This is useful for plotting cracks as a vector field (paraview etc.).
    void setCrackVector(FloatArray cv) { crackVector = cv; }

    void sete0(double e0) { var_e0 = e0; }
    void setgf(double gf) { var_gf = gf; }
    double givee0() { return var_e0; }
    double givegf() { return var_gf; }

    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

#ifdef supplementary_info
    void setCrackWidth(double src) { crackWidth = src; }
    double giveCrackWidth(void) { return crackWidth; }
    void setResidualTensileStrength(double src) { residTensileStrength = src; }
    double giveResidualTensileStrength(void) { return residTensileStrength; }
#endif


    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    // definition
    virtual const char *giveClassName() const { return "MPSDamMaterialStatus"; }

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

    bool timeDepFracturing;
    double fib_s;
    double fib_fcm28;
    bool isotropic;

    /// dummy Young's modulus
    double E;

    /// Maximum limit on omega. The purpose is elimination of a too compliant material which may cause convergence problems. Set to something like 0.99 if needed.
    double maxOmega;

    /// Equivalent strain at stress peak (or a similar parameter).
    double const_e0;

    /**
     * Determines the softening -> corresponds to the initial fracture energy. For a linear law, it is the area
     * under the stress/strain curve. For an exponential law, it is the area bounded by the elastic range
     * and a tangent to the softening part of the curve at the peak stress. For a bilinear law,
     * gf corresponds to area bounded by elasticity and the first linear softening line projected to zero stress.
     */
    double const_gf;

    /** Type characterizing the formula for the damage law. For example, linear softening can be specified
     *   with fracturing strain or crack opening.
     */
    enum SofteningType { ST_Exponential_Cohesive_Crack, ST_Linear_Cohesive_Crack, ST_Disable_Damage };

    /// Check possible snap back flag
    int checkSnapBack;

    /// Parameter specifying the type of softening (damage law).
    SofteningType softType;

    /// Method used for evaluation of characteristic element size
    ElementCharSizeMethod ecsMethod;

    /// 28-day value of tensile strength. Used only with "timedepfracturing"
    double ft28;
    /// 28-day value of fracture energy. Used only with "timedepfracturing"
    double gf28;

public:
    MPSDamMaterial(int n, Domain *d);
    virtual ~MPSDamMaterial() { }

    virtual int hasNonLinearBehaviour() { return 1; }

    virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual const char *giveInputRecordName() const { return _IFT_MPSDamMaterial_Name; }
    virtual const char *giveClassName() const { return "MPSDamMaterial"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);

    double givee0(GaussPoint *gp);
    double givegf(GaussPoint *gp);

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
     * @param[out] omega Contains result.
     * @param kappa Equivalent strain measure.
     * @param strain Total strain in full form.
     * @param gp Integration point.
     */
    void computeDamage(double &omega, double kappa, GaussPoint *gp);

    /**
     * computes the value of damage parameter omega,
     * based on a given value of equivalent strain,
     * using iterations to achieve objectivity,
     * based on the crack band concept (effective element size used)
     * @param[out] omega Contains the resulting damage.
     * @param kappa Equivalent strain measure.
     * @param gp Integration point.
     */
    void computeDamageForCohesiveCrack(double &omega, double kappa, GaussPoint *gp);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;


    virtual double computeTensileStrength(double equivalentTime);
    virtual double computeFractureEnergy(double equivalentTime);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);

    virtual void givePlaneStressStiffMtrx(FloatMatrix &answer,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);

    virtual void givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);

    virtual void give1dStressStiffMtrx(FloatMatrix &answer,
                                       MatResponseMode mode,
                                       GaussPoint *gp,
                                       TimeStep *tStep);


protected:

    //virtual double giveEModulus(GaussPoint *gp, TimeStep *tStep);
};
} // end namespace oofem
#endif // mpsdammat_h
