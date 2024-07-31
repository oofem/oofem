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

#ifndef ConcreteDPM2_h
#define ConcreteDPM2_h

#include "sm/Materials/structuralmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "cltypes.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "mathfem.h"

#define CDPM2_TOL 1.e-6
#define keep_track_of_dissipated_energy
///@name Input fields for ConcreteDPM2
//@{
#define _IFT_ConcreteDPM2_Name "con2dpm"
#define _IFT_ConcreteDPM2_fc "fc"
#define _IFT_ConcreteDPM2_ft "ft"
#define _IFT_ConcreteDPM2_ecc "ecc"
#define _IFT_ConcreteDPM2_kinit "kinit"
#define _IFT_ConcreteDPM2_ahard "ahard"
#define _IFT_ConcreteDPM2_bhard "bhard"
#define _IFT_ConcreteDPM2_chard "chard"
#define _IFT_ConcreteDPM2_dhard "dhard"
#define _IFT_ConcreteDPM2_dilation "dilation"
#define _IFT_ConcreteDPM2_asoft "asoft"
#define _IFT_ConcreteDPM2_hp "hp"
#define _IFT_ConcreteDPM2_yieldtol "yieldtol"
#define _IFT_ConcreteDPM2_newtoniter "newtoniter"
#define _IFT_ConcreteDPM2_wf "wf"
#define _IFT_ConcreteDPM2_efc "efc"
#define _IFT_ConcreteDPM2_softeningType "stype"
#define _IFT_ConcreteDPM2_ftOne "ft1"
#define _IFT_ConcreteDPM2_wfOne "wf1"
#define _IFT_ConcreteDPM2_strengthratetype "sratetype"
#define _IFT_ConcreteDPM2_energyratetype "eratetype"
#define _IFT_ConcreteDPM2_deltatime "deltat"
#define _IFT_ConcreteDPM2_helem "helem"
#define _IFT_ConcreteDPM2_damflag "damflag"
//@}

namespace oofem {
/**
 * This class implements the material status associated to ConcreteDPM2.
 * Main article is "CDPM2: A damage-plasticity approach to modelling the failure of concrete"
 * International Journal of Solids and Structures, Volume 50, Issue 24, November 2013, Pages 3805-3816
 * See also https://petergrassl.com/Research/Constitutive/index.html for more info.
 * @author Peter Grassl, Dimitrios Xenos
 */
class ConcreteDPM2Status : public StructuralMaterialStatus
{
public:
    /// Values of history variable state_flag.
    enum state_flag_values {
        ConcreteDPM2_Elastic,
        ConcreteDPM2_Unloading,
        ConcreteDPM2_Plastic,
        ConcreteDPM2_Damage,
        ConcreteDPM2_PlasticDamage,
        ConcreteDPM2_VertexCompression,
        ConcreteDPM2_VertexTension,
        ConcreteDPM2_VertexCompressionDamage,
        ConcreteDPM2_VertexTensionDamage
    };



protected:

    /// @name History variables of the plasticity model
    //@{
    FloatArrayF< 6 >plasticStrain;
    FloatArrayF< 6 >tempPlasticStrain;

    FloatArrayF< 6 >reducedStrain;
    FloatArrayF< 6 >tempReducedStrain;

    FloatArrayF< 6 >effectiveStress;
    FloatArrayF< 6 >tempEffectiveStress;
    //@}

    /// @name Hardening variable
    //@{
    double kappaP = 0.;
    double tempKappaP = 0.;
    //@}

    double kappaPPeak = 0.;

    double deltaLambda = 0.;

    double le = 0.;

    double alpha = 0.;
    double tempAlpha = 0.;

    double equivStrain = 0.;
    double tempEquivStrain = 0.;

    double equivStrainTension = 0.;
    double tempEquivStrainTension = 0.;

    double equivStrainCompression = 0.;
    double tempEquivStrainCompression = 0.;

    double kappaDTension = 0.;
    double tempKappaDTension = 0.;

    double kappaDCompression = 0.;
    double tempKappaDCompression = 0.;

    double kappaDTensionOne = 0.;
    double tempKappaDTensionOne = 0.;

    double kappaDCompressionOne = 0.;
    double tempKappaDCompressionOne = 0.;

    double kappaDTensionTwo = 0.;
    double tempKappaDTensionTwo = 0.;

    double kappaDCompressionTwo = 0.;
    double tempKappaDCompressionTwo = 0.;

    double damageTension = 0.;
    double tempDamageTension = 0.;

    double damageCompression = 0.;
    double tempDamageCompression = 0.;

    double deltaEquivStrain = 0.;

    double rateFactor = 1.;
    double tempRateFactor = 0.;

    /// Strains that are used for calculation of strain rates
    double rateStrain = 0.;
    double tempRateStrain = 0.;

    /// Indicates the state (i.e. elastic, unloading, plastic, damage, vertex) of the Gauss point
    int state_flag = ConcreteDPM2Status::ConcreteDPM2_Elastic;
    int temp_state_flag = ConcreteDPM2Status::ConcreteDPM2_Elastic;


#ifdef keep_track_of_dissipated_energy
    /// Density of total work done by stresses on strain increments.
    double stressWork = 0.;
    /// Non-equilibrated density of total work done by stresses on strain increments.
    double tempStressWork = 0.;
    /// Density of dissipated work.
    double dissWork = 0.;
    /// Non-equilibrated density of dissipated work.
    double tempDissWork = 0.;
#endif

public:
    /// Constructor
    ConcreteDPM2Status(GaussPoint *gp);

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;
    void printOutputAt(FILE *file, TimeStep *tStep) const override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;
    const char *giveClassName() const override { return "ConcreteDPM2Status"; }


    // Inline functions for access to state variables

    /**
     * Get the reduced strain vector from the material status.
     * @return Strain vector.
     */
    const FloatArrayF< 6 > &giveReducedStrain() const { return reducedStrain; }

    /**
     * Get the reduced strain vector from the material status.
     * @return Strain vector.
     */
    const FloatArrayF< 6 > &giveTempReducedStrain() const { return tempReducedStrain; }


    /**
     * Get the effective stress vector from the material status.
     * @return effective stress vector.
     */
    const FloatArrayF< 6 > &giveTempEffectiveStress() const { return tempEffectiveStress; }

    /**
     * Get the plastic strain vector from the material status.
     * @return Strain vector.
     */
    const FloatArrayF< 6 > &givePlasticStrain() const { return plasticStrain; }

    /**
     * Get the deviatoric plastic strain norm from the material status.
     * @return DeviatoricPlasticStrainNorm.
     */
    double giveDeviatoricPlasticStrainNorm() const
    {
        auto dev = StructuralMaterial::computeDeviator(plasticStrain);
        return sqrt( .5 * ( 2. * dev [ 0 ] * dev [ 0 ] + 2. * dev [ 1 ] * dev [ 1 ] + 2. * dev [ 2 ] * dev [ 2 ] +
                            dev [ 3 ] * dev [ 3 ] + dev [ 4 ] * dev [ 4 ] + dev [ 5 ] * dev [ 5 ] ) );
    }

    /**
     * Get the volumetric plastic strain from the material status.
     * @return volumetricPlasticStrainNorm.
     */
    double giveVolumetricPlasticStrain() const
    {
        return 1. / 3. * ( plasticStrain [ 0 ] + plasticStrain [ 1 ] + plasticStrain [ 2 ] );
    }

    /**
     * Get the hardening variable of the plasticity model.
     * @return The hardening variable of the plasticity model.
     */
    double giveKappaP() const
    { return kappaP; }

    /**
     * Get the hardening variable of the damage model from the
     * material status.
     * @return Hardening variable kappaD.
     */
    double giveKappaDTensionOne() const
    { return kappaDTensionOne; }

    /**
     * Get the compression hardening variable one of the damage model from the
     * material status.
     * @return Hardening variable kappaDCompressionOne.
     */
    double giveKappaDCompressionOne() const
    { return kappaDCompressionOne; }


    /**
     * Get the tension hardening variable two of the damage model from the
     * material status.
     * @return Hardening variable kappaDTensionTwo.
     */
    double giveKappaDTensionTwo() const
    { return kappaDTensionTwo; }


    /**
     * Get the compression hardening variable two of the damage model from the
     * material status.
     * @return Hardening variable kappaDCompressionTwo.
     */
    double giveKappaDCompressionTwo() const
    { return kappaDCompressionTwo; }


    /**
     * Get the equivalent strain from the
     * material status.
     * @return Equivalent strain equivStrain.
     */
    double giveEquivStrain() const
    { return equivStrain; }

    /**
     * Get the tension equivalent strain from the
     * material status.
     * @return Equivalent strain equivStrainTension.
     */
    double giveEquivStrainTension() const
    { return equivStrainTension; }


    /**
     * Get the compression equivalent strain from the
     * material status.
     * @return Equivalent strain equivStrainCompression.
     */
    double giveEquivStrainCompression() const
    { return equivStrainCompression; }

    /**
     * Get the tension damage variable of the damage model from the
     * material status.
     * @return Tension damage variable damageTension.
     */
    double giveDamageTension() const
    { return damageTension; }

    /**
     * Get the compressive damage variable of the damage model from the
     * material status.
     * @return Compressive damage variable damageCompression.
     */
    double giveDamageCompression() const
    { return damageCompression; }

    /**
     * Get the rate factor of the damage model from the
     * material status.
     * @return rate factor rateFactor.
     */
    double giveRateFactor() const
    { return rateFactor; }

    /**
     * Get the temp variable of the damage model from the
     * material status.
     * @return Damage variable damage.
     */
    double giveTempRateFactor() const
    { return tempRateFactor; }


    double giveRateStrain() const
    { return rateStrain; }

    void letTempRateStrainBe(double v)
    { tempRateStrain = v; }


    void letTempAlphaBe(double v)
    { tempAlpha = v; }

    /**
     * Get the state flag from the material status.
     * @return State flag (i.e. elastic, unloading, yielding, vertex case yielding)
     */
    int giveStateFlag() const
    { return state_flag; }


    // giveTemp:

    // Functions used to access the temp variables.
    /**
     * Get the temp value of the full plastic strain vector from the material status.
     * @return Temp value of plastic strain vector.
     */
    const FloatArrayF< 6 > &giveTempPlasticStrain() const { return tempPlasticStrain; }

    /**
     *  Get the temp value of the volumetric plastic strain in plane stress
     */
    double giveTempVolumetricPlasticStrain() const
    { return 1. / 3. * ( tempPlasticStrain [ 0 ] + tempPlasticStrain [ 1 ] + tempPlasticStrain [ 2 ] ); }

    /**
     * Get the temp value of the hardening variable of the plasticity model
     * from the material status.
     * @return Temp value of hardening variable kappaP.
     */
    double giveTempKappaP() const
    { return tempKappaP; }

    double giveDeltaLambda() const
    { return deltaLambda; }

    /**
     * Get the temp value of the hardening variable of the damage model
     * from the material status.
     * @return Temp value of the damage variable damage.
     */
    double giveKappaDTension() const
    { return kappaDTension; }

    /**
     * Get value of alpha from the status.
     * @return value of alpha.
     */
    double giveAlpha() const
    { return alpha; }

    /**
     * Get value of temp alpha from the status.
     * @return value of temp alpha.
     */
    double giveTempAlpha() const
    { return tempAlpha; }



    /**
     * Get the temp value of the hardening variable of the damage model
     * from the material status.
     * @return Temp value of the damage variable damage.
     */
    double giveKappaDCompression() const
    { return kappaDCompression; }

    /**
     * Get the temp value of the hardening variable of the damage model
     * from the material status.
     * @return Temp value of the damage variable damage.
     */
    double giveTempDamageTension() const
    { return tempDamageTension; }

    /**
     * Get the temp value of the hardening variable of the damage model
     * from the material status.
     * @return Temp value of the damage variable damage.
     */
    double giveTempDamageCompression() const
    { return tempDamageCompression; }

    /**
     * Get the temp value of the hardening variable of the damage model
     * from the material status.
     * @return Temp value of the damage variable damage.
     */
    double giveDeltaEquivStrain() const
    { return deltaEquivStrain; }

    /**
     * Get the temp value of the state flag from the material status.
     * @return Temp value of the state flag (i.e. elastic, unloading,
     * yielding, vertex case yielding).
     */
    int giveTempStateFlag() const
    { return temp_state_flag; }

    // letTemp...be :
    // Functions used by the material to assign a new value to a temp variable.
    /**
     * Assign the temp value of deviatoric plastic strain.
     * @param v New temp value of deviatoric plastic strain
     */
    void letTempPlasticStrainBe(const FloatArrayF< 6 > &v)
    { tempPlasticStrain = v; }


    void letTempReducedStrainBe(const FloatArrayF< 6 > &v)
    { tempReducedStrain = v; }

    void letTempEffectiveStressBe(const FloatArrayF< 6 > &v)
    { tempEffectiveStress = v; }


    /**
     * Assign the temp value of the hardening variable of the plasticity model.
     * @param v New temp value of the hardening variable
     */
    void letTempKappaPBe(double v)
    { tempKappaP = v; }

    void letDeltaLambdaBe(double v)
    { deltaLambda = v; }

    /**
     * Assign the temp value of the rate factor of the damage model.
     * @param v New temp value of the damage variable
     */
    void letTempKappaDTensionBe(double v)
    { tempKappaDTension = v; }

    /**
     * Assign the temp value of the rate factor of the damage model.
     * @param v New temp value of the damage variable
     */
    void letTempKappaDCompressionBe(double v)
    { tempKappaDCompression = v; }

    /**
     * Assign the temp value of the hardening variable of the damage model.
     * @param v New temp value of the hardening variable
     */
    void letTempKappaDTensionOneBe(double v)
    { tempKappaDTensionOne = v; }

    /**
     * Assign the temp value of the hardening variable of the damage model.
     * @param v New temp value of the hardening variable
     */
    void letTempKappaDCompressionOneBe(double v)
    { tempKappaDCompressionOne = v; }

    /**
     * Assign the temp value of the second tension hardening variable of the damage model.
     * @param v New temp value of the second tension hardening variable
     */
    void letTempKappaDTensionTwoBe(double v)
    { tempKappaDTensionTwo = v; }

    /**
     * Assign the temp value of the second compression hardening variable of the damage model.
     * @param v New temp value of the second compression hardening variable
     */
    void letTempKappaDCompressionTwoBe(double v)
    { tempKappaDCompressionTwo = v; }

    /**
     * Assign the temp value of the tensile damage variable of the damage model.
     * @param v New temp value of the tensile damage variable
     */
    void letTempDamageTensionBe(double v)
    { tempDamageTension = v; }

    /**
     * Assign the temp value of the compressive damage variable of the damage model.
     * @param v New temp value of the compressive damage variable
     */
    void letTempDamageCompressionBe(double v)
    { tempDamageCompression = v; }

    /**
     * Assign the temp value of the rate factor of the damage model.
     * @param v New temp value of the damage variable
     */
    void letTempRateFactorBe(double v)
    { tempRateFactor = v; }

    /**
     * Assign the temp value of the rate factor of the damage model.
     * @param v New temp value of the damage variable
     */
    void letTempEquivStrainBe(double v)
    { tempEquivStrain = v; }

    /**
     * Assign the temp value of the rate factor of the damage model.
     * @param v New temp value of the damage variable
     */
    void letTempEquivStrainTensionBe(double v)
    { tempEquivStrainTension = v; }

    /**
     * Assign the temp value of the rate factor of the damage model.
     * @param v New temp value of the damage variable
     */
    void letTempEquivStrainCompressionBe(double v)
    { tempEquivStrainCompression = v; }

    /**
     *  Gives the characteristic length.
     */
    double giveLe() const { return le; }

    /**
     *  Sets the characteristic length.
     */
    void setLe(double ls)
    { le = ls; }
    /**
     * Assign the temp value of the state flag.
     * @param v New temp value of the state flag (i.e. elastic, unloading, yielding,
     * vertex case yielding).
     */
    void letTempStateFlagBe(const int v)
    { temp_state_flag = v; }

    void letKappaPPeakBe(double kappa)
    { kappaPPeak = kappa; }
#ifdef keep_track_of_dissipated_energy
    /// Returns the density of total work of stress on strain increments.
    double giveStressWork() { return stressWork; }
    /// Returns the temp density of total work of stress on strain increments.
    double giveTempStressWork() { return tempStressWork; }
    /// Sets the density of total work of stress on strain increments to given value.
    void setTempStressWork(double w) { tempStressWork = w; }
    /// Returns the density of dissipated work.
    double giveDissWork() { return dissWork; }
    /// Returns the density of temp dissipated work.
    double giveTempDissWork() { return tempDissWork; }
    /// Sets the density of dissipated work to given value.
    void setTempDissWork(double w) { tempDissWork = w; }
    /**
     * Computes the increment of total stress work and of dissipated work
     * (gf is the dissipation density per unit volume at complete failure,
     * it is needed only to determine which extremely small dissipation
     * can be set to zero to get clean results, but parameter gf can be
     * set to zero if not available).
     */
    void computeWork(GaussPoint *gp, double ft);
#endif
};


//   **************************************************
//   *** CLASS CONCRETE DAMAGE PLASTICITY MODEL 2   ***
//   **************************************************

/**
 * This class contains the combination of a local plasticity model for concrete with a local isotropic damage model.
 * This is an extension of concretedpm2. The yield surface of the plasticity model is based on the extension of the Menetrey and Willam yield criterion.
 * The flow rule is nonassociated. The evolution laws of the hardening variables depend on the stress state.
 * The plasticity model describes only hardening and perfect plasticity.
 * It is based on the effective stress. The damage parameter of the isotropic damage model is based on the total volumetric strain.
 * An exponential softening law is implemented.
 *
 * @author Peter Grassl, Dimitrios Xenos
 */
class ConcreteDPM2 : public StructuralMaterial
{
public:

    enum ConcreteDPM2_ReturnResult {
        RR_Unknown,
        RR_NotConverged,
        RR_Converged
    };

    enum ConcreteDPM2_ReturnType {
        RT_Unknown,
        RT_Regular,
        RT_Tension,
        RT_Compression,
        RT_Auxiliary
    };


protected:
    /// Parameters of the yield surface of the plasticity model. fc is the uniaxial compressive strength, ft the uniaxial tensile strength and ecc controls the out of roundness of the deviatoric section.
    double fc = 0., ft = 0., ecc = 0.;

    /** Parameter which controls the type of damage that is used together with plasticity
     * 0: tensile and compression damage based on split of stress state
     */
    int damageFlag = 0;

    double e0 = 0.;

    /// Parameter of the ductilityMeasure of the plasticity model.
    double AHard = 0.;
    /// Parameter of the ductilityMeasure of the plasticity model.
    double BHard = 0.;
    /// Parameter of the ductilityMeasure of the plasticity model.
    double CHard = 0.;
    /// Parameter of the ductilityMeasure of the plasticity model.
    double DHard = 0.;

    /// Hardening modulus.
    double hardeningModulus = 0.;

    /// Parameter of the ductilityMeasure of the damage model.
    double ASoft = 0.;

    /// Parameter of the hardening law of the plasticity model.
    double yieldHardPrimePeak = 0.;

    /// Parameter of the hardening law of the plasticity model.
    double yieldHardInitial = 0.;

    /// Control parameter for te volumetric plastic flow of the plastic potential.
    double dilationConst = 0.;

    /// Friction parameter of the yield surface.
    double m = 0.;

    /// Dilation parameter of the plastic potential.
    double mQ = 0.;

    /// Element size (to be used in fracture energy approach (crack band).
    double helem = 0.;

    /// Pointer for linear elastic material.
    IsotropicLinearElasticMaterial linearElasticMaterial;

    /// Elastic Young's modulus.
    double eM = 0.;
    /// Elastic shear modulus.
    double gM = 0.;
    /// Elastic bulk modulus.
    double kM = 0.;
    /// Elastic poisson's ration.
    double nu = 0.;

    /// Control parameter for the exponential softening law.
    double efCompression = 0.;

    /// Control parameter for the linear/bilinear softening law in tension.
    double wf = 0.;

    /// Control parameter for the bilinear softening law in tension.
    double wfOne = 0.;

    /// Control parameter for the bilinear softening law.
    double ftOne = 0.;

    /// yield tolerance for the plasticity model.
    double yieldTol = 0.;

    /// yield tolerance for the damage model.
    double yieldTolDamage = 0.;

    /// Maximum number of iterations for stress return.
    int newtonIter = 0;

    /// Type of softening function used.
    int softeningType = 0;

    /// Input parameter which simulates a loading rate. Only for debugging purposes.
    double deltaTime = 0.;

    /** Type of strength strain rate dependence used.
     * 0 = no strain rate (default)
     * 1 = Model Code 2010 initial branch of strain rate effect for strength
     * 2 = Model Code 2010 initial and second branch of strain rate effect for strength
     */
    int strengthRateType = 0;

    /** Type of energy strain rate dependence used if strengthRateType >0 .
     * 0 = mod. CEB strain rate effect for strength with constant fracture energy
     * 1 = CEB strain rate effect for strength and linear for fracture energy
     * 2 = mod. CEB strain rate effect for strength and squared for fracture energy
     *
     */
    int energyRateType = 0;

public:
    /// Constructor
    ConcreteDPM2(int n, Domain *d);

    void initializeFrom(InputRecord &ir) override;

    const char *giveClassName() const override { return "ConcreteDPM2"; }
    const char *giveInputRecordName() const override { return _IFT_ConcreteDPM2_Name; }

    FloatArrayF< 6 >giveRealStressVector_3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) const override;

    bool hasMaterialModeCapability(MaterialMode mode) const override;

    /**
     * Perform stress return of the plasticity model and compute history variables.
     * @param gp Gauss point.
     * @param D stiffness matrix
     * @param strain Strain vector of this Gauss point.
     */
    FloatArrayF< 6 >performPlasticityReturn(GaussPoint *gp,
                                            const FloatMatrixF< 6, 6 > &D,
                                            const FloatArrayF< 6 > &strain) const;

    /**
     * Check if the trial stress state falls within the vertex region of the plasticity model at the apex of triaxial extension or triaxial compression.
     * @returns true for vertex case and false if regular stress return can be used.
     * @param answer Volumetric apex stress.
     * @param sig Volumetric stress.
     * @param tempKappa Hardening variable.
     * @param gp Gauss point
     */
    void checkForVertexCase(double &answer,
                            ConcreteDPM2_ReturnType &returnType,
                            double sig,
                            double tempKappa,
                            GaussPoint *gp) const;

    /**
     * Perform regular stress return for the plasticity model, i.e. if the trial stress state does not lie in the vertex region.
     * @param stress Stress vector which is computed.
     * @param kappaP Initial guess for kappa P (i.e. previous kappa)
     * @param gp Gauss point.
     * @param theta Load angle of trial stress (remains constant throughout return).
     */
    double performRegularReturn(FloatArrayF< 6 > &stress,
                                ConcreteDPM2_ReturnResult &returnResult,
                                ConcreteDPM2_ReturnType &returnType,
                                double kappaP,
                                GaussPoint *gp,
                                double theta) const;

    /**
     * Compute jacobian for 2D(plane strain) and 3d cases
     * @param sig volumetric strain
     * @param rho deviatoric
     * @param tempKappa plastic strain
     * @param deltaLambda plastic multiplier
     * @param gp Gauss point
     */
    FloatMatrixF< 4, 4 >computeJacobian(double sig,
                                        double rho,
                                        double theta,
                                        double tempKappa,
                                        double deltaLambda,
                                        GaussPoint *gp) const;

    /**
     * Perform stress return for vertex case of the plasticity model, i.e. if the trial stress state lies within the vertex region.
     * @param stress Stress vector of this Gauss point.
     * @param apexStress Volumetric stress at the apex of the yield surface.
     * @param theta Lode angle.
     * @param tempKappaP temporary cummulative plastic strain
     * @param gp Gauss point.
     * @returns updated temporary cummulative plastic strain
     */
    double performVertexReturn(FloatArrayF< 6 > &stress,
                               ConcreteDPM2_ReturnResult &returnResult,
                               ConcreteDPM2_ReturnType &returnType,
                               double apexStress,
                               double tempKappaP,
                               GaussPoint *gp) const;

    /**
     * Compute the yield value based on stress and hardening variable.
     * @param sig Volumetric stress.
     * @param rho Length of the deviatoric stress.
     * @param theta Lode angle of the stress state.
     * @param tempKappa Hardening variable.
     * @returns Yield value.
     */
    double computeYieldValue(double sig,
                             double rho,
                             double theta,
                             double tempKappa) const;

    /**
     * Compute the value of the hardening function based on the hardening variable.
     * @param tempKappa Hardening variable.
     * @returns Value of the hardening function.
     */
    double computeHardeningOne(double tempKappa) const;

    /** Compute the derivative of the hardening function based on the
     *  hardening parameter.
     *  @param tempKappa Hardening variable.
     *  @return Derivative of the hardening function.
     */
    double computeHardeningOnePrime(double tempKappa) const;

    /**
     * Compute the value of the hardening function based on the hardening variable.
     * @param tempKappa Hardening variable.
     * @return Value of the hardening function.
     */
    double computeHardeningTwo(double tempKappa) const;

    /**
     * Compute the derivative of the hardening function based on the
     * hardening parameter.
     * @param tempKappa Hardening variable.
     * @return Derivative of the hardening function.
     */
    double computeHardeningTwoPrime(double tempKappa) const;

    /** Compute the derivative of the yield surface with respect to the hardening
     *  variable based on the stress state and the hardening variable
     *  @param sig Volumetric stress.
     *  @param rho Deviatoric length.
     *  @param theta Lode angle.
     *  @param tempKappa Hardening variable.
     *  @return Derivative of the yield surface.
     */
    double computeDFDKappa(double sig,
                           double rho,
                           double theta,
                           double tempKappa) const;

    /**
     * Compute the derivative of kappa with respect of delta lambda based on the stress state and the hardening variable.
     * @param sig Volumetric stress.
     * @param rho Length of the deviatoric stress.
     * @param theta Lode angle
     * @param tempKappa Hardening variable.
     * @return Derivative of kappa with respect to delta lambda.
     */
    double computeDKappaDDeltaLambda(double sig, double rho, double theta, double tempKappa) const;


    /**
     * Compute the ductility measure based on the stress state.
     * @param sig Volumetric stress.
     * @param rho Length of the deviatoric strength.
     * @param theta Lode angle of stress state.
     * @return Ductility measure.
     */
    virtual double computeDuctilityMeasure(double sig,
                                           double rho,
                                           double theta) const;

    /**
     * Compute derivative the ductility measure with respect to  the stress state.
     * @param answer array of the derivative of the ductility measure with respect to volumetric and deviatoric stress
     * @param sig Volumetric stress.
     * @param rho Length of the deviatoric strength.
     * @param theta Lode angle.
     * @param tempKappa plastic strain
     */
    FloatArrayF< 2 >computeDDuctilityMeasureDInv(double sig,
                                                 double rho,
                                                 double theta,
                                                 double tempKappa) const;

    /**
     * Compute derivative the palstic potential function with respect to  the stress state.
     * @param answer array of the derivative of the plastic potential with respect to volumetric and deviatoric stress
     * @param sig volumetric stress.
     * @param rho deviatoric stress.
     * @param tempKappa plastic strain
     */
    FloatArrayF< 2 >computeDGDInv(double sig,
                                  double rho,
                                  double tempKappa) const;

    /**
     * This function computes the ratio of the volumetric and deviatoric component
     * of the flow direction. It is used within the vertex return to check,
     * if the vertex return is admissible.
     */
    double computeRatioPotential(double sig,
                                 double rho,
                                 double tempKappa) const;

    /**
     * This function computes the rate factor which is used to take into account the strain rate dependence of the material.
     */
    double computeRateFactor(double alpha,
                             double timeFactor,
                             GaussPoint *gp,
                             TimeStep *deltaTime) const;

    /**
     * Computes the second derivative of the plastic potential with respect to the
     * invariants sig and rho are computed.
     */
    FloatMatrixF< 2, 2 >computeDDGDDInv(double sig,
                                        double rho,
                                        double tempKappa) const;


    /**
     * Computes the second derivative of the plastic potential with respect to the
     * stress.
     */
    FloatMatrixF< 6, 6 >computeDDGDDStress(const FloatArrayF< 6 > &stress,
                                           const double tempKappa) const;


    /**
     * Computes the derivative of cos theta with respect to the
     * stress.
     */
    FloatArrayF< 6 >computeDCosThetaDStress(const FloatArrayF< 6 > &stress) const;

    /**
     * Computes the second derivative of cos theta with respect to the
     * stress.
     */
    FloatMatrixF< 6, 6 >computeDDCosThetaDDStress(const FloatArrayF< 6 > &stress) const;

    /**
     * The mixed derivative of the plastic potential with respect
     * to the invariants and the hardening parameter are determined.
     */
    FloatArrayF< 2 >computeDDGDInvDKappa(double sig,
                                         double rho,
                                         double tempKappa) const;


    /**
     * Computes the mixed derivative of the hardening parameter kappa with
     * respect to the plastic multiplier delta Lambda and the invariants sig
     * and rho.
     */
    FloatArrayF< 2 >computeDDKappaDDeltaLambdaDInv(double sig,
                                                   double rho,
                                                   double theta,
                                                   double tempKappa) const;

    /**
     * Computes the derivative of evolution law of hardening variable with respect to the
     * stress.
     */
    FloatArrayF< 6 >computeDDKappaDDeltaLambdaDStress(const FloatArrayF< 6 > &stress, double tempKappa) const;


    /**
     * Computes the second mixed derivative of plastic potential with respect to the
     * stress and kappa.
     */
    FloatArrayF< 6 >computeDDGDStressDKappa(const FloatArrayF< 6 > &stress, double tempKappa) const;

    /**
     * Computes the derivative of the evolution law of the hardening parameter kappa with respect to the hardening variable kappa.
     */
    double computeDDKappaDDeltaLambdaDKappa(double sig, double rho, double theta, double tempKappa) const;


    /**
     * Computes the derivative of the yield surface with respect to the
     * invariants sig and rho.
     */
    FloatArrayF< 2 >computeDFDInv(double sig,
                                  double rho,
                                  double theta,
                                  double tempKappa) const;

    /**
     * Computes the derivative of the yield surface with respect to the
     * stress.
     */
    FloatArrayF< 6 >computeDFDStress(const FloatArrayF< 6 > &stress,
                                     const double tempKappa) const;

    /**
     * Computes the derivative of the plastic potential with respect to the
     * stress.
     */
    FloatArrayF< 6 >computeDGDStress(const FloatArrayF< 6 > &stress,
                                     const double tempKappa) const;

    /**
     * Computes full Jacobian used for the algorithmic tangent stiffness
     */
    FloatMatrixF< 8, 8 >computeFullJacobian(const FloatArrayF< 6 > &stress,
                                            const double deltaLambda,
                                            GaussPoint *gp,
                                            TimeStep *atTime,
                                            const double tempKappa) const;

    /// Compute tempKappa.
    double computeTempKappa(double kappaInitial,
                            double sigTrial,
                            double rhoTrial,
                            double sig) const;

    /// Compute damage parameters
    FloatArrayF< 2 >computeDamage(const FloatArrayF< 6 > &strain, const FloatMatrixF< 6, 6 > &D, double timeFactor, GaussPoint *gp, TimeStep *tStep, double alpha, const FloatArrayF< 6 > &effectiveStress) const;


    /// Check for un- and reloading in the damage part
    int checkForUnAndReloading(double &tempEquivStrain,
                               double &minEquivStrain,
                               const FloatMatrixF< 6, 6 > &D,
                               GaussPoint *gp) const;

    /// Compute alpha for rate effect
    double computeAlpha(FloatArrayF< 6 > &effectiveStressTension, FloatArrayF< 6 > &effectiveStressCompression, const FloatArrayF< 6 > &effectiveStress) const;

    /// Compute damage parameter in tension.
    virtual double computeDamageParamTension(double equivStrain, double kappaOne, double kappaTwo, double le, double omegaOld, double rateFactor) const;

    /// Compute damage parameter in compression.
    virtual double computeDamageParamCompression(double equivStrain, double kappaOne, double kappaTwo, double omegaOld, double rateFactor) const;

    /// Compute equivalent strain value for tension.
    double computeDeltaPlasticStrainNormTension(double tempKappaD, double kappaD, GaussPoint *gp) const;

    /// Compute equivalent strain value for compression.
    double computeDeltaPlasticStrainNormCompression(double tempAlpha, double tempKappaD, double kappaD, GaussPoint *gp, const double rho) const;

    /// Compute the base equivalent strain value.
    virtual double computeEquivalentStrain(double sig, double rho, double theta) const;

    /// Compute the ductility measure for the damage model.
    double computeDuctilityMeasureDamage(GaussPoint *gp, const double sig, const double rho) const;

    /**
     * Initialize the characteristic length, if damage is not yet activated
     * Set the increase factor for the strain rate dependence
     */
    void initDamaged(double kappa,
                     const FloatArrayF< 6 > &strain,
                     GaussPoint *gp) const;

    /// Compute the Haigh-Westergaard coordinates.
    void computeCoordinates(const FloatArrayF< 6 > &stress, double &sig, double &rho, double &theta) const;


    /// Assign state flag.
    void assignStateFlag(GaussPoint *gp) const;

    /// Computes the derivative of rho with respect to the stress.
    FloatArrayF< 6 >computeDRhoDStress(const FloatArrayF< 6 > &stress) const;

    /// Computes the derivative of function r with respect to cos theta
    double computeDRDCosTheta(const double theta, const double ecc) const;

    /// Computes the second derivative of function r with respect to cos theta
    double computeDDRDDCosTheta(const double theta, const double ecc) const;


    /// Computes the derivative of sig with respect to the stress.
    FloatArrayF< 6 >computeDSigDStress() const;

    /// Computes the second derivative of rho with the respect to the stress.
    FloatMatrixF< 6, 6 >computeDDRhoDDStress(const FloatArrayF< 6 > &stress) const;

    FloatMatrixF< 6, 6 >give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    /// Compute the 3d secant stiffness matrix.
    FloatMatrixF< 6, 6 >compute3dSecantStiffness(GaussPoint *gp, TimeStep *tStep) const;

    /// Compute the 3d tangent stiffness matrix.
    FloatMatrixF< 6, 6 >compute3dTangentStiffness(GaussPoint *gp, TimeStep *tStep) const;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;


protected:
    MaterialStatus *CreateStatus(GaussPoint *gp) const override;
};
} //end namespace oofem
#endif
