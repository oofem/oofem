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

#include "../sm/Materials/structuralmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "cltypes.h"
#include "../sm/Materials/structuralms.h"
#include "Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"

#define DYNCON_TOL 1.e-6
#define keep_track_of_dissipated_energy
///@name Input fields for ConcreteDPM2
//@{
#define _IFT_ConcreteDPM2_Name "con2dpm"
#define _IFT_ConcreteDPM2_fc "fc"
#define _IFT_ConcreteDPM2_fcZero "fczero"
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
#define _IFT_ConcreteDPM2_rateFlag "rateflag"
#define _IFT_ConcreteDPM2_timeFactor "timefactor"
#define _IFT_ConcreteDPM2_helem "helem"
#define _IFT_ConcreteDPM2_isoflag "isoflag"
//@}

namespace oofem {
/**
 * This class implements the material status associated to ConcreteDPM2.
 * @author Peter Grassl
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
    FloatArray plasticStrain;
    FloatArray tempPlasticStrain;

    double dFDKappa;
    double deltaLambda;
    //@}

    /// @name Hardening variable
    //@{
    double kappaP;
    double tempKappaP;
    //@}

    double kappaPPeak;

    double le;

    double alpha;
    double tempAlpha;

    double equivStrain;
    double tempEquivStrain;

    double equivStrainTension;
    double tempEquivStrainTension;

    double equivStrainCompression;
    double tempEquivStrainCompression;

    double kappaDTension;
    double tempKappaDTension;

    double kappaDCompression;
    double tempKappaDCompression;

    double kappaDTensionOne;
    double tempKappaDTensionOne;

    double kappaDCompressionOne;
    double tempKappaDCompressionOne;

    double kappaDTensionTwo;
    double tempKappaDTensionTwo;

    double kappaDCompressionTwo;
    double tempKappaDCompressionTwo;

    double damageTension;
    double tempDamageTension;

    double damageCompression;
    double tempDamageCompression;

    double deltaEquivStrain;

    double rateFactor;
    double tempRateFactor;

    /// Strains that are used for calculation of strain rates
    double rateStrain;
    double tempRateStrain;

    /// Indicates the state (i.e. elastic, unloading, plastic, damage, vertex) of the Gauss point
    int state_flag;
    int temp_state_flag;
#ifdef keep_track_of_dissipated_energy
    /// Density of total work done by stresses on strain increments.
    double stressWork;
    /// Non-equilibrated density of total work done by stresses on strain increments.
    double tempStressWork;
    /// Density of dissipated work.
    double dissWork;
    /// Non-equilibrated density of dissipated work.
    double tempDissWork;
#endif

public:
    /// Constructor
    ConcreteDPM2Status(int n, Domain *d, GaussPoint *gp);

    /// Destructor
    virtual ~ConcreteDPM2Status();

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual contextIOResultType saveContext(DataStream &, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &, ContextMode mode, void *obj = NULL);
    virtual const char *giveClassName() const { return "ConcreteDPM2Status"; }

    // Inline functions for access to state variables

    /**
     * Get the plastic strain deviator from the material status.
     * @return Plastic strain deviator.
     */
    const FloatArray &givePlasticStrain() const { return plasticStrain; }


    /**
     * Get the deviatoric plastic strain norm from the material status.
     * @return DeviatoricPlasticStrainNorm.
     */
    double giveDeviatoricPlasticStrainNorm()
    {
        FloatArray dev;
        StructuralMaterial :: computeDeviatoricVolumetricSplit(dev, plasticStrain);
        return sqrt( .5 * ( 2. * dev [ 0 ] * dev [ 0 ] + 2. * dev [ 1 ] * dev [ 1 ] + 2. * dev [ 2 ] * dev [ 2 ] +
                    dev [ 3 ] * dev [ 3 ] + dev [ 4 ] * dev [ 4 ] + dev [ 5 ] * dev [ 5 ] ) );
    }

    /**
     * Get the volumetric plastic strain from the material status.
     * @return volumetricPlasticStrainNorm.
     */
    double giveVolumetricPlasticStrain() const
    {
        return 1. / 3. * ( plasticStrain(0) + plasticStrain(1) + plasticStrain(2) );
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
    const FloatArray &giveTempPlasticStrain() const { return tempPlasticStrain; }

    /**
     *  Get the temp value of the volumetric plastic strain in plane stress
     */
    double giveTempVolumetricPlasticStrain() const
    { return 1. / 3. * ( tempPlasticStrain(0) + tempPlasticStrain(1) + tempPlasticStrain(2) ); }


    /**
     * Get the temp value of the hardening variable of the plasticity model
     * from the material status.
     * @return Temp value of hardening variable kappaP.
     */
    double giveTempKappaP() const
    { return tempKappaP; }


    /**
     * Get the temp value of the hardening variable of the damage model
     * from the material status.
     * @return Temp value of the damage variable damage.
     */
    double giveKappaDTension() const
    { return kappaDTension; }

    double giveAlpha() const
    { return alpha; }


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
    void letTempPlasticStrainBe(const FloatArray &v)
    { tempPlasticStrain = v; }

    /**
     * Assign the temp value of the rate factor of the damage model.
     * @param v New temp value of the damage variable
     */
    void letDeltaLambdaBe(double v)
    { deltaLambda = v; }

    /**
     * Assign the temp value of the hardening variable of the plasticity model.
     * @param v New temp value of the hardening variable
     */
    void letTempKappaPBe(double v)
    { tempKappaP = v; }

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
    double giveLe() { return le; }

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


//   ********************************
//   *** CLASS DYNAMIC CONCRETE   ***
//   ********************************

/**
 * This class contains the combination of a local plasticity model for concrete with a local isotropic damage model.
 * This is an extension of concretedpm2. The yield surface of the plasticity model is based on the extension of the Menetrey and Willam yield criterion.
 * The flow rule is nonassociated. The evolution laws of the hardening variables depend on the stress state.
 * The plasticity model describes only hardening and perfect plasticity.
 * It is based on the effective stress. The damage parameter of the isotropic damage model is based on the total volumetric strain.
 * An exponential softening law is implemented.
 */
class ConcreteDPM2 : public StructuralMaterial
{
public:

protected:
    enum ConcreteDPM2_ReturnType { RT_Regular, RT_Tension, RT_Compression, RT_Auxiliary };
    ConcreteDPM2_ReturnType returnType;

    enum ConcreteDPM2_ReturnResult { RR_NotConverged, RR_Converged };
    ConcreteDPM2_ReturnResult returnResult;

    /// Parameters of the yield surface of the plasticity model. fc is the uniaxial compressive strength, ft the uniaxial tensile strength and ecc controls the out of roundness of the deviatoric section.
    double fc, ft, ecc;

    int isotropicFlag;

    double e0;

    /// Parameter of the ductilityMeasure of the plasticity model.
    double AHard;
    /// Parameter of the ductilityMeasure of the plasticity model.
    double BHard;
    /// Parameter of the ductilityMeasure of the plasticity model.
    double CHard;
    /// Parameter of the ductilityMeasure of the plasticity model.
    double DHard;

    /// Hardening modulus.
    double hardeningModulus;

    /// Parameter of the ductilityMeasure of the damage model.
    double ASoft;

    /// Parameter of the hardening law of the plasticity model.
    double yieldHardPrimePeak;

    /// Parameter of the hardening law of the plasticity model.
    double yieldHardInitial;

    /// Control parameter for te volumetric plastic flow of the plastic potential.
    double dilationConst;

    ///@todo These values should not be stored by the material model itself. As temp values in the status they would be OK, but this will be thread unsafe (and makes the model into a complete spagetti code) / Mikael
#if 1
    /// Volumetric stress.
    double sig;

    /// Length of the deviatoric stress.
    double rho;

    /// Lode angle of the trial stress..
    double thetaTrial;
#endif

    /// Friction parameter of the yield surface.
    double m;

    /// Dilation parameter of the plastic potential.
    double mQ;

    /// Element size (to be used in fracture energy approach (crack band).
    double helem;

    /// Pointer for linear elastic material.
    LinearElasticMaterial *linearElasticMaterial;

    /// Elastic Young's modulus.
    double eM;
    /// Elastic shear modulus.
    double gM;
    /// Elastic bulk modulus.
    double kM;
    /// Elastic poisson's ration.
    double nu;

    /// Control parameter for the exponential softening law.
    double efCompression;

    /// Control parameter for the linear/bilinear softening law in tension.
    double wf;

    /// Control parameter for the bilinear softening law in tension.
    double wfOne;

    /// Control parameter for the bilinear softening law.
    double ftOne;

    /// yield tolerance for the plasticity model.
    double yieldTol;

    /// yield tolerance for the damage model.
    double yieldTolDamage;

    /// Maximum number of iterations for stress return.
    int newtonIter;

    /// Type of softening function used.
    int softeningType;

    /// Input parameter which simulates a loading rate. Only for debugging purposes.
    double timeFactor;

    /// This parameter is needed for the rate dependence. It should be read in if rate dependence is considered.
    double fcZero;

    /// Flag which signals if strainRate effects should be considered.
    int strainRateFlag;


public:
    /// Constructor
    ConcreteDPM2(int n, Domain *d);
    /// Destructor
    virtual ~ConcreteDPM2();
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual const char *giveClassName() const { return "ConcreteDPM2"; }
    virtual const char *giveInputRecordName() const { return _IFT_ConcreteDPM2_Name; }

    LinearElasticMaterial *giveLinearElasticMaterial()
    { return linearElasticMaterial; }

    virtual void giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *tStep);

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &strainVector, TimeStep *tStep);

    /**
     * Perform stress return of the plasticity model and compute history variables.
     * @param gp Gauss point.
     * @param strain Strain vector of this Gauss point.
     */
    void performPlasticityReturn(GaussPoint *gp,
                                 const FloatMatrix &D,
                                 const FloatArray &strain);

    /**
     * Check if the trial stress state falls within the vertex region of the plasticity model at the apex of triaxial extension or triaxial compression.
     * @returns true for vertex case and false if regular stress return can be used.
     * @param answer Volumetric apex stress.
     * @param sig Volumetric stress.
     * @param tempKappa Hardening variable.
     */
    bool checkForVertexCase(double &answer,
                            double sig,
                            double tempKappa,
                            bool mode1d);

    /**
     * Perform regular stress return for the plasticity model, i.e. if the trial stress state does not lie in the vertex region.
     * @param stress Stress vector which is computed.
     * @param kappaP Initial guess for kappa P (i.e. previous kappa)
     * @param gp Gauss point.
     */
    double performRegularReturn(FloatArray &stress,
                                double kappaP,
                                GaussPoint *gp);

    /**
     * Compute jacobian for 1D case
     * @param totalsigma stress value
     * @param tempKappa plastic strain
     * @param deltaLambda plastic multiplier
     * @param gp Gauss point
     */

    void compute1dJacobian(FloatMatrix &answer,
                           double totalsigma,
                           double tempKappa,
                           double deltaLambda,
                           GaussPoint *gp);
    /**
     * Compute jacobian for 2D(plane strain) and 3d cases
     * @param sig volumetric strain
     * @param rho deviatoric
     * @param tempKappa plastic strain
     * @param deltaLambda plastic multiplier
     * @param gp Gauss point
     */
    void computeJacobian(FloatMatrix &answer,
                         double sig,
                         double rho,
                         double tempKappa,
                         double deltaLambda,
                         GaussPoint *gp);

    /**
     * Perform stress return for vertex case of the plasticity model, i.e. if the trial stress state lies within the vertex region.
     * @param stress Stress vector of this Gauss point.
     * @param apexStress Volumetric stress at the apex of the yield surface.
     * @param tempKappaP temporary cummulative plastic strain
     * @param gp Gauss point.
     * @returns updated temporary cummulative plastic strain
     */
    double performVertexReturn(FloatArray &stress,
                               double apexStress,
                               double tempKappaP,
                               GaussPoint *gp);

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
     *  @param tempKappa Hardening variable.
     *  @return Derivative of the yield surface.
     */
    double computeDFDKappa(double sig,
                           double rho,
                           double tempKappa,
                           bool mode1d);


    /**
     * Compute the derivative of kappa with respect of delta lambda based on the stress state and the hardening variable.
     * @param sig Volumetric stress.
     * @param rho Length of the deviatoric stress.
     * @param tempKappa Hardening variable.
     * @return Derivative of kappa with respect to delta lambda.
     */
    double computeDKappaDDeltaLambda(double sig, double rho, double tempKappa);
    double computeDKappaDDeltaLambda1d(double sig, double tempKappa);


    /**
     * Compute the ductility measure based on the stress state.
     * @param sig Volumetric stress.
     * @param rho Length of the deviatoric strength.
     * @param theta Lode angle of stress state.
     * @return Ductility measure.
     */
    virtual double computeDuctilityMeasure(double sig,
                                           double rho,
                                           double theta);



    /**
     * Compute derivative the ductility measure with respect to  the stress state.
     * @param answer array of the derivative of the ductility measure with respect to volumetric and deviatoric stress
     * @param sig Volumetric stress.
     * @param rho Length of the deviatoric strength.
     * @param tempKappa plastic strain
     */
    void computeDDuctilityMeasureDInv(FloatArray &answer,
                                      double sig,
                                      double rho,
                                      double tempKappa);
    /**
     * Compute derivative the ductility measure with respect to  the stress state.
     * @return The derivative of the ductility measure with respect to stress
     * @param sig stress.
     * @param tempKappa plastic strain
     */
    double computeDDuctilityMeasureDInv1d(double sigma, double tempKappa); //Dimitris change 1d implementation

    /**
     * Compute derivative the palstic potential function with respect to  the stress state.
     * @param answer array of the derivative of the plastic potential with respect to volumetric and deviatoric stress
     * @param sig volumetric stress.
     * @param rho deviatoric stress.
     * @param tempKappa plastic strain
     */
    void computeDGDInv(FloatArray &answer,
                       double sig,
                       double rho,
                       double tempKappa);
    /**
     * Compute derivative the palstic potential function with respect to  the stress state.
     * @param return The derivative of the plastic potential with respect to stress
     * @param sig stress.
     * @param tempKappa plastic strain
     */
    double computeDGDInv1d(double sig, double tempKappa);

    /**
     * This function computes the ratio of the volumetric and deviatoric component
     * of the flow direction. It is used within the vertex return to check,
     * if the vertex return is admissible.
     */
    double computeRatioPotential(double sig,
                                 double tempKappa);

    /**
     * This function computes the rate factor which is used to take into account the strain rate dependence of the material.
     */
    double computeRateFactor(double alpha,
                             double timeFactor,
                             GaussPoint *gp,
                             TimeStep *deltaTime);

    /**
     * Here, the second derivative of the plastic potential with respect to the
     * invariants sig and rho are computed.
     */
    void computeDDGDDInv(FloatMatrix &answer,
                         double sig,
                         double rho,
                         double tempKappa);
    /**
     * Here, the second derivative of the plastic potential with respect to the
     * invariants sig and rho are computed.
     */
    double computeDDGDDInv1d(double sigma, double tempKappa);
    /**
     * Here, the mixed derivative of the plastic potential with respect
     * to the invariants and the hardening parameter are determined.
     */
    void computeDDGDInvDKappa(FloatArray &answer,
                              double sig,
                              double rho,
                              double tempKappa);

    double computeDDGDInvDKappa1d(double sigma, double tempKappa);
    /**
     * Computes the mixed derivative of the hardening parameter kappa with
     * respect to the plastic multiplier delta Lambda and the invariants sig
     * and rho.
     */
    void computeDDKappaDDeltaLambdaDInv(FloatArray &answer,
                                        double sig,
                                        double rho,
                                        double tempKappa);

    double computeDDKappaDDeltaLambdaDInv1d(double sigma, double tempKappa);
    /**
     * Computes the derivative of the evolution law of the hardening parameter kappa with respect to the hardening variable kappa.
     */
    double computeDDKappaDDeltaLambdaDKappa(double sig, double rho, double tempKappa);
    double computeDDKappaDDeltaLambdaDKappa1d(double sig, double tempKappa);


    /**
     * Computes the derivative of the yield surface with respect to the
     * invariants sig and rho.
     */
    void computeDFDInv(FloatArray &answer,
                       double sig,
                       double rho,
                       double tempKappa) const;
    double computeDFDInv1d(double sigma, double tempKappa) const;
    /**
     * Compute tempKappa.
     */
    double computeTempKappa(double kappaInitial,
                            double sigTrial,
                            double rhoTrial,
                            double sig);


    /**
     * Compute damage parameters
     */
    void  computeDamage(FloatArray &answer, const FloatArray &strain, const FloatMatrix &D, double timeFactor, GaussPoint *gp, TimeStep *tStep, double alpha);


    /**
     * Check for un- and reloading in the damage part
     */
    int checkForUnAndReloading(double &tempEquivStrain,
                               double &minEquivStrain,
                               const FloatMatrix &D,
                               GaussPoint *gp);

    double computeAlpha(FloatArray &effectiveStressTension, FloatArray &effectiveStressCompression, FloatArray &effectiveStress);

    /// Compute damage parameter.
    virtual double computeDamageParamTension(double equivStrain, double kappaOne, double kappaTwo, double le, double omegaOld);

    virtual double computeDamageParamCompression(double equivStrain, double kappaOne, double kappaTwo, double omegaOld);

    /// Compute equivalent strain value.
    double computeDeltaPlasticStrainNormTension(double tempKappaD, double kappaD, GaussPoint *gp);

    double computeDeltaPlasticStrainNormCompression(double tempAlpha, double tempKappaD, double kappaD, GaussPoint *gp);

    virtual double computeEquivalentStrain(double sig, double rho, double theta);


    /// Compute the ductility measure for the damage model.
    double computeDuctilityMeasureDamage(const FloatArray &strain, GaussPoint *gp);

    /**
     * Initialize the characteristic length, if damage is not yet activated
     * Set the increase factor for the strain rate dependence
     */
    void initDamaged(double kappa,
                     const FloatArray &strain,
                     GaussPoint *gp);


    /// Compute the trial coordinates.
    void computeTrialCoordinates(const FloatArray &stress, double &sig, double &rho, double &theta);


    /// Assign state flag.
    void assignStateFlag(GaussPoint *gp);

    /// Computes the derivative of rho with respect to the stress.
    void computeDRhoDStress(FloatArray &answer,
                            const FloatArray &stress) const;

    /// Computes the derivative of sig with respect to the stress.
    void computeDSigDStress(FloatArray &answer) const;

    /// Computes the seconfd derivative of rho with the respect to the stress.
    void computeDDRhoDDStress(FloatMatrix &answer,
                              const FloatArray &stress) const;

    /// Computes the derivative of costheta with respect to the stress.
    void computeDCosThetaDStress(FloatArray &answer,
                                 const FloatArray &stress) const;

    /// Compute the derivative of R with respect to costheta.
    double computeDRDCosTheta(double theta, double ecc) const;

    virtual void give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    void compute3dSecantStiffness(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep);



    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual int giveIPValue(FloatArray &answer,
                            GaussPoint *gp,
                            InternalStateType type,
                            TimeStep *tStep);

protected:
    MaterialStatus *CreateStatus(GaussPoint *gp) const;
};
} //end namespace oofem
#endif
