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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef ConcreteDPM2_h
#define ConcreteDPM2_h

#include "structuralmaterial.h"
#include "gausspnt.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"
#include "structuralms.h"
#include "strainvector.h"
#include "stressvector.h"
#include "isolinearelasticmaterial.h"

#define DYNCON_TOL 1.e-6

namespace oofem {
/**
 * This class implements the material status associated to ConcreteDPM2.
 * @author Peter Grassl
 */
class ConcreteDPM2Status : public StructuralMaterialStatus
{
public:
    /// values of history variable state_flag
    enum state_flag_values { ConcreteDPM2_Elastic, ConcreteDPM2_Unloading, ConcreteDPM2_Plastic, ConcreteDPM2_Damage, ConcreteDPM2_PlasticDamage,  ConcreteDPM2_VertexCompression, ConcreteDPM2_VertexTension, ConcreteDPM2_VertexCompressionDamage, ConcreteDPM2_VertexTensionDamage };
protected:
    /// @name History variables of the plasticity model
    //@{
    StrainVector plasticStrain;
    StrainVector tempPlasticStrain;

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

    double damage;
    double tempDamage;

    double deltaEquivStrain;

    double rateFactor;
    double tempRateFactor;

    /// Strains that are used for calculation of strain rates
    double rateStrain;
    double tempRateStrain;

    /// Indicates the state (i.e. elastic, unloading, plastic, damage, vertex) of the Gauss point
    int state_flag;
    int temp_state_flag;


public:
    /// Constructor
    ConcreteDPM2Status(int n, Domain *d, GaussPoint *gp);

    /// Destructor
    virtual ~ConcreteDPM2Status();

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual contextIOResultType saveContext(DataStream *, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *, ContextMode mode, void *obj = NULL);
    virtual const char *giveClassName() const { return "ConcreteDPM2Status"; }
    virtual classType giveClassID() const { return ConcreteDPM2StatusClass; }

    // Inline functions for access to state variables
    // give:
    // Functions used to access the value of a internal variable.
    /**
     *   Get the full plastic strain vector from the material status.
     *  @param answer Plastic strain vector.
     */
    void  giveFullPlasticStrainVector(StrainVector &answer) const
    {
        StrainVector plasticStrainVector(_Unknown);
        givePlasticStrain(plasticStrainVector);
        plasticStrainVector.convertToFullForm(answer);
    }
    /**
     *   Get the plastic strain deviator from the material status.
     *  @param answer Plastic strain deviator.
     */
    void  givePlasticStrain(StrainVector &answer) const
    { answer = plasticStrain; }


    /**
     *  Get the deviatoric plastic strain norm from the material status.
     *  @return Deviatoric plasticStrainNorm.
     */
    double giveDeviatoricPlasticStrainNorm()
    {
        StrainVector deviatoricPlasticStrain( gp->giveMaterialMode() );
        double volumetricPlasticStrain;
        plasticStrain.computeDeviatoricVolumetricSplit(deviatoricPlasticStrain,
                                                       volumetricPlasticStrain);
        return deviatoricPlasticStrain.computeStrainNorm();
    }


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
     * Get the hardening variable of the damage model from the
     * material status.
     * @return Hardening variable kappaD.
     */
    double giveKappaDCompressionOne() const
    { return kappaDCompressionOne; }


    /**
     * Get the hardening variable of the damage model from the
     * material status.
     * @return Hardening variable kappaD.
     */
    double giveKappaDTensionTwo() const
    { return kappaDTensionTwo; }


    /**
     * Get the hardening variable of the damage model from the
     * material status.
     * @return Hardening variable kappaD.
     */
    double giveKappaDCompressionTwo() const
    { return kappaDCompressionTwo; }

    double giveEquivStrain() const
    { return equivStrain; }

    double giveEquivStrainTension() const
    { return equivStrainTension; }

    double giveEquivStrainCompression() const
    { return equivStrainCompression; }

    /**
     * Get the damage variable of the damage model from the
     * material status.
     * @return Damage variable damage.
     */
    double giveDamageTension() const
    { return damageTension; }


    /**
     * Get the damage variable of the damage model from the
     * material status.
     * @return Damage variable damage.
     */
    double giveDamageCompression() const
    { return damageCompression; }


    /**
     * Get the damage variable of the damage model from the
     * material status.
     * @return Damage variable damage.
     */
    double giveDamage() const
    { return damage; }


    /**
     * Get the damage variable of the damage model from the
     * material status.
     * @return Damage variable damage.
     */
    double giveRateFactor() const
    { return rateFactor; }

    double giveTempRateFactor() const
    { return tempRateFactor; }


    double giveRateStrain() const
    { return rateStrain; }

    void letTempRateStrainBe(const double v)
    { tempRateStrain = v; }

    void letTempAlphaBe(const double v)
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
     * @param answer Temp value of plastic strain vector.
     */
    void giveTempPlasticStrain(StrainVector &answer) const
    { answer = tempPlasticStrain; }



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
    double giveTempDamage() const
    { return tempDamage; }


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
    void letTempPlasticStrainBe(const StrainVector &v)
    { tempPlasticStrain = v; }

    /**
     * Assign the temp value of the rate factor of the damage model.
     * @param v New temp value of the damage variable
     */
    void letDeltaLambdaBe(const double v)
    { deltaLambda = v; }

    /**
     * Assign the temp value of the hardening variable of the plasticity model.
     * @param v New temp value of the hardening variable
     */
    void letTempKappaPBe(const double v)
    { tempKappaP = v; }

    /**
     * Assign the temp value of the rate factor of the damage model.
     * @param v New temp value of the damage variable
     */
    void letTempKappaDTensionBe(const double v)
    { tempKappaDTension = v; }

    /**
     * Assign the temp value of the rate factor of the damage model.
     * @param v New temp value of the damage variable
     */
    void letTempKappaDCompressionBe(const double v)
    { tempKappaDCompression = v; }

    /**
     * Assign the temp value of the hardening variable of the damage model.
     * @param v New temp value of the hardening variable
     */
    void letTempKappaDTensionOneBe(const double v)
    { tempKappaDTensionOne = v; }

    /**
     * Assign the temp value of the hardening variable of the damage model.
     * @param v New temp value of the hardening variable
     */
    void letTempKappaDCompressionOneBe(const double v)
    { tempKappaDCompressionOne = v; }

    /**
     * Assign the temp value of the second tension hardening variable of the damage model.
     * @param v New temp value of the second tension hardening variable
     */
    void letTempKappaDTensionTwoBe(const double v)
    { tempKappaDTensionTwo = v; }

    /**
     * Assign the temp value of the second compression hardening variable of the damage model.
     * @param v New temp value of the second compression hardening variable
     */
    void letTempKappaDCompressionTwoBe(const double v)
    { tempKappaDCompressionTwo = v; }

    /**
     * Assign the temp value of the tensile damage variable of the damage model.
     * @param v New temp value of the tensile damage variable
     */
    void letTempDamageTensionBe(const double v)
    { tempDamageTension = v; }

    /**
     * Assign the temp value of the compressive damage variable of the damage model.
     * @param v New temp value of the compressive damage variable
     */
    void letTempDamageCompressionBe(const double v)
    { tempDamageCompression = v; }

    /**
     * Assign the temp value of the rate factor of the damage model.
     * @param v New temp value of the damage variable
     */
    void letTempRateFactorBe(const double v)
    { tempRateFactor = v; }

    /**
     * Assign the temp value of the rate factor of the damage model.
     * @param v New temp value of the damage variable
     */
    void letTempEquivStrainBe(const double v)
    { tempEquivStrain = v; }

    /**
     * Assign the temp value of the rate factor of the damage model.
     * @param v New temp value of the damage variable
     */
    void letTempEquivStrainTensionBe(const double v)
    { tempEquivStrainTension = v; }

    /**
     * Assign the temp value of the rate factor of the damage model.
     * @param v New temp value of the damage variable
     */
    void letTempEquivStrainCompressionBe(const double v)
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
};


//   ********************************
//   *** CLASS DYNAMIC CONCRETE   ***
//   ********************************

#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"
#include "structuralmaterial.h"
#include "isolinearelasticmaterial.h"
#include "strainvector.h"
#include "stressvector.h"

/**
 * This class contains the combination of a local plasticity model for concrete with a local isotropic damage model. This is an extension of concretedpm2. The yield surface of the plasticity model is based on the extension of the Menetrey and Willam yield criterion. The flow rule is nonassociated. The evolution laws of the hardening variables depend on the stress state. The plasticity model describes only hardening and perfect plasticity. It is based on the effective stress. The damage parameter of the isotropic damage model is based on the total volumetric strain. An exponential softening law is implemented.
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
    double ASoft, BSoft;

    /// Parameter of the hardening law of the plasticity model.
    double yieldHardPrimePeak;

    /// Parameter of the hardening law of the plasticity model.
    double yieldHardInitial;

    /// Control parameter for te volumetric plastic flow of the plastic potential.
    double dilationConst;

    /// Volumetric stress.
    double sig;

    /// Length of the deviatoric stress.
    double rho;

    /// Lode angle of the trial stress..
    double thetaTrial;

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
    double wf;

    /// Control parameter for the exponential softening law.
    double wfOne;

    /// Control parameter for the exponential softening law.
    double ftOne;

    /// yield tolerance for the plasticity model.
    double yieldTol;

    /// Maximum number of iterations for stress return.
    int newtonIter;

    /// Material mode for convenient access.
    MaterialMode matMode;

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

    virtual int hasMaterialModeCapability(MaterialMode mMode);

    virtual const char *giveClassName() const { return "ConcreteDPM2"; }
    virtual classType giveClassID() const { return ConcreteDPM2Class; }

    LinearElasticMaterial *giveLinearElasticMaterial()
    { return linearElasticMaterial; }

    int computeInverseOfJacobian(FloatMatrix &answer, const FloatMatrix &src);

    /**
     * Computes real macro stress in corresponding macro integration point for
     * given macroscopic strain vector.
     */
    virtual void giveRealStressVector(FloatArray &answer,
                                      MatResponseForm form,
                                      GaussPoint *gp,
                                      const FloatArray &strainVector,
                                      TimeStep *atTime);

    /** 
     * Perform stress return of the plasticity model and compute history variables.
     * @param gp Gauss point.
     * @param strain Strain vector of this Gauss point.
     */
    void performPlasticityReturn(GaussPoint *gp,
                                 StrainVector &strain);

    /**
     * Check if the trial stress state falls within the vertex region of the plasticity model at the apex of triaxial extension or triaxial compression.
     * @returns true for vertex case and false if regular stress return can be used.
     * @param answer Volumetric apex stress.
     * @param sig Volumetric stress.
     * @param tempKappa Hardening variable.
     */

    bool checkForVertexCase(double answer,
                            const double sig,
                            const double tempKappa);

    /**
     * Perform regular stress return for the plasticity model, i.e. if the trial stress state does not lie in the vertex region.
     * @param stress Stress vector which is computed.
     * @param gp Gauss point.
     */
    double performRegularReturn(StressVector &stress,
                                GaussPoint *gp);


    void computeJacobian(FloatMatrix &answer,
                         const double sig,
                         const double rho,
                         const double tempKappa,
                         const double deltaLambda,
                         GaussPoint *gp);

    /**
     * Perform stress return for vertex case of the plasticity model, i.e. if the trial stress state lies within the vertex region.
     * @param strain Strain vector of this Gauss point.
     * @param apexStress Volumetric stress at the apex of the yield surface.
     * @param gp Gauss point.
     */
    void performVertexReturn(StressVector &stress,
                             double apexStress,
                             GaussPoint *gp);

    /**
     * Compute the yield value based on stress and hardening variable.
     * @param sig Volumetric stress.
     * @param rho Length of the deviatoric stress.
     * @param theta Lode angle of the stress state.
     * @param tempKappa Hardening variable.
     * @returns Yield value.
     */
    double computeYieldValue(const double sig,
                             const double rho,
                             const double theta,
                             const double tempKappa) const;

    /**
     * Compute the value of the hardening function based on the hardening variable.
     * @param tempKappa Hardening variable.
     * @returns Value of the hardening function.
     */
    double computeHardeningOne(const double tempKappa) const;

    /** Compute the derivative of the hardening function based on the
     *  hardening parameter.
     *  @param tempKappa Hardening variable.
     *  @return Derivative of the hardening function.
     */
    double computeHardeningOnePrime(const double tempKappa) const;

    /**
     * Compute the value of the hardening function based on the hardening variable.
     * @param tempKappa Hardening variable.
     * @return Value of the hardening function.
     */
    double computeHardeningTwo(const double tempKappa) const;

    /**
     * Compute the derivative of the hardening function based on the
     * hardening parameter.
     * @param tempKappa Hardening variable.
     * @return Derivative of the hardening function.
     */
    double computeHardeningTwoPrime(const double tempKappa) const;

    /** Compute the derivative of the yield surface with respect to the hardening
     *  variable based on the stress state and the hardening variable
     *  @param sig Volumetric stress.
     *  @param rho Deviatoric length.
     *  @param tempKappa Hardening variable.
     *  @return Derivative of the yield surface.
     */
    double computeDFDKappa(const double sig,
                           const double rho,
                           const double tempKappa);


    /**
     * Compute the derivative of kappa with respect of delta lambda based on the stress state and the hardening variable.
     * @param sig Volumetric stress.
     * @param rho Length of the deviatoric stress.
     * @param tempKappa Hardening variable.
     * @return Derivative of kappa with respect to delta lambda.
     */
    double computeDKappaDDeltaLambda(const double sig,
                                     const double rho,
                                     const double tempKappa);




    /**
     * Compute the ductility measure based on the stress state.
     * @param sig Volumetric stress.
     * @param rho Length of the deviatoric strength.
     * @param theta Lode angle of stress state.
     * @return Ductility measure.
     */
    virtual double computeDuctilityMeasure(const double sig,
                                           const double rho,
					   const double theta);


    /// Compute the first derivative of the ductility measure with respect to the invariants sig and rho based on the stress state and the hardening parameter.
    void computeDDuctilityMeasureDInv(FloatArray &answer,
                                      const double sig,
                                      const double rho,
                                      const double tempKappa);

    /** 
     * Here, the first derivative of the plastic potential with respect
     * to the invariants sig and rho are computed
     */
    void computeDGDInv(FloatArray &answer,
                       const double sig,
                       const double rho,
                       const double tempKappa);

    /**
     * This function computes the ratio of the volumetric and deviatoric component
     * of the flow direction. It is used within the vertex return to check,
     * if the vertex return is admissible.
     */
    double computeRatioPotential(const double sig,
                                 const double tempKappa);

    /**
     * This function computes the rate factor which is used to take into account the strain rate dependence of the material.
     */
    double computeRateFactor(const double alpha,
                             const double timeFactor,
                             GaussPoint *gp,
                             TimeStep *deltaTime);

    /**
     * Here, the second derivative of the plastic potential with respect to the
     * invariants sig and rho are computed.
     */
    void computeDDGDDInv(FloatMatrix &answer,
                         const double sig,
                         const double rho,
                         const double tempKappa);

    /**
     * Here, the mixed derivative of the plastic potential with respect
     * to the invariants and the hardening parameter are determined.
     */
    void computeDDGDInvDKappa(FloatArray &answer,
                              const double sig,
                              const double rho,
                              const double tempKappa);

    /**
     * Computes the mixed derivative of the hardening parameter kappa with
     * respect to the plastic multiplier delta Lambda and the invariants sig
     * and rho.
     */
    void computeDDKappaDDeltaLambdaDInv(FloatArray &answer,
                                        const double sig,
                                        const double rho,
                                        const double tempKappa);

    /**
     * Computes the derivative of the evolution law of the hardening parameter kappa with respect to the hardening variable kappa.
     */
    double computeDDKappaDDeltaLambdaDKappa(const double sig,
                                            const double rho,
                                            const double tempKappa);



    /**
     * Computes the derivative of the yield surface with respect to the
     * invariants sig and rho.
     */
    void computeDFDInv(FloatArray &answer,
                       const double sig,
                       const double rho,
                       const double tempKappa) const;

    /**
     * Compute tempKappa.
     */
    double computeTempKappa(const double kappaInitial,
                            const double sigTrial,
                            const double rhoTrial,
                            const double sig);


    /**
     * Perform stress return for the damage model, i.e. if the trial stress state does not violate the plasticity surface.
     */
    void  computeDamage(FloatArray &answer, const StrainVector &strain, const double timeFactor, GaussPoint *gp, TimeStep *atTime, const double alpha);


    double computeAlpha(StressVector &effectiveStressTension, StressVector &effectiveStressCompression, StressVector &effectiveStress);

    /// Compute damage parameter.
    virtual double computeDamageParam(double equivStrain, double kappaOne, double kappaTwo, GaussPoint *gp);

    /// Compute equivalent strain value.
    double computeDeltaPlasticStrainNormTension(double tempKappaD, double kappaD, GaussPoint *gp);

    double computeDeltaPlasticStrainNormCompression(const double tempAlpha, double tempKappaD, double kappaD, GaussPoint *gp);

    virtual void computeEquivalentStrain(double &kappaD, const StrainVector &strain, GaussPoint *gp, TimeStep *atTime);

    /// Compute the ductility measure for the damage model.
    double computeDuctilityMeasureDamage(const StrainVector &strain, GaussPoint *gp);

    /**
     * Initialize the characteristic length, if damage is not yet activated
     * Set the increase factor for the strain rate dependence
     */
    void initDamaged(double kappa,
                     const StrainVector &strain,
                     GaussPoint *gp);


    /// Compute the trial coordinates.
    void computeTrialCoordinates(const StressVector &stress, double &sig, double &rho, double &theta);


    /// Assign state flag.
    void assignStateFlag(GaussPoint *gp);

    /// Computes the derivative of rho with respect to the stress.
    void computeDRhoDStress(FloatArray &answer,
                            const StressVector &stress) const;

    /// Computes the derivative of sig with respect to the stress.
    void computeDSigDStress(FloatArray &answer) const;

    /// Computes the seconfd derivative of rho with the respect to the stress.
    void computeDDRhoDDStress(FloatMatrix &answer,
                              const StressVector &stress) const;

    /// Computes the derivative of costheta with respect to the stress.
    void computeDCosThetaDStress(FloatArray &answer,
                                 const StressVector &stress) const;

    /// Compute the derivative of R with respect to costheta.
    double computeDRDCosTheta(const double theta, const double ecc) const;

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form,
                                               MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    void computeSecantStiffness(FloatMatrix &answer, MatResponseForm form, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);



    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual int giveIPValue(FloatArray &answer,
                            GaussPoint *gp,
                            InternalStateType type,
                            TimeStep *atTime);

    virtual int giveIPValueSize(InternalStateType type,
                                GaussPoint *gp);

    virtual int giveIntVarCompFullIndx(IntArray &answer,
                                       InternalStateType type,
                                       MaterialMode mmode);

    virtual InternalStateValueType giveIPValueType(InternalStateType type);

protected:
    MaterialStatus *CreateStatus(GaussPoint *gp) const;
};
} //end namespace oofem
#endif
