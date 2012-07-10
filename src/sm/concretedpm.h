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

#ifndef ConcreteDPM_h
#define ConcreteDPM_h

#include "structuralmaterial.h"
#include "gausspnt.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"
#include "structuralms.h"
#include "strainvector.h"
#include "stressvector.h"
#include "isolinearelasticmaterial.h"

namespace oofem {
/**
 * This class implements the material status associated to ConcreteDPM.
 * @author Peter Grassl
 */
// new approach to size-dependent adjustment of damage evolution
// (can be deactivated by commenting out the following line)
#define SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT

class ConcreteDPMStatus : public StructuralMaterialStatus
{
public:
    /// Values of history variable state_flag
    enum state_flag_values {
        ConcreteDPM_Elastic,
        ConcreteDPM_Unloading,
        ConcreteDPM_Plastic,
        ConcreteDPM_Damage,
        ConcreteDPM_PlasticDamage,
        ConcreteDPM_VertexCompression,
        ConcreteDPM_VertexTension,
        ConcreteDPM_VertexCompressionDamage,
        ConcreteDPM_VertexTensionDamage
    };

protected:
    /// @name History variables of the plasticity model
    //@{
    StrainVector plasticStrain;
    StrainVector tempPlasticStrain;

    double tempVolumetricPlasticStrain;

    double dFDKappa;
    double deltaLambda;
    //@}

    /// @name Hardening variable
    //@{
    double kappaP;
    double tempKappaP;
    //@}

    double le;

    /// @name History variables of the damage model
    //@{
    double equivStrain;
    double tempEquivStrain;

    double kappaD;
    double tempKappaD;

    double damage;
    double tempDamage;

    double deltaEquivStrain;
    //@}

    /// @name Indicates the state (i.e. elastic, unloading, plastic, damage, vertex) of the Gauss point
    //@{
    int state_flag;
    int temp_state_flag;
    //@}

#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
    ///  @name History variable of the modified size-dependent adjustment
    /// (indicating value of omega*ft/E+kappaD at the onset of localization)
    double epsloc;
#endif

public:
    /// Constructor.
    ConcreteDPMStatus(int n, Domain *d, GaussPoint *gp);

    /// Destructor.
    virtual ~ConcreteDPMStatus();

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    int setIPValue(const FloatArray value, InternalStateType type);

    virtual void restoreConsistency();
    virtual const char *giveClassName() const { return "ConcreteDPMStatus"; }
    virtual classType giveClassID() const { return ConcreteDPMStatusClass; }

    /**
     * Get the full plastic strain vector from the material status.
     * @param answer Plastic strain vector.
     */
    void giveFullPlasticStrainVector(StrainVector &answer) const
    {
        StrainVector plasticStrainVector(_Unknown);
        givePlasticStrain(plasticStrainVector);
        plasticStrainVector.convertToFullForm(answer);
    }
    /**
     * Get the plastic strain deviator from the material status.
     * @param answer Plastic strain deviator.
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
    double giveKappaP() const { return kappaP; }

    /**
     * Get the hardening variable of the damage model from the
     * material status.
     * @return Hardening variable kappaD.
     */
    double giveKappaD() const { return kappaD; }

    /**
     * Get the equivalent strain from the
     * material status.
     * @return Equivalent strain.
     */
    double giveEquivStrain() const { return equivStrain; }

#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
    /**
     * Get the value of omega*ft/E at the expected onset of localization
     * (defined by negative second-order work).
     * @returns Variable epsloc.
     */
    double giveEpsLoc() const { return epsloc; }
#endif

    /**
     * Get the damage variable of the damage model from the
     * material status.
     * @return Damage variable damage.
     */
    double giveDamage() const { return damage; }

    /**
     * Get the state flag from the material status.
     * @return State flag (i.e. elastic, unloading, yielding, vertex case yielding)
     */
    int giveStateFlag() const { return state_flag; }

    /**
     * Get the temp value of the full plastic strain vector from the material status.
     * @param answer Temp value of plastic strain vector.
     */
    void giveTempPlasticStrain(StrainVector &answer) const { answer = tempPlasticStrain; }

    /**
     * Get the temp value of the volumetric plastic strain in plane stress.
     * @return Temp value of volumetric plastic strain.
     */
    double giveTempVolumetricPlasticStrain() const { return tempVolumetricPlasticStrain; }


    /**
     * Get the temp value of the hardening variable of the plasticity model
     * from the material status.
     * @return Temp value of hardening variable kappaP.
     */
    double giveTempKappaP() const { return tempKappaP; }

    /**
     * Get the temp value of the hardening variable of the damage model
     * from the material status.
     * @return Temp value of hardening variable kappaD.
     */
    double giveTempKappaD() const { return tempKappaD; }

    /**
     * Get the temp value of the hardening variable of the damage model
     * from the material status.
     * @return Temp value of the damage variable damage.
     */
    double giveTempDamage() const { return tempDamage; }

    /**
     * Get the temp value of the hardening variable of the damage model
     * from the material status.
     * @return Temp value of the damage variable damage.
     */
    double giveDeltaEquivStrain() const { return deltaEquivStrain; }


    /**
     * Get the temp value of the state flag from the material status.
     * @return The temp value of the state flag (i.e. elastic, unloading,
     * yielding, vertex case yielding).
     */
    int giveTempStateFlag() const { return temp_state_flag; }

    /**
     * Assign the temp value of deviatoric plastic strain.
     * @param v New temp value of deviatoric plastic strain.
     */
    void letTempPlasticStrainBe(const StrainVector &v) { tempPlasticStrain = v; }

    /**
     * Assign the value of deviatoric plastic strain.
     * @param v New temp value of deviatoric plastic strain.
     */
    void letDeltaLambdaBe(const double v) { deltaLambda = v; }

    /**
     *  Assign the temp value of the volumetric
     *  plastic strain in plane stress
     */
    void  letTempVolumetricPlasticStrainBe(const double v)
    { tempVolumetricPlasticStrain = v; }

    /**
     * Assign the temp value of the hardening variable of the plasticity model.
     * @param v New temp value of the hardening variable
     */
    void  letTempKappaPBe(const double v)
    { tempKappaP = v; }


    /**
     * Assign the temp value of the hardening variable of the damage model.
     * @param v New temp value of the hardening variable.
     */
    void letTempKappaDBe(const double v) { tempKappaD = v; }

    /**
     * Assign the temp value of the hardening variable of the damage model.
     * @param v New temp value of the hardening variable.
     */
    void letTempEquivStrainBe(const double v) { tempEquivStrain = v; }

    /**
     * Assign the temp value of the damage variable of the damage model.
     * @param v New temp value of the damage variable.
     */
    void letTempDamageBe(const double v) { tempDamage = v; }

    /**
     * Assign the temp value of the damage variable of the damage model.
     * @param v New temp value of the damage variable.
     */
    void letDeltaEquivStrainBe(const double v) { deltaEquivStrain = v; }

    /**
     * Gives the characteristic length.
     */
    double giveLe() { return le; }

    /**
     * Sets the characteristic length.
     * @param ls New characteristic length.
     */
    void setLe(double ls) { le = ls; }

    /**
     * Assign the temp value of the state flag.
     * @param v New temp value of the state flag (i.e. elastic, unloading, yielding,
     * vertex case yielding).
     */
    void letTempStateFlagBe(const int v) { temp_state_flag = v; }
};


/**
 * This class contains the combination of a local plasticity model for concrete with a local isotropic damage model.
 * The yield surface of the plasticity model is based on the extension of the Menetrey and Willam yield criterion.
 * The flow rule is nonassociated. The evolution laws of the hardening variables depend on the stress state.
 * The plasticity model describes only hardening and perfect plasticity. It is based on h effective stress.
 * The damage parameter of the isotropic damage model is based on the total volumetric strain.
 * An exponential softening law is implemented.
 *
 * @author Peter Grassl
 */
class ConcreteDPM : public StructuralMaterial
{
protected:
    enum Concrete_VertexType { VT_Regular, VT_Tension, VT_Compression };
    Concrete_VertexType vertexType;

    /**
     * Parameters of the yield surface of the plasticity model. fc is the uniaxial compressive strength, ft the uniaxial tensile strength and ecc controls the out of roundness of the deviatoric section.
     */
    double fc, ft, ecc;

    /// Parameter of the ductilityMeasure of the plasticity model.
    double AHard;
    /// Parameter of the ductilityMeasure of the plasticity model.
    double BHard;
    /// Parameter of the ductilityMeasure of the plasticity model.
    double CHard;
    /// Parameter of the ductilityMeasure of the plasticity model.
    double DHard;

    /// Parameter of the ductilityMeasure of the damage model.
    double ASoft;

    /// Parameter of the hardening law of the plasticity model.
    double yieldHardInitial;

    /// Control parameter for te volumetric plastic flow of the plastic potential.
    double dilationConst;

    /// Plastic multiplier of the plasticity model.
    double deltaLambda;

    /// the volumetric stress.
    double sig;
    /// The length of the deviatoric stress.
    double rho;

    /// The lode angle of the trial stress.
    double thetaTrial;

    /// The friction parameter of the yield surface.
    double m;

    /// The dilation parameter of the plastic potential.
    double mQ;

    /// Element size (to be used in fracture energy approach (crack band).
    double helem;

    /// Pointer for linear elastic material
    LinearElasticMaterial *linearElasticMaterial;

    /// Elastic Young's modulus.
    double eM;
    /// Elastic shear modulus.
    double gM;
    /// Elastic bulk modulus.
    double kM;
    /// Elastic Poisson's ration.
    double nu;

    /// Hardening variable of plasticity model.
    double kappaP;
    double tempKappaP;

    /// Hardening variable of damage model.
    double kappaD;
    double tempKappaD;

    /// Damage variable of damage model.
    double damage;
    double tempDamage;

    /// Control parameter for the exponential softening law.
    double ef;

    /// Yield tolerance for the plasticity model.
    double yieldTol;

    /// Maximum number of iterations for stress return.
    int newtonIter;

    /// Material mode for convenient access.
    MaterialMode matMode;

    /// Stress and its deviatoric part.
    StressVector effectiveStress;

#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
    /// Material parameter of the size-dependent adjustment
    /// (reference element size)
    double href;
#endif

public:
    /// Constructor
    ConcreteDPM(int n, Domain *d);
    /// Destructor
    virtual ~ConcreteDPM();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int hasMaterialModeCapability(MaterialMode mMode);

    virtual const char *giveClassName() const { return "ConcreteDPM"; }
    virtual classType giveClassID() const { return ConcreteDPMClass; }

    virtual ConcreteDPMStatus *giveStatus(GaussPoint *gp) const
    { return static_cast< ConcreteDPMStatus * >( this->Material :: giveStatus(gp) ); }

    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    virtual void giveRealStressVector(FloatArray &answer,
                                      MatResponseForm form,
                                      GaussPoint *gp,
                                      const FloatArray &reducedStrain,
                                      TimeStep *tStep);


    /**
     * @param gp Gauss point.
     * @param strain Strain vector of this Gauss point.
     */
    void performPlasticityReturn(GaussPoint *gp, StrainVector &strain);

    /**
     * Check if the trial stress state falls within the vertex region of the plasticity model at the apex of triaxial extension or triaxial compression.
     * @returns true for vertex case and false if regular stress return can be used.
     * @param answer The volumetric apex stress.
     * @param sig The volumetric stress.
     * @param tempKappa The hardening variable.
     */

    bool checkForVertexCase(double answer,
                            const double sig,
                            const double tempKappa);

    /**
     * Perform regular stress return for the plasticity model, i.e. if the trial stress state does not lie in the vertex region.
     * @param stress Stress vector which is computed.
     * @param gp Gauss point.
     */
    void performRegularReturn(StressVector &stress,
                              GaussPoint *gp);


    /**
     * Perform stress return for vertex case of the plasticity model, i.e. if the trial stress state lies within the vertex region.
     * @param stress Stress vector of this Gauss point.
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
     * @return Yield value.
     */
    double computeYieldValue(const double sig,
                             const double rho,
                             const double theta,
                             const double tempKappa) const;

    /**
     * Compute the value of the hardening function based on the hardening
     * variable.
     * @param tempKappa Hardening variable.
     * @return Value of the hardening function.
     */
    double computeHardeningOne(const double tempKappa) const;

    /**
     * Compute the derivative of the hardening function based on the
     * hardening parameter.
     * @param tempKappa Hardening variable.
     * @return The derivative of the hardening function.
     */
    double computeHardeningOnePrime(const double tempKappa) const;


    /**
     * Compute the derivative of the yield surface with respect to the hardening
     * variable based on the stress state and the hardening variable.
     * @param sig Volumetric stress.
     * @param rho Deviatoric length.
     * @param tempKappa Hardening variable.
     * @return The derivative of the yield surface.
     */
    double computeDFDKappa(const double sig,
                           const double rho,
                           const double tempKappa);

    /**
     * Compute the derivative of kappa with respect of delta
     * lambda based on the stress state and the hardening variable.
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
     * @returns Ductility measure.
     */
    virtual double computeDuctilityMeasure(const double sig,
                                           const double rho,
					   const double theta);


    /**
     * Computes the first derivative of the ductility measure with respect to
     * the invariants sig and rho based on the stress state and the hardening parameter.
     */
    void computeDDuctilityMeasureDInv(FloatArray &answer,
                                      const double sig,
                                      const double rho,
                                      const double tempKappa);

    /**
     * This matrix is the core of the closest point projection and collects
     * the derivative of the flow rule and the hardening parameters.
     */
    void computeAMatrix(FloatMatrix &answer,
                        const double sig,
                        const double rho,
                        const double tempKappa);

    /**
     * Here, the first derivative of the plastic potential with respect
     * to the invariants sig and rho are computed.
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
     * Computes the derivative of the evolution law of the hardening
     * parameter kappa with respect to the hardening variable kappa.
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
     * Compute temporary kappa.
     */
    double computeTempKappa(const double kappaInitial,
                            const double sigTrial,
                            const double rhoTrial,
                            const double sig);


    /**
     * Perform stress return for the damage model, i.e. if the trial stress state does not violate the plasticity surface.
     * @param strain Strain.
     * @param gp Gauss point.
     * @param tStep Time step.
     * @return Damage.
     */
    double computeDamage(const StrainVector &strain, GaussPoint *gp, TimeStep *tStep);



    /// Compute damage parameter.
    virtual double computeDamageParam(double kappa, GaussPoint *gp);

    /// Compute the damage-driving variable from given damage.
    double computeInverseDamage(double dam, GaussPoint *gp);

    /// Compute equivalent strain value.
    virtual void computeEquivalentStrain(double &kappaD, const StrainVector &elasticStrain, GaussPoint *gp, TimeStep *tStep);


    /// Compute the ductility measure for the damage model.
    double computeDuctilityMeasureDamage(const StrainVector &strain, GaussPoint *gp);

    /**
     * Initialize the characteristic length, if damage is not yet activated.
     */
    void initDamaged(double kappa,
                     const StrainVector &elasticStrain,
                     GaussPoint *gp);


    /// Compute the trial coordinates.
    void computeTrialCoordinates(const StressVector &stress);

    /// Assign state flag.
    void assignStateFlag(GaussPoint *gp);

    /// Computes the derivative of rho with respect to the stress.
    void computeDRhoDStress(FloatArray &answer, const StressVector &stress) const;

    /// Computes the derivative of sig with respect to the stress.
    void computeDSigDStress(FloatArray &answer) const;

    /// Computes the second derivative of rho with the respect to the stress.
    void computeDDRhoDDStress(FloatMatrix &answer, const StressVector &stress) const;

    /// Computes the derivative of costheta with respect to the stress.
    void computeDCosThetaDStress(FloatArray &answer, const StressVector &stress) const;

    /// Compute the derivative of R with respect to costheta.
    double computeDRDCosTheta(const double theta, const double ecc) const;

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer, MatResponseForm form,
                                               MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);


    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual int setIPValue(const FloatArray value, GaussPoint *aGaussPoint, InternalStateType type);

    virtual int giveIPValue(FloatArray &answer,
                            GaussPoint *gp,
                            InternalStateType type,
                            TimeStep *tStep);

    virtual int giveIPValueSize(InternalStateType type,
                                GaussPoint *gp);

    virtual int giveIntVarCompFullIndx(IntArray &answer,
                                       InternalStateType type,
                                       MaterialMode mmode);

    virtual InternalStateValueType giveIPValueType(InternalStateType type);

protected:
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;
};
} // end namespace oofem
#endif
