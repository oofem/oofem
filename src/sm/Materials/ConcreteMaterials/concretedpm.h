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

#ifndef ConcreteDPM_h
#define ConcreteDPM_h

#include "sm/Materials/structuralmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "cltypes.h"
#include "sm/Materials/structuralms.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "gausspoint.h"
#include "mathfem.h"

///@name Input fields for ConcreteDPM
//@{
#define _IFT_ConcreteDPM_Name "concretedpm"
#define _IFT_ConcreteDPM_fc "fc"
#define _IFT_ConcreteDPM_ft "ft"
#define _IFT_ConcreteDPM_ecc "ecc"
#define _IFT_ConcreteDPM_kinit "kinit"
#define _IFT_ConcreteDPM_ahard "ahard"
#define _IFT_ConcreteDPM_bhard "bhard"
#define _IFT_ConcreteDPM_chard "chard"
#define _IFT_ConcreteDPM_dhard "dhard"
#define _IFT_ConcreteDPM_asoft "asoft"
#define _IFT_ConcreteDPM_dilation "dilation"
#define _IFT_ConcreteDPM_yieldtol "yieldtol"
#define _IFT_ConcreteDPM_newtoniter "newtoniter"
#define _IFT_ConcreteDPM_wf "wf"
#define _IFT_ConcreteDPM_gf "gf"
#define _IFT_ConcreteDPM_href "href"
#define _IFT_ConcreteDPM_helem "helem"
//@}

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

    enum Concrete_VertexType {
        VT_Regular,
        VT_Tension,
        VT_Compression
    };

protected:
    /// @name History variables of the plasticity model
    //@{
    FloatArrayF< 6 >plasticStrain;
    FloatArrayF< 6 >tempPlasticStrain;

    double tempVolumetricPlasticStrain = 0.;

    double dFDKappa = 0.;
    double deltaLambda = 0.;
    //@}

    /// @name Hardening variable
    //@{
    double kappaP = 0.;
    double tempKappaP = 0.;
    //@}

    double le = 0.;

    /// @name History variables of the damage model
    //@{
    double equivStrain = 0.;
    double tempEquivStrain = 0.;

    double kappaD = 0.;
    double tempKappaD = 0.;

    double damage = 0.;
    double tempDamage = 0.;

    double deltaEquivStrain = 0.;
    //@}

    /// @name Indicates the state (i.e. elastic, unloading, plastic, damage, vertex) of the Gauss point
    //@{
    int state_flag = ConcreteDPM_Elastic;
    int temp_state_flag = ConcreteDPM_Elastic;
    //@}

    int vertexType = ConcreteDPMStatus::VT_Regular;
    int tempVertexType = ConcreteDPMStatus::VT_Regular;


#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
    ///  @name History variable of the modified size-dependent adjustment
    /// (indicating value of omega*ft/E+kappaD at the onset of localization)
    double epsloc = -1.;
    double tempEpsloc = -1.;
#endif

public:
    /// Constructor.
    ConcreteDPMStatus(GaussPoint *gp);

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;
    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    int setIPValue(const FloatArray &value, InternalStateType type);

    const char *giveClassName() const override { return "ConcreteDPMStatus"; }

    /**
     * Get the plastic strain deviator from the material status.
     * @return Plastic strain deviator.
     */
    const FloatArrayF< 6 > &givePlasticStrain() const { return plasticStrain; }

    /**
     * Get the deviatoric plastic strain norm from the material status.
     * @return Deviatoric plasticStrainNorm.
     */
    double giveDeviatoricPlasticStrainNorm()
    {
        auto dev = StructuralMaterial::computeDeviator(plasticStrain);
        return sqrt(.5 * ( 2. * dev [ 0 ] * dev [ 0 ] + 2. * dev [ 1 ] * dev [ 1 ] + 2. * dev [ 2 ] * dev [ 2 ] +
                           dev [ 3 ] * dev [ 3 ] + dev [ 4 ] * dev [ 4 ] + dev [ 5 ] * dev [ 5 ] ) );
    }

    double giveVolumetricPlasticStrain() const
    {
        return 1. / 3. * ( plasticStrain [ 0 ] + plasticStrain [ 1 ] + plasticStrain [ 2 ] );
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

    /**
     * History variable of the modified size-dependent adjustment
     * Assign the temp value of the damage variable of the damage model.
     * @param v New temp value of the damage variable.
     */
    void letTempEpslocBe(double v) { tempEpsloc = v; }

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
     * @return Temp value of plastic strain vector.
     */
    const FloatArrayF< 6 > &giveTempPlasticStrain() const { return tempPlasticStrain; }

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

    int giveTempVertexType() const { return tempVertexType; }

    /**
     * Assign the temp value of deviatoric plastic strain.
     * @param v New temp value of deviatoric plastic strain.
     */
    void letTempPlasticStrainBe(const FloatArrayF< 6 > &v) { tempPlasticStrain = v; }
    /**
     * Assigns converged value (only to be used when restoring consistency)
     */
    void letPlasticStrainBe(const FloatArrayF< 6 > &v) { plasticStrain = v; }

    /**
     * Assign the value of deviatoric plastic strain.
     * @param v New temp value of deviatoric plastic strain.
     */
    void letDeltaLambdaBe(double v) { deltaLambda = v; }

    /**
     *  Assign the temp value of the volumetric
     *  plastic strain in plane stress
     */
    void letTempVolumetricPlasticStrainBe(double v)
    { tempVolumetricPlasticStrain = v; }

    /**
     * Assign the temp value of the hardening variable of the plasticity model.
     * @param v New temp value of the hardening variable
     */
    void letTempKappaPBe(double v)
    { tempKappaP = v; }


    /**
     * Assign the temp value of the hardening variable of the damage model.
     * @param v New temp value of the hardening variable.
     */
    void letTempKappaDBe(double v) { tempKappaD = v; }
    /**
     * Assigns converged value (only to be used when restoring consistency)
     */
    void letKappaDBe(double v) { kappaD = v; }
    /**
     * Assign the temp value of the hardening variable of the damage model.
     * @param v New temp value of the hardening variable.
     */
    void letTempEquivStrainBe(double v) { tempEquivStrain = v; }
    /**
     * Assigns converged value (only to be used when restoring consistency)
     */
    void letEquivStrainBe(double v) { equivStrain = v; }

    /**
     * Assign the temp value of the damage variable of the damage model.
     * @param v New temp value of the damage variable.
     */
    void letTempDamageBe(double v) { tempDamage = v; }

    /**
     * Assign the temp value of the damage variable of the damage model.
     * @param v New temp value of the damage variable.
     */
    void letDeltaEquivStrainBe(double v) { deltaEquivStrain = v; }

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
    void letTempStateFlagBe(int v) { temp_state_flag = v; }

    void letTempVertexTypeBe(const int type) { tempVertexType = type; }
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
    //    enum Concrete_VertexType { VT_Regular, VT_Tension, VT_Compression };
    //  mutable Concrete_VertexType vertexType; /// FIXME: This must be removed, material models must never have mutable state.

    /// Linear elastic material
    IsotropicLinearElasticMaterial linearElasticMaterial;

    /**
     * Parameters of the yield surface of the plasticity model. fc is the uniaxial compressive strength, ft the uniaxial tensile strength and ecc controls the out of roundness of the deviatoric section.
     */
    double fc = 0., ft = 0., ecc = 0.;

    /// Parameter of the ductilityMeasure of the plasticity model.
    double AHard = 0.;
    /// Parameter of the ductilityMeasure of the plasticity model.
    double BHard = 0.;
    /// Parameter of the ductilityMeasure of the plasticity model.
    double CHard = 0.;
    /// Parameter of the ductilityMeasure of the plasticity model.
    double DHard = 0.;

    /// Parameter of the ductilityMeasure of the damage model.
    double ASoft = 0.;

    /// Parameter of the hardening law of the plasticity model.
    double yieldHardInitial = 0.;

    /// Control parameter for te volumetric plastic flow of the plastic potential.
    double dilationConst = 0.;

    /// Plastic multiplier of the plasticity model.
    //double deltaLambda = 0.;

    /// the volumetric stress.
    //double sig = 0.;
    /// The length of the deviatoric stress.
    //double rho = 0.;

    /// The lode angle of the trial stress.
    //    double thetaTrial = 0.;

    /// The friction parameter of the yield surface.
    double m = 0.;

    /// The dilation parameter of the plastic potential.
    double mQ = 0.;

    /// Element size (to be used in fracture energy approach (crack band).
    double helem = 0.;

    /// Elastic Young's modulus.
    double eM = 0.;
    /// Elastic shear modulus.
    double gM = 0.;
    /// Elastic bulk modulus.
    double kM = 0.;
    /// Elastic Poisson's ration.
    double nu = 0.;

    /// Hardening variable of plasticity model.
    //double kappaP = 0.;
    //double tempKappaP = 0.;

    /// Hardening variable of damage model.
    //double kappaD = 0.;
    //double tempKappaD = 0.;

    /// Damage variable of damage model.
    //double damage = 0.;
    //double tempDamage = 0.;

    /// Control parameter for the exponential softening law.
    double ef = 0.;

    /// Yield tolerance for the plasticity model.
    double yieldTol = 0.;

    /// Maximum number of iterations for stress return.
    int newtonIter = 0;

    /// Stress and its deviatoric part.
    //FloatArray effectiveStress;

#ifdef SOPHISTICATED_SIZEDEPENDENT_ADJUSTMENT
    /// Material parameter of the size-dependent adjustment
    /// (reference element size)
    double href = 0.;
#endif

public:
    /// Constructor
    ConcreteDPM(int n, Domain *d);

    void initializeFrom(InputRecord &ir) override;

    const char *giveClassName() const override { return "ConcreteDPM"; }
    const char *giveInputRecordName() const override { return _IFT_ConcreteDPM_Name; }

    ConcreteDPMStatus *giveConcreteDPMStatus(GaussPoint *gp) const
    { return static_cast< ConcreteDPMStatus * >( this->Material::giveStatus(gp) ); }

    FloatArrayF< 6 >giveRealStressVector_3d(const FloatArrayF< 6 > &strain, GaussPoint *gp,
                                            TimeStep *tStep) const override;

    /**
     * @param gp Gauss point.
     * @param strain Strain vector of this Gauss point.
     */
    void performPlasticityReturn(GaussPoint *gp, FloatArrayF< 6 > &strain) const;

    /**
     * Check if the trial stress state falls within the vertex region of the plasticity model at the apex of triaxial extension or triaxial compression.
     * @returns true for vertex case and false if regular stress return can be used.
     * @param answer The volumetric apex stress.
     * @param sig The volumetric stress.
     * @param tempKappa The hardening variable.
     */

    void checkForVertexCase(double &answer,
                            double sig,
                            double tempKappa,
                            GaussPoint *gp) const;

    /**
     * Perform regular stress return for the plasticity model, i.e. if the trial stress state does not lie in the vertex region.
     * @param stress Stress vector which is computed.
     * @param gp Gauss point.
     * @param theta Lode angle of the stress state.
     */
    void performRegularReturn(FloatArrayF< 6 > &stress,
                              GaussPoint *gp,
                              double theta) const;


    /**
     * Perform stress return for vertex case of the plasticity model, i.e. if the trial stress state lies within the vertex region.
     * @param stress Stress vector of this Gauss point.
     * @param apexStress Volumetric stress at the apex of the yield surface.
     * @param gp Gauss point.
     */
    void performVertexReturn(FloatArrayF< 6 > &stress,
                             double apexStress,
                             GaussPoint *gp) const;

    /**
     * Compute the yield value based on stress and hardening variable.
     * @param sig Volumetric stress.
     * @param rho Length of the deviatoric stress.
     * @param theta Lode angle of the stress state.
     * @param tempKappa Hardening variable.
     * @return Yield value.
     */
    double computeYieldValue(double sig,
                             double rho,
                             double theta,
                             double tempKappa) const;

    /**
     * Compute the value of the hardening function based on the hardening
     * variable.
     * @param tempKappa Hardening variable.
     * @return Value of the hardening function.
     */
    double computeHardeningOne(double tempKappa) const;

    /**
     * Compute the derivative of the hardening function based on the
     * hardening parameter.
     * @param tempKappa Hardening variable.
     * @return The derivative of the hardening function.
     */
    double computeHardeningOnePrime(double tempKappa) const;


    /**
     * Compute the derivative of the yield surface with respect to the hardening
     * variable based on the stress state and the hardening variable.
     * @param sig Volumetric stress.
     * @param rho Deviatoric length.
     * @param theta Lode angle of the stress state.
     * @param tempKappa Hardening variable.
     * @return The derivative of the yield surface.
     */
    double computeDFDKappa(double sig,
                           double rho,
                           double theta,
                           double tempKappa) const;

    /**
     * Compute the derivative of kappa with respect of delta
     * lambda based on the stress state and the hardening variable.
     * @param sig Volumetric stress.
     * @param rho Length of the deviatoric stress.
     * @param theta Lode angle of the stress state.
     * @param tempKappa Hardening variable.
     * @return Derivative of kappa with respect to delta lambda.
     */
    double computeDKappaDDeltaLambda(double sig,
                                     double rho,
                                     double theta,
                                     double tempKappa) const;

    /**
     * Compute the ductility measure based on the stress state.
     * @param sig Volumetric stress.
     * @param rho Length of the deviatoric strength.
     * @param theta Lode angle of stress state.
     * @returns Ductility measure.
     */
    virtual double computeDuctilityMeasure(double sig,
                                           double rho,
                                           double theta) const;


    /**
     * Computes the first derivative of the ductility measure with respect to
     * the invariants sig and rho based on the stress state and the hardening parameter.
     */
    FloatArrayF< 2 >computeDDuctilityMeasureDInv(double sig,
                                                 double rho,
                                                 double theta,
                                                 double tempKappa) const;

    /**
     * This matrix is the core of the closest point projection and collects
     * the derivative of the flow rule and the hardening parameters.
     */
    FloatMatrixF< 3, 3 >computeAMatrix(double sig,
                                       double rho,
                                       double theta,
                                       double tempKappa,
                                       double deltaLambda) const;

    /**
     * Here, the first derivative of the plastic potential with respect
     * to the invariants sig and rho are computed.
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
                                 double tempKappa) const;

    /**
     * Here, the second derivative of the plastic potential with respect to the
     * invariants sig and rho are computed.
     */
    FloatMatrixF< 2, 2 >computeDDGDDInv(double sig,
                                        double rho,
                                        double tempKappa) const;

    /**
     * Here, the mixed derivative of the plastic potential with respect
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
     * Computes the derivative of the evolution law of the hardening
     * parameter kappa with respect to the hardening variable kappa.
     */
    double computeDDKappaDDeltaLambdaDKappa(double sig,
                                            double rho,
                                            double theta,
                                            double tempKappa) const;


    /**
     * Computes the derivative of the yield surface with respect to the
     * invariants sig and rho.
     */
    FloatArrayF< 2 >computeDFDInv(double sig,
                                  double rho,
                                  double theta,
                                  double tempKappa) const;

    /**
     * Compute temporary kappa.
     */
    double computeTempKappa(double kappaInitial,
                            double sigTrial,
                            double rhoTrial,
                            double sig) const;

    /**
     * Perform stress return for the damage model, i.e. if the trial stress state does not violate the plasticity surface.
     * @param strain Strain.
     * @param gp Gauss point.
     * @param tStep Time step.
     * @return Damage.
     */
    std::pair< double, double >computeDamage(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) const;

    /// Compute damage parameter.
    virtual double computeDamageParam(double kappa, GaussPoint *gp) const;

    /// Compute the damage-driving variable from given damage.
    double computeInverseDamage(double dam, GaussPoint *gp) const;

    /// Compute equivalent strain value.
    virtual double computeEquivalentStrain(const FloatArrayF< 6 > &elasticStrain, GaussPoint *gp, TimeStep *tStep) const;

    /// Compute the ductility measure for the damage model.
    double computeDuctilityMeasureDamage(const FloatArrayF< 6 > &strain, GaussPoint *gp) const;

    /**
     * Initialize the characteristic length, if damage is not yet activated.
     */
    void initDamaged(double kappa,
                     const FloatArrayF< 6 > &elasticStrain,
                     GaussPoint *gp) const;


    /// Compute the trial coordinates.
    std::tuple< double, double, double >computeTrialCoordinates(const FloatArrayF< 6 > &stress, GaussPoint *gp) const;

    /// Assign state flag.
    void assignStateFlag(GaussPoint *gp) const;

    /// Computes the derivative of rho with respect to the stress.
    FloatArrayF< 6 >computeDRhoDStress(const FloatArrayF< 6 > &stress) const;

    /// Computes the derivative of sig with respect to the stress.
    void computeDSigDStress(FloatArrayF< 6 > &answer) const;

    /// Computes the second derivative of rho with the respect to the stress.
    FloatMatrixF< 6, 6 >computeDDRhoDDStress(const FloatArrayF< 6 > &stress) const;

    /// Computes the derivative of costheta with respect to the stress.
    FloatArrayF< 6 >computeDCosThetaDStress(const FloatArrayF< 6 > &stress) const;

    /// Compute the derivative of R with respect to costheta.
    double computeDRDCosTheta(double theta, double ecc) const;

    FloatMatrixF< 6, 6 >give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    int setIPValue(const FloatArray &value, GaussPoint *gp, InternalStateType type) override;

    int giveIPValue(FloatArray &answer,
                    GaussPoint *gp,
                    InternalStateType type,
                    TimeStep *tStep) override;

    void restoreConsistency(GaussPoint *gp) override;
    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;


protected:
    MaterialStatus *CreateStatus(GaussPoint *gp) const override;
};
} // end namespace oofem
#endif
