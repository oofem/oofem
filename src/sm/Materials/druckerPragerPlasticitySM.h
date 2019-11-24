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

#ifndef druckerpragerplasticitysm_h
#define druckerpragerplasticitysm_h

#include "floatarray.h"
#include "floatmatrix.h"

#include "sm/Materials/structuralms.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/isolinearelasticmaterial.h"

///@name Input fields for DruckerPragerPlasticitySM
//@{
#define _IFT_DruckerPragerPlasticitySM_Name "druckerprager"
#define _IFT_DruckerPragerPlasticitySM_iys "iys" ///< Initial yield stress under pure shear
#define _IFT_DruckerPragerPlasticitySM_alpha "alpha" ///< Friction coefficient
#define _IFT_DruckerPragerPlasticitySM_alphapsi "alphapsi" ///< Dilatancy coefficient
#define _IFT_DruckerPragerPlasticitySM_ht "ht"
#define _IFT_DruckerPragerPlasticitySM_hm "hm"
#define _IFT_DruckerPragerPlasticitySM_kc "kc"
#define _IFT_DruckerPragerPlasticitySM_lys "lys"
#define _IFT_DruckerPragerPlasticitySM_yieldtol "yieldtol"
#define _IFT_DruckerPragerPlasticitySM_newtoniter "newtoniter"
//@}

namespace oofem {
/**
 * This class implements the material status associated to DruckerPragerPlasticitySM.
 * Tracks volumetric and deviatoric plastic strain and hardening.
 * @author Simon Rolshoven
 */
class DruckerPragerPlasticitySMStatus : public StructuralMaterialStatus
{
public:
    /// Values of history variable state_flag.
    enum state_flag_values { DP_Elastic, DP_Unloading, DP_Yielding, DP_Vertex };

protected:
    /// Volumetric plastic strain.
    double volumetricPlasticStrain = 0.;
    double tempVolumetricPlasticStrain = 0.;

    /// Deviatoric of plastic strain.
    FloatArrayF<6> plasticStrainDeviator;
    FloatArrayF<6> tempPlasticStrainDeviator;

    /// Hardening variable.
    double kappa = 0.;
    double tempKappa = 0.;

    /// Indicates the state (i.e. elastic, yielding, vertex, unloading) of the Gauss point
    int state_flag = DruckerPragerPlasticitySMStatus :: DP_Elastic;
    int temp_state_flag = DruckerPragerPlasticitySMStatus :: DP_Elastic;

public:
    /// Constructor
    DruckerPragerPlasticitySMStatus(GaussPoint * gp);

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;
    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    const char *giveClassName() const override { return "DruckerPragerPlasticitySMStatus"; }

    /**
     * Get the full plastic strain vector from the material status.
     * @param answer Plastic strain vector.
     */
    FloatArrayF<6> givePlasticStrainVector() const
    {
        auto plasticStrain = plasticStrainDeviator;
        plasticStrain[0] += volumetricPlasticStrain;
        plasticStrain[1] += volumetricPlasticStrain;
        plasticStrain[2] += volumetricPlasticStrain;
        return plasticStrain;
    }
    /**
     * Get the plastic strain deviator from the material status.
     * @return Plastic strain deviator.
     */
    const FloatArrayF<6> &givePlasticStrainDeviator() const { return plasticStrainDeviator; }
    /**
     * Get the volumetric plastic strain from the material status.
     * @return Volumetric plastic strain.
     */
    double giveVolumetricPlasticStrain() const { return volumetricPlasticStrain; }
    /**
     * Get the hardening variable from the material status.
     * @return hardening variable kappa
     */
    double giveKappa() const { return kappa; }
    /**
     * Get the state flag from the material status.
     * @return State flag (i.e. elastic, unloading, yielding, vertex case yielding)
     */
    int giveStateFlag() const { return state_flag; }

    /**
     * Get the temp value of the full plastic strain vector from the material status.
     * @param answer Temp value of plastic strain vector.
     */
    FloatArrayF<6> giveTempPlasticStrainVector() const
    {
        auto plasticStrain = tempPlasticStrainDeviator;
        plasticStrain[0] += tempVolumetricPlasticStrain;
        plasticStrain[1] += tempVolumetricPlasticStrain;
        plasticStrain[2] += tempVolumetricPlasticStrain;
        return plasticStrain;
    }
    /**
     * Get the temp value of the plastic strain deviator from the material status.
     * @param answer Temp value of plastic strain deviator.
     */
    const FloatArrayF<6> &giveTempPlasticStrainDeviator() const { return tempPlasticStrainDeviator; }
    /**
     * Get the temp value of the volumetric strain deviator from the material status.
     * @return Temp value of volumetric plastic strain
     */
    double giveTempVolumetricPlasticStrain() const { return tempVolumetricPlasticStrain; }
    /**
     * Get the temp value of the hardening variable from the material status.
     * @return Temp value of hardening variable kappa.
     */
    double giveTempKappa() const { return tempKappa; }
    /**
     * Get the temp value of the state flag from the material status.
     * @return Temp value of the state flag (i.e. elastic, unloading, yielding, vertex case yielding)
     */
    int giveTempStateFlag() const { return temp_state_flag; }

    /**
     * Assign the temp value of deviatoric plastic strain.
     * @param v New temp value of deviatoric plastic strain.
     */
    void letTempPlasticStrainDeviatorBe(const FloatArrayF<6> &v) { tempPlasticStrainDeviator = v; }
    /**
     * Assign the temp value of volumetric plastic strain.
     * @param v New temp value of volumetric plastic strain.
     */
    void letTempVolumetricPlasticStrainBe(double v) { tempVolumetricPlasticStrain = v; }
    /**
     * Assign the temp value of the hardening variable.
     * @param v New temp value of the hardening variable.
     */
    void letTempKappaBe(double v) { tempKappa = v; }
    /**
     * Assign the temp value of the state flag.
     * @param v New temp value of the state flag (i.e. elastic, unloading, yielding, vertex case yielding).
     */
    void letTempStateFlagBe(int v) { temp_state_flag = v; }
};

/**
 * This class implements a (local) nonassociated plasticity model based on the Drucker-Prager yield criterion with hardening and softening.
 * @author Simon Rolshoven
 */
class DruckerPragerPlasticitySM : public StructuralMaterial
{
public:
    /// Values of history variable state_flag.
    enum state_flag_values { DP_Elastic, DP_Unloading, DP_Yielding, DP_Vertex };

protected:
    /**
     * Controls the hardening function in the yield stress:
     * 1: linear hardening/softening with cutoff at zero stress.
     * 2: exponential hardening/softening to limitYieldStress.
     */
    int hardeningType = 1;
    /// Parameter of the exponential laws.
    double kappaC = 0.;
    /// Hardening modulus normalized with the elastic modulus, parameter of the linear hardening/softening law.
    double hardeningModulus = 0.;
    /// Parameter of the exponential hardening law.
    double limitYieldStress = 0.;
    /// Parameter of all three laws, this is the initial value of the yield stress in pure shear.
    double initialYieldStress = 0.;
    /// Friction coefficient, parameter of the yield criterion.
    double alpha = 0.;
    /// Dilatancy coefficient, parameter of the flow rule.
    double alphaPsi = 0.;
    /// Scalar factor between rate of plastic multiplier and rate of hardening variable.
    double kFactor = 0.;

    /// Associated linear elastic material.
    IsotropicLinearElasticMaterial LEMaterial;

    /// Yield tolerance.
    double yieldTol = 0.;
    /// Maximum number of iterations for stress return.
    int newtonIter = 0;

public:
    /// Constructor
    DruckerPragerPlasticitySM(int n, Domain * d);

    void initializeFrom(InputRecord &ir) override;

    const char *giveClassName() const override { return "DruckerPragerPlasticitySM"; }
    const char *giveInputRecordName() const override { return _IFT_DruckerPragerPlasticitySM_Name; }

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp,
                                           TimeStep *tStep) const override;

    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;

    /**
     * Perform a standard local stress return using the function computeYieldValue at the specified Gauss point.
     * This function computes the history variables, i.e. the temp variable of plastic strain, hardening variable, state flag, and the temp stress, and stores them in the temp-status.
     * @param gp Gauss point.
     * @param strain Strain vector of this Gauss point.
     */
    void performLocalStressReturn(GaussPoint *gp, const FloatArrayF<6> &strain) const;
    /**
     * Check if the trial stress state falls within the vertex region.
     * @param eM Elasticity modulus.
     * @param gM Shear modulus.
     * @param kM Bulk modulus.
     * @return True for vertex case and false if regular stress return has to be used.
     */
    bool checkForVertexCase(double eM, double gM, double kM, double trialStressJTwo, double volumetricStress, double tempKappa) const;
    /**
     * Perform stress return for regular case, i.e. if the trial stress state does not lie within the vertex region.
     * @param eM Elasticity modulus.
     * @param gM Shear modulus.
     * @param kM Bulk modulus.
     */
    void performRegularReturn(double eM, double gM, double kM, double trialStressJTwo, FloatArrayF<6> &stressDeviator, double &volumetricStress, double &tempKappa) const;
    /**
     * Perform stress return for vertex case, i.e. if the trial stress state lies within the vertex region.
     * @param eM Elasticity modulus.
     * @param gM Shear modulus.
     * @param kM Bulk modulus.
     */
    void performVertexReturn(double eM, double gM, double kM, double trialStressJTwo, FloatArrayF<6> &stressDeviator, double &volumetricStress, double &tempKappa, double volumetricElasticTrialStrain, double kappa) const;
    /**
     * Compute the yield value based on stress and hardening variable.
     * @param meanStress 1/3 of trace of sigma.
     * @param JTwo Second deviatoric invariant.
     * @param kappa Hardening variable.
     * @param eM Elasticity modulus.
     * @return Yield value.
     */
    double computeYieldValue(double meanStress,
                             double JTwo,
                             double kappa,
                             double eM) const;
    /**
     * Compute the current yield stress in pure shear of the Drucker-Prager model according to the used hardening law. The yield stress is tauY in f(sigma, kappa) = F(sigma) - tauY(kappa).
     * @param kappa Hardening variable.
     * @param eM Elasticity modulus.
     * @returns Yield stress in pure shear.
     */
    virtual double computeYieldStressInShear(double kappa, double eM) const;

    /**
     * Compute derivative of yield stress with respect to the hardening variable kappa.
     * @param kappa Hardening variable.
     * @param eM Elasticity modulus.
     * @return Derivative of yield stress with respect to kappa.
     */
    virtual double computeYieldStressPrime(double kappa, double eM) const;

    /**
     * Compute and give back algorithmic stiffness matrix for the regular case (no vertex).
     * @param mode Material reponse mode.
     * @param gp Gauss point.
     * @param tStep Time step.
     * @return Consistent stiffness matrix.
     */
    FloatMatrixF<6,6> giveRegAlgorithmicStiffMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const;
    /**
     * Compute consistent stiffness matrix for the vertex case.
     * @param mode Material reponse mode.
     * @param gp Gauss point.
     * @param tStep Time step.
     * @return Consistent stiffness matrix.
     */
    FloatMatrixF<6,6> giveVertexAlgorithmicStiffMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const;

    int giveIPValue(FloatArray &answer,
                    GaussPoint *gp,
                    InternalStateType type,
                    TimeStep *tStep) override;

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return false; }

    FloatArrayF<6> giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const override
    {
        return LEMaterial.giveThermalDilatationVector(gp, tStep);
    }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    double predictRelativeComputationalCost(GaussPoint *gp) override;
    double predictRelativeRedistributionCost(GaussPoint *gp) override { return 1.0; }
};
} // end namespace oofem
#endif // druckerpragerplasticitysm_h
