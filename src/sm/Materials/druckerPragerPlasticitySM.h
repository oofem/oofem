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

#include "../sm/Materials/structuralms.h"
#include "../sm/Materials/structuralmaterial.h"
#include "Materials/isolinearelasticmaterial.h"

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
    double volumetricPlasticStrain;
    double tempVolumetricPlasticStrain;

    /// Deviatoric of plastic strain.
    FloatArray plasticStrainDeviator;
    FloatArray tempPlasticStrainDeviator;

    /// Hardening variable.
    double kappa;
    double tempKappa;

    /// Indicates the state (i.e. elastic, yielding, vertex, unloading) of the Gauss point
    int state_flag;
    int temp_state_flag;

public:
    /// Constructor
    DruckerPragerPlasticitySMStatus(int n, Domain * d, GaussPoint * gp);

    /// Destructor
    virtual ~DruckerPragerPlasticitySMStatus();

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    virtual const char *giveClassName() const { return "DruckerPragerPlasticitySMStatus"; }

    /**
     * Get the full plastic strain vector from the material status.
     * @param answer Plastic strain vector.
     */
    void  givePlasticStrainVector(FloatArray &answer) const
    {
        answer = plasticStrainDeviator;
        answer[0] += volumetricPlasticStrain;
        answer[1] += volumetricPlasticStrain;
        answer[2] += volumetricPlasticStrain;
    }
    /**
     * Get the plastic strain deviator from the material status.
     * @return Plastic strain deviator.
     */
    const FloatArray &givePlasticStrainDeviator() const { return plasticStrainDeviator; }
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
    void giveTempPlasticStrainVector(FloatArray &answer) const
    {
        answer = tempPlasticStrainDeviator;
        answer[0] += tempVolumetricPlasticStrain;
        answer[1] += tempVolumetricPlasticStrain;
        answer[2] += tempVolumetricPlasticStrain;
    }
    /**
     * Get the temp value of the plastic strain deviator from the material status.
     * @param answer Temp value of plastic strain deviator.
     */
    const FloatArray &giveTempPlasticStrainDeviator() const { return tempPlasticStrainDeviator; }
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
    void letTempPlasticStrainDeviatorBe(const FloatArray &v) { tempPlasticStrainDeviator = v; }
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
    int hardeningType;
    /// Parameter of the exponential laws.
    double kappaC;
    /// Hardening modulus normalized with the elastic modulus, parameter of the linear hardening/softening law.
    double hardeningModulus;
    /// Parameter of the exponential hardening law.
    double limitYieldStress;
    /// Parameter of all three laws, this is the initial value of the yield stress in pure shear.
    double initialYieldStress;
    /// Friction coefficient, parameter of the yield criterion.
    double alpha;
    /// Dilatancy coefficient, parameter of the flow rule.
    double alphaPsi;
    /// Scalar factor between rate of plastic multiplier and rate of hardening variable.
    double kFactor;

    /// Associated linear elastic material.
    IsotropicLinearElasticMaterial *LEMaterial;

    /// Yield tolerance.
    double yieldTol;
    /// Maximum number of iterations for stress return.
    int newtonIter;

public:
    /// Constructor
    DruckerPragerPlasticitySM(int n, Domain * d);
    /// Destructor
    virtual ~DruckerPragerPlasticitySM();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual const char *giveClassName() const { return "DruckerPragerPlasticitySM"; }
    virtual const char *giveInputRecordName() const { return _IFT_DruckerPragerPlasticitySM_Name; }

    virtual void giveRealStressVector_3d(FloatArray &answer,
                                      GaussPoint *gp,
                                      const FloatArray &strainVector,
                                      TimeStep *tStep);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep);

    /**
     * Perform a standard local stress return using the function computeYieldValue at the specified Gauss point.
     * This function computes the history variables, i.e. the temp variable of plastic strain, hardening variable, state flag, and the temp stress, and stores them in the temp-status.
     * @param gp Gauss point.
     * @param strain Strain vector of this Gauss point.
     */
    void performLocalStressReturn(GaussPoint *gp, const FloatArray &strain);
    /**
     * Check if the trial stress state falls within the vertex region.
     * @param eM Elasticity modulus.
     * @param gM Shear modulus.
     * @param kM Bulk modulus.
     * @return True for vertex case and false if regular stress return has to be used.
     */
    bool checkForVertexCase(double eM, double gM, double kM, double trialStressJTwo, double volumetricStress, double tempKappa);
    /**
     * Perform stress return for regular case, i.e. if the trial stress state does not lie within the vertex region.
     * @param eM Elasticity modulus.
     * @param gM Shear modulus.
     * @param kM Bulk modulus.
     */
    void performRegularReturn(double eM, double gM, double kM, double trialStressJTwo, FloatArray &stressDeviator, double &volumetricStress, double &tempKappa);
    /**
     * Perform stress return for vertex case, i.e. if the trial stress state lies within the vertex region.
     * @param eM Elasticity modulus.
     * @param gM Shear modulus.
     * @param kM Bulk modulus.
     */
    void performVertexReturn(double eM, double gM, double kM, double trialStressJTwo, FloatArray &stressDeviator, double &volumetricStress, double &tempKappa, double volumetricElasticTrialStrain, double kappa);
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
     * @param answer Consistent stiffness matrix.
     * @param mode Material reponse mode.
     * @param gp Gauss point.
     * @param tStep Time step.
     */
    void giveRegAlgorithmicStiffMatrix(FloatMatrix &answer,
                                       MatResponseMode mode,
                                       GaussPoint *gp,
                                       TimeStep *tStep);
    /**
     * Compute consistent stiffness matrix for the vertex case.
     * @param answer Consistent stiffness matrix.
     * @param mode Material reponse mode.
     * @param gp Gauss point.
     * @param tStep Time step.
     */
    void giveVertexAlgorithmicStiffMatrix(FloatMatrix &answer,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer,
                            GaussPoint *gp,
                            InternalStateType type,
                            TimeStep *tStep);

    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
    {
        LEMaterial->giveThermalDilatationVector(answer, gp, tStep);
    }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual double predictRelativeComputationalCost(GaussPoint *gp);
    virtual double predictRelativeRedistributionCost(GaussPoint *gp) { return 1.0; }
};
} // end namespace oofem
#endif // druckerpragerplasticitysm_h
