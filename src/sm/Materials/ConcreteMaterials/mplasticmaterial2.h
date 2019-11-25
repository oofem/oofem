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

#ifndef mplasticmaterial2_h
#define mplasticmaterial2_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/linearelasticmaterial.h"
#include "dictionary.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "sm/Materials/structuralms.h"

#include <vector>
#include <set>

namespace oofem {
class GaussPoint;

/**
 * This class implements associated Material Status to MPlasticMaterial.
 * It is attribute of matStatusDictionary at every GaussPoint, for which this material
 * is active.
 *
 * Description:
 * Idea used there is that we have variables
 * describing:
 * -# state at previous equilibrium state (variables without temp)
 * -# state during searching new equilibrium (variables with temp)
 *    when we start search new state from previous equilibrium one we copy
 *    non-tem variables into temp ones. And after we reach new equilibrium
 *    (now described by temp variables) we copy tem-var into non-tepm ones
 *    (see function updateYourself).
 */
class MPlasticMaterial2Status : public StructuralMaterialStatus
{
public:
    enum state_flag_values { PM_Elastic, PM_Yielding, PM_Unloading };

protected:
    /// Plastic strain vector.
    FloatArray plasticStrainVector;
    FloatArray tempPlasticStrainVector;

    /**
     * Strain space hardening variables, e.g. cumulative plastic strain which can be used
     * in yield stress evaluation or damage.
     */
    FloatArray strainSpaceHardeningVarsVector;
    FloatArray tempStrainSpaceHardeningVarsVector;

    /// Yield function status indicator.
    state_flag_values state_flag = MPlasticMaterial2Status :: PM_Elastic;
    state_flag_values temp_state_flag = MPlasticMaterial2Status :: PM_Elastic;

    /// Isotropic damage variables
    double damage = 0., tempDamage = 0.;

    /// Consistency parameter values (needed for algorithmic stiffness).
    FloatArray gamma, tempGamma;
    /// Active set of yield functions (needed for algorithmic stiffness).
    IntArray activeConditionMap, tempActiveConditionMap;

public:
    MPlasticMaterial2Status(GaussPoint * g, int statusSize);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    /// Returns the equilibrated strain vector.
    const FloatArray & givePlasticStrainVector() const { return plasticStrainVector; }
    /// Returns the actual (temp) strain vector.
    const FloatArray & giveTempPlasticStrainVector() const { return tempPlasticStrainVector; }
    /// Returns the equilibrated hardening variable vector.
    const FloatArray & giveStrainSpaceHardeningVars() const { return strainSpaceHardeningVarsVector; }
    /// Returns the actual (temp) hardening variable vector.
    const FloatArray & giveTempStrainSpaceHardeningVarsVector() const { return tempStrainSpaceHardeningVarsVector; }

    void letPlasticStrainVectorBe(FloatArray v) { plasticStrainVector = std :: move(v); }
    void letTempPlasticStrainVectorBe(FloatArray v) { tempPlasticStrainVector = std :: move(v); }
    void letTempStrainSpaceHardeningVarsVectorBe(FloatArray v) { tempStrainSpaceHardeningVarsVector = std :: move(v); }
    void letStrainSpaceHardeningVarsVectorBe(FloatArray v) { strainSpaceHardeningVarsVector = std :: move(v); }

    void letTempDamageBe(double v) { tempDamage = v; }
    double giveDamage() { return damage; }
    double giveTempDamage() { return tempDamage; }

    state_flag_values giveStateFlag() { return state_flag; }
    state_flag_values giveTempStateFlag() { return temp_state_flag; }
    void letTempStateFlagBe(state_flag_values v) { temp_state_flag = v; }

    const IntArray &giveTempActiveConditionMap() { return tempActiveConditionMap; }
    void setTempActiveConditionMap(IntArray v) { tempActiveConditionMap = std :: move(v); }
    const FloatArray &giveTempGamma() { return tempGamma; }
    void setTempGamma(FloatArray v) { tempGamma = std :: move(v); }

    const char *giveClassName() const override { return "MPlasticMaterial2Status"; }
};

/**
 * This class represents a base class for non-associated multisurface plasticity.
 * The Multisurface plasticity is characterized by the following:
 * Let @f$\sigma, \varepsilon @f$, and @f$\varepsilon^p @f$ be the stress, total strain, and plastic strain vectors, respectively.
 * It is assumed that the total strain is decomposed into reversible elastic and irreversible plastic parts
 * @f$\varepsilon=\varepsilon^e+\varepsilon^p @f$.
 * The elastic response is characterized in terms of elastic constitutive matrix @f$ D^e @f$ as @f$ \sigma=D^e(\varepsilon-\varepsilon^e) @f$
 * As long as the stress remains inside the elastic domain, the deformation process is purely elastic and the
 * plastic strain does not change.
 *
 * It is assumed that the elastic domain, denoted as @f$ IE @f$ is bounded by a composite yield surface. It is defined as
 * @f[
 * IE=\{(\sigma,\kappa)|f_i(\sigma,\kappa)<0,\; i\in\{1,\cdots,m\}\}
 * @f]
 * where @f$ f_i(\sigma,\kappa) @f$ are @f$ m\ge1 @f$ yield functions intersecting in a possibly non-smooth fashion. The
 * vector @f$ \kappa @f$ contains internal variables controlling the evolution of yield surfaces (amount of hardening or softening).
 * The evolution of plastic strain @f$ \varepsilon^p @f$ is expressed in Koiter's form. Assuming the non-associated plasticity, this reads
 * @f[
 * \label{epe}
 * \varepsilon^p=\sum^{m}_{i=1} \lambda^i \partial_{\sigma}g_i(\sigma,\kappa)
 * @f]
 * where @f$ g_i @f$ are plastic potential functions. The @f$ \lambda^i @f$ are referred as plastic consistency parameters, which satisfy the following Kuhn-Tucker conditions
 * @f[
 * \label{ktc}
 * \lambda^i\ge0,\;f_i\le0,\;{\rm and}\ \lambda^i f_i=0
 * @f]
 * These conditions imply that in the elastic regime the yield function must remain negative and the rate of the plastic multiplier is zero
 * (plastic strain remains constant) while in the plastic regime the yield function must be equal to zero (stress remains on the surface) and the rate of the plastic multiplier is positive.
 * The evolution of vector of internal hardening/softening variables @f$ \kappa @f$  is expressed in terms of a general
 * hardening/softening law of the form
 * @f[
 * \dot{\kappa} = \dot{\kappa}(\sigma, \lambda)
 * @f]
 * where @f$ \lambda @f$ is the vector of plastic consistency parameters @f$ \lambda_i @f$.
 *
 */
class MPlasticMaterial2 : public StructuralMaterial
{
protected:
    /// Reference to bulk (undamaged) material.
    LinearElasticMaterial *linearElasticMaterial;
    /// Number of yield surfaces.
    int nsurf;
    /// Protected type to determine the return mapping algorithm.
    enum ReturnMappingAlgoType { mpm_ClosestPoint, mpm_CuttingPlane } rmType;
    /// Type that allows to distinguish between yield function and loading function.
    enum functType { yieldFunction, loadFunction };
    enum plastType { associatedPT, nonassociatedPT } plType;
    /// Flag indicating whether iterative update of a set of active yield conditions takes place.
    bool iterativeUpdateOfActiveConds;
    /// Set for keeping record of generated populations of active yield conditions during return.
    mutable std :: set< long >populationSet; // FIXME: Race condition when threaded.

public:
    MPlasticMaterial2(int n, Domain * d);
    virtual ~MPlasticMaterial2();

    bool hasMaterialModeCapability(MaterialMode mode) const override;
    const char *giveClassName() const override { return "MPlasticMaterial2"; }

    /// Returns reference to undamaged (bulk) material
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return true; }

    FloatArrayF<6> giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const override
    {
        double alpha = this->linearElasticMaterial->give(tAlpha, gp);
        return {alpha, alpha, alpha, 0., 0., 0.};
    }

    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const override;

    void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                              const FloatArray &strain, TimeStep *tStep) override;

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<MPlasticMaterial2*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<4> giveRealStressVector_PlaneStrain(const FloatArrayF<4> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<MPlasticMaterial2*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<3> giveRealStressVector_PlaneStress(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<MPlasticMaterial2*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<1> giveRealStressVector_1d(const FloatArrayF<1> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<MPlasticMaterial2*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }

    virtual double computeDamage(GaussPoint *gp, const FloatArray &strainSpaceHardeningVariables, TimeStep *tStep);

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    // auxiliary functions
    virtual int giveSizeOfFullHardeningVarsVector() const { return 0; }
    virtual int giveSizeOfReducedHardeningVarsVector(GaussPoint *) const { return 0; }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

protected:

    virtual int giveMaxNumberOfActiveYieldConds(GaussPoint *gp) const = 0;
    void closestPointReturn(FloatArray &answer, IntArray &activeConditionMap, FloatArray &gamma,
                            GaussPoint *gp,
                            const FloatArray &totalStrain, FloatArray &plasticStrainR,
                            FloatArray &strainSpaceHardeningVariables, TimeStep *tStep) const;

    void cuttingPlaneReturn(FloatArray &answer, IntArray &activeConditionMap, FloatArray &gamma,
                            GaussPoint *gp,
                            const FloatArray &totalStrain, FloatArray &plasticStrainR,
                            FloatArray &strainSpaceHardeningVariables, TimeStep *tStep) const;

    // add here some auxiliary functions if needed
    /* void  computeGradientVector (FloatArray& answer, functType ftype, int isurf, GaussPoint* gp, const FloatArray& fullStressVector,
     * const FloatArray& fullStressSpaceHardeningVars);*/
    void computeResidualVector(FloatArray &answer, GaussPoint *gp, const FloatArray &gamma,
                               const IntArray &activeConditionMap, const FloatArray &plasticStrainVectorR,
                               const FloatArray &strainSpaceHardeningVariables, std :: vector< FloatArray > &gradVec) const;
    virtual FloatMatrix giveConsistentStiffnessMatrix(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const;

    virtual FloatMatrix giveElastoPlasticStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const;

    void computeAlgorithmicModuli(FloatMatrix &answer,
                                  GaussPoint *gp, const FloatMatrix &elasticModuliInverse,
                                  const FloatArray &gamma, const IntArray &activeConditionMap,
                                  const FloatArray &fullStressVector,
                                  const FloatArray &strainSpaceHardeningVariables) const;

    /*void  computeDiagModuli(FloatMatrix& answer,
     *                      GaussPoint *gp, FloatMatrix &elasticModuliInverse,
     *                      FloatMatrix &hardeningModuliInverse);*/

    /// Computes the value of yield function.
    virtual double computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector,
                                       const FloatArray &strainSpaceHardeningVariables) const = 0;

    /// Computes the stress gradient of yield/loading function (df/d_sigma).
    virtual void computeStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &stressVector,
                                             const FloatArray &strainSpaceHardeningVariables) const = 0;
    void computeReducedStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &stressVector,
                                            const FloatArray &strainSpaceHardeningVariables) const;

    /**
     * Computes the increment of strain-space hardening variables.
     * @param answer Result.
     * @param gp Gauss point to compute at.
     * @param stress Updated stress (corresponds to newly reached state).
     * @param dlambda Increment of consistency parameters.
     * @param dplasticStrain Actual plastic strain increment.
     * @param activeConditionMap Array of active yield conditions.
     */
    virtual void computeStrainHardeningVarsIncrement(FloatArray &answer, GaussPoint *gp,
                                                     const FloatArray &stress, const FloatArray &dlambda,
                                                     const FloatArray &dplasticStrain, const IntArray &activeConditionMap) const = 0;
    /**
     * Computes the derivative of yield/loading function with respect to @f$ \kappa @f$ vector
     */
    virtual void computeKGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, FloatArray &fullStressVector,
                                        const FloatArray &strainSpaceHardeningVariables) const = 0;

    /**
     * Computes derivative of @f$ \kappa @f$ vector with respect to stress
     */
    virtual void computeReducedHardeningVarsSigmaGradient(FloatMatrix &answer, GaussPoint *gp, const IntArray &activeConditionMap,
                                                          const FloatArray &fullStressVector,
                                                          const FloatArray &strainSpaceHardeningVars,
                                                          const FloatArray &gamma) const = 0;
    /// computes derivative of @f$ \kappa @f$ vector with respect to lambda vector
    virtual void computeReducedHardeningVarsLamGradient(FloatMatrix &answer, GaussPoint *gp, int actSurf,
                                                        const IntArray &activeConditionMap,
                                                        const FloatArray &fullStressVector,
                                                        const FloatArray &strainSpaceHardeningVars,
                                                        const FloatArray &gamma) const = 0;
    /**
     * Indicates, whether receiver model has hardening/softening behavior or behaves according to perfect plasticity theory.
     */
    virtual int hasHardening() const = 0;
    /* virtual void  computeReducedGradientMatrix (FloatMatrix& answer, int isurf,
     *                                          GaussPoint *gp,
     *                                          const FloatArray& stressVector,
     *                                          const FloatArray& stressSpaceHardeningVars) = 0;*/
    /// Computes second derivative of loading function with respect to stress.
    virtual void computeReducedSSGradientMatrix(FloatMatrix &gradientMatrix,  int i, GaussPoint *gp, const FloatArray &fullStressVector,
                                                const FloatArray &strainSpaceHardeningVariables) const = 0;
    /// Computes second derivative of loading function with respect to stress and hardening vars
    virtual void computeReducedSKGradientMatrix(FloatMatrix &gradientMatrix,  int i, GaussPoint *gp, const FloatArray &fullStressVector,
                                                const FloatArray &strainSpaceHardeningVariables) const = 0;

    /**
     * Computes full-space trial stress increment (elastic).
     * @param answer Contains result in full space.
     * @param gp Integration point.
     * @param strainIncrement Strain increment vector.
     * @param tStep Solution step.
     */
    virtual void computeTrialStressIncrement(FloatArray &answer, GaussPoint *gp,
                                             const FloatArray &strainIncrement, TimeStep *tStep) const;
    virtual void computeReducedElasticModuli(FloatMatrix &answer, GaussPoint *gp,
                                             TimeStep *tStep) const;
    //virtual void compute3dElasticModuli(FloatMatrix& answer, GaussPoint *gp,
    //                                    TimeStep *tStep) = 0;

    // next functions overloaded rom structural material level
    FloatMatrixF<3,3> givePlaneStressStiffMtrx(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<4,4> givePlaneStrainStiffMtrx(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<1,1> give1dStressStiffMtrx(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<2,2> give2dBeamLayerStiffMtrx(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<5,5> givePlateLayerStiffMtrx(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> giveFiberStiffMtrx(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const override;

protected:
    long getPopulationSignature(IntArray &mask) const;
    int testPopulation(long pop) const;
    void clearPopulationSet() const;
    void addNewPopulation(IntArray &mask) const;
    int getNewPopulation(IntArray &result, IntArray &candidateMask, int degree, int size) const;
};
} // end namespace oofem
#endif // mplasticmaterial2_h
