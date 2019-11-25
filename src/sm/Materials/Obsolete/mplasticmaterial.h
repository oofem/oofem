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

#ifndef mplasticmaterial_h
#define mplasticmaterial_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/linearelasticmaterial.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "sm/Materials/structuralms.h"

#include <vector>

namespace oofem {
class GaussPoint;

/**
 * This class implements associated Material Status to MPlasticMaterial.
 * It is attribute of matStatusDictionary at every GaussPoint, for which this material
 * is active.
 *
 * Idea used there is that we have variables
 * describing:
 * -# state at previous equilibrium state (variables without temp)
 * -# state during searching new equilibrium (variables with temp)
 *    when we start search new state from previous equilibrium one we copy
 *    non-temp variables into temp ones. And after we reach new equilibrium
 *    (now described by temp variables) we copy temp-var into non-temp ones
 *    (see function updateYourself).
 */
class MPlasticMaterialStatus : public StructuralMaterialStatus
{
public:
    enum state_flag_values { PM_Elastic, PM_Yielding, PM_Unloading };

protected:
    /// Plastic strain vector.
    FloatArray plasticStrainVector;
    FloatArray tempPlasticStrainVector;

    /// Strain space hardening variables.
    FloatArray strainSpaceHardeningVarsVector;
    FloatArray tempStrainSpaceHardeningVarsVector;

    /// Yield function status indicator.
    int state_flag = MPlasticMaterialStatus :: PM_Elastic;
    int temp_state_flag = MPlasticMaterialStatus :: PM_Elastic;

    /// Consistency parameter values (needed for algorithmic stiffness).
    FloatArray gamma, tempGamma;
    /// Active set of yield functions (needed for algorithmic stiffness).
    IntArray activeConditionMap, tempActiveConditionMap;

public:
    MPlasticMaterialStatus(GaussPoint * g, int statusSize);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    /// Returns the equilibrated strain vector.
    const FloatArray &givePlasticStrainVector() const { return plasticStrainVector; }
    /// Returns the actual (temp) strain vector.
    const FloatArray &giveTempPlasticStrainVector() const { return tempPlasticStrainVector; }
    /// Returns the equilibrated hardening variable vector.
    const FloatArray &giveStrainSpaceHardeningVars() const
    { return strainSpaceHardeningVarsVector; }
    /// Returns the actual (temp) hardening variable vector.
    const FloatArray &givetempStrainSpaceHardeningVarsVector() const
    { return tempStrainSpaceHardeningVarsVector; }

    void letPlasticStrainVectorBe(FloatArray v)
    { plasticStrainVector = std :: move(v); }
    void letTempPlasticStrainVectorBe(FloatArray v)
    { tempPlasticStrainVector = std :: move(v); }
    void letTempStrainSpaceHardeningVarsVectorBe(FloatArray v)
    { tempStrainSpaceHardeningVarsVector = std :: move(v); }
    void letStrainSpaceHardeningVarsVectorBe(FloatArray v)
    { strainSpaceHardeningVarsVector = std :: move(v); }

    int giveStateFlag() { return state_flag; }
    int giveTempStateFlag() { return temp_state_flag; }
    void letTempStateFlagBe(int v) { temp_state_flag = v; }

    const IntArray &giveTempActiveConditionMap() { return tempActiveConditionMap; }
    void setTempActiveConditionMap(IntArray v) { tempActiveConditionMap = std :: move(v); }
    const FloatArray &giveTempGamma() { return tempGamma; }
    void setTempGamma(FloatArray v) { tempGamma = std :: move(v); }

    // definition
    const char *giveClassName() const override { return "MPlasticMaterialStatus"; }
};

/**
 * This class implements a general plastic material.
 * It is assumed to be a base class for many material
 * models based on different yield conditions and hardening laws.
 * Sress return mapping algorithm  is based on general
 * return mapping algorithm, with following assumptions
 * - associative flow rule
 * - general hardening law is supported
 * - Kunt-Tucker conditions apply.
 */
class MPlasticMaterial : public StructuralMaterial
{
protected:
    /// Reference to bulk (undamaged) material.
    LinearElasticMaterial *linearElasticMaterial = nullptr;
    /// Number of yield surfaces.
    int nsurf = 0;
    /// Protected type to determine the return mapping algorithm.
    enum ReturnMappingAlgoType { mpm_ClosestPoint, mpm_CuttingPlane } rmType = mpm_ClosestPoint;
    /// Type that allows to distinguish between yield function and loading function.
    enum functType { yieldFunction, loadFunction };
    enum plastType { associatedPT, nonassociatedPT } plType = associatedPT;

public:
    MPlasticMaterial(int n, Domain * d);
    virtual ~MPlasticMaterial();

    // identification and auxiliary functions
    bool hasMaterialModeCapability(MaterialMode mode) const override;
    const char *giveClassName() const override { return "MPlasticMaterial"; }

    /// Returns reference to undamaged (bulk) material.
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    /**
     * Returns true if stiffness matrix of receiver is symmetric.
     * Default implementation returns true.
     */
    bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) const override { return true; }

    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const override;


    void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                              const FloatArray &, TimeStep *tStep) override;

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<MPlasticMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<4> giveRealStressVector_PlaneStrain(const FloatArrayF<4> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<MPlasticMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<3> giveRealStressVector_PlaneStress(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<MPlasticMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<1> giveRealStressVector_1d(const FloatArrayF<1> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<MPlasticMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    // auxiliary functions
    virtual int giveSizeOfFullHardeningVarsVector() const { return 0; }
    virtual int giveSizeOfReducedHardeningVarsVector(GaussPoint *) const { return 0; }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

protected:
    void closestPointReturn(FloatArray &answer, IntArray &activeConditionMap, FloatArray &gamma,
                            GaussPoint *gp,
                            const FloatArray &totalStrain, FloatArray &plasticStrainR,
                            FloatArray &strainSpaceHardeningVariables, TimeStep *tStep);

    void cuttingPlaneReturn(FloatArray &answer, IntArray &activeConditionMap, FloatArray &gamma,
                            GaussPoint *gp,
                            const FloatArray &totalStrain, FloatArray &plasticStrainR,
                            FloatArray &strainSpaceHardeningVariables, TimeStep *tStep);

    void computeGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &fullStressVector,
                               const FloatArray &fullStressSpaceHardeningVars) const;
    void computeResidualVector(FloatArray &answer, GaussPoint *gp, const FloatArray &gamma,
                               const IntArray &activeConditionMap, const FloatArray &plasticStrainVectorR,
                               const FloatArray &strainSpaceHardeningVariables, std :: vector< FloatArray > &gradVec);
    virtual FloatMatrix giveConsistentStiffnessMatrix(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const;

    virtual void giveElastoPlasticStiffnessMatrix(FloatMatrix &answer,
                                                  MatResponseMode mode,
                                                  GaussPoint *gp,
                                                  TimeStep *tStep) const;

    void computeAlgorithmicModuli(FloatMatrix &answer,
                                  GaussPoint *gp, const FloatMatrix &elasticModuliInverse,
                                  const FloatMatrix &hardeningModuliInverse,
                                  const FloatArray &gamma, const IntArray &activeConditionMap,
                                  const FloatArray &fullStressVector,
                                  const FloatArray &fullStressSpaceHardeningVars) const;
    void computeDiagModuli(FloatMatrix &answer,
                           GaussPoint *gp, FloatMatrix &elasticModuliInverse,
                           FloatMatrix &hardeningModuliInverse) const ;

    virtual void computeStressSpaceHardeningVars(FloatArray &answer, GaussPoint *gp,
                                                 const FloatArray &strainSpaceHardeningVariables) const = 0;
    virtual double computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector,
                                       const FloatArray &stressSpaceHardeningVars) const = 0;
    virtual void computeHardeningReducedModuli(FloatMatrix &answer,
                                               GaussPoint *gp,
                                               const FloatArray &strainSpaceHardeningVariables,
                                               TimeStep *tStep) const = 0;
    virtual void computeStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &stressVector,
                                             const FloatArray &stressSpaceHardeningVars) const = 0;
    virtual void computeStressSpaceHardeningVarsReducedGradient(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp,
                                                                const FloatArray &stressVector,
                                                                const FloatArray &stressSpaceHardeningVars) const = 0;
    virtual int hasHardening() const { return 0; }
    virtual void computeReducedGradientMatrix(FloatMatrix &answer, int isurf,
                                              GaussPoint *gp,
                                              const FloatArray &stressVector,
                                              const FloatArray &stressSpaceHardeningVars) const = 0;

    virtual void computeTrialStressIncrement(FloatArray &answer, GaussPoint *gp,
                                             const FloatArray &strainIncrement, TimeStep *tStep) const;
    virtual void computeReducedElasticModuli(FloatMatrix &answer, GaussPoint *gp,
                                             TimeStep *tStep) const;
    //virtual void compute3dElasticModuli(FloatMatrix& answer, GaussPoint *gp,
    //                                    TimeStep *tStep) = 0;

    // next functions overloaded from structural material level
    FloatMatrixF<3,3> givePlaneStressStiffMtrx(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<4,4> givePlaneStrainStiffMtrx(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<1,1> give1dStressStiffMtrx(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<2,2> give2dBeamLayerStiffMtrx(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<5,5> givePlateLayerStiffMtrx(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF<3,3> giveFiberStiffMtrx(MatResponseMode, GaussPoint *gp, TimeStep *tStep) const override;
};
} // end namespace oofem
#endif // mplasticmaterial_h
