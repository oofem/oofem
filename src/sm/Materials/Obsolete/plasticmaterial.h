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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef plasticmaterial_h
#define plasticmaterial_h

#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/linearelasticmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "sm/Materials/structuralms.h"
#include "matstatmapperint.h"

namespace oofem {
class GaussPoint;

enum state_flag_values { PM_Elastic, PM_Yielding, PM_Unloading };
/**
 * This class implements associated Material Status to PlasticMaterial.
 */
class PlasticMaterialStatus : public StructuralMaterialStatus
{
protected:
    /// Plastic strain vector (reduced form).
    FloatArray plasticStrainVector;
    FloatArray tempPlasticStrainVector;

    /// Strain space hardening variables.
    FloatArray strainSpaceHardeningVarsVector;
    FloatArray tempStrainSpaceHardeningVarsVector;

    /// Yield function status indicator.
    int state_flag = PM_Elastic;
    int temp_state_flag = PM_Elastic;

    /// Plastic consistency parameter.
    double gamma = 0., temp_gamma = 0.;

public:
    PlasticMaterialStatus(GaussPoint * g, int statusSize);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    const FloatArray &givePlasticStrainVector() const { return plasticStrainVector; }
    const FloatArray &giveTempPlasticStrainVector() const { return tempPlasticStrainVector; }
    const FloatArray &giveStrainSpaceHardeningVars() const { return strainSpaceHardeningVarsVector; }
    const FloatArray &givetempStrainSpaceHardeningVarsVector() const { return tempStrainSpaceHardeningVarsVector; }

    void letPlasticStrainVectorBe(FloatArray v) { plasticStrainVector = std :: move(v); }
    void letTempPlasticStrainVectorBe(FloatArray v) { tempPlasticStrainVector = std :: move(v); }
    void letTempStrainSpaceHardeningVarsVectorBe(FloatArray v) { tempStrainSpaceHardeningVarsVector = std :: move(v); }
    void letStrainSpaceHardeningVarsVectorBe(FloatArray v) { strainSpaceHardeningVarsVector = std :: move(v); }

    int giveStateFlag() const { return state_flag; }
    int giveTempStateFlag() const { return temp_state_flag; }
    double givePlasticConsistencyPrameter() const { return gamma; }
    double giveTempPlasticConsistencyPrameter() const { return temp_gamma; }
    void letTempPlasticConsistencyPrameterBe(double v) { gamma = v; }

    void letTempStateFlagBe(int v) { temp_state_flag = v; }

    const char *giveClassName() const override { return "PlasticMaterialStatus"; }

    /// Functions for MaterialStatusMapperInterface
    void copyStateVariables(const MaterialStatus &iStatus) override;
    void addStateVariables(const MaterialStatus &iStatus) override;
};

/**
 * This class implements a general plastic material.
 * It is assumed to be a base class for many material
 * models based on different yield conditions and hardening laws.
 * Stress return mapping algorithm  is based on general
 * return mapping algorithm, with following assumptions
 * - associative flow rule
 * - general hardening law is supported
 * - Kunt-Tucker conditions apply.
 */
class PlasticMaterial : public StructuralMaterial
{
protected:
    /// Reference to bulk (undamaged) material
    LinearElasticMaterial *linearElasticMaterial = nullptr;

public:
    PlasticMaterial(int n, Domain * d);
    virtual ~PlasticMaterial();

    // identification and auxiliary functions
    bool hasMaterialModeCapability(MaterialMode mode) const override;
    const char *giveClassName() const override { return "PlasticMaterial"; }

    /// Returns reference to undamaged (bulk) material.
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    FloatMatrixF<6,6> give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                              const FloatArray &reducedStrain, TimeStep *tStep) const override;

    FloatArrayF<6> giveRealStressVector_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<PlasticMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<4> giveRealStressVector_PlaneStrain(const FloatArrayF<4> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<PlasticMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<3> giveRealStressVector_PlaneStress(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<PlasticMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<1> giveRealStressVector_1d(const FloatArrayF<1> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<PlasticMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<2> giveRealStressVector_2dBeamLayer(const FloatArrayF<2> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<PlasticMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<5> giveRealStressVector_PlateLayer(const FloatArrayF<5> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<PlasticMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF<3> giveRealStressVector_Fiber(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast<PlasticMaterial*>(this)->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }


    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;

protected:
    // add here some auxiliary functions if needed
    FloatArray *ComputeGradientVector(GaussPoint *gp, FloatArray *fullStressVector,
                                      FloatArray *fullStressSpaceHardeningVars) const;
    FloatArray *ComputeResidualVector(GaussPoint *gp, double Gamma,
                                      FloatArray *plasticStrainVectorR,
                                      FloatArray *strainSpaceHardeningVariables,
                                      FloatArray *gradientVectorR) const;
    virtual FloatMatrix giveConsistentStiffnessMatrix(MatResponseMode,
                                               GaussPoint *gp,
                                               TimeStep *tStep) const;

    void  computeConsistentModuli(FloatMatrix &answer,
                                  GaussPoint *gp, FloatMatrix &elasticModuliInverse,
                                  FloatMatrix &hardeningModuliInverse,
                                  double Gamma, const FloatArray &fullStressVector,
                                  const FloatArray &fullStressSpaceHardeningVars) const;
    void  computeDiagModuli(FloatMatrix &answer,
                            GaussPoint *gp, FloatMatrix &elasticModuliInverse,
                            FloatMatrix &hardeningModuliInverse) const;

    virtual FloatArray *ComputeStressSpaceHardeningVars(GaussPoint *gp,
                                                        FloatArray *strainSpaceHardeningVariables) const
    { return NULL; }
    virtual double computeYieldValueAt(GaussPoint *gp, FloatArray *stressVector,
                                       FloatArray *stressSpaceHardeningVars) const { return 0.; }
    virtual void   computeHardeningReducedModuli(FloatMatrix &answer,
                                                 GaussPoint *gp,
                                                 FloatArray *strainSpaceHardeningVariables,
                                                 TimeStep *tStep) const = 0;
    virtual FloatArray *ComputeStressGradient(GaussPoint *gp, FloatArray *stressVector,
                                              FloatArray *stressSpaceHardeningVars) const { return NULL; }
    virtual FloatArray *ComputeStressSpaceHardeningVarsReducedGradient(GaussPoint *gp,
                                                                       FloatArray *stressVector,
                                                                       FloatArray *stressSpaceHardeningVars) const
    { return NULL; }
    virtual int hasHardening() const { return 0; }
    virtual void  computeReducedGradientMatrix(FloatMatrix &answer,
                                               GaussPoint *gp,
                                               const FloatArray &stressVector,
                                               const FloatArray &stressSpaceHardeningVars) const = 0;

    virtual void computeTrialStressIncrement(FloatArray &answer, GaussPoint *gp,
                                             const FloatArray &strainIncrement, TimeStep *tStep) const;
    virtual void computeReducedElasticModuli(FloatMatrix &answer, GaussPoint *gp,
                                             TimeStep *tStep) const;
    virtual void compute3dElasticModuli(FloatMatrix &answer, GaussPoint *gp,
                                        TimeStep *tStep) const = 0;

    // auxiliary functions
    virtual int giveSizeOfFullHardeningVarsVector() const { return 0; }
    virtual int giveSizeOfReducedHardeningVarsVector(GaussPoint *) const { return 0; }

    friend class PlasticMaterialStatus;

    // next functions overloaded rom structural material level
    FloatMatrixF<3,3> givePlaneStressStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<4,4> givePlaneStrainStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<1,1> give1dStressStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<2,2> give2dBeamLayerStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<5,5> givePlateLayerStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;

    FloatMatrixF<3,3> giveFiberStiffMtrx(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
};
} // end namespace oofem
#endif // plasticmaterial_h
