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

#ifndef mplasticmaterial_h
#define mplasticmaterial_h

#include "structuralmaterial.h"
#include "linearelasticmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"

#include "structuralms.h"
#ifndef __MAKEDEPEND
 #include <vector>
#endif

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
    int state_flag;
    int temp_state_flag;

    /// Consistency parameter values (needed for algorithmic stiffness).
    FloatArray gamma, tempGamma;
    /// Active set of yield functions (needed for algorithmic stiffness).
    IntArray activeConditionMap, tempActiveConditionMap;

public:
    MPlasticMaterialStatus(int n, Domain *d, GaussPoint *g);
    virtual ~MPlasticMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    /// Returns the equilibrated strain vector.
    void givePlasticStrainVector(FloatArray &answer) const { answer = plasticStrainVector; }
    /// Returns the actual (temp) strain vector.
    void giveTempPlasticStrainVector(FloatArray &answer) const { answer = tempPlasticStrainVector; }
    /// Returns the equilibrated hardening variable vector.
    void giveStrainSpaceHardeningVars(FloatArray &answer) const
    { answer = strainSpaceHardeningVarsVector; }
    /// Returns the actual (temp) hardening variable vector.
    void givetempStrainSpaceHardeningVarsVector(FloatArray &answer) const
    { answer = tempStrainSpaceHardeningVarsVector; }

    void letPlasticStrainVectorBe(const FloatArray &v)
    { plasticStrainVector = v; }
    void letTempPlasticStrainVectorBe(const FloatArray &v)
    { tempPlasticStrainVector = v; }
    void letTempStrainSpaceHardeningVarsVectorBe(const FloatArray &v)
    { tempStrainSpaceHardeningVarsVector = v; }
    void letStrainSpaceHardeningVarsVectorBe(const FloatArray &v)
    { strainSpaceHardeningVarsVector = v; }

    int giveStateFlag() { return state_flag; }
    int giveTempStateFlag() { return temp_state_flag; }
    void letTempStateFlagBe(int v) { temp_state_flag = v; }

    void giveTempActiveConditionMap(IntArray &answer) { answer = tempActiveConditionMap; }
    void setTempActiveConditionMap(const IntArray &v) { tempActiveConditionMap = v; }
    void giveTempGamma(FloatArray &answer) { answer = tempGamma; }
    void setTempGamma(const FloatArray &v) { tempGamma = v; }

    // definition
    virtual const char *giveClassName() const { return "MPlasticMaterialStatus"; }
    virtual classType giveClassID() const { return MPlasticMaterialStatusClass; }
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
    LinearElasticMaterial *linearElasticMaterial;
    /// Number of yield surfaces.
    int nsurf;
    /// Protected type to determine the return mapping algorithm.
    enum ReturnMappingAlgoType { mpm_ClosestPoint, mpm_CuttingPlane } rmType;
    /// Type that allows to distinguish between yield function and loading function.
    enum functType { yieldFunction, loadFunction };
    enum plastType { associatedPT, nonassociatedPT } plType;

public:
    MPlasticMaterial(int n, Domain *d);
    ~MPlasticMaterial();

    // identification and auxiliary functions
    virtual int hasNonLinearBehaviour() { return 1; }
    virtual int hasMaterialModeCapability(MaterialMode mode);
    virtual const char *giveClassName() const { return "MPlasticMaterial"; }
    virtual classType giveClassID() const { return MPlasticMaterialClass; }

    /// Returns reference to undamaged (bulk) material.
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    /**
     * Returns true if stiffness matrix of receiver is symmetric.
     * Default implementation returns true.
     */
    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return true; }

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix & answer,
                                               MatResponseForm, MatResponseMode,
                                               GaussPoint * gp,
                                               TimeStep * atTime);


    virtual void giveRealStressVector(FloatArray & answer,  MatResponseForm, GaussPoint *,
                              const FloatArray &, TimeStep *);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);

    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);

    // auxiliary functions
    virtual int giveSizeOfFullHardeningVarsVector()  { return 0; }
    virtual int giveSizeOfReducedHardeningVarsVector(GaussPoint *)  { return 0; }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
    void closestPointReturn(FloatArray &answer, IntArray &activeConditionMap, FloatArray &gamma,
                            MatResponseForm form, GaussPoint *gp,
                            const FloatArray &totalStrain, FloatArray &plasticStrainR,
                            FloatArray &strainSpaceHardeningVariables, TimeStep *atTime);

    void cuttingPlaneReturn(FloatArray &answer, IntArray &activeConditionMap, FloatArray &gamma,
                            MatResponseForm form, GaussPoint *gp,
                            const FloatArray &totalStrain, FloatArray &plasticStrainR,
                            FloatArray &strainSpaceHardeningVariables, TimeStep *atTime);

    void computeGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &fullStressVector,
                               const FloatArray &fullStressSpaceHardeningVars);
    void computeResidualVector(FloatArray &answer, GaussPoint *gp, const FloatArray &gamma,
                               const IntArray &activeConditionMap, const FloatArray &plasticStrainVectorR,
                               const FloatArray &strainSpaceHardeningVariables, std :: vector< FloatArray > &gradVec);
    virtual void giveConsistentStiffnessMatrix(FloatMatrix & answer,
                                               MatResponseForm,
                                               MatResponseMode,
                                               GaussPoint * gp,
                                               TimeStep * atTime);

    virtual void giveElastoPlasticStiffnessMatrix(FloatMatrix &answer,
                                                  MatResponseForm form,
                                                  MatResponseMode mode,
                                                  GaussPoint *gp,
                                                  TimeStep *atTime);

    void computeAlgorithmicModuli(FloatMatrix &answer,
                                  GaussPoint *gp, const FloatMatrix &elasticModuliInverse,
                                  const FloatMatrix &hardeningModuliInverse,
                                  const FloatArray &gamma, const IntArray &activeConditionMap,
                                  const FloatArray &fullStressVector,
                                  const FloatArray &fullStressSpaceHardeningVars);
    void computeDiagModuli(FloatMatrix &answer,
                           GaussPoint *gp, FloatMatrix &elasticModuliInverse,
                           FloatMatrix &hardeningModuliInverse);

    virtual void computeStressSpaceHardeningVars(FloatArray &answer, GaussPoint *gp,
                                                 const FloatArray &strainSpaceHardeningVariables) = 0;
    virtual double computeYieldValueAt(GaussPoint *gp, int isurf, const FloatArray &stressVector,
                                       const FloatArray &stressSpaceHardeningVars) = 0;
    virtual void computeHardeningReducedModuli(FloatMatrix &answer,
                                               GaussPoint *gp,
                                               const FloatArray &strainSpaceHardeningVariables,
                                               TimeStep *atTime) = 0;
    virtual void computeStressGradientVector(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp, const FloatArray &stressVector,
                                             const FloatArray &stressSpaceHardeningVars) = 0;
    virtual void computeStressSpaceHardeningVarsReducedGradient(FloatArray &answer, functType ftype, int isurf, GaussPoint *gp,
                                                                const FloatArray &stressVector,
                                                                const FloatArray &stressSpaceHardeningVars) = 0;
    virtual int hasHardening() { return 0; }
    virtual void computeReducedGradientMatrix(FloatMatrix &answer, int isurf,
                                              GaussPoint *gp,
                                              const FloatArray &stressVector,
                                              const FloatArray &stressSpaceHardeningVars) = 0;

    virtual void computeTrialStressIncrement(FloatArray &answer, GaussPoint *gp,
                                             const FloatArray &strainIncrement, TimeStep *atTime);
    virtual void computeReducedElasticModuli(FloatMatrix &answer, GaussPoint *gp,
                                             TimeStep *atTime);
    //virtual void compute3dElasticModuli(FloatMatrix& answer, GaussPoint *gp,
    //                                    TimeStep *atTime) = 0;

    // next functions overloaded from structural material level
    virtual void givePlaneStressStiffMtrx(FloatMatrix & answer,
                                          MatResponseForm, MatResponseMode,
                                          GaussPoint * gp,
                                          TimeStep * atTime);
    virtual void givePlaneStrainStiffMtrx(FloatMatrix & answer,
                                          MatResponseForm, MatResponseMode,
                                          GaussPoint * gp,
                                          TimeStep * atTime);
    virtual void give1dStressStiffMtrx(FloatMatrix & answer,
                                       MatResponseForm, MatResponseMode,
                                       GaussPoint * gp,
                                       TimeStep * atTime);
    virtual void give2dBeamLayerStiffMtrx(FloatMatrix & answer,
                                          MatResponseForm, MatResponseMode,
                                          GaussPoint * gp,
                                          TimeStep * atTime);
    virtual void give2dPlateLayerStiffMtrx(FloatMatrix & answer,
                                           MatResponseForm, MatResponseMode,
                                           GaussPoint * gp,
                                           TimeStep * atTime);

    virtual void give1dFiberStiffMtrx(FloatMatrix & answer,
                                      MatResponseForm, MatResponseMode, GaussPoint * gp,
                                      TimeStep * atTime);

    virtual void give3dShellLayerStiffMtrx(FloatMatrix & answer,
                                           MatResponseForm, MatResponseMode,
                                           GaussPoint * gp,
                                           TimeStep * atTime);
};
} // end namespace oofem
#endif // mplasticmaterial_h
