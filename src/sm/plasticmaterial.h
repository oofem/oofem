/* $Header: /home/cvs/bp/oofem/sm/src/plasticmaterial.h,v 1.5 2003/04/06 14:08:31 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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


//   ****************************************
//   *** CLASS PLASTIC  MATERIAL ************
//   ****************************************

#ifndef plasticmaterial_h
#define plasticmaterial_h

#include "structuralmaterial.h"
#include "linearelasticmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"

#include "structuralms.h"

namespace oofem {

class GaussPoint;

// material contant's keys for give()
//#define pscm_Ee 300
//#define pscm_Et 301
//#define pscm_Gf 302
//#define pscm_Beta 303
//#define pscm_G  304
//#define pscm_Ft 305

enum state_flag_values { PM_Elastic, PM_Yielding, PM_Unloading };

class PlasticMaterialStatus : public StructuralMaterialStatus
{
    /*
     * This class implements associated Material Status to PlasticMaterial.
     * It is atribute of matStatusDictionary at every GaussPoint, for which this material
     * is active.
     *
     * DESCRIPTION:
     * Idea used there is that we have variables
     * describing:
     *   1) state at previous equilibrium state (variables without temp)
     * 2) state during searching new equilibrium (variables with temp)
     * when we start search new state from previous equilibrium one we copy
     * non-tem variables into temp ones. And after we reach new equilibrium
     * (now decribed by temp variables) we copy tem-var into non-tepm ones
     * (see function updateYourself).
     *
     * variables description:
     *
     *
     * TASK:
     *
     *
     *
     */
protected:

    // add internal variables here

    // plastic strain vector
    FloatArray plasticStrainVector; // reduced form
    FloatArray tempPlasticStrainVector;

    // strain space hardening variables
    FloatArray strainSpaceHardeningVarsVector;
    FloatArray tempStrainSpaceHardeningVarsVector;

    // yield function status indicator
    int state_flag;
    int temp_state_flag;

    // plastic consistency parameter
    double gamma, temp_gamma;

public:
    PlasticMaterialStatus(int n, Domain *d, GaussPoint *g);
    ~PlasticMaterialStatus();

    void   printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *); // update after new equilibrium state reached

    // saves current context(state) into stream
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    // add internal variables access functions and other auxiliary functions here
    void givePlasticStrainVector(FloatArray &answer) const { answer = plasticStrainVector; }
    void giveTempPlasticStrainVector(FloatArray &answer) const { answer = tempPlasticStrainVector; }
    void giveStrainSpaceHardeningVars(FloatArray &answer) const
    { answer = strainSpaceHardeningVarsVector; }
    void givetempStrainSpaceHardeningVarsVector(FloatArray &answer) const
    { answer = tempStrainSpaceHardeningVarsVector; }

    void         letPlasticStrainVectorBe(const FloatArray &v)
    { plasticStrainVector = v; }
    void         letTempPlasticStrainVectorBe(const FloatArray &v)
    { tempPlasticStrainVector = v; }
    void         letTempStrainSpaceHardeningVarsVectorBe(const FloatArray &v)
    { tempStrainSpaceHardeningVarsVector = v; }
    void         letStrainSpaceHardeningVarsVectorBe(const FloatArray &v)
    { strainSpaceHardeningVarsVector = v; }

    int giveStateFlag() { return state_flag; }
    int giveTempStateFlag() { return temp_state_flag; }
    double givePlasticConsistencyPrameter() { return gamma; }
    double giveTempPlasticConsistencyPrameter() { return temp_gamma; }
    void letTempPlasticConsistencyPrameterBe(double v) { gamma = v; }

    void letTempStateFlagBe(int v) { temp_state_flag = v; }


    // definition
    const char *giveClassName() const { return "PlasticMaterialStatus"; }
    classType             giveClassID() const
    { return PerfectlyPlasticMaterialStatusClass; }
};

class PlasticMaterial : public StructuralMaterial
{
    /*
     * This class implements a general plastic material.
     * It is assumed to be a base class for many material
     * models based on different yield conditions and hardening laws.
     * Sress return mapping algorithm  is based on general
     * return mapping algorithm, with following assumptions
     * - associative flow rule
     * - general hardening law is supported
     * - Kunt-Tucker conditions apply.
     *
     * DESCRIPTION
     *
     * TASK
     *
     */

protected:
    // add common (same for every gauss point) material parameters here

    /// Reference to bulk (undamaged) material
    LinearElasticMaterial *linearElasticMaterial;

public:

    PlasticMaterial(int n, Domain *d);
    ~PlasticMaterial();

    // identification and auxiliary functions
    int hasNonLinearBehaviour()   { return 1; }
    int hasMaterialModeCapability(MaterialMode mode);
    const char *giveClassName() const { return "PlasticMaterial"; }
    classType giveClassID()         const { return PerfectlyPlasticMaterialClass; }

    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    /// Returns reference to undamaged (bulk) material
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    // non-standard - returns time independent material constant
    // double   give (int) ;

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix & answer,
                                               MatResponseForm, MatResponseMode,
                                               GaussPoint * gp,
                                               TimeStep * atTime);


    void giveRealStressVector(FloatArray & answer,  MatResponseForm, GaussPoint *,
                              const FloatArray &, TimeStep *);

    /**
     * Returns the integration point corresponding value in Reduced form.
     * @param answer contain corresponding ip value, zero sized if not available
     * @param aGaussPoint integration point
     * @param type determines the type of internal variable
     * @param type determines the type of internal variable
     * @returns nonzero if ok, zero if var not supported
     */
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    /**
     * Returns the mask of reduced indexes of Internal Variable component .
     * @param answer mask of Full VectorSize, with components beeing the indexes to reduced form vectors.
     * @param type determines the internal variable requested (physical meaning)
     * @returns nonzero if ok or error is generated for unknown mat mode.
     */
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);


    /**
     * Returns the type of internal variable (scalar, vector, tensor,...).
     * @param type determines the type of internal variable
     * @returns type of internal variable
     */
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    /**
     * Returns the corresponding integration point  value size in Reduced form.
     * @param type determines the type of internal variable
     * @returns var size, zero if var not supported
     */
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);

#ifdef __OOFEG
#endif

    MaterialStatus *CreateStatus(GaussPoint *gp) const;


protected:

    // add here some auxiliary functions if needed
    FloatArray *ComputeGradientVector(GaussPoint *gp, FloatArray *fullStressVector,
                                      FloatArray *fullStressSpaceHardeningVars);
    FloatArray *ComputeResidualVector(GaussPoint *gp, double Gamma,
                                      FloatArray *plasticStrainVectorR,
                                      FloatArray *strainSpaceHardeningVariables,
                                      FloatArray *gradientVectorR);
    virtual void giveConsistentStiffnessMatrix(FloatMatrix & answer,
                                               MatResponseForm,
                                               MatResponseMode,
                                               GaussPoint * gp,
                                               TimeStep * atTime);

    void  computeConsistentModuli(FloatMatrix &answer,
                                  GaussPoint *gp, FloatMatrix &elasticModuliInverse,
                                  FloatMatrix &hardeningModuliInverse,
                                  double Gamma, const FloatArray &fullStressVector,
                                  const FloatArray &fullStressSpaceHardeningVars);
    void  computeDiagModuli(FloatMatrix &answer,
                            GaussPoint *gp, FloatMatrix &elasticModuliInverse,
                            FloatMatrix &hardeningModuliInverse);

    virtual FloatArray *ComputeStressSpaceHardeningVars(GaussPoint *gp,
                                                        FloatArray *strainSpaceHardeningVariables)
    { return NULL; }
    virtual double computeYieldValueAt(GaussPoint *gp, FloatArray *stressVector,
                                       FloatArray *stressSpaceHardeningVars) { return 0.; }
    virtual void   computeHardeningReducedModuli(FloatMatrix &answer,
                                                 GaussPoint *gp,
                                                 FloatArray *strainSpaceHardeningVariables,
                                                 TimeStep *atTime) = 0;
    virtual FloatArray *ComputeStressGradient(GaussPoint *gp, FloatArray *stressVector,
                                              FloatArray *stressSpaceHardeningVars) { return NULL; }
    virtual FloatArray *ComputeStressSpaceHardeningVarsReducedGradient(GaussPoint *gp,
                                                                       FloatArray *stressVector,
                                                                       FloatArray *stressSpaceHardeningVars)
    { return NULL; }
    virtual int hasHardening() { return 0; }
    virtual void  computeReducedGradientMatrix(FloatMatrix &answer,
                                               GaussPoint *gp,
                                               const FloatArray &stressVector,
                                               const FloatArray &stressSpaceHardeningVars) = 0;

    virtual void computeTrialStressIncrement(FloatArray &answer, GaussPoint *gp,
                                             const FloatArray &strainIncrement, TimeStep *atTime);
    virtual void computeReducedElasticModuli(FloatMatrix &answer, GaussPoint *gp,
                                             TimeStep *atTime);
    virtual void compute3dElasticModuli(FloatMatrix &answer, GaussPoint *gp,
                                        TimeStep *atTime) = 0;

    // auxiliary functions
    virtual int         giveSizeOfFullHardeningVarsVector()  { return 0; }
    virtual int         giveSizeOfReducedHardeningVarsVector(GaussPoint *)  { return 0; }

    friend class PlasticMaterialStatus;



    // next functions overloaded rom structural material level
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
#endif // plasticmaterial_h
