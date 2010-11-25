/* $Header: /home/cvs/bp/oofem/sm/src/perfectlyplasticmaterial.h,v 1.5 2003/04/06 14:08:31 bp Exp $ */
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

//   ********************************************************
//   *** CLASS LINEAR ELACSTIC PERFECTLY PLASTIC MATERIAL ***
//   ********************************************************

#ifndef perfectlyplasticmaterial_h
#define perfectlyplasticmaterial_h

#include "structuralmaterial.h"
#include "linearelasticmaterial.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"

#include "structuralms.h"

namespace oofem {
class GaussPoint;


class PerfectlyPlasticMaterialStatus : public StructuralMaterialStatus
{
    /*
     * This class implements associated Material Status to PerfectlyPlasticMaterial.
     * It is atribute of matStatusDictionary at every GaussPoint, for which this material
     * is active.
     * DESCRIPTION:
     * This class contains only yield_flag to inform, whether
     * gp is yielding or is in elastic state.
     * temp_yield_flag is used to indicate yield_flag during a process when we try
     * to find a equilibrium state. After This State is reached we update the Yield_Flag
     * accordingly.
     * When lastly reched state is eqilibrium state we call updateYourself(), so
     * yield_flag containing last reched status of yielding is
     * copied into yield_flag indicating yield status at previously
     * reched equilibrium state .
     * TASK:
     *   setting the yield_flag
     * Printing to stdin
     * Printing to file
     */

protected:

    FloatArray plasticStrainVector; // reduced form to save memory
    FloatArray plasticStrainIncrementVector;
    int yield_flag;
    int temp_yield_flag;

public:
    PerfectlyPlasticMaterialStatus(int n, Domain *d, GaussPoint *g);
    ~PerfectlyPlasticMaterialStatus();
    void   printOutputAt(FILE *file, TimeStep *tStep);

    int    setTempYieldFlag(int i) { return temp_yield_flag = i; }
    int    giveTempYieldFlag()     { return temp_yield_flag; }
    int    giveYieldFlag() { return yield_flag; }
    // saves current context(state) into stream
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual void initTempStatus();
    void   updateYourself(TimeStep *tNow);

    void givePlasticStrainVector(FloatArray &answer) const { answer =  plasticStrainVector; }
    void givePlasticStrainIncrementVector(FloatArray &answer) const { answer = plasticStrainIncrementVector; }
    void         letPlasticStrainVectorBe(FloatArray &v)
    { plasticStrainVector = v; }
    void         letPlasticStrainIncrementVectorBe(FloatArray &v)
    { plasticStrainIncrementVector = v; }

    // definition
    const char *giveClassName() const { return "PerfectlyPlasticMaterialStatus"; }
    classType             giveClassID() const
    { return PerfectlyPlasticMaterialStatusClass; }
};

class PerfectlyPlasticMaterial : public StructuralMaterial
{
    /*
     * This class implements a perfectly plastic material in a finite element problem. A material
     * is an attribute of a domain. It is usually also attribute of many elements.
     * DESCRIPTION
     * This class is mainly abstract, no new methods or data added. This class is
     * mainly intended to be a base class for perfectly plastic materials.
     * (for perfectly plastic material a linear elastic behaviour is asumed
     * before yielding surface is beeing reached, and also unloading is assumed to
     * be linear elastic).
     * This class takes care about possible failure and fracture and hardening,
     * even if it is not required
     * for perfectly-plastic material, but this class is assumed to be a parent
     * class for more sophisticated materials where these phenomena should appear.
     * TASK
     * - Returning standard material stiffness and flexibility marices for 3d-case.
     * according to current state determined by using data stored
     * in Gausspoint.
     * - Returning a material property (method 'give'). Only for non-standard elements.
     * - Returning real stress state vector(tensor) at gauss point for 3d - case.
     */

protected:

    int yieldCriteria;
    int loadingCriteria;
    LinearElasticMaterial *linearElasticMaterial;
    // ComputePlasticStiffnessAt for description
public:

    PerfectlyPlasticMaterial(int n, Domain *d) : StructuralMaterial(n, d)
    {
        yieldCriteria = 0;
        loadingCriteria = 0;
        linearElasticMaterial = NULL;
    }

    ~PerfectlyPlasticMaterial() { delete linearElasticMaterial; }

    void giveRealStressVector(FloatArray & answer, MatResponseForm, GaussPoint *,
                              const FloatArray &, TimeStep *);

    // uptates MatStatus to newly reched (equilibrium) state
    virtual void updateYourself(GaussPoint *gp, TimeStep *);

    virtual void  updateIfFailure(GaussPoint *gp,
                                  FloatArray *stressVector3d,
                                  FloatArray *PlasticStrainVector3d) { }
    // updateIfFailure .. can be moved to parent classes
    //
    // identification and auxiliary functions
    int hasNonLinearBehaviour()   { return 1; }
    int  hasMaterialModeCapability(MaterialMode mode);
    const char *giveClassName() const { return "PerfectlyPlasticMaterial"; }
    classType giveClassID()         const { return PerfectlyPlasticMaterialClass; }
    IRResultType initializeFrom(InputRecord *ir);
    void     printYourself();
    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);


    // non-standard  - returns time independent material constant
    double   give(int, GaussPoint *);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix & answer,
                                               MatResponseForm, MatResponseMode,
                                               GaussPoint * gp,
                                               TimeStep * atTime);
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

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
    // two functions used to initialize and updating temporary variables in
    // gps status. These variables are used to controll process, when
    // we try to find equlibrium state

    virtual void giveMaterialStiffnessMatrix(FloatMatrix & answer,
                                             MatResponseMode, GaussPoint * gp,
                                             TimeStep * atTime);
    virtual void giveEffectiveMaterialStiffnessMatrix(FloatMatrix & answer,
                                                      MatResponseMode, GaussPoint * gp,
                                                      TimeStep * atTime);
    // Give3dMaterialStiffnessMatrix should return 3d material stiffness matrix
    // taking into account possible failure or fracture of material

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
    virtual void give3dShellLayerStiffMtrx(FloatMatrix & answer,
                                           MatResponseForm, MatResponseMode,
                                           GaussPoint * gp,
                                           TimeStep * atTime);




    // Give3dMaterialStiffnessMatrix should return 3d material stiffness matrix
    // taking into account possible failure or fracture of material
    //
    void computeTrialStressIncrement(FloatArray &answer, GaussPoint *gp,
                                     const FloatArray &strainIncrement, TimeStep *);
    void computePlasticStiffnessAt(FloatMatrix &answer,
                                   GaussPoint *gp,
                                   FloatArray *currentStressVector,
                                   FloatArray *currentPlasticStrainVector,
                                   FloatArray *strainIncrement3d,
                                   TimeStep *atTime,
                                   double &lambda);
    FloatArray *GiveStressCorrectionBackToYieldSurface(GaussPoint *gp,
                                                       FloatArray *stressVector3d,
                                                       FloatArray *plasticVector3d);
    //
    // yield(YC-like functions) and loading(LC-like functions) criteria specific section
    //

    virtual double      computeYCValueAt(GaussPoint *, FloatArray *, FloatArray *)
    { return -1.0; }
    virtual FloatArray *GiveYCStressGradient(GaussPoint *, FloatArray *, FloatArray *)
    { return NULL; }
    virtual FloatArray *GiveLCStressGradient(GaussPoint *, FloatArray *, FloatArray *)
    { return NULL; }
    virtual FloatArray *GiveYCPlasticStrainGradient(GaussPoint *, FloatArray *, FloatArray *)
    { return NULL; }
    virtual FloatArray *GiveLCPlasticStrainGradient(GaussPoint *, FloatArray *, FloatArray *)
    { return NULL; }
    // virtual FloatArray* GiveHardeningGradient (GaussPoint*,FloatArray*,FloatArray*)
    // {return NULL;}
    virtual void        updateTempYC(GaussPoint *, FloatArray *, FloatArray *) { }
    virtual void        updateTempLC(GaussPoint *, FloatArray *, FloatArray *) { }
    // update during computation

    // virtual int  updateYieldStatus (GaussPoint* gp, FloatArray* strainIncrementIn3d);
};
} // end namespace oofem
#endif // perfectlyplasticmaterial_h
