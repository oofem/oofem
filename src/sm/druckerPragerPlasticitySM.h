/* $Header: /home/cvs/bp/oofem/sm/src/druckerPragerPlasticitySM.h,v 1.1 2003/04/06 14:08:30 bp Exp $ */
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
/*
 *  Author: Simon Rolshoven
 */

// file : DruckerPragerPlasticitySM.h
//   ********************************************************************
//   *** CLASS DRUCKER-PRAGER PLASTICITY STRUCTURAL MATERIAL STATUS   ***
//   ********************************************************************

#ifndef druckerpragerplasticitysm_h
#define druckerpragerplasticitysm_h

#include "flotarry.h"
#include "flotmtrx.h"
#include "cltypes.h"
#include "structuralms.h"
#include "strainvector.h"
#include "stressvector.h"
#include "structuralmaterial.h"
#include "isolinearelasticmaterial.h"


class DruckerPragerPlasticitySMStatus : public StructuralMaterialStatus
{
    /*
     *       This class implements the material status associated to DruckerPragerPlasticitySM.
     *
     *       At every GaussPoint for which this material is active, it is an attribute
     *       of matStatusDictionary.
     *
     *       DESCRIPTION:
     *       Idea used there is that we have variables describing:
     *
     *       temp-variables:
     *       at the beginning of an iteration for material equilibrium,
     *       these variables carry the previous equilibrium state.
     *       When material equilibrium is achieved,
     *       they carry the values associated to this state. If it corresponds
     *       to a new equilibrium, i.e. if convergence in the force balance is satisfied,
     *       they are copied to the
     *
     *       non-temp-variables:
     *       they always carry the state corresponding to the latest global equilibrium.
     */


public:
    /// values of history variable state_flag
    enum state_flag_values { DP_Elastic, DP_Unloading, DP_Yielding, DP_Vertex };

protected:
    /// volumetric plastic strain
    double volumetricPlasticStrain;
    double tempVolumetricPlasticStrain;

    /// deviator of plastic strain
    StrainVector plasticStrainDeviator;
    StrainVector tempPlasticStrainDeviator;

    /// hardening variable
    double kappa;
    double tempKappa;

    /// Indicates the state (i.e. elastic, yielding, vertex, unloading) of the Gauss point
    int state_flag;
    int temp_state_flag;

public:
    /// constructor
    DruckerPragerPlasticitySMStatus(int, Domain *, GaussPoint *);

    /// destructor
    ~DruckerPragerPlasticitySMStatus();

    void  initTempStatus();
    void  updateYourself(TimeStep *);
    void  printOutputAt(FILE *file, TimeStep *tStep);
    contextIOResultType  saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType  restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    const char *giveClassName() const
    { return "DruckerPragerPlasticitySMStatus"; }
    classType  giveClassID() const
    { return DruckerPragerPlasticitySMStatusClass; }

    // Inline functions for access to state variables
    // give:
    // Functions used to access the value of a internal variable.
    /**
     *      Get the full plastic strain vector from the material status.
     *      @param answer plastic strain vector.
     */
    void  giveFullPlasticStrainVector(StrainVector &answer) const
    {
        StrainVector plasticStrainVector(_3dMat);
        plasticStrainDeviator.computeDeviatoricVolumetricSum(plasticStrainVector,
                                                             volumetricPlasticStrain);
        plasticStrainVector.convertToFullForm(answer);
    }
    /**
     *      Get the plastic strain deviator from the material status.
     *      @param answer plastic strain deviator.
     */
    void  givePlasticStrainDeviator(StrainVector &answer) const
    { answer = plasticStrainDeviator; }
    /**
     *      Get the volumetric strain deviator from the material status.
     *      @returns volumetric plastic strain
     */
    double giveVolumetricPlasticStrain() const
    { return volumetricPlasticStrain; }
    /**
     *      Get the hardening variable from the material status.
     *      @returns hardening variable kappa
     */
    double  giveKappa() const
    { return kappa; }
    /**
     *      Get the state flag from the material status.
     *      @returns state flag (i.e. elastic, unloading, yielding, vertex case yielding)
     */
    int giveStateFlag() const
    { return state_flag; }

    // giveTemp:
    // Functions used to access the temp variables.
    /**
     *      Get the temp value of the full plastic strain vector from the material status.
     *      @param answer temp value of plastic strain vector.
     */
    void  giveTempPlasticStrainVector(StrainVector &answer) const
    {
        tempPlasticStrainDeviator.computeDeviatoricVolumetricSum(answer,
                                                                 tempVolumetricPlasticStrain);
    }
    /**
     *      Get the temp value of the plastic strain deviator from the material status.
     *      @param answer temp value of plastic strain deviator.
     */
    void  giveTempPlasticStrainDeviator(StrainVector &answer) const
    { answer = tempPlasticStrainDeviator; }
    /**
     *      Get the temp value of the volumetric strain deviator from the material status.
     *      @returns temp value of volumetric plastic strain
     */
    double giveTempVolumetricPlasticStrain() const
    { return tempVolumetricPlasticStrain; }
    /**
     *      Get the temp value of the hardening variable from the material status.
     *      @returns temp value ofhardening variable kappa
     */
    double  giveTempKappa() const
    { return tempKappa; }
    /**
     *      Get the temp value of the state flag from the material status.
     *      @returns the temp value of the state flag (i.e. elastic, unloading, yielding, vertex case yielding)
     */
    int  giveTempStateFlag() const
    { return temp_state_flag; }

    // letTemp...be :
    // Functions used by the material to assign a new value to a temp variable.
    /**
     *      Assign the temp value of deviatoric plastic strain.
     *      @v new temp value of deviatoric plastic strain
     */
    void  letTempPlasticStrainDeviatorBe(const StrainVector &v)
    { tempPlasticStrainDeviator = v; }
    /**
     *      Assign the temp value of volumetric plastic strain.
     *      @v new temp value of volumetric plastic strain
     */
    void  letTempVolumetricPlasticStrainBe(const double v)
    { tempVolumetricPlasticStrain = v; }
    /**
     *      Assign the temp value of the hardening variable.
     *      @v new temp value of the hardening variable
     */
    void  letTempKappaBe(const double v)
    { tempKappa = v; }
    /**
     *      Assign the temp value of the state flag.
     *      @v new temp value of the state flag (i.e. elastic, unloading, yielding, vertex case yielding)
     */
    void  letTempStateFlagBe(const int v)
    { temp_state_flag = v; }
};

//   *************************************************************
//   *** CLASS DRUCKER-PRAGER PLASTICITY STRUCTURAL MATERIAL   ***
//   *************************************************************



/**
 * This class implements a (local) nonassociated plasticity model based on the Drucker-Prager yield criterion with hardening and softening.
 */
class DruckerPragerPlasticitySM : public StructuralMaterial
{
public:
    /// values of history variable state_flag
    enum state_flag_values { DP_Elastic, DP_Unloading, DP_Yielding, DP_Vertex };

protected:
    /**
     *      Controls the hardening function in the yield stress:
     *      1: linear hardening/softening with cutoff at zero stress.
     *      2: exponential hardening/softening to limitYieldStress.
     */
    int hardeningType;
    /// parameter of the exponential laws
    double kappaC;
    /// hardening modulus normalized with the elastic modulus, parameter of the linear hardening/softening law
    double hardeningModulus;
    /// parameter of the exponential hardening law
    double limitYieldStress;
    /// parameter of all three laws, this is the initial value of the yield stress in pure shear.
    double initialYieldStress;
    /// friction coefficient, parameter of the yield criterion
    double alpha;
    /// dilatancy coefficient, parameter of the flow rule
    double alphaPsi;
    /// scalar factor between rate of plastic multiplier and rate of hardening variable
    double kFactor;

    /// pointer for linear elastic material
    IsotropicLinearElasticMaterial *LEMaterial;
    /// elastic material constants stored here for convenient access
    double eM, gM, kM, nu;

    /// hardening variable
    double kappa;
    double tempKappa;

    /// volumetric stress, i.e. 1/3 of the trace of sigma
    double volumetricStress;
    /// volumetric part of elastic trial strain
    double volumetricElasticTrialStrain;
    /// J2-invariant of elastic trial stress
    double trialStressJTwo;
    /// deviatoric part of stress
    StressVector stressDeviator;

    /// yield tolerance
    double yieldTol;

    /// Maximum number of iterations for stress return.
    int newtonIter;

public:
    /// constructor
    DruckerPragerPlasticitySM(int n, Domain *d);
    /// destructor
    ~DruckerPragerPlasticitySM();
    IRResultType initializeFrom(InputRecord *ir);
    int  hasMaterialModeCapability(MaterialMode mMode);
    const char *giveClassName() const
    { return "DruckerPragerPlasticitySM"; }
    classType giveClassID() const
    { return DruckerPragerPlasticitySMClass; }
    void  giveRealStressVector(FloatArray &answer,
                               MatResponseForm form,
                               GaussPoint *gp,
                               const FloatArray &strainVector,
                               TimeStep *atTime);


    /// Perform local stress return and compute history variables.
    /**
     *      Perform a standard local stress return using the function computeYieldValue at the specified Gauss point. This function computes the history variables, i.e. the temp variable of plastic strain, hardening variable, state flag, and the temp stress, and stores them in the temp-status.
     *      @param gp Gauss point
     *      @param strain strain vector of this Gauss point
     */
    void  performLocalStressReturn(GaussPoint *gp,
                                   const StrainVector &strain);
    /// Check for vertex stress return.
    /**
     *      Check if the trial stress state falls within the vertex region.
     *      @returns true for vertex case and false if regular stress return has to be used.
     */
    bool checkForVertexCase();
    /// Perform stress return for regular case.
    /**
     *      Perform stress return for regular case, i.e. if the trial stress state does not lie within the vertex region.
     */
    void  performRegularReturn();
    /// Perform stress return for vertex case.
    /**
     *      Perform stress return for vertex case, i.e. if the trial stress state lies within the vertex region.
     */
    void  performVertexReturn();
    /// Compute yield value.
    /**
     *      Compute the yield value based on stress and hardening variable.
     *      @param volumetricStress 1/3 of trace of sigma
     *      @param JTwo second deviatoric invariant
     *      @param kappa hardening variable
     *      @returns yield value
     */
    double  computeYieldValue(const double volumetricStress,
                              const double JTwo,
                              const double kappa) const;
    /// Compute current yield stress in pure shear
    /**
     *      Compute the current yield stress in pure shear of the Drucker-Prager model according to the used hardening law. The yield stress is tauY in f(sigma, kappa) = F(sigma) - tauY(kappa).
     *      @param kappa hardening variable
     *      @returns yield stress in pure shear
     */
    virtual double  computeYieldStressInShear(const double kappa) const;

    /// Compute derivative of yield stress with respect to kappa
    /**
     *      Compute derivative of yield stress with respect to the hardening variable kappa.
     *      @param kappa hardening variable
     *      @returns derivative of yield stress with respect to kappa
     */
    virtual double  computeYieldStressPrime(const double kappa) const;
    void give3dMaterialStiffnessMatrix(FloatMatrix &, MatResponseForm,
                                       MatResponseMode, GaussPoint *, TimeStep *);


    /**
     *      Compute and give back algorithmic stiffness matrix for the regular case (no vertex).
     *      @param answer consistent stiffness matrix
     *      @param form material response form
     *      @param mode material reponse mode
     *      @param gp Gauss point
     *      @param atTime time step
     */
    void giveRegAlgorithmicStiffMatrix(FloatMatrix &answer,
                                       MatResponseForm form,
                                       MatResponseMode mode,
                                       GaussPoint *gp,
                                       TimeStep *atTime);
    /**
     *      Compute consistent stiffness matrix for the vertex case.
     *      @param answer consistent stiffness matrix
     *      @param form material response form
     *      @param mode material reponse mode
     *      @param gp Gauss point
     *      @param atTime time step
     */
    void giveVertexAlgorithmicStiffMatrix(FloatMatrix &answer,
                                          MatResponseForm form,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *atTime);

    /**
     *      Returns the value of a state variable at the specified integration point in reduced form.
     *      @param answer contains corresponding ip value, zero sized if not available
     *      @param gp integration point
     *      @param type determines the type of internal variable
     *      @param atTime time step
     *      @returns nonzero if ok, zero if var not supported
     */
    int giveIPValue(FloatArray &answer,
                    GaussPoint *gp,
                    InternalStateType type,
                    TimeStep *atTime);

    /**
     *      Returns the size of a state variable at the specified integration point in reduced form.
     *      @param type determines the type of internal variable
     *      @param gp Gauss point
     *      @returns var size, zero if var not supported
     */
    int giveIPValueSize(InternalStateType type,
                        GaussPoint *gp);

    /**
     *      Returns the mask of reduced indexes of Internal Variable component.
     *      @param answer mask of Full VectorSize, with components being the indexes to reduced form vectors.
     *      @param type determines the internal variable requested (physical meaning)
     *      @param mmode material mode
     *      @returns nonzero if ok or error is generated for unknown mat mode.
     */
    int giveIntVarCompFullIndx(IntArray &answer,
                               InternalStateType type,
                               MaterialMode mmode);

    /**
     *      Returns the type of internal variable (scalar, vector, tensor,...).
     *      @param type determines the type of internal variable
     *      @returns type of internal variable
     */
    InternalStateValueType giveIPValueType(InternalStateType type);
    /**
     * Returns true if stiffness matrix of receiver is symmetric
     * Default implementation returns true.
     */
    virtual bool isCharacteristicMtrxSymmetric(MatResponseMode rMode) { return false; }

    /**
     * Returns a vector of coefficients of thermal dilatation in direction
     * of each material principal (local) axis.
     * @param answer vector of thermal dilatation coefficients
     * @param gp integration point
     * @param tStep time step (most models are able to respond only when atTime is current time step)
     */
    void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *atTime)
    {
        LEMaterial->giveThermalDilatationVector(answer, gp, atTime);
    }

    MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
};

#endif // druckerpragerplasticitysm_h

