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
 *               Copyright (C) 1993 - 2009   Borek Patzak
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

//   ***********************************
//   *** CLASS Rheologic Chain Model ***
//   ***********************************

#ifndef rheochm_h
#define rheochm_h

#include "femcmpnn.h"
#include "structuralmaterial.h"
#include "linearelasticmaterial.h"
#include "flotarry.h"
#include "flotmtrx.h"

#include "matconst.h"
#include "structuralelement.h"
#include "structuralms.h"

namespace oofem {

#define MNC_NPOINTS 30
#define TIME_DIFF   1.e-10


class RheoChainMaterialStatus : public StructuralMaterialStatus
{
    /*
     * This class implements associated Material Status to RheoChainMaterial.
     * It is an attribute of matStatusDictionary at every GaussPoint 
     * for which this material
     * is active. Isotropic linear viscoelastic material is assumed.
     * DESCRIPTION:
     * Idea used there is that we have variables
     * describing:
     * 1) state at previous equilibrium state (variables without temp)
     * 2) state during searching new equilibrium (variables with temp)
     * when we start search new state from previous equilibrium one we copy
     * non-tem variables into temp ones. And after we reach new equilibrium
     * (now decribed by temp variables) we copy temp-var into non-temp ones
     * (see function updateYourself).
     *
     */

protected:
    /// number of units in the chain 
    int nUnits;
    /// hidden (internal) variables, the meaning of which depends on the type of chain
    FloatArray **hiddenVars;
    /* total shrinkage strain (needed only when the shrinkage evolution 
       is described in the incremental form)
    */
    FloatArray shrinkageStrain;

public:
    RheoChainMaterialStatus(int n, Domain *d, GaussPoint *g, int nunits);
    ~RheoChainMaterialStatus();
    void printOutputAt(FILE *file, TimeStep *tStep);

    FloatArray *giveHiddenVarsVector(int i) { return hiddenVars [ i - 1 ]; }
    FloatArray *letHiddenVarsVectorBe(int i, FloatArray *);

    FloatArray *giveShrinkageStrainVector() { return & shrinkageStrain; }
    void        setShrinkageStrainVector(const FloatArray &src) { shrinkageStrain = src; }

    /// initialize the status
    virtual void initTempStatus();
    /// update after new equilibrium state reached
    virtual void updateYourself(TimeStep *); 

    /// save current context (state) into a stream
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /// restore current context (state) from a stream
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    // definition
    const char *giveClassName() const { return "RheoChainMaterialStatus"; }
    classType             giveClassID() const
    { return RheoChainMaterialStatusClass; }
};


//=================================================================================

class RheoChainMaterial : public StructuralMaterial
{
    /*
     * This class implements a rheologic chain model 
     * describing a viscoelastic material.
     * It serves as the parent class for Maxwell and Kelvin chains.
     *
     */

protected:

    /// number of (Maxwell or Kelvin) units in the rheologic chain 
    int nUnits;
    /// physical age of the material at simulation time = 0
    double relMatAge;
    /// Poisson's ratio (assumed to be constant, unaffected by creep)
    double nu;
    /// incremental modulus (links the stress increment to the strain increment)
    double Einc;
    /// time for which the partial moduli of individual units have been evaluated
    double EparValTime;

    /// time from which the model should give a good approximation. Optional field. Default value is 0.1 [day].
    double begOfTimeOfInterest; // local one or taken from e-model
    /// time (age???) up to which the model should give a good approximation 
    double endOfTimeOfInterest; // local one or taken from e-model
    // associated linearElasticMaterial, with E = 1;
    LinearElasticMaterial *linearElasticMaterial;
    /// partial moduli of individual units
    FloatArray EparVal;
    //FloatArray relaxationTimes;
    /// characteristic times of individual units (relaxation or retardation times)
    FloatArray charTimes;
    /// times at which the errors are evaluated if the least-square method is used
    FloatArray discreteTimeScale;

    /** scaling factor transforming the simulation time units into days
      * (gives the number of simulation time units in one day,
      *  e.g. 86400 if the simulation works with seconds as the time units) */
    double timeFactor;


public:
    RheoChainMaterial(int n, Domain *d);
    ~RheoChainMaterial();

    /// evaluation of standard material stiffness matrices
    virtual void  giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime);

    /// evaluation of stress from a given strain and internal variables
    virtual void giveRealStressVector(FloatArray & answer,  MatResponseForm, GaussPoint *,
                                      const FloatArray &, TimeStep *);

    // returns a FloatArray(3) of coefficients of thermal dilation in direction
    // of each (local) axis given by the principal axes of the material
    //
    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *, TimeStep *)
    { answer.resize(0); }

    /// evaluation of the incremental modulus
    virtual double giveEModulus(GaussPoint *gp, TimeStep *atTime){return 0.0;}

    /// evaluation of the moduli of individual units
    virtual void computeCharCoefficients(FloatArray &answer, GaussPoint *gp, double){};

    /// update of MatStatus to the newly reached (equilibrium) state
    virtual void updateYourself(GaussPoint *gp, TimeStep *){};

    // identification and auxiliary functions
    virtual int hasNonLinearBehaviour()   { return 0; }
    virtual int hasMaterialModeCapability(MaterialMode mode);
    const char *giveClassName() const { return "RheoChainMaterial"; }
    classType giveClassID()         const { return RheoChainMaterialClass; }
    IRResultType initializeFrom(InputRecord *ir);
    void     printYourself();

    // store & restore context functions
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);


    /// evaluation of material stiffness matrix in 3D
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix & answer,
                                               MatResponseForm, MatResponseMode,
                                               GaussPoint * gp,
                                               TimeStep * atTime);

    /**
     * Computes, for the given integration point, 
     * the strain vector induced by stress-independent
     * internal processes in the material.
     * Default implementation takes into account only temperature-induced strains.
     * Overloaded to take into account shrinkage and creep strains.
     * @param answer returned strain vector
     * @param gp integration point
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     * @param determines response mode (Total or incremental)
     */
    virtual void computeStressIndependentStrainVector(FloatArray &answer,
                                                      GaussPoint *gp, TimeStep *stepN, ValueModeType mode);
    /**
     * Computes, for the given integration point, 
     * the strain vector induced by stress-independent shrinkage 
     * @param answer returned strain vector
     * @param form material response form
     * @param gp integration point
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     * @param determines response mode (Total or incremental)
     */
    virtual void  giveShrinkageStrainVector(FloatArray &answer,
                                            MatResponseForm form,
                                            GaussPoint *gp,
                                            TimeStep *atTime,
                                            ValueModeType mode)
    { answer.resize(0); }

    // Note: must take LoadResponseMode into account
    /**
     * Computes, for the given integration point, 
     * the strain vector induced by the stress history (typically creep strain)
     * @param answer computed strains
     * @param form material response form
     * @param gp integration point
     * @param atTime time step (most models are able to respond only when atTime is the current time step)
     * @param mode determines response mode
     */
    virtual void  giveEigenStrainVector(FloatArray &answer, MatResponseForm form,
                                        GaussPoint *gp, TimeStep *atTime, ValueModeType mode){};

#ifdef __OOFEG
#endif

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;



protected:

    /** if only incremental shrinkage strain formulation is provided, then total shrinkage strain must be tracked
     * in status in order to be able to compute total value. */
    virtual int  hasIncrementalShrinkageFormulation() { return 0; }

    void         generateLogTimeScale(FloatArray &answer, double from, double to, int nsteps,
                                      int fromIncluded = 0);
    const FloatArray &giveDiscreteTimes();
    /// evaluation of the creep compliance function
    virtual double computeCreepFunction(GaussPoint *gp, double ofAge, double atTime) = 0;

    /// evaluation of the relaxation function at given times
    void computeDiscreteRelaxationFunction(FloatArray &answer, GaussPoint *gp,
                                                   const FloatArray &atTimes,
                                                   double t0, double tr);

    /// evaluation of elastic compliance matrix for unit Young's modulus
    void giveUnitComplianceMatrix(FloatMatrix & answer, MatResponseForm,
                                   GaussPoint * gp, TimeStep * tStep);
    /// evaluation of elastic stiffness matrix for unit Young's modulus
    void giveUnitStiffnessMatrix(FloatMatrix & answer,
                                      MatResponseForm, GaussPoint * gp, TimeStep * tStep);

    /// update of partial moduli of individual chain units
    void         updateEparModuli(GaussPoint *gp, double atTime);

    /// access to partial modulus of a given unit
    double       giveEparModulus(int iChain);

    /// evaluation of characteristic times
    virtual void    computeCharTimes();

    /// access to the characteristic time of a given unit
    double       giveCharTime(int);

    /// exponent to be used with the char time of a given unit, usually = 1.0
    virtual double giveCharTimeExponent(int i) { return 1.0; }

    /// access to the underlying linear elastic material with unit Young's modulus
    LinearElasticMaterial *giveLinearElasticMaterial();

    /// access to the time up to which the response should be accurate
    double giveEndOfTimeOfInterest();

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
    /**
     * Computes, for the given integration point, 
     * the strain vector induced by stress-independent
     * internal processes in the material.
     * Takes into account only temperature and shrinkage-induced strains.
     * @param answer returned strain vector
     * @param gp integration point
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     * @param determines response mode (Total or incremental)
     */
    void computeTrueStressIndependentStrainVector(FloatArray &answer, GaussPoint *gp,
                                                  TimeStep *stepN, ValueModeType mode);
};

} // end namespace oofem
#endif // rheochm_h
