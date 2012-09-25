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

#ifndef rheochm_h
#define rheochm_h

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

/**
 * This class implements associated Material Status to RheoChainMaterial.
 */
class RheoChainMaterialStatus : public StructuralMaterialStatus
{
protected:
    /// Number of units in the chain.
    int nUnits;
    /// Hidden (internal) variables, the meaning of which depends on the type of chain.
    FloatArray **hiddenVars;
    /**
     * Total shrinkage strain (needed only when the shrinkage evolution
     * is described in the incremental form).
     */
    FloatArray shrinkageStrain;

public:
    RheoChainMaterialStatus(int n, Domain *d, GaussPoint *g, int nunits);
    virtual ~RheoChainMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    FloatArray *giveHiddenVarsVector(int i) { return hiddenVars [ i - 1 ]; }
    FloatArray *letHiddenVarsVectorBe(int i, FloatArray *);

    FloatArray *giveShrinkageStrainVector() { return & shrinkageStrain; }
    void setShrinkageStrainVector(const FloatArray &src) { shrinkageStrain = src; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    // definition
    virtual const char *giveClassName() const { return "RheoChainMaterialStatus"; }
    virtual classType giveClassID() const { return RheoChainMaterialStatusClass; }
};


/**
 * This class implements a rheologic chain model
 * describing a viscoelastic material.
 * It serves as the parent class for Maxwell and Kelvin chains.
 */
class RheoChainMaterial : public StructuralMaterial
{
protected:
    /// Number of (Maxwell or Kelvin) units in the rheologic chain.
    int nUnits;
    /// Physical age of the material at simulation time = 0.
    double relMatAge;
    /// Poisson's ratio (assumed to be constant, unaffected by creep).
    double nu;
    /// Time for which the partial moduli of individual units have been evaluated.
    double EparValTime;

    /// Time from which the model should give a good approximation. Optional field. Default value is 0.1 [day].
    double begOfTimeOfInterest; // local one or taken from e-model
    /// Time (age???) up to which the model should give a good approximation.
    double endOfTimeOfInterest; // local one or taken from e-model
    /// Associated linearElasticMaterial, with E = 1.
    LinearElasticMaterial *linearElasticMaterial;
    /// Partial moduli of individual units.
    FloatArray EparVal;
    //FloatArray relaxationTimes;
    /// Characteristic times of individual units (relaxation or retardation times).
    FloatArray charTimes;
    /// Times at which the errors are evaluated if the least-square method is used.
    FloatArray discreteTimeScale;

    /**
     * Scaling factor transforming the simulation time units into days
     * (gives the number of simulation time units in one day,
     *  e.g. 86400 if the simulation works with seconds as the time units)
     */
    double timeFactor;

public:
    RheoChainMaterial(int n, Domain *d);
    virtual ~RheoChainMaterial();

    virtual void giveCharacteristicMatrix(FloatMatrix &answer,
                                          MatResponseForm form,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);

    virtual void giveRealStressVector(FloatArray &answer,  MatResponseForm form, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
    { answer.resize(0); }

    /// Evaluation of the incremental modulus.
    virtual double giveEModulus(GaussPoint *gp, TimeStep *atTime) { return 0.0; }

    /// Evaluation of the moduli of individual units.
    virtual void computeCharCoefficients(FloatArray &answer, GaussPoint *gp, double atTime) {};

    /// Update of MatStatus to the newly reached (equilibrium) state.
    virtual void updateYourself(GaussPoint *gp, TimeStep *tStep) {};

    // identification and auxiliary functions
    virtual int hasNonLinearBehaviour() { return 0; }
    virtual int hasMaterialModeCapability(MaterialMode mode);
    virtual const char *giveClassName() const { return "RheoChainMaterial"; }
    virtual classType giveClassID() const { return RheoChainMaterialClass; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    // store & restore context functions
    virtual contextIOResultType saveIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp);
    virtual contextIOResultType restoreIPContext(DataStream *stream, ContextMode mode, GaussPoint *gp);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseForm form, MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);

    virtual void computeStressIndependentStrainVector(FloatArray &answer,
                                                      GaussPoint *gp, TimeStep *tStep, ValueModeType mode);

    /**
     * Computes, for the given integration point,
     * the strain vector induced by stress-independent shrinkage.
     * @param answer Returned strain vector.
     * @param form Material response form.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     * @param mode Determines response mode (Total or incremental).
     */
    virtual void giveShrinkageStrainVector(FloatArray &answer,
                                           MatResponseForm form,
                                           GaussPoint *gp,
                                           TimeStep *tStep,
                                           ValueModeType mode)
    { answer.resize(0); }

    // Note: must take LoadResponseMode into account
    /**
     * Computes, for the given integration point,
     * the strain vector induced by the stress history (typically creep strain).
     * @param answer Computed strains.
     * @param form Material response form.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is the current time step).
     * @param mode Determines response mode.
     */
    virtual void giveEigenStrainVector(FloatArray &answer, MatResponseForm form,
                                       GaussPoint *gp, TimeStep *tStep, ValueModeType mode) {};

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
    /**
     * If only incremental shrinkage strain formulation is provided, then total shrinkage strain must be tracked
     * in status in order to be able to compute total value.
     */
    virtual int hasIncrementalShrinkageFormulation() { return 0; }

    /**
     * Generates discrete times starting from time "from" to time "to"
     * uniformly distributed in log time scale.
     * The time interval (to-from) is divided to nsteps intervals.
     * We return times starting from ("from" + first increment)
     * @param from Starting time
     * @param to End time
     * @param nsteps Number of discrete steps.
     */
    void generateLogTimeScale(FloatArray &answer, double from, double to, int nsteps);
    const FloatArray &giveDiscreteTimes();
    /// Evaluation of the creep compliance function.
    virtual double computeCreepFunction(GaussPoint *gp, double ofAge, double atTime) = 0;

    /**
     * Evaluation of the relaxation function at given times.
     * This functions solves numerically an integral equation of the form
     * @f[
     * \varepsilon(t) = \int_{0}^{t} J(t, \tau) \mathrm{d}\sigma(\tau) + \varepsilon_n(t)
     * @f]
     * where @f$ \varepsilon_n(t) @f$ is stress-independent deformation, for the case
     * where @f$ \varepsilon(t) = 1 @f$ is kept at constant value in time.
     *
     * @param t0 Age of material when load is applied.
     * @param tr Age of material when relaxation has begun ???
     * @param atTimes At which times the relaxation function will be evaluated.
     * @warning atTimes should be uniformly distributed in log time scale and relatively dense (100 intervals) in order to achieve a reasonable accuracy.
     */
    void computeDiscreteRelaxationFunction(FloatArray &answer, GaussPoint *gp,
                                           const FloatArray &atTimes,
                                           double t0, double tr);

    /// Evaluation of elastic compliance matrix for unit Young's modulus.
    void giveUnitComplianceMatrix(FloatMatrix &answer, MatResponseForm form,
                                  GaussPoint *gp, TimeStep *tStep);
    /// Evaluation of elastic stiffness matrix for unit Young's modulus.
    void giveUnitStiffnessMatrix(FloatMatrix &answer,
                                 MatResponseForm form, GaussPoint *gp, TimeStep *tStep);

    /// Update of partial moduli of individual chain units
    void updateEparModuli(GaussPoint *gp, double atTime);

    /// Access to partial modulus of a given unit
    double giveEparModulus(int iChain);

    /// Evaluation of characteristic times
    virtual void computeCharTimes();

    /// Access to the characteristic time of a given unit
    double giveCharTime(int);

    /// Exponent to be used with the char time of a given unit, usually = 1.0
    virtual double giveCharTimeExponent(int i) { return 1.0; }

    /// Access to the underlying linear elastic material with unit Young's modulus
    LinearElasticMaterial *giveLinearElasticMaterial();

    /// Access to the time up to which the response should be accurate
    double giveEndOfTimeOfInterest();

    virtual void givePlaneStressStiffMtrx(FloatMatrix &answer,
                                          MatResponseForm, MatResponseMode mmode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);
    virtual void givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                          MatResponseForm, MatResponseMode mmode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);
    virtual void give1dStressStiffMtrx(FloatMatrix &answer,
                                       MatResponseForm, MatResponseMode mmode,
                                       GaussPoint *gp,
                                       TimeStep *tStep);
    virtual void give2dBeamLayerStiffMtrx(FloatMatrix &answer,
                                          MatResponseForm, MatResponseMode mmode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);
    virtual void give2dPlateLayerStiffMtrx(FloatMatrix &answer,
                                           MatResponseForm, MatResponseMode mmode,
                                           GaussPoint *gp,
                                           TimeStep *tStep);
    virtual void give3dShellLayerStiffMtrx(FloatMatrix &answer,
                                           MatResponseForm, MatResponseMode mmode,
                                           GaussPoint *gp,
                                           TimeStep *tStep);
    /**
     * Computes, for the given integration point,
     * the strain vector induced by stress-independent
     * internal processes in the material.
     * Takes into account only temperature and shrinkage-induced strains.
     * @param answer Returned strain vector.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     * @param mode Determines response mode (Total or incremental).
     */
    void computeTrueStressIndependentStrainVector(FloatArray &answer, GaussPoint *gp,
                                                  TimeStep *tStep, ValueModeType mode);
};
} // end namespace oofem
#endif // rheochm_h
