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

#ifndef rheochm_h
#define rheochm_h

/// thermal, shrinkage (drying & autogenous) and creep strains are stored for export
#define keep_track_of_strains

#include "../sm/Materials/structuralmaterial.h"
#include "Materials/linearelasticmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"

#include "matconst.h"
#include "../sm/Elements/structuralelement.h"
#include "../sm/Materials/structuralms.h"

///@name Input fields for RheoChainMaterial
//@{
#define _IFT_RheoChainMaterial_n "n"
#define _IFT_RheoChainMaterial_alphaOne "a1"
#define _IFT_RheoChainMaterial_alphaTwo "a2"
#define _IFT_RheoChainMaterial_lattice "lattice"
#define _IFT_RheoChainMaterial_relmatage "relmatage"
#define _IFT_RheoChainMaterial_begoftimeofinterest "begoftimeofinterest"
#define _IFT_RheoChainMaterial_endoftimeofinterest "endoftimeofinterest"
#define _IFT_RheoChainMaterial_timefactor "timefactor"
//@}

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
    std :: vector< FloatArray >hiddenVars;
    std :: vector< FloatArray >tempHiddenVars;

    /**
     * Total shrinkage strain (needed only when the shrinkage evolution
     * is described in the incremental form).
     */
    FloatArray shrinkageStrain;

#ifdef keep_track_of_strains
    double thermalStrain;
    double tempThermalStrain;
#endif

public:
    RheoChainMaterialStatus(int n, Domain *d, GaussPoint *g, int nunits);
    virtual ~RheoChainMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual const FloatArray &giveViscoelasticStressVector() const { return stressVector; }

    FloatArray &giveHiddenVarsVector(int i) { return hiddenVars [ i - 1 ]; }
    FloatArray &giveTempHiddenVarsVector(int i) { return tempHiddenVars [ i - 1 ]; }
    FloatArray *letHiddenVarsVectorBe(int i, FloatArray *);
    void letTempHiddenVarsVectorBe(int i, FloatArray &valueArray);

    FloatArray *giveShrinkageStrainVector() { return & shrinkageStrain; }
    void setShrinkageStrainVector(FloatArray src) { shrinkageStrain = std :: move(src); }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);


#ifdef keep_track_of_strains
    void setTempThermalStrain(double src) { tempThermalStrain = src; }
    double giveTempThermalStrain(void) { return tempThermalStrain; }
    double giveThermalStrain(void) { return thermalStrain; }
#endif

    // definition
    virtual const char *giveClassName() const { return "RheoChainMaterialStatus"; }
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

    bool lattice;
    /// Poisson's ratio (assumed to be constant, unaffected by creep).
    double nu;
    /// Parameters for the lattice model
    double alphaOne, alphaTwo;
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

    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_2dBeamLayer(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_PlateLayer(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }

    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
    { answer.clear(); }

    /// Evaluation of the incremental modulus.
    virtual double giveEModulus(GaussPoint *gp, TimeStep *tStep) = 0;

    /// Evaluation of the moduli of individual units.
    virtual void computeCharCoefficients(FloatArray &answer, double tStep) = 0;

    // identification and auxiliary functions
    virtual int hasNonLinearBehaviour() { return 0; }
    virtual int hasMaterialModeCapability(MaterialMode mode);
    virtual const char *giveClassName() const { return "RheoChainMaterial"; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    // store & restore context functions
    virtual contextIOResultType saveIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp);
    virtual contextIOResultType restoreIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);

    virtual void givePlaneStressStiffMtrx(FloatMatrix &answer,
                                          MatResponseMode mmode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);
    virtual void givePlaneStrainStiffMtrx(FloatMatrix &answer,
                                          MatResponseMode mmode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);
    virtual void give1dStressStiffMtrx(FloatMatrix &answer,
                                       MatResponseMode mmode,
                                       GaussPoint *gp,
                                       TimeStep *tStep);

    virtual void give2dLatticeStiffMtrx(FloatMatrix &answer,
                                        MatResponseMode mmode,
                                        GaussPoint *gp,
                                        TimeStep *tStep);

    virtual void give3dLatticeStiffMtrx(FloatMatrix &answer,
                                        MatResponseMode mmode,
                                        GaussPoint *gp,
                                        TimeStep *tStep);

    virtual void computeStressIndependentStrainVector(FloatArray &answer,
                                                      GaussPoint *gp, TimeStep *tStep, ValueModeType mode);

    /**
     * Computes, for the given integration point,
     * the strain vector induced by stress-independent shrinkage.
     * @param answer Returned strain vector.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     * @param mode Determines response mode (Total or incremental).
     */
    virtual void giveShrinkageStrainVector(FloatArray &answer,
                                           GaussPoint *gp,
                                           TimeStep *tStep,
                                           ValueModeType mode)
    { answer.clear(); }

    // Note: must take LoadResponseMode into account
    /**
     * Computes, for the given integration point,
     * the strain vector induced by the stress history (typically creep strain).
     * @param answer Computed strains.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is the current time step).
     * @param mode Determines response mode.
     */
    virtual void giveEigenStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) { }

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    double giveAlphaOne() const { return this->alphaOne; }
    double giveAlphaTwo() const { return this->alphaTwo; }

    /// Evaluation of the creep compliance function at time t when loading is acting from time t_prime
    virtual double computeCreepFunction(double t, double t_prime) = 0;

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
     * @param answer Resulting array of discrete times.
     * @param from Starting time
     * @param to End time
     * @param nsteps Number of discrete steps.
     */
    static void generateLogTimeScale(FloatArray &answer, double from, double to, int nsteps);
    const FloatArray &giveDiscreteTimes();

    /**
     * Evaluation of the relaxation function at given times.
     * This functions solves numerically an integral equation of the form
     * @f[
     * \varepsilon(t) = \int_{0}^{t} J(t, \tau) \mathrm{d}\sigma(\tau) + \varepsilon_n(t)
     * @f]
     * where @f$ \varepsilon_n(t) @f$ is stress-independent deformation, for the case
     * where @f$ \varepsilon(t) = 1 @f$ is kept at constant value in time.
     *
     * @param[out] answer Array with evaluated relaxation function.
     * @param t0 Age of material when load is applied.
     * @param tr Age of material when relaxation has begun ???
     * @param tSteps At which times the relaxation function will be evaluated.
     * @warning tSteps should be uniformly distributed in log time scale and relatively dense (100 intervals) in order to achieve a reasonable accuracy.
     */
    void computeDiscreteRelaxationFunction(FloatArray &answer, const FloatArray &tSteps, double t0, double tr);

    /// Evaluation of elastic compliance matrix for unit Young's modulus.
    void giveUnitComplianceMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep);
    /// Evaluation of elastic stiffness matrix for unit Young's modulus.
    void giveUnitStiffnessMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep);

    /// Update of partial moduli of individual chain units
    virtual void updateEparModuli(double tStep);

    /// Access to partial modulus of a given unit
    double giveEparModulus(int iChain);

    /// Evaluation of characteristic times
    virtual void computeCharTimes();

    /// Access to the characteristic time of a given unit
    double giveCharTime(int) const;

    /// Exponent to be used with the char time of a given unit, usually = 1.0
    virtual double giveCharTimeExponent(int i) { return 1.0; }

    /// Access to the underlying linear elastic material with unit Young's modulus
    LinearElasticMaterial *giveLinearElasticMaterial();

    /// Access to the time up to which the response should be accurate
    double giveEndOfTimeOfInterest();

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
