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

#include "sm/Materials/structuralmaterial.h"
//#include "sm/Materials/linearelasticmaterial.h"
#include "floatarray.h"
#include "floatmatrix.h"

#include "matconst.h"
#include "sm/Elements/structuralelement.h"
#include "sm/Materials/structuralms.h"

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
#define _IFT_RheoChainMaterial_talpha "talpha"
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
    int nUnits = 0;
    /// Hidden (internal) variables, the meaning of which depends on the type of chain.
    std :: vector< FloatArray >hiddenVars;
    std :: vector< FloatArray >tempHiddenVars;

    /**
     * Total shrinkage strain (needed only when the shrinkage evolution
     * is described in the incremental form).
     */
    FloatArray shrinkageStrain;

    /** Internal variable with the meaning of the solution time
     * essential when giveRealStressVector of the viscoelastic material is called
     * more than once in one time step and at one time step.
     */
    double currentTime = 0.;

#ifdef keep_track_of_strains
    double thermalStrain = 0.;
    double tempThermalStrain = 0.;
#endif

public:
    RheoChainMaterialStatus(GaussPoint *g, int nunits);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    virtual const FloatArray &giveViscoelasticStressVector() const { return stressVector; }

    FloatArray &giveHiddenVarsVector(int i) { return hiddenVars [ i - 1 ]; }
    FloatArray &giveTempHiddenVarsVector(int i) { return tempHiddenVars [ i - 1 ]; }
    FloatArray *letHiddenVarsVectorBe(int i, FloatArray *);
    void letTempHiddenVarsVectorBe(int i, FloatArray &valueArray);

    FloatArray *giveShrinkageStrainVector() { return & shrinkageStrain; }
    void setShrinkageStrainVector(FloatArray src) { shrinkageStrain = std :: move(src); }

    /// Returns current time - see explanation near initTempStatus in giveRealStressVector
    double giveCurrentTime() { return currentTime; }
    /// Stores current time
    void setCurrentTime(double src) { currentTime = src; }

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;


#ifdef keep_track_of_strains
    void setTempThermalStrain(double src) { tempThermalStrain = src; }
    double giveTempThermalStrain(void) { return tempThermalStrain; }
    double giveThermalStrain(void) { return thermalStrain; }
#endif

    // definition
    const char *giveClassName() const override { return "RheoChainMaterialStatus"; }
};


/**
 * This class implements a rheologic chain model
 * describing a viscoelastic material.
 * It serves as the parent class for Maxwell and Kelvin chains.
 */
class RheoChainMaterial : public StructuralMaterial
{
protected:
    /// thermal dilatation coeff.
    double talpha = 0.;
    /// Number of (Maxwell or Kelvin) units in the rheologic chain.
    int nUnits = 0;
    /// Physical age of the material at castingTime
    double relMatAge = 0.;

    bool lattice = false;
    /// Poisson's ratio (assumed to be constant, unaffected by creep).
    double nu = 0.;
    /// Parameters for the lattice model
    double alphaOne = 0., alphaTwo = 0.;
    /// Time for which the partial moduli of individual units have been evaluated.
    mutable double EparValTime = -1.;

    /// Time from which the model should give a good approximation. Optional field. Default value is 0.1 [day].
    double begOfTimeOfInterest = 0.; // local one or taken from e-model
    /// Time (age???) up to which the model should give a good approximation.
    double endOfTimeOfInterest = 0.; // local one or taken from e-model
    /// Associated linearElasticMaterial, with E = 1.
    //    LinearElasticMaterial *linearElasticMaterial = nullptr;
    StructuralMaterial *linearElasticMaterial = nullptr;

    /// Partial moduli of individual units.
    mutable FloatArray EparVal;
    //FloatArray relaxationTimes;
    /// Characteristic times of individual units (relaxation or retardation times).
    mutable FloatArray charTimes;
    /// Times at which the errors are evaluated if the least-square method is used.
    FloatArray discreteTimeScale;

    /**
     * Scaling factor transforming the simulation time units into days
     * (gives the number of simulation time units in one day,
     *  e.g. 86400 if the simulation works with seconds as the time units)
     */
    double timeFactor = 0.;

public:
    RheoChainMaterial(int n, Domain *d);
    virtual ~RheoChainMaterial();

    void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                              const FloatArray &reducedStrain, TimeStep *tStep) override;

    FloatArrayF< 6 >giveRealStressVector_3d(const FloatArrayF< 6 > &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast< RheoChainMaterial * >( this )->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF< 4 >giveRealStressVector_PlaneStrain(const FloatArrayF< 4 > &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast< RheoChainMaterial * >( this )->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF< 3 >giveRealStressVector_PlaneStress(const FloatArrayF< 3 > &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast< RheoChainMaterial * >( this )->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF< 1 >giveRealStressVector_1d(const FloatArrayF< 1 > &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast< RheoChainMaterial * >( this )->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF< 2 >giveRealStressVector_2dBeamLayer(const FloatArrayF< 2 > &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast< RheoChainMaterial * >( this )->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }
    FloatArrayF< 5 >giveRealStressVector_PlateLayer(const FloatArrayF< 5 > &strain, GaussPoint *gp, TimeStep *tStep) const override
    {
        FloatArray answer;
        const_cast< RheoChainMaterial * >( this )->giveRealStressVector(answer, gp, strain, tStep);
        return answer;
    }

    FloatArrayF< 6 >giveThermalDilatationVector(GaussPoint *gp, TimeStep *tStep) const override;

    /*    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
     *    { answer.clear(); }*/

    /// Evaluation of the incremental modulus.
    virtual double giveEModulus(GaussPoint *gp, TimeStep *tStep) const = 0;

    /*    virtual double giveIncrementalModulus(GaussPoint *gp, TimeStep *tStep) {
     * if ( (tStep->giveIntrinsicTime() < this->castingTime) && ( this->zeroStiffness > 0. ) ) {
     *  return this->zeroStiffness;
     * } else {
     *  return this->giveEModulus(gp, tStep);
     * }
     * }*/

    /// Evaluation of the moduli of individual units.
    virtual FloatArray computeCharCoefficients(double tPrime, GaussPoint *gp, TimeStep *tStep) const = 0;

    // identification and auxiliary functions
    bool hasMaterialModeCapability(MaterialMode mode) const override;

    bool hasCastingTimeSupport() const override { return true; }

    const char *giveClassName() const override { return "RheoChainMaterial"; }
    void initializeFrom(InputRecord &ir) override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    // store & restore context functions
    void saveIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp) override;
    void restoreIPContext(DataStream &stream, ContextMode mode, GaussPoint *gp) override;

    /**
     * Computes the stiffness matrix for giveRealStressVector of receiver in given integration point, respecting its history.
     * The algorithm should use temporary or equilibrium  history variables stored in integration point status
     * to compute and return required result.
     * @param answer Contains result.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    /*    void giveStiffnessMatrix(FloatMatrix &answer,
     *                               MatResponseMode mode,
     *                               GaussPoint *gp,
     *                               TimeStep *tStep) override;
     */

    FloatMatrixF< 6, 6 >give3dMaterialStiffnessMatrix(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF< 3, 3 >givePlaneStressStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF< 4, 4 >givePlaneStrainStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;
    FloatMatrixF< 1, 1 >give1dStressStiffMtrx(MatResponseMode mmode, GaussPoint *gp, TimeStep *tStep) const override;

    // maybe not needed afterall?
    //	FloatMatrixF< 3, 3 >give2dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const;

    //	FloatMatrixF< 6, 6 >give3dLatticeStiffnessMatrix(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const;

    FloatArray computeStressIndependentStrainVector(GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const override;

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
                                           ValueModeType mode) const
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
    virtual void giveEigenStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const { }

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    double giveAlphaOne() const { return this->alphaOne; }
    double giveAlphaTwo() const { return this->alphaTwo; }

    /// Returns Poisson's ratio.
    double givePoissonsRatio() const { return nu; }

    /// Evaluation of the creep compliance function at time t when loading is acting from time t_prime
    virtual double computeCreepFunction(double t, double t_prime, GaussPoint *gp, TimeStep *tStep) const = 0;

    /// Extended meaning: returns true if the material is cast (target time > casting time) or the precasing time mat is defined
    bool isActivated(TimeStep *tStep) const override {
        if ( this->preCastingTimeMat > 0 ) {
            return true;
        } else {
            return Material :: isActivated(tStep);
        }
    }

    /// By default returns equivalent time in the middle of the time step
    virtual double giveEquivalentTime(GaussPoint *gp, TimeStep *tStep) const
    { return ( tStep->giveTargetTime() - tStep->giveTimeIncrement() / 2 ); }


protected:
    /**
     * If only incremental shrinkage strain formulation is provided, then total shrinkage strain must be tracked
     * in status in order to be able to compute total value.
     */
    virtual bool hasIncrementalShrinkageFormulation() const { return false; }

    /**
     * Generates discrete times starting from time "from" to time "to"
     * uniformly distributed in log time scale.
     * The time interval (to-from) is divided to nsteps intervals.
     * We return times starting from ("from" + first increment)
     * @param from Starting time
     * @param to End time
     * @param nsteps Number of discrete steps.
     * @return Resulting array of discrete times.
     */
    static FloatArray generateLogTimeScale(double from, double to, int nsteps);
    const FloatArray &giveDiscreteTimes() const;

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
    void computeDiscreteRelaxationFunction(FloatArray &answer, const FloatArray &tSteps, double t0, double tr, GaussPoint *gp, TimeStep *tStep) const;

    /// Evaluation of elastic compliance matrix for unit Young's modulus.
    void giveUnitComplianceMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep) const;
    /// Evaluation of elastic stiffness matrix for unit Young's modulus.
    void giveUnitStiffnessMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep) const;

    /// Update of partial moduli of individual chain units
    virtual void updateEparModuli(double tPrime, GaussPoint *gp, TimeStep *tStep) const;

    /// Access to partial modulus of a given unit
    double giveEparModulus(int iChain) const;

    /// Evaluation of characteristic times
    virtual void computeCharTimes();

    /// Access to the characteristic time of a given unit
    double giveCharTime(int) const;

    /// Exponent to be used with the char time of a given unit, usually = 1.0
    virtual double giveCharTimeExponent(int i) const { return 1.0; }

    /// Access to the underlying linear elastic material with unit Young's modulus
    //    LinearElasticMaterial *giveLinearElasticMaterial();
    StructuralMaterial *giveLinearElasticMaterial();

    /// Access to the time up to which the response should be accurate
    double giveEndOfTimeOfInterest();

    /**
     * Computes, for the given integration point,
     * internal processes in the material.
     * Takes into account only temperature and shrinkage-induced strains.
     * @param answer Returned strain vector.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     * @param mode Determines response mode (Total or incremental).
     */
    void computeTrueStressIndependentStrainVector(FloatArray &answer, GaussPoint *gp,
                                                  TimeStep *tStep, ValueModeType mode) const;
};
} // end namespace oofem
#endif // rheochm_h
