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
 *               Copyright (C) 1993 - 2015   Borek Patzak
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

//   ***********************************************************************************
//   *** CLASS Material Model based on Microprestress Solidification Theory & Status ***
//   ***********************************************************************************

#ifndef mps_h
#define mps_h

#include "kelvinChSolM.h"

///@name Input fields for MPSMaterial
//@{
#define _IFT_MPSMaterial_Name "mps"
#define _IFT_MPSMaterial_mode "mode"
#define _IFT_MPSMaterial_coupledanalysistype "coupledanalysistype"
#define _IFT_MPSMaterial_fc "fc"
#define _IFT_MPSMaterial_cc "cc"
#define _IFT_MPSMaterial_wc "w/c"
#define _IFT_MPSMaterial_ac "a/c"
#define _IFT_MPSMaterial_q1 "q1"
#define _IFT_MPSMaterial_q2 "q2"
#define _IFT_MPSMaterial_q3 "q3"
#define _IFT_MPSMaterial_q4 "q4"
#define _IFT_MPSMaterial_lambda0 "lambda0"
// in input record t0 replaced by relMatAge
//#define _IFT_MPSMaterial_t0 "t0"
#define _IFT_MPSMaterial_ksh "ksh"
#define _IFT_MPSMaterial_qetor "qetor"
#define _IFT_MPSMaterial_qrtor "qrtor"
#define _IFT_MPSMaterial_qstor "qstor"
#define _IFT_MPSMaterial_alphae "alphae"
#define _IFT_MPSMaterial_alphar "alphar"
#define _IFT_MPSMaterial_alphas "alphas"
#define _IFT_MPSMaterial_mus "mus"
#define _IFT_MPSMaterial_ktm "ktm"
#define _IFT_MPSMaterial_ktc "ktc"
#define _IFT_MPSMaterial_stiffnessfactor "stiffnessfactor"
#define _IFT_MPSMaterial_p "p"
#define _IFT_MPSMaterial_p_tilde "p_tilde"
#define _IFT_MPSMaterial_k3 "k3"
#define _IFT_MPSMaterial_sh_a "sh_a"
#define _IFT_MPSMaterial_sh_hC "sh_hc"
#define _IFT_MPSMaterial_sh_n "sh_n"
#define _IFT_MPSMaterial_factor_ksh_h "factor_ksh_h"
#define _IFT_MPSMaterial_factor_ksh_fh "factor_ksh_f(h)"
#define _IFT_MPSMaterial_timedependent_ksh "timedependent_ksh"
#define _IFT_MPSMaterial_alpha_as "alpha_as"
#define _IFT_MPSMaterial_eps_cas0 "eps_cas0"
#define _IFT_MPSMaterial_B4_eps_au_infty "b4_eps_au_infty"
#define _IFT_MPSMaterial_B4_tau_au "b4_tau_au"
#define _IFT_MPSMaterial_B4_alpha "b4_alpha"
#define _IFT_MPSMaterial_B4_r_t "b4_r_t"
#define _IFT_MPSMaterial_B4_cem_type "b4_cem_type"
#define _IFT_MPSMaterial_temperInCelsius "temperincelsius"
#define _IFT_MPSMaterial_hydrationTimescaleTF "hydrationtimescaletf"
#define _IFT_MPSMaterial_autoShrinkageTF "autoshrinkagetf"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to MPSMaterial,
 * which corresponds to a model for humidity- and temperature-dependent
 * creep of concrete according to the microprestress-solidification theory.
 * At room temperature and 100% relative humidity, it reduces to model B3
 * for basic creep of concrete.
 */
class MPSMaterialStatus : public KelvinChainSolidMaterialStatus
{
protected:
    /// Values of humidity and temperature in a particular GP and their increment
    double hum = -1.;
    double hum_increment = -1.;
    double T = -1.;
    double T_increment = -1.;
    double T_max = 0.;
    /// Hidden variable - equivalent time: necessary to compute solidified volume
    double equivalentTime = 0.;
    double equivalentTimeTemp = 0.;
    double flowTermViscosity = 0.;
    double flowTermViscosityTemp = -1.;
    /// flag for Emodulus - true if modulus has been already computed in the current time step
    bool storedEmodulusFlag = false;
    double storedEmodulus = -1.;
#ifdef keep_track_of_strains
    double dryingShrinkageStrain = 0.;
    double tempDryingShrinkageStrain = 0.;
    double autogenousShrinkageStrain = 0.;
    double tempAutogenousShrinkageStrain = 0.;
    FloatArray creepStrain;
    FloatArray creepStrainIncrement;
#endif


public:
    MPSMaterialStatus(GaussPoint *g, int nunits);

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    /// Returns relative humidity
    double giveHum() { return hum; }
    /// Stores relative humidity
    void setHum(double src) { hum = src; }
    /// Returns relative humidity increment
    double giveHumIncrement() { return hum_increment; }
    /// Stores relative humidity increment
    void setHumIncrement(double src) { hum_increment = src; }

    /// Returns temperature
    double giveT() { return T; }
    /// Stores temperature
    void setT(double src) { T = src; }

    /// Returns temperature increment
    double giveTIncrement() { return T_increment; }
    /// Stores temperature increment
    void setTIncrement(double src) { T_increment = src; }

    /// Returns previously maximum reached temperature
    double giveTmax() { return T_max; }
    /// Stores maximum reached temperature
    void setTmax(double src) { T_max = src; }

    /// Returns equivalent time
    double giveEquivalentTime() { return equivalentTime; }
    /// Stores equivalent time
    void setEquivalentTimeTemp(double src) { equivalentTimeTemp = src; }

    /// Returns viscosity of the flow term (associated with q4 and microprestress evolution)
    double giveFlowTermViscosity() { return flowTermViscosity; }
    double giveFlowTermViscosityTemp() { return flowTermViscosityTemp; }
    void setFlowTermViscosityTemp(double src) { flowTermViscosityTemp = src; }

    /// Returns Emodulus if computed previously in the same tStep
    void storeEmodulus(double src) { storedEmodulus = src; }
    void setEmodulusFlag(bool src) { storedEmodulusFlag = src; }
    /// Returns Emodulus if computed previously in the same tStep
    double giveStoredEmodulus(void) { return storedEmodulus; }
    bool giveStoredEmodulusFlag(void) { return storedEmodulusFlag; }

#ifdef keep_track_of_strains
    void setTempDryingShrinkageStrain(double src) { tempDryingShrinkageStrain = src; }
    double giveTempDryingShrinkageStrain(void) { return tempDryingShrinkageStrain; }
    double giveDryingShrinkageStrain(void) { return dryingShrinkageStrain; }

    void setTempAutogenousShrinkageStrain(double src) { tempAutogenousShrinkageStrain = src; }
    double giveTempAutogenousShrinkageStrain(void) { return tempAutogenousShrinkageStrain; }
    double giveAutogenousShrinkageStrain(void) { return autogenousShrinkageStrain; }

    void setCreepStrainIncrement(FloatArray src) { creepStrainIncrement  = std :: move(src); }
    const FloatArray &giveCreepStrain() const { return creepStrain; }
#endif

    // definition
    const char *giveClassName() const override { return "MPSMaterialStatus"; }
};


/**
 * This class implements the extended B3 model for concrete creep and shrinkage
 * based on the microprestress-solidification theory.
 * The implementation exploits a solidifying Kelvin chain.
 * Creep is affected by variable temperature and humidity.
 * Comparing to other material models for concrete creep MPS material model is
 * unit independent (except parameters corresponding to concrete strength
 * and composition: fc, c, wc, ac ).
 * In the input record is is necessary to specify time factor equal to 1!!!
 */
class MPSMaterial : public KelvinChainSolidMaterial
{
protected:
    /// age when temperature or humidity starts to change
    double t0 = 0.;
    /// compliances of the B3 model
    double q1 = 0., q2 = 0., q3 = 0., q4 = 0.;
    /// constant equal to one day in time units of analysis (eg. 86400 if the analysis runs in seconds)
    double lambda0 = 0.;


    enum coupledAnalysisType { Basic, MPS_full, MPS_humidity, MPS_temperature } CoupledAnalysis;

    mutable double EspringVal = 0.; // elastic modulus of the aging spring (first member of Kelvin chain if retardation spectrum is used)

    /// additional parameters for sorption isotherm (used to compute relative humidity from water content)
    //double w_h, n, a; //constant (obtained from experiments) A [Pedersen, 1990]

    // MPS theory parameters
    /// proportionality parameter between change of humidity and shrinkage
    double kSh = 0.;
    /// fluidity parameter used in viscosity evolution equation
    double muS = 0., k3 = 0.;
    /// kTm replaces ln(h) on RHS of the differential equation describing evolution of MPS
    double kTm = 0.;
    /// parameter reducing creep effects of thermal cycling (replaces kTm in such case)
    double kTc = 0.;
    /// parameter reducing creep effects of thermal cycling
    double ct = 0.;
    /// reference room temperature for MPS algorithm [K]
    double roomTemperature = 0.;
    /// activation energies
    double QEtoR = 0., QRtoR = 0., QStoR = 0.; //[K]
    /// parameters that control the effect of humidity on rates of hydration, creep and microprestress relaxation
    double alphaE = 0., alphaR = 0., alphaS = 0.; //[-]
    /// exponent in the microprestress/viscosity governing equation
    double p = 0.;
    /// parameters for nonlinear shrinkage function
    double sh_a = 0., sh_hC = 0., sh_n = 0.;
    /// flag activating equivalent-time-dependence of shrinkage coefficient
    bool timeDependent_ksh = false;
    /// nonlinear shrinkage function - user-defined piecewise-linear function
    FloatArray ksh_h, ksh_fh;
    /// parameter for autogenous shrinkage according to fib MC 2010
    double eps_cas0 = 0.;
    /// parameters for autogenous shrinkage according to B4 model
    double b4_eps_au_infty = 0., b4_tau_au = 0., b4_alpha = 0., b4_r_t = 0.;

    /// scaling factor 1. for Pa, 1.e6 for MPa - only for empirical formulas - q1-q4 and ft and gf
    double stiffnessFactor = 0.;

    /// 0 for Kelvin, 273.15 for Celsius
    double temperScaleDifference = 0.;

    /**
     * Further scaling of creep, shrinkage, tensile strength etc. with regards to this time function.
     * Let us assume two hydrating cements with different kinetics. If the fast cement is a reference one,
     * the simulation of the slow one is carried out on the reference cement and scaling by 
     * this time function, which is then ≤1. Results from isothermal calorimetry are typically used for
     * defining this function.
     */
    int hydrationTimescaleTF;

    /**
     * Possibility to prescribe the evolution of autogenous shrinkage strain by an auxiliary time function, the function is evaluated according to the current value of the equivalenet age, i.e. maturity.
     */
    int autoShrinkageTF;



public:
    MPSMaterial(int n, Domain *d) : KelvinChainSolidMaterial(n, d) { }

    const char *giveInputRecordName() const override { return _IFT_MPSMaterial_Name; }
    const char *giveClassName() const override { return "MPSMaterial"; }

    void initializeFrom(InputRecord &ir) override;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    void giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) override;
    //void updateYourself(GaussPoint *gp, TimeStep *tStep) override;

    void giveShrinkageStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const override;

    /// Evaluation of the basic creep compliance function - can be used to compute elastic modulus in derived damage material
    double computeCreepFunction(double t, double t_prime, GaussPoint *gp, TimeStep *tStep) const override;

    double giveEquivalentTime(GaussPoint *gp, TimeStep *tStep) const override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

protected:
    void predictParametersFrom(double, double, double, double);

    void computeCharTimes() override;

    /// Evaluation of characteristic moduli of the non-aging Kelvin chain
    FloatArray computeCharCoefficients(double tPrime, GaussPoint *gp, TimeStep *tStep) const override;

    double giveEModulus(GaussPoint *gp, TimeStep *tStep) const override;

    double computeSolidifiedVolume(GaussPoint *gp, TimeStep *tStep) const override;

    double computeBetaMu(GaussPoint *gp, TimeStep *tStep, int Mu) const override;
    double computeLambdaMu(GaussPoint *gp, TimeStep *tStep, int Mu) const override;

    /// Evaluation of the flow term viscosity
    double computeFlowTermViscosity(GaussPoint *gp, TimeStep *tStep) const;

    /// Returns initial value of the flow term viscosity
    double giveInitViscosity(TimeStep *tStep) const;

    void giveEigenStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const override;

    bool hasIncrementalShrinkageFormulation() const override { return true; }

    /// Evaluation of the shrinkageStrainVector - shrinkage is fully dependent on humidity rate in given GP
    void computePointShrinkageStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) const;

    /// Evaluation of the autogenousShrinkageStrainVector according to fib MC 2010 - autogenous shrinkage is fully dependent on the equivalent age at given GP
    void computeFibAutogenousShrinkageStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) const;

    /// Evaluation of the autogenousShrinkageStrainVector according to Bazant's B4 model. In the model the evolution depends on temperature adjusted age, here on equivalent age (additional humidity influence)
    void computeB4AutogenousShrinkageStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) const;

   /// Evaluation of the autogenousShrinkageStrainVector given by an auxiliary time function (autoShrinkageTF parameter). The time scale in that time function corresponds to the equivalent material age.
   void computeAutogenousShrinkageDefinedByTF(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) const;

    //double inverse_sorption_isotherm(double w) const;

    /// Gives value of humidity at given GP and timestep
    /// option = 0 ... beginning of the time step
    /// option = 1 ... end of the time step
    /// option = 2 ... average values
    /// option = 3 ... incremental values
    double giveHumidity(GaussPoint *gp, TimeStep *tStep, int option) const;

    /// Gives value of temperature at given GP and timestep
    /// option = 0 ... beginning of the time step
    /// option = 1 ... end of the time step
    /// option = 2 ... average values
    /// option = 3 ... incremental values
    double giveTemperature(GaussPoint *gp, TimeStep *tStep, int option) const;

    /// Evaluation of the factor transforming real time to reduced time (effect on the flow term)
    /// option = 0 ... beginning of the time step
    /// option = 1 ... end of the time step
    /// option = 2 ... average value
    double computePsiR(GaussPoint *gp, TimeStep *tStep, int option) const;

    /// Evaluation of the factor transforming real time to reduced time (effect on the evolution of microprestress)
    double computePsiS(GaussPoint *gp, TimeStep *tStep) const;

    /// Evaluation of the factor transforming real time to equivalent time (effect on the solidified volume)
    double computePsiE(GaussPoint *gp, TimeStep *tStep) const;

    /// Computes equivalent time at given time step and GP.
    /// If option == 0, equivalentTime is evaluated in the middle of the time step (to determine solidified ratio).
    /// If option == 1, equivalentTime is evaluated at the end of the time step. (for updating).
    double computeEquivalentTime(GaussPoint *gp, TimeStep *tStep, int option) const;

    friend class RankineMPSmat;
};
} // end namespace oofem
#endif // mps_h
