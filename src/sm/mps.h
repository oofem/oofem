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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef mps_h
#define mps_h

#include "kelvinChSolM.h"

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
    /// Values of humidity and temperature in a particular GP and their increment.
    double hum;
    double hum_increment;
    double T;
    double T_increment;
    double T_max;
    /// Hidden variable - equivalent time: necessary to compute solidified volume.
    double equivalentTime;
    double equivalentTimeTemp;
    double flowTermViscosity;
    double flowTermViscosityTemp;

public:
    MPSMaterialStatus(int n, Domain *d, GaussPoint *g, int nunits);
    ~MPSMaterialStatus() {}

    virtual void updateYourself(TimeStep *tStep);

    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    /// Returns the humidity.
    double giveHum() { return hum; }
    void setHum(const double src) { hum = src; }

    double giveHumIncrement() { return hum_increment; }
    void setHumIncrement(const double src) { hum_increment = src; }

    /// Returns the temperature.
    double giveT() { return T; }
    void setT(const double src) { T = src; }

    /// Returns the temperature increment.
    double giveTIncrement() { return T_increment; }
    void setTIncrement(const double src) { T_increment = src; }

    double giveTmax() { return T_max; }
    void setTmax(const double src) { T_max = src; }

    double giveEquivalentTime() { return equivalentTime; }
    void setEquivalentTime(const double src) { equivalentTimeTemp = src; }

    double giveFlowTermViscosity() { return flowTermViscosity; }
    double giveFlowTermViscosityTemp() { return flowTermViscosityTemp; }
    void setFlowTermViscosityTemp(const double src) { flowTermViscosityTemp = src; }

    // definition
    const char *giveClassName() const { return "MPSMaterialStatus"; }
    classType giveClassID() const { return MPSMaterialStatusClass; }
};


/**
 * This class implements the extended B3 model for concrete creep and shrinkage
 * based on the microprestress-solidification theory.
 * The implementation exploits a solidifying Kelvin chain.
 * Creep is affected by variable temperature and humidity.
 */
class MPSMaterial : public KelvinChainSolidMaterial
{
protected:
    double talpha; ///< Thermal dilatation coeff.
    double t0; ///< Age when temperature or humidity starts to change
    double q1, q2, q3, q4; ///< Predicted data.

    enum coupledAnalysisType { Basic, MPS } CoupledAnalysis;

    double EspringVal; ///< Elastic modulus of the aging spring (first member of Kelvin chain if retardation spectrum is used).

    // additional parameters for sorption isotherm (used to compute relative humidity from water content)
    double w_h, n, a; //constant (obtained from experiments) A [Pedersen, 1990]

    // MPS theory parameters
    /// Proportionality parameter between change of humidity and shrinkage.
    double kSh;
    /// Fluidity parameter used in viscosity evolution equation.
    double muS;
    /// kappaT replaces ln(h) on RHS of the differential equation describing evolution of MPS.
    double kappaT;
    /// Parameter reducing creep effects of thermal cycling.
    double cyclicTparam;
    /// Reference room temperature for MPS algorithm [K].
    double roomTemperature;
    /// Activation energies [K].
    double QEtoR, QRtoR, QStoR;
    /// Parameters that control the effect of humidity on rates of hydration, creep and microprestress relaxation.
    double alphaE, alphaR, alphaS;

public:
    MPSMaterial(int n, Domain *d) : KelvinChainSolidMaterial(n, d) { }
    ~MPSMaterial() { }

    const char *giveClassName()  const { return "MPSMaterial"; }
    classType giveClassID() const { return MPSMaterialClass; }

    IRResultType initializeFrom(InputRecord *ir);

    void updateYourself(GaussPoint *gp, TimeStep *tStep);

    void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual void giveShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
                                           GaussPoint *gp, TimeStep *tStep, ValueModeType mode);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
    void predictParametersFrom(double, double, double, double);

    /// Evaluation of characteristic times.
    virtual void computeCharTimes();

    /// Evaluation of characteristic moduli of the non-aging Kelvin chain.
    virtual void computeCharCoefficients(FloatArray &answer, GaussPoint *gp, double);

    /// Evaluation of the incremental modulus.
    virtual double giveEModulus(GaussPoint *gp, TimeStep *tStep);

    /// Evaluation of the relative volume of the solidified material.
    virtual double computeSolidifiedVolume(GaussPoint *gp, TimeStep *tStep);

    /// Factors for exponential algorithm.
    virtual double computeBetaMu(GaussPoint *gp, TimeStep *tStep, double Mu);
    virtual double computeLambdaMu(GaussPoint *gp, TimeStep *tStep, double Mu);

    /// Evaluation of the flow term viscosity.
    double computeFlowTermViscosity(GaussPoint *gp, TimeStep *tStep);

    /// Returns initial value of the flow term viscosity.
    double giveInitViscosity(TimeStep *tStep);

    virtual void  giveEigenStrainVector(FloatArray &answer, MatResponseForm form,
                                        GaussPoint *gp, TimeStep *tStep, ValueModeType mode);

    virtual int hasIncrementalShrinkageFormulation() { return 1; }


    /// Evaluation of the shrinkageStrainVector - shrinkage is fully dependent on humidity rate in given GP.
    void computePointShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
                                           GaussPoint *gp, TimeStep *tStep);


    double inverse_sorption_isotherm(double w);

    /**
     * Gives value of humidity at given GP and time step.
     * option = 0 ... beginning of the time step.
     * option = 1 ... end of the time step.
     * option = 2 ... average values.
     * option = 3 ... incremental values.
     */
    double giveHumidity(GaussPoint *gp, TimeStep *tStep, int option);

    /**
     * Gives value of temperature at given GP and time step.
     * option = 0 ... beginning of the time step.
     * option = 1 ... end of the time step.
     * option = 2 ... average values.
     * option = 3 ... incremental values.
     */
    double giveTemperature(GaussPoint *gp, TimeStep *tStep, int option);

    /**
     * Evaluation of the factor transforming real time to reduced time (effect on the flow term).
     * option = 0 ... beginning of the time step.
     * option = 1 ... end of the time step.
     * option = 2 ... average value.
     */
    double computePsiR(GaussPoint *gp, TimeStep *tStep, double option);

    /**
     * Evaluation of the factor transforming real time to reduced time (effect on the evolution of microprestress).
     */
    double computePsiS(GaussPoint *gp, TimeStep *tStep);

    /**
     * Evaluation of the factor transforming real time to equivalent time (effect on the solidified volume).
     */
    double computePsiE(GaussPoint *gp, TimeStep *tStep);

    /**
     * Computes equivalent time at given time step and GP.
     * If option == 0, equivalentTime is evaluated in the middle of the time step (to determine solidified ratio).
     * If option == 1, equivalentTime is evaluated at the end of the time step. (for updating).
     */
    double computeEquivalentTime(GaussPoint *gp, TimeStep *tStep, int option);
};
} // end namespace oofem
#endif // mps_h
