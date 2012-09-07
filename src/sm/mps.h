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

//   ***********************************************************************************
//   *** CLASS Material Model based no Microprestress Solidification Theory & Status ***
//   ***********************************************************************************

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
    /// Values of humidity and temperature in a particular GP and their increment
    double hum;
    double hum_increment;
    double T;
    double T_increment;
    double T_max;
    /// Hidden variable - equivalent time: necessary to compute solidified volume
    double equivalentTime;
    double equivalentTimeTemp;
    double flowTermViscosity;
    double flowTermViscosityTemp;

public:
    MPSMaterialStatus(int n, Domain *d, GaussPoint *g, int nunits);
    virtual ~MPSMaterialStatus() {}

    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    /// Returns relative humidity
    double giveHum(void) { return hum; }
    /// Stores relative humidity
    void setHum(const double src) { hum = src; }
    /// Returns relative humidity increment
    double giveHumIncrement(void) { return hum_increment; }
    /// Stores relative humidity increment
    void setHumIncrement(const double src) { hum_increment = src; }

    /// Returns temperature
    double giveT(void) { return T; }
    /// Stores temperature
    void setT(const double src) { T = src; }

    /// Returns temperature increment
    double giveTIncrement(void) { return T_increment; }
    /// Stores temperature increment
    void setTIncrement(const double src) { T_increment = src; }

    /// Returns previously maximum reached temperature
    double giveTmax(void) { return T_max; }
    /// Stores maximum reached temperature
    void setTmax(const double src) { T_max = src; }

    /// Returns equivalent time
    double giveEquivalentTime(void) { return equivalentTime; }
    /// Stores equivalent time
    void setEquivalentTime(const double src) { equivalentTimeTemp = src; }

    /// Returns viscosity of the flow term (associated with q4 and microprestress evolution)
    double giveFlowTermViscosity(void) { return flowTermViscosity; }
    double giveFlowTermViscosityTemp(void) { return flowTermViscosityTemp; }
    void setFlowTermViscosityTemp(const double src) { flowTermViscosityTemp = src; }

    // definition
    virtual const char *giveClassName() const { return "MPSMaterialStatus"; }
    virtual classType giveClassID() const { return MPSMaterialStatusClass; }
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

    /// thermal dilatation coeff.
    double talpha;
    /// age when temperature or humidity starts to change
    double t0;
    /// compliances of the B3 model
    double q1, q2, q3, q4;
    /// constant equal to one day in time units of analysis (eg. 86400 if the analysis runs in seconds)
    double lambda0;


    enum coupledAnalysisType { Basic, MPS } CoupledAnalysis;

    double EspringVal; // elastic modulus of the aging spring (first member of Kelvin chain if retardation spectrum is used)

    /// additional parameters for sorption isotherm (used to compute relative humidity from water content)
    double w_h, n, a; //constant (obtained from experiments) A [Pedersen, 1990]

    // MPS theory parameters
    /// proportionality parameter between change of humidity and shrinkage
    double kSh;
    /// fluidity parameter used in viscosity evolution equation
    double muS;
    /// kappaT replaces ln(h) on RHS of the differential equation describing evolution of MPS
    double kappaT;
    /// parameter reducing creep effects of thermal cycling
    double ct;
    /// reference room temperature for MPS algorithm [K]
    double roomTemperature;
    /// activation energies
    double QEtoR, QRtoR, QStoR; //[K]
    /// parameters that control the effect of humidity on rates of hydration, creep and microprestress relaxation
    double alphaE, alphaR, alphaS; //[-]

public:
    MPSMaterial(int n, Domain *d) : KelvinChainSolidMaterial(n, d) { }
    virtual ~MPSMaterial() { }

    virtual const char *giveClassName()  const { return "MPSMaterial"; }
    virtual classType giveClassID()          const { return MPSMaterialClass; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void updateYourself(GaussPoint *gp, TimeStep *tStep);

    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual void giveShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
                                           GaussPoint *gp, TimeStep *atTime, ValueModeType mode);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
    void predictParametersFrom(double, double, double, double, double);

    virtual void computeCharTimes();

    /// Evaluation of characteristic moduli of the non-aging Kelvin chain
    virtual void computeCharCoefficients(FloatArray &answer, GaussPoint *gp, double);

    virtual double giveEModulus(GaussPoint *gp, TimeStep *atTime);

    virtual double computeSolidifiedVolume(GaussPoint *gp, TimeStep *atTime);

    virtual double computeBetaMu(GaussPoint *gp, TimeStep *atTime, int Mu);
    virtual double computeLambdaMu(GaussPoint *gp, TimeStep *atTime, int Mu);

    /// Evaluation of the flow term viscosity
    double computeFlowTermViscosity(GaussPoint *gp, TimeStep *atTime);

    /// Returns initial value of the flow term viscosity
    double giveInitViscosity(TimeStep *atTime);

    virtual void  giveEigenStrainVector(FloatArray &answer, MatResponseForm form,
                                        GaussPoint *gp, TimeStep *atTime, ValueModeType mode);

    virtual int hasIncrementalShrinkageFormulation() { return 1; }


    /// Evaluation of the shrinkageStrainVector - shrinkage is fully dependent on humidity rate in given GP
    void computePointShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
                                           GaussPoint *gp, TimeStep *atTime);


    double inverse_sorption_isotherm(double w);

    /// Gives value of humidity at given GP and timestep
    /// option = 0 ... beginning of the time step
    /// option = 1 ... end of the time step
    /// option = 2 ... average values
    /// option = 3 ... incremental values
    double giveHumidity(GaussPoint *gp, TimeStep *atTime, int option);

    /// Gives value of temperature at given GP and timestep
    /// option = 0 ... beginning of the time step
    /// option = 1 ... end of the time step
    /// option = 2 ... average values
    /// option = 3 ... incremental values
    double giveTemperature(GaussPoint *gp, TimeStep *atTime, int option);

    /// Evaluation of the factor transforming real time to reduced time (effect on the flow term)
    /// option = 0 ... beginning of the time step
    /// option = 1 ... end of the time step
    /// option = 2 ... average value
    double computePsiR(GaussPoint *gp, TimeStep *atTime, int option);

    /// Evaluation of the factor transforming real time to reduced time (effect on the evolution of microprestress)
    double computePsiS(GaussPoint *gp, TimeStep *atTime);

    /// Evaluation of the factor transforming real time to equivalent time (effect on the solidified volume)
    double computePsiE(GaussPoint *gp, TimeStep *atTime);

    /// Computes equivalent time at given time step and GP.
    /// If option == 0, equivalentTime is evaluated in the middle of the time step (to determine solidified ratio).
    /// If option == 1, equivalentTime is evaluated at the end of the time step. (for updating).
    double computeEquivalentTime(GaussPoint *gp, TimeStep *atTime, int option);

    friend class RankineMPSmat;
};
} // end namespace oofem
#endif // mps_h
