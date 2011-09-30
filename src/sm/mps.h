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
 *               Copyright (C) 1993 - 2010   Borek Patzak
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
/**********     B3SolidMaterialStatus - HUMIDITY ****************************************/

class MPSMaterialStatus : public KelvinChainSolidMaterialStatus
{
    /**
     * This class implements associated Material Status to MPSMaterial,
     * which corresponds to a model for humidity- and temperature-dependent
     * creep of concrete according to the microprestress-solidification theory.
     * At room temperature and 100% relative humidity, it reduces to model B3
     * for basic creep of concrete.
     */

protected:
    /// values of humidity and temperature in a particular GP and their increment
    double hum;
    double hum_increment;
    double T;
    double T_increment;
    double T_max;
    /// hidden variable - equivalent time: necessary to compute solidified volume
    double equivalentTime;
    double equivalentTimeTemp;
    double flowTermViscosity;
    double flowTermViscosityTemp;

public:
    MPSMaterialStatus(int n, Domain *d, GaussPoint *g, int nunits);
    ~MPSMaterialStatus() {}
    //void printOutputAt(FILE *file, TimeStep *tStep);

    /// update after new equilibrium state reached
    virtual void updateYourself(TimeStep *);

    /// save current context (state) into a stream
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /// restore current context (state) from a stream
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    /// HUMIDITY
    double giveHum(void) { return hum; }
    void setHum(const double src) { hum = src; }

    double giveHumIncrement(void) { return hum_increment; }
    void setHumIncrement(const double src) { hum_increment = src; }

    /// TEMPERATURE
    double giveT(void) { return T; }
    void setT(const double src) { T = src; }

    double giveTIncrement(void) { return T_increment; }
    void setTIncrement(const double src) { T_increment = src; }

    double giveTmax(void) { return T_max; }
    void setTmax(const double src) { T_max = src; }

    /// HIDDEN VARIABLES
    double giveEquivalentTime(void) { return equivalentTime; }
    void setEquivalentTime(const double src) { equivalentTimeTemp = src; }

    double giveFlowTermViscosity(void) { return flowTermViscosity; }
    double giveFlowTermViscosityTemp(void) { return flowTermViscosityTemp; }
    void setFlowTermViscosityTemp(const double src) { flowTermViscosityTemp = src; }

    // definition
    const char *giveClassName() const { return "MPSMaterialStatus"; }
    classType             giveClassID() const
    { return MPSMaterialStatusClass; }
};



//=================================================================================


class MPSMaterial : public KelvinChainSolidMaterial
{
    /**
     * This class implements the extended B3 model for concrete creep and shrinkage
     * based on the microprestress-solidification theory.
     * The implementation exploits a solidifying Kelvin chain.
     * Creep is affected by variable temperature and humidity.
     */
protected:

    double talpha;   // thermal dilatation coeff.
    double t0; // age when temperature or humidity starts to change
    double q1, q2, q3, q4; // predicted data

    enum coupledAnalysisType { Basic, MPS } CoupledAnalysis;

    double EspringVal; // elastic modulus of the aging spring (first member of Kelvin chain if retardation spectrum is used)

    // additional parameters for sorption isotherm (used to compute relative humidity from water content)
    double w_h, n, a; //constant (obtained from experiments) A [Pedersen, 1990]

    // MPS theory parameters
    /// proportionality parameter between change of humidity and shrinkage
    double kSh;
    /// fluidity parameter used in viscosity evolution equation
    double muS;
    /// kappaT replaces ln(h) on RHS of the differential equation describing evolution of MPS
    double kappaT;
    /// parameter reducing creep effects of thermal cycling
    double cyclicTparam;
    /// reference room temperature for MPS algorithm [K]
    double roomTemperature;
    /// activation energies
    double QEtoR, QRtoR, QStoR; //[K]
    /// parameters that control the effect of humidity on rates of hydration, creep and microprestress relaxation
    double alphaE, alphaR, alphaS; //[-]

public:
    MPSMaterial(int n, Domain *d) : KelvinChainSolidMaterial(n, d) { }
    ~MPSMaterial() { }

    const char *giveClassName()  const { return "MPSMaterial"; }
    classType giveClassID()          const { return MPSMaterialClass; }

    IRResultType initializeFrom(InputRecord *ir);

    /// updates MatStatus to the newly reached (equilibrium) state
    void updateYourself(GaussPoint *gp, TimeStep *);

    /**
     * Returns a vector of coefficients of thermal dilatation in direction
     * of each material principal (local) axis.
     * @param answer vector of thermal dilatation coefficients
     * @param gp integration point
     * @param tStep time step (most models are able to respond only when atTime is current time step)
     */
    void giveThermalDilatationVector(FloatArray &answer, GaussPoint *, TimeStep *);

    virtual void giveShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
                                           GaussPoint *gp, TimeStep *atTime, ValueModeType mode);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;



protected:
    void  predictParametersFrom(double, double, double, double);

    /// evaluation of characteristic times
    virtual void    computeCharTimes();

    /// evaluation of characteristic moduli of the non-aging Kelvin chain
    virtual void         computeCharCoefficients(FloatArray &answer, GaussPoint *gp, double);

    /// evaluation of the incremental modulus
    virtual double       giveEModulus(GaussPoint *gp, TimeStep *atTime);

    /// evaluation of the relative volume of the solidified material
    virtual double computeSolidifiedVolume(GaussPoint *gp, TimeStep *atTime);

    /// factors for exponential algorithm
    virtual double computeBetaMu(GaussPoint *gp, TimeStep *atTime, double Mu);
    virtual double computeLambdaMu(GaussPoint *gp, TimeStep *atTime, double Mu);

    /// evaluation of the flow term viscosity
    double computeFlowTermViscosity(GaussPoint *gp, TimeStep *atTime);

    /// returns initial value of the flow term viscosity
    double giveInitViscosity(TimeStep *atTime);

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
                                        GaussPoint *gp, TimeStep *atTime, ValueModeType mode);

    /** if only incremental shrinkage strain formulation is provided, then total shrinkage strain must be tracked
     * in status in order to be able to compute total value. */
    virtual int  hasIncrementalShrinkageFormulation() { return 1; }


    /// evaluation of the shrinkageStrainVector - shrinkage is fully dependent on humidity rate in given GP
    void computePointShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
                                           GaussPoint *gp, TimeStep *atTime);


    double inverse_sorption_isotherm(double w);

    /// gives value of humidity at given GP and timestep
    /// option = 0 ... beginning of the time step
    /// option = 1 ... end of the time step
    /// option = 2 ... average values
    /// option = 3 ... incremental values
    double giveHumidity(GaussPoint *gp, TimeStep *atTime, int option);

    /// gives value of temperature at given GP and timestep
    /// option = 0 ... beginning of the time step
    /// option = 1 ... end of the time step
    /// option = 2 ... average values
    /// option = 3 ... incremental values
    double giveTemperature(GaussPoint *gp, TimeStep *atTime, int option);

    /// evaluation of the factor transforming real time to reduced time (effect on the flow term)
    /// option = 0 ... beginning of the time step
    /// option = 1 ... end of the time step
    /// option = 2 ... average value
    double computePsiR(GaussPoint *gp, TimeStep *atTime, double option);

    /// evaluation of the factor transforming real time to reduced time (effect on the evolution of microprestress)
    double computePsiS(GaussPoint *gp, TimeStep *atTime);

    /// evaluation of the factor transforming real time to equivalent time (effect on the solidified volume)
    double computePsiE(GaussPoint *gp, TimeStep *atTime);

    /// computes equivalent time at given time step and GP.
    /// If option == 0, equivalentTime is evaluated in the middle of the time step (to determine solidified ratio).
    /// If option == 1, equivalentTime is evaluated at the end of the time step. (for updating).
    double computeEquivalentTime(GaussPoint *gp, TimeStep *atTime, int option);
};
} // end namespace oofem
#endif // mps_h
