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

#ifndef b3solidmat_h
#define b3solidmat_h

#include "kelvinChM.h"

///@name Input fields for B3Material
//@{
#define _IFT_B3SolidMaterial_Name "b3solidmat"
#define _IFT_B3SolidMaterial_emodulimode "emodulimode"
#define _IFT_B3SolidMaterial_microprestress "microprestress"
#define _IFT_B3SolidMaterial_c0 "c0"
#define _IFT_B3SolidMaterial_c1 "c1"
#define _IFT_B3SolidMaterial_ksh "ksh"
#define _IFT_B3SolidMaterial_finalhumidity "finalhumidity"
#define _IFT_B3SolidMaterial_initialhumidity "initialhumidity"
#define _IFT_B3SolidMaterial_lambda0 "lambda0"
//@}

namespace oofem {
/**
 * This class implements associated Material Status to B3SolidMaterial.
 */
class B3SolidMaterialStatus : public KelvinChainMaterialStatus
{
protected:
    /// Microprestresses
    double microprestress_old;
    double microprestress_new;

public:
    B3SolidMaterialStatus(int n, Domain *d, GaussPoint *g, int nunits);
    virtual ~B3SolidMaterialStatus() { }

    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    double giveMPS() const { return microprestress_old; }
    void setMPS(double src) { microprestress_new = src; }

    // definition
    virtual const char *giveClassName() const { return "B3SolidMaterialStatus"; }
};


/**
 * This class implements the B3 model for concrete creep and shrinkage
 * based on the solidification theory.
 * The implementation exploits a solidifying Kelvin chain.
 */
class B3SolidMaterial : public KelvinChainMaterial
{
protected:
    double t0;
    double w, E28, q1, q2, q3, q4, q5; // predicted data

    /// constant equal to one day in time units of analysis (eg. 86400 if the analysis runs in seconds)
    double lambda0;

    enum b3ShModeType { B3_NoShrinkage, B3_AverageShrinkage, B3_PointShrinkage } shMode;
    /// Additional parameters for average cross section shrinkage
    double EpsSinf, kt, ks, vs, hum;
    /// Additional parameters for free shrinkage at material point
    double es0, r, rprime, at;
    // Additional parameters for sorption isotherm (used to compute relative humidity from water content)
    double w_h;    ///< Constant water content (obtained from experiments) w_h [Pedersen, 1990]
    double n;      ///< Constant-exponent (obtained from experiments) n [Pedersen, 1990]
    double a;      ///< Constant (obtained from experiments) A [Pedersen, 1990]
    double talpha; ///< Thermal dilatation coeff.
    double EspringVal; ///< elastic modulus of the aging spring (first member of Kelvin chain if retardation spectrum is used)
    /**
     * If 0, analysis of retardation spectrum is used for evaluation of Kelvin units moduli (default).
     * If 1, least-squares method is used for evaluation of Kelvin units moduli.
     */
    int EmoduliMode;
    /**
     * If 1, computation exploiting Microprestress solidification theory is done.
     * Default value is 0 = without external fields it can be used for basic creep.
     */
    int MicroPrestress;
    double c0;  ///< MPS constant c0 [MPa^-1 * day^-1]
    double c1;  ///< MPS constant c1 (=C1*R*T/M)
    double tS0; ///< MPS tS0 - necessary for the initial value of microprestress (age when the load is applied)
    double kSh; ///< MPS shrinkage parameter. Either this or inithum and finalhum must be given in input record


public:
    B3SolidMaterial(int n, Domain *d) : KelvinChainMaterial(n, d) {
        shMode = B3_NoShrinkage;
    }
    virtual ~B3SolidMaterial() { }

    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);
    //virtual void updateYourself(GaussPoint *gp, TimeStep *tStep);

    virtual void giveShrinkageStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode);

    virtual const char *giveClassName() const { return "B3SolidMaterial"; }
    virtual const char *giveInputRecordName() const { return _IFT_B3SolidMaterial_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *, TimeStep *);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    /// Evaluation of the compliance function of the non-aging solidifying constituent.
    virtual double computeCreepFunction(double t, double t_prime);

protected:
    virtual int hasIncrementalShrinkageFormulation() { return 1; }

    void computeTotalAverageShrinkageStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    /// Evaluation of the shrinkageStrainVector. Shrinkage is fully dependent on humidity rate in given GP
    void computePointShrinkageStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    void predictParametersFrom(double, double, double, double, double, double, double);

    /// Evaluation of the relative volume of the solidified material.
    double computeSolidifiedVolume(TimeStep *tStep);

    /// Evaluation of the flow term viscosity.
    double computeFlowTermViscosity(GaussPoint *gp, TimeStep *tStep);

    double inverse_sorption_isotherm(double w);

    /// Evaluation of characteristic moduli of the non-aging Kelvin chain.
    virtual void computeCharCoefficients(FloatArray &answer, double tStep);

    /// Update of partial moduli of individual chain units
    virtual void updateEparModuli(double tStep);

    virtual void computeCharTimes();

    virtual double giveEModulus(GaussPoint *gp, TimeStep *tStep);

    virtual void giveEigenStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode);

    /**
     * Computes microprestress at given time step and GP.
     * @param gp Gauss point to compute at.
     * @param tStep Time step to compute for.
     * @param option If 0, microprestress is evaluated in the middle of the time step (used for stiffnesses).
     * If 1, MPS is evaluated at the end of the time step. (Used for updating).
     */
    double computeMicroPrestress(GaussPoint *gp, TimeStep *tStep, int option);

    /// Computes initial value of the MicroPrestress
    double giveInitMicroPrestress(void);

    /// Computes relative humidity at given time step and GP
    double giveHumidity(GaussPoint *gp, TimeStep *tStep);

    /// Computes relative humidity increment at given time step and GP
    double giveHumidityIncrement(GaussPoint *gp, TimeStep *tStep);
};
} // end namespace oofem
#endif // b3solidmat_h
