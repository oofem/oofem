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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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
    double microprestress_old = 0.;
    double microprestress_new = 0.;

public:
    B3SolidMaterialStatus(GaussPoint *g, int nunits);

    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    double giveMPS() const { return microprestress_old; }
    void setMPS(double src) { microprestress_new = src; }

    // definition
    const char *giveClassName() const override { return "B3SolidMaterialStatus"; }
};


/**
 * This class implements the B3 model for concrete creep and shrinkage
 * based on the solidification theory.
 * The implementation exploits a solidifying Kelvin chain.
 */
class B3SolidMaterial : public KelvinChainMaterial
{
protected:
    double t0 = 0.;
    double w = 0., E28 = 0., q1 = 0., q2 = 0., q3 = 0., q4 = 0., q5 = 0.; // predicted data

    /// constant equal to one day in time units of analysis (eg. 86400 if the analysis runs in seconds)
    double lambda0 = 0.;

    enum b3ShModeType { B3_NoShrinkage, B3_AverageShrinkage, B3_PointShrinkage } shMode = B3_NoShrinkage;
    /// Additional parameters for average cross section shrinkage
    double EpsSinf = 0., kt = 0., ks = 0., vs = 0., hum = 0.;
    /// Additional parameters for free shrinkage at material point
    double es0 = 0., r = 0., rprime = 0., at = 0.;
    // Additional parameters for sorption isotherm (used to compute relative humidity from water content)
    double w_h = 0.;    ///< Constant water content (obtained from experiments) w_h [Pedersen, 1990]
    double n = 0.;      ///< Constant-exponent (obtained from experiments) n [Pedersen, 1990]
    double a = 0.;      ///< Constant (obtained from experiments) A [Pedersen, 1990]
    mutable double EspringVal = 0.; ///< elastic modulus of the aging spring (first member of Kelvin chain if retardation spectrum is used)
    /**
     * If 0, analysis of retardation spectrum is used for evaluation of Kelvin units moduli (default).
     * If 1, least-squares method is used for evaluation of Kelvin units moduli.
     */
    int EmoduliMode = 0;
    /**
     * If 1, computation exploiting Microprestress solidification theory is done.
     * Default value is 0 = without external fields it can be used for basic creep.
     */
    int MicroPrestress = 0;
    double c0 = 0.;  ///< MPS constant c0 [MPa^-1 * day^-1]
    double c1 = 0.;  ///< MPS constant c1 (=C1*R*T/M)
    double tS0 = 0.; ///< MPS tS0 - necessary for the initial value of microprestress (age when the load is applied)
    double kSh = 0.; ///< MPS shrinkage parameter. Either this or inithum and finalhum must be given in input record


public:
    B3SolidMaterial(int n, Domain *d) : KelvinChainMaterial(n, d) {}

    void giveRealStressVector(FloatArray &answer, GaussPoint *gp,
                              const FloatArray &reducedStrain, TimeStep *tStep) const override;

    void giveShrinkageStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const override;

    const char *giveClassName() const override { return "B3SolidMaterial"; }
    const char *giveInputRecordName() const override { return _IFT_B3SolidMaterial_Name; }
    void initializeFrom(InputRecord &ir) override;

    std::unique_ptr<MaterialStatus> CreateStatus(GaussPoint *gp) const override;

    /// Evaluation of the compliance function of the non-aging solidifying constituent.
    double computeCreepFunction(double t, double t_prime, GaussPoint *gp, TimeStep *tStep) const override;

protected:
    bool hasIncrementalShrinkageFormulation() const override { return true; }

    void computeTotalAverageShrinkageStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) const;

    /// Evaluation of the shrinkageStrainVector. Shrinkage is fully dependent on humidity rate in given GP
    void computePointShrinkageStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) const;

    void predictParametersFrom(double, double, double, double, double, double, double);

    /// Evaluation of the relative volume of the solidified material.
    double computeSolidifiedVolume(TimeStep *tStep) const;

    /// Evaluation of the flow term viscosity.
    double computeFlowTermViscosity(GaussPoint *gp, TimeStep *tStep) const;

    double inverse_sorption_isotherm(double w) const;

    /// Evaluation of characteristic moduli of the non-aging Kelvin chain.
    FloatArray computeCharCoefficients(double tPrime, GaussPoint *gp, TimeStep *tStep) const override;

    /// Update of partial moduli of individual chain units
    void updateEparModuli(double tPrime, GaussPoint *gp, TimeStep *tStep) const override;

    void computeCharTimes() override;

    double giveEModulus(GaussPoint *gp, TimeStep *tStep) const override;

    void giveEigenStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const override;

    /**
     * Computes microprestress at given time step and GP.
     * @param gp Gauss point to compute at.
     * @param tStep Time step to compute for.
     * @param option If 0, microprestress is evaluated in the middle of the time step (used for stiffnesses).
     * If 1, MPS is evaluated at the end of the time step. (Used for updating).
     */
    double computeMicroPrestress(GaussPoint *gp, TimeStep *tStep, int option) const;

    /// Computes initial value of the MicroPrestress
    double giveInitMicroPrestress() const;

    /// Computes relative humidity at given time step and GP
    double giveHumidity(GaussPoint *gp, TimeStep *tStep) const;

    /// Computes relative humidity increment at given time step and GP
    double giveHumidityIncrement(GaussPoint *gp, TimeStep *tStep) const;
};
} // end namespace oofem
#endif // b3solidmat_h
