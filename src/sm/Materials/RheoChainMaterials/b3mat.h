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

#ifndef b3mat_h
#define b3mat_h

#include "maxwellChM.h"

///@name Input fields for B3Material
//@{
#define _IFT_B3Material_Name "b3mat"
#define _IFT_B3Material_mode "mode"
#define _IFT_B3Material_shmode "shmode"
#define _IFT_B3Material_fc "fc"
#define _IFT_B3Material_cc "cc"
#define _IFT_B3Material_wc "w/c"
#define _IFT_B3Material_ac "a/c"
#define _IFT_B3Material_t0 "t0"
#define _IFT_B3Material_es0 "es0"
#define _IFT_B3Material_r "r"
#define _IFT_B3Material_rprime "rprime"
#define _IFT_B3Material_at "at"
#define _IFT_B3Material_wh "w_h"
#define _IFT_B3Material_ncoeff "ncoeff"
#define _IFT_B3Material_a "a"
#define _IFT_B3Material_alpha1 "alpha1"
#define _IFT_B3Material_alpha2 "alpha2"
#define _IFT_B3Material_ks "ks"
#define _IFT_B3Material_hum "hum"
#define _IFT_B3Material_vs "vs"
#define _IFT_B3Material_q1 "q1"
#define _IFT_B3Material_q2 "q2"
#define _IFT_B3Material_q3 "q3"
#define _IFT_B3Material_q4 "q4"
#define _IFT_B3Material_q5 "q5"
#define _IFT_B3Material_kt "kt"
#define _IFT_B3Material_EpsSinf "epssinf"
//@}

namespace oofem {
/**
 * This class implements the B3 model for concrete creep and shrinkage.
 */
class B3Material : public MaxwellChainMaterial
{
protected:
    double t0 = 0.; ///< Age when drying begins (in days)

    double w = 0., E28 = 0.;
    double q1 = 0., q2 = 0., q3 = 0., q4 = 0., q5 = 0.; // predicted data
    enum b3ShModeType { B3_NoShrinkage, B3_AverageShrinkage, B3_PointShrinkage } shMode = B3_NoShrinkage;
    ///@name Additional parameters for average cross section shrinkage.
    //@{
    double EpsSinf = 0., kt = 0., ks = 0., vs = 0., hum = 0.;
    //@}
    ///@name Additional parameters for free shrinkage at material point.
    //@{
    double es0 = 0., r = 0., rprime = 0., at = 0.;
    //@}
    ///@name Additional parameters for absorption isotherm (used to compute relative humidity from water content)
    //@{
    double w_h = 0.;      ///< Constant water content (obtained from experiments) w_h [Pedersen, 1990]
    double n = 0.;        ///< Constant-exponent (obtained from experiments) n [Pedersen, 1990]
    double a = 0.;        ///< Constant (obtained from experiments) A [Pedersen, 1990]
    //@}
public:
    B3Material(int n, Domain *d) : MaxwellChainMaterial(n, d) {}

    void giveShrinkageStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const override;

    const char *giveClassName() const override { return "B3Material"; }
    const char *giveInputRecordName() const override { return _IFT_B3Material_Name; }
    void initializeFrom(InputRecord &ir) override;

    double computeCreepFunction(double t, double t_prime, GaussPoint *gp, TimeStep *tStep) const override;

protected:
    bool hasIncrementalShrinkageFormulation() const override { return true; }

    virtual void computeTotalAverageShrinkageStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) const;
    /// Free shrinkage at material point, requires staggered analysis.
    virtual void computeShrinkageStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const;
    void predictParametersFrom(double, double, double, double, double, double, double);

    /**
     * Function calculates relative humidity from water content (inverse relation form sorption isotherm).
     * Relative humidity (phi) is from range 0.2 - 0.98 !!!
     * Sorption isotherm by C. R. Pedersen (1990), Combined heat and moisture transfer in building constructions,
     * PhD-thesis, Technical University of Denmark, Lingby.
     * @param w Water content (kg/kg).
     */
    double inverse_sorption_isotherm(double w) const;
};
} // end namespace oofem
#endif // b3mat_h
