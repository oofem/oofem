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

#ifndef b3mat_h
#define b3mat_h

#include "maxwellChM.h"

namespace oofem {
/**
 * This class implements the B3 model for concrete creep and shrinkage.
 */
class B3Material : public MaxwellChainMaterial
{
protected:
    double t0; ///< Age when drying begins (in days)

    double w, E28;
    double q1, q2, q3, q4, q5; // predicted data
    enum b3ShModeType { B3_NoShrinkage, B3_AverageShrinkage, B3_PointShrinkage } shMode;
    ///@name Additional parameters for average cross section shrinkage.
    //@{
    double EpsSinf, kt, ks, vs, hum;
    //@}
    ///@name Additional parameters for free shrinkage at material point.
    //@{
    double es0, r, rprime, at;
    //@}
    ///@name Additional parameters for absorption isotherm (used to compute relative humidity from water content)
    //@{
    double w_h;      ///< Constant water content (obtained from experiments) w_h [Pedersen, 1990]
    double n;        ///< Constant-exponent (obtained from experiments) n [Pedersen, 1990]
    double a;        ///< Constant (obtained from experiments) A [Pedersen, 1990]
    double talpha;   ///< Thermal dilatation coeff.
    //@}
public:
    B3Material(int n, Domain *d) : MaxwellChainMaterial(n, d) { shMode = B3_NoShrinkage; }
    virtual ~B3Material() { }

    virtual void giveShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
                                           GaussPoint *gp, TimeStep *atTime, ValueModeType mode);

    virtual const char *giveClassName() const { return "B3Material"; }
    virtual classType giveClassID() const { return B3MaterialClass; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

protected:
    virtual int hasIncrementalShrinkageFormulation() { return 1; }

    virtual void computeTotalAverageShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
                                                  GaussPoint *gp, TimeStep *atTime);
    /// Free shrinkage at material point, requires staggered analysis.
    virtual void computeShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
                                      GaussPoint *gp, TimeStep *atTime, ValueModeType mode);
    void predictParametersFrom(double, double, double, double, double, double, double);
    virtual double computeCreepFunction(GaussPoint *gp, double atTime, double ofAge);

    /**
     * Function calculates relative humidity from water content (inverse relation form sorption isotherm).
     * Relative humidity (phi) is from range 0.2 - 0.98 !!!
     * Sorption isotherm by C. R. Pedersen (1990), Combined heat and moisture transfer in building constructions,
     * PhD-thesis, Technical University of Denmark, Lingby.
     * @param w Water content (kg/kg).
     */
    double inverse_sorption_isotherm(double w);
};
} // end namespace oofem
#endif // b3mat_h
