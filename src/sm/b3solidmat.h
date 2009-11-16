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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   *********************************************
//   *** CLASS B3 Solidification Material Model **
//   *********************************************

#ifndef b3solidmat_h
#define b3solidmat_h

#include "kelvinChM.h"

class B3SolidMaterial : public KelvinChainMaterial
{
    /*
     * This class implements the B3 model for concrete creep and shrinkage
     * based on the solidification theory.
     * The implementation exploits a solidifying Kelvin chain.
     *
     */
protected:
    double t0;
    double w, E28, q1, q2, q3, q4, q5; // predicted data
    enum b3ShModeType { B3_NoShrinkage, B3_AverageShrinkage, B3_PointShrinkage } shMode;
    /// additional parameters for average cross section shrinkage
    double EpsSinf, kt, ks, vs, hum;
    /// additional parameters for free shrinkage at material point
    double es0, r, rprime, at;
    // additional parameters for sorption isotherm (used to compute relative humidity from water content)
    double w_h;       //constant water content (obtained from experiments) w_h [Pedersen, 1990]
    double n;         //constant-exponent (obtained from experiments) n [Pedersen, 1990]
    double a;         //constant (obtained from experiments) A [Pedersen, 1990]
    double talpha;   // thermal dilatation coeff.
    double EspringVal; // elastic modulus of the aging spring (first member of Kelvin chain if retardation spectrum is used)
    int EmoduliMode;  // 0 = analysis of retardation spectrum is used for evaluation of Kelvin units moduli (default)
					// 1 = least-squares method is used for evaluation of Kelvin units moduli
public:
    B3SolidMaterial(int n, Domain *d) : KelvinChainMaterial(n, d) { shMode = B3_NoShrinkage; }
    ~B3SolidMaterial() { }

    virtual void giveShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
                                           GaussPoint *gp, TimeStep *atTime, ValueModeType mode);

    const char *giveClassName()  const { return "B3SolidMaterial"; }
    classType giveClassID()          const { return B3SolidMaterialClass; }
    IRResultType initializeFrom(InputRecord *ir);

    /**
     * Returns a vector of coefficients of thermal dilatation in direction
     * of each material principal (local) axis.
     * @param answer vector of thermal dilatation coefficients
     * @param gp integration point
     * @param tStep time step (most models are able to respond only when atTime is current time step)
     */
    void giveThermalDilatationVector(FloatArray &answer, GaussPoint *, TimeStep *);

protected:

    /** if only incremental shrinkage strain formulation is provided, then total shrinkage strain must be tracked
     * in status in order to be able to compute total value. */
    virtual int  hasIncrementalShrinkageFormulation() { return 1; }

    void computeTotalAverageShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
                                                  GaussPoint *gp, TimeStep *atTime);
    void computeShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
                                      GaussPoint *gp, TimeStep *atTime, ValueModeType mode);
    void  predictParametersFrom(double, double, double, double, double, double, double);

    /// evaluation of the compliance function of the non-aging solidifying constituent
    double computeNonAgingCreepFunction(GaussPoint *gp, double loadDuration);

    /// evaluation of the relative volume of the solidified material
    double computeSolidifiedVolume(GaussPoint *gp, double atAge);

    /// evaluation of the compliance function
    virtual double  computeCreepFunction(GaussPoint *gp, double atTime, double ofAge);

    double inverse_sorption_isotherm(double w);

    /// evaluation of characteristic moduli of the non-aging Kelvin chain
    virtual void         computeCharCoefficients(FloatArray &answer, GaussPoint *gp, double);

    /// evaluation of characteristic times
    virtual void    computeCharTimes();

    /// evaluation of the incremental modulus
    virtual double       giveEModulus(GaussPoint *gp, TimeStep *atTime);

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

    void computeMicroPrestress(GaussPoint *gp, TimeStep *atTime); //PH docasna fce na vyzkouseni nacitani vlhkosti a teploty
};

// Note: There is no associated material status - everything is handled by KelvinChainMaterialStatus

#endif // b3solidmat_h
