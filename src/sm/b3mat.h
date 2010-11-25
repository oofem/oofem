/* $Header: /home/cvs/bp/oofem/sm/src/b3mat.h,v 1.4 2003/04/06 14:08:30 bp Exp $ */
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
//   *** CLASS RHEOLOGIC B3 Material Model
//   *********************************************
#ifndef b3mat_h
#define b3mat_h


#include "maxwellChM.h"

namespace oofem {
class B3Material : public MaxwellChainMaterial
{
    /*
     * This class implements the B3 model for concrete creep and shrinkage.
     *
     * DESCRIPTION
     * TASK
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
public:
    B3Material(int n, Domain *d) : MaxwellChainMaterial(n, d) { shMode = B3_NoShrinkage; }
    ~B3Material() { }


    virtual void giveShrinkageStrainVector(FloatArray &answer, MatResponseForm form,
                                           GaussPoint *gp, TimeStep *atTime, ValueModeType mode);

    const char *giveClassName()  const { return "B3Material"; }
    classType giveClassID()          const { return B3MaterialClass; }
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
    virtual double  computeCreepFunction(GaussPoint *gp, double atTime, double ofAge);

    double inverse_sorption_isotherm(double w);
};
} // end namespace oofem
#endif // b3mat_h
