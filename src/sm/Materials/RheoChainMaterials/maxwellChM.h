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

#ifndef maxwellchm_h
#define maxwellchm_h

#include "rheoChM.h"

namespace oofem {
/**
 * This class implements associated Material Status to MaxwellChainMaterial.
 */
class MaxwellChainMaterialStatus : public RheoChainMaterialStatus
{
public:
    MaxwellChainMaterialStatus(GaussPoint *g, int nunits);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    // definition
    const char *giveClassName() const override { return "MaxwellChainMaterialStatus"; }
};



/**
 * This class implements an aging Maxwell chain model
 * describing a viscoelastic material.
 */
class MaxwellChainMaterial : public RheoChainMaterial
{
public:
    MaxwellChainMaterial(int n, Domain *d);

    // overload thesse function such that computation of hidden vars can be done after the computation of stress
    void giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) override;
    void computeHiddenVars(GaussPoint *gp, TimeStep *tStep);

    // identification and auxiliary functions
    const char *giveClassName() const override { return "MaxwellChainMaterial"; }
    void initializeFrom(InputRecord &ir) override;

    void giveShrinkageStrainVector(FloatArray &answer,
                                   GaussPoint *gp,
                                   TimeStep *tStep,
                                   ValueModeType mode) const override
    { answer.clear(); }

    void giveEigenStrainVector(FloatArray &answer,
                               GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

protected:
    bool hasIncrementalShrinkageFormulation() const override { return false; }
    /**
     * This function computes the moduli of individual Maxwell units
     * such that the corresponding Dirichlet series gives the best
     * approximation of the actual relaxation function.
     *
     * The optimal moduli are obtained using the least-square method,
     * i.e. by minimizing the following functional (tStep = t_0):
     * @f[
     * F=\sum^{k}_{r=1} \left[ \sum^{N}_{\mu=1} E_m(t_0) \exp^{-(t_r-t_0)/\tau_{\mu} - \bar{R}(t_r, t_0)} \right]^2 = min
     * @f]
     *
     * @param[out] answer Array with coefficients
     * @param tStep Age of material when load is applied ???
     */
    FloatArray computeCharCoefficients(double tPrime, GaussPoint *gp, TimeStep *tStep) const override;

    double giveEModulus(GaussPoint *gp, TimeStep *tStep) const override;
    //    LinearElasticMaterial *giveLinearElasticMaterial();
};
} // end namespace oofem
#endif // maxwellchm_h
