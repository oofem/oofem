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

#ifndef kelvinchm_h
#define kelvinchm_h

#include "rheoChM.h"

namespace oofem {
/**
 * This class implements associated Material Status to KelvinChainMaterial.
 */
class KelvinChainMaterialStatus : public RheoChainMaterialStatus
{
public:
    KelvinChainMaterialStatus(GaussPoint *g, int nunits);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    // definition
    const char *giveClassName() const override { return "KelvinChainMaterialStatus"; }
};



/**
 * This class implements a solidifying Kelvin chain model describing a viscoelastic material.
 */
class KelvinChainMaterial : public RheoChainMaterial
{
public:
    KelvinChainMaterial(int n, Domain *d);

    // identification and auxiliary functions
    const char *giveClassName() const override { return "KelvinChainMaterial"; }

    void initializeFrom(InputRecord &ir) override;

    void giveShrinkageStrainVector(FloatArray &answer,
                                   GaussPoint *gp,
                                   TimeStep *tStep,
                                   ValueModeType mode) const override
    { answer.clear(); }

    void giveEigenStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    void giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) override;

    void computeHiddenVars(GaussPoint *gp, TimeStep *tStep);

protected:
    bool hasIncrementalShrinkageFormulation() const override { return false; }

    FloatArray computeCharCoefficients(double tPrime, GaussPoint *gp, TimeStep *tStep) const override;

    double giveEModulus(GaussPoint *gp, TimeStep *tStep) const override;

    //    LinearElasticMaterial *giveLinearElasticMaterial();
};
} // end namespace oofem
#endif // kelvinchm_h
