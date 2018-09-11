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

#ifndef kelvinchsol_h
#define kelvinchsol_h

#include "rheoChM.h"

namespace oofem {
/**
 * This class implements associated Material Status to KelvinChainSolidMaterial,
 * which corresponds to a solidifying Kelvin chain model (framework for creep with aging).
 */
class KelvinChainSolidMaterialStatus : public RheoChainMaterialStatus
{
public:
    KelvinChainSolidMaterialStatus(int n, Domain * d, GaussPoint * g, int nunits);
    virtual ~KelvinChainSolidMaterialStatus() { }
    void printOutputAt(FILE *file, TimeStep *tStep) override;

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    // definition
    const char *giveClassName() const override { return "KelvinChainSolidMaterialStatus"; }
};


/**
 * This class implements a solidifying Kelvin chain model
 * describing a viscoelastic material.
 */
class KelvinChainSolidMaterial : public RheoChainMaterial
{
public:
    KelvinChainSolidMaterial(int n, Domain * d);
    virtual ~KelvinChainSolidMaterial() { }

    void giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep) override;
    void computeHiddenVars(GaussPoint *gp, TimeStep *tStep);

    // identification and auxiliary functions
    const char *giveClassName() const override { return "KelvinChainSolidMaterial"; }
    IRResultType initializeFrom(InputRecord *ir) override;

    void  giveShrinkageStrainVector(FloatArray &answer,
                                    GaussPoint *gp,
                                    TimeStep *tStep,
                                    ValueModeType mode) override
    { answer.clear(); }

    void giveEigenStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override;

    /// Evaluation of the creep compliance function - function useless here
    double computeCreepFunction(double ofAge, double tPrime, GaussPoint *gp, TimeStep *tStep) override;

protected:
    int hasIncrementalShrinkageFormulation() override { return 0; }

    double giveEModulus(GaussPoint *gp, TimeStep *tStep) override;

    /// Evaluation of the relative volume of the solidified material
    virtual double computeSolidifiedVolume(GaussPoint *gp, TimeStep *tStep) = 0;

    /// factors for exponential algorithm
    virtual double computeBetaMu(GaussPoint *gp, TimeStep *tStep, int Mu);
    virtual double computeLambdaMu(GaussPoint *gp, TimeStep *tStep, int Mu);
};
} // end namespace oofem
#endif // kelvinchsol_h
