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
    KelvinChainMaterialStatus(int n, Domain * d, GaussPoint * g, int nunits);
    virtual ~KelvinChainMaterialStatus() { }

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    // definition
    virtual const char *giveClassName() const { return "KelvinChainMaterialStatus"; }
};



/**
 * This class implements a solidifying Kelvin chain model describing a viscoelastic material.
 */
class KelvinChainMaterial : public RheoChainMaterial
{
public:
    KelvinChainMaterial(int n, Domain * d);
    virtual ~KelvinChainMaterial() { }

    // identification and auxiliary functions
    virtual int hasNonLinearBehaviour() { return 0; }
    virtual const char *giveClassName() const { return "KelvinChainMaterial"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void  giveShrinkageStrainVector(FloatArray &answer,
                                            GaussPoint *gp,
                                            TimeStep *tStep,
                                            ValueModeType mode)
    { answer.clear(); }

    virtual void  giveEigenStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, ValueModeType mode);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

    virtual void giveRealStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep);

    void computeHiddenVars(GaussPoint *gp, TimeStep *tStep);

protected:
    virtual int hasIncrementalShrinkageFormulation() { return 0; }

    virtual void computeCharCoefficients(FloatArray &answer, double tStep);

    virtual double giveEModulus(GaussPoint *gp, TimeStep *tStep);

    LinearElasticMaterial *giveLinearElasticMaterial();
};
} // end namespace oofem
#endif // kelvinchm_h
