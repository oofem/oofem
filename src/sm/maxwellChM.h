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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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
    MaxwellChainMaterialStatus(int n, Domain *d, GaussPoint *g, int nunits);
    virtual ~MaxwellChainMaterialStatus() {}

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    // definition
    virtual const char *giveClassName() const { return "MaxwellChainMaterialStatus"; }
    virtual classType giveClassID() const { return MaxwellChainMaterialStatusClass; }
};



/**
 * This class implements an aging Maxwell chain model
 * describing a viscoelastic material.
 */
class MaxwellChainMaterial : public RheoChainMaterial
{
public:
    MaxwellChainMaterial(int n, Domain *d);
    virtual ~MaxwellChainMaterial() {}

    virtual void updateYourself(GaussPoint *gp, TimeStep *tStep);

    // identification and auxiliary functions
    virtual int hasNonLinearBehaviour()   { return 0; }
    virtual const char *giveClassName() const { return "MaxwellChainMaterial"; }
    virtual classType giveClassID() const { return MaxwellChainMaterialClass; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    // store & restore context functions
    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual void giveShrinkageStrainVector(FloatArray &answer,
                                           MatResponseForm form,
                                           GaussPoint *gp,
                                           TimeStep *tStep,
                                           ValueModeType mode)
    { answer.resize(0); }

    virtual void giveEigenStrainVector(FloatArray &answer, MatResponseForm form,
                                       GaussPoint *gp, TimeStep *tStep, ValueModeType mode);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:
    virtual int hasIncrementalShrinkageFormulation() { return 0; }
    virtual double computeCreepFunction(GaussPoint *gp, double ofAge, double atTime) = 0;
    virtual void computeCharCoefficients(FloatArray &answer, GaussPoint *gp, double atTime);

    virtual double giveEModulus(GaussPoint *gp, TimeStep *atTime);
    LinearElasticMaterial *giveLinearElasticMaterial();
};
} // end namespace oofem
#endif // maxwellchm_h
