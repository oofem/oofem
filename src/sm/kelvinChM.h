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
 *               Copyright (C) 1993 - 2009   Borek Patzak
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

//   ********************************
//   *** CLASS Kelvin Chain Model ***
//   ********************************

#ifndef kelvinchm_h
#define kelvinchm_h

#include "rheoChM.h"

namespace oofem {
class KelvinChainMaterialStatus : public RheoChainMaterialStatus
{
    /*
     * This class implements associated Material Status to KelvinChainMaterial.
     * It is an attribute of matStatusDictionary at every GaussPoint
     * for which this material
     * DESCRIPTION:
     * Idea used there is that we have variables
     * describing:
     * 1) state at previous equilibrium state (variables without temp)
     * 2) state during searching new equilibrium (variables with temp)
     * when we start search new state from previous equilibrium one we copy
     * non-tem variables into temp ones. And after we reach new equilibrium
     * (now decribed by temp variables) we copy temp-var into non-temp ones
     * (see function updateYourself).
     *
     * variables description:
     *
     *
     * TASK:
     *
     */

protected:

public:
    KelvinChainMaterialStatus(int n, Domain *d, GaussPoint *g, int nunits);
    ~KelvinChainMaterialStatus() {}
    void printOutputAt(FILE *file, TimeStep *tStep);

    /// initialize the status
    virtual void initTempStatus();
    /// update after new equilibrium state reached
    virtual void updateYourself(TimeStep *); // update after new equilibrium state reached

    /// save current context (state) into a stream
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /// restore current context (state) from a stream
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    // definition
    const char *giveClassName() const { return "KelvinChainMaterialStatus"; }
    classType             giveClassID() const
    { return KelvinChainMaterialStatusClass; }
};



//=================================================================================


class KelvinChainMaterial : public RheoChainMaterial
{
    /*
     * This class implements a solidifying Kelvin chain model
     * describing a viscoelastic material.
     *
     * DESCRIPTION
     * TASK
     */

protected:

public:
    KelvinChainMaterial(int n, Domain *d);
    ~KelvinChainMaterial() {}

    // updates MatStatus to the newly reached (equilibrium) state
    void updateYourself(GaussPoint *gp, TimeStep *);

    // identification and auxiliary functions
    virtual int hasNonLinearBehaviour()   { return 0; }
    const char *giveClassName() const { return "KelvinChainMaterial"; }
    classType giveClassID()         const { return KelvinChainMaterialClass; }
    IRResultType initializeFrom(InputRecord *ir);
    void     printYourself();

    // store & restore context functions
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    /**
     * Computes, for the given integration point,
     * the strain vector induced by stress-independent shrinkage
     * @param answer returned strain vector
     * @param form material response form
     * @param gp integration point
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     * @param determines response mode (Total or incremental)
     */

    virtual void  giveShrinkageStrainVector(FloatArray &answer,
                                            MatResponseForm form,
                                            GaussPoint *gp,
                                            TimeStep *atTime,
                                            ValueModeType mode)
    { answer.resize(0); }


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

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;


protected:

    /** if only incremental shrinkage strain formulation is provided, then total shrinkage strain must be tracked
     * in status in order to be able to compute total value. */
    virtual int  hasIncrementalShrinkageFormulation() { return 0; }

    /// evaluation of the creep compliance function
    virtual double  computeCreepFunction(GaussPoint *gp, double ofAge, double atTime) = 0;

    virtual void         computeCharCoefficients(FloatArray &answer, GaussPoint *gp, double);

    /// evaluation of the incremental modulus
    virtual double       giveEModulus(GaussPoint *gp, TimeStep *atTime);

    LinearElasticMaterial *giveLinearElasticMaterial();

    /// this is a modified function from its parent - now it is possible to get a timescale for times < 0.1
    void         generateLogTimeScale(FloatArray &answer, double from, double to, int nsteps,
                                      int fromIncluded = 0);
    const FloatArray &giveDiscreteTimes();
};
} // end namespace oofem
#endif // kelvinchm_h
