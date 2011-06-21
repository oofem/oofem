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

//   ********************************************
//   *** CLASS Solidifying Kelvin Chain Model ***
//   ********************************************

#ifndef kelvinchsol_h
#define kelvinchsol_h

#include "rheoChM.h"

namespace oofem {
class KelvinChainSolidMaterialStatus : public RheoChainMaterialStatus
{
    /**
     * This class implements associated Material Status to KelvinChainSolidMaterial,
     * which corresponds to a solidifying Kelvin chain model (framework for creep with aging).
     **/

protected:

public:
    KelvinChainSolidMaterialStatus(int n, Domain *d, GaussPoint *g, int nunits);
    ~KelvinChainSolidMaterialStatus() {}
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
    const char *giveClassName() const { return "KelvinChainSolidMaterialStatus"; }
    classType             giveClassID() const
    { return KelvinChainSolidMaterialStatusClass; }
};



//=================================================================================


class KelvinChainSolidMaterial : public RheoChainMaterial
{
    /**
     * This class implements a solidifying Kelvin chain model
     * describing a viscoelastic material.
     */

protected:

public:
    KelvinChainSolidMaterial(int n, Domain *d);
    ~KelvinChainSolidMaterial() {}

    // updates MatStatus to the newly reached (equilibrium) state
    void updateYourself(GaussPoint *gp, TimeStep *);

    // identification and auxiliary functions
    virtual int hasNonLinearBehaviour()   { return 0; }
    const char *giveClassName() const { return "KelvinChainSolidMaterial"; }
    classType giveClassID()         const { return KelvinChainSolidMaterialClass; }
    IRResultType initializeFrom(InputRecord *ir);

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

    /// evaluation of the creep compliance function - function useless here
    virtual double  computeCreepFunction(GaussPoint *gp, double ofAge, double atTime);

    /// evaluation of the incremental modulus
    virtual double       giveEModulus(GaussPoint *gp, TimeStep *atTime);

    /// evaluation of the relative volume of the solidified material
    virtual double computeSolidifiedVolume(GaussPoint *gp, TimeStep *atTime) = 0;

    /// factors for exponential algorithm
    virtual double computeBetaMu(GaussPoint *gp, TimeStep *atTime, double Mu);
    virtual double computeLambdaMu(GaussPoint *gp, TimeStep *atTime, double Mu);

    LinearElasticMaterial *giveLinearElasticMaterial();
};
} // end namespace oofem
#endif // kelvinchsol_h
