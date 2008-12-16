/* $Header: /home/cvs/bp/oofem/tm/src/transportmaterial.h,v 1.1 2003/04/14 16:01:40 bp Exp $ */
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

#ifndef twofluidmaterial_h
#define twofluidmaterial_h

#include "fluiddynamicmaterial.h"
#include "flotarry.h"
#include "flotmtrx.h"

#include "matconst.h"
#include "structuralelement.h"
#include "matstatus.h"

class GaussPoint;

/**
 * Material coupling the behaviour of two particular materials based on
 * rule of mixture. The weighting factor is VOF fraction.
 */
class TwoFluidMaterial : public FluidDynamicMaterial
{
protected:
    int masterMat;
    int slaveMaterial [ 2 ];
public:

    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n material number
     * @param d domain to which new material will belong
     */
    TwoFluidMaterial(int n, Domain *d) : FluidDynamicMaterial(n, d) { masterMat = 0; }
    /// Destructor.
    ~TwoFluidMaterial()                { }

    /**
     * Computes the characteristic matrix of receiver in given integration point, respecting its history.
     * The algorithm should use temporary or equlibrium  history variables stored in integration point status
     * to compute and return required result.
     * @param answer contains result
     * @param form material response form
     * @param mode  material response mode
     * @param gp integration point
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    virtual void  giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime) { }

    /**
     * Computes the characteristic value of receiver in given integration point, respecting its history.
     * The algorithm should use temporary or equlibrium  history variables stored in integration point status
     * to compute and return required result.
     * @param mode material response mode
     * @param gp integration point
     * @param atTime time step (most models are able to respond only when atTime is current time step)
     */
    virtual double  giveCharacteristicValue(MatResponseMode mode,
                                            GaussPoint *gp,
                                            TimeStep *atTime);

    /**
     * Computes devatoric stress vector from given strain
     */
    virtual void computeDeviatoricStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep);

    /**
     * Computes Deviatoric stiffness (derivative of deviatoric stress tensor with respect to strain)
     */
    virtual void giveDeviatoricStiffnessMatrix(FloatMatrix & answer, MatResponseMode, GaussPoint * gp,
                                               TimeStep * atTime);
    /**
     * Returns the value of material property 'aProperty'. Property must be identified
     * by unique int id.
     * @param aProperty id of peroperty requested
     * @param gp integration point
     * @return property value
     */
    virtual double   give(int aProperty, GaussPoint* gp);
    /**
     * Initializes receiver acording to object description stored in input record.
     * The density of material is read into property dictionary (keyword 'd')
     */
    IRResultType initializeFrom(InputRecord *ir);
    /** Setups the input record string of receiver
     * @param str string to be filled by input record
     * @param keyword print record keyword (default true)
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    /**
     * Tests, if material supports material mode.
     * @param mode required material mode
     * @return nonzero if supported, zero otherwise
     */
    virtual int hasMaterialModeCapability(MaterialMode mode);
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "TwoFluidMaterial"; }
    /// Returns classType id of receiver.
    classType giveClassID()         const { return TwoFluidMaterialClass; }

    /** Allows programmer to test some internal data, before computation begins.
     *  For example, one may use this function, to ensure that element has material with
     *  required capabilities is assigned to element. This must be done after all
     *  mesh components are instanciated.
     *  @return nonzero if receiver check is o.k. */
    virtual int    checkConsistency();
#ifdef __OOFEG
#endif

    /**
     * Creates new copy of associated status and inserts it into given integration point.
     * @param gp Integration point where newly created status will be stored.
     * @return reference to new status.
     */
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const;

protected:

    FluidDynamicMaterial *giveMaterial(int i) const { return static_cast< FluidDynamicMaterial * >( domain->giveMaterial(slaveMaterial [ i ]) ); }
    double giveTempVOF(GaussPoint *gp);
};


#endif // twofluidmaterial_h





