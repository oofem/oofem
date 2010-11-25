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

#ifndef binghamfluid2_h
#define binghamfluid2_h

#include "fluiddynamicmaterial.h"
#include "flotarry.h"
#include "flotmtrx.h"

#include "matconst.h"
#include "structuralelement.h"
#include "matstatus.h"

namespace oofem {
class GaussPoint;

/**
 * Class representing material status for Bingham material
 */
class BinghamFluidMaterial2Status : public FluidDynamicMaterialStatus
{
protected:

    /// magnitude of deviatoric strains
    double devStrainMagnitude, temp_devStrainMagnitude;
    /// magnitude of deviatoric stresses
    double devStressMagnitude, temp_devStressMagnitude;
    /// devaitoric stresses and strains
    FloatArray deviatoricStrainVector, temp_deviatoricStrainVector;  // reduced form

public:
    /// Constructor - creates new BinghamFluidMaterial2Status with number n, belonging to domain d and IntegrationPoint g.
    BinghamFluidMaterial2Status(int n, Domain *d, GaussPoint *g);
    /// Destructor
    ~BinghamFluidMaterial2Status() { }

    /// Print receiver's output to given stream.
    void   printOutputAt(FILE *, TimeStep *);

    /**
     * Initializes the temporary internal variables (stresss and strains vectors),
     * describing the current state according to
     * previously reached equilibrium internal variables.
     */
    virtual void initTempStatus();
    /**
     * Update equilibrium history variables according to temp-variables.
     * Invoked, after new equilibrium state has been reached.
     */
    virtual void updateYourself(TimeStep *); // update after new equilibrium state reached


    /**
     * Stores context of receiver into given stream (the equilibriun stress and strains vectors are stored).
     * Generally, only non-temp internal history variables should be stored.
     * @param stream stream where to write data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores context of receiver from given stream (the equilibriun stress and strains vectors are restored).
     * @param stream stream where to read data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    // saves current context(state) into stream

    const double giveTempDevStressMagnitude() { return temp_devStressMagnitude; }
    const double giveTempDevStrainMagnitude() { return temp_devStrainMagnitude; }
    const double giveDevStressMagnitude() { return devStressMagnitude; }
    const double giveDevStrainMagnitude() { return devStrainMagnitude; }

    void letTempDevStrainMagnitudeBe(double _val) { temp_devStrainMagnitude = _val; }
    void letTempDevStressMagnitudeBe(double _val) { temp_devStressMagnitude = _val; }

    const FloatArray &giveDeviatoricStrainVector()         { return deviatoricStrainVector; }
    const FloatArray &giveTempDeviatoricStrainVector()     { return temp_deviatoricStrainVector; }
    void         letTempDeviatoricStrainVectorBe(const FloatArray &v)
    { temp_deviatoricStrainVector = v; }


    /// Returns "BinghamFluidMaterial2Status" string - class name of the receiver.
    const char *giveClassName() const { return "BinghamFluidMaterialStatus"; }
    /// Returns TransportMaterialStatusClass - classType id of receiver.
    classType                giveClassID() const
    { return BinghamFluidMaterialStatusClass; }
};



/**
 * Constitutive model of Bingham fluid for concentrated suspensions and pastes.
 * This is the simplest two-constant model, with yield stress and viscosity as parameters.
 */
class BinghamFluidMaterial2 : public FluidDynamicMaterial
{
protected:
    /// viscosity
    double mu_0;
    /// yield stress
    double tau_0;
    double tau_c;
    double mu_inf;
public:

    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n material number
     * @param d domain to which new material will belong
     */
    BinghamFluidMaterial2(int n, Domain *d) : FluidDynamicMaterial(n, d) { mu_inf = 1.e6; }
    /// Destructor.
    ~BinghamFluidMaterial2()                { }

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
    virtual double   give(int aProperty, GaussPoint *);
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
    const char *giveClassName() const { return "BinghamFluidMaterial2"; }
    /// Returns classType id of receiver.
    classType giveClassID()         const { return BinghamFluidMaterialClass; }

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
    double  computeActualViscosity(double Tau, double shearRate);
    double  computeDevStrainMagnitude(MaterialMode mmode, const FloatArray &epsd);
    double  computeDevStressMagnitude(MaterialMode mmode, const FloatArray &sigd);
    void    computeDeviatoricStrain(FloatArray &answer, const FloatArray &eps, MaterialMode mmode);
    void    computeDeviatoricStress(FloatArray &answer, const FloatArray &deps,
                                    double _nu, MaterialMode mmode);

    void    __debug(GaussPoint *gp, TimeStep *atTime);
};
} // end namespace oofem
#endif // binghamfluid2_h
