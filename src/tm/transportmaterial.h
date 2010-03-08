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

#ifndef transportmaterial_h
#define transportmaterial_h

#include "material.h"
#include "dictionr.h"
#include "flotarry.h"
#include "flotmtrx.h"

#include "matconst.h"
#include "structuralelement.h"
#include "matstatus.h"
#include "transportelement.h"

namespace oofem {
class GaussPoint;

/**
 * This class implements a transport material status information. It is atribute of
 * gaussPoint. This is only an abstract class, for every instance of material class
 * there should be specialized derived class, which handles are history variables.
 * It only adds attributes common to all "transport problem" material models - the
 * state value vectors (both the temp and equlibrium) containing the state values
 * in associated integration point. The corresponding services
 * for accessing, setting, initializing and updating these attributes are provided.
 *
 * For general description of material status, and its role @see MaterialStatus class.
 */
class TransportMaterialStatus : public MaterialStatus
{
protected:

    /// Equilibrated state vector in reduced form. The physical meaning corresponds to temperature, concentration etc.
    FloatArray stateVector;
    /// Temporary state vector in a reduced form, used mainly in a nonlinear analysis
    FloatArray tempStateVector;

public:
    /// Constructor - creates new TransportMaterialStatus with number n, belonging to domain d and IntegrationPoint g.
    TransportMaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    ~TransportMaterialStatus() { }

    /// Print receiver's output to given stream.
    void   printOutputAt(FILE *, TimeStep *);

    /**
     * Initializes temporary internal variables (state vectors).
     * Assign previously reached equilibrium internal variables to them.
     */
    virtual void initTempStatus();
    /**
     * Update equilibrium history variables according to temp-variables.
     * Invoked, after new equilibrium state has been reached.
     */
    virtual void updateYourself(TimeStep *); // update after new equilibrium state reached


    /**
     * Stores context of receiver into given stream (the equilibrated state variable).
     * Generally, only non-temp internal history variables should be stored.
     * @param stream stream where to write data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores context of receiver from given stream (the equilibrated state variable).
     * @param stream stream where to read data
     * @param mode determines ammount of info required in stream (state, definition,...)
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    // saves current context(state) into stream

    /// Returns the const pointer to receiver's stateVector
    const FloatArray &giveStateVector()         { return stateVector; }
    /// Returns the const pointer to receiver's tempStateVector
    const FloatArray &giveTempStateVector()      { return tempStateVector; }
    /// Assigns tempStateVector from a given vector v
    void         letTempStateVectorBe(const FloatArray &v)
    { tempStateVector = v; }

    /// Returns "TransportMaterialStatus" - class name of the receiver.
    const char *giveClassName() const { return "TransportMaterialStatus"; }
    /// Returns TransportMaterialStatusClass - classType id of receiver.
    classType                giveClassID() const
    { return TransportMaterialStatusClass; }
};


/**
 * Abstract base class for all constitutive models for transport problems. It declares common  services provided
 * by all structural material models. The implementation of these services is partly left on derived classes,
 * which will implement constitutive model dependent part.
 * Some general purpose services are implemented on this level. For details, how to store
 * material model related history variables in integration points, see base class \ref Material documentation.
 *
 * The constitutive model can in general support several material modes (plane stress, plane strain ,... modes).
 * Its capabilities can be examined using hasMaterialModeCapability  service.
 * It is generally assumed, that results obtained from constitutive model services are according to
 * valid material mode. This mode is determined from integration point, which is compulsory parameter of all material
 * services.
 */
class TransportMaterial : public Material
{
protected:
public:

    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n material number
     * @param d domain to which new material will belong
     */
    TransportMaterial(int n, Domain *d) : Material(n, d) { }
    /// Destructor.
    ~TransportMaterial()                { }

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
                                           TimeStep *atTime) = 0;

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
                                            TimeStep *atTime) = 0;
    /**
     * Updates internal state of material according to new state vector.
     * @param stateVec new state vector
     * @param gp integration point
     * @param tStep solution step
     */
    virtual void updateInternalState(const FloatArray &stateVec, GaussPoint *gp, TimeStep *);
    /**
     * Returns nonzero if receiver genarets internal source of state variable(s), zero otherwise.
     */
    virtual int hasInternalSource() { return 0; }
    /**
     * Computes the internal source vector of receiver.
     * @param val contains response
     * @param gp integration point
     * @param atTime solution step
     * @param mode determines response mode
     */
    virtual void computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
    { val.resize(0); }
    /**
     * Returns positive value of humidity if implemented and enabled in derived material, -1 otherwise.
     */
    virtual double giveHumidity(GaussPoint *gp, ValueModeType mode) { return -1.0; }

    /**
     * Request material extension.
     * @param ext material extension tested
     * @return nonzero if implemented
     */
    virtual int testMaterialExtension(MaterialExtension ext) { return ( ( ext == Material_TransportCapability ) ? 1 : 0 ); }
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "TransportMaterial"; }
    /// Returns classType id of receiver.
    classType giveClassID()         const { return TransportMaterialClass; }

    // post-processing
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);

#ifdef __OOFEG
#endif
};
} // end namespace oofem
#endif // transportmaterial_h
