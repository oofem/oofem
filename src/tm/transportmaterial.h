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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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
#include "flotarry.h"
#include "flotmtrx.h"

#include "matconst.h"
#include "structuralelement.h"
#include "matstatus.h"
#include "transportelement.h"

namespace oofem {

/**
 * This class implements a transport material status information. It is attribute of
 * a Gauss point. This is only an abstract class, for every instance of material class
 * there should be specialized derived class, which handles are history variables.
 * It only adds attributes common to all "transport problem" material models - the
 * state value vectors (both the temp and equilibrium) containing the state values
 * in associated integration point. The corresponding services
 * for accessing, setting, initializing and updating these attributes are provided.
 *
 * @see MaterialStatus For general description of material status, and its role.
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
    virtual ~TransportMaterialStatus() { }

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    /// Returns the const pointer to receiver's state vector.
    const FloatArray &giveStateVector() { return stateVector; }
    /// Returns the const pointer to receiver's temporary state vector.
    const FloatArray &giveTempStateVector() { return tempStateVector; }
    /// Assigns temporary state vector from a given vector v.
    void letTempStateVectorBe(const FloatArray &v) { tempStateVector = v; }

    virtual const char *giveClassName() const { return "TransportMaterialStatus"; }
    virtual classType giveClassID() const { return TransportMaterialStatusClass; }
};


/**
 * Abstract base class for all constitutive models for transport problems. It declares common  services provided
 * by all structural material models. The implementation of these services is partly left on derived classes,
 * which will implement constitutive model dependent part.
 * Some general purpose services are implemented on this level. For details, how to store
 * material model related history variables in integration points, see base class @ref Material documentation.
 *
 * The constitutive model can in general support several material modes (plane stress, plane strain ,... modes).
 * Its capabilities can be examined using hasMaterialModeCapability  service.
 * It is generally assumed, that results obtained from constitutive model services are according to
 * valid material mode. This mode is determined from integration point, which is compulsory parameter of all material
 * services.
 */
class TransportMaterial : public Material
{
public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    TransportMaterial(int n, Domain *d) : Material(n, d) { }
    /// Destructor.
    virtual ~TransportMaterial() { }

    virtual void giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep) {};

    virtual void giveCharacteristicMatrix(FloatMatrix &answer,
                                          MatResponseForm form,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *atTime) = 0;

    virtual double giveCharacteristicValue(MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime) = 0;

    /**
     * Updates internal state of material according to new state vector.
     * @param state New state vector.
     * @param gp Integration point.
     * @param tStep Solution step.
     */
    virtual void updateInternalState(const FloatArray &state, GaussPoint *gp, TimeStep *tStep);

    /**
     * Returns nonzero if receiver generates internal source of state variable(s), zero otherwise.
     */
    virtual int hasInternalSource() { return 0; }
    /**
     * Computes the internal source vector of receiver.
     * @param val Contains response.
     * @param gp Integration point.
     * @param atTime Solution step.
     * @param mode Determines response mode.
     */
    virtual void computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *atTime, ValueModeType mode)
    { val.resize(0); }
    /**
     * Returns positive value of humidity if implemented and enabled in derived material, -1 otherwise.
     */
    virtual double giveHumidity(GaussPoint *gp, ValueModeType mode) { return -1.0; }

    virtual const char *giveClassName() const { return "TransportMaterial"; }
    virtual classType giveClassID() const { return TransportMaterialClass; }

    // post-processing
    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);
};
} // end namespace oofem
#endif // transportmaterial_h
