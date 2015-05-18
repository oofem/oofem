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

#ifndef transportmaterial_h
#define transportmaterial_h

#include "material.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "matconst.h"
#include "matstatus.h"
#include "transportelement.h"

namespace oofem {
/**
 * This class implements a transport material status information.
 * When needed, new materials should specialized a derived class from this base.
 * It is attribute of a Gauss point.
 *
 * @see MaterialStatus For general description of material status, and its role.
 */
class TransportMaterialStatus : public MaterialStatus
{
protected:
    FloatArray temp_field; ///< Vector containing the last used field.
    FloatArray temp_gradient; ///< Vector containing the last used gradient.
    FloatArray temp_flux; ///< Vector containing the last computed flux.

    FloatArray field; ///< Vector containing the last equilibrated field. The physical meaning corresponds to temperature, concentration etc.
    FloatArray gradient; ///< Vector containing the last equilibrated gradient. It is the spatial gradient of the field.
    FloatArray flux; ///< Vector containing the last equilibrated flux. The physical meaning corresponds to energy flux, mass flow, etc.

    /// A scalar containing maturity (integration of temperature over time)
    double maturity;

public:
    /// Constructor - creates new TransportMaterialStatus with number n, belonging to domain d and IntegrationPoint g.
    TransportMaterialStatus(int n, Domain * d, GaussPoint * g);
    /// Destructor
    virtual ~TransportMaterialStatus() { }

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    ///@todo REMOVE THESE:
    /// Assigns temporary state vector from a given vector v.
    void letTempStateVectorBe(FloatArray v) { temp_field = std :: move(v); }

    virtual const char *giveClassName() const { return "TransportMaterialStatus"; }

    /// Set gradient.
    void setTempGradient(FloatArray grad);
    /// Set field.
    void setTempField(FloatArray newField);
    /// Set flux.
    void setTempFlux(FloatArray w);

    /// Return last gradient.
    const FloatArray &giveGradient() { return gradient; }
    /// Return last field.
    const FloatArray &giveField() { return field; }
    /// Returns last flux.
    const FloatArray &giveFlux() { return flux; }

    /// Return last gradient.
    const FloatArray &giveTempGradient() { return temp_gradient; }
    /// Return last field.
    const FloatArray &giveTempField() { return temp_field; }
    /// Returns last flux.
    const FloatArray &giveTempFlux() { return temp_flux; }
    /// Returns maturity.
    double giveMaturity() { return maturity; }
};


/**
 * Abstract base class for all constitutive models for transport problems. It declares common services provided
 * by all structural material models. The implementation of these services is partly left on derived classes,
 * which will implement constitutive model dependent part.
 * Some general purpose services are implemented on this level. For details, how to store
 * material model related history variables in integration points, see base class @ref Material documentation.
 */
class TransportMaterial : public Material
{
public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    TransportMaterial(int n, Domain * d) : Material(n, d) { }
    /// Destructor.
    virtual ~TransportMaterial() { }

    /**
     * Returns the flux for the field and its gradient.
     * @todo { Should the field variable just be a scalar? This might change when we rethink the coupled-fields approach.
     * Now its either just [temperature], or [temperature, concentration] so to cover both cases there is a floatarray. }
     * @param answer The flux.
     * @param gp Gauss point.
     * @param grad Gradient of the primary field, usually the main input.
     * @param field The value of the field itself.
     * @param tStep Active time step.
     */
    virtual void giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep) = 0;

    /**
     * Computes characteristic matrix of receiver in given integration point.
     * The algorithm should use temporary or equilibrium  history variables stored in integration point status
     * to compute and return required result.
     * @param answer Contains result.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual void giveCharacteristicMatrix(FloatMatrix &answer,
                                          MatResponseMode mode,
                                          GaussPoint *gp,
                                          TimeStep *tStep) = 0;

    /**
     * Computes the characteristic value of receiver in given integration point, respecting its history.
     * The algorithm should use temporary or equilibrium  history variables stored in integration point status
     * to compute and return required result.
     * @param mode Material response mode.
     * @param gp Integration point.
     * @param tStep Time step (most models are able to respond only when tStep is current time step).
     */
    virtual double giveCharacteristicValue(MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *tStep) = 0;

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
     * @param tStep Solution step.
     * @param mode Determines response mode.
     */
    virtual void computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *tStep, ValueModeType mode)
    { val.clear(); }
    /**
     * Returns positive value of humidity if implemented and enabled in derived material, -1 otherwise.
     */
    virtual double giveHumidity(GaussPoint *gp, ValueModeType mode) { return -1.0; }

    // post-processing
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new TransportMaterialStatus(1, domain, gp); }
};
} // end namespace oofem
#endif // transportmaterial_h
