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
#include "floatarrayf.h"
#include "floatmatrix.h"
#include "matconst.h"
#include "matstatus.h"
#include "tm/Elements/transportelement.h"

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
    double field = 0.; ///< General field (temperature, concentration, etc.).
    FloatArrayF<3> gradient; ///< General gradient.
    FloatArrayF<3> flux; ///< General flux (energy flux, mass flow, etc.)

    double temp_field = 0.; ///< Temp. Primary field.
    FloatArrayF<3> temp_gradient; ///< Temp. Gradient
    FloatArrayF<3> temp_flux; ///< Vector containing the last computed flux.

    double maturity = 0.; ///< A scalar containing maturity (integration of temperature over time)

public:
    TransportMaterialStatus(GaussPoint * g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    const char *giveClassName() const override { return "TransportMaterialStatus"; }

    /// Set gradient.
    void setTempGradient(const FloatArrayF<3> &newGradient) { temp_gradient = newGradient; }
    /// Set field.
    void setTempField(double newField) { temp_field = newField; }
    /// Set flux.
    void setTempFlux(const FloatArrayF<3> &newFlux) { temp_flux = newFlux; }

    /// Return last gradient.
    const FloatArrayF<3> &giveGradient() const { return gradient; }
    /// Return last field.
    double giveField() const { return field; }
    /// Returns last flux.
    const FloatArrayF<3> &giveFlux() const { return flux; }

    /// Return last gradient.
    const FloatArrayF<3> &giveTempGradient() const { return temp_gradient; }
    /// Return last field.
    double giveTempField() const { return temp_field; }
    /// Returns last flux.
    const FloatArrayF<3> &giveTempFlux() const { return temp_flux; }
    /// Returns maturity.
    double giveMaturity() const { return maturity; }
};

/**
 * The temperature is stored im the general "field" value, and this adds the additional humidity field.
 */
class HeMoTransportMaterialStatus : public MaterialStatus
{
protected:
    double temperature = 0.; ///< Temperature.
    FloatArrayF<3> t_gradient; ///< Temperature gradient.
    FloatArrayF<3> t_flux; ///< Heat flux.

    double humidity = 0.; ///< Humidity.
    FloatArrayF<3> h_gradient; ///< Humidity gradient.
    FloatArrayF<3> h_flux; ///< Humidity flux.

    double temp_temperature = 0.; ///< Temp temperature.
    FloatArrayF<3> temp_t_gradient; ///< Temp temperature gradient.
    FloatArrayF<3> temp_t_flux; ///< Temp heat flux.

    double temp_humidity = 0.; ///< Temp humidity.
    FloatArrayF<3> temp_h_gradient; ///< Temp humidity gradient.
    FloatArrayF<3> temp_h_flux; ///< Temp humidity flux.

public:
    HeMoTransportMaterialStatus(GaussPoint * g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void initTempStatus() override;
    void updateYourself(TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    const char *giveClassName() const override { return "TransportMaterialStatus"; }

    /// Set gradient.
    void setTempTemperatureGradient(const FloatArrayF<3> &newGradient) { temp_t_gradient = newGradient; }
    /// Set field.
    void setTempTemperature(double newField) { temp_temperature = newField; }
    /// Set flux.
    void setTempHeatFlux(const FloatArrayF<3> &newFlux) { temp_t_flux = newFlux; }

    /// Return last gradient.
    const FloatArrayF<3> &giveTemperatureGradient() const { return t_gradient; }
    /// Return last field.
    double giveTemperature() const { return temperature; }
    /// Returns last flux.
    const FloatArrayF<3> &giveHeatFlux() const { return t_flux; }

    /// Return last gradient.
    const FloatArrayF<3> &giveTempTemperatureGradient() const { return temp_t_gradient; }
    /// Return last field.
    double giveTempTemperature() const { return temp_temperature; }
    /// Returns last flux.
    const FloatArrayF<3> &giveTempHeatFlux() const { return temp_t_flux; }
    

    /// Set gradient.
    void setTempHumidityGradient(const FloatArrayF<3> &newGradient) { temp_h_gradient = newGradient; }
    /// Set field.
    void setTempHumidity(double newField) { temp_humidity = newField; }
    /// Set flux.
    void setTempHumidityFlux(const FloatArrayF<3> &newFlux) { temp_h_flux = newFlux; }

    /// Return last gradient.
    const FloatArrayF<3> &giveHumidityGradient() const { return h_gradient; }
    /// Return last field.
    double giveHumidity() const { return humidity; }
    /// Returns last flux.
    const FloatArrayF<3> &giveHumidityFlux() const { return h_flux; }

    /// Return last gradient.
    const FloatArrayF<3> &giveTempHumidityGradient() const { return temp_h_gradient; }
    /// Return last field.
    double giveTempHumidity() const { return temp_humidity; }
    /// Returns last flux.
    const FloatArrayF<3> &giveTempHumidityFlux() const { return temp_h_flux; }
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
    virtual void giveFluxVector(FloatArray &answer, GaussPoint *gp, const FloatArray &grad, const FloatArray &field, TimeStep *tStep) const;
    virtual FloatArrayF<3> computeFlux3D(const FloatArrayF<3> &grad, double field, GaussPoint *gp, TimeStep *tStep) const;
    FloatArrayF<2> computeFlux2D(const FloatArrayF<2> &grad, double field, GaussPoint *gp, TimeStep *tStep) const;
    FloatArrayF<1> computeFlux1D(const FloatArrayF<1> &grad, double field, GaussPoint *gp, TimeStep *tStep) const;

    // HeMo:
    virtual std::pair<FloatArrayF<3>, FloatArrayF<3>> computeHeMoFlux3D(const FloatArrayF<3> &grad_t, const FloatArrayF<3> &grad_w, double t, double h, GaussPoint *gp, TimeStep *tStep) const;
    std::pair<FloatArrayF<2>, FloatArrayF<2>> computeHeMoFlux2D(const FloatArrayF<2> &grad_t, const FloatArrayF<2> &grad_w, double t, double h, GaussPoint *gp, TimeStep *tStep) const;
    std::pair<FloatArrayF<1>, FloatArrayF<1>> computeHeMoFlux1D(const FloatArrayF<1> &grad_t, const FloatArrayF<1> &grad_w, double t, double h, GaussPoint *gp, TimeStep *tStep) const;

    //virtual FloatMatrixF<3,3> computeHeMoTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const;
    //virtual FloatMatrixF<2,2> computeHeMoTangent2D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const;
    //virtual FloatMatrixF<1,1> computeHeMoTangent1D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const;

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
                                          TimeStep *tStep) const;

    virtual FloatMatrixF<3,3> computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const = 0;
    FloatMatrixF<2,2> computeTangent2D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const;
    FloatMatrixF<1,1> computeTangent1D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const;

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
                                           TimeStep *tStep) const = 0;

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
    virtual bool hasInternalSource() const { return false; }
    /**
     * Computes the internal source vector of receiver.
     * @param val Contains response.
     * @param gp Integration point.
     * @param tStep Solution step.
     * @param mode Determines response mode.
     */
    virtual void computeInternalSourceVector(FloatArray &val, GaussPoint *gp, TimeStep *tStep, ValueModeType mode) const
    { val.clear(); }
    /**
     * Returns positive value of humidity if implemented and enabled in derived material, -1 otherwise.
     */
    virtual double giveHumidity(GaussPoint *gp, ValueModeType mode) const { return -1.0; }

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;

    MaterialStatus *CreateStatus(GaussPoint *gp) const override { return new TransportMaterialStatus(gp); }
};
} // end namespace oofem
#endif // transportmaterial_h
