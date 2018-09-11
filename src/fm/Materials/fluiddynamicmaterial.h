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

#ifndef fluiddynamicmaterial_h
#define fluiddynamicmaterial_h

#include "material.h"
#include "matstatus.h"
#include "floatarray.h"


namespace oofem {
/**
 * This class implements a transport material status information. It is an attribute of
 * a Gauss point. This is only an abstract class, for every instance of material class
 * there should be specialized derived class, which handles are history variables.
 * It only adds attributes common to all "transport problem" material models - the
 * state value vectors (both the temporary and equilibrium) containing the state values
 * in associated integration point. The corresponding services
 * for accessing, setting, initializing and updating these attributes are provided.
 *
 * @see MaterialStatus General description of material status and its role.
 */
class FluidDynamicMaterialStatus : public MaterialStatus
{
protected:
    /// Stress vector in reduced form.
    FloatArray deviatoricStressVector;
    /// Strain vector in reduced form.
    FloatArray deviatoricStrainRateVector;

public:
    /// Constructor - creates new TransportMaterialStatus with number n, belonging to domain d and integration point g.
    FluidDynamicMaterialStatus(int n, Domain * d, GaussPoint * g);
    /// Destructor.
    virtual ~FluidDynamicMaterialStatus() { }

    void printOutputAt(FILE *file, TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    /**
     * Gives the deviatoric stress.
     */
    const FloatArray &giveDeviatoricStressVector() { return deviatoricStressVector; }
    const FloatArray &giveDeviatoricStrainRateVector() { return deviatoricStrainRateVector; }
    /**
     * Sets the deviatoric stress.
     */
    void letDeviatoricStressVectorBe(FloatArray v) { deviatoricStressVector = std :: move(v); }
    void letDeviatoricStrainRateVectorBe(FloatArray v) { deviatoricStrainRateVector = std :: move(v); }
};


/**
 * Abstract base class for all fluid materials. The fluid materials can have compressible response, though most will be incompressible.
 * The models are characterized by having deviatoric strain rate, and pressure as input at each integration point.
 */
class FluidDynamicMaterial : public Material
{
public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    FluidDynamicMaterial(int n, Domain * d) : Material(n, d) { }
    /// Destructor.
    virtual ~FluidDynamicMaterial() { }

    /**
     * Computes the deviatoric stress vector and volumetric strain rate from given deviatoric strain and pressure.
     * Defaults to incompressible function.
     * @param stress_dev Deviatoric stress.
     * @param epsp_vol Volumetric strain-rate (minus the volumetric part of the eps argument).
     * @param gp Integration point.
     * @param eps Strain-rate.
     * @param pressure Pressure.
     * @param tStep Time step.
     */
    virtual void computeDeviatoricStress3D(FloatArray &stress_dev, double &epsp_vol, GaussPoint *gp, const FloatArray &eps, double pressure, TimeStep *tStep);
    virtual void computeDeviatoricStress2D(FloatArray &stress_dev, double &epsp_vol, GaussPoint *gp, const FloatArray &eps, double pressure, TimeStep *tStep);
    /**
     * Computes the deviatoric stress vector from given strain.
     * @param answer Deviatoric stress.
     * @param gp Integration point.
     * @param eps Strain-rate.
     * @param tStep Time step.
     */
    virtual void computeDeviatoricStress3D(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep) = 0;
    virtual void computeDeviatoricStress2D(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep);
    virtual void computeDeviatoricStressAxi(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep);

    /**
     * Computes the deviatoric stiffness; @f$ \frac{\partial\sigma_{\mathrm{dev}}}{\partial \epsilon_{\mathrm{dev}}}@f$.
     * @param answer Stiffness matrix.
     * @param mode Mode of result.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual void computeTangent3D(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    virtual void computeTangent2D(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void computeTangentAxi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    /**
     * Computes the 4 tangents of the compressible material response in 3D.
     * @param dsdd @f$ \frac{\partial\sigma_{\mathrm{dev}}}{\partial \epsilon_{\mathrm{dev}}}@f$.
     * @param dsdp @f$ \frac{\partial\sigma_{\mathrm{dev}}}{\partial p}@f$
     * @param dedd @f$ \frac{\partial\epsilon_{\mathrm{vol}}}{\partial\epsilon_{\mathrm{dev}}}@f$
     * @param dedp @f$ \frac{\partial\epsilon_{\mathrm{vol}}}{\partial p}@f$
     * @param mode Mode of result.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual void computeTangents3D(FloatMatrix &dsdd, FloatArray &dsdp, FloatArray &dedd, double &dedp,
                                         MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void computeTangents2D(FloatMatrix &dsdd, FloatArray &dsdp, FloatArray &dedd, double &dedp,
                                         MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    /**
     * Gives the effective viscosity for the given integration point.
     * @param gp Gauss point of interest.
     * @param tStep Time step.
     * @return The effective viscosity in the point.
     */
    virtual double giveEffectiveViscosity(GaussPoint *gp, TimeStep *tStep) = 0;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
};
} // end namespace oofem
#endif // fluiddynamicmaterial_h
