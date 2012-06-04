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

#ifndef fluiddynamicmaterial_h
#define fluiddynamicmaterial_h

#include "material.h"
#include "matstatus.h"

#include "flotarry.h"
#include "flotmtrx.h"


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
    /// Equilibrated state vector in reduced form.
    FloatArray deviatoricStressVector;

public:
    /// Constructor - creates new TransportMaterialStatus with number n, belonging to domain d and integration point g.
    FluidDynamicMaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor.
    virtual ~FluidDynamicMaterialStatus() { }

    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    /**
     * Gives the deviatoric stress.
     */
    const FloatArray &giveDeviatoricStressVector() { return deviatoricStressVector; }
    /**
     * Sets the deviatoric stress.
     */
    void letTempDeviatoricStressVectorBe(const FloatArray &v) { deviatoricStressVector = v; }

    virtual const char *giveClassName() const { return "FluidDynamicMaterialStatus"; }
    virtual classType giveClassID() const { return FluidDynamicMaterialStatusClass; }
};


/**
 * Abstract base class for all constitutive models for transport problems. It declares common  services provided
 * by all structural material models. The implementation of these services is partly left on derived classes,
 * which will implement constitutive model dependent part.
 * Some general purpose services are implemented on this level. For details, how to store
 * material model related history variables in integration points, see base class @ref Material documentation.
 *
 * The constitutive model can in general support several material modes (plane stress, plane strain ,... modes).
 * Its capabilities can be examined using hasMaterialModeCapability service.
 * It is generally assumed, that results obtained from constitutive model services are according to
 * valid material mode. This mode is determined from integration point, which is compulsory parameter of all material
 * services.
 */
class FluidDynamicMaterial : public Material
{
public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    FluidDynamicMaterial(int n, Domain *d) : Material(n, d) { }
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
    virtual void computeDeviatoricStressVector(FloatArray &stress_dev, double &epsp_vol, GaussPoint *gp, const FloatArray &eps, double pressure, TimeStep *tStep);
    /**
     * Computes the deviatoric stress vector from given strain.
     * @param answer Deviatoric stress.
     * @param gp Integration point.
     * @param eps Strain-rate.
     * @param tStep Time step.
     */
    virtual void computeDeviatoricStressVector(FloatArray &answer, GaussPoint *gp, const FloatArray &eps, TimeStep *tStep) = 0;


    /**
     * Computes the deviatoric stiffness; @f$ \frac{\partial\sigma_{\mathrm{dev}}}{\partial \epsilon_{\mathrm{dev}}}@f$.
     * @param answer Stiffness matrix.
     * @param mode Mode of result.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual void giveDeviatoricStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;
    /**
     * Computes the tangent @f$ \frac{\partial\sigma_{\mathrm{dev}}}{\partial p}@f$.
     * @param answer Tangent vector.
     * @param mode Mode of result.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual void giveDeviatoricPressureStiffness(FloatArray &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    /**
     * Computes the tangent @f$ \frac{\partial\epsilon_{\mathrm{vol}}}{\partial\epsilon_{\mathrm{dev}}}@f$.
     * @param answer Tangent vector.
     * @param mode Mode of result.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual void giveVolumetricDeviatoricStiffness(FloatArray &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    /**
     * Computes the tangent @f$ \frac{\partial\epsilon_{\mathrm{vol}}}{\partial p}@f$.
     * @param answer Tangent.
     * @param mode Mode of result.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual void giveVolumetricPressureStiffness(double &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    /**
     * Updates internal state of material according to new state vector.
     * @param vec New state vector.
     * @param gp Integration point.
     * @param tStep Solution step.
     */
    virtual void updateInternalState(const FloatArray &vec, GaussPoint *gp, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);

    virtual const char *giveClassName() const { return "FluidDynamicMaterial"; }
    virtual classType giveClassID() const { return FluidDynamicMaterialClass; }
};
} // end namespace oofem
#endif // fluiddynamicmaterial_h
