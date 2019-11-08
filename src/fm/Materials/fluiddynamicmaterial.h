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
#include "floatarrayf.h"
#include "floatmatrixf.h"

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
    FloatArrayF<6> deviatoricStressVector;
    /// Strain vector in reduced form.
    FloatArrayF<6> deviatoricStrainRateVector;

public:
    /// Constructor - creates new TransportMaterialStatus with number n, belonging to domain d and integration point g.
    FluidDynamicMaterialStatus(GaussPoint * g);

    void printOutputAt(FILE *file, TimeStep *tStep) const override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    /**
     * Gives the deviatoric stress.
     */
    const FloatArrayF<6> &giveDeviatoricStressVector() const { return deviatoricStressVector; }
    const FloatArrayF<6> &giveDeviatoricStrainRateVector() const { return deviatoricStrainRateVector; }
    /**
     * Sets the deviatoric stress.
     */
    void letDeviatoricStressVectorBe(const FloatArrayF<6> &v) { deviatoricStressVector = v; }
    void letDeviatoricStrainRateVectorBe(const FloatArrayF<6> &v) { deviatoricStrainRateVector = v; }

    const char *giveClassName() const override { return "FluidDynamicMaterialStatus"; }
};


/**
 * Abstract base class for all fluid materials. The fluid materials can have compressible response, though most will be incompressible.
 * The models are characterized by having deviatoric strain rate, and pressure as input at each integration point.
 */
class FluidDynamicMaterial : public Material
{
public:
    template<int N> struct Tangents {
        FloatMatrixF<N,N> dsdd;
        FloatArrayF<N> dsdp;
        FloatArrayF<N> dedd;
        double dedp;
    };

    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n Material number.
     * @param d Domain to which new material will belong.
     */
    FluidDynamicMaterial(int n, Domain * d) : Material(n, d) { }

    /**
     * Computes the deviatoric stress vector and volumetric strain rate from given deviatoric strain and pressure.
     * Defaults to incompressible function.
     * @param gp Integration point.
     * @param eps Strain-rate.
     * @param pressure Pressure.
     * @param tStep Time step.
     * @return Deviatoric stress and volumetric strain-rate (minus the volumetric part of the eps argument).
     */
    virtual std::pair<FloatArrayF<6>, double> computeDeviatoricStress3D(const FloatArrayF<6> &eps, double pressure,
                                                                        GaussPoint *gp, TimeStep *tStep) const;
    virtual std::pair<FloatArrayF<3>, double> computeDeviatoricStress2D(const FloatArrayF<3> &eps, double pressure,
                                                                        GaussPoint *gp, TimeStep *tStep) const;
    /**
     * Computes the deviatoric stress vector from given strain.
     * @param answer Deviatoric stress.
     * @param gp Integration point.
     * @param eps Strain-rate.
     * @param tStep Time step.
     */
    virtual FloatArrayF<6> computeDeviatoricStress3D(const FloatArrayF<6> &eps, GaussPoint *gp, TimeStep *tStep) const = 0;
    virtual FloatArrayF<3> computeDeviatoricStress2D(const FloatArrayF<3> &eps, GaussPoint *gp, TimeStep *tStep) const;
    virtual FloatArrayF<4> computeDeviatoricStressAxi(const FloatArrayF<4> &eps, GaussPoint *gp, TimeStep *tStep) const;

    /**
     * Computes the deviatoric stiffness; @f$ \frac{\partial\sigma_{\mathrm{dev}}}{\partial \epsilon_{\mathrm{dev}}}@f$.
     * @param answer Stiffness matrix.
     * @param mode Mode of result.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual FloatMatrixF<6,6> computeTangent3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const = 0;
    virtual FloatMatrixF<3,3> computeTangent2D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const;
    virtual FloatMatrixF<4,4> computeTangentAxi(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const;

    /**
     * Computes the 4 tangents of the compressible material response in 3D.
     * - dsdd @f$ \frac{\partial\sigma_{\mathrm{dev}}}{\partial \epsilon_{\mathrm{dev}}}@f$.
     * - dsdp @f$ \frac{\partial\sigma_{\mathrm{dev}}}{\partial p}@f$
     * - dedd @f$ \frac{\partial\epsilon_{\mathrm{vol}}}{\partial\epsilon_{\mathrm{dev}}}@f$
     * - dedp @f$ \frac{\partial\epsilon_{\mathrm{vol}}}{\partial p}@f$
     * @param mode Mode of result.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual Tangents<6> computeTangents3D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const;
    virtual Tangents<3> computeTangents2D(MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) const;

    /**
     * Gives the effective viscosity for the given integration point.
     * @param gp Gauss point of interest.
     * @param tStep Time step.
     * @return The effective viscosity in the point.
     */
    virtual double giveEffectiveViscosity(GaussPoint *gp, TimeStep *tStep) const = 0;

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
};
} // end namespace oofem
#endif // fluiddynamicmaterial_h
