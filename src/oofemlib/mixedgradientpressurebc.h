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

#ifndef mixedgradientpressurebc_h
#define mixedgradientpressurebc_h

#include "activebc.h"
#include "valuemodetype.h"

namespace oofem {

/**
 * General class for boundary condition that prolongates macroscopic fields to incompressible flow.
 * Both the deviatoric strain rate, and the pressure is controlled. 
 *
 * @see MixedGradientPressureDirichlet For implementation of Dirichlet type
 * @see MixedGradientPressureNeumann For implementation of Neumann type
 * 
 * @author Mikael Öhman
 */
class MixedGradientPressureBC : public ActiveBoundaryCondition
{
public:
    /**
     * Creates boundary condition with given number, belonging to given domain.
     * @param n Boundary condition number.
     * @param d Domain to which new object will belongs.
     */
    MixedGradientPressureBC(int n, Domain *d) : ActiveBoundaryCondition(n,d) { }

    /// Destructor
    virtual ~MixedGradientPressureBC() { };

    /// Not relevant for this boundary condition.
    virtual bcType giveType() const { return UnknownBT; }

    /**
     * Initializes receiver according to object description stored in input record.
     * The input record contains two fields;
     * - devGradient \#columns { d_11 d_22 ... d_21 ... } (required)
     * - pressure p (required)
     * The gradient should be deviatoric and in Voigt notation.
     * The prescribed tensor's columns must be equal to the size of the center coordinates.
     * The size of the center coordinates must be equal to the size of the coordinates in the applied nodes.
     */
    virtual IRResultType initializeFrom(InputRecord *ir);

    /**
     * Computes the size (including pores) by surface integral over the domain
     */
    double domainSize();

    /**
     * Computes the homogenized fields through sensitivity analysis.
     * @param[out] stressDev Computes the homogenized deviatoric stress.
     * @param[out] vol Computes the homogenized volumetric gradient.
     * @param eid Equation ID that fields belong to.
     * @param tStep Time step for which field to obtain.
     */
    virtual void computeFields(FloatArray &sigmaDev, double &vol, EquationID eid, TimeStep *tStep) = 0;

    /**
     * Computes the macroscopic tangents through sensitivity analysis.
     * @param[out] Ed Tangent @f$ \frac{\partial \sigma_{\mathrm{dev}}}{\partial d_{\mathrm{dev}}} @f$.
     * @param[out] Ep Tangent @f$ \frac{\partial \sigma_{\mathrm{dev}}}{\partial p} @f$.
     * @param[out] Cd Tangent @f$ \frac{\partial d_{\mathrm{vol}}}{\partial d_{\mathrm{dev}}} @f$.
     * @param[out] Cp Tangent @f$ \frac{\partial d_{\mathrm{vol}}}{\partial p} @f$.
     * @param eid Equation ID for which the tangents belong.
     * @param tStep Time step for the tangents.
     */
    virtual void computeTangents(FloatMatrix &Ed, FloatArray &Ep, FloatArray &Cd, double &Cp, EquationID eid, TimeStep *tStep) = 0;

    /**
     * Set prescribed pressure.
     * @param p New prescribed pressure.
     */
    virtual void setPrescribedPressure(double p) = 0;

    /**
     * Sets the prescribed tensor from the matrix from given Voigt notation.
     * Assumes use of double values (gamma) for off-diagonal, usually the way for strain in Voigt form.
     * @param ddev Vector in Voigt format.
     */
    virtual void setPrescribedDeviatoricGradientFromVoigt(const FloatArray &ddev) = 0;

    virtual const char *giveClassName() const { return "MixedGradientPressureBC"; }
    virtual classType giveClassID() const { return MixedGradientPressureBCClass; }
};
} // end namespace oofem

#endif // mixedgradientpressurebc_h

