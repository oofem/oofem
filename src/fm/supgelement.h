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

#ifndef supgelement_h
#define supgelement_h

#include "fmelement.h"
#include "matresponsemode.h"

///@name Input fields for SUPG elements
//@{
#define _IFT_SUPGElement_bsides "bsides"
#define _IFT_SUPGElement_bcodes "bcodes"
//@}

namespace oofem {
class TimeStep;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;
class Load;

/**
 * General stabilized SUPG/PSPG element for CFD analysis.
 */
class SUPGElement : public FMElement
{
protected:
    /// Array of boundary sides.
    IntArray boundarySides;
    /// Boundary sides codes.
    IntArray boundaryCodes;

    /**
     * Stabilization coefficients, updated for each solution step
     * in updateStabilizationCoeffs()
     */
    double t_supg, t_pspg, t_lsic;

public:
    SUPGElement(int n, Domain * aDomain);
    virtual ~SUPGElement();

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep);
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep);
    virtual void updateStabilizationCoeffs(TimeStep *tStep) { }
    virtual void updateElementForNewInterfacePosition(TimeStep *tStep) { }

    /**
     * Computes acceleration terms (generalized mass matrix with stabilization terms) for momentum balance equations(s).
     */
    virtual void computeAccelerationTerm_MB(FloatMatrix &answer, TimeStep *tStep) = 0;
    /**
     * Computes nonlinear advection terms for momentum balance equations(s).
     */
    virtual void computeAdvectionTerm_MB(FloatArray &answer, TimeStep *tStep) = 0;
    /**
     * Computes the derivative of advection terms for momentum balance equations(s)
     * with respect to nodal velocities.
     */
    virtual void computeAdvectionDerivativeTerm_MB(FloatMatrix &answer, TimeStep *tStep) = 0;
    /**
     * Computes diffusion terms for momentum balance equations(s).
     */
    virtual void computeDiffusionTerm_MB(FloatArray &answer, TimeStep *tStep) = 0;
    /**
     * Computes the derivative of diffusion terms for momentum balance equations(s)
     * with respect to nodal velocities.
     */
    virtual void computeDiffusionDerivativeTerm_MB(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep) = 0;
    /**
     * Computes pressure terms for momentum balance equations(s).
     */
    virtual void computePressureTerm_MB(FloatMatrix &answer, TimeStep *tStep) = 0;
    /**
     * Computes Lhs terms due to boundary conditions - pressure.
     */
    virtual void computeBCLhsPressureTerm_MC(FloatMatrix &answer, TimeStep *tStep);
    /**
     * Computes SLIC stabilization term for momentum balance equation(s).
     */
    virtual void computeLSICStabilizationTerm_MB(FloatMatrix &answer, TimeStep *tStep) = 0;
    /**
     * Computes the linear advection term for mass conservation equation.
     */
    virtual void computeLinearAdvectionTerm_MC(FloatMatrix &answer, TimeStep *tStep) = 0;
    /**
     * Computes advection terms for mass conservation equation.
     */
    virtual void computeAdvectionTerm_MC(FloatArray &answer, TimeStep *tStep) = 0;
    /**
     * Computes the derivative of advection terms for mass conservation equation
     * with respect to nodal velocities.
     */
    virtual void computeAdvectionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *tStep) = 0;
    /**
     * Computes diffusion derivative terms for mass conservation equation.
     */
    virtual void computeDiffusionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *tStep) = 0;
    /**
     * Computes diffusion terms for mass conservation equation.
     */
    virtual void computeDiffusionTerm_MC(FloatArray &answer, TimeStep *tStep) = 0;
    /**
     * Computes acceleration terms for mass conservation equation.
     */
    virtual void computeAccelerationTerm_MC(FloatMatrix &answer, TimeStep *tStep) = 0;
    /**
     * Computes pressure terms for mass conservation equation.
     */
    virtual void computePressureTerm_MC(FloatMatrix &answer, TimeStep *tStep) = 0;
    /**
     * Computes Lhs terms due to boundary conditions - velocity.
     */
    virtual void computeBCLhsTerm_MB(FloatMatrix &answer, TimeStep *tStep);
    /**
     * Computes Lhs terms due to boundary conditions - pressure.
     */
    virtual void computeBCLhsPressureTerm_MB(FloatMatrix &answer, TimeStep *tStep);
    /**
     * Computes Rhs terms due to boundary conditions.
     */
    virtual void computeBCRhsTerm_MB(FloatArray &answer, TimeStep *tStep) = 0;
    /**
     * Computes Rhs terms due to boundary conditions.
     */
    virtual void computeBCRhsTerm_MC(FloatArray &answer, TimeStep *tStep) = 0;
    /**
     * Computes Lhs term due to applied slip with friction bc.
     */
    virtual void computeSlipWithFrictionBCTerm_MB(FloatMatrix &answer, Load *load, int side, TimeStep *tStep) {
        OOFEM_WARNING("computeSlipWithFrictionBCTerm_MB not implemented");
        answer.clear();
    }
    /**
     * Computes Lhs contribution due to applied Penetration bc.
     */
    virtual void computePenetrationWithResistanceBCTerm_MB(FloatMatrix &answer, Load *load, int side, TimeStep *tStep) {
        OOFEM_WARNING("computePenetrationWithResistanceBCTerm_MB not implemented");
        answer.clear();
    }
    /**
     * Computes Lhs contribution due to outflow BC.
     */
    virtual void computeOutFlowBCTerm_MB(FloatMatrix &answer, int side, TimeStep *tStep) {
        OOFEM_WARNING("computeOutFlowBCTerm_MB not implemented");
        answer.clear();
    }


    virtual void computeHomogenizedReinforceTerm_MB(FloatMatrix &answer,  Load *load, TimeStep *tStep) {
        OOFEM_WARNING("computeHomogenizedReinforceTerm_MB");
        answer.clear();
    }
    virtual void computeHomogenizedReinforceTerm_MC(FloatMatrix &answer,  Load *load, TimeStep *tStep) {
        OOFEM_WARNING("computeHomogenizedReinforceTerm_MB");
        answer.clear();
    }


    /// Computes the critical time increment.
    virtual double computeCriticalTimeStep(TimeStep *tStep) = 0;

    // time step termination
    virtual void updateInternalState(TimeStep *tStep);
    virtual int checkConsistency();

#ifdef __OOFEG
    int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                int node, TimeStep *tStep);
#endif

    virtual void giveLocalVelocityDofMap(IntArray &map) { }
    virtual void giveLocalPressureDofMap(IntArray &map) { }

protected:
    virtual void computeDeviatoricStrain(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) = 0;
    virtual void computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);
};
} // end namespace oofem
#endif // supgelement_h
