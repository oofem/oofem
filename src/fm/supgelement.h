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

#ifndef supgelement_h
#define supgelement_h

#include "fmelement.h"
#include "fluiddynamicmaterial.h"

namespace oofem {
class TimeStep;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * General stabilized SUPG/PSPG element for CFD analysis.
 */
class SUPGElement : public FMElement
{
protected:
    /**
     * Stabilization coefficients, updated for each solution step
     * in updateStabilizationCoeffs()
     */
    double t_supg, t_pspg, t_lsic;

public:
    SUPGElement(int n, Domain *aDomain);
    virtual ~SUPGElement();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep);
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep);
    virtual double giveCharacteristicValue(CharType type, TimeStep *tStep);
    virtual void updateStabilizationCoeffs(TimeStep *tStep) { }
    virtual void updateElementForNewInterfacePosition(TimeStep *tStep) { }

    /**
     * Computes acceleration terms (generalized mass matrix with stabilization terms) for momentum balance equations(s).
     */
    virtual void computeAccelerationTerm_MB(FloatMatrix &answer, TimeStep *atTime) = 0;
    /**
     * Computes nonlinear advection terms for momentum balance equations(s).
     */
    virtual void computeAdvectionTerm_MB(FloatArray &answer, TimeStep *atTime) = 0;
    /**
     * Computes the derivative of advection terms for momentum balance equations(s)
     * with respect to nodal velocities.
     */
    virtual void computeAdvectionDerivativeTerm_MB(FloatMatrix &answer, TimeStep *atTime) = 0;
    /**
     * Computes diffusion terms for momentum balance equations(s).
     */
    virtual void computeDiffusionTerm_MB(FloatArray &answer, TimeStep *atTime) = 0;
    /**
     * Computes the derivative of diffusion terms for momentum balance equations(s)
     * with respect to nodal velocities.
     */
    virtual void computeDiffusionDerivativeTerm_MB(FloatMatrix &answer, MatResponseMode mode, TimeStep *atTime) = 0;
    /**
     * Computes pressure terms for momentum balance equations(s).
     */
    virtual void computePressureTerm_MB(FloatMatrix &answer, TimeStep *atTime) = 0;
    /**
     * Computes Lhs terms due to boundary conditions - pressure.
     */
    virtual void computeBCLhsPressureTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    /**
     * Computes SLIC stabilization term for momentum balance equation(s).
     */
    virtual void computeLSICStabilizationTerm_MB(FloatMatrix &answer, TimeStep *atTime) = 0;
    /**
     * Computes the linear advection term for mass conservation equation.
     */
    virtual void computeLinearAdvectionTerm_MC(FloatMatrix &answer, TimeStep *atTime) = 0;
    /**
     * Computes advection terms for mass conservation equation.
     */
    virtual void computeAdvectionTerm_MC(FloatArray &answer, TimeStep *atTime) = 0;
    /**
     * Computes the derivative of advection terms for mass conservation equation
     * with respect to nodal velocities.
     */
    virtual void computeAdvectionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime) = 0;
    /**
     * Computes diffusion derivative terms for mass conservation equation.
     */
    virtual void computeDiffusionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime) = 0;
    /**
     * Computes diffusion terms for mass conservation equation.
     */
    virtual void computeDiffusionTerm_MC(FloatArray &answer, TimeStep *atTime) = 0;
    /**
     * Computes acceleration terms for mass conservation equation.
     */
    virtual void computeAccelerationTerm_MC(FloatMatrix &answer, TimeStep *atTime) = 0;
    /**
     * Computes pressure terms for mass conservation equation.
     */
    virtual void computePressureTerm_MC(FloatMatrix &answer, TimeStep *atTime) = 0;
    /**
     * Computes Lhs terms due to boundary conditions - velocity.
     */
    virtual void computeBCLhsTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    /**
     * Computes Lhs terms due to boundary conditions - pressure.
     */
    virtual void computeBCLhsPressureTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    /**
     * Computes Rhs terms due to boundary conditions.
     */
    virtual void computeBCRhsTerm_MB(FloatArray &answer, TimeStep *atTime) = 0;
    /**
     * Computes Rhs terms due to boundary conditions.
     */
    virtual void computeBCRhsTerm_MC(FloatArray &answer, TimeStep *atTime) = 0;
    /**
     * Computes Lhs term due to applied slip with friction bc.
     */
    virtual void computeSlipWithFrictionBCTerm_MB(FloatMatrix &answer, Load *load, int side, TimeStep *atTime) {
      _warning("computeSlipWithFrictionBCTerm_MB not implemented");
      answer.resize(0,0);
    }
    /**
     * Computes Lhs contribution due to applied Penetration bc.
     */
    virtual void computePenetrationWithResistanceBCTerm_MB(FloatMatrix &answer, Load *load, int side, TimeStep *atTime){
      _warning("computePenetrationWithResistanceBCTerm_MB not implemented");
      answer.resize(0,0);
    }
    /**
     * Computes Lhs contribution due to outflow BC.
     */
    virtual void computeOutFlowBCTerm_MB(FloatMatrix &answer, int side, TimeStep *atTime){
      _warning("computeOutFlowBCTerm_MB not implemented");
      answer.resize(0,0);
    }
    

    virtual void computeHomogenizedReinforceTerm_MB(FloatMatrix &answer,  Load * load, TimeStep *atTime){
      _warning("computeHomogenizedReinforceTerm_MB");
      answer.resize(0,0); 
    }
    virtual void computeHomogenizedReinforceTerm_MC(FloatMatrix &answer,  Load * load, TimeStep *atTime){
      _warning("computeHomogenizedReinforceTerm_MB");
      answer.resize(0,0); 
    }


    /// Computes the critical time increment.
    virtual double computeCriticalTimeStep(TimeStep *tStep) = 0;

    // time step termination
    virtual void updateInternalState(TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual int checkConsistency();

    // definition
    virtual const char *giveClassName() const { return "SUPGElement"; }
    virtual classType giveClassID() const { return SUPGElementClass; }

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *gp);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type);

#ifdef __OOFEG
    int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                int node, TimeStep *atTime);
    // Graphics output
    //void drawYourself(oofegGraphicContext&);
    //virtual void drawRawGeometry(oofegGraphicContext&) {}
    //virtual void drawDeformedGeometry(oofegGraphicContext&, UnknownType) {}
#endif

protected:
    virtual void giveLocalVelocityDofMap(IntArray &map) {}
    virtual void giveLocalPressureDofMap(IntArray &map) {}

    virtual void computeDeviatoricStrain(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) = 0;
    virtual void computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);
};
} // end namespace oofem
#endif // supgelement_h
