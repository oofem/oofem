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

//   **************************************************************************************
//   *** CLASS GENERAL ELEMENT CLASS FOR FLUID DYNAMIC PROBLEMS SOLVED WITH PFEM METHOD ***
//   **************************************************************************************

#ifndef pfemelement_h
#define pfemelement_h


#include "fmelement.h"
#include "femcmpnn.h"
#include "domain.h"
#include "floatmatrix.h"
#include "material.h"

#include "primaryfield.h"

namespace oofem {
class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * This abstract class represent a general base element class for
 * fluid dynamic problems solved using PFEM algorithm.
 */
class PFEMElement : public FMElement
{
public:
protected:

public:
    /// Constructor
    PFEMElement(int, Domain *);
    // Destructor
    ~PFEMElement();

    ///Initializes receiver acording to object description stored in input record.
    IRResultType initializeFrom(InputRecord *ir);

    // characteristic  matrix
    void giveCharacteristicMatrix(FloatMatrix &answer, CharType, TimeStep *);
    void giveCharacteristicVector(FloatArray &answer, CharType, ValueModeType, TimeStep *);
    virtual double giveCharacteristicValue(CharType, TimeStep *);

    /** Calculates consistent mass matrix */
    virtual void          computeConsistentMassMtrx(FloatMatrix &answer, TimeStep *) = 0;
    /** Calculates diagonal mass matrix in form of an FloatArray*/
    virtual void          computeDiagonalMassMtrx(FloatArray &answer, TimeStep *) = 0;
    /** Calculates diagonal mass matrix */
    virtual void          computeDiagonalMassMtrx(FloatMatrix &answer, TimeStep *) = 0;
    /** Calculates inverted stabilization diagonal mass matrix*/
    virtual void computeInvertedStabilizationDiagonalMassMtrx(FloatMatrix &answer, TimeStep *atTime) = 0;

    /// Calculates critical time step
    virtual double        computeCriticalTimeStep(TimeStep *tStep) = 0;
    // time step termination
    /**
     * Updates element state corresponding to newly reached solution.
     * It computes stress vector in each element integration point (to ensure that data in integration point's
     * statuses are valid).
     * @param tStep finished time step
     */
    void                  updateInternalState(TimeStep *);
    void                  printOutputAt(FILE *, TimeStep *);
    virtual int           checkConsistency();

    /// returns the interpolation for velocity
    virtual FEInterpolation *giveVelocityInterpolation() = 0;
    /// returns the interpolation for the pressure
    virtual FEInterpolation *givePressureInterpolation() = 0;

    // definition
    const char *giveClassName() const { return "PFEMElement"; }
    classType                giveClassID() const { return PFEMElementClass; }

    /**
     * Returns the mask of reduced indexes of Internal Variable component .
     * @param answer mask of Full VectorSize, with components beeing the indexes to reduced form vectors.
     * @param type determines the internal variable requested (physical meaning)
     * @returns nonzero if ok or error is generated for unknown mat mode.
     */
    //virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type);

#ifdef __OOFEG
    int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                int node, TimeStep *atTime);
    //
    // Graphics output
    //
    //void          drawYourself (oofegGraphicContext&);
    //virtual void  drawRawGeometry (oofegGraphicContext&) {}
    //virtual void  drawDeformedGeometry(oofegGraphicContext&, UnknownType) {}
#endif

protected:
    /// Computes deviatoric stress vector in given integration point and solution step from given total strain vector
    virtual void computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *) = 0;

    /// Calculates the element shape funciton matrix for velocity degrees of freedom
    virtual void computeNuMatrix(FloatMatrix &answer, GaussPoint *gp) = 0;
    /// Calculates the element shape funciton matrix for pressure degrees of freedom
    virtual void computeNpMatrix(FloatMatrix &answer, GaussPoint *gp) = 0;
    /// Calculates the element shape function derivative matrix
    virtual void computeBMatrix(FloatMatrix &answer, GaussPoint *gp) = 0;
    /// Calculates the stiffness matrix
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *atTime) = 0; //K
    /// Calculates the velocity laplacian matrix
    virtual void computeVelocityLaplacianMatrix(FloatMatrix &answer, TimeStep *atTime) = 0;
    /// Calculates the pressure laplacian matrix
    virtual void computePressureLaplacianMatrix(FloatMatrix &answer, TimeStep *atTime) = 0;
    /// Calculates the stabilized pressure laplacian matrix
    virtual void computeStabilizedLaplacianMatrix(FloatMatrix &answer, TimeStep *atTime) = 0; //L_tau
    /// Calculates the pressure gradient matrix
    virtual void computeGradientMatrix(FloatMatrix &answer, TimeStep *atTime) = 0; //G
    /// Calculates the velocity divergence matrix
    virtual void computeDivergenceMatrix(FloatMatrix &answer, TimeStep *atTime) = 0; //D
    /// Calculates the stabilized gradient matrix
    virtual void computeStabilizationGradientMatrix(FloatMatrix &answer, TimeStep *atTime) = 0; //Q
    /// Calculates the stabilized mass matrix
    virtual void computeStabilizationMassMatrix(FloatMatrix &answer, TimeStep *atTime) = 0; //M_hat
    /// Calculates the force vector
    virtual void computeForceVector(FloatArray &answer, TimeStep *atTime) = 0; //F
    /// Calculates the so called substitution matrix
    virtual void computePFEMSubstitutionMatrix(FloatMatrix &answer,  TimeStep *atTime) = 0; //S
    /// Calculates the prescribed velocity vector for the right hand side of the pressure equation
    virtual void computePrescribedRhsVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode) = 0;
	/**
	 * Calculates the prescribed pressure vector for the right hand side of the velocity equation
	 * Implementation not finished yet. At the moment not needed due to improved boundary condition setting in alpha shape
	 */ 
	virtual void computePrescribedPressureRhsVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode) = 0;
};
} // end namespace oofem
#endif // pfemelement_h







