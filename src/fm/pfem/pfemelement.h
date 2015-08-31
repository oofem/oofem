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

//   ********************************************************************************
//   *** GENERAL ELEMENT CLASS FOR FLUID DYNAMIC PROBLEMS SOLVED WITH PFEM METHOD ***
//   ********************************************************************************

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
 *
 * @author David Krybus
 */
class PFEMElement : public FMElement
{
public:
    /// Constructor
    PFEMElement(int, Domain *);
    /// Destructor
    ~PFEMElement();

    ///Initializes receiver acording to object description stored in input record.
    IRResultType initializeFrom(InputRecord *ir);

    // characteristic  matrix
    void giveCharacteristicMatrix(FloatMatrix & answer, CharType, TimeStep *);
    void giveCharacteristicVector(FloatArray & answer, CharType, ValueModeType, TimeStep *);

    /** Calculates diagonal mass matrix as vector*/
    virtual void computeDiagonalMassMtrx(FloatArray &answer, TimeStep *) = 0;
    /** Calculates diagonal mass matrix */
    virtual void computeDiagonalMassMtrx(FloatMatrix &answer, TimeStep *) = 0;
    /// Calculates critical time step
    virtual double computeCriticalTimeStep(TimeStep *tStep) = 0;

    virtual void updateInternalState(TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual int checkConsistency();

    /// Returns the interpolation for velocity
    virtual FEInterpolation *giveVelocityInterpolation() = 0;
    /// Returns the interpolation for the pressure
    virtual FEInterpolation *givePressureInterpolation() = 0;

    // definition
    virtual const char *giveClassName() const { return "PFEMElement"; }

    /// Returns mask of velocity Dofs
    virtual const IntArray &giveVelocityDofMask() const = 0;
    /// Returns mask of pressure Dofs
    virtual const IntArray &givePressureDofMask() const = 0;

#ifdef __OOFEG
    virtual int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                        int node, TimeStep *atTime);
    //
    // Graphics output
    //
    //virtual void drawYourself(oofegGraphicContext&);
    //virtual void drawRawGeometry(oofegGraphicContext&) {}
    //virtual void drawDeformedGeometry(oofegGraphicContext&, UnknownType) {}
#endif

    /// Computes deviatoric stress vector in given integration point and solution step from given total strain vector
    virtual void computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *) = 0;
    /// Calculates the divergence of the deviatoric stress
    virtual void computeDeviatoricStressDivergence(FloatArray &answer, TimeStep *atTime) = 0;
    /// Calculates the element shape function derivative matrix
    virtual void computeBMatrix(FloatMatrix &answer, GaussPoint *gp) = 0;
    /// Calculates the stiffness matrix
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode mode, TimeStep *atTime) = 0; //K
    /// Calculates the pressure laplacian matrix
    virtual void computePressureLaplacianMatrix(FloatMatrix &answer, TimeStep *atTime) = 0;
    /// Calculates the pressure gradient matrix
    virtual void computeGradientMatrix(FloatMatrix &answer, TimeStep *atTime) = 0; //G
    /// Calculates the velocity divergence matrix
    virtual void computeDivergenceMatrix(FloatMatrix &answer, TimeStep *atTime) = 0; //D
    /// Calculates the force vector
    virtual void computeForceVector(FloatArray &answer, TimeStep *atTime) = 0; //F
    /// Calculates the prescribed velocity vector for the right hand side of the pressure equation
    virtual void computePrescribedRhsVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode) = 0;
};
} // end namespace oofem
#endif // pfemelement_h
