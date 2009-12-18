/* $Header: /home/cvs/bp/oofem/tm/src/transportelement.h,v 1.3 2003/04/23 14:22:15 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   ****************************************************************************************
//   *** CLASS GENERAL ELEMENT CLASS FOR FLUID DYNAMIC PROBLEMS SOLVED WITH CBS ALGORITHM ***
//   ****************************************************************************************

#ifndef cbselement_h
#define cbselement_h


#include "fmelement.h"
#include "femcmpnn.h"
#include "domain.h"
#include "flotmtrx.h"

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
 * fluid dynamic problems solved using CBS algorithm.
 */
class CBSElement : public FMElement
{
public:
protected:

public:
    // constructor
    CBSElement(int, Domain *);
    ~CBSElement();                        // destructor

    ///Initializes receiver acording to object description stored in input record.
    IRResultType initializeFrom(InputRecord *ir);

    // characteristic  matrix
    void giveCharacteristicMatrix(FloatMatrix & answer, CharType, TimeStep *);
    void giveCharacteristicVector(FloatArray & answer, CharType, ValueModeType, TimeStep *);
    virtual double giveCharacteristicValue(CharType, TimeStep *);
    /**
     * Returns local vector of prescribed unknowns. Local vector of prescribed unknowns is
     * extracted from nodal (and side - if they hold unknowns) boundary conditions.
     * @param u    Identifies mode of unknown (eg. total values or velocity of unknown).
     * @param stepN Time step, when vector of prescribed unknowns is requested.
     * @param answer Local vector of prescribed unknowns. If unknown is not prescibed,
     * zero value is placed on position of free dof.
     */
    ///void                  computeVectorOfPrescribed (EquationID ut, ValueModeType type, TimeStep* stepN, FloatArray& answer) ;
    /** Calculates consistent mass matrix */
    virtual void          computeConsistentMassMtrx(FloatMatrix &answer, TimeStep *) = 0;
    /** Calculates diagonal mass matrix */
    virtual void          computeDiagonalMassMtrx(FloatArray &answer, TimeStep *) = 0;
    /** Calculates rhs due to prescribed (*) velocities for (*) velocities */
    virtual void computePrescribedTermsI(FloatArray & answer, ValueModeType, TimeStep *);
    /** Calculates rhs due to prescribed density/pressure */
    //virtual void          computePrescribedTermsII (FloatArray& answer, ValueModeType, TimeStep*);
    /** Calculates convection component for (*) velocities */
    virtual void          computeConvectionTermsI(FloatArray &answer, TimeStep *) = 0;
    /** Calculates contribution from diffusion terms for (*) velocities */
    virtual void          computeDiffusionTermsI(FloatArray &answer, TimeStep *) = 0;
    /// computes velocity terms on RHS for density equation
    virtual void          computeDensityRhsVelocityTerms(FloatArray &answer, TimeStep *tStep) = 0;
    /// computes pressure terms on RHS for density equation
    virtual void          computeDensityRhsPressureTerms(FloatArray &answer, TimeStep *tStep) = 0;
    /// computes prescribed pressure due to applied tractions
    virtual void          computePrescribedTractionPressure(FloatArray &answer, TimeStep *tStep) = 0;
    /// computes number of edges/sides with prescribed traction contributing to node with prescribed pressure
    virtual void          computeNumberOfNodalPrescribedTractionPressureContributions(FloatArray &answer, TimeStep *tStep) = 0;
    /// calculates the pressure LHS
    virtual void          computePressureLhs(FloatMatrix &answer, TimeStep *tStep) = 0;
    /// calculates the RHS of velocity correction step
    virtual void          computeCorrectionRhs(FloatArray &answer, TimeStep *tStep) = 0;
    /// calculates critical time step
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

    // definition
    const char *giveClassName() const { return "CBSElement"; }
    classType                giveClassID() const { return CBSElementClass; }

    /**
     * Returns the mask of reduced indexes of Internal Variable component .
     * @param answer mask of Full VectorSize, with components beeing the indexes to reduced form vectors.
     * @param type determines the internal variable requested (physical meaning)
     * @returns nonzero if ok or error is generated for unknown mat mode.
     */
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type);

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
};

} // end namespace oofem
#endif // cbselement_h
