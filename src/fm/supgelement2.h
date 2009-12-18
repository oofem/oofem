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

//   *******************************************************************
//   *** CLASS GENERAL STABILIZED SUPG/PSPG ELEMENT for CFD Analysis ***
//   *******************************************************************

#ifndef supgelement2_h
#define supgelement2_h


#include "supgelement.h"
#include "femcmpnn.h"
#include "domain.h"
#include "flotmtrx.h"

#include "primaryfield.h"
#include "fluiddynamicmaterial.h"

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
 * fluid dynamic problems. It is derived from SUPGElement base,
 * but it provides general algorithms how to obtain characteristic
 * vectors&matrices by means of numerical integration.
 * It is therefore suitable for high-order elements where exact
 * integration would be difficult.
 */
class SUPGElement2 : public SUPGElement
{
public:
protected:
    /**
     * stabilization coefficients, updated for each solution step
     * in updateStabilizationCoeffs()
     */
    //double t_supg, t_pspg, t_lsic;

public:
    // constructor
    SUPGElement2(int, Domain *);
    ~SUPGElement2();                        // destructor

    ///Initializes receiver acording to object description stored in input record.
    IRResultType initializeFrom(InputRecord *ir);

    // characteristic  matrix
    void giveCharacteristicMatrix(FloatMatrix & answer, CharType, TimeStep *);
    void giveCharacteristicVector(FloatArray & answer, CharType, ValueModeType, TimeStep *);
    virtual double giveCharacteristicValue(CharType, TimeStep *);
    virtual void     updateStabilizationCoeffs(TimeStep *);
    virtual void     updateElementForNewInterfacePosition(TimeStep *) { }

    /**
     * Computes acceleration terms (generalized mass matrix with stabilization terms ) for momentum balance equations(s)
     */
    virtual void computeAccelerationTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    /**
     * Computes nonlinear advection terms for momentum balance equations(s)
     */
    virtual void computeAdvectionTerm_MB(FloatArray &answer, TimeStep *atTime);
    /**
     * Computes the derivative of advection terms for momentum balance equations(s)
     * with respect to nodal velocities
     */
    virtual void computeAdvectionDerivativeTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    /**
     *  Computes diffusion terms for momentum balance equations(s)
     */
    virtual void computeDiffusionTerm_MB(FloatArray &answer, TimeStep *atTime);
    /** Computes the derivative of diffusion terms for momentum balance equations(s)
     *  with respect to nodal velocities
     */
    virtual void computeDiffusionDerivativeTerm_MB(FloatMatrix &answer, MatResponseMode mode, TimeStep *atTime);
    /** Computes pressure terms for momentum balance equations(s) */
    virtual void computePressureTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    /** Computes SLIC stabilization term for momentum balance equation(s) */
    virtual void computeLSICStabilizationTerm_MB(FloatMatrix &answer, TimeStep *atTime);
    /** Computes the linear advection term for mass conservation equation
     */
    virtual void computeLinearAdvectionTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    /**
     * Computes advection terms for mass conservation equation
     */
    virtual void computeAdvectionTerm_MC(FloatArray &answer, TimeStep *atTime);
    /** Computes the derivative of advection terms for mass conservation equation
     *  with respect to nodal velocities
     */
    virtual void computeAdvectionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    /**
     * Computes diffusion terms for mass conservation equation
     */
    virtual void computeDiffusionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    virtual void computeDiffusionTerm_MC(FloatArray &answer, TimeStep *atTime);
    /**
     * Computes acceleration terms for mass conservation equation
     */
    virtual void  computeAccelerationTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    /**
     * Computes pressure terms for mass conservation equation
     */
    virtual void computePressureTerm_MC(FloatMatrix &answer, TimeStep *atTime);
    /**
     * Computes Rhs terms due to boundary conditions
     */
    virtual void  computeBCRhsTerm_MB(FloatArray &answer, TimeStep *atTime);
    /**
     * Computes Rhs terms due to boundary conditions
     */
    virtual void  computeBCRhsTerm_MC(FloatArray &answer, TimeStep *atTime);


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
    const char *giveClassName() const { return "SUPGElement2"; }
    classType                giveClassID() const { return SUPGElementClass; }

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
    virtual void computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *);
    virtual void computeNuMatrix(FloatMatrix &answer, GaussPoint *gp) = 0;
    virtual void computeUDotGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *atTime) = 0;
    virtual void computeBMatrix(FloatMatrix &anwer, GaussPoint *gp) = 0;
    virtual void computeDivUMatrix(FloatMatrix &answer, GaussPoint *gp) = 0;
    virtual void computeNpMatrix(FloatMatrix &answer, GaussPoint *gp) = 0;
    virtual void computeGradPMatrix(FloatMatrix &answer, GaussPoint *gp) = 0;
    virtual int  giveNumberOfSpatialDimensions() = 0;



    //virtual void computeEdgeNuMatrix (FloatMatrix& answer, GaussPoint* gp) = 0;
    //virtual double computeEdgeVolumeAround (GaussPoint*, int iedge) = 0;
    //virtual void computeEdgeIpGlobalCoords (FloatArray& answer, GaussPoint* gp, int iedge) = 0;
    //virtual int computeLoadGToLRotationMtrx (FloatMatrix& answer) = 0;
    //virtual int computeLoadLEToLRotationMatrix (FloatMatrix& answer, int iedge, GaussPoint* gp) = 0;
    //virtual void giveEdgeUDofMapping (IntArray& answer, int iedge) = 0;
    //virtual int  giveUApproxOrder () = 0;
    virtual void computeEdgeLoadVector_MB(FloatArray &answer, Load *load, int id, TimeStep *stepN);
    virtual void computeSurfaceLoadVector_MB(FloatArray &answer, Load *load, int id, TimeStep *stepN);
    virtual void computeEdgeLoadVector_MC(FloatArray &answer, Load *load, int id, TimeStep *stepN);
    virtual void computeSurfaceLoadVector_MC(FloatArray &answer, Load *load, int id, TimeStep *stepN);
};

} // end namespace oofem
#endif // supgelement2_h
