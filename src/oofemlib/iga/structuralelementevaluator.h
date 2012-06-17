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

#ifndef structuralelementevaluator_h
#define structuralelementevaluator_h

#include "iga.h"
#include "matresponsemode.h"

namespace oofem {
/**
 * This class represent a new concept on how to define elements.
 * Traditionally, Elements are derived from problem specific class (structural element, for example)
 * define their interpolation and implement methods to evaluate shape function matrix, geometrical matrix, etc.
 * This new concept, here represented by  StructuralElementEvaluator and derived classes allows to
 * define all problem specific methods (including shape function matrix, geometrical matrix evaluation) to be defined
 * once for all elements of the same type (plane stress elements, space3d elements, etc) just relying on
 * services of finite element interpolation classes.
 * Definition of particular element is then done simply by deriving if from evaluator and providing interpolation.
 *
 * StructuralElementEvaluator - base class of all structural elements
 * Individual elements supposed to be derived from StructuralElementEvaluator and IGAElement
 * @todo{Class is missing much documentation.}
 */
class StructuralElementEvaluator
{
protected:
    FloatMatrix rotationMatrix;
    /// Flag indicating if transformation matrix has been already computed
    int rotationMatrixDefined;

    StructuralElementEvaluator();
    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep);
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep);

    /**
     * Returns the integration rule for mass matrices, if relevant.
     * @return Number of integration points for mass matrix.
     */
    virtual IntegrationRule* giveMassMtrxIntegrationRule() { return NULL; }
    /**
     * Returns mask indicating, which unknowns (their type and ordering is the same as
     * element unknown vector) participate in mass matrix integration.
     * Nonzero value at i-th position
     * indicates that corresponding row in interpolation matrix N will participate in
     * mass matrix integration (typically only displacements are taken into account).
     * @param answer Integration mask, if zero sized, all unknowns participate. This is default.
     */
    virtual void giveMassMtrxIntegrationMask(IntArray &answer) { answer.resize(0); }
    /**
     * Computes lumped mass matrix of receiver. Default implementation returns lumped consistent mass matrix.
     * Then returns lumped mass transformed into nodal coordinate system.
     * The lumping procedure zeroes all off-diagonal members and zeroes also all diagonal members
     * corresponding to non-displacement DOFs. Such diagonal matrix is then rescaled, to preserve
     * the element mass.
     * Requires the computeNmatrixAt and giveMassMtrxIntegrationgMask services to be implemented.
     * @param answer Lumped mass matrix.
     * @param tStep Time step.
     */
    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    /**
     * Computes consistent mass matrix of receiver using numerical integration over element volume.
     * Mass matrix is computed as @f$ M = \int_V N^{\mathrm{T}} \rho N dV @f$, where @f$ N @f$ is displacement approximation matrix.
     * The number of necessary integration points  is determined using this->giveNumberOfIPForMassMtrxIntegration
     * service. Only selected degrees of freedom participate in integration of mass matrix. This is described
     * using dof mass integration mask. This mask is obtained from this->giveMassMtrxIntegrationgMask service.
     * The nonzero mask value at i-th position indicates that i-th element DOF participates in mass matrix
     * computation. The result is in element local coordinate system.
     * @param answer Mass matrix.
     * @param tStep Time step.
     * @param mass Total mass of receiver.
     */
    virtual void computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass);

protected:
    virtual Element *giveElement() = 0;

    /**
     * Computes the matrix for which the unknown field is obtained, typically [N1, 0, N2, 0, ...; 0, N1, 0, N2, ...].
     */
    virtual void computeNMatrixAt(FloatMatrix &answer, GaussPoint *gp) = 0;
    virtual void computeBMatrixAt(FloatMatrix &answer, GaussPoint *gp) = 0;
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual double computeVolumeAround(GaussPoint *gp) { return 0.; }
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, bool useUpdatedGpRecord = false);
    void computeVectorOf(EquationID type, ValueModeType u,
                         TimeStep *stepN, FloatArray &answer) {
        this->giveElement()->computeVectorOf(type, u, stepN, answer);
    }
    void computeVectorOf(PrimaryField &field, ValueModeType u, TimeStep *stepN, FloatArray &answer) {
        this->giveElement()->computeVectorOf(field, u, stepN, answer);
    }
    void computeVectorOfPrescribed(EquationID ut, ValueModeType type, TimeStep *stepN, FloatArray &answer) {
        this->giveElement()->computeVectorOfPrescribed(ut, type, stepN, answer);
    }
    bool isActivated(TimeStep *atTime) { return true; }
    void updateInternalState(TimeStep *stepN);
    void computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);
    void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);
    /*
     * Optimized version, allowing to pass element displacements as parameter.
     * Standard version has a huge performance leak; in typical IGA element the element vector is VERY large
     * and its querying for each point take more time than strain evaluation. And this has to be done for each
     * integration point. This optimized version allows to assemble displacement vector only once (for all IP)
     * and pass this vector as parameter
     */
    void computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN, FloatArray &u);
    void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN, FloatArray &u);

    /*
     * Assembles the code numbers of given integration element (sub-patch)
     * This is done by obtaining list of nonzero shape functions and
     * by collecting the code numbers of nodes corresponding to these
     * shape functions
     * @returns returns nonzero if integration rule code numbers differ from element code numbers
     */
    //virtual int giveIntegrationElementCodeNumbers(IntArray &answer, Element *elem,
    //                                              IntegrationRule *ie, EquationID ut);

    /**
     * Assembles the local element code numbers of given integration element (sub-patch)
     * This is done by obtaining list of nonzero shape functions and
     * by collecting the code numbers of nodes corresponding to these
     * shape functions
     * @return Nonzero if integration rule code numbers differ from element code numbers
     */
    virtual int giveIntegrationElementLocalCodeNumbers(IntArray &answer, Element *elem,
                                                       IntegrationRule *ie, EquationID ut);
#ifdef __OOFEG
    friend void drawIGAPatchDeformedGeometry(Element * elem, StructuralElementEvaluator * se, oofegGraphicContext & gc, UnknownType);
#endif
    public:
    void elem(int arg1, EquationID arg2, IntArray arg3);
};
} // end namespace oofem
#endif //structuralelementevaluator_h
