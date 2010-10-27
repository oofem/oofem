/* $Header: /home/cvs/bp/oofem/oofemlib/src/element.h,v 1.27 2003/04/06 14:08:24 bp Exp $ */
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

#ifndef structuralelementevaluator_h
#define structuralelementevaluator_h

#include "iga.h"
#include "matresponsemode.h"

namespace oofem {

/**
 * 
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
 *
 */
class StructuralElementEvaluator
{
protected:
    FloatMatrix *rotationMatrix;
    /// Flag indicating if tranformation matrix has been already computed
    int rotationMatrixDefined;

    StructuralElementEvaluator();
    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep);
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep) {
        if ( type == ElementNonForceLoadVector ) {
            this->computeNonForceLoadVector(answer, tStep, mode);
        } else {
            answer.resize(0);
        }
    }

protected:
    virtual Element *giveElement() = 0;
    virtual void computeNMatrixAt(FloatMatrix &answer, GaussPoint *gp) = 0;
    virtual void computeBMatrixAt(FloatMatrix &answer, GaussPoint *gp) = 0;
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual double computeVolumeAround(GaussPoint *gp) { return 0.; }
    void  computeNonForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode);
    void  computeBcLoadVectorAt(FloatArray &answer, TimeStep *, ValueModeType mode);
    virtual void giveInternalForcesVector(FloatArray &answer,
                                          TimeStep *, int useUpdatedGpRecord = 0) {
        answer.resize(0);
    }
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
    bool   isActivated(TimeStep *atTime) { return true; }
    void   updateInternalState(TimeStep *stepN);
    void   computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);
    void   computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);
    /* Optimized version, allowing to pass element displacents as parameter.
     * Standart version has a huge performance leak; in typical iga element the element vector is VERY large
     * and its querying for each point take more time than strain evaluation. And this has to be done for each
     * integration point. This optimized version allows to assemble displacement vector only once (for all IP)
     * and pass this vector as parameter
     */
    void   computeStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN, FloatArray &u);
    void   computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *stepN, FloatArray &u);

    /**
     * Updates rotation matrix r(l)=T r(g*) between  local and global coordinate system
     * taking into account also possible local - coordinate system in some elements
     * nodes.
     * Default implementation uses \ref computeGtoLRotationMatrix and
     \ref computeGNDofRotationMatrix  services to compute result.
     * Default implementation uses cached rotation matrix in
     * rotationMatrix attribute, so rotation matrix is computed only once.
     * @return nonzero if transformation is necessary.
     */
    virtual int updateRotationMatrix() ;
    /**
     * Returns transformation matrix for DOFs from global coordinate system
     * to local coordinate system in nodes (i.e. r(n)=T r(g)) if mode == _toNodalCS.
     * If mode == _toGlobalCS, the transformation from local nodal cs to
     * global cs in node is returned. If no trasformation is
     * necessary sets answer to empty mtrx and returns zero value.
     * @return nonzero if transformation is necessary, zero otherwise.
     */
    virtual int  computeGNDofRotationMatrix(FloatMatrix &answer, DofManTransfType mode);
    /**
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
     * @returns returns nonzero if integration rule code numbers differ from element code numbers
     */
    virtual int giveIntegrationElementLocalCodeNumbers(IntArray &answer, Element *elem,
                                                       IntegrationRule *ie, EquationID ut);
#ifdef __OOFEG
    friend void drawIGAPatchDeformedGeometry(Element * elem, StructuralElementEvaluator * se, oofegGraphicContext & gc, UnknownType);
#endif
};

} // end namespace oofem
#endif //structuralelementevaluator_h
