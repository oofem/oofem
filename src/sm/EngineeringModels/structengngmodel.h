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

#ifndef structengngmodel_h
#define structengngmodel_h

#include "engngm.h"
#include "statecountertype.h"
#include "floatarray.h"

namespace oofem {
class StructuralElement;

/// Assembles the internal forces, without updating the strain.
///@todo The need for this is just due to some other design choices. 
class LastEquilibratedInternalForceAssembler : public InternalForceAssembler
{
    virtual void vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const;
};

/**
 * Callback class for assembling linearized thermal "loads", useful for computing initial guesses.
 * @author Mikael Öhman
 */
class LinearizedDilationForceAssembler : public VectorAssembler
{
public:
    virtual void vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const;
};

/**
 * Callback class for assembling initial stress matrices
 * @author Mikael Öhman
 */
class InitialStressMatrixAssembler : public MatrixAssembler
{
public:
    virtual void matrixFromElement(FloatMatrix &mat, Element &element, TimeStep *tStep) const;
};


/**
 * This class implements extension of EngngModel for structural models.
 * Its purpose is to declare and implement general methods for computing reaction forces.
 */
class StructuralEngngModel : public EngngModel
{
protected:
    /**
     * Contains last time stamp of internal variable update.
     * This update is made via various services
     * (like those for computing real internal forces or updating the internal state).
     */
    StateCounterType internalVarUpdateStamp;

    /// Norm of nodal internal forces evaluated on element by element basis (squared)
    FloatArray internalForcesEBENorm;
    /**
     * Computes and prints reaction forces, computed from nodal internal forces. Assumes, that real
     * stresses corresponding to reached state are already computed (uses giveInternalForcesVector
     * structural element service with useUpdatedGpRecord = 1 parameter). Only the dof managers selected for
     * output (OutputManager) are handled.
     * @see StructuralElement::giveInternalForcesVector
     * @see OutputManager
     * @param tStep Time step.
     * @param id Domain number.
     */
    void printReactionForces(TimeStep *tStep, int id, FILE *out);

    /**
     * Computes the contribution external loading to reaction forces in given domain. Default implementations adds the
     * contribution from computeElementLoadReactionContribution and computeElementLoadReactionContribution methods.
     * @param reactions Contains the computed contributions.
     * @param tStep Solution step.
     * @param di Domain number.
     */
    virtual void computeExternalLoadReactionContribution(FloatArray &reactions, TimeStep *tStep, int di);
    /**
     * Evaluates the nodal representation of internal forces by assembling contributions from individual elements.
     * @param answer Vector of nodal internal forces.
     * @param normFlag True if element by element norm of internal forces (internalForcesEBENorm) is to be computed.
     * @param di Domain number.
     * @param tStep Solution step.
     */
    virtual void giveInternalForces(FloatArray &answer, bool normFlag, int di, TimeStep *tStep);

    /**
     * Updates nodal values
     * (calls also this->updateDofUnknownsDictionary for updating dofs unknowns dictionaries
     * if model supports changes of static system). The element internal state update is also forced using
     * updateInternalState service.
     * @param tStep Solution step.
     */
    void updateInternalState(TimeStep *tStep);

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

public:
    /// Creates new StructuralEngngModel with number i, associated to domain d.
    StructuralEngngModel(int i, EngngModel * _master = NULL);
    /// Destructor.
    virtual ~StructuralEngngModel();

    virtual void updateYourself(TimeStep *tStep);

    virtual int checkConsistency();

    /**
     * Computes reaction forces. The implementation assumes, that real
     * stresses corresponding to reached state are already computed (uses giveInternalForcesVector
     * structural element service with useUpdatedGpRecord = 1 parameter).
     * To be safe, this method should be called after convergence has been reached, eq.,
     * after engngModel->updateYourself() has been called.
     * @param answer Reactions, the ordering of individual values follows numbering of prescribed equations.
     * @param tStep Time step.
     * @param di Domain number.
     */
    void computeReaction(FloatArray &answer, TimeStep *tStep, int di);

    /**
     * Terminates the solution of time step. Default implementation calls prinOutput() service and if specified,
     * context of whole domain is stored and output for given time step is printed.
     */
    virtual void terminate(TimeStep *tStep);
    
    /**
     * Builds the reaction force table. For each prescribed equation number it will find
     * corresponding node and dof number. The entries in the restrDofMans, restrDofs, and eqn
     * arrays are sorted with increasing dofman number and with increasing dof number as
     * a second minor criterion.
     * @param restrDofMans Contains numbers of restrained Dofmanagers, with size equal to total number of prescribed equations.
     * @param restrDofs Contains IDs of restrained Dofs, with size equal to total number of prescribed equations.
     * @param eqn Contains the corresponding restrained equation numbers.
     * @param tStep Time step.
     * @param di Domain number.
     */
    void buildReactionTable(IntArray &restrDofMans, IntArray &restrDofs, IntArray &eqn, TimeStep *tStep, int di);

#ifdef __OOFEG
    /**
     * Shows the sparse structure of required matrix, type == 1 stiffness.
     */
    virtual void showSparseMtrxStructure(int type, oofegGraphicContext &gc, TimeStep *tStep);
#endif
};
} // end namespace oofem
#endif // structengngmodel_h
