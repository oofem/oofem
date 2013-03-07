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

#ifndef structengngmodel_h
#define structengngmodel_h

#include "engngm.h"

#ifdef __PARALLEL_MODE
 #include "problemcomm.h"
 #include "processcomm.h"
#endif

namespace oofem {
class StructuralElement;

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

#ifdef __PARALLEL_MODE
    /// Common Communicator buffer.
    CommunicatorBuff *commBuff;
    /// Communicator.
    ProblemCommunicator *communicator;

    /// Flag indicating if nonlocal extension active.
    int nonlocalExt;
    /// NonLocal Communicator. Necessary when nonlocal constitutive models are used.
    ProblemCommunicator *nonlocCommunicator;
    /// Message tags
    enum { InternalForcesExchangeTag, MassExchangeTag, LoadExchangeTag, ReactionExchangeTag, RemoteElementExchangeTag };
    /**
     * Packing function for vector values of DofManagers. Packs vector values of shared/remote DofManagers
     * into send communication buffer of given process communicator.
     * @param processComm Task communicator.
     * @param src Source vector.
     * @param prescribedEquations True if src uses prescribed equations numbering.
     * @return Nonzero if successful.
     */
    int packDofManagers(FloatArray *src, ProcessCommunicator &processComm, bool prescribedEquations);
    /**
     * Unpacking function for vector values of DofManagers . Unpacks vector of shared/remote DofManagers
     * from  receive communication buffer of given process communicator.
     * @param processComm Task communicator.
     * @param dest Destination vector.
     * @param prescribedEquations True if src uses prescribed equations numbering.
     * @return Nonzero if successful.
     */
    int unpackDofManagers(FloatArray *dest, ProcessCommunicator &processComm, bool prescribedEquations);
#endif

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
    void printReactionForces(TimeStep *tStep, int id);
    /**
     * Builds the reaction force table. For each prescribed equation number it will find
     * corresponding node and dof number. The entries in the restrDofMans, restrDofs, and eqn
     * arrays are sorted with increasing dofman number and with increasing dof number as
     * a second minor criterion.
     * @param restrDofMans Contains numbers of restrained Dofmanagers, with size equal to total number of prescribed equations.
     * @param restrDofs Contains numbers of restrained Dofs, with size equal to total number of prescribed equations.
     * @param eqn Contains the corresponding restrained equation numbers.
     * @param tStep Time step.
     * @param di Domain number.
     */
    void buildReactionTable(IntArray &restrDofMans, IntArray &restrDofs, IntArray &eqn, TimeStep *tStep, int di);

    /**
     * Computes the contribution of internal element forces to reaction forces in given domain.
     * @param reactions Contains the computed contributions.
     * @param tStep Solution step.
     * @param di Domain number.
     */
    void computeInternalForceReactionContribution(FloatArray &reactions, TimeStep *tStep, int di);
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
public:
    /// Creates new StructuralEngngModel with number i, associated to domain d.
    StructuralEngngModel(int i, EngngModel *_master = NULL) : EngngModel(i, _master),
        internalVarUpdateStamp(0), internalForcesEBENorm()
    { }
    /// Destructor.
    virtual ~StructuralEngngModel();

    // identification
    virtual const char *giveClassName() const { return "StructuralEngngModel"; }
    virtual classType giveClassID() const { return StructuralEngngModelClass; }

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
     * Exchanges necessary remote DofManagers data.
     * @param answer Array with collected values.
     * @param ExchangeTag Exchange tag used by communicator.
     * @return Nonzero if successful.
     */
    int updateSharedDofManagers(FloatArray &answer, int ExchangeTag);
    /**
     * Exchanges necessary remote prescribed DofManagers data.
     * @param answer Array with collected values.
     * @param ExchangeTag Exchange tag used by communicator.
     * @return Nonzero if successful.
     */
    int updateSharedPrescribedDofManagers(FloatArray &answer, int ExchangeTag);

    void initializeCommMaps(bool forceInit = false);

#ifdef __PARALLEL_MODE
    /**
     * Exchanges necessary remote element data with remote partitions. The receiver's nonlocalExt flag must be set.
     * Uses receiver nonlocCommunicator to perform the task using packRemoteElementData and unpackRemoteElementData
     * receiver's services.
     * @param ExchangeTag Exchange tag used by communicator.
     * @return Nonzero if successful.
     */
    int exchangeRemoteElementData(int ExchangeTag);
    /**
     * Packs data of local element to be received by their remote counterpart on remote partitions.
     * Remote elements are introduced when nonlocal constitutive models are used, in order to
     * allow local averaging procedure (remote elements, which are involved in averaging on local partition are
     * mirrored on this local partition) instead of implementing inefficient fine-grain communication.
     * Remote element data are exchanged only if necessary and once for all of them.
     * Current implementation calls packUnknowns service for all elements listed in
     * given process communicator send map.
     * @param processComm Corresponding process communicator.
     * @return Nonzero if successful.
     */
    int packRemoteElementData(ProcessCommunicator &processComm);
    /**
     * Unpacks data for remote elements (which are mirrors of remote partition's local elements).
     * Remote elements are introduced when nonlocal constitutive models are used, in order to
     * allow local averaging procedure (remote elements, which are involved in averaging on local partition are
     * mirrored on this local partition) instead of implementing inefficient fine-grain communication.
     * Remote element data are exchanged only if necessary and once for all of them.
     * Current implementation calls unpackAndUpdateUnknowns service for all elements listed in
     * given process communicator receive map.
     * @param processComm Corresponding process communicator.
     * @return Nonzero if successful.
     */
    int unpackRemoteElementData(ProcessCommunicator &processComm);
    /**
     * Packing function for vector values of DofManagers. Packs vector values of shared/remote DofManagers
     * into send communication buffer of given process communicator.
     * @param processComm Task communicator.
     * @param src Source vector.
     * @return Nonzero if successful.
     */
    int packDofManagers(FloatArray *src, ProcessCommunicator &processComm);
    /**
     * Packing function for vector values of DofManagers. Packs vector values of shared/remote prescribed DofManagers
     * into send communication buffer of given process communicator.
     * @param processComm Task communicator.
     * @param src Source vector.
     * @return Nonzero if successful.
     */
    int packPrescribedDofManagers(FloatArray *src, ProcessCommunicator &processComm);
    /**
     * Unpacking function for vector values of DofManagers . Unpacks vector of shared/remote DofManagers
     * from  receive communication buffer of given process communicator.
     * @param processComm Task communicator.
     * @param dest Destination vector.
     * @return Nonzero if successful.
     */
    int unpackDofManagers(FloatArray *dest, ProcessCommunicator &processComm);
    /**
     * Unpacking function for vector values of DofManagers . Unpacks vector of shared/remote prescribed DofManagers
     * from  receive communication buffer of given process communicator.
     * @param processComm Task communicator.
     * @param dest Destination vector.
     * @return Nonzero if successful.
     */
    int unpackPrescribedDofManagers(FloatArray *dest, ProcessCommunicator &processComm);

    ProblemCommunicator *giveProblemCommunicator(EngngModelCommType t) {
        if ( t == PC_default ) { return communicator; } else if ( t == PC_nonlocal ) { return nonlocCommunicator; } else { return NULL; }
    }

#endif
#ifdef __PETSC_MODULE
    /**
     * Creates PETSc contexts. Must be implemented by derived classes since the governing equation type is required
     * for context creation.
     */
    virtual void initPetscContexts();
#endif

#ifdef __OOFEG
    /**
     * Shows the sparse structure of required matrix, type == 1 stiffness.
     */
    void showSparseMtrxStructure(int type, oofegGraphicContext &context, TimeStep *atTime);
#endif
};
} // end namespace oofem
#endif // structengngmodel_h
