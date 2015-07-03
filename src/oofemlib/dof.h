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

#ifndef dof_h
#define dof_h

#include <cstdio>

#include "oofemcfg.h"
#include "floatarray.h"
#include "intarray.h"
#include "unknowntype.h"
#include "valuemodetype.h"
#include "dofiditem.h"
#include "doftype.h"
#include "contextioresulttype.h"

namespace oofem {
class DataStream;
class Dictionary;
class PrimaryField;
class Domain;
class DofManager;
class TimeStep;
class BoundaryCondition;
class InitialCondition;
class UnknownNumberingScheme;
class EntityRenumberingFunctor;

/**
 * Abstract class Dof represents Degree Of Freedom in finite element mesh.
 * DOFs are possessed by DofManagers (i.e, nodes, sides or whatever) and
 * one DOF belongs to only one DofManager.
 * Dof maintain its related equation or prescribed equation number.
 * This equation number is usually assigned by EngngModel, however, several
 * numbering schemes can exists (see giveEquationNumber and similar services).
 *
 * It maintains also its physical meaning and reference to
 * related DofManager (reference to DofManager which possess particular DOF).
 * To describe physical meaning of particular Dof, special enum type DofId has
 * been introduced (see cltypes.h). This type is more descriptive than
 * UnknownType, which determines physical meaning for unknowns generally
 * (displacement or temperature). DofId type must distinguish between
 * dofs having displacement unknown, but in different directions, because
 * only some of these may be required by elements.
 *
 * Dof can be subjected to boundary (BC) or initial (IC) condition. Method for
 * obtaining corresponding DOF unknown value is provided. If no IC condition
 * has been given, zero value IC is assumed otherwise when needed.
 *
 * Dof class generally supports changes of static system during computation.
 * This feature generally leads to equation renumbering. Then because equation number
 * associated to dof may change, it may become extremely complicated to ask EngngModel
 * for unknown from previous time step (because equation number may have been changed).
 * To overcome this problem, derived class will implement so called unknown dictionary,
 * which is updated after finishing each time step and where unknowns for particular
 * dof are stored. Dof then uses this dictionary for requests for unknowns instead of
 * asking EngngModel for unknowns. Unknowns in dof dictionary are updated by EngngModel
 * automatically (if EngngModel supports changes of static system) after finishing time
 * step.
 */
class OOFEM_EXPORT Dof
{
protected:
    /// Link to related DofManager.
    DofManager *dofManager;
    /// Physical meaning of DOF.
    DofIDItem dofID;

public:

    /**
     * Constructor. Creates DOF with number i, belonging to DofManager aNode and with
     * physical meaning described by id.
     * @param aNode DofManager which possess DOF.
     * @param id Physical meaning type.
     * @see cltypes.h, DofIDItem type
     */
    Dof(DofManager * aNode, DofIDItem id = Undef);
    /// Destructor.
    virtual ~Dof() { }

    /// Returns the type of the receiver.
    virtual dofType giveDofType() = 0;
    /// Returns class name of the receiver.
    virtual const char *giveClassName() const { return "Dof"; }

    /// @return Index of possessing DofManager.
    int giveDofManNumber() const;

    /// @return Associated DofManager.
    DofManager *giveDofManager() const { return dofManager; }

    int giveDofManGlobalNumber() const;

    /**
     * Returns value of boundary condition of dof if it is prescribed.
     * Use hasBc service to determine, if boundary condition is active.
     * The physical meaning of BC is determined by corresponding DOF.
     * @param mode Unknown char type (if total or incremental value is returned).
     * @param tStep Time step.
     * @return Prescribed value of unknown or zero if not prescribed.
     */
    virtual double giveBcValue(ValueModeType mode, TimeStep *tStep);
    /**
     * Returns equation number of receiver for given equation numbering scheme.
     * @param s Numbering scheme used to obtain equation number.
     * @return Equation number for DOF.
     */
    int giveEquationNumber(const UnknownNumberingScheme &s);
    /**
     * Returns equation number of receiver, usually assigned by emodel.
     * If Dof has active BC, returned equation number
     * is zero. After initializing Dof by calling constructor, Dof has no equation
     * number assigned. When firstly invoked, this function asks EngngModel object
     * for next equation prescribed equation number (this will increase also total number of equation
     * at EngngModel level). Note: By asking nodal code numbers or element code numbers
     * when initializing code numbers in EngngMode, designer should alter equation
     * numbering strategy.
     * @return Equation number of DOF.
     */
    virtual int __giveEquationNumber() const = 0;
    /**
     * Returns equation number of receiver. If Dof has active BC, returned equation number
     * is zero. After initializing Dof by calling constructor, Dof has no equation
     * number assigned. When firstly invoked, this function asks EngngModel object
     * for next equation prescribed equation number (this will increase also total number of equation
     * at EngngModel level). Note: By asking nodal code numbers or element code numbers
     * when initializing code numbers in EngngMode, designer should alter equation
     * numbering strategy.
     *
     * For slave dofs (dependent on other primary dofs) the array of master equation numbers is returned.
     * @param[out] masterEqNumbers Equation numbers of master DOFs for receiver.
     * @param[in] s Numbering scheme used enumeration of equations.
     */
    virtual void giveEquationNumbers(IntArray &masterEqNumbers, const UnknownNumberingScheme &s);
    /**
     * As giveEquationNumbers but for dof IDs.
     * @param[out] masterDofIDs Dof IDs of master DOFs for receiver.
     */
    virtual void giveDofIDs(IntArray &masterDofIDs);

    /**
     * Returns prescribed equation number of receiver. If Dof has inactive BC,
     * returned prescribed equation number is zero.
     * If Dof has active BC, then the corresponding prescribed equation number is returned.
     * is zero. After initializing Dof by calling constructor, Dof has no prescribed equation
     * number assigned. When firstly invoked, this function asks EngngModel object
     * for next equation or prescribed equation number (this will increase also total number of equation
     * at EngngModel level). Note: By asking nodal code numbers or element code numbers
     * when initializing code numbers in EngngMode, designer should alter equation
     * numbering strategy.
     * @return Prescribed equation number.
     */
    virtual int __givePrescribedEquationNumber() = 0;

    /**
     * Asks EngngModel for new equation number. Necessary for EngngModels supporting
     * changes of static system during solution. Then it is necessary to force
     * equation renumbering after finishing each time step.
     * @param tStep Time step determining the time.
     * @see Dof::updateUnknownsDictionary
     * @see EngngModel::requiresUnknownsDictionaryUpdate
     */
    virtual int askNewEquationNumber(TimeStep *tStep) = 0;
    /**
     * The key method of class Dof. Returns the value of the unknown of the receiver
     * at given time step. Unknown is characterized by its physical meaning (i.g., displacement)
     * an by its mode (e.g., value of displacement, velocity of displacement or acceleration of
     * displacement). UnknownType of requested unknown must be same as UnknownType of Dof.
     * @param mode Mode of unknown (e.g, total value, velocity or acceleration of unknown).
     * @param tStep Time step when unknown is requested. See documentation of particular EngngModel
     * class for valid tStep values (most implementation can return only values for current
     * and possibly for previous time step).
     * @return Value of unknown. If activeBC exist then returns value prescribed by BC. If tStep is time step
     * when IC apply, returns value given by this IC.
     */
    virtual double giveUnknown(ValueModeType mode, TimeStep *tStep) = 0;
    /**
     * The key method of class Dof. Returns the value of the unknown of the receiver
     * at given time step associated to given field.
     * @param field Field used to provide values.
     * @param mode Mode of unknown.
     * @param tStep Time step when unknown is requested. See documentation of particular EngngModel
     * class for valid tStep values (most implementation can return only values for current
     * and possibly for previous time step).
     * @return Value of unknown. If activeBC exist then returns value prescribed by BC. If tStep is time step
     * when IC apply, returns value given by this IC.
     */
    virtual double giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *tStep) = 0;
    /**
     * The key method of class Dof. Returns the value of the unknown of the receiver
     * at given time step associated to given field. For primary dof it returns is associated unknown value,
     * for slave dofs it returns an array of master values (in recursive way).
     * @param masterUnknowns Values of master unknowns for receiver.
     * @param mode Value mode for unknowns.
     * @param tStep Time step for when unknowns are requested.
     */
    virtual void giveUnknowns(FloatArray &masterUnknowns, ValueModeType mode, TimeStep *tStep);
    /**
     * The key method of class Dof. Returns the value of the unknown of the receiver
     * at given time step associated to given field. For primary dof it returns is associated unknown value,
     * for slave dofs it returns an array of master values (in recursive way).
     * @param masterUnknowns
     * @param field The field to pick unknowns from.
     * @param mode Value mode for unknowns.
     * @param tStep Time step for when unknowns are requested.
     */
    virtual void giveUnknowns(FloatArray &masterUnknowns, PrimaryField &field, ValueModeType mode, TimeStep *tStep);

    /**
     * Computes dof transformation array, which describes the dependence of receiver value on values of master dofs.
     * For primary dof, this transformation is unity, however, for slave DOFs, this array contains weights, which are multiplied by
     * corresponding master DOF values to obtain slave value.
     * @param masterContribs Contributions of master dofs for receiver.
     */
    virtual void computeDofTransformation(FloatArray &masterContribs);
    /**
     * @return Number of primary dofs, on which receiver value depends on (even recursively).
     */
    virtual int giveNumberOfPrimaryMasterDofs() { return 1; }

    /**
     * Test if Dof has active boundary condition.
     * @param tStep Time when test is evaluated.
     * @return True if active BC exists, false otherwise.
     */
    virtual bool hasBc(TimeStep *tStep) = 0;
    /**
     * Test if Dof has initial condition.
     * @return True if IC exists, false otherwise.
     */
    virtual bool hasIc() = 0;
    /**
     * Test if Dof has initial condition of required ValueModeType.
     * @param u Type of required IC.
     * @return True if IC exists, false otherwise.
     * @see ValueModeType
     */
    virtual bool hasIcOn(ValueModeType u) = 0;
    /**
     * Returns DofID value of receiver, which determines type of
     * of unknown connected to receiver (e.g., u-displacement, v-displacement, ...).
     */
    DofIDItem giveDofID() { return dofID; }
    /**
     * Sets the ID of receiver.
     */
    void setDofID(DofIDItem id) { this->dofID = id; }
    /**
     * Tests if receiver is primary DOF. Dof is primary if it possess or directly represent
     * certain DOF. If it is linked somehow (rigid arm, doubled node) to other DOF(s) then it is not
     * primary DOF.
     * @return True if receiver is primary DOF, false otherwise (default).
     */
    virtual bool isPrimaryDof() { return false; }
    /**
     * Returns the id of associated boundary condition, if there is any.
     * Used only for printing purposes. In general, id could not be used
     * to decide whether BC is active. Use appropriate services instead.
     * @return Id of associated Boundary condition, zero otherwise
     */
    virtual int giveBcId() = 0;
    /**
     * Returns the id of associated initial condition, if there is any.
     * Used only for printing purposes. In general, id could not be used
     * to decide whether IC is active. Use appropriate services instead.
     * @return Id of associated initial condition, zero otherwise
     */
    virtual int giveIcId() = 0;

    /**
     * @return Array of master DofManagers to which the receiver is linked.
     */
    virtual void giveMasterDofManArray(IntArray &answer);
    /**
     * Local renumbering support. For some tasks (parallel load balancing, for example) it is necessary to
     * renumber the entities. The various FEM components (such as nodes or elements) typically contain
     * links to other entities in terms of their local numbers, etc. This service allows to update
     * these relations to reflect updated numbering. The renumbering function is passed, which is supposed
     * to return an updated number of specified entity type based on old number.
     * @param f Function that converts old to new equation number.
     */
    virtual void updateLocalNumbering(EntityRenumberingFunctor &f) { }

    /**
     * Prints Dof output (it prints value of unknown related to dof at given timeStep).
     * The format of output depends on analysis type.
     * Called from corresponding e-model.
     */
    virtual void printSingleOutputAt(FILE *file, TimeStep *tStep, char ch, ValueModeType mode, double scale = 1.0);
    /**
     * Prints Dof output (it prints value of unknown related to dof at given timeStep).
     * The format of output depends on analysis type.
     * Called from corresponding e-model.
     */
    virtual void printMultipleOutputAt(FILE *File, TimeStep *tStep, char *ch, ValueModeType *mode, int nite);
    void printSingleOutputWithAdditionAt(FILE *File, TimeStep *tStep, char ch, ValueModeType mode, double addend);
    /// Prints the receiver state on stdout.
    virtual void printYourself();
    /**
     * Updates receiver after finishing time step.
     * @param tStep Finished time step.
     */
    virtual void updateYourself(TimeStep *tStep) { }
    /**
     * Abstract function, allowing Dof to store its unknowns in its own private
     * dictionary. Dof then uses this dictionary instead of forwarding the requests to
     * EngngModel (with equationNUmber as parameter). If EngngModel does not support changes
     * of static system (see EngngModel::requiresUnknownsDictionaryUpdate method), the dof
     * forwards the requests for its unknowns to EngngModel, where unknowns are naturally kept.
     * This is possible, because dof equation number is same during whole solution.
     * But when changes of static system are allowed, several problem arise. For example
     * by solving simple  incremental static with allowed static changes, the incremental displacement
     * vector of structure can not be added to total displacement vector of structure, because
     * equation numbers may have changed, and one can not simply add these vector to obtain new
     * total displacement vector, because incompatible displacement will be added.
     * To solve this problem, unknown dictionary at dof level has been assumed. Dof then keeps
     * its unknowns in its own private dictionary.
     * After computing increment of solution, engngModel updates for each dof its unknowns  in its
     * dictionary (using updateUnknownsDictionary function). For aforementioned example
     * engngModel updates incremental values but also total value by asking dof for previous total
     * value (dof will use its dictionary, does not asks back EngngModel) adds corresponding increment
     * and updates total value in dictionary.
     * In fact on EngngModel level only incremental solution is stored, but total values are
     * always stored in dofs dictionaries.
     * Implementation is not provided, only interface declared. Children must implement this method.
     * @param tStep Time step when unknowns are updated. In current version it is unused parameter.
     * It is EngngModel responsibility to update values, and values stored in dictionary
     * are always related to timeStep when they were lastly updated.
     * @param mode Mode of stored unknown.
     * @param dofValue Value of unknown. Old value will generally be lost.
     * @see EngngModel::requiresUnknownsDictionaryUpdate
     */
    virtual void updateUnknownsDictionary(TimeStep *tStep, ValueModeType mode, double dofValue) { }
    /**
     * Access dictionary value, if not present zero is returned.
     * @param tStep Time step.
     * @param mode Mode of value.
     * @return The requested value.
     */
    virtual double giveUnknownsDictionaryValue(TimeStep *tStep, ValueModeType mode) { return 0; }

    /// Returns string for prepending output (used by error reporting macros).
    std :: string errorInfo(const char *func) const;
    /// Stores receiver state to output stream.
    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    /// Restores the receiver state previously written in stream.
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    /// Overwrites the boundary condition id (0-inactive BC), intended for specific purposes such as coupling of bc's in multiscale simulations
    virtual void setBcId(int bcId) { }
    /// Overwrites the initial condition id (0-inactive IC)
    virtual void setIcId(int icId) { }

    /**
     * Sets a specific equation number to receiver.
     * @param equationNumber New equation number.
     */
    virtual void setEquationNumber(int equationNumber) { }

    /**
     * Sets the dictionary of unknowns for receiver.
     * @param unknowns New dictionary of unknowns.
     */
    virtual void setUnknowns(Dictionary *unknowns) { }
    /**
     * Receives the dictionary of unknowns in receiver.
     * @return Dictionary of unknowns for DOF.
     */
    virtual Dictionary *giveUnknowns() { return NULL; }

    /**
     * Gives number for equation, negative for prescribed equations.
     * @return Equation number of receiver.
     */
    virtual int giveEqn() { return 0; }

    /**
     * Packs specific  DOF Manager's dofs unknowns into communication buffer.
     * If dof is slave, then no packing is done, this is maintained by master. This requires master
     * be available at same partition as slave.
     * @param buff Communication buffer to pack data.
     * @param mode Mode of unknown (e.g, total value, velocity or acceleration of unknown).
     * @param tStep Time step when unknown requested. See documentation of particular EngngModel
     * class for valid tStep values (most implementations can return only values for current
     * and possibly for previous time step).
     * @return Nonzero if successful.
     * @todo Remove this? It is not inherited by MasterDof. Is this leftovers? / Mikael
     */
    virtual int packUnknowns(DataStream &buff, ValueModeType mode, TimeStep *tStep)
    { return 1; }
    /**
     * Unpacks DOF unknown from communication buffer and updates unknown if necessary.
     * Unknown is always updated using EngngModel::updateUnknownComponent, if DOFManager
     * to which receiver belongs has DofManager_shared dofManagerParallelMode type.
     * Unknown is unpacked and stored in unknowns dictionary, if DOFManager
     * to which receiver belongs has DofManager_remote dofManagerParallelMode type.
     * There is no reason for invoking this service if DOFManager has DofManager_local mode.
     * If do is slave, then no unpacking and updating is done. This is left on master, which must be
     * available on same partition.
     * @param buff Buffer containing packed message.
     * @param mode Mode of unknown (e.g, total value, velocity or acceleration of unknown).
     * @param tStep Time step when unknown requested. See documentation of particular EngngModel
     * class for valid tStep values (most implementations can return only values for current
     * and possibly for previous time step).
     * @return Nonzero if successful.
     * @todo Remove this? It is not inherited by MasterDof. Is this leftovers? / Mikael
     */
    virtual int unpackAndUpdateUnknown(DataStream &buff,
                                       ValueModeType mode, TimeStep *tStep) { return 1; }

protected:
    /**
     * Returns boundary condition of dof if it is prescribed.
     * @return NULL if no BC applied, otherwise pointer to corresponding BC.
     */
    virtual BoundaryCondition *giveBc() { return NULL; }
    /**
     * Returns initial condition of dof if it is prescribed.
     * @return NULL if no IC applied, otherwise pointer to corresponding IC.
     */
    virtual InitialCondition *giveIc() { return NULL; }

    friend class SimpleSlaveDof;
};
} // end namespace oofem
#endif // dof_h
