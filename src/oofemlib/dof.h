/* $Header: /home/cvs/bp/oofem/oofemlib/src/dof.h,v 1.13.4.1 2004/04/05 15:19:43 bp Exp $ */
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



//   *****************
//   *** CLASS DOF ***
//   *****************


#ifndef dof_h
#define dof_h

#include "compiler.h"
#include "dictionr.h"
#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <string.h>
#endif

#include "error.h"
#include "primaryfield.h"
#include "classtype.h"
#include "unknowntype.h"
#include "equationid.h"
#include "valuemodetype.h"
#include "dofiditem.h"
#include "contextioresulttype.h"
#include "entityrenumberingscheme.h"

#ifdef __PARALLEL_MODE
class CommunicationBuffer;
#endif

class Domain;
class DofManager;
class TimeStep;
class BoundaryCondition;
class InitialCondition;

/**
 * Abstract class Dof represents Degree Of Freedom in finite element mesh.
 * DOFs are possesed by DofManagers (i.e, nodes, sides or whatever) and
 * one DOF belongs to only one DofManager.
 * Dof maintain its related equation or prescribed equation number,
 * physical meaning and reference to
 * related DofManager (reference to DofManager which possess particular DOF).
 * To describe physical meaning of particular Dof, special enum type "DofId" has
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
 * associated to dof may change, it may become extremly complicated to ask EngngModel
 * for unknown from previous time step (because equation number may have been changed).
 * To overcome this problem, derived class will implement so called uknown dictionary,
 * which is updated after finishing each time step and where unknowns for particular
 * dof are stored. Dof then uses this dictionary for requests for unknowns instead of
 * asking EngngModel for unknowns. Unknowns in dof dictionary are updated by EngngModel
 * automatically (if EngngModel supports changes of static system) after finishing time
 * step.
 */
class Dof
{
    /*
     * This class implements an abstract class for  nodal (or whatever i.e. side)
     * degree of freedom. A dof is usually attribute of one node.
     * DESCRIPTION
     * 'number' and 'node' are used for reading/writing data in the data file.
     * 'DofID' is parameter determining physical meaning of receiver.
     * This parameter is also used in member function giveUnknownType, which returns
     * UnknownType type according to DofID parameter.
     *
     * I don't know whether to implement following feature (now  not implemented)
     * 'unknowns' and 'pastUnknowns' are the dictionaries where the dof stores
     * its unknowns (e.g., the displacement 'd', the velocity 'v' and the acce-
     * leration 'a'), at the current time step and at the previous one.
     *
     * TASKS
     * - equation numbering, in method 'giveEquationNumber'and 'givePrescribedEquationNumber';
     * - managing its b.c. and its i.c., if any (methods 'hasBc', 'giveBc', etc);
     * - managing its unknowns. This includes retrieving the associated solution
     *   from the Engng. System , or from receiver's dictionary
     * (based on emodel->requiresNodeUnknowsDictionaryUpdate() function,
     * which determines whether to use dictionary or ask unknowns values from
     * emodel)
     * - managing physical meaning of dof (dofID variable)
     *
     * REMARKS
     * - class Dof is not a subclass of FEMComponent : a dof belongs to a single
     *   node, not to the domain ;
     * - class Dof is not restricted to structural analysis problems. Unknowns
     *   may also be pressures, temperatures, etc.
     * - method give returns unknown value quantity according to ValueModeType parameter,
     *  UnknownType parameter is used to check whether pysical meaning of
     * unknown corresponds.
     *
     */

protected:
    /// Dof number.
    int number;
    /// Link to related DofManager.
    DofManager *dofManager;
    /// Physical meaning of DOF.
    DofIDItem dofID;
    /*      Dictionary*  unknowns ;
     * Dictionary*  pastUnknowns ; */

public:

    /**
     * Constructor. Creates DOF with number i, belonging to DofManager aNode and with
     * physical meaning described by id.
     * @param i DOF number.
     * @param aNode DofManager which possess DOF.
     * @param id Physical meaning type.
     * @see cltypes.h, DofID type
     */
    Dof(int i, DofManager *aNode, DofID id = Undef);   // constructor
    /// Destructor.
    virtual ~Dof()   { }  // destructor.

    /// Returns class name of the receiver.
    virtual const char *giveClassName() const { return "Dof"; }
    /// Returns classType id of receiver.
    virtual classType           giveClassID() const { return DofClass; }

    /// Returns receiver number.
    int giveNumber() const { return number; }

    int giveDofManNumber() const;
#ifdef __PARALLEL_MODE
    int giveDofManGlobalNumber() const;
#endif
    /**
     * Returns value of boundary condition of dof if it is precsribed.
     * Use hasBc service to determine, if boundary condition is active.
     * The physical meaning of Bc is determined by corresponding DOF.
     * @param mode unknown char type (if total or incremental value is returned)
     * @param tStep time step
     * @return prescribed value of unknown or zero if not prescribed
     */
    virtual double  giveBcValue(ValueModeType mode, TimeStep *tStep);
    /**
     * Returns array with values of boundary condition of dof if it is precsribed.
     * For primary dof it has only one value, for slave dof the value of corresponding
     * masters dofs boundary condition is assembled and returned.
     * Use hasBc service to determine, if boundary condition is active.
     * The physical meaning of Bc is determined by corresponding DOF.
     */
    virtual void giveBcValues(FloatArray &masterBcValues, ValueModeType mode, TimeStep *stepN)
    {
        masterBcValues.resize(1);
        masterBcValues.at(1) = this->giveBcValue(mode, stepN);
    }

    //  virtual double  giveBcValue (UnknownType type, ValueModeType mode, TimeStep* tStep) ;
    /**
     * Returns equation number of receiver. If Dof has active BC, returned equation number
     * is zero. After initializing Dof by calling constructor, Dof has no equation
     * number assigned. When firstly invoked, this function asks EngngModel object
     * for next equation prescribed equation number (this will increase also total number of equation
     * at EngngModel level). Note: By asking nodal code numbers or element code numbers
     * when initializing code numbers in EngngMode, designer should alter equation
     * numbering strategy.
     */
    virtual int                 giveEquationNumber()   = 0;
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
     */
    virtual void giveEquationNumbers(IntArray &masterEqNumbers)
    {
        masterEqNumbers.resize(1);
        masterEqNumbers.at(1) = this->giveEquationNumber();
    }

    /**
     * Returns prescribed equation number of receiver. If Dof has inactive BC,
     * returned prescribed equation number is zero.
     * If Dof has active BC, then the corresponding  prescribed equation number is returned.
     * is zero. After initializing Dof by calling constructor, Dof has no prescribed equation
     * number assigned. When firstly invoked, this function asks EngngModel object
     * for next equation or prescribed equation number (this will increase also total number of equation
     * at EngngModel level). Note: By asking nodal code numbers or element code numbers
     * when initializing code numbers in EngngMode, designer should alter equation
     * numbering strategy.
     */
    virtual int                 givePrescribedEquationNumber()   = 0;
    /**
     * Returns prescribed equation number of receiver. If Dof has inactive BC,
     * returned prescribed equation number is zero.
     * If Dof has active BC, then the corresponding  prescribed equation number is returned.
     * is zero. After initializing Dof by calling constructor, Dof has no prescribed equation
     * number assigned. When firstly invoked, this function asks EngngModel object
     * for next equation or prescribed equation number (this will increase also total number of equation
     * at EngngModel level). Note: By asking nodal code numbers or element code numbers
     * when initializing code numbers in EngngMode, designer should alter equation
     * numbering strategy.
     *
     * For slave dofs (dependent on other primary dofs) the array of master equation numbers is returned.
     */

    virtual void givePrescribedEquationNumbers(IntArray &masterEqNumbers)
    {
        masterEqNumbers.resize(1);
        masterEqNumbers.at(1) = this->givePrescribedEquationNumber();
    }

    /**
     * Asks EngngModel for new equation number. Necessary for EngngModels supporting
     * changes of static system during solution. Then it is necessary to force
     * equation renumbering after finishing each time step.
     * @param tStep time step determining the time
     * @see Dof::updateUnknownsDictionary function.
     * @see EngngModel::requiresUnknownsDictionaryUpdate.
     */
    virtual int                 askNewEquationNumber(TimeStep *tStep) = 0;
    //      double              givePastUnknown (char,TimeStep*) ;
    //      double              giveUnknown (char,TimeStep*) ;
    /**
     * The key method of class Dof. Returns the value of the unknown of the receiver
     * at given time step. Unknown is characterized by its physical meaning (i.g., displacement)
     * an by its mode (e.g., value of displacement, velocity of displacement or acceleration of
     * displacement). UnknownType of requested unknown must be same as UnknownType of Dof.
     * @param type physical meaning of  unknown.
     * @param mode mode of unknown (e.g, total value, velocity or acceleration of unknown).
     * @stepN time step when unknown requested. See documentation of particular EngngModel
     * class for valid StepN values (most implementaion can return only values for current
     * and possibly for previous time step).
     * @return unknown value, if activeBC exist then returns value prescribed by BC. If stepN is time step
     * when IC apply, returns value given by this IC.
     */
    virtual double              giveUnknown(EquationID type, ValueModeType mode, TimeStep *stepN)  = 0;
    /**
     * The key method of class Dof. Returns the value of the unknown of the receiver
     * at given time step associated to given field.
     * @param field field used to provide values.
     * @stepN time step when unknown requested. See documentation of particular EngngModel
     * class for valid StepN values (most implementaion can return only values for current
     * and possibly for previous time step).
     * @return unknown value, if activeBC exist then returns value prescribed by BC. If stepN is time step
     * when IC apply, returns value given by this IC.
     */
    virtual double              giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *stepN)  = 0;
    /**
     * The key method of class Dof. Returns the value of the unknown of the receiver
     * at given time step associated to given field. For primary dof it returns is associated unknown value,
     * for slave dofs it returns an array of master values (in recursive way).
     */
    virtual void giveUnknowns(FloatArray &masterUnknowns, EquationID type, ValueModeType mode, TimeStep *stepN)
    {
        masterUnknowns.resize(1);
        masterUnknowns.at(1) = this->giveUnknown(type, mode, stepN);
    }
    /**
     * The key method of class Dof. Returns the value of the unknown of the receiver
     * at given time step associated to given field. For primary dof it returns is associated unknown value,
     * for slave dofs it returns an array of master values (in recursive way).
     */
    virtual void giveUnknowns(FloatArray &masterUnknowns, PrimaryField &field, ValueModeType mode, TimeStep *stepN)
    {
        masterUnknowns.resize(1);
        masterUnknowns.at(1) = this->giveUnknown(field, mode, stepN);
    }

    /**
     * Computes dof transformation array, which describes the dependence of receiver value on values of master dofs.
     * For primary dof, this thransformation is unity, hovewer, for slave DOFs, thisarray contains weights, which are multiplied by
     * corresponding master DOF values to obtain slave value.
     */
    virtual void computeDofTransformation(FloatArray &masterContribs)
    {
        masterContribs.resize(1);
        masterContribs.at(1) = 1.0;
    }
    /**
     * Returns number of primary dofs, on which receiver value depends on (even recursivelly)
     */
    virtual int giveNumberOfPrimaryMasterDofs() { return 1; }

    /**
     * Test if Dof has active boundary condition.
     * @param tStep time when test is evaluated.
     * @return nonzero if active BC exists, zero otherwise.
     */
    virtual int                 hasBc(TimeStep *tStep)  = 0;
    /**
     * Test if Dof has initial condition.
     * @return nonzero if IC exists, zero otherwise.
     */
    virtual int                 hasIc()  = 0;
    /**
     * Test if Dof has initial condition of required ValueModeType.
     * @param u type of required IC
     * @return nonzero if IC exists, zero otherwise.
     * @see ValueModeType.
     */
    virtual int                 hasIcOn(ValueModeType u)  = 0;
    /**
     * Returns DofID value of receiver, which determines type of
     * of unknown connected to receiver (e.g., u-displacement, v-displacement, ...).
     */
    DofIDItem               giveDofID()  { return dofID; }
    /**
     * Returns char representation of DofID value of receiver, which determines physical meaning
     * of unknown connected to receiver. Usefull only for printing. More conviniently,
     * one should use giveDofID function.
     * @see Dof::giveDofID.
     */
    char *giveDofIDName(char *s);
    /**
     * Return CharType value of unknown related to receiver. (e.g., displacement).
     */
    UnknownType            giveUnknownType();
    /**
     * Tests if receiver is primary DOF. Dof is primary if it posses or directly represent
     * certain DOF. If it is linked somehow (rigid arm, doubled node) to other DOF(s) then it is not
     * primary DOF.
     * @return nonzero if receiver is primary DOF, zero otherwise (default).
     */
    virtual int            isPrimaryDof() { return 0; }
    /** Returns the id of associated boundary condition, if there is any.
     * Used only for printing purposes. In general, id could not be used
     * to decide whether bc is active. Use appropriate services instead.
     * @param id of associated Boundary condition, zero otherwise
     */
    virtual int giveBcId() = 0;
    /** Returns the id of associated initial condition, if there is any.
     * Used only for printing purposes. In general, id could not be used
     * to decide whether bc is active. Use appropriate services instead.
     * @param id of associated initial condition, zero otherwise
     */
    virtual int giveIcId() = 0;

    /**
     * Returns an array of master DofManagers  to which the recever is linked
     */
    virtual void giveMasterDofManArray(IntArray &answer) { answer.resize(0); } // termitovo
    /**
     * Local renumbering support. For some tasks (parallel load balancing, for example) it is necessary to
     * renumber the entities. The various fem components (such as nodes or elements) typically contain
     * links to other entities in terms of their local numbers, etc. This service allows to update
     * these relations to reflext updated numbering. The renumbering funciton is passed, which is supposed
     * to return an updated number of specified entyty type based on old number.
     */
    virtual void updateLocalNumbering( EntityRenumberingFunctor &f ) {}

    /**
     * Prints Dof output (it prints value of unknown related to dof at given timeStep).
     * The format of output depends on analysis type.
     * Called from corresponding e-model.
     */
    virtual void        printSingleOutputAt(FILE *, TimeStep *stepN, char ch, EquationID type, ValueModeType mode, double scale = 1.0);
    /**
     * Prints Dof output (it prints value of unknown related to dof at given timeStep).
     * The format of output depends on analysis type.
     * Called from corresponding e-model.
     */
    virtual void        printMultipleOutputAt(FILE *File, TimeStep *stepN, char *ch,
                                              EquationID type, ValueModeType *mode, int nite);

    /// Prints the receiver state on stdin.
    virtual void       printYourself();
    /// Updates receiver after finishing time step.
    virtual void       updateYourself(TimeStep *) { }
    /**
     * Abstract function, allowing Dof to store its unknowns in its own private
     * dictionary. Dof then uses this dictionary instead of forwarding the requests to
     * EngngModel (with equationNUmber as parameter). If EngngModel does not support changes
     * of static system (see EngngModel::requiresUnknownsDictionaryUpdate method), the dof
     * frowards the requests for its unknowns to EngngModel, where unknowns are naturaly kept.
     * This is posible, because dof equation number is same during whole solution.
     * But when changes of static system are allowed, several problem arise. For example
     * by solving simple  incremental static with allowed static changes, the incremetal displacement
     * vector of structure can not be added to total displacement vector of structure, because
     * equation numbers may have changed, and one can not simply add these vector to obtain new
     * total displacement vector, because uncompatible displacement will be added.
     * To solve this problem, uknown dictionary at dof level has been assumed. Dof then keeps
     * its unknowns in its onw private dictionary.
     * After computing increment of solution, engngModel updates for each dof its unknowns  in its
     * dictionary (using updateUnknownsDictionary function). For aforementioned example
     * engngModel updates incremental values but also total value by asking dof for previous total
     * value (dof will use its dictionary, does not asks back EngngModel) adds corresponding increment
     * and updates total value in dictionary.
     * In fact on EngngModel level only incremental solution is stored, but total values are
     * always stored in dofs dictionaries.
     * Implementaion is not provided, only interface declared. Children must implement this method.
     * @param tStep time step when unknowns are updated. In current version it is unused parameter.
     * It is EngngModel responsibility to update values, and values stored in dictionary
     * are always related to timeStep when they were lastly updated.
     * @param type identifies type of unknown. It is not possible to store values of different
     * UnknownType types then  UnknownType type of receiver.
     * @param mode mode of stored unknown.
     * @param dofValue value of unknown. Old value will generaly be lost.
     * @see EngngModel::requiresUnknownsDictionaryUpdate
     *
     */
    virtual void       updateUnknownsDictionary(TimeStep *tStep, EquationID type,
                                                ValueModeType mode, double dofValue) { }
    /** access dictionary value, if not present zero is returned */
    virtual void       giveUnknownsDictionaryValue(TimeStep *tStep, EquationID type,
                                                   ValueModeType mode, double &dofValue) { }

    /// prints simple error message and exits
    void error(const char *file, int line, const char *format, ...);
    /// Stores receiver state to output stream.
    virtual contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /// Restores the receiver state previously written in stream.
    virtual contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual void setEquationNumber(int equationNumber) {}; // rch
    virtual void setUnknowns(Dictionary *unknowns) {}; // rch
    virtual Dictionary *giveUnknowns() { return NULL; } // rch
    virtual int giveEqn() {return 0;}

#ifdef __PARALLEL_MODE
    /**
     * Packs specific  DOF Manager's dofs unknowns into communication buffer.
     * If dof is slave, then no packing is done, this is maintained by master. This requires master
     * be available at same partition as slave.
     * @param buff communication buffer to pack data.
     * @param type physical meaning of  unknown.
     * @param mode mode of unknown (e.g, total value, velocity or acceleration of unknown).
     * @stepN time step when unknown requested. See documentation of particular EngngModel
     * class for valid StepN values (most implementaion can return only values for current
     * and possibly for previous time step).
     * @return nonzero if succesfull
     */
    virtual int packUnknowns(CommunicationBuffer &buff, EquationID type, ValueModeType mode, TimeStep *stepN)
    { return 1; }
    /**
     * Unpacks DOF unknown from communication buffer and updates unknown if necessary.
     * Unknown is always updated using EngngModel::updateUnknownComponent, if DOFManager
     * to which receiver belongs has DofManager_shared dofManagerParallelMode type.
     * Unknown is unpacked and stored in unknowns dictionary, if DOFManager
     * to which receiver belongs has DofManager_remote dofManagerParallelMode type.
     * There is no reason for invoking this service if  DOFManager has DofManager_local mode.
     * If do is slave, then no unpacking and updating is done. This is left on master, which must be
     * available on same partition.
     * @see Dof::unpackAndUpdateDOFsUnknown for description
     * @param buff buffer containing packed message
     * @param type physical meaning of  unknown.
     * @param mode mode of unknown (e.g, total value, velocity or acceleration of unknown).
     * @stepN time step when unknown requested. See documentation of particular EngngModel
     * class for valid StepN values (most implementaion can return only values for current
     * and possibly for previous time step).
     * @return nonzero if succesfull
     */
    virtual int unpackAndUpdateUnknown(CommunicationBuffer &buff, EquationID type,
                                       ValueModeType mode, TimeStep *stepN) { return 1; }
#endif

protected:
    /**
     * Returns boundary condition of dof if it is precsribed.
     * @return returns NULL if no BC applied, otherwise pointer to correcpondig BC.
     */
    virtual BoundaryCondition *giveBc() { return NULL; }
    /**
     * Returns initial condition of dof if it is precsribed.
     * @return returns NULL if no IC applied, otherwise pointer to correcpondig IC.
     */
    virtual InitialCondition *giveIc() { return NULL; }

    friend class SimpleSlaveDof;
};

#endif // dof_h
