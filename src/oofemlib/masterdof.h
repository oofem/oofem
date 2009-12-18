/* $Header: /home/cvs/bp/oofem/oofemlib/src/masterdof.h,v 1.12 2003/04/06 14:08:25 bp Exp $ */
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



//   ************************
//   *** CLASS MASTER DOF ***
//   ************************


#ifndef masterdof_h
#define masterdof_h

#include "dof.h"
#include "compiler.h"
#include "dictionr.h"
#ifndef __MAKEDEPEND
 #include <stdio.h>
 #include <string.h>
#endif

#include "error.h"
#include "dofmanager.h"

#ifdef __PARALLEL_MODE
 #include "combuff.h"
#endif

namespace oofem {

class Domain;
class DofManager;
class TimeStep;
class BoundaryCondition;
class InitialCondition;

/**
 * Class representing "master" degree of freedom. Master is degree of freedom, which has
 * its related unknown and corresponding equation number.
 */
class MasterDof  : public Dof
{
    /*
     * This class implements a nodal degree of freedom. A dof is usually attri-
     * bute of one node.
     * DESCRIPTION
     * 'number' and 'node' are used for reading/writing data in the data file.
     * 'equationNumber' keeps the number of the associated equation in the linear
     * system (value>0) or the associated prescribed equation number if the dof
     * is subjected to a boundary condition (b.c.) (value<0).
     * 'bc' is the number of the b.c. the dof is subjected to, if any.
     * 'DofID' is parameter determining physical meaning of receiver.
     * This parameter is also used in member function giveUnknownType, which returns
     * CharType type according to DofID parameter.
     *
     * I don't know whether to implement following feature (now  not implemented)
     * 'unknowns' and 'pastUnknowns' are the dictionaries where the dof stores
     * its unknowns (e.g., the displacement 'd', the velocity 'v' and the acce-
     * leration 'a'), at the current time step and at the previous one.
     *
     * TASKS
     * - equation numbering, in method 'giveEquationNumber' and 'givePrescribedEquationNumber' ;
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
     *  CharType parameter is used to check whether pysical meaning of
     * unknown corresponds.
     *
     */

protected:
    /// corresponding equation number (positive value) or prescribed equation number (negative value)
    int equationNumber;
    /// Boundary condition number associated to dof.
    int bc;
    /// Initial condition number associated to dof.
    int ic;
    /// Unknowns dictionary to support changes of static system.
    Dictionary *unknowns;
    /*      Dictionary*  unknowns ;
     *      Dictionary*  pastUnknowns ; */

public:
    /**
     * Constructor. Creates master dof with number i, belonging to DofManager aNode and with
     * physical meaning described by id.
     * @param i DOF number.
     * @param aNode DofManager which possess DOF.
     * @param nbc number of associated boundary condition, zero if none.
     * @param nic number of associated initial condition, zero if none.
     * @param id Physical meaning type.
     * @see cltypes.h, DofID type
     */
    MasterDof(int i, DofManager *aNode, int nbc, int nic, DofID id);     // constructor
    MasterDof(int i, DofManager *aNode, DofID id = Undef);
    /// Destructor.
    ~MasterDof()   { delete unknowns; /*delete unknowns ; delete pastUnknowns ;*/ }      // destructor.
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "MasterDof"; }
    /// Returns classType id of receiver.
    classType           giveClassID() const { return MasterDofClass; }
    /**
     * Returns assigned equation number of receiver. If Dof has active BC, returned equation number
     * is zero. After initializing Dof by calling constructor, Dof has no equation
     * number assigned. When firstly invoked, this function asks EngngModel object
     * for next equation prescribed equation number (this will increase also total number of equation
     * at EngngModel level). Note: By asking nodal code numbers or element code numbers
     * when initializing code numbers in EngngMode, designer should alter equation
     * numbering strategy.
     */
    virtual int                 __giveEquationNumber();
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
    virtual int                 __givePrescribedEquationNumber();
    /**
     * Asks EngngModel for new equation number. Necessary for EngngModels supporting
     * changes of static system during solution. Then it is necessary to force
     * equation renumbering after finishing each time step.
     * @param tStep time step determining the time
     * @see Dof::updateUnknownsDictionary function.
     * @see EngngModel::requiresUnknownsDictionaryUpdate.
     */
    int                 askNewEquationNumber(TimeStep *tStep);
    //      double              givePastUnknown (char,TimeStep*) ;
    //      double              giveUnknown (char,TimeStep*) ;
    /**
     * Returns the value of the unknown associated with the receiver
     * at given time step. Unknown is characterized by its physical meaning (i.g., displacement)
     * an by its mode (e.g., value of displacement, velocity of displacement or acceleration of
     * displacement). CharType of requested unknown must be same as CharType of Dof.
     * @param type physical meaning of  unknown.
     * @param mode mode of unknown (e.g, total value, velocity or acceleration of unknown).
     * @stepN time step when unknown requested. See documentation of particular EngngModel
     * class for valid StepN values (most implementaion can return only values for current
     * and possibly for previous time step).
     * @return unknown value, if activeBC exist then returns value prescribed by BC. If stepN is time step
     * when IC apply, returns value given by this IC.
     */
    double giveUnknown(EquationID, ValueModeType, TimeStep *);
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
    double giveUnknown(PrimaryField & field, ValueModeType, TimeStep * stepN);
    /**
     * Test if Dof has active boundary condition.
     * @param tStep time when test is evaluated.
     * @return nonzero if active BC exists, zero otherwise.
     */
    int                 hasBc(TimeStep *tStep);
    /**
     * Test if Dof has initial condition.
     * @return nonzero if IC exists, zero otherwise.
     */
    int                 hasIc();
    /**
     * Test if Dof has initial condition of required ValueModeType.
     * @param u type of required IC
     * @return nonzero if IC exists, zero otherwise.
     * @see ValueModeType.
     */
    int hasIcOn(ValueModeType);
    // DofID               giveDofID () {return dofID;}
    // char*              giveDofIDName (char* s);
    // CharType            giveUnknownType ();

    /**
     * Tests if receiver is primary DOF. Dof is primary if it posses or directly represent
     * certain DOF. If it is linked somehow (rigid arm, doubled node) to other DOF(s) then it is not
     * primary DOF.
     * @returns nonzero indicating primary DOF.
     */
    int                 isPrimaryDof() { return 1; }
    /** Returns the id of associated boundary condition, if there is any.
     * Used only for printing purposes. In general, id culd not be used
     * to decide whether bc is active. Use appropriate services instead.
     * @param id of associated Boubdaray condition, zero otherwise
     */
    int giveBcId();
    /** Returns the id of associated initial condition, if there is any.
     * Used only for printing purposes. In general, id could not be used
     * to decide whether bc is active. Use appropriate services instead.
     * @param id of associated initial condition, zero otherwise
     */
    int giveIcId();

    /// Prints the receiver state on stdin.
    void                printYourself();
    /// Updates receiver after finishing time step.
    void                updateYourself(TimeStep *);
    /**
     * Updates entry in unknowns dictionary for given unknown (determined form type and mode parameters)
     * to given value. Dof dictionary is used only if EngngModel supports changes of static system.
     * When dof is asked for unknown, uses values stored in its private dictionary rather than requesting
     * EngngModel for this value. But EngngModel must correctly update these values after finishing each time step.
     * If EngngModel does not support changes of static system
     * (recognized calling  EngngModel::requiresUnknownsDictionaryUpdate method), the dof
     * frowards the requests for its unknowns to EngngModel, where unknowns are naturaly kept.
     * @see Dof::updateUnknownsDictionary for detailed explanation
     * @see EngngModel::requiresUnknownsDictionaryUpdate
     */
    void                updateUnknownsDictionary(TimeStep *tStep, EquationID type,
                                                 ValueModeType mode, double dofValue);
    /** access dictionary value, if not present zero is returned */
    virtual void       giveUnknownsDictionaryValue(TimeStep *tStep, EquationID type,
                                                   ValueModeType mode, double &dofValue);

    /**
     * Stores receiver state to output stream.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    contextIOResultType    saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores the receiver state previously written in stream.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    contextIOResultType    restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    
    void setBcId(int bcId) { this->bc = bcId; }
    void setEquationNumber(int equationNumber) { this->equationNumber = equationNumber; } // rch
    void setUnknowns(Dictionary *unknowns) { this->unknowns = unknowns; } // rch
    Dictionary *giveUnknowns() { return this->unknowns; } // rch
    int giveEqn() {return equationNumber;}

#ifdef __PARALLEL_MODE
    /**
     * Packs specific  DOF Manager's dofs unknowns into communication buffer.
     * @param buff communication buffer to pack data.
     * @param type physical meaning of  unknown.
     * @param mode mode of unknown (e.g, total value, velocity or acceleration of unknown).
     * @stepN time step when unknown requested. See documentation of particular EngngModel
     * class for valid StepN values (most implementaion can return only values for current
     * and possibly for previous time step).
     * @return nonzero if succesfull
     */
    virtual int packUnknowns(CommunicationBuffer &buff, EquationID type, ValueModeType mode, TimeStep *stepN);
    /**
     * Unpacks DOF unknown from communication buffer and updates unknown if necessary.
     * Unknown is always updated using EngngModel::updateUnknownComponent, if DOFManager
     * to which receiver belongs has DofManager_shared dofManagerParallelMode type.
     * Unknown is unpacked and stored in unknowns dictionary, if DOFManager
     * to which receiver belongs has DofManager_remote dofManagerParallelMode type.
     * There is no reason for invoking this service if  DOFManager has DofManager_local mode.
     * If dof belonging to shared or remote DofManager, engng model unknowns are updated to
     * accommodate remote contribution or "prescribed" remote values.
     * The unknown dictionary is not updated in this case, even if unknown dictionary should be used,
     * This is engng model job to update all unknowns dictionaries.
     * @param buff buffer containing packed message
     * @param type physical meaning of  unknown.
     * @param mode mode of unknown (e.g, total value, velocity or acceleration of unknown).
     * @stepN time step when unknown requested. See documentation of particular EngngModel
     * class for valid StepN values (most implementaion can return only values for current
     * and possibly for previous time step).
     * @return nonzero if succesfull
     * @see Dof::unpackAndUpdateDOFsUnknown for description
     */
    virtual int unpackAndUpdateUnknown(CommunicationBuffer &buff, EquationID type,
                                       ValueModeType mode, TimeStep *stepN);
#endif

protected:
    /**
     * Returns boundary condition of dof if it is prescribed.
     * @return returns NULL if no BC applied, otherwise pointer to correcpondig BC.
     */
    BoundaryCondition *giveBc();
    /**
     * Returns initial condition of dof if it is prescribed.
     * @return returns NULL if no IC applied, otherwise pointer to correcpondig IC.
     */
    InitialCondition *giveIc();
};

} // end namespace oofem
#endif // masterdof_h
