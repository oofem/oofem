/* $Header: /home/cvs/bp/oofem/oofemlib/src/slavedof.h,v 1.11 2003/04/06 14:08:25 bp Exp $ */
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
//   *** CLASS SLAVE DOF ***
//   ************************


#ifndef simpleslavedof_h
#define simpleslavedof_h

#include "dof.h"
#include "compiler.h"
#include "dictionr.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#include <string.h>
#endif

#include "error.h"

class Domain;
class DofManager;
class TimeStep;
class BoundaryCondition;
class InitialCondition;

/**
 * Class representing "slave" degree of freedom. Slave dof is linked to some master dof
 * (link slave-slave is not allowed). Slaves can be used to implement duplicated joints
 * (by specifying some dof in node as masters and some as slaves linked to dofs in
 * other node (with same or possibly different coordinates).
 * Slave dof share the same equation number with master. Almost all operation as request
 * for BC or IC or equation number are simply forwarded to master.
 * Other functions (which change internal state) like updateYourself, updateUnknownsDictionary,
 * askNewEquationNumber or context storing/restoring functions are empty functions,
 * relying on fact that same funtion will be called also for master.
 * From this point of wiev, one can see slave dof as link to other dof.
 *
 * It is important to ensure (on input) that both slave dofManager and master dofManager
 * are using the same local coordinate system. In future releases, this can be checked
 * using checkConsistency function, where this check could be performed.
 */
class SimpleSlaveDof : public Dof
{
    /*
     * This class implements a "slave" nodal degree of freedom.
     * A dof is usually attribute of one node.
     * DESCRIPTION
     *
     * Slave dof is special dof - connected to some master dof on other node (side), with
     * same equation number and same boundary and initial conditions. Typically, slave dof's node
     * and master node are sharing the same space position. This will allow to have multiple dof's
     * (i.e. displacement in x-dir, or rotations) in such node.
     *
     * It is simple wrapper, it translates queries to master.
     * So it does not needs to have its equation number, boundary or initial conditions -
     * all is queried from master. Also unknown dictionary is not necessary -
     * this is done in  master .
     *
     * 'number' and 'node' are used for reading/writing data in the data file.
     * 'DofID' is parameter determining physical meaning of receiver.
     * This parameter is also used in member function giveUnknownType, which returns
     * CharType type according to DofID parameter.
     *
     *
     * TASKS
     * - equation numbering, in method 'giveEquationNumber' ;
     * - managing its b.c. and its i.c., if any (methods 'hasBc', 'giveBc', etc);
     * - managing its unknowns. This includes retrieving the associated solution
     * from the Engng. System , or from receiver's dictionary
     * (based on emodel->requiresNodeUnknowsDictionaryUpdate() function,
     * which determines whether to use dictionary or ask unknowns values from
     * emodel)
     * - managing physical meaning of dof (dofID variable)
     *
     * REMARKS
     * - class Dof is not a subclass of FEMComponent : a dof belongs to a single
     * node, not to the domain ;
     * - class Dof is not restricted to structural analysis problems. Unknowns
     * may also be pressures, temperatures, etc.
     * - method give returns unknown value quantity according to ValueModeType parameter,
     * CharType parameter is used to check whether pysical meaning of
     * unknown corresponds.
     *
     */

private:
    /// number of DofManager containing master dof (Master DofManager)
    int masterDofMngr;
    /// number of master dof in master dofManager.
    int masterDofIndx;
    /*      Dictionary*  unknowns ;
     *      Dictionary*  pastUnknowns ; */

public:
    /**
     * Constructor. Creates slave dof vith number i, belonging to aNode dof manager.
     * Slave will be linked to master dof with id type belonging to dofManager with
     * number given in master variable.
     * @param i dof number
     * @param aNode receiver will belong to aNode dof manager
     * @param master number of dofManager which contain master dof
     * @param id DofID of master dof (and slave too).
     */
    SimpleSlaveDof(int i, DofManager *aNode, int master, DofID id);  // constructor
    SimpleSlaveDof(int i, DofManager *aNode, DofID id = Undef);   // constructor

    /// Destructor.
    ~SimpleSlaveDof()   { }   // destructor.

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "SimpleSlaveDof"; }
    /// Returns classType id of receiver.
    classType           giveClassID() const { return SimpleSlaveDofClass; }
    /**
     * Returns equation number corresponding to receiver.
     * Slave simply forwards this mesage to master.
     * @return equation number, if active BC exists, returns zero
     */
    int                 giveEquationNumber();
    /**
     * Returns prescribed equation number corresponding to receiver.
     * Slave simply forwards this mesage to master.
     * @return prescribed equation number, if active BC exists, returns zero
     */
    int                 givePrescribedEquationNumber();
    /**
     * Asks new equation number. Empty function (master is assumed to receive same message).
     */
    int                 askNewEquationNumber(TimeStep *tStep) { return 1; }
    /**
     * Returns the value of the unknown associated with the receiver
     * at given time step. Slave simply forwards this message to master dof.
     * @see MasterDof::giveUnknown function
     */
    double giveUnknown(EquationID, ValueModeType, TimeStep *);
    /**
     * Returns the value of the unknown associated to given field of the receiver
     * at given time step. Slave simply forwards this message to master dof.
     * @see MasterDof::giveUnknown function
     */
    double giveUnknown(PrimaryField & field, ValueModeType, TimeStep * stepN);
    /**
     * Returns boundary condition of dof if it is precsribed.
     * Slave simply forwards this message to master.
     * @return returns NULL if no BC applied, otherwise pointer to correcpondig BC.
     */
    int                 hasBc(TimeStep *tStep);
    /**
     * Slave simply forwards this message to master.
     * @see MasterDof::hasIc
     */

    int                 hasIc();
    /**
     * Slave simply forwards this message to master.
     * @see MasterDof::hasIc
     */
    int hasIcOn(ValueModeType);
    /** Returns the id of associated boundary condition, if there is any.
     * Used only for printing purposes. In general, id could not be used
     * to decide whether bc is active. Use appropriate services instead.
     * @param id of associated Boundary condition, zero otherwise
     */
    int giveBcId () ;
    /** Returns the id of associated initial condition, if there is any.
     * Used only for printing purposes. In general, id could not be used
     * to decide whether bc is active. Use appropriate services instead.
     * @param id of associated initial condition, zero otherwise
     */
    int giveIcId () ;


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

    /// Returns Master Dof Manager Number.
    int giveMasterDofManagerNum() const { return masterDofMngr; }
    /// Sets masterDofMngr
    void setMasterDofManagerNum(int i) { masterDofMngr = i; }
    /// Returns number of master dof in master dofManager.
    int giveMasterDofIndx() const { return masterDofIndx; }

    /**
     * Local renumbering support. For some tasks (parallel load balancing, for example) it is necessary to
     * renumber the entities. The various fem components (such as nodes or elements) typically contain
     * links to other entities in terms of their local numbers, etc. This service allows to update
     * these relations to reflext updated numbering. The renumbering funciton is passed, which is supposed
     * to return an updated number of specified entyty type based on old number.
     */
    virtual void updateLocalNumbering( EntityRenumberingFunctor &f );


private:
    /// returns reference to master dof.
    Dof *giveMasterDof();


protected:
    /**
     * Returns boundary condition of dof if it is precsribed.
     * Slave simply forwards this mesage to master.
     * @return returns NULL if no BC applied, otherwise pointer to correcpondig BC.
     */
    BoundaryCondition *giveBc();
    /**
     * Returns initial condition of dof if it is precsribed.
     * Slave forwards this message to master.
     * @return returns NULL if no IC applied, otherwise pointer to correcpondig IC.
     */
    InitialCondition *giveIc();
};


#endif // simpleslavedof_h
