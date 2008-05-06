/* $Header: /home/cvs/bp/oofem/oofemlib/src/rigidarmslavedof.h,v 1.13 2003/04/06 14:08:25 bp Exp $ */
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


//   *********************************
//   *** CLASS RIGID ARM SLAVE DOF ***
//   *********************************


#ifndef rigidarmslavedof_h
#define rigidarmslavedof_h

#include "dof.h"
#include "compiler.h"
#include "dictionr.h"
#ifndef __MAKEDEPEND
#include <stdio.h>
#include <string.h>
#endif

#include "error.h"
#include "flotmtrx.h"

class Domain;
class DofManager;
class TimeStep;
class BoundaryCondition;
class InitialCondition;

/**
 * Class representing "rigid arm slave" degree of freedom. This dof is linked to some master dof
 * (link slave-slave is not allowed) using rigid arm.
 * Implemented using nodal transformation. The Rigid Arm Slave dof represent dof, which is directly
 * related to master dofs. Therefore the rigid arm slave's equation number is undefined.
 * Similarly, rigid arm slave cannot have own boundary or initial conditions - these are
 * entirely detrmined using master boundary or initial conditions.
 *
 * The transformation for DOFs and load is not ortogonal - the inverse transformation can
 * not be constructed by transposition. Because of time consuming inversion, methods
 * can generally compute both transformations for dofs as well as loads.
 */
class RigidArmSlaveDof : public Dof
{
    /*
     * This class implements a "rigid arm slave" nodal degree of freedom.
     * A dof is usually attribute of one node.
     * DESCRIPTION
     *
     * Rigid Arm Slave dof is special dof - connected using rigid arm to some master dof on other node (side),
     * with given equation number and boundary and initial conditions.
     * Implemented using nodal transformation. The Rigid Arm Slave dof has no dofs, its unknowns are
     * related to master dofs.
     *
     */


public:
    enum RigidArmSlaveDofTransfType {
        _toMaster,
        _toSlave
    };

private:
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
    RigidArmSlaveDof(int i, DofManager *aNode, DofID id);        // constructor
    /// Destructor.
    ~RigidArmSlaveDof()   { }       // destructor.

    /// Returns class name of the receiver.
    const char *giveClassName() const { return "RigidArmSlaveDof"; }
    /// Returns classType id of receiver.
    classType           giveClassID() const { return RigidArmSlaveDofClass; }
    /**
     * Returns equation number corresponding to receiver.
     * Rigid Arm Slave have equation number undefined.
     * Usually single dof in node connected using rigid arm is
     * contributing to several master dofs (diplacement to displacement and rotations in master).
     * @return prints error msg and exits.
     */
    int                 giveEquationNumber() { _error("giveEquationNumber: undefined");
                                               return 0; }
    /**
     * Returns equation number corresponding to receiver.
     * Rigid Arm Slave have equation number undefined.
     * Usually single dof in node connected using rigid arm is
     * contributing to several master dofs (diplacement to displacement and rotations in master).
     * @return prints error msg and exits.
     */
    int                 givePrescribedEquationNumber() { _error("givePrescribedEquationNumber: undefined");
                                                         return 0; }
    /**
     * Asks new equation number. Empty function (master is assumed to receive same message).
     */
    int                 askNewEquationNumber(TimeStep *tStep) { return 1; }
    /**
     * Returns the value of the unknown associated with the receiver
     * at given time step. Slave simply asks corresponding master dof and
     * returns master result. Standart element services have to transform
     * global unknown vector transform into their local c.s before using it
     * (when computing strain vector by \eps=Br, for example,
     * where B is element geometrical matrix). This transformation should contain also
     * nodal to global coordinate system transformation. So, this specialized
     * standard method for unknown query returns the corresponding master DOF value.
     * @see MasterDof::giveUnknown function
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
    double              giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *stepN);
    /**
     * Returns the value of the unknown associated with the receiver
     * at given time step. Slave simply asks necessary master dofs and
     * computes the results.
     * @see MasterDof::giveUnknown function
     */
    double giveLocalUnknown(EquationID, ValueModeType, TimeStep *);
    /**
     * Returns boundary condition of dof if it is precsribed.
     * RA Slave can not be subjected to bc - it is a mapping of master dofs.
     * @return returns NULL if no BC applied, otherwise pointer to correcpondig BC.
     */
    int                 hasBc(TimeStep *tStep);
    /**
     * RA Slave can not be subjected to initial condition - it is mapping of master dofs.
     * @see MasterDof::hasIc
     */
    int                 hasIc();
    /**
     * RigidArmSlaveDof can not be subjected to ic - it is only mapping to master.
     * @see MasterDof::hasIc
     */
    int hasIcOn(ValueModeType);
    /** Returns the id of associated boundary condition, if there is any.
     * Used only for printing purposes. In general, id culd not be used
     * to decide whether bc is active. Use appropriate services instead.
     * @param id of associated Boubdaray condition, zero otherwise
     */
    int giveBcIdValue();

    /**
     * Stores receiver state to output stream.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    contextIOResultType    saveContext(FILE *stream, void *obj = NULL);
    /**
     * Restores the receiver state previously written in stream.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    contextIOResultType    restoreContext(FILE *stream, void *obj = NULL);

    /**
     * Computes DOF transformation mtrx of receiver. This transformation include
     * transformation to master dofs (in this case, all master dofs are considered to contribute).
     * This transformation is not orthogonal.
     * @param answer contains result.
     * @param ttype determines the "direction" of transformation.
     */
    void    computeDofTransformation(FloatMatrix &answer, RigidArmSlaveDofTransfType ttype);
    /**
     * Computes LOAD transformation mtrx of receiver. This transformation include
     * transformation to master dofs (in this case, all master dofs are considered to contribute).
     * This transformation is not orthogonal.
     * @param answer contains result.
     * @param ttype determines the "direction" of transformation.
     */
    void    computeLoadTransformation(FloatMatrix &answer, RigidArmSlaveDofTransfType ttype);
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
                                              EquationID type, IntArray &mode, int nite);
};

#endif // rigidarmslavedof_h
