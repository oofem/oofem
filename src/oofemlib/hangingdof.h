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

/*
 * Contribution by Ladislav Svoboda
 */



#ifndef hangingdof_h
#define hangingdof_h

#include "dof.h"
#include "compiler.h"
#include "dictionr.h"
#include "cltypes.h"
#include "error.h"
#include "flotmtrx.h"

#ifndef __MAKEDEPEND
#include <stdio.h>
#include <string.h>
#endif

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

class HangingDof : public Dof
{
public:
    HangingDof(int i, DofManager *aNode, DofID id);
    ~HangingDof()   { }

    double giveUnknown(EquationID type, ValueModeType mode, TimeStep *stepN);
    double giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *stepN);
    double giveMasterUnknowns(FloatArray &masterUnknowns);

    /**
     * Returns the value of the unknown associated with the receiver
     * at given time step. Slave simply asks necessary master dofs and
     * computes the results.
     * @see MasterDof::giveUnknown function
     */
    double giveLocalUnknown(EquationID, ValueModeType, TimeStep *) { _error("HangingDof :: giveLocalUnknown: local coordinate system doesn't exist");
                                                                     return 0.0; }

    /**
     * Returns equation number corresponding to receiver.
     * Rigid Arm Slave have equation number undefined.
     * Usually single dof in node connected using rigid arm is
     * contributing to several master dofs (diplacement to displacement and rotations in master).
     * @return prints error msg and exits.
     */
    int giveEquationNumber(void) { _error("giveEquationNumber: undefined");
                                   return 0; }

    /**
     * Returns equation number corresponding to receiver.
     * Rigid Arm Slave have equation number undefined.
     * Usually single dof in node connected using rigid arm is
     * contributing to several master dofs (diplacement to displacement and rotations in master).
     * @return prints error msg and exits.
     */
    int givePrescribedEquationNumber() {
        //  _error ("givePrescribedEquationNumber: undefined");                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        //  v structengngmodel.C na line 195 ja tato fce volana
        // budes si muset rozmyslet Borku co ma tato fce vracet

        return 0;
    }

    /**
     * Asks new equation number. Empty function (master is assumed to receive same message).
     */
    int askNewEquationNumber(TimeStep *tStep) { return 1; }

    /**
     * Returns boundary condition of dof if it is precsribed.
     * HangingDof can not be subjected to bc, it is only mapping to master
     * @return returns NULL if no BC applied, otherwise pointer to correcpondig BC.
     */
    int hasBc(TimeStep *tStep) { return 0; }

    /**
     * Returns initial condition of dof if it is precsribed.
     * HangingDof can not be subjected to ic, it is only mapping to master
     * @see MasterDof::hasIc
     */
    int hasIc() { return 0; }

    /**
     * RigidArmSlaveDof can not be subjected to ic - it is only mapping to master.
     * @see MasterDof::hasIc
     */
    int hasIcOn(ValueModeType) { return 0; }

    /**
     * Returns the id of associated boundary condition, if there is any.
     * Used only for printing purposes. In general, id culd not be used
     * to decide whether bc is active. Use appropriate services instead.
     * @param id of associated Boubdaray condition, zero otherwise
     */
    int giveBcIdValue() { return 0; }

    /**
     * Stores receiver state to output stream.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    contextIOResultType saveContext(FILE *stream, void *obj = NULL) { return CIO_OK; }

    /**
     * Restores the receiver state previously written in stream.
     * @exception throws an ContextIOERR exception if error encountered.
     */
    contextIOResultType restoreContext(FILE *stream, void *obj = NULL) { return CIO_OK; }

    /**
     * Returns class name of the receiver.
     */
    const char *giveClassName() const { return "HangingDof"; }

    /**
     * Returns classType id of receiver.
     */
    classType giveClassID() const { return HangingDofClass; }
};

#endif // hangingdof_h
