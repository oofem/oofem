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

#ifndef slavedof_h
#define slavedof_h

#include "dof.h"

namespace oofem {
/**
 * Class representing "slave" degree of freedom. This dof is generally linked to some master dofs
 * with various weights - contributions - combination (link slave-slave is allowed).
 * The slave dof is fully dependent upon master dofs, therefore it's equation number is undefined.
 * It can have own boundary conditions - it is not yet implemented.
 */
class SlaveDof : public Dof
{
protected:
    /// Count of master DofManagers.
    int countOfMasterDofs;
    /// Count of master DofManagers.
    int countOfPrimaryMasterDofs;
    /// Array of master  DofManagers.
    IntArray masterDofMans;
    /// Array of master dofIDs.
    IntArray dofIDs;
    /// Vector of master contribution coefficients.
    FloatArray masterContribution;

public:
    /**
     * Constructor. Creates slave dof with number n, belonging to aNode dof manager.
     * @param n Dof number.
     * @param aNode Node receiver will belong to.
     * @param id DofID of slave dof.
     */
    SlaveDof(int n, DofManager *aNode, DofIDItem id = Undef);
    /// Destructor.
    virtual ~SlaveDof(void) { }

    void initialize(int cntOfMstrDfMngr, const IntArray &masterNodes, const IntArray *mstrDofID, const FloatArray &mstrContribution);
    virtual int giveNumberOfPrimaryMasterDofs();
    virtual void giveMasterDofManArray(IntArray &answer);
    virtual void giveUnknowns(FloatArray &masterUnknowns, EquationID type, ValueModeType mode, TimeStep *stepN);
    virtual void giveUnknowns(FloatArray &masterUnknowns, PrimaryField &field, ValueModeType mode, TimeStep *stepN);
    virtual void computeDofTransformation(FloatArray &primaryMasterContribs);
    virtual void giveEquationNumbers(IntArray &masterEqNumbers, const UnknownNumberingScheme &s);

    /**
     * Returns the value of the unknown associated with the receiver at given time step.
     * Slave simply asks vector of corresponding master dofs and own transformation
     * vector and returns result as dot product of these vectors. Standard element
     * services have to transform global unknown vector transform into their local c.s
     * before using it (when computing strain vector by @f$\epsilon = B\cdot r @f$, for example,
     * where B is element geometrical matrix). This transformation should contain also
     * nodal to global coordinate system transformation. So, this specialized
     * standard method for unknown query returns the corresponding master DOF value.
     * @see MasterDof::giveUnknown
     */
    virtual double giveUnknown(EquationID type, ValueModeType mode, TimeStep *stepN);
    virtual double giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *stepN);

    /**
     * Returns equation number corresponding to receiver.
     * Rigid Arm Slave have equation number undefined.
     * Usually single dof in node connected using rigid arm is
     * contributing to several master dofs (displacement to displacement and rotations in master).
     * @return Prints error message and exits.
     */
    virtual int __giveEquationNumber() const {
        _error("giveEquationNumber: undefined");
        return 0;
    }

    /**
     * Returns equation number corresponding to receiver.
     * Rigid Arm Slave have equation number undefined.
     * Usually single dof in node connected using rigid arm is
     * contributing to several master dofs (displacement to displacement and rotations in master).
     * @return Prints error message and exits.
     */
    virtual int __givePrescribedEquationNumber() {
        _error("givePrescribedEquationNumber: undefined");
        return 0;
    }
    /**
     * Asks new equation number. Empty function (master is assumed to receive same message).
     */
    virtual int askNewEquationNumber(TimeStep *tStep) { return 1; }

    /**
     * Returns boundary condition of dof if it is prescribed.
     * HangingDof can not be subjected to BC, it is only mapping to master
     * @return NULL if no BC applied, otherwise pointer to corresponding BC.
     */
    virtual bool hasBc(TimeStep *tStep) { return false; }

    /**
     * Returns initial condition of dof if it is prescribed.
     * HangingDof can not be subjected to IC, it is only mapping to master
     * @see MasterDof::hasIc
     */
    virtual bool hasIc() { return false; }

    /**
     * RigidArmSlaveDof can not be subjected to IC - it is only mapping to master.
     * @see MasterDof::hasIc
     */
    virtual bool hasIcOn(ValueModeType) { return false; }

    virtual int giveBcId() { return 0; }
    virtual int giveIcId() { return 0; }

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    virtual const char *giveClassName() const { return "SlaveDof"; }
    virtual classType giveClassID() const { return SlaveDofClass; }

    virtual void updateLocalNumbering(EntityRenumberingFunctor &f);

protected:
    inline Dof *giveMasterDof(int i);
};
} // end namespace oofem
#endif // slavedof_h
