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

#ifndef simpleslavedof_h
#define simpleslavedof_h

#include "dof.h"
#include "error.h"

namespace oofem {
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
 * relying on fact that same function will be called also for master.
 * From this point of view, one can see slave dof as link to other dof.
 *
 * It is important to ensure (on input) that both slave dofManager and master dofManager
 * are using the same local coordinate system. In future releases, this can be checked
 * using checkConsistency function, where this check could be performed.
 *
 * Slave dof is special dof - connected to some master dof on other node (side), with
 * same equation number and same boundary and initial conditions. Typically, slave dof's node
 * and master node are sharing the same space position. This will allow to have multiple dof's
 * (i.e. displacement in x-dir, or rotations) in such node.
 *
 */
class OOFEM_EXPORT SimpleSlaveDof : public Dof
{
private:
    /// Number of DofManager containing master dof (Master DofManager)
    int masterDofMngr;
    /// Number of master dof in master dofManager.
    mutable int masterDofIndx;

public:
    /**
     * Constructor. Creates slave dof vith number i, belonging to aNode dof manager.
     * Slave will be linked to master dof with id type belonging to dofManager with
     * number given in master variable.
     * @param aNode receiver will belong to aNode dof manager
     * @param master number of dofManager which contain master dof
     * @param id DofID of master dof (and slave too).
     */
    SimpleSlaveDof(DofManager * aNode, int master, DofIDItem id);
    SimpleSlaveDof(DofManager * aNode, DofIDItem id = Undef);
    /// Destructor.
    virtual ~SimpleSlaveDof() { }

    dofType giveDofType() override { return DT_simpleSlave; }
    const char *giveClassName() const override { return "SimpleSlaveDof"; }
    /**
     * Returns equation number corresponding to receiver.
     * Slave simply forwards this message to master.
     * @return Equation number, if active BC exists, returns zero
     */
    int __giveEquationNumber() const override;
    /**
     * Returns prescribed equation number corresponding to receiver.
     * Slave simply forwards this message to master.
     * @return Prescribed equation number, if active BC exists, returns zero
     */
    int __givePrescribedEquationNumber() override;
    /**
     * Asks new equation number. Empty function (master is assumed to receive same message).
     */
    int askNewEquationNumber(TimeStep *tStep) override { return 1; }
    double giveUnknown(ValueModeType, TimeStep *) override;
    double giveUnknown(PrimaryField &field, ValueModeType, TimeStep *tStep) override;
    bool hasBc(TimeStep *tStep) override;
    bool hasIc() override;
    bool hasIcOn(ValueModeType) override;
    int giveBcId() override;
    int giveIcId() override;
    double giveBcValue(ValueModeType mode, TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    /// Returns Master Dof Manager Number.
    int giveMasterDofManagerNum() const { return masterDofMngr; }
    void giveMasterDofManArray(IntArray &answer) override {
        answer.resize(1);
        answer.at(1) = masterDofMngr;
    }
    /// Sets master dof manager
    void setMasterDofManagerNum(int i) { masterDofMngr = i; }
    /// Returns number of master dof in master dofManager.
    int giveMasterDofIndx() const { return masterDofIndx; }

    void updateLocalNumbering(EntityRenumberingFunctor &f) override;

private:
    /// Returns reference to master dof.
    Dof *giveMasterDof() const;

protected:
    BoundaryCondition *giveBc() override;
    InitialCondition *giveIc() override;
};
} // end namespace oofem
#endif // simpleslavedof_h
