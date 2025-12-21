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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#ifndef activedof_h
#define activedof_h

#include "dof.h"

namespace oofem {
class ActiveBoundaryCondition;

/**
 * Class representing "slave" degree of freedom with an active boundary condition.
 * It is similar to SlaveDof, but its actual value is controlled by an active boundary condition.
 * The code is simple, the functions just pass on the evaluation to the corresponding active boundary condition.
 */
class OOFEM_EXPORT ActiveDof : public Dof
{
protected:
    /// Corresponding equation number (positive value) or prescribed equation number (negative value).
    int equationNumber;
    /// Boundary condition number.
    int bc;
    /// Active boundary condition number.
    ActiveBoundaryCondition *activeBC;

public:
    /**
     * Constructor. Creates slave dof with number n, belonging to aNode dof manager.
     * @param aNode Node receiver will belong to.
     * @param id DofID of slave dof.
     * @param bc Boundary condition dof belongs to.
     */
    ActiveDof(DofManager * aNode, DofIDItem id = Undef, int bc = 0);
    /// Destructor.
    virtual ~ActiveDof() { }

    void initialize(int cntOfMstrDfMngr, const IntArray &masterNodes, const IntArray *mstrDofID, const FloatArray &mstrContribution);
    int giveNumberOfPrimaryMasterDofs() override;
    bool isPrimaryDof() override;
    int giveNumberOfMasterDofs();
    void giveMasterDofManArray(IntArray &answer) override;
    void giveUnknowns(FloatArray &masterUnknowns, ValueModeType mode, TimeStep *tStep) override;
    void giveUnknowns(FloatArray &masterUnknowns, PrimaryField &field, ValueModeType mode, TimeStep *tStep) override;
    void computeDofTransformation(FloatArray &primaryMasterContribs) override;
    void giveEquationNumbers(IntArray &masterEqNumbers, const UnknownNumberingScheme &s) override;
    void giveDofIDs(IntArray &masterDofIDs) override;

    double giveUnknown(ValueModeType mode, TimeStep *tStep) override;
    double giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *tStep) override;

    void saveContext(DataStream &stream, ContextMode mode) override;
    void restoreContext(DataStream &stream, ContextMode mode) override;

    dofType giveDofType() override { return DT_active; }
    const char *giveClassName() const override { return "ActiveDof"; }

    void updateLocalNumbering(EntityRenumberingFunctor &f) override;

    int __giveEquationNumber() const override;
    int __givePrescribedEquationNumber() override;
    int askNewEquationNumber(TimeStep *tStep) override;
    bool hasBc(TimeStep *tStep) override;
    int giveBcId() override;
    void setBcId(int bcId) override;
    double giveBcValue(ValueModeType mode, TimeStep *tStep) override;

    bool hasIcOn(ValueModeType type) override;
    InitialCondition *giveIc() override;
    bool hasIc() override;
    int giveIcId() override;

    ActiveBoundaryCondition *giveActiveBoundaryCondition();

protected:
    inline Dof *giveMasterDof(int i);
};
} // end namespace oofem
#endif // slavedof_h
