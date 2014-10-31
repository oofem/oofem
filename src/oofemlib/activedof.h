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

    virtual void initialize(int cntOfMstrDfMngr, const IntArray &masterNodes, const IntArray *mstrDofID, const FloatArray &mstrContribution);
    virtual int giveNumberOfPrimaryMasterDofs();
    virtual bool isPrimaryDof();
    int giveNumberOfMasterDofs();
    virtual void giveMasterDofManArray(IntArray &answer);
    virtual void giveUnknowns(FloatArray &masterUnknowns, ValueModeType mode, TimeStep *tStep);
    virtual void giveUnknowns(FloatArray &masterUnknowns, PrimaryField &field, ValueModeType mode, TimeStep *tStep);
    virtual void computeDofTransformation(FloatArray &primaryMasterContribs);
    virtual void giveEquationNumbers(IntArray &masterEqNumbers, const UnknownNumberingScheme &s);
    virtual void giveDofIDs(IntArray &masterDofIDs);

    virtual double giveUnknown(ValueModeType mode, TimeStep *tStep);
    virtual double giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *tStep);

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    virtual dofType giveDofType() { return DT_active; }
    virtual const char *giveClassName() const { return "ActiveDof"; }

    virtual void updateLocalNumbering(EntityRenumberingFunctor &f);

    virtual int __giveEquationNumber() const;
    virtual int __givePrescribedEquationNumber();
    virtual int askNewEquationNumber(TimeStep *tStep);
    virtual bool hasBc(TimeStep *tStep);
    virtual int giveBcId();
    virtual void setBcId(int bcId);
    virtual double giveBcValue(ValueModeType mode, TimeStep *tStep);

    virtual bool hasIc(TimeStep *tStep);
    virtual bool hasIcOn(ValueModeType type);
    virtual InitialCondition *giveIc();
    virtual bool hasIc();
    virtual int giveIcId();

    ActiveBoundaryCondition *giveActiveBoundaryCondition();

protected:
    inline Dof *giveMasterDof(int i);
};
} // end namespace oofem
#endif // slavedof_h
