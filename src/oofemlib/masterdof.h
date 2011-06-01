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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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
 *
 * This class implements a nodal degree of freedom. A dof is usually attribute of one node.

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
 * its unknowns (e.g., the displacement 'd', the velocity 'v' and the acceleration 'a'),
 * at the current time step and at the previous one.
 *
 * TASKS
 * - Equation numbering, in method 'giveEquationNumber' and 'givePrescribedEquationNumber' ;
 * - Managing its b.c. and its i.c., if any (methods 'hasBc', 'giveBc', etc);
 * - Managing its unknowns. This includes retrieving the associated solution
 *   from the engineering model, or from receiver's dictionary
 *   (based on emodel->requiresNodeUnknowsDictionaryUpdate() function,
 *   which determines whether to use dictionary or ask unknowns values from
 *   emodel)
 * - Managing physical meaning of dof (dofID variable)
 *
 * REMARKS
 * - Class Dof is not a subclass of FEMComponent : a dof belongs to a single
 *   node, not to the domain ;
 * - Class Dof is not restricted to structural analysis problems. Unknowns
 *   may also be pressures, temperatures, etc.
 * - Method give returns unknown value quantity according to ValueModeType parameter,
 *   CharType parameter is used to check whether physical meaning of
 *   unknown corresponds.
 *
 */
class MasterDof  : public Dof
{
protected:
    /// Corresponding equation number (positive value) or prescribed equation number (negative value).
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
     * @see cltypes.h, DofIDItem type 
     */
    MasterDof(int i, DofManager *aNode, int nbc, int nic, DofIDItem id);
    MasterDof(int i, DofManager *aNode, DofIDItem id = Undef);
    /// Destructor.
    ~MasterDof()   { delete unknowns; /*delete pastUnknowns ;*/ }

    const char *giveClassName() const { return "MasterDof"; }
    classType giveClassID() const { return MasterDofClass; }

    virtual int __giveEquationNumber() const ;

    virtual int __givePrescribedEquationNumber();

    int askNewEquationNumber(TimeStep *tStep);

    double giveUnknown(EquationID eid, ValueModeType mode, TimeStep *stepN);
    double giveUnknown(PrimaryField & field, ValueModeType, TimeStep *stepN);

    bool hasBc(TimeStep *tStep);
    bool hasIc();
    bool hasIcOn(ValueModeType);

    bool isPrimaryDof() { return true; }

    int giveBcId();
    int giveIcId();

    void printYourself();
    void updateYourself(TimeStep *);

    void updateUnknownsDictionary(TimeStep *tStep, EquationID type, ValueModeType mode, double dofValue);

    virtual void giveUnknownsDictionaryValue(TimeStep *tStep, EquationID type, ValueModeType mode, double &dofValue);

    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

    void setBcId(int bcId) { this->bc = bcId; }
    void setEquationNumber(int equationNumber) { this->equationNumber = equationNumber; }
    void setUnknowns(Dictionary *unknowns) { this->unknowns = unknowns; }
    Dictionary *giveUnknowns() { return this->unknowns; }
    int giveEqn() { return equationNumber; }

#ifdef __PARALLEL_MODE
    virtual int packUnknowns(CommunicationBuffer &buff, EquationID type, ValueModeType mode, TimeStep *stepN);
    virtual int unpackAndUpdateUnknown(CommunicationBuffer &buff, EquationID type, ValueModeType mode, TimeStep *stepN);
#endif

protected:
    BoundaryCondition *giveBc();
    InitialCondition *giveIc();
};
} // end namespace oofem
#endif // masterdof_h
