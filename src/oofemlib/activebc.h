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

#ifndef activebc_h
#define activebc_h

#include "generalbc.h"
#include "alist.h"
#include "intarray.h"
#include "equationid.h"
#include "chartype.h"
#include "valuemodetype.h"
#include "error.h"

namespace oofem {
class SparseMtrx;
class UnknownNumberingScheme;
class ActiveDof;
class Dof;
class PrimaryField;

/**
 * Abstract base class for all active boundary conditions.
 * Design of active boundary conditions are subject to change.
 */
class ActiveBoundaryCondition : public GeneralBoundaryCondition
{
public:
    /**
     * Constructor. Creates boundary an active condition with given number, belonging to given domain.
     * @param n Boundary condition number.
     * @param d Domain to which new object will belongs.
     */
    ActiveBoundaryCondition(int n, Domain *d) : GeneralBoundaryCondition(n, d) { }
    /// Destructor.
    virtual ~ActiveBoundaryCondition() { }

    IRResultType initializeFrom(InputRecord *ir)
    {
        GeneralBoundaryCondition :: initializeFrom(ir);

        const char *__proc = "initializeFrom";
        IRResultType result;
        IntArray tempA, tempB, tempC;
        IR_GIVE_OPTIONAL_FIELD(ir, tempA, IFT_ActiveBoundaryCondition_elements, "elements");
        for (int i = 0; i < tempA.giveSize(); ++i) {
            this->addElement(tempA(i));
        }
        IR_GIVE_OPTIONAL_FIELD(ir, tempB, IFT_ActiveBoundaryCondition_elementSides, "elementsides");
        for (int i = 0; i < tempB.giveSize()/2; ++i) {
            this->addElementSide(tempB(i*2),tempB(i*2+1));
        }
        IR_GIVE_OPTIONAL_FIELD(ir, tempC, IFT_ActiveBoundaryCondition_dofManagers, "dofmans");
        for (int i = 0; i < tempC.giveSize(); ++i) {
            this->addDofman(tempC(i));
        }

        return IRRT_OK;
    }

    ///  @name Methods supporting classical input files.
    //@{
    /**
     * Adds element for active boundary condition.
     * @param elem Element number.
     */
    virtual void addElement(int elem) { OOFEM_ERROR2("%s :: addElement - Not supported", giveClassName()); }

    /**
     * Adds element for active boundary condition.
     * @param elem Element number.
     * @param side Side number.
     */
    virtual void addElementSide(int elem, int side) { OOFEM_ERROR2("%s :: addElement - Not supported", giveClassName()); }

    /**
     * Adds dof manager for active boundary condition.
     * @param dman Dof manager number.
     */
    virtual void addDofman(int dman) { OOFEM_ERROR2("%s :: addElement - Not supported", giveClassName()); }
    //@}

    /**
     * Assembles B.C. contributions to specified matrix.
     * @param[in,out] answer Matrix to assemble to.
     * @param tStep Active time step.
     * @param tStep Active time step.
     * @param eid Equation ID.
     * @param type Type of matrix to assemble.
     * @param r_s Row numbering scheme.
     * @param c_s Column numbering scheme.
     * @param domain Domain to assemble from.
     */
    virtual void assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, Domain *domain) {};

    /**
     * Assembles B.C. contributions to specified vector.
     * @param[in,out] answer Vector to assemble to.
     * @param tStep Active time step.
     * @param eid Equation ID.
     * @param type Type of matrix to assemble.
     * @param mode Mode of value.
     * @param s Numbering scheme.
     * @param domain Domain to assemble from.
     * @return Equivalent of the sum of the element norm (squared) of assembled vector.
     */
    virtual double assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                  CharType type, ValueModeType mode,
                                  const UnknownNumberingScheme &s, Domain *domain) { return 0.0; };

    /**
     * Gives a list of location arrays that will be assembled.
     * This should only be used to construct zero structure in sparse matrices.
     * @param rows List of location arrays for r_s.
     * @param cols List of location arrays for c_s.
     * @param eid Equation ID.
     * @param type Type of matrix to assemble.
     * @param r_s Row numbering scheme.
     * @param c_s Column numbering scheme.
     * @param domain Domain to assemble from.
     */
    virtual void giveLocationArrays(AList<IntArray> &rows, AList<IntArray> &cols, EquationID eid, CharType type,
                                    const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, Domain *domain) {};


    /// @name Functions related to boundary conditions which have to deal with special DOFs.
    //@{
    /**
     * NOT ACTUALLY USED YET due to some design difficulties.
     * Checks to see if active boundary condition requires special DOFs.
     * @return True if ActiveDof should be created.
     */
    virtual bool requiresActiveDofs() { return false; }
    /**
     * Checks to see if the dof is a primary DOF.
     * @return True if ActiveDof is a primary DOF.
     */
    virtual bool isPrimaryDof(ActiveDof *dof) { return true; }
    /**
     * Returns the prescribed value of a dof (if any).
     */
    virtual double giveBcValue(ActiveDof *dof, ValueModeType mode, TimeStep *tStep) { return 0.0; }
    /**
     * Returns the prescribed value of a dof (if any).
     * @return True if dof is prescribed.
     */
    virtual bool hasBc(ActiveDof *dof, TimeStep *tStep) { return false; }
    /**
     * Allows for active boundary conditions to handle their own special DOF.
     * @param dof Active dof belonging to receiver.
     * @return Number of primary master DOFs.
     */
    virtual int giveNumberOfMasterDofs(ActiveDof *dof)
    {
        OOFEM_ERROR2("%s :: giveNumberOfPrimaryMasterDofs - Not supported by bc.", giveClassName());
        return 0;
    }
    /**
     * Give the pointer to master dof belonging to active DOF.
     * @param dof Active dof belonging to receiver.
     * @param mdof Local master dof number.
     * @return Master dof.
     */
    virtual Dof *giveMasterDof(ActiveDof *dof, int mdof)
    {
        OOFEM_ERROR2("%s :: giveMasterDof - Not supported by bc.", giveClassName());
        return NULL;
    }
    virtual void computeDofTransformation(ActiveDof *dof, FloatArray &masterContribs)
    {
        OOFEM_ERROR2("%s :: computeDofTransformation - Not supported by bc.", giveClassName());
    }
    /**
     * Computes the value of the dof.
     * @param field Field to take value from.
     * @param mode Mode of unknown value.
     * @param tStep Time step.
     * @return Value of dof.
     */
    virtual double giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *tStep, ActiveDof *dof)
    {
        OOFEM_ERROR2("%s :: giveUnknown - Not supported by bc.", giveClassName());
        return 0.0;
    }
    /**
     * Computes the value of the dof.
     * @param eid Equation ID for the unknown value.
     * @param mode Mode of unknown value.
     * @param tStep Time step.
     * @return Value of dof.
     */
    virtual double giveUnknown(EquationID type, ValueModeType mode, TimeStep *tStep, ActiveDof *dof)
    {
        OOFEM_ERROR2("%s :: giveUnknown - Not supported by bc.", giveClassName());
        return 0.0;
    }
    //@}

    virtual classType giveClassID() const { return ActiveBoundaryConditionClass; }
    virtual const char *giveClassName() const { return "ActiveBoundaryCondition"; }
};
} // end namespace oofem
#endif // activebc_h

